// $(CC) -O3 -march=native -DMAIN -o dozeu dozeu.c
// #define DEBUG
// #define DZ_PRINT_VECTOR
/**
 * @file dozeu.h
 * @brief SIMD X-drop DP for read-to-graph alignment
 *
 * @author Hajime Suzuki
 * @date 18/03/20 at Kyoto
 * @license MIT
 *
 * @datail
 * minimum requirements:
 *   SSE4.1 (Core 2 / Bobcat or later)
 *   { gcc, clang, icc } x { Linux, Mac OS }
 * input sequence encoding to be either of the three (not required to be terminated with '\0'):
 *   1. {'A','C','G','T','U','a','c','g','t','u'}      (ASCII nucleotide)
 *   2. { 0, 1, 2, 3 }                                 (2-bit nucleotide)
 *   3. { 0, 1, 2, 3, ..., 19 }                        (integer protein)
 *
 * References:
 *   1. BLAST X-drop DP
 *     https://github.com/elucify/blast-docs/wiki/The-Developer's-Guide-to-BLAST (the most datailed document that describes the original X-drop DP algorithm)
 *   2. Myers bit vector (describes vertical (horizontal) tiling)
 *     Gene Myers, A fast bit-vector algorithm for approximate string matching based on dynamic programming, JACM (1999)
 */
#ifdef __cplusplus
extern "C" {
#endif

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* _MM_FROUND_CUR_DIRECTION is defined in smmintrin.h (the include
   file for SSE4.1), so if it's already defined we assume that the
   SSE4.1 API is available and we can skip including <smmintrin.h>.
   This allows Dozeu to use alternative implementations (e.g.,
   <https://github.com/nemequ/simde>) without depending on them. */
#if !defined(_MM_FROUND_CUR_DIRECTION)
#include <smmintrin.h>
#endif

#ifndef DZ_CIGAR_OP
#  define DZ_CIGAR_OP				0x04030201
#endif
#ifndef dz_cmp_max
#  define dz_cmp_max(x, y)			( (x) > (y) )
#endif

/* default encoding: ascii */
#ifndef DZ_PROTEIN
#  if !defined(DZ_NUCL_ASCII) && !defined(DZ_NUCL_2BIT)
#    define DZ_NUCL_ASCII
#  endif

/* define DZ_N_AS_UNMATCHING_BASE to penalize Ns, otherwise scores for (x, N) and (N, x) are always set zero */
// #define DZ_N_AS_UNMATCHING_BASE
enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, qN = 0x02, qS = 0x02
	#else
		rN = 0x90, qN = 0x90, qS = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};
#define dz_pair_score(_self, _q, _r, _i)	( (_self)->matrix[((_r) | (_q)->arr[(_i)]) & 0x1f] - DZ_SCORE_OFS )
#define dz_pair_eq(_self, _q, _r, _i)		( (uint32_t)((_q)->arr[(_i)]) == ((uint32_t)(_r)<<2) )
#else

/* protein */
#  ifndef DZ_MAT_SIZE
#    define DZ_MAT_SIZE				( 32 )
#  endif
#define dz_pair_score(_self, _q, _r, _i)	( (int8_t)((_q)->arr[(_r) * (_q)->blen * L + (_i)]) - DZ_SCORE_OFS )
#define dz_pair_eq(_self, _q, _r, _i)		( (uint32_t)((_q)->q[(_i) - 1] & 0x1f) == (uint32_t)(_r) )
#endif


/* utils */
#define dz_pp_cat_intl(x, y)		x##y
#define dz_pp_cat(x, y)				dz_pp_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_pp_cat(_st_, __LINE__)[(expr) ? 1 : -1]
#define dz_trap()					{ *((volatile uint8_t *)NULL); }
#ifdef DZ_PROTEIN
	dz_static_assert(DZ_MAT_SIZE <= 32);
#endif

#if (defined(DEBUG) || defined(UNITTEST)) && !defined(__cplusplus)
#  include "log.h"
#  define UNITTEST_ALIAS_MAIN		0
#  define UNITTEST_UNIQUE_ID		3213
#  include "unittest.h"
unittest_config( "dozeu" );
unittest() { debug("hello"); }
#else
#  define unittest(...)				static void dz_pp_cat(dz_unused_, __LINE__)(void)
#  define ut_assert(...)			;
#  define debug(...)				;
#  define trap()					;
#endif


#if defined(DZ_NUCL_ASCII)
#  define DZ_UNITTEST_INDEX			0
#elif defined(DZ_NUCL_2BIT)
#  define DZ_UNITTEST_INDEX			1
#elif defined(DZ_PROTEIN)
#  define DZ_UNITTEST_INDEX			2
#endif
#define dz_ut_sel(a, b, c)			( (DZ_UNITTEST_INDEX == 0) ? (a) : ((DZ_UNITTEST_INDEX == 1) ? (b) : (c)) )

/* vectorize */
#ifndef __x86_64__
#  error "x86_64 is required"
#endif
#define __dz_vectorize			/* follow the compiler options */

/* inlining (FIXME: add appropriate __force_inline flag) */
#define __dz_force_inline			inline
#ifdef __cplusplus
#  define dz_unused(x)
#else
#  define dz_unused(x)				(void)(x)
#endif
#define dz_likely(x)				__builtin_expect(!!(x), 1)
#define dz_unlikely(x)				__builtin_expect(!!(x), 0)
#define dz_roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )
#define dz_rounddown(x, base)		( (x) & ~((base) - 1) )
#define dz_max2(x, y)				( (x) < (y) ? (y) : (x) )
#define dz_min2(x, y)				( (x) < (y) ? (x) : (y) )
#define dz_inside(x, y, z)			( ((uint64_t)(y) - (uint64_t)(x)) < ((uint64_t)(z) - (uint64_t)(x)) )
#define dz_loadu_u64(p)				({ uint8_t const *_p = (uint8_t const *)(p); *((uint64_t const *)_p); })
#define dz_storeu_u64(p, e)			{ uint8_t *_p = (uint8_t *)(p); *((uint64_t *)(_p)) = (e); }

#define dz_add_ofs(_x)				( (uint16_t)(_x) ^ (uint16_t)0x8000 )
#define dz_rm_ofs(_x)				( (int16_t)((uint16_t)(_x) ^ (uint16_t)0x8000) )


#define dz_is_all_zero(x)			( _mm_test_all_zeros((x), (x)) == 1 )

#define DZ_MEM_MARGIN_SIZE			( 256 )
#define DZ_MEM_ALIGN_SIZE			( 16 )
#define DZ_MEM_INIT_SIZE			( 16 * 1024 * 1024 )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_CELL_MARGIN				( 32 )
#define DZ_CELL_MARGINED_MIN		( DZ_CELL_MIN + DZ_CELL_MARGIN )
#define DZ_CELL_MARGINED_MAX		( DZ_CELL_MAX - DZ_CELL_MARGIN )

#define DZ_SCORE_OFS				( 64 )


/* query; preconverted query sequence; blen = roundup(qlen, L) / L; array must have 16-byte-length margin at the tail */
struct dz_query_s { uint64_t blen; char const *q; int16_t bonus[2 * sizeof(__m128i) / sizeof(int16_t)]; uint8_t arr[]; };
dz_static_assert(sizeof(struct dz_query_s) % sizeof(__m128i) == 0);

/* node (reference) */
struct dz_node_s { int32_t id, len; uint8_t const *ptr; };
dz_static_assert(sizeof(struct dz_node_s) % sizeof(__m128i) == 0);

/* DP matrix structures */
struct dz_swgv_s { __m128i e, f, s; };						/* followed by dz_cap_s */
struct dz_range_s { uint32_t spos, epos; };					/* placed just after every score vector to indicate the length */

struct dz_head_s {
	struct dz_range_s r;
	uint32_t rch, n_forefronts;
};
struct dz_cap_s {											/* followed by dz_forefront_s; spos and epos are shared to forefront_s */
	struct dz_range_s r;
	uint32_t rch; int32_t rrem;
};
struct dz_forefront_s {
	struct dz_range_s r;
	uint32_t rid;
	int32_t rlen;
	uint32_t rsum, rcnt; int32_t max, inc;
	struct dz_query_s const *query;
	struct dz_cap_s const *mcap;
};
dz_static_assert(sizeof(struct dz_swgv_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_cap_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_forefront_s) % sizeof(__m128i) == 0);
#define dz_swgv(_p)					( (struct dz_swgv_s *)(_p) )
#define dz_cswgv(_p)				( (struct dz_swgv_s const *)(_p) )
#define dz_range(_p)				( (struct dz_range_s *)(_p) )
#define dz_crange(_p)				( (struct dz_range_s const *)(_p) )
#define dz_head(_p)					( (struct dz_head_s *)(_p) )
#define dz_chead(_p)				( (struct dz_head_s const *)(_p) )
#define dz_cap(_p)					( (struct dz_cap_s *)(_p) )
#define dz_ccap(_p)					( (struct dz_cap_s const *)(_p) )
#define dz_forefront(_p)			( (struct dz_forefront_s *)(_p) )
#define dz_cff(_p)					( (struct dz_forefront_s const *)(_p) )

/* alignment path */
struct dz_path_span_s {
	uint32_t id;
	uint32_t offset;
};
struct dz_alignment_s {
	struct dz_path_span_s const *span;
	uint8_t const *path;
	uint32_t span_length, path_length, ref_length, query_length;
	int32_t rrem, score;
	uint32_t mismatch_count, match_count, ins_count, del_count;
};

/* context (constants and working buffers) */
struct dz_mem_block_s { struct dz_mem_block_s *next; size_t size; };
struct dz_stack_s { struct dz_mem_block_s *curr; uint8_t *top, *end; uint64_t _pad[3]; };
struct dz_mem_s { struct dz_mem_block_s blk; struct dz_stack_s stack; };
#define dz_mem_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )

struct dz_s {
	int8_t matrix[32];					/* DNA score matrix */
	uint16_t iiv[8], iev[8];
	uint16_t div[8], dev[8];			/* deletion: open and extend */
	uint16_t xt, bonus, max_gap_len, _pad[9];
	struct dz_forefront_s const *root;
	int8_t protein_matrix[];			/* allocated at the tail if not DNA matrix */
};
dz_static_assert(sizeof(struct dz_s) % sizeof(__m128i) == 0);
#define dz_mem(_self)				( (struct dz_mem_s *)(_self) - 1 )

#define dz_root(_self)				( (struct dz_forefront_s const **)(&_self->root) )
#define dz_is_terminated(_ff)		( dz_cff(_ff)->r.spos >= dz_cff(_ff)->r.epos )
#define dz_gt(_ff)					( dz_cff(_ff)->inc > 0 )
#define dz_geq(_ff)					( dz_cff(_ff)->mcap != NULL )

#ifdef DZ_PRINT_VECTOR
#define print_vector(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	dz_rm_ofs(_mm_extract_epi16(v, 7)), \
	dz_rm_ofs(_mm_extract_epi16(v, 6)), \
	dz_rm_ofs(_mm_extract_epi16(v, 5)), \
	dz_rm_ofs(_mm_extract_epi16(v, 4)), \
	dz_rm_ofs(_mm_extract_epi16(v, 3)), \
	dz_rm_ofs(_mm_extract_epi16(v, 2)), \
	dz_rm_ofs(_mm_extract_epi16(v, 1)), \
	dz_rm_ofs(_mm_extract_epi16(v, 0))); \
}
#else
#define print_vector(v) ;
#endif

/**
 * @fn dz_malloc, dz_free
 * @brief aligned and margined malloc and free
 */
static __dz_force_inline
void *dz_malloc(
	size_t size)
{
	void *ptr = NULL;

	/* roundup to align boundary, add margin at the head and forefront */
	size = dz_roundup(size, DZ_MEM_ALIGN_SIZE);
	if(posix_memalign(&ptr, DZ_MEM_ALIGN_SIZE, size + 2 * DZ_MEM_MARGIN_SIZE) != 0) {
		debug("posix_memalign failed");
		dz_trap(); return(NULL);
	}
	debug("posix_memalign(%p), size(%lu)", ptr, size);
	return((void *)((uint8_t *)ptr + DZ_MEM_MARGIN_SIZE));
}
static __dz_force_inline
void dz_free(
	void *ptr)
{
	debug("free(%p)", (uint8_t const *)ptr - DZ_MEM_MARGIN_SIZE);
	free((void *)((uint8_t *)ptr - DZ_MEM_MARGIN_SIZE));
	return;
}

unittest() {
	size_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(size_t i = 0; i < sizeof(size) / sizeof(size_t); i++) {
		uint8_t *p = (uint8_t *)dz_malloc(size[i]);
		ut_assert(p != NULL);
		volatile uint8_t a = p[-32], b = p[0], c = p[size[i]], d = p[size[i] + 32];
		dz_unused(a); dz_unused(b); dz_unused(c); dz_unused(d);
		memset(p, 0, size[i]);
		dz_free(p);
	}
}

/**
 * @fn dz_mem_init, dz_mem_destroy, dz_mem_add_stack, dz_mem_malloc, dz_mem_flush
 * @brief stack chain
 */
static __dz_force_inline
void dz_mem_flush(
	struct dz_mem_s *mem)
{
	mem->stack.curr = &mem->blk;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem + 1), DZ_MEM_ALIGN_SIZE);
	mem->stack.end = (uint8_t *)mem + mem->blk.size - DZ_MEM_MARGIN_SIZE;
	return;
}
static __dz_force_inline
struct dz_mem_s *dz_mem_init(
	size_t size)
{
	size = dz_min2(DZ_MEM_INIT_SIZE, sizeof(struct dz_mem_s) + dz_roundup(size, DZ_MEM_ALIGN_SIZE));
	struct dz_mem_s *mem = (struct dz_mem_s *)dz_malloc(size);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init mem object then stack pointers */
	mem->blk = (struct dz_mem_block_s){ .next = NULL, .size = size };
	dz_mem_flush(mem);
	return(mem);
}
static __dz_force_inline
void dz_mem_destroy(
	struct dz_mem_s *mem)
{
	struct dz_mem_block_s *blk = mem->blk.next;
	debug("cleanup memory chain, blk(%p)", blk);
	while(blk != NULL) {
		struct dz_mem_block_s *next = blk->next;
		debug("free blk(%p), next(%p)", blk, next);
		dz_free(blk); blk = next;
	}
	dz_free(mem);
	return;
}
static __dz_force_inline
uint64_t dz_mem_add_stack(
	struct dz_mem_s *mem,
	size_t size)
{
	debug("add_stack, ptr(%p)", mem->stack.curr->next);
	if(mem->stack.curr->next == NULL) {
		/* current stack is the forefront of the memory block chain, add new block */
		size = dz_max2(
			size + dz_roundup(sizeof(struct dz_mem_block_s), DZ_MEM_ALIGN_SIZE),
			2 * mem->stack.curr->size
		);
		struct dz_mem_block_s *blk = (struct dz_mem_block_s *)dz_malloc(size);
		debug("malloc called, blk(%p)", blk);
		if(blk == NULL) { return(1); }

		/* link new node to the forefront of the current chain */
		mem->stack.curr->next = blk;

		/* initialize the new memory block */
		blk->next = NULL;
		blk->size = size;
	}
	/* follow a forward link and init stack pointers */
	mem->stack.curr = mem->stack.curr->next;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem->stack.curr + 1), DZ_MEM_ALIGN_SIZE);
	mem->stack.end = (uint8_t *)mem->stack.curr + mem->stack.curr->size - DZ_MEM_MARGIN_SIZE;
	return(0);
}
static __dz_force_inline
void *dz_mem_malloc(
	struct dz_mem_s *mem,
	size_t size)
{
	if(dz_mem_stack_rem(mem) < 4096) { dz_mem_add_stack(mem, 0); }
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += dz_roundup(size, sizeof(__m128i));
	return(ptr);
}

unittest() {
	size_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(size_t i = 0; i < sizeof(size) / sizeof(size_t); i++) {
		struct dz_mem_s *mem = dz_mem_init(size[i]);
		uint8_t *p = (uint8_t *)dz_mem_malloc(mem, size[i]);
		ut_assert(p != NULL);
		memset(p, 0, size[i]);
		dz_mem_destroy(mem);
	}
}

/**
 * vector update macros
 */
#define _calc_next_size(_sp, _ep, _nt) ({ \
	size_t forefront_arr_size = dz_roundup(sizeof(struct dz_forefront_s *) * (_nt), sizeof(__m128i)); \
	size_t est_column_size = 2 * ((_ep) - (_sp)) * sizeof(struct dz_swgv_s); \
	size_t next_req = forefront_arr_size + est_column_size + sizeof(struct dz_cap_s); \
	/* debug("est_column_size(%lu), next_req(%lu)", est_column_size, next_req); */ \
	next_req; \
})
#define _init_cap(_adj, _rch, _forefronts, _n_forefronts) ({ \
	/* push forefront pointers */ \
	size_t forefront_arr_size = dz_roundup(sizeof(struct dz_forefront_s *) * (_n_forefronts), sizeof(__m128i)); \
	struct dz_forefront_s const **dst = (struct dz_forefront_s const **)(dz_mem(self)->stack.top + forefront_arr_size); \
	struct dz_forefront_s const **src = (struct dz_forefront_s const **)(_forefronts); \
	for(size_t i = 0; i < (_n_forefronts); i++) { dst[-((int64_t)(_n_forefronts)) + i] = src[i]; } \
	/* push head-cap info */ \
	struct dz_head_s *_head = (struct dz_head_s *)dst; \
	(_head)->r.spos = (_adj);					/* save merging adjustment */ \
	(_head)->r.epos = 0;						/* head marked as zero */ \
	(_head)->rch = (_rch);						/* rch for the first column */ \
	(_head)->n_forefronts = (_n_forefronts);	/* record n_forefronts */ \
	debug("create head cap(%p), n_forefronts(%lu)", _head, (uint64_t)(_n_forefronts)); \
	dz_cap(_head); \
})
#define _begin_column_head(_spos, _epos, _adj, _forefronts, _n_forefronts) ({ \
	/* calculate sizes */ \
	size_t next_req = _calc_next_size(_spos, _epos, _n_forefronts); \
	/* allocate from heap */ \
	if(dz_mem_stack_rem(dz_mem(self)) < next_req) { dz_mem_add_stack(dz_mem(self), 0); } \
	/* push head-cap */ \
	struct dz_cap_s *cap = _init_cap(_adj, 0xff, _forefronts, _n_forefronts); \
	/* return array pointer */ \
	(struct dz_swgv_s *)(cap + 1) - (_spos); \
})
#define _begin_column(_w, _rch, _rlen) ({ \
	/* push cap info */ \
	struct dz_cap_s *cap = dz_cap(dz_mem(self)->stack.top); \
	/* calculate sizes */ \
	size_t next_req = _calc_next_size((_w).r.spos, (_w).r.epos, 0); \
	/* allocate from heap */ \
	cap->rch = (_rch);							/* record rch for the next column */ \
	cap->rrem = (_rlen);						/* record rlen for use in traceback */ \
	if(dz_likely(dz_mem_stack_rem(dz_mem(self)) < next_req)) { \
		dz_mem_add_stack(dz_mem(self), 0); \
		cap = _init_cap(0, _rch, &cap, 1); \
	} \
	debug("create column(%p), [%u, %u), span(%u), rrem(%ld), max(%d), inc(%d)", cap, (_w).r.spos, (_w).r.epos, (_w).r.epos - (_w).r.spos, (_rlen), (_w).max, (_w).inc); \
	/* return array pointer */ \
	(struct dz_swgv_s *)(cap + 1) - (_w).r.spos; \
})
#define _end_column(_p, _spos, _epos) ({ \
	/* write back the stack pointer and return a cap */ \
	struct dz_range_s *r = dz_range(&dz_swgv(_p)[(_epos)]); \
	debug("create range(%p), [%u, %u)", r, (_spos), (_epos)); \
	dz_mem(self)->stack.top = (uint8_t *)r; \
	r->spos = (_spos); r->epos = (_epos); \
	(struct dz_cap_s *)r; \
})
#define _end_matrix(_p, _wp, _rrem) ({ \
	/* create forefront object */ \
	struct dz_forefront_s *forefront = dz_forefront(&(dz_swgv(_p))[(_wp)->r.epos]); \
	/* if((_wp)->r.spos >= (_wp)->r.epos) { (_wp)->r.epos = 0; } */ \
	forefront->rid = (_wp)->rid; \
	forefront->rlen = (_wp)->rlen; \
	forefront->rsum = (_wp)->rsum + ((_wp)->rlen < 0 ? -1 : 1) * ((_wp)->rlen - (_rrem)); \
	forefront->rcnt = (_wp)->rcnt + 1; \
	forefront->max = (_wp)->max + (_wp)->inc; \
	forefront->inc = (_wp)->inc; \
	forefront->query = (_wp)->query; \
	forefront->mcap = (_wp)->mcap; \
	debug("create forefront(%p), [%u, %u), max(%d), inc(%d), rlen(%d), rrem(%d)", forefront, (_wp)->r.spos, (_wp)->r.epos, (_wp)->max, (_wp)->inc, (int32_t)(_wp)->rlen, (int32_t)(_rrem)); \
	/* write back stack pointer */ \
	dz_mem(self)->stack.top = (uint8_t *)(forefront + 1); \
	/* return the forefront pointer */ \
	(struct dz_forefront_s const *)forefront; \
})
/**
 * @macro _calc_bonus
 * @brief add bonus to cells on the tail row; bonus constant vectors are always placed just before parr
 */
#ifdef DZ_FULL_LENGTH_BONUS
#define _init_bonus(_query) \
	uint32_t blim = (_query)->blen - 1; \
	uint8_t const *pbonus = (_query)->arr;

#define _add_bonus(_i, _v)			( _mm_add_epi16((_v), _mm_load_si128(&((__m128i const *)pbonus)[-2 + ((_i) == blim)])) )
#else
#define _init_bonus(_query)			;
#define _add_bonus(_i, _v)			( (_v) )
#endif

#if defined(DZ_NUCL_ASCII)
#define _init_rch(_query, _rt, _rrem) \
	_init_bonus(_query); \
	uint32_t rch = conv[_rt[-_rrem] & 0x0f]; \
	/* debug("rch(%c, %u, %x)", _rt[-_rrem], rch, rch); */ \
	uint8_t const *parr = (_query)->arr; \
	__m128i const rv = _mm_set1_epi8(rch);

#define _calc_score_profile(_i) ({ \
	__m128i qv = _mm_loadl_epi64((__m128i const *)&parr[(_i) * L]); \
	__m128i sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)self->matrix), _mm_or_si128(rv, qv))); \
	/* print_vector(_mm_cvtepi8_epi16(rv)); print_vector(_mm_cvtepi8_epi16(qv)); */ \
	sc; \
})
#elif defined(DZ_NUCL_2BIT)
#define _init_rch(_query, _rt, _rrem) \
	_init_bonus(_query); \
	uint32_t rch = dir < 0 ? _rt[-_rrem] : (_rt[-_rrem] ^ 0x03); \
	/* debug("rch(%c, %u, %x)", _rt[-_rrem], rch, rch); */ \
	uint8_t const *parr = (_query)->arr; \
	__m128i const rv = _mm_set1_epi8(rch);

#define _calc_score_profile(_i) ({ \
	__m128i qv = _mm_loadl_epi64((__m128i const *)&parr[(_i) * L]); \
	__m128i sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)self->matrix), _mm_or_si128(rv, qv))); \
	sc; \
})
#else /* DZ_PROTEIN */
#define _init_rch(_query, _rt, _rrem) \
	_init_bonus(_query); \
	uint32_t rch = _rt[-_rrem] & 0x1f; \
	/* debug("rch(%c, %u, %x)", _rt[-_rrem], rch, rch); */ \
	int8_t const *parr = (int8_t const *)&(_query)->arr[rch * (_query)->blen * L];

#define _calc_score_profile(_i) ({ \
	__m128i sc = _mm_cvtepi8_epi16(_mm_loadl_epi64((__m128i const *)&parr[(_i) * L])); \
	sc; \
})
#endif

#define _load_vector(_p) \
	__m128i e = _mm_load_si128((__m128i const *)(&dz_cswgv(_p)->e)); /* print_vector(e); */ \
	__m128i s = _mm_load_si128((__m128i const *)(&dz_cswgv(_p)->s)); /* print_vector(s); */
#define _update_vector(_p) { \
	__m128i sc = _calc_score_profile(_p); \
	__m128i te = _mm_subs_epu16(_mm_max_epu16(e, _mm_subs_epu16(s, div)), dev); \
	/* print_vector(_mm_alignr_epi8(s, ps, 14)); print_vector(sc); */ \
	__m128i ts = _mm_max_epu16(te, _mm_subs_epu16(_mm_adds_epu16(sc, _mm_alignr_epi8(s, ps, 14)), ofs)); ps = s; \
	__m128i tf = _mm_max_epu16(_mm_subs_epu16(ts, iiv), _mm_subs_epu16(_mm_srli_si128(f, 14), iev1)); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 2), iev1)); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 4), iev2)); \
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 8), iev4)); \
	ts = _mm_max_epu16(ts, tf); \
	maxv = _mm_max_epu16(maxv, _add_bonus(_p, ts)); \
	/* print_vector(te); print_vector(_add_bonus(_p, ts)); print_vector(tf); print_vector(maxv);*/ \
	e = te; f = tf; s = ts; \
}
#define _store_vector(_p) { \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->e), e); \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->f), f); \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->s), s); \
}
#define _hmax_vector(_v) ({ \
	__m128i _t = _mm_max_epu16(_v, _mm_srli_si128(_v, 8)); \
	_t = _mm_max_epu16(_t, _mm_srli_si128(_t, 4)); \
	_t = _mm_max_epu16(_t, _mm_srli_si128(_t, 2)); \
	dz_rm_ofs(_mm_extract_epi16(_t, 0)); \
})
#define _test_xdrop(_s, _xtv) ({ \
	__m128i xtest = _mm_subs_epu16(_s, _xtv); \
	/* print_vector(_s); print_vector(_xtv); */ \
	dz_is_all_zero(xtest); \
})

/**
 * @fn dz_init, dz_destroy
 */
static __dz_vectorize
struct dz_s *dz_init_intl(
	int8_t const *score_matrix,		/* match award in positive, mismatch penalty in negative. s(A,A) at [0], s(A,C) at [1], ... s(T,T) at [15] where s(ref_base, query_base) is a score function */
	uint16_t ins_open,				/* gap penalties in positive */
	uint16_t ins_extend,
	uint16_t del_open,				/* gap penalties in positive */
	uint16_t del_extend,
	uint64_t max_gap_len,			/* as X-drop threshold */
	uint16_t full_length_bonus)		/* end-to-end mapping bonus; only activated when compiled with -DDZ_FULL_LENGTH_BONUS */
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	uint16_t const gi = dz_max2(ins_open, del_open), ge = dz_max2(ins_extend, del_extend);
	/*
	static uint8_t const transpose[16] __attribute__(( aligned(16) )) = {
		0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15
	};
	*/

	struct dz_mem_s *mem = dz_mem_init(DZ_MEM_INIT_SIZE);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}
	#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
		struct dz_s *self = (struct dz_s *)dz_mem_malloc(mem, sizeof(struct dz_s));

		/* constants */
		__m128i const tmat = _mm_loadu_si128((__m128i const *)score_matrix);
		__m128i const ofs = _mm_set1_epi8(DZ_SCORE_OFS);
		_mm_store_si128((__m128i *)&self->matrix[0], _mm_add_epi8(tmat, ofs));
		_mm_store_si128((__m128i *)&self->matrix[16], _mm_setzero_si128());
	#else
		struct dz_s *self = (struct dz_s *)dz_mem_malloc(mem, sizeof(struct dz_s) + DZ_MAT_SIZE * DZ_MAT_SIZE + 2 * sizeof(__m128i));

		/* clear the first matrix field for protein */
		_mm_store_si128((__m128i *)&self->matrix[0], _mm_setzero_si128());
		_mm_store_si128((__m128i *)&self->matrix[16], _mm_setzero_si128());

		/* transpose */
		for(uint64_t i = 0; i < DZ_MAT_SIZE; i++) {
			for(uint64_t j = 0; j < DZ_MAT_SIZE; j++) {
				self->protein_matrix[j * DZ_MAT_SIZE + i] = score_matrix[i * DZ_MAT_SIZE + j] + DZ_SCORE_OFS;
			}
		}
	#endif

	/* insertion penalties */
	__m128i const iiv = _mm_set1_epi16(ins_open);
	__m128i const iev = _mm_set1_epi16(ins_extend);
	_mm_store_si128((__m128i *)self->iiv, iiv);
	_mm_store_si128((__m128i *)self->iev, iev);

	/* deletion penalties */
	__m128i const div = _mm_set1_epi16(del_open);
	__m128i const dev = _mm_set1_epi16(del_extend);
	_mm_store_si128((__m128i *)self->div, div);
	_mm_store_si128((__m128i *)self->dev, dev);

	self->xt = gi + ge * max_gap_len;			/* X-drop threshold */
	self->bonus = full_length_bonus;
	self->max_gap_len = max_gap_len;			/* save raw value */
	debug("gi(%u), ge(%u), xdrop_threshold(%u), full_length_bonus(%u), max_gap_len(%lu)", gi, ge, self->xt, self->bonus, self->max_gap_len);

	/* create root head */
	struct dz_cap_s *cap = (struct dz_cap_s *)dz_mem_malloc(mem, sizeof(struct dz_cap_s));
	_mm_store_si128((__m128i *)cap, _mm_setzero_si128());

	/* calc vector length; query = NULL for the first (root) column */
	max_gap_len = dz_roundup(max_gap_len, L);
	struct dz_forefront_s w = { { 0, (uint32_t)(max_gap_len / L) }, 0, 0, 0, 0, 0 /* dz_add_ofs(0) */, 0, NULL, NULL };
	struct dz_forefront_s *a = &w;

	/* malloc the first column */
	struct dz_swgv_s *dp = _begin_column_head(0, max_gap_len / L, 0, &a, 0);

	/* fill the root (the leftmost) column; first init vectors */
	__m128i s = _mm_setr_epi16(
		dz_add_ofs(0),
		dz_add_ofs(-(ins_open +     ins_extend)),
		dz_add_ofs(-(ins_open + 2 * ins_extend)),
		dz_add_ofs(-(ins_open + 3 * ins_extend)),
		dz_add_ofs(-(ins_open + 4 * ins_extend)),
		dz_add_ofs(-(ins_open + 5 * ins_extend)),
		dz_add_ofs(-(ins_open + 6 * ins_extend)),
		dz_add_ofs(-(ins_open + 7 * ins_extend))
	);
	__m128i const e = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));
	__m128i const xtv = _mm_set1_epi16(dz_add_ofs(-self->xt));

	/* until the X-drop test fails on all the cells in a vector */
	for(size_t p = 0; p < max_gap_len / L; p++) {
		__m128i const f = s;
		if(dz_unlikely(_test_xdrop(s, xtv))) { debug("p(%lu)", p); w.r.epos = p; break; }
		_store_vector(&dp[p]); print_vector(s);
		if(p == 0) {
			/* reset s vector */
			s = _mm_setr_epi16(
				dz_add_ofs(-ins_open),
				dz_add_ofs(-(ins_open +     ins_extend)),
				dz_add_ofs(-(ins_open + 2 * ins_extend)),
				dz_add_ofs(-(ins_open + 3 * ins_extend)),
				dz_add_ofs(-(ins_open + 4 * ins_extend)),
				dz_add_ofs(-(ins_open + 5 * ins_extend)),
				dz_add_ofs(-(ins_open + 6 * ins_extend)),
				dz_add_ofs(-(ins_open + 7 * ins_extend))
			);
		}
		s = _mm_subs_epu16(s, _mm_slli_epi16(iev, 3));		/* iev[i] * 8 */
	}

	/* done; create and return a forefront object */
	_end_column(dp, w.r.spos, w.r.epos);
	self->root = _end_matrix(dp, &w, 0);
	debug("self(%p), mem(%p), root(%p)", self, dz_mem(self), self->root);
	return(self);
}

#ifdef DZ_FULL_LENGTH_BONUS
static __dz_vectorize
struct dz_s *dz_init(
	int8_t const *score_matrix,
	uint16_t gap_open,
	uint16_t gap_extend,
	uint64_t max_gap_len,
	uint16_t full_length_bonus)
{
	return(dz_init_intl(score_matrix, gap_open, gap_extend, gap_open, gap_extend, max_gap_len, full_length_bonus));
}
static __dz_vectorize
struct dz_s *dz_initx(
	int8_t const *score_matrix,
	uint16_t ins_open,
	uint16_t ins_extend,
	uint16_t del_open,
	uint16_t del_extend,
	uint64_t max_gap_len,
	uint16_t full_length_bonus)
{
	return(dz_init_intl(score_matrix, ins_open, ins_extend, del_open, del_extend, max_gap_len, full_length_bonus));
}

#else
static __dz_vectorize
struct dz_s *dz_init(
	int8_t const *score_matrix,
	uint16_t gap_open,
	uint16_t gap_extend,
	uint64_t max_gap_len)
{
	return(dz_init_intl(score_matrix, gap_open, gap_extend, gap_open, gap_extend, max_gap_len, 0));
}
static __dz_vectorize
struct dz_s *dz_initx(
	int8_t const *score_matrix,
	uint16_t ins_open,
	uint16_t ins_extend,
	uint16_t del_open,
	uint16_t del_extend,
	uint64_t max_gap_len)
{
	return(dz_init_intl(score_matrix, ins_open, ins_extend, del_open, del_extend, max_gap_len, 0));
}

#endif

static __dz_vectorize
void dz_destroy(
	struct dz_s *self)
{
	debug("self(%p)", self);
	if(self == NULL) { return; }
	dz_mem_destroy(dz_mem(self));
	return;
}
static __dz_vectorize
void dz_flush(
	struct dz_s *self)
{
	dz_mem_flush(dz_mem(self));
	dz_mem(self)->stack.top = (uint8_t *)(self->root + 1);
	return;
}

#if defined(DZ_NUCL_ASCII)
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
	"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
	"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
	"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
	"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
	"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
	"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
	"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

static int8_t const dz_unittest_score_matrix[16] = {
	2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2
};

#elif defined(DZ_NUCL_2BIT)
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"\x0\x2\x1\x3\x3\x3\x3\x1\x0\x3\x3\x1\x3\x2\x0\x1\x3\x2\x1\x0\x0\x1\x2\x2\x2\x1\x0\x0\x3\x0\x3\x2\x3\x1\x3"
	"\x1\x3\x2\x3\x2\x3\x2\x2\x0\x3\x3\x0\x0\x0\x0\x0\x0\x0\x2\x0\x2\x3\x2\x3\x1\x3\x2\x0\x3\x0\x2\x1\x0\x2\x1"
	"\x3\x3\x1\x3\x2\x0\x0\x1\x3\x2\x2\x3\x3\x0\x1\x1\x3\x2\x1\x1\x2\x3\x2\x0\x2\x3\x0\x0\x0\x3\x3\x0\x0\x0\x0"
	"\x3\x3\x3\x3\x0\x3\x3\x2\x0\x1\x3\x3\x0\x2\x2\x3\x1\x0\x1\x3\x0\x0\x0\x3\x0\x1\x3\x3\x3\x0\x0\x1\x1\x0\x0"
	"\x3\x0\x3\x0\x2\x2\x1\x0\x3\x0\x2\x1\x2\x1\x0\x1\x0\x2\x0\x1\x0\x2\x0\x3\x0\x0\x0\x0\x0\x3\x3\x0\x1\x0\x2"
	"\x0\x2\x3\x0\x1\x0\x1\x0\x0\x1\x0\x3\x1\x1\x0\x3\x2\x0\x0\x0\x1\x2\x1\x0\x3\x3\x0\x2\x1\x0\x1\x1\x0\x1\x1"
	"\x0\x3\x3\x0\x1\x1\x0\x1\x1\x0\x1\x1\x0\x3\x1\x0\x1\x1\x0\x3\x3\x0\x1\x1\x0\x1\x0\x2\x2\x3\x0\x0\x1\x2\x2"
	"\x3\x2\x1\x2\x2\x2\x1\x3\x2\x0\x1\x2\x1\x2\x3\x0\x1\x0\x2\x2\x0\x0\x0\x1\x0\x1\x0\x2\x0\x0\x0\x0\x0\x0\x2"
	"\x1\x1\x1\x2\x1\x0\x1\x1\x3\x2\x0\x1\x0\x2\x3\x2\x1\x2\x2\x2\x1\x3\x3\x3\x3\x3\x3\x3\x3\x3\x1\x2\x0\x1\x1"
	"\x0\x0\x0\x2\x2\x3\x0\x0\x1\x2\x0\x2\x2\x3\x0\x0\x1\x0\x0\x1\x1\x0\x3\x2\x1\x2\x0\x2\x3\x2\x3\x3\x2\x0\x0"
	"\x2\x3\x3\x1\x2\x2\x1\x2\x2\x3\x0\x1\x0\x3\x1\x0\x2\x3\x2\x2\x1\x0\x0\x0\x3\x2\x1\x0\x2\x0\x0\x1\x2\x3\x3"
	"\x3\x3\x1\x3\x2\x1\x2\x3\x2\x3\x3\x2\x1\x1\x2\x0\x3\x0\x3\x3\x1\x3\x2\x2\x0\x0\x0\x2\x1\x0\x0\x3\x2\x1\x1"
	"\x0\x2\x2\x1\x0\x2\x2\x2\x2\x1\x0\x2\x2\x3\x2\x2\x1\x1\x0\x1\x1\x2\x3\x1\x1\x3\x1\x3\x1\x3\x2\x1\x1\x1\x1"
	"\x1\x2\x1\x1\x0\x0\x0\x0\x3\x1\x0\x1\x1\x0\x0\x1\x1\x0\x1\x1\x3\x2\x2\x3\x2\x2\x1\x2\x0\x3\x2\x0\x3\x3\x2"
	"\x0\x0\x0\x0\x0\x0\x1\x1\x0\x3\x3\x0\x2\x1\x2\x2\x1\x1\x0\x2\x2\x0\x3\x2\x1\x3\x3\x3\x0\x1\x1\x1\x0\x0\x3"
	"\x0\x3\x1\x0\x2\x1\x2\x0\x3\x2\x1\x1\x2\x0\x0\x1\x2\x3\x0\x3\x3\x3\x3\x3\x2\x1\x1\x2\x0\x0\x1\x3\x3\x3\x3";

static int8_t const dz_unittest_score_matrix[16] = {
	2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2
};

#elif defined(DZ_PROTEIN)
static size_t const dz_unittest_query_length = 511;
static char const *dz_unittest_query =
	"MATLVQTGKAKQLTLLGFFAITASMVMAVYEYPTFATSGFSLVFFLLLGGILWFIPVGLC"
	"AAEMATVDGWEEGGVFAWVSNTLGPRWGFAAISFGYLQIAIGFIPMLYFVLGALSYILKW"
	"PALNEDPITKTIAALIILWALALTQFGGTKYTARIAKVGFFAGILLPAFILIALAAIYLH"
	"SGAPVAIEMDSKTFFPDFSKVGTLVVFVAFILSYMGVEASATHVNEMSNPGRDYPLAMLL"
	"LMVAAICLSSVGGLSIAMVIPGNEINLSAGVMQTFTVLMSHVAPEIEWTVRVISALLLLG"
	"VLAEIASWIVGPSRGMYVTAQKNLLPAAFAKMNKNGVPVTLVISQLVITSIALIILTNTG"
	"GGNNMSFLIALALTVVIYLCAYFMLFIGYIVLVLKHPDLKRTFNIPGGKGVKLVVAIVGL"
	"LTSIMAFIVSFLPPDNIQGDSTDMYVELLVVSFLVVLALPFILYAVHDRKGKANTGVTLE"
	"PINSQNAPKGHFFLHPRARSPHYIVMNDKKH";

/* BLOSUM62 */
dz_static_assert(DZ_MAT_SIZE == 32);
static int8_t const dz_unittest_score_matrix[DZ_MAT_SIZE * DZ_MAT_SIZE] = {
	 1, -4, -4, -4, -4, -4, -4, -4, -4, -4,  0, -4, -4, -4, -4,  0, -4, -4, -4, -4, -4,  0, -4, -4, -4, -4, -4,  0,  0,  0,  0,  0,
	-4,  4, -2,  0, -2, -1, -2,  0, -2, -1,  0, -1, -1, -1, -2,  0, -1, -1, -1,  1,  0,  0,  0, -3,  0, -2, -1,  0,  0,  0,  0,  0,
	-4, -2,  4, -3,  4,  1, -3, -1,  0, -3,  0,  0, -4, -3,  3,  0, -2,  0, -1,  0, -1,  0, -3, -4, -1, -3,  1,  0,  0,  0,  0,  0,
	-4,  0, -3,  9, -3, -4, -2, -3, -3, -1,  0, -3, -1, -1, -3,  0, -3, -3, -3, -1, -1,  0, -1, -2, -2, -2, -3,  0,  0,  0,  0,  0,
	-4, -2,  4, -3,  6,  2, -3, -1, -1, -3,  0, -1, -4, -3,  1,  0, -1,  0, -2,  0, -1,  0, -3, -4, -1, -3,  1,  0,  0,  0,  0,  0,
	-4, -1,  1, -4,  2,  5, -3, -2,  0, -3,  0,  1, -3, -2,  0,  0, -1,  2,  0,  0, -1,  0, -2, -3, -1, -2,  4,  0,  0,  0,  0,  0,
	-4, -2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3,  0, -4, -3, -3, -2, -2,  0, -1,  1, -1,  3, -3,  0,  0,  0,  0,  0,
	-4,  0, -1, -3, -1, -2, -3,  6, -2, -4,  0, -2, -4, -3,  0,  0, -2, -2, -2,  0, -2,  0, -3, -2, -1, -3, -2,  0,  0,  0,  0,  0,
	-4, -2,  0, -3, -1,  0, -1, -2,  8, -3,  0, -1, -3, -2,  1,  0, -2,  0,  0, -1, -2,  0, -3, -2, -1,  2,  0,  0,  0,  0,  0,  0,
	-4, -1, -3, -1, -3, -3,  0, -4, -3,  4,  0, -3,  2,  1, -3,  0, -3, -3, -3, -2, -1,  0,  3, -3, -1, -1, -3,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4, -1,  0, -3, -1,  1, -3, -2, -1, -3,  0,  5, -2, -1,  0,  0, -1,  1,  2,  0, -1,  0, -2, -3, -1, -2,  1,  0,  0,  0,  0,  0,
	-4, -1, -4, -1, -4, -3,  0, -4, -3,  2,  0, -2,  4,  2, -3,  0, -3, -2, -2, -2, -1,  0,  1, -2, -1, -1, -3,  0,  0,  0,  0,  0,
	-4, -1, -3, -1, -3, -2,  0, -3, -2,  1,  0, -1,  2,  5, -2,  0, -2,  0, -1, -1, -1,  0,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,
	-4, -2,  3, -3,  1,  0, -3,  0,  1, -3,  0,  0, -3, -2,  6,  0, -2,  0,  0,  1,  0,  0, -3, -4, -1, -2,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4, -1, -2, -3, -1, -1, -4, -2, -2, -3,  0, -1, -3, -2, -2,  0,  7, -1, -2, -1, -1,  0, -2, -4, -2, -3, -1,  0,  0,  0,  0,  0,
	-4, -1,  0, -3,  0,  2, -3, -2,  0, -3,  0,  1, -2,  0,  0,  0, -1,  5,  1,  0, -1,  0, -2, -2, -1, -1,  3,  0,  0,  0,  0,  0,
	-4, -1, -1, -3, -2,  0, -3, -2,  0, -3,  0,  2, -2, -1,  0,  0, -2,  1,  5, -1, -1,  0, -3, -3, -1, -2,  0,  0,  0,  0,  0,  0,
	-4,  1,  0, -1,  0,  0, -2,  0, -1, -2,  0,  0, -2, -1,  1,  0, -1,  0, -1,  4,  1,  0, -2, -3,  0, -2,  0,  0,  0,  0,  0,  0,
	-4,  0, -1, -1, -1, -1, -2, -2, -2, -1,  0, -1, -1, -1,  0,  0, -1, -1, -1,  1,  5,  0,  0, -2,  0, -2, -1,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-4,  0, -3, -1, -3, -2, -1, -3, -3,  3,  0, -2,  1,  1, -3,  0, -2, -2, -3, -2,  0,  0,  4, -3, -1, -1, -2,  0,  0,  0,  0,  0,
	-4, -3, -4, -2, -4, -3,  1, -2, -2, -3,  0, -3, -2, -1, -4,  0, -4, -2, -3, -3, -2,  0, -3, 11, -2,  2, -3,  0,  0,  0,  0,  0,
	-4,  0, -1, -2, -1, -1, -1, -1, -1, -1,  0, -1, -1, -1, -1,  0, -2, -1, -1,  0,  0,  0, -1, -2, -1, -1, -1,  0,  0,  0,  0,  0,
	-4, -2, -3, -2, -3, -2,  3, -3,  2, -1,  0, -2, -1, -1, -2,  0, -3, -1, -2, -2, -2,  0, -1,  2, -1,  7, -2,  0,  0,  0,  0,  0,
	-4, -1,  1, -3,  1,  4, -3, -2,  0, -3,  0,  1, -3, -1,  0,  0, -1,  3,  0,  0, -1,  0, -2, -3, -1, -2,  4,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

#endif

#ifdef DZ_FULL_LENGTH_BONUS
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1, 20, 0
#else
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1, 20
#endif

unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);
	ut_assert(dz->root != NULL);
	dz_destroy(dz);
}

/**
 * @fn dz_pack_query, dz_pack_query_reverse
 * packed sequence object is allocated inside the context so freed by dz_flush
 */
#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)

static __dz_vectorize
struct dz_query_s *dz_pack_query_forward(
	struct dz_s *self,
	char const *query, size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
	/*
	 * tentative support of full-length bonus; precaclulated bonus score is saved at q->bonus[0] and q->bonus[1]; [0] for non-tail vector and [1] for tail vector
	 */
	q->bonus[L + (qlen % L)] = self->bonus;

	/*
	 * ASCII to 2-bit conversion table (shifted by two bits to be used as the upper (row) shuffle index)
	 * @ABC_DEFG_HIJK_LMNO
	 * PQRS_TUVW_XYZ
	 */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		#ifdef DZ_NUCL_ASCII
			/* ['A' & 0x0f] = A, ['C' & 0x0f] = C, ['G' & 0x0f] = G, ['T' & 0x0f] = T, ['U' & 0x0f] = T, ['N' & 0x0f] = N */
			0, qA, 0, qC, qT, qU, 0, qG, 0, 0, 0, 0, 0, 0, qN, 0
		#else /* DZ_NUCL_2BIT */
			qA, qC, qG, qT, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
		#endif
	};
	__m128i pv = _mm_set1_epi8((int8_t)qN);
	__m128i const fv = _mm_set1_epi8(0x0f);														/* conversion mask */
	__m128i const cv = _mm_load_si128((__m128i const *)conv);									/* conversion table */

	/* until the end of the query sequence */
	for(size_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
		__m128i tv = _mm_shuffle_epi8(cv, _mm_and_si128(qv, fv));
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
	}

	/* continue the same conversion on the remainings */
	// _mm_store_si128((__m128i *)&q->arr[dz_rounddown(qlen, sizeof(__m128i))], _mm_srli_si128(pv, 15));
	q->arr[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
	for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[(uint8_t)query[i] & 0x0f];
	}
	for(size_t i = qlen; i < dz_roundup(qlen + 1, sizeof(__m128i)); i++) {
		q->arr[i + 1] = qS;
	}

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(q);
}
static __dz_vectorize
struct dz_query_s *dz_pack_query_reverse(
	struct dz_s *self,
	char const *query, size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
	q->bonus[L + (qlen % L)] = self->bonus;

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		#ifdef DZ_NUCL_ASCII
			0, qT, 0, qG, qA, qA, 0, qC, 0, 0, 0, 0, 0, 0, qN, 0
		#else /* DZ_NUCL_2BIT */
			qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
		#endif
	};
	static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
		15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	};
	__m128i pv = _mm_set1_epi8((int8_t)qN);
	__m128i const fv = _mm_set1_epi8(0x0f);
	__m128i const cv = _mm_load_si128((__m128i const *)conv);
	__m128i const rv = _mm_load_si128((__m128i const *)rev);

	/* until the end of the query sequence */
	for(size_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
		__m128i tv = _mm_shuffle_epi8(_mm_shuffle_epi8(cv, _mm_and_si128(qv, fv)), rv);
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
	}

	/* continue the same conversion on the remainings */
	q->arr[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
	for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[(uint8_t)query[qlen - 1 - i] & 0x0f];
	}
	for(size_t i = qlen; i < dz_roundup(qlen + 1, sizeof(__m128i)); i++) {
		q->arr[i + 1] = qS;
	}

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(q);
}

#else /* DZ_PROTEIN */
static __dz_vectorize
struct dz_query_s *dz_pack_query_forward(
	struct dz_s *self,
	char const *query, size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + DZ_MAT_SIZE * dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
	q->bonus[L + (qlen % L)] = self->bonus;

	for(uint64_t j = 0; j < DZ_MAT_SIZE; j++) {
		size_t const clen = dz_roundup(qlen + 1, L);
		int8_t *a = (int8_t *)&q->arr[j * clen];
		a[0] = 0;

		int8_t const *conv = &self->protein_matrix[j * DZ_MAT_SIZE];
		__m128i pv = _mm_setzero_si128();
		__m128i const fv = _mm_set1_epi8(0x0f);													/* conversion mask */
		__m128i const lcv = _mm_loadu_si128((__m128i const *)&conv[0]);							/* conversion table */
		__m128i const hcv = _mm_loadu_si128((__m128i const *)&conv[16]);						/* conversion table */

		/* until the end of the query sequence */
		for(size_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
			__m128i const _qv = _mm_and_si128(qv, fv);
			__m128i lv = _mm_shuffle_epi8(lcv, _qv), hv = _mm_shuffle_epi8(hcv, _qv);
			__m128i tv = _mm_blendv_epi8(lv, hv, _mm_slli_epi32(qv, 3));
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
		}

		/* continue the same conversion on the remainings */
		a[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
		for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
			a[i + 1] = conv[(uint8_t)query[i] & 0x1f];
		}
		for(size_t i = qlen; i < dz_roundup(qlen + 1, L); i++) {
			a[i + 1] = 0;				/* as -DZ_SCORE_OFS */
		}
	}
	debug("qlen(%lu), q(%s)", qlen, query);
	return(q);
}
static __dz_vectorize
struct dz_query_s *dz_pack_query_reverse(
	struct dz_s *self,
	char const *query, size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + DZ_MAT_SIZE * dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
	q->bonus[L + (qlen % L)] = self->bonus;

	for(uint64_t j = 0; j < DZ_MAT_SIZE; j++) {
		size_t const clen = dz_roundup(qlen + 1, L);
		int8_t *a = (int8_t *)&q->arr[j * clen];
		a[0] = 0;

		static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
			15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
		};
		int8_t const *conv = &self->protein_matrix[j * DZ_MAT_SIZE];
		__m128i pv = _mm_setzero_si128();
		__m128i const fv = _mm_set1_epi8(0x0f);
		__m128i const lcv = _mm_loadu_si128((__m128i const *)&conv[0]);							/* conversion table */
		__m128i const hcv = _mm_loadu_si128((__m128i const *)&conv[16]);						/* conversion table */
		__m128i const rv = _mm_load_si128((__m128i const *)rev);

		/* toward the head of the query sequence */
		for(size_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
			__m128i const qv = _mm_shuffle_epi8(_mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]), rv);
			__m128i const _qv = _mm_and_si128(qv, fv);
			__m128i lv = _mm_shuffle_epi8(lcv, _qv), hv = _mm_shuffle_epi8(hcv, _qv);
			__m128i tv = _mm_blendv_epi8(lv, hv, _mm_slli_epi32(qv, 3));
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
		}

		/* continue the same conversion on the remainings */
		a[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
		for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
			a[i + 1] = conv[(uint8_t)query[qlen - 1 - i] & 0x0f];
		}
		for(size_t i = qlen; i < dz_roundup(qlen + 1, L); i++) {
			a[i + 1] = 0;				/* as -DZ_SCORE_OFS */
		}
	}

	debug("qlen(%lu), q(%s)", qlen, query);
	return(q);
}
#endif
static __dz_vectorize
struct dz_query_s *dz_pack_query(
	struct dz_s *self,
	char const *query, int64_t qlen)
{
	if(qlen >= 0) {
		return(dz_pack_query_forward(self, query, (size_t)qlen));
	} else {
		return(dz_pack_query_reverse(self, &query[qlen], (size_t)-qlen));
	}
}


unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_unused(q);
	dz_destroy(dz);
}

#define _merge_column(w, adj, forefronts, n_forefronts, query, init_s) ({ \
	for(size_t i = 0; i < n_forefronts; i++) { \
		/* update max and pos */ \
		w.r.spos = dz_min2(w.r.spos, forefronts[i]->r.spos); \
		w.r.epos = dz_max2(w.r.epos, forefronts[i]->r.epos); \
		w.rsum = dz_max2(w.rsum, forefronts[i]->rsum); \
		w.rcnt = dz_max2(w.rcnt, forefronts[i]->rcnt); \
		if(init_s != 0) { w.max = dz_max2(w.max, forefronts[i]->max); } \
		debug("i(%lu, %p), [%u, %u), inc(%d), max(%d)", i, forefronts[i], forefronts[i]->r.spos, forefronts[i]->r.epos, forefronts[i]->inc, forefronts[i]->max); \
	} \
	w.r.spos = w.r.spos > INT32_MAX ? 0 : w.r.spos; \
	w.r.epos = dz_min2(w.r.epos, query->blen);						/* make sure epos and p-iteration index is no larger than blen */ \
	debug("start extension [%u, %u), inc(%d), max(%d), stack(%p, %p), rlen(%d)", w.r.spos, w.r.epos, w.inc, w.max, dz_mem(self)->stack.top, dz_mem(self)->stack.end, (int32_t)rlen); \
	if(init_s != 0) { \
		for(size_t i = 0; i < n_forefronts; i++) { \
			adj[i] = w.max - (forefronts[i]->max - forefronts[i]->inc);	/* base = max - inc */ \
			debug("i(%lu), adj(%lu)", i, adj[i]); \
		} \
	} \
	/* once; merge incoming vectors */ \
	struct dz_swgv_s *cdp = _begin_column_head(w.r.spos, w.r.epos, w.max, forefronts, n_forefronts);	/* allocate memory for the first column */ \
	for(uint64_t p = w.r.spos; p < w.r.epos; p++) { \
		/* memset(cdp, 0xff, sizeof(struct dz_swgv_s) * (w.r.epos - w.r.spos)); */ \
		__m128i const e = _mm_set1_epi16(dz_add_ofs(INT16_MIN)), f = e, s = e; _store_vector(&cdp[p]); \
	} \
	/* paste the last vectors */ \
	for(size_t i = 0; i < n_forefronts; i++) { \
		struct dz_swgv_s const *tdp = (struct dz_swgv_s const *)forefronts[i] - forefronts[i]->r.epos; \
		__m128i const adjv = _mm_set1_epi16(init_s == 0 ? 0 : adj[i]); \
		for(uint64_t p = forefronts[i]->r.spos; p < forefronts[i]->r.epos; p++) { \
			/* adjust offset */ \
			__m128i e = _mm_subs_epu16(_mm_load_si128(&tdp[p].e), adjv); \
			__m128i f = _mm_subs_epu16(_mm_load_si128(&tdp[p].f), adjv); \
			__m128i s = _mm_subs_epu16(_mm_load_si128(&tdp[p].s), adjv); \
			/* read-max-write */ \
			_mm_store_si128(&cdp[p].e, _mm_max_epu16(e, _mm_load_si128(&cdp[p].e))); \
			_mm_store_si128(&cdp[p].f, _mm_max_epu16(f, _mm_load_si128(&cdp[p].f))); \
			_mm_store_si128(&cdp[p].s, _mm_max_epu16(s, _mm_load_si128(&cdp[p].s))); \
			print_vector(s); \
		} \
	} \
	_end_column(cdp, w.r.spos, w.r.epos); \
	cdp; \
})
#define _fill_column(w, pdp, query, rt, rrem, init_s) ({ \
	/* load the base */ \
	_init_rch(query, rt, rrem); \
	struct dz_swgv_s *cdp = _begin_column(w, rch, rrem); \
	/* init vectors */ \
	__m128i f = minv; \
	__m128i ps = _mm_set1_epi16(dz_add_ofs(init_s)); \
	__m128i maxv = _mm_set1_epi16(dz_add_ofs(INT16_MIN)); \
	__m128i const xtv = _mm_set1_epi16(dz_add_ofs(w.inc - self->xt));	/* next offset == current max thus X-drop threshold is always -xt */ \
	__m128i const ofs = _mm_set1_epi16(DZ_SCORE_OFS); \
	/* until the bottommost vertically placed band... */ \
	uint32_t sspos = w.r.spos;					/* save spos on the stack */ \
	for(uint64_t p = w.r.spos; p < w.r.epos; p++) { \
		_load_vector(&pdp[p]); _update_vector(p); \
		if(dz_unlikely(_test_xdrop(s, xtv))) {	/* mark _unlikely to move out of the core loop */ \
			/* drop either leading or trailing vector, skip the forefront extension when the forefront is clipped */ \
			if(p == w.r.spos) { w.r.spos++; cdp--; continue; } else { w.r.epos = p; goto dz_pp_cat(_forefront_, __LINE__); } \
		} \
		_store_vector(&cdp[p]); \
	} \
	/* if reached the forefront of the query sequence, finish the extension */ \
	if(w.r.epos < query->blen) { \
		/* forefront extension; clip the column length if too long */ \
		__m128i e = minv, s = minv; _update_vector(w.r.epos); \
		do { \
			if(_test_xdrop(s, xtv)) { break; } \
			_store_vector(&cdp[w.r.epos]); w.r.epos++; \
			f = _mm_subs_epu16(f, iev8); s = _mm_subs_epu16(s, iev8); \
		} while(w.r.epos < query->blen); \
	} \
dz_pp_cat(_forefront_, __LINE__):; \
	/* create cap object that contains [spos, epos) range (for use in the traceback routine) */ \
	struct dz_cap_s *cap = _end_column(cdp, w.r.spos, w.r.epos); \
	int32_t inc = _hmax_vector(maxv); \
	if(dz_cmp_max(inc, w.inc)) { w.inc = inc; w.mcap = cap; }/* update max; the actual score (absolute score accumulated from the origin) is expressed as max + inc; 2 x cmov */ \
	/* FIXME: rescue overflow */ \
	/* debug("update inc(%u), max(%u, %u, %p), cap(%p)", w.inc, w.max, w.max + w.inc, w.mcap, cap); */ \
	cdp; \
})

/**
 * @fn dz_extend_intl
 */
static __dz_vectorize
struct dz_forefront_s const *dz_extend_intl(
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid,
	uint16_t init_s)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	if(n_forefronts == 0) { return(NULL); }										/* invalid */
	if(rlen == 0 && n_forefronts == 1) { return(forefronts[0]); }				/* no need to merge incoming vectors */

	#ifdef DZ_NUCL_ASCII
		/*
		 * load constants
		 * @ABC_DEFG_HIJK_LMNO
		 * PQRS_TUVW_XYZ
		 */
		static uint8_t const conv_fr[32] __attribute__(( aligned(16) )) = {
			/* ['A' & 0x0f] = A, ['C' & 0x0f] = C, ['G' & 0x0f] = G, ['T' & 0x0f] = T, ['U' & 0x0f] = T, ['N' & 0x0f] = N */
			0, rA, 0, rC, rT, rU, 0, rG, 0, 0, 0, 0, 0, 0, rN, 0,
			0, rT, 0, rG, rA, rA, 0, rC, 0, 0, 0, 0, 0, 0, rN, 0
		};

		/* create conversion table (accessed by indirect %rbp) */
		uint8_t conv[16] __attribute__(( aligned(16) ));
		_mm_store_si128((__m128i *)conv, _mm_load_si128((__m128i const *)&conv_fr[rlen > 0 ? 0 : 16]));
	#endif

	/* constant */
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));

	/* insertion penalties */
	__m128i const iiv = _mm_load_si128((__m128i const *)self->iiv);
	__m128i const iev1 = _mm_load_si128((__m128i const *)self->iev);
	__m128i const iev2 = _mm_add_epi16(iev1, iev1);
	__m128i const iev4 = _mm_slli_epi16(iev1, 2);
	__m128i const iev8 = _mm_slli_epi16(iev1, 3);

	/* deletion penalties */
	__m128i const div = _mm_load_si128((__m128i const *)self->div);
	__m128i const dev = _mm_load_si128((__m128i const *)self->dev);


	struct dz_forefront_s w = { { UINT32_MAX, 0 }, 0, 0, 0, 0, 0, 0, NULL, NULL };	/* uint32_t spos, epos, max, inc; struct dz_query_s const *query; struct dz_cap_s const *cap; */ \
	w.rlen = rlen;
	w.rid = rid;
	w.query = query;
	debug("max(%d), inc(%d), forefront[0](%p)", w.max, w.inc, forefronts[0]);

	/* first iterate over the incoming edge objects to get the current max */
	uint64_t adj[n_forefronts];							/* keep variable length array out of statement expression to avoid a bug of icc */
	struct dz_swgv_s *pdp = _merge_column(w, adj, forefronts, n_forefronts, query, init_s);

	/* fetch the first base */
	int64_t rrem = rlen, dir = rlen < 0 ? 1 : -1;
	uint8_t const *rt = (uint8_t const *)&ref[rrem - (rlen < 0)];

	while(rrem != 0 && w.r.spos < w.r.epos) {
		pdp = _fill_column(w, pdp, query, rt, rrem, init_s); rrem += dir;
		debug("rrem(%d), [%u, %u)", rrem, w.r.spos, w.r.epos);
	}
	return(_end_matrix(pdp, &w, rrem));					/* update max vector (pdp always points at the last vector) */
}

/**
 * @fn dz_extend, dz_scan
 * @brief dz_extend for semi-global alignment (head gaps are penalized for the both directions),
 * and dz_scan for searching anchor (head gap is only penalized for the query side).
 *
 * NOTE: DZ_N_AS_UNMATCHING_BASE is not recommended when dz_scan is used
 */
static __dz_vectorize
struct dz_forefront_s const *dz_extend(
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid)
{
	return(dz_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, INT16_MIN));
}
static __dz_vectorize
struct dz_forefront_s const *dz_scan(
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid)
{
	debug("dz_scan called");
	return(dz_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, 0));
}


#undef _end_matrix
#undef _begin_column_head
#undef _begin_column
#undef _end_column
#undef _load_vector
#undef _update_vector
#undef _store_vector
#undef _hmax_vector
#undef _test_xdrop

#undef _merge_column
#undef _fill_column

/* short, exact matching sequences */
unittest( "extend.base" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 3; trial++) {
		struct dz_query_s const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			:              10
		);
		struct dz_forefront_s const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront == *dz_root(dz));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 5), "max(%u), cell(%u)", forefront->max, dz_ut_sel(2, 2, 5));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 9), "max(%u), cell(%u)", forefront->max, dz_ut_sel(4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28), "max(%u), cell(%u)", forefront->max, dz_ut_sel(9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 34), "max(%u), cell(%u)", forefront->max, dz_ut_sel(11, 11, 34));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 39), "max(%u), cell(%u)", forefront->max, dz_ut_sel(13, 13, 39));

		(void)forefront;
	}
	dz_destroy(dz);
}
#ifndef DZ_PROTEIN
unittest( "extend.base.revcomp" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	char revcomp[dz_unittest_query_length];
	for(size_t i = 0; i < dz_unittest_query_length; i++) {
		switch(dz_unittest_query[dz_unittest_query_length - i - 1]) {
			case 'A': revcomp[i] = 'T'; break;
			case 'C': revcomp[i] = 'G'; break;
			case 'G': revcomp[i] = 'C'; break;
			case 'T': revcomp[i] = 'A'; break;
			case 0: revcomp[i] = 3; break;
			case 1: revcomp[i] = 2; break;
			case 2: revcomp[i] = 1; break;
			case 3: revcomp[i] = 0; break;
		}
	}
	for(size_t trial = 0; trial < 3; trial++) {
		size_t length = trial == 0 ? dz_unittest_query_length
					  : trial == 1 ? 3
					  :              10;
		struct dz_query_s const *q = dz_pack_query_reverse(dz, &revcomp[dz_unittest_query_length - length], length);
		struct dz_forefront_s const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront == *dz_root(dz));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 5));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"AAAATCT"[7], &"\x0\x0\x0\x0\x3\x1\x3"[7], "MASLVQT"), -7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 34));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"GAAAATCT"[8], &"\x2\x0\x0\x0\x0\x3\x1\x3"[8], "MASLVQTG"), -8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 28));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 39));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"TGAAAATCT"[9], &"\x3\x2\x0\x0\x0\x0\x3\x1\x3"[9], "MASLVQTGK"), -9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 28));

		(void)forefront;
	}
	dz_destroy(dz);
}
#endif

/* a small graph */
unittest( "extend.small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 4; trial++) {
		struct dz_query_s const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			: trial == 2 ? 4
			:              7
		);
		struct dz_forefront_s const *forefronts[5] = { NULL };

		/*
		 * AG---TTTT------CTGA
		 *   \ /    \    /
		 *    C      CATT
		 */
		forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 1);
		ut_assert(forefronts[0] != NULL && forefronts[0]->max == dz_ut_sel(4, 4, 9));

		forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "T"), 1, 2);
		ut_assert(forefronts[1] != NULL && forefronts[1]->max == dz_ut_sel(6, 6, 14));
		if(trial == 1) { continue; }

		forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "LVQT"), 4, 3);
		if(trial == 2) {
			ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(8, 8, 18));
			continue;
		}
		ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(14, 14, 32));
		if(trial == 3) { continue; }

		forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATT", "\x1\x0\x3\x3", "CKAK"), 4, 4);
		ut_assert(forefronts[3] != NULL && forefronts[3]->max == dz_ut_sel(22, 22, 43));

		forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "QLTL"), 4, 5);
		ut_assert(forefronts[4] != NULL && forefronts[4]->max == dz_ut_sel(30, 30, 61));
	}
	dz_destroy(dz);
}

/**
 * @fn dz_calc_max_rpos
 */
static __dz_vectorize
int64_t dz_calc_max_rpos(
	struct dz_s *self,
	struct dz_forefront_s const *forefront)
{
	(void)self;
	struct dz_cap_s const *pcap = forefront->mcap;
	int32_t rpos = (pcap == dz_ccap(forefront)) ? 0 : pcap->rrem;
	return((int64_t)rpos);									/* signed expansion */
}

/**
 * @fn dz_calc_max_qpos
 */
static __dz_vectorize
uint64_t dz_calc_max_qpos(
	struct dz_s *self,
	struct dz_forefront_s const *forefront)
{
	(void)self;
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	#define _dp(_cap)					( (struct dz_swgv_s const *)(_cap) - (_cap)->r.epos )

	_init_bonus(forefront->query);
	struct dz_cap_s const *pcap = forefront->mcap;
	__m128i const maxv = _mm_set1_epi16(dz_add_ofs(forefront->inc));
	for(uint64_t p = pcap->r.spos; p < pcap->r.epos; p++) {
		__m128i const s = _add_bonus(p, _mm_load_si128(&_dp(pcap)[p].s)); print_vector(s);
		uint64_t eq = _mm_movemask_epi8(_mm_cmpeq_epi16(s, maxv));
		if(eq != 0) {
			/* tzcntq is faster but avoid using it b/c it requires relatively newer archs */
			uint64_t zcnt = (eq - 1) & (~eq & 0x5555);		/* subq, andnq, andq; chain length == 2 */
			zcnt += zcnt>>2; zcnt &= 0x3333;				/* shrq, addq, andq; chain length might be shorter */
			zcnt += zcnt>>4; zcnt &= 0x0f0f;
			zcnt += zcnt>>8; zcnt &= 0x00ff;
			debug("found, eq(%lx), zcnt(%lu), idx(%lu)", eq, zcnt, p * L + zcnt);
			return(p * L + zcnt);
		}
	}
	return(0);

	#undef _dp
}

/**
 * @fn dz_calc_max_pos
 * @brief rpos (signed int) in upper 32bit, qpos (unsigned int) in lower 32bit
 */
static __dz_vectorize
uint64_t dz_calc_max_pos(
	struct dz_s *self,
	struct dz_forefront_s const *forefront)
{
	struct dz_cap_s const *pcap = forefront->mcap;
	if(pcap == NULL) { debug("pcap == NULL, rlen(%d)", (int32_t)forefront->rlen); return(((uint64_t)forefront->rlen)<<32); }							/* is root */
	debug("forefront(%p), pcap(%p), rrem(%d), rlen(%d)", forefront, pcap, (int32_t)pcap->rrem, (int32_t)forefront->rlen);
	return((dz_calc_max_rpos(self, forefront)<<32) | dz_calc_max_qpos(self, forefront));
}

/* traceback macros */
#define _is_head(_cap)				( (_cap)->r.epos == 0 )
#define _dp(_cap)					( dz_cswgv(_cap) - (_cap)->r.epos )
#define _vector_idx(_idx)			( (_idx) / L )
#define _cell_idx(_idx)				( (_idx) & (L - 1) )
#define _s(_l, _cap, _idx)			( ((uint16_t const *)(&_dp(_cap)[_vector_idx(_idx)]._l))[_cell_idx(_idx)] )
#define _unwind_cap(_c)				( dz_ccap(dz_cswgv(_c) - (_c)->r.epos + (_c)->r.spos) - 1 )
#define _push_span(_id)				{ span->offset = path_base - path; span--; span->id = (_id); }
#define _idx_asc(_i)				( -((int64_t)n_forefronts) + (_i) )
#define _idx_dsc(_i)				( -1LL - (_i) )
#define _load_prev_cap(_l, _score, _idx) ({ \
	while(1) { \
		cap = pcap; pcap = _unwind_cap(cap); \
		debug("idx(%lu), cap(%p), pcap(%p), [%u, %u)", idx, cap, pcap, pcap->r.spos, pcap->r.epos); \
		/* escape if not head */ \
		if(dz_likely(!_is_head(pcap))) { break; } \
		/* escape if root */ \
		debug("head, n_forefronts(%u)", dz_chead(pcap)->n_forefronts); \
		if(pcap->rrem == 0) { debug("reached root"); goto _trace_tail; } \
		/* just load pcap if internal bridge */ \
		struct dz_forefront_s const **forefront_arr = (struct dz_forefront_s const **)pcap; \
		if(pcap->rch != 0xff) { debug("internal bridge, pcap(%p)", dz_ccap(forefront_arr[-1])); pcap = dz_ccap(forefront_arr[-1]); break; } \
		/* merging vector; load contents to find an edge */ \
		uint16_t prev_score = _s(_l, cap, idx) + pcap->r.spos;		/* adj */ \
		size_t n_forefronts = dz_chead(pcap)->n_forefronts; \
		debug("merging vector, %s(%d), n_forefronts(%lu)", #_l, prev_score, n_forefronts); \
		for(size_t i = 0; i < n_forefronts; i++) { \
			struct dz_forefront_s const *ff = forefront_arr[_idx(i)]; \
			pcap = dz_ccap(ff); \
			debug("i(%lu), %s(%d, %d), max(%d), inc(%d), [%u, %u)", i, #_l, prev_score, (uint16_t)(_s(_l, pcap, idx) + (ff->max - ff->inc)), ff->max, ff->inc, ff->r.spos, ff->r.epos); \
			/* adj[i] = w.max - (ffs[i]->max - ffs[i]->inc); base = max - inc */ \
			if(!dz_inside(pcap->r.spos, _vector_idx(idx), pcap->r.epos)) { debug("out of range(%u), [%u, %u)", _vector_idx(idx), pcap->r.spos, pcap->r.epos); continue; } \
			if(prev_score == (uint16_t)(_s(_l, pcap, idx) + (ff->max - ff->inc))) { debug("found, i(%lu), id(%u), prev_score(%d), score(%d), max(%d), inc(%d)", i, dz_cff(pcap)->rid, prev_score, _s(_l, pcap, idx), ff->max, ff->inc); break; } \
		} \
		_push_span(dz_cff(pcap)->rid);								/* push segment info */ \
		_score = _s(_l, pcap, idx); \
		*drp++ = ':'; *dqp++ = ':'; \
	} \
	/* return the reference-side base */ \
	ref_length++; pcap->rch; \
})

/**
 * @fn dz_trace
 */
static __dz_vectorize
struct dz_alignment_s *dz_trace(
	struct dz_s *self,
	struct dz_forefront_s const *forefront)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	if(forefront->mcap == NULL) { return(NULL); }

	/* detect pos */
	uint64_t idx = dz_calc_max_qpos(self, forefront);							/* vector index, cell index */
	uint64_t ref_length = 0, query_length = idx;

	/* allocate aln object */
	size_t aln_size = (sizeof(struct dz_alignment_s)
		+ (forefront->rcnt + 6) * sizeof(struct dz_path_span_s)
		+ dz_roundup(forefront->rsum + idx + 1, 8) * sizeof(uint8_t)			/* +1 for tail '\0' */
	);
	struct dz_alignment_s *aln = (struct dz_alignment_s *)dz_mem_malloc(dz_mem(self), aln_size);
	struct dz_path_span_s *span = (struct dz_path_span_s *)(aln + 1) + forefront->rcnt + 4, *span_base = span;
	uint8_t *path = ((uint8_t *)(span + 2)) + forefront->rsum + idx + 1, *path_base = --path;
	_push_span(forefront->rid); *path = '\0';									/* make sure readable as C string */

	/* load max column pointers */
	struct dz_cap_s const *pcap = forefront->mcap, *cap = NULL;
	struct dz_query_s const *query = forefront->query;
	uint16_t score = _s(s, pcap, idx), cnt[4] = { 0 };

	uint8_t debug_ref[1024], debug_query[1024];
	uint8_t *drp = debug_ref, *dqp = debug_query;

	/* traceback loop */
	#define _debug(_label) { \
		debug("test %s, idx(%lu), rbase(%d, %c), qbase(%c), c(%d), e(%d), f(%d), s(%d), score(%d, %lu), diag(%d)", \
			#_label, idx, rch, '-', '-', /*"ACGTNNNNNNNNNNNN"[rch & 0xff], "ANNNCNNNGNNNTNNN"[query->arr[idx]],*/ \
			dz_rm_ofs(score), dz_rm_ofs(_s(e, cap, idx)), dz_rm_ofs(_s(f, cap, idx)), dz_rm_ofs(_s(s, pcap, idx - 1)), \
			dz_pair_score(self, query, rch, idx), (uint64_t)dz_pair_eq(self, query, rch, idx), \
			(score == (_s(s, pcap, idx - 1) + dz_pair_score(self, query, rch, idx)))); \
	}
	#define _match(_idx) { \
		if(dz_inside(pcap->r.spos, _vector_idx(idx - 1), pcap->r.epos) \
		&& score == (uint16_t)(_s(s, pcap, idx - 1) + dz_pair_score(self, query, rch, idx))) { \
			uint64_t eq = dz_pair_eq(self, query, rch, idx); \
			*--path = DZ_CIGAR_OP>>(eq<<3); cnt[eq]++; \
			score = (uint16_t)_s(s, pcap, idx - 1); idx--; rch = _load_prev_cap(s, score, _idx); \
			continue; \
		} \
	}
	#define _ins(_idx) { \
		if(_vector_idx(idx - 1) >= cap->r.spos && score == (uint16_t)_s(f, cap, idx)) { \
			_debug(I); \
			while(_vector_idx(idx - 1) >= cap->r.spos && score != (uint16_t)(_s(s, cap, idx - 1) - self->iev[0] - self->iiv[0])) { \
				*--path = (DZ_CIGAR_OP>>16) & 0xff; cnt[2]++; score = _s(f, cap, idx - 1); idx--; _debug(I); \
			} \
			*--path = (DZ_CIGAR_OP>>16) & 0xff; cnt[2]++; score = _s(s, cap, idx - 1); idx--; \
			continue; \
		} \
	}
	#define _del(_idx) { \
		if(dz_inside(pcap->r.spos, _vector_idx(idx), pcap->r.epos) && score == (uint16_t)_s(e, cap, idx)) { \
			_debug(D); \
			while(dz_inside(pcap->r.spos, _vector_idx(idx), pcap->r.epos) && score == (uint16_t)(_s(e, pcap, idx) - self->iev[0])) { \
				*--path = (DZ_CIGAR_OP>>24) & 0xff; cnt[3]++; score = _s(e, pcap, idx); _load_prev_cap(e, score, _idx); _debug(D); \
			} \
			*--path = (DZ_CIGAR_OP>>24) & 0xff; cnt[3]++; score = _s(s, pcap, idx); rch = _load_prev_cap(s, score, _idx); \
			continue; \
		} \
	}

	uint32_t rch = _load_prev_cap(s, score, _idx_asc);
	while(1) { _debug(M); _match(_idx_asc); _ins(_idx_asc); _del(_idx_asc); dz_trap(); }
_trace_tail:;
	span++;				/* adjust root span */

	/* finish alignment object */
	uint32_t span_length = span_base - span, path_length = path_base - path;

	assert((uintptr_t)(aln + 1) <= (uintptr_t)span);
	assert((uintptr_t)(&span[span_length + 1]) <= (uintptr_t)path);

	aln->span = span;
	aln->path = path;
	aln->span_length = span_length;
	aln->path_length = path_length;
	aln->ref_length = ref_length;
	aln->query_length = query_length;
	aln->rrem = dz_calc_max_rpos(self, forefront);
	aln->score = forefront->max;
	aln->mismatch_count = cnt[0];
	aln->match_count = cnt[1];
	aln->ins_count = cnt[2];
	aln->del_count = cnt[3];

	/* finish offsets */
	span[0].offset = 0;
	for(size_t i = 1; i < (size_t)aln->span_length; i++) {
		span[i].offset = path_length - span[i].offset;
	}
	span[span_length].offset = path_length;
	return(aln);
}

#undef _is_head
#undef _dp
#undef _s
#undef _vector_idx
#undef _cell_idx
#undef _unwind_cap
#undef _push_span
#undef _idx_asc
#undef _idx_dsc
#undef _load_prev_cap

/* short, exact matching sequences */
unittest( "trace" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefront = NULL;
	struct dz_alignment_s *aln = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, "A", 1, 1);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "GC", 2, 2);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "TTTT", 4, 3);
	aln = dz_trace(dz, forefront);

	(void)aln;
	dz_destroy(dz);
}

#ifdef __cplusplus
};	/* extern "C" { */
#endif

/**
 * end of dozeu.h
 */
