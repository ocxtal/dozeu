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
#include <x86intrin.h>

#ifndef DZ_CIGAR_OP
#  define DZ_CIGAR_OP				0x03040201
#endif
#define DZ_CIGAR_INTL				( (((uint32_t)DZ_CIGAR_OP)>>16) | (((uint32_t)DZ_CIGAR_OP)<<16) )

#ifndef dz_cmp_max
#  define dz_cmp_max(x, y)			( (x) > (y) )
#endif

/* default encoding: ascii */
#if !defined(DZ_PROTEIN) && !defined(DZ_NUCL_ASCII) && !defined(DZ_NUCL_2BIT) && !defined(DZ_NUCL_4BIT)
#  define DZ_NUCL_ASCII
#endif


#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
#define DZ_MAT_SIZE					( 4 )		/* 4 x 4 */
/*
 * define DZ_N_AS_UNMATCHING_BASE to penalize Ns, otherwise scores for (x, N) and (N, x) are always set zero
 */
// #define DZ_N_AS_UNMATCHING_BASE

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

/* matcher; bit[4] for invalid flag */
#define dz_pair_score(_profile, _q, _r, _i)		( (_profile)->matrix[((_r) | (_q)->arr[(_i)]) & 0x1f] - DZ_SCORE_OFS )
#define dz_pair_eq(_profile, _q, _r, _i)		( (uint32_t)((_q)->arr[(_i)]) == ((uint32_t)(_r)<<2) )

#elif defined(DZ_NUCL_4BIT)

#ifndef DZ_MAT_SIZE
#  define DZ_MAT_SIZE					( 16 )		/* 16 x 16 */
#endif

/* 4bit nucleotide (with ambiguous bases) */
enum dz_alphabet_reference {
	rA = 0x01, rC = 0x02, rG = 0x04, rT = 0x08, rU = 0x08,

	rR = rA | rG, rY = rC | rT, rS = rG | rC,
	rW = rA | rT, rK = rG | rT, rM = rA | rC,

	rB = rC | rG | rT, rD = rA | rG | rT,
	rH = rA | rC | rT, rV = rA | rC | rG,

	rN = 0x00
};

enum dz_alphabet_query {
	qA = 0x01, qC = 0x02, qG = 0x04, qT = 0x08, qU = 0x08,

	qR = qA | qG, qY = qC | qT, qS = qG | qC,
	qW = qA | qT, qK = qG | qT, qM = qA | qC,

	qB = qC | qG | qT, qD = qA | qG | qT,
	qH = qA | qC | qT, qV = qA | qC | qG,

	qN = 0x00, qX = 0x00
};

/* matcher; bit[4] for invalid flag */
#define dz_pair_score(_profile, _q, _r, _i)		( (_profile)->matrix[(_r) * DZ_MAT_SIZE + (_q)->arr[(_i)]] - DZ_SCORE_OFS )
#define dz_pair_eq(_profile, _q, _r, _i)		( (uint32_t)((_q)->arr[(_i)]) == (uint32_t)(_r) )

#elif defined(DZ_PROTEIN)

/* protein */
#ifndef DZ_MAT_SIZE
#  define DZ_MAT_SIZE				( 32 )		/* 32 x 32 */
#endif
#define dz_pair_score(_profile, _q, _r, _i)		( (int8_t)((_q)->arr[(_r) * (_q)->blen * DZ_L + (_i)]) - DZ_SCORE_OFS )
#define dz_pair_eq(_profile, _q, _r, _i)		( (uint32_t)((_q)->q[(_i) - 1] & 0x1f) == (uint32_t)(_r) )
#endif


/* utils */
#define dz_pp_cat_intl(x, y)		x##y
#define dz_pp_cat(x, y)				dz_pp_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_pp_cat(_st_, __LINE__)[(expr) ? 1 : -1]
#define dz_trap()					{ *((volatile uint8_t *)NULL); }
#ifdef DZ_PROTEIN
	dz_static_assert(DZ_MAT_SIZE <= 32);
#endif

#if (defined(DEBUG) || (defined(UNITTEST) && UNITTEST != 0)) && !defined(__cplusplus)
#  include "log.h"
#  if !defined(UNITTEST_UNIQUE_ID)
#    define UNITTEST_ALIAS_MAIN		0
#    define UNITTEST_UNIQUE_ID		3213
#  endif
#  include "unittest.h"
unittest_config( .name = "dozeu" );
unittest() { debug("hello"); }
#else
#  ifndef UNITTEST_H_INCLUDED
#    define unittest(...)				static void dz_pp_cat(dz_unused_, __LINE__)(void)
#    define ut_assert(...)			;
#    define debug(...)				;
#    define trap()					;
#  endif
#endif


#if defined(DZ_NUCL_ASCII)
#  define DZ_UNITTEST_INDEX			0
#elif defined(DZ_NUCL_2BIT)
#  define DZ_UNITTEST_INDEX			1
#elif defined(DZ_NUCL_4BIT)
#  define DZ_UNITTEST_INDEX			2
#elif defined(DZ_PROTEIN)
#  define DZ_UNITTEST_INDEX			3
#endif
#define dz_ut_sel(a, b, c, d)		( (DZ_UNITTEST_INDEX == 0) ? (a) : (DZ_UNITTEST_INDEX == 1) ? (b) : (DZ_UNITTEST_INDEX == 2) ? (c) : (d) )

/* vectorize */
#ifndef __x86_64__
#  error "x86_64 is required"
#endif
#ifndef __SSE4_1__
#  warning "SSE4.1 is automatically enabled in dozeu.h, please check compatibility to the system."
#  define __dz_vectorize			__attribute__(( target( "sse4.1" ) )) inline
#else
#  define __dz_vectorize			/* follow the compiler options */
#endif

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

#define dz_add_ofs(_x)				( (int16_t)((uint16_t)(_x) ^ (uint16_t)0x8000) )
#define dz_rm_ofs(_x)				( (int16_t)((uint16_t)(_x) ^ (uint16_t)0x8000) )

#define dz_add_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) + (size_t)(_ofs)) )
#define dz_sub_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) - (size_t)(_ofs)) )


/* print_vector for debugging */
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
#define print_vector_raw(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	_mm_extract_epi16(v, 7), \
	_mm_extract_epi16(v, 6), \
	_mm_extract_epi16(v, 5), \
	_mm_extract_epi16(v, 4), \
	_mm_extract_epi16(v, 3), \
	_mm_extract_epi16(v, 2), \
	_mm_extract_epi16(v, 1), \
	_mm_extract_epi16(v, 0)); \
}

#define print_vector_8(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	_mm_extract_epi8(v, 8 + 7), \
	_mm_extract_epi8(v, 8 + 6), \
	_mm_extract_epi8(v, 8 + 5), \
	_mm_extract_epi8(v, 8 + 4), \
	_mm_extract_epi8(v, 8 + 3), \
	_mm_extract_epi8(v, 8 + 2), \
	_mm_extract_epi8(v, 8 + 1), \
	_mm_extract_epi8(v, 8 + 0), \
	_mm_extract_epi8(v, 7), \
	_mm_extract_epi8(v, 6), \
	_mm_extract_epi8(v, 5), \
	_mm_extract_epi8(v, 4), \
	_mm_extract_epi8(v, 3), \
	_mm_extract_epi8(v, 2), \
	_mm_extract_epi8(v, 1), \
	_mm_extract_epi8(v, 0)); \
}


#else
#define print_vector(v) ;
#define print_vector_raw(v) ;

#endif



#ifdef __SSE4_1__
#  define dz_is_all_zero(x)			( _mm_test_all_zeros((x), (x)) == 1 )
#else
// #  define dz_is_all_zero(x)			( _mm_movemask_epi8((x)) == 0 )
#endif

#define DZ_MEM_MARGIN_SIZE			( 256 )
#define DZ_MEM_ALIGN_SIZE			( 16 )
#define DZ_MEM_INIT_SIZE			( 16 * 1024 * 1024 )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_CELL_MARGIN				( 32 )
#define DZ_CELL_MARGINED_MIN		( DZ_CELL_MIN + DZ_CELL_MARGIN )
#define DZ_CELL_MARGINED_MAX		( DZ_CELL_MAX - DZ_CELL_MARGIN )
#define DZ_L						( sizeof(__m128i) / sizeof(uint16_t) )

#define DZ_SCORE_OFS				( 64 )
#define DZ_HEAD_RCH					( 0x40 )
#define DZ_ROOT_RCH					( 0x80 )


/* query; preconverted query sequence; blen = roundup(qlen, DZ_L) / DZ_L; array must have 16-byte-length margin at the tail */
typedef struct dz_query_s {
	uint64_t blen;		/* #vectors */
	char const *q;		/* query string (encoded and aligned) */
	int16_t bonus[2 * DZ_L];
	uint8_t arr[];
} dz_query_t;
dz_static_assert(sizeof(dz_query_t) % sizeof(__m128i) == 0);

/* full-length bonus */
static __dz_vectorize
__m128i dz_query_add_bonus(dz_query_t const *query, __m128i v, size_t p)
{
	#ifdef DZ_FULL_LENGTH_BONUS
		__m128i const *pbonus = &((__m128i const *)query->arr)[-2];
		__m128i bv = _mm_load_si128(&pbonus[p == query->blen - 1]);
		v = _mm_add_epi16(v, bv);

	#else
		dz_unused(query);
		dz_unused(p);

	#endif
	return(v);
}


/* node (reference; reserved for realigner implementation) */
typedef struct dz_node_s {
	int32_t id, len;
	uint8_t const *ptr;
} dz_node_t;
dz_static_assert(sizeof(dz_node_t) % sizeof(__m128i) == 0);



/* internal memory allocator */
typedef struct dz_arena_block_s {
	struct dz_arena_block_s *next;
	size_t size;
} dz_arena_block_t;

typedef struct dz_stack_s {
	dz_arena_block_t *curr;
	uint8_t *top, *end;
	uint64_t _pad[3];
} dz_stack_t;

typedef struct dz_arena_s {
	dz_arena_block_t blk;
	dz_stack_t stack;
} dz_arena_t;
#define dz_arena_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )


/**
 * @fn dz_malloc, dz_free
 * @brief aligned and margined malloc and free
 */
static __dz_force_inline
void *dz_malloc(size_t size)
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
void dz_free(void *ptr)
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
 * @fn dz_arena_init, dz_arena_destroy, dz_arena_add_stack, dz_arena_malloc, dz_arena_flush
 * @brief stack chain
 */
static __dz_force_inline
void dz_arena_flush(dz_arena_t *mem)
{
	mem->stack.curr = &mem->blk;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem + 1), DZ_MEM_ALIGN_SIZE);
	mem->stack.end = (uint8_t *)mem + mem->blk.size - DZ_MEM_MARGIN_SIZE;
	return;
}

static __dz_force_inline
dz_arena_t *dz_arena_init(size_t size)
{
	size = dz_min2(DZ_MEM_INIT_SIZE, sizeof(dz_arena_t) + dz_roundup(size, DZ_MEM_ALIGN_SIZE));
	dz_arena_t *mem = (dz_arena_t *)dz_malloc(size);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init mem object then stack pointers */
	mem->blk = (dz_arena_block_t){ .next = NULL, .size = size };
	dz_arena_flush(mem);
	return(mem);
}

static __dz_force_inline
void dz_arena_destroy(dz_arena_t *mem)
{
	dz_arena_block_t *blk = mem->blk.next;
	debug("cleanup memory chain, blk(%p)", blk);
	while(blk != NULL) {
		dz_arena_block_t *next = blk->next;
		debug("free blk(%p), next(%p)", blk, next);
		dz_free(blk); blk = next;
	}
	dz_free(mem);
	return;
}
static __dz_force_inline
uint64_t dz_arena_add_stack(dz_arena_t *mem, size_t size)
{
	debug("add_stack, ptr(%p)", mem->stack.curr->next);
	if(mem->stack.curr->next == NULL) {
		/* current stack is the forefront of the memory block chain, add new block */
		size = dz_max2(
			size + dz_roundup(sizeof(dz_arena_block_t), DZ_MEM_ALIGN_SIZE),
			2 * mem->stack.curr->size
		);
		dz_arena_block_t *blk = (dz_arena_block_t *)dz_malloc(size);
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
void *dz_arena_malloc(dz_arena_t *mem, size_t size)
{
	if(dz_arena_stack_rem(mem) < 4096) { dz_arena_add_stack(mem, 0); }
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += dz_roundup(size, sizeof(__m128i));
	return(ptr);
}

unittest() {
	size_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(size_t i = 0; i < sizeof(size) / sizeof(size_t); i++) {
		dz_arena_t *mem = dz_arena_init(size[i]);
		uint8_t *p = (uint8_t *)dz_arena_malloc(mem, size[i]);
		ut_assert(p != NULL);
		memset(p, 0, size[i]);
		dz_arena_destroy(mem);
	}
}

static __dz_force_inline
uint8_t *dz_reserve_stack(dz_arena_t *mem, size_t size)
{
	/* allocate from heap */
	if(dz_arena_stack_rem(mem) < size) {
		dz_arena_add_stack(mem, 0);
	}
	return(mem->stack.top);
}

static __dz_force_inline
void dz_save_stack(dz_arena_t *mem, void *ptr)
{
	mem->stack.top = (uint8_t *)ptr;
	return;
}



/* score profile object */
typedef struct dz_profile_s {
	uint16_t iiv[8], iev[8];
	uint16_t div[8], dev[8];	/* deletion: open and extend */

	/* root column */
	struct dz_forefront_s const *root;

	/* constants */
	uint16_t xt, bonus;
	uint16_t max_gap_len;
	uint16_t size;

	int8_t matrix[];
} dz_profile_t;
dz_static_assert(sizeof(dz_profile_t) % sizeof(__m128i) == 0);
// #define dz_arena(_self)				( (dz_arena_t *)(_self) - 1 )


/* DP matrix internal objects */

enum dz_layer {
	DZ_E_MATRIX = 0,
	DZ_F_MATRIX = 1,
	DZ_S_MATRIX = 2
};

/* score matrix block; followed by dz_cap_s */
typedef struct dz_swgv_s {
	__m128i e, f, s;	/* horizontal gap, vertical gap, and score */
} dz_swgv_t;
dz_static_assert(sizeof(dz_swgv_t) % sizeof(__m128i) == 0);
#define dz_swgv(_p)					( (dz_swgv_t *)(_p) )
#define dz_cswgv(_p)				( (dz_swgv_t const *)(_p) )


static __dz_vectorize
dz_swgv_t dz_load_swgv(dz_swgv_t const *p)
{
	dz_swgv_t v;
	v.e = _mm_loadu_si128(&p->e);
	v.f = _mm_loadu_si128(&p->f);
	v.s = _mm_loadu_si128(&p->s);
	return(v);
}

static __dz_vectorize
dz_swgv_t *dz_store_swgv(dz_swgv_t *p, dz_swgv_t v)
{
	_mm_storeu_si128(&p->e, v.e);
	_mm_storeu_si128(&p->f, v.f);
	_mm_storeu_si128(&p->s, v.s);
	return(p);
}

static __dz_vectorize
dz_swgv_t dz_subs_swgv(dz_swgv_t v, __m128i ofs)
{
	v.e = _mm_subs_epu16(v.e, ofs);
	v.f = _mm_subs_epu16(v.f, ofs);
	v.s = _mm_subs_epu16(v.s, ofs);
	return(v);
}

static __dz_vectorize
dz_swgv_t dz_max_swgv(dz_swgv_t v, dz_swgv_t w)
{
	return((dz_swgv_t){
		.e = _mm_max_epu16(v.e, w.e),
		.f = _mm_max_epu16(v.f, w.f),
		.s = _mm_max_epu16(v.s, w.s)
	});
}



/* sequence fetcher comes here */
typedef struct {
	char const *ptr;
	uint32_t len;
	uint32_t id;

	int32_t dir;
	int32_t init_s;
} dz_ref_t;

typedef struct {
	uint32_t id;
	int32_t len;
} dz_ref_save_t;

typedef struct {
	uint32_t ch;		/* fetched char (after conversion) */
	int32_t rem;		/* remaining reference length */
} dz_ref_state_t;

typedef struct {
	/* encoded dz_ref_t */
	uint8_t const *rt;
	int32_t dir;
	uint16_t pad2, init_s;

	/* encoded reference base */
	#if defined(DZ_NUCL_ASCII)
		__m128i rv, matrix;
		uint8_t conv[16];

	#elif defined(DZ_NUCL_2BIT)
		__m128i rv, matrix;

	#elif defined(DZ_NUCL_4BIT)
		__m128i matrix;
		int8_t const *score_matrix;
		uint64_t conv;

	#elif defined(DZ_PROTEIN)
		int8_t const *parr;
		uint64_t pad1;
	#endif

	/* query */
	dz_query_t const *query;

	dz_ref_save_t rsave;
	dz_ref_state_t rstate;
} dz_ref_fetcher_t;
dz_static_assert((offsetof(dz_ref_fetcher_t, rt) % DZ_L) == 0);

static __dz_vectorize
dz_ref_fetcher_t dz_init_fetcher(dz_profile_t const *profile, dz_ref_t const *ref, dz_query_t const *query)
{
	dz_ref_fetcher_t fetcher;

	#if defined(DZ_NUCL_ASCII)
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
		_mm_store_si128((__m128i *)&fetcher.conv[0],
			_mm_load_si128((__m128i const *)&conv_fr[ref->dir < 0 ? 16 : 0])
		);
		fetcher.matrix = _mm_loadu_si128((__m128i const *)profile->matrix);

	#elif defined(DZ_NUCL_2BIT)
		fetcher.matrix = _mm_loadu_si128((__m128i const *)profile->matrix);

	#elif defined(DZ_NUCL_4BIT)
		fetcher.score_matrix = profile->matrix;
		fetcher.conv = ref->dir < 0 ? 0xf7b3d591e6a2c480 : 0xfedcba9876543210;

	#else
		dz_unused(profile);

	#endif

	/* save query pointer */
	fetcher.query = query;

	/* load length and pointer */
	size_t tail_pos = ref->dir < 0 ? (-1LL - ref->len) : ref->len;
	fetcher.rt     = (uint8_t const *)&ref->ptr[tail_pos];
	fetcher.dir    = ref->dir;
	fetcher.init_s = dz_add_ofs(ref->init_s);

	/* save metadata */
	fetcher.rsave.id  = ref->id;
	fetcher.rsave.len = ref->dir * ref->len;

	/* init state */
	fetcher.rstate.ch  = DZ_HEAD_RCH;
	fetcher.rstate.rem = ref->dir * ref->len;
	return(fetcher);
}

static __dz_vectorize
uint64_t dz_fetch_next(dz_ref_fetcher_t *fetcher)
{
	if(fetcher->rstate.rem == 0) {
		return(1);
	}

	#if defined(DZ_NUCL_ASCII)
		uint32_t const c = fetcher->rt[-fetcher->rstate.rem];
		uint32_t const e = fetcher->conv[c & 0x0f];

		fetcher->rstate.ch = e;
		fetcher->rv = _mm_set1_epi8(e);
		debug("ch(%c, %x), rem(%d), dir(%d)", c, e, fetcher->rstate.rem, fetcher->dir);

	#elif defined(DZ_NUCL_2BIT)
		uint32_t const c = fetcher->rt[-fetcher->rstate.rem];
		uint32_t const e = c ^ ((fetcher->dir>>1) & 0x03);

		fetcher->rstate.ch = e;
		fetcher->rv = _mm_set1_epi8(e);

	#elif defined(DZ_NUCL_4BIT)
		uint32_t const c =  fetcher->rt[-fetcher->rstate.rem];
		uint32_t const e = (fetcher->conv>>(4 * c)) & 0x0f;
		fetcher->rstate.ch = e;
		fetcher->matrix = _mm_load_si128((__m128i const *)&fetcher->score_matrix[e * DZ_MAT_SIZE]);
		debug("ch(%x, %x), rem(%d), dir(%d)", c, e, fetcher->rstate.rem, fetcher->dir);
		// print_vector_8(fetcher->matrix);

	#elif defined(DZ_PROTEIN)
		uint32_t const c = fetcher->rt[-fetcher->rstate.rem] & 0x1f;
		fetcher->rstate.ch = c;
		fetcher->parr = (int8_t const *)&fetcher->query->arr[c * fetcher->query->blen * DZ_L];
		debug("ch(%x), rem(%d), dir(%d)", c, fetcher->rstate.rem, fetcher->dir);

	#endif

	fetcher->rstate.rem -= fetcher->dir;
	return(0);			/* 1 if reached tail */
}

static __dz_vectorize
__m128i dz_calc_score_profile(dz_ref_fetcher_t *fetcher, size_t p)
{
	#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
		__m128i const qv = _mm_loadl_epi64(
			(__m128i const *)&fetcher->query->arr[p * DZ_L]
		);
		__m128i const sc = _mm_shuffle_epi8(
			fetcher->matrix,
			_mm_or_si128(fetcher->rv, qv)
		);
		__m128i const v = _mm_cvtepi8_epi16(sc);

		// print_vector_raw(_mm_cvtepi8_epi16(qv));
		// print_vector_raw(_mm_cvtepi8_epi16(sc));
		// print_vector_raw(v);

	#elif defined(DZ_NUCL_4BIT)
		__m128i const qv = _mm_loadl_epi64(
			(__m128i const *)&fetcher->query->arr[p * DZ_L]
		);
		__m128i const sc = _mm_shuffle_epi8(fetcher->matrix, qv);
		__m128i const v = _mm_cvtepi8_epi16(sc);

		print_vector_raw(_mm_cvtepi8_epi16(qv));
		print_vector_raw(_mm_cvtepi8_epi16(sc));
		print_vector_raw(v);

	#elif defined(DZ_PROTEIN)
		__m128i const v = _mm_cvtepi8_epi16(
			_mm_loadl_epi64((__m128i const *)&fetcher->parr[p * DZ_L])
		);
	#endif
	return(v);
}

static __dz_vectorize
size_t dz_finalize_fetcher(dz_ref_fetcher_t const *fetcher)
{
	int64_t dir = fetcher->dir;
	size_t fwd  = dir * (fetcher->rsave.len - fetcher->rstate.rem);
	return(fwd);		/* forwarded columns */
}


/* scan */
static __dz_vectorize
__m128i dz_fetch_init_vector(dz_ref_fetcher_t const *fetcher)
{
	__m128i v = _mm_load_si128((__m128i const *)&(&fetcher->init_s)[-7]);
	return(v);
}

static __dz_vectorize
uint64_t dz_is_scan(dz_ref_fetcher_t const *fetcher)
{
	return(fetcher->init_s == (uint16_t)dz_add_ofs(0));
}


/* DP matrix internal objects */

/* placed just after every score vector to indicate the length */
typedef struct dz_range_s {
	uint32_t spos, epos;
} dz_range_t;
#define dz_range(_p)				( (dz_range_t *)(_p) )
#define dz_crange(_p)				( (dz_range_t const *)(_p) )

typedef struct {
	uint32_t column;	/* column count */
	uint32_t section;	/* section count */
} dz_ref_cnt_t;

typedef struct {
	int32_t score, inc;
	struct dz_cap_s const *cap;
} dz_max_t;

/* as working buffer */
/* typedef */
typedef struct dz_state_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;
	dz_max_t max;
} dz_state_t;
#define dz_state(_p)				( (dz_state_t *)(_p) )
#define dz_cstate(_p)				( (dz_state_t const *)(_p) )
#define dz_pstate(_p)				( (dz_state_t **)(_p) )
#define dz_cpstate(_p)				( (dz_state_t const **)(_p) )


/* struct dz_forefront_s is alias of dz_state_t */
typedef struct dz_forefront_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;
	int32_t max, inc;
	struct dz_cap_s const *cap;
} dz_forefront_t;
dz_static_assert(sizeof(dz_forefront_t) == sizeof(dz_state_t));
dz_static_assert(offsetof(dz_forefront_t, max) == offsetof(dz_state_t, max.score));

#define dz_cff(_p)					( (dz_forefront_t const *)(_p) )
#define dz_cpff(_p)					( (dz_forefront_t const **)(_p) )

/* control */
#define dz_profile_root(_profile)	( dz_cpff(&(_profile)->root) )
#define dz_root(_self)				( dz_profile_root((_self)->profile) )


static __dz_vectorize
dz_state_t dz_merge_state(dz_state_t const **ff, size_t fcnt)
{
	/* load constant */
	static int32_t const neg[4] __attribute__(( aligned(16) )) = { -1, 0, 0, 0 };
	__m128i const nv = _mm_load_si128((__m128i const *)neg);

	/* working registers */
	__m128i mv0 = nv;
	__m128i mv1 = _mm_setzero_si128();

	for(size_t i = 0; i < fcnt; i++) {
		/* spos, epos, ccnt, scnt */
		__m128i const v0 = _mm_loadu_si128((__m128i const *)&ff[i]->range.spos);
		mv0 = _mm_max_epu32(mv0, _mm_xor_si128(v0, nv));

		/* max */
		__m128i const v1 = _mm_loadu_si128((__m128i const *)&ff[i]->max.score);
		mv1 = _mm_max_epu32(mv1, v1);

		debug("i(%zu), fcnt(%zu), ff(%p), range(%u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
			i, fcnt, ff[i], ff[i]->range.spos, ff[i]->range.epos, ff[i]->cnt.column, ff[i]->cnt.section, ff[i]->max.score, ff[i]->max.inc
		);
	}

	/* save */
	dz_state_t state __attribute__(( aligned(16) ));

	/* spos, epos, column, section */
	_mm_store_si128((__m128i *)&state.range.spos,
		_mm_xor_si128(mv0, nv)
	);

	/* score, inc = 0, cap = NULL */
	_mm_store_si128((__m128i *)&state.max.score,
		_mm_and_si128(mv1, nv)
	);
	debug("merged state, range(%u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
		state.range.spos, state.range.epos, state.cnt.column, state.cnt.section, state.max.score, state.max.inc
	);
	return(state);
}

static __dz_vectorize
void dz_finalize_state(dz_state_t *state, size_t cols)
{
	state->cnt.column += cols;
	state->cnt.section++;
	state->max.score += state->max.inc;
	debug("score(%d), inc(%d)", state->max.score, state->max.inc);
	return;
}


/* working buffer */
typedef struct {
	/* constants */
	__m128i minv;		/* dz_add_ofs(DZ_CELL_MIN) */
	__m128i ofsv;

	/* working vectors */
	__m128i iiv;		/* gap (insertion) open penalty */
	__m128i iev1, iev2, iev4, iev8;		/* insertion extension penalty x { 1, 2, 4, 8 } for vertical chaining */
	__m128i div, dev;	/* gap (deletion) open and extension penalty */

	/* constants (contd; X-drop threshold) */
	uint16_t xt, _pad[3];

	/* mem stack */
	dz_arena_t *mem;

	/* incoming matrices */
	struct {
		dz_state_t const **ptr;
		size_t cnt;
	} incoming;

	/* sequence (query and reference) info */
	dz_ref_fetcher_t fetcher;

	/* states */
	dz_state_t state;

	/* sequence info (continued) */
	dz_query_t const *query;
} dz_work_t;

static __dz_vectorize
void dz_load_state(dz_work_t *w, dz_state_t const **ff, size_t fcnt)
{
	/* save pointer */
	w->incoming.ptr = ff;
	w->incoming.cnt = fcnt;

	/* merge state (before merging vectors) */
	w->state = dz_merge_state(ff, fcnt);
	return;
}

static __dz_vectorize
void dz_load_profile(dz_work_t *w, dz_profile_t const *profile)
{
	/* constant */
	w->minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));
	w->ofsv = _mm_set1_epi16(DZ_SCORE_OFS);

	/* insertion penalties */
	__m128i const iev = _mm_load_si128((__m128i const *)profile->iev);
	w->iiv  = _mm_load_si128((__m128i const *)profile->iiv);
	w->iev1 = iev;
	w->iev2 = _mm_add_epi16(iev, iev);
	w->iev4 = _mm_slli_epi16(iev, 2);
	w->iev8 = _mm_slli_epi16(iev, 3);

	/* deletion penalties */
	w->div  = _mm_load_si128((__m128i const *)profile->div);
	w->dev  = _mm_load_si128((__m128i const *)profile->dev);

	/* X-drop threshold */
	w->xt = profile->xt;
	return;
}



/* DP matrix scaffolds */

/* placed at the head */
typedef struct dz_head_s {
	uint32_t magic;		/* 0xff for merged head */
	uint32_t fcnt;		/* #merged vectors (n_forefronts) */
	int32_t adj;		/* always positive */
	uint32_t zero;		/* always zero */
} dz_head_t;
#define dz_head(_p)					( (dz_head_t *)(_p) )
#define dz_chead(_p)				( (dz_head_t const *)(_p) )

/* followed by dz_tail_t; spos and epos are shared to forefront_s */
typedef struct dz_cap_s {
	dz_ref_state_t rstate;
	dz_range_t range;	/* column range */
} dz_cap_t;
dz_static_assert(sizeof(dz_cap_t) % sizeof(__m128i) == 0);
#define dz_cap(_p)					( (dz_cap_t *)(_p) )
#define dz_ccap(_p)					( (dz_cap_t const *)(_p) )

/* dz_head_t and dz_cap_t are compatibile for head two elements */
dz_static_assert(offsetof(dz_head_t, magic) == offsetof(dz_cap_t, rstate.ch));
dz_static_assert(offsetof(dz_head_t, zero)  == offsetof(dz_cap_t, range.epos));

/* DP matrix object (exported) */
typedef struct dz_tail_s {
	dz_ref_save_t rsave;
	dz_state_t state;
	dz_query_t const *query;
} dz_tail_t;
dz_static_assert(sizeof(dz_tail_t) % sizeof(__m128i) == 0);
#define dz_tail(_p)			( (dz_tail_t *)(_p) )
#define dz_ctail(_p)		( (dz_tail_t const *)(_p) )


static __dz_force_inline
size_t dz_calc_column_size(size_t spos, size_t epos)
{
	size_t const span = epos - spos;
	size_t const column_size = sizeof(dz_swgv_t) * (span + dz_max2(span, 8));
	size_t const size = (
		  sizeof(dz_cap_t)	/* header size */
		+ column_size		/* estimated column size */
		+ sizeof(dz_state_t)	/* tail cap */
	);
	return(size);
}

static __dz_vectorize
dz_swgv_t *dz_init_root_head(dz_head_t *head)
{
	head->magic = DZ_ROOT_RCH | DZ_HEAD_RCH;
	head->fcnt  = 0;		/* is_root */
	head->adj   = 0;
	head->zero  = 0;		/* is_head */
	return(dz_swgv(head + 1));
}

static __dz_vectorize
uint64_t dz_is_root(dz_cap_t const *cap)
{
	/* test for root tail */
	// return(cap->range.epos == 0);
	return((dz_chead(cap)->magic & DZ_ROOT_RCH) != 0);
}

static __dz_vectorize
dz_cap_t *dz_slice_head(dz_work_t *w)
{
	/* forefront pointer array size */
	size_t const fsize = dz_roundup(
		sizeof(dz_state_t *) * w->incoming.cnt,
		sizeof(__m128i)
	);

	/* column size (including head and internal caps) */
	size_t const csize = dz_calc_column_size(w->state.range.spos, w->state.range.epos);

	/* allocate mem */
	dz_head_t *head = dz_head(dz_add_ptr(dz_reserve_stack(w->mem, fsize + csize), fsize));

	/* fill magic numbers */
	head->magic = DZ_HEAD_RCH;
	head->zero  = 0;		/* overlaps range.epos */

	/* push head-cap info */
	head->fcnt  = w->incoming.cnt;
	head->adj   = w->state.max.score;	/* save merging adjustment */
	return(dz_cap(head));
}

static __dz_vectorize
uint64_t dz_is_head(dz_cap_t const *cap)
{
	return((dz_chead(cap)->magic & DZ_HEAD_RCH) != 0);
}

static __dz_vectorize
void dz_save_incoming(dz_cap_t *head, dz_state_t const **ff, size_t fcnt)
{
	dz_state_t const **q = dz_cpstate(head) - fcnt;
	for(size_t i = 0; i < fcnt; i++) {
		q[i] = ff[i];
		debug("i(%zu), fcnt(%zu), ff(%p)", i, fcnt, ff[i]);
	}
	return;
}


/* ping-pong pointer conversion between dz_cap_t and dz_swgv_t */
static __dz_vectorize
dz_swgv_t *dz_slice_column(dz_cap_t *prev_cap, size_t spos)
{
	dz_swgv_t *col = dz_swgv(prev_cap + 1);
	return(col - spos);
}

static __dz_vectorize
dz_cap_t *dz_cap_column(dz_swgv_t *col, size_t epos)
{
	/* put cap at the tail of current column */
	return(dz_cap(&col[epos]));
}

static __dz_vectorize
dz_swgv_t const *dz_restore_column(dz_cap_t const *cap)
{
	return(dz_cswgv(cap) - cap->range.epos);
}

static __dz_vectorize
dz_cap_t const *dz_unwind_cap(dz_cap_t const *cap)
{
	dz_swgv_t const *col = dz_restore_column(cap);
	debug("cap(%p), range(%u, %u), col(%p), prev_cap(%p)", cap, cap->range.spos, cap->range.epos, col, dz_ccap(col + cap->range.spos) - 1);
	return(dz_ccap(col + cap->range.spos) - 1);
}


static __dz_vectorize
dz_state_t const *dz_slice_tail(dz_work_t *w, dz_swgv_t *col)
{
	dz_tail_t *tail = dz_tail(dz_cap_column(col, w->state.range.epos));
	debug("col(%p), range(%u, %u), tail(%p, %p)", col, w->state.range.spos, w->state.range.epos, tail, &tail->state);

	/* just copy */
	tail->rsave = w->fetcher.rsave;
	tail->state = w->state;
	tail->query = w->fetcher.query;

	/* write back stack pointer */
	dz_save_stack(w->mem, tail + 1);
	return(&tail->state);
}

static __dz_vectorize
dz_tail_t const *dz_restore_tail(dz_state_t const *ff)
{
	return(dz_ctail(dz_sub_ptr(ff, offsetof(dz_tail_t, state))));
}

static __dz_vectorize
dz_swgv_t const *dz_restore_tail_column(dz_state_t const *ff)
{
	dz_tail_t const *tail = dz_restore_tail(ff);
	debug("tail(%p, %p), range(%u, %u), col(%p)", tail, ff, ff->range.spos, ff->range.epos, dz_cswgv(tail) - tail->state.range.epos);
	return(dz_cswgv(tail) - tail->state.range.epos);
}


/* internal link (ilink) */
static __dz_vectorize
dz_swgv_t *dz_slice_ilink(dz_work_t *w, dz_cap_t const *prev_cap, uint8_t *ptr)
{
	/* save link */
	dz_state_t const **ff = ((dz_state_t const **)ptr) + 2;
	ff[-1] = dz_state(&prev_cap->range);

	/* save head info */
	dz_head_t *head = dz_head(ff);
	head->magic = w->fetcher.rstate.ch;
	head->zero  = 0;

	head->fcnt  = 1;		/* #incoming vectors == 1 */
	head->adj   = 0;
	return(dz_swgv(head + 1));
}

static __dz_vectorize
uint64_t dz_is_ilink(dz_cap_t const *cap)
{
	return(cap->rstate.ch != DZ_HEAD_RCH);
}

static __dz_vectorize
dz_cap_t const *dz_rewind_ilink(dz_cap_t const *cap)
{
	dz_state_t const **ff = dz_cpstate(cap);
	return(dz_ccap(dz_restore_tail(ff[-1])));
}

static __dz_vectorize
dz_swgv_t *dz_slice_cap_core(dz_work_t *w, dz_swgv_t *prev_col)
{
	dz_cap_t *cap = dz_cap_column(prev_col, w->state.range.epos);
	debug("prev_col(%p), cap(%p), range(%u, %u)", prev_col, cap, w->state.range.spos, w->state.range.epos);

	/* save states */
	cap->rstate = w->fetcher.rstate;
	cap->range  = w->state.range;
	dz_save_stack(w->mem, cap + 1);

	/* add stack if memory starved */
	size_t size = 2 * sizeof(dz_state_t *);
	size += dz_calc_column_size(w->state.range.spos, w->state.range.epos);

	/* allocate memory from stack */
	uint8_t *ptr = dz_reserve_stack(w->mem, size);
	debug("cap(%p, %p)", cap, dz_cap(ptr) - 1);
	if(dz_likely(dz_cap(ptr) == cap + 1)) { return(dz_swgv(cap + 1)); }

	/* new stack allocated */
	return(dz_slice_ilink(w, cap, ptr));
}

static __dz_vectorize
dz_swgv_t *dz_slice_cap(dz_work_t *w, dz_swgv_t *prev_col)
{
	dz_swgv_t *col = dz_slice_cap_core(w, prev_col);
	return(col - w->state.range.spos);
}


/* merge incoming vectors */
static __dz_vectorize
void dz_merge_init_col(dz_work_t *w, dz_swgv_t *col)
{
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(INT16_MIN));
	dz_swgv_t const v = { minv, minv, minv };

	for(uint64_t p = w->state.range.spos; p < w->state.range.epos; p++) {
		dz_store_swgv(&col[p], v);
	}
	return;
}

static __dz_vectorize
int16_t dz_merge_calc_adj(dz_work_t *w, dz_state_t const *ff)
{
	return(dz_is_scan(&w->fetcher)
		? 0
		: (w->state.max.score - ff->max.score + ff->max.inc)
	);
}

static __dz_vectorize
void dz_merge_fold_col(dz_work_t *w, dz_swgv_t *col)
{
	/* squash max if scan */
	if(dz_is_scan(&w->fetcher)) { w->state.max.score = 0; }

	/* for each incoming vectors */
	for(size_t i = 0; i < w->incoming.cnt; i++) {
		dz_state_t const *ff = w->incoming.ptr[i];

		int16_t const adj  = dz_merge_calc_adj(w, ff);
		__m128i const adjv = _mm_set1_epi16(adj);
		debug("adj(%d)", adj);

		dz_swgv_t const *prev_col = dz_restore_tail_column(ff);
		for(uint64_t p = ff->range.spos; p < ff->range.epos; p++) {
			/* adjust offset */
			dz_swgv_t const v = dz_load_swgv(&prev_col[p]);
			dz_swgv_t const x = dz_load_swgv(&col[p]);
			dz_swgv_t const u = dz_max_swgv(x, dz_subs_swgv(v, adjv));
			print_vector(u.s);
			dz_store_swgv(&col[p], u);
		}
	}
	return;
}

static __dz_vectorize
dz_swgv_t *dz_merge_column(dz_work_t *w, dz_cap_t *prev_cap)
{
	/* prev_cap is head cap */
	dz_swgv_t *col = dz_slice_column(prev_cap, w->state.range.spos);

	dz_merge_init_col(w, col);
	dz_merge_fold_col(w, col);
	return(col);
}


/* fill in */
typedef struct {
	__m128i f, e, s;
	__m128i ps;			/* prev s */

	/* max tracker */
	__m128i maxv, xtv;

	dz_swgv_t *prev_col, *col;
} dz_fill_work_t;

static __dz_vectorize
void dz_fill_work_init(dz_work_t *w, dz_fill_work_t *fw)
{
	fw->f    = w->minv;
	fw->e    = w->minv;
	fw->s    = w->minv;
	fw->ps   = dz_fetch_init_vector(&w->fetcher);
	fw->maxv = w->minv;
	fw->xtv  = _mm_set1_epi16(dz_add_ofs(w->state.max.inc - w->xt));	/* next offset == current max thus X-drop threshold is always -xt */
	return;
}

static __dz_vectorize
void dz_fill_load_vector(dz_fill_work_t *fw, size_t p)
{
	fw->e = _mm_load_si128((__m128i const *)&fw->prev_col[p].e);
	fw->s = _mm_load_si128((__m128i const *)&fw->prev_col[p].s);

	debug("col(%p, %p), p(%zu)", fw->prev_col, fw->col, p);
	// print_vector(fw->e); print_vector(fw->s);
	return;
}

static __dz_vectorize
void dz_fill_update_vector(dz_work_t *w, dz_fill_work_t *fw, size_t p)
{
	__m128i const sc = dz_calc_score_profile(&w->fetcher, p);

	/* E[i, j] = max{ E[i - 1, j], S[i - 1, j] - Gi } - Ge */
	__m128i const tte = _mm_max_epu16(fw->e, _mm_subs_epu16(fw->s, w->div));
	__m128i const te = _mm_subs_epu16(tte, w->dev);
	/* print_vector(_mm_alignr_epi8(s, ps, 14)); print_vector(sc); */

	/* U[i, j] = max{ E[i, j], S[i - 1, j - 1] + sc(i, j) } */
	__m128i const tts = _mm_adds_epu16(sc, _mm_alignr_epi8(fw->s, fw->ps, 14));
	__m128i const ts = _mm_max_epu16(te, _mm_subs_epu16(tts, w->ofsv));
	fw->ps = fw->s;

	/* fold F[i, j] = max{ U[i, j] - Gi, F[i, j - 1] - Ge } */
	__m128i tf = _mm_max_epu16(
		_mm_subs_epu16(ts, w->iiv),
		_mm_subs_epu16(_mm_srli_si128(fw->f, 14), w->iev1)
	);
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 2), w->iev1));
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 4), w->iev2));
	tf = _mm_max_epu16(tf, _mm_subs_epu16(_mm_slli_si128(tf, 8), w->iev4));

	/* S[i, j] = max{ U[i, j], F[i, j] } */
	__m128i const us = _mm_max_epu16(ts, tf);

	/* update max */
	__m128i const bs = dz_query_add_bonus(w->fetcher.query, us, p);
	fw->maxv = _mm_max_epu16(fw->maxv, bs);
	print_vector(us); print_vector(bs); print_vector(fw->maxv);

	/* done */
	fw->f = tf;
	fw->e = te;
	fw->s = us;
	return;
}

static __dz_vectorize
void dz_fill_store_vector(dz_fill_work_t *fw, size_t p)
{
	_mm_store_si128((__m128i *)&fw->col[p].e, fw->e);
	_mm_store_si128((__m128i *)&fw->col[p].f, fw->f);
	_mm_store_si128((__m128i *)&fw->col[p].s, fw->s);

	// debug("col(%p), p(%zu)", fw->col, p);
	// print_vector(fw->e);
	// print_vector(fw->f);
	// print_vector(fw->s);
	return;
}

static __dz_vectorize
int16_t dz_fill_fold_max(__m128i v)
{
	__m128i t = _mm_max_epu16(v, _mm_srli_si128(v, 8));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 2));
	return(dz_rm_ofs(_mm_extract_epi16(t, 0)));
}

static __dz_vectorize
uint64_t dz_fill_test_xdrop(dz_fill_work_t *fw)
{
	__m128i xtest = _mm_subs_epu16(fw->s, fw->xtv);

	/* print_vector(_s); print_vector(_xtv); */
	return(dz_is_all_zero(xtest));
}

static __dz_vectorize
uint64_t dz_fill_is_bottom(dz_work_t *w)
{
	return(w->state.range.epos >= w->fetcher.query->blen);
}

static __dz_vectorize
uint64_t dz_fill_column_body(dz_work_t *w, dz_fill_work_t *fw)
{
	for(uint64_t p = w->state.range.spos; p < w->state.range.epos; p++) {
		/* load previous DP matrix vectors */
		dz_fill_load_vector(fw, p);

		/* load score profile and update vectors */
		dz_fill_update_vector(w, fw, p);

		/* save if at least one cell survived */
		if(dz_likely(!dz_fill_test_xdrop(fw))) {
			dz_fill_store_vector(fw, p);
			continue;
		}

		/* vector removed; clip head */
		if(p == w->state.range.spos) {
			w->state.range.spos++;
			fw->col--;
			continue;
		}

		/* terminate */
		w->state.range.epos = p;
		return(1);
	}
	return(dz_fill_is_bottom(w));
}

static __dz_vectorize
void dz_fill_column_tail(dz_work_t *w, dz_fill_work_t *fw)
{
	/* forefront extension; clip the column length if too long */
	fw->e = w->minv;
	fw->s = w->minv;

	/* update the last vector */
	dz_fill_update_vector(w, fw, w->state.range.epos);

	/* until all the cells drop */
	do {
		if(dz_fill_test_xdrop(fw)) { break; }
		dz_fill_store_vector(fw, w->state.range.epos);

		fw->f = _mm_subs_epu16(fw->f, w->iev8);
		fw->s = _mm_subs_epu16(fw->s, w->iev8);
	} while(++w->state.range.epos < w->fetcher.query->blen);
	return;
}

static __dz_vectorize
void dz_fill_update_max(dz_work_t *w, dz_fill_work_t *fw)
{
	/* update max; the actual score (absolute score accumulated from the origin) is expressed as max + inc; 2 x cmov */
	int32_t inc = dz_fill_fold_max(fw->maxv);		/* without offset */
	if(dz_cmp_max(inc, w->state.max.inc)) {
		w->state.max.inc = inc;
		w->state.max.cap = dz_cap_column(fw->col, w->state.range.epos);
	}
	debug("inc(%d), cap(%p, %p)", inc, dz_cap_column(fw->col, w->state.range.epos), w->state.max.cap);
	return;
}

static __dz_vectorize
dz_swgv_t *dz_fill_column(dz_work_t *w, dz_swgv_t *prev_col)
{
	dz_fill_work_t fw;

	/* slice new cap */
	fw.prev_col = prev_col;
	fw.col = dz_slice_cap(w, prev_col);

	/* init vectors */
	dz_fill_work_init(w, &fw);

	/* if reached the forefront of the query sequence, finish the extension */
	if(dz_fill_column_body(w, &fw)) {
		dz_fill_column_tail(w, &fw);
	}

	/* create cap object that contains [spos, epos) range (for use in the traceback routine) */
	dz_fill_update_max(w, &fw);
	return(fw.col);	/* FIXME: rescue overflow */
}

static __dz_vectorize
uint64_t dz_is_end(dz_work_t *w)
{
	return((w->state.range.epos - w->state.range.spos) == 0);
}

static __dz_vectorize
dz_state_t const *dz_extend_core(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	dz_query_t const *query,
	dz_ref_t const *ref,
	dz_state_t const **ff,
	size_t fcnt)
{
	if(fcnt == 0) { return(NULL); }		/* invalid */

	/* no extension for len == 0 */
	if(ref->len == 0 && fcnt == 1) {
		return(ff[0]);					/* no need to merge incoming vectors */
	}

	/* init working buffer */
	dz_work_t w __attribute__(( aligned(16) ));
	w.mem = mem;
	w.fetcher = dz_init_fetcher(profile, ref, query);
	dz_load_state(&w, ff, fcnt);		/* iterate over the incoming edge objects to get the current max and range */

	/* load constants */
	dz_load_profile(&w, profile);

	/* save head info */
	dz_cap_t *cap = dz_slice_head(&w);	/* range obtained in dz_load_state */
	dz_save_incoming(cap, ff, fcnt);

	/* merge incoming vectors */
	dz_swgv_t *col = dz_merge_column(&w, cap);
	debug("merged, cap(%p), col(%p)", cap, col);

	/* until X-drop */
	while(!dz_is_end(&w)) {
		if(dz_fetch_next(&w.fetcher)) { break; }

		/* forward one */
		col = dz_fill_column(&w, col);
	}

	/* calc #columns filled */
	size_t const cols = dz_finalize_fetcher(&w.fetcher);

	/* update max and range */
	dz_finalize_state(&w.state, cols);

	/* save them in tail object */
	return(dz_slice_tail(&w, col));
}



/* profile constructor */

typedef struct dz_score_conf_s {
	int8_t const *score_matrix;		/* dimension is determined by macro at compile time */

	/* gap penalties */
	uint16_t ins_open;				/* gap penalties in positive */
	uint16_t ins_extend;
	uint16_t del_open;				/* gap penalties in positive */
	uint16_t del_extend;
	uint16_t full_length_bonus;		/* end-to-end mapping bonus; only activated when compiled with -DDZ_FULL_LENGTH_BONUS */

	/* X-drop threshold */
	uint64_t max_gap_len;			/* as X-drop threshold */
} dz_score_conf_t;

/* external (custom) allocator */
typedef void *(*dz_malloc_t)(void *ctx, size_t size);
typedef struct dz_allocator_s {
	void *ctx;
	dz_malloc_t fp;
} dz_allocator_t;


#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
static __dz_vectorize
dz_profile_t *dz_alloc_profile(dz_allocator_t *alloc)
{
	dz_profile_t *profile = (dz_profile_t *)alloc->fp(alloc->ctx,
		  sizeof(dz_profile_t)
		+ 2 * DZ_MAT_SIZE * DZ_MAT_SIZE
		+ sizeof(__m128i)
	);
	return(profile);
}

static __dz_vectorize
void dz_init_score_matrix(
	dz_profile_t *profile,
	int8_t const *score_matrix)
{
	/* constants */
	__m128i const tmat = _mm_loadu_si128((__m128i const *)score_matrix);
	__m128i const ofs = _mm_set1_epi8(DZ_SCORE_OFS);
	_mm_store_si128((__m128i *)&profile->matrix[0],
		_mm_add_epi8(tmat, ofs)
	);

	/* invalid bases */
	_mm_store_si128((__m128i *)&profile->matrix[16],
		_mm_setzero_si128()
	);
	return;
}

#elif defined(DZ_NUCL_4BIT) || defined(DZ_PROTEIN)
static __dz_vectorize
dz_profile_t *dz_alloc_profile(dz_allocator_t *alloc)
{
	dz_profile_t *self = (dz_profile_t *)alloc->fp(alloc->ctx,
		  sizeof(dz_profile_t)
		+ DZ_MAT_SIZE * DZ_MAT_SIZE			/* 16 x 16 */
		+ 2 * sizeof(__m128i)
	);
	return(self);
}

static __dz_vectorize
void dz_init_score_matrix(
	dz_profile_t *profile,
	int8_t const *score_matrix)
{
	/* transpose */
	for(size_t i = 0; i < DZ_MAT_SIZE; i++) {
		for(size_t j = 0; j < DZ_MAT_SIZE; j++) {
			profile->matrix[j * DZ_MAT_SIZE + i] = score_matrix[i * DZ_MAT_SIZE + j] + DZ_SCORE_OFS;
		}
	}
	return;
}

#endif


static __dz_vectorize
void dz_init_gap_penalties(
	dz_profile_t *profile,
	dz_score_conf_t const *conf)
{
	/* insertion penalties */
	__m128i const iiv = _mm_set1_epi16(conf->ins_open);
	__m128i const iev = _mm_set1_epi16(conf->ins_extend);
	_mm_store_si128((__m128i *)profile->iiv, iiv);
	_mm_store_si128((__m128i *)profile->iev, iev);

	/* deletion penalties */
	__m128i const div = _mm_set1_epi16(conf->del_open);
	__m128i const dev = _mm_set1_epi16(conf->del_extend);
	_mm_store_si128((__m128i *)profile->div, div);
	_mm_store_si128((__m128i *)profile->dev, dev);

	/* X-drop threshold */
	uint16_t const gi = dz_max2(conf->ins_open,   conf->del_open);
	uint16_t const ge = dz_max2(conf->ins_extend, conf->del_extend);
	profile->xt = gi + ge * conf->max_gap_len;	/* X-drop threshold */
	profile->bonus = conf->full_length_bonus;
	profile->max_gap_len = conf->max_gap_len;	/* save raw value */
	debug("gi(%u), ge(%u), xdrop_threshold(%u), full_length_bonus(%u), max_gap_len(%u)", gi, ge, profile->xt, profile->bonus, profile->max_gap_len);

	return;
}

static __dz_vectorize
dz_head_t *dz_slice_root_head(dz_allocator_t *alloc, size_t blen)
{
	/* malloc the first column */
	size_t size = dz_calc_column_size(0, blen + 1);
	return(dz_head(alloc->fp(alloc->ctx, size)));
}

static __dz_vectorize
size_t dz_init_root_column(dz_profile_t const *profile, dz_score_conf_t const *conf, dz_swgv_t *col, size_t blen)
{
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(DZ_CELL_MIN));

	/* fill the root (the leftmost) column; first init vectors */
	uint16_t const ins_open = conf->ins_open;
	uint16_t const ins_extend = conf->ins_extend;

	__m128i const hv = _mm_setr_epi16(
		dz_add_ofs(0),
		dz_add_ofs(-(ins_open +     ins_extend)),
		dz_add_ofs(-(ins_open + 2 * ins_extend)),
		dz_add_ofs(-(ins_open + 3 * ins_extend)),
		dz_add_ofs(-(ins_open + 4 * ins_extend)),
		dz_add_ofs(-(ins_open + 5 * ins_extend)),
		dz_add_ofs(-(ins_open + 6 * ins_extend)),
		dz_add_ofs(-(ins_open + 7 * ins_extend))
	);
	__m128i const bv = _mm_setr_epi16(
		dz_add_ofs(-ins_open),
		dz_add_ofs(-(ins_open +     ins_extend)),
		dz_add_ofs(-(ins_open + 2 * ins_extend)),
		dz_add_ofs(-(ins_open + 3 * ins_extend)),
		dz_add_ofs(-(ins_open + 4 * ins_extend)),
		dz_add_ofs(-(ins_open + 5 * ins_extend)),
		dz_add_ofs(-(ins_open + 6 * ins_extend)),
		dz_add_ofs(-(ins_open + 7 * ins_extend))
	);
	__m128i const iev8 = _mm_slli_epi16(_mm_load_si128((__m128i const *)profile->iev), 3);

	dz_fill_work_t fw;
	fw.f = minv;
	fw.e = minv;
	fw.s = hv;
	fw.ps   = _mm_setzero_si128();
	fw.maxv = _mm_setzero_si128();
	fw.xtv  = _mm_set1_epi16(dz_add_ofs(-profile->xt));
	fw.col  = col;

	/* until the X-drop test fails on all the cells in a vector */
	for(size_t p = 0; p < blen; p++) {
		fw.f = fw.s;
		print_vector(fw.s);

		/* terminate when X-drop condition hold on all the cells */
		if(dz_fill_test_xdrop(&fw)) {
			return(p);
		}
		dz_fill_store_vector(&fw, p);

		/* reset s vector */
		if(p == 0) { fw.s = bv; }
		fw.s = _mm_subs_epu16(fw.s, iev8);		/* iev[i] * 8 */
	}
	return(blen);
}

static __dz_vectorize
dz_state_t const *dz_init_root_cap(dz_tail_t *tail, size_t epos)
{
	/* ref_save_t: dummy */
	dz_ref_save_t rsave;
	rsave.id = 0;
	rsave.len = 0;
	tail->rsave = rsave;

	/* initialize state */
	dz_state_t state;
	state.range.spos  = 0;
	state.range.epos  = epos;
	state.cnt.column  = 0;
	state.cnt.section = 0;
	state.max.score   = 0;
	state.max.inc     = 0;
	state.max.cap     = NULL;
	tail->state = state;

	/* dummy */
	tail->query = NULL;
	return(&tail->state);
}

static __dz_vectorize
void dz_init_root(
	dz_profile_t *profile,
	dz_score_conf_t const *conf,
	dz_allocator_t *alloc)
{
	/* calc vector length; query = NULL for the first (root) column */
	size_t const max_gap_len = dz_roundup(conf->max_gap_len, DZ_L);
	size_t const blen = max_gap_len / DZ_L;

	/* fill root head */
	dz_head_t *head = dz_slice_root_head(alloc, blen);
	dz_swgv_t *col = dz_init_root_head(head);

	/* DP matrix */
	size_t epos = dz_init_root_column(profile, conf, col, blen);

	/* tail object */
	dz_tail_t *tail = dz_tail(dz_cap_column(col, epos));
	dz_state_t const *root = dz_init_root_cap(tail, epos);
	debug("col(%p), range(%u, %u), tail(%p, %p)", col, root->range.spos, root->range.epos, tail, root);

	profile->root = dz_cff(root);
	return;
}

static __dz_vectorize
dz_profile_t *dz_init_profile(
	dz_allocator_t *alloc,
	dz_score_conf_t const *conf)
{
	/* allocate score priflie object and fill the root column */
	dz_profile_t *profile = dz_alloc_profile(alloc);
	dz_init_score_matrix(profile, conf->score_matrix);
	dz_init_gap_penalties(profile, conf);

	/* allocate first column and fill cells */
	dz_init_root(profile, conf, alloc);
	debug("profile(%p), root(%p)", profile, profile->root);
	return(profile);
}


/* dz_init */

/* dz_t: wrapper context */
typedef struct dz_s {
	dz_arena_t *mem;
	dz_profile_t *profile;
} dz_t;


static __dz_vectorize
dz_t *dz_initx(
	int8_t const *score_matrix,
	uint16_t ins_open,
	uint16_t ins_extend,
	uint16_t del_open,
	uint16_t del_extend,
	uint64_t max_gap_len,
	uint16_t full_length_bonus)
{
	dz_arena_t *mem = dz_arena_init(DZ_MEM_INIT_SIZE);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* allocate wrapper context */
	dz_t *self = dz_arena_malloc(mem, sizeof(dz_t));
	self->mem = mem;

	/* init profile */
	dz_allocator_t alloc;
	alloc.ctx = (void *)mem;	/* ctx */
	alloc.fp  = (dz_malloc_t)dz_arena_malloc;	/* malloc */

	dz_score_conf_t conf;
	conf.score_matrix = score_matrix;
	conf.ins_open     = ins_open;
	conf.ins_extend   = ins_extend;
	conf.del_open     = del_open;
	conf.del_extend   = del_extend;
	conf.full_length_bonus = full_length_bonus;
	conf.max_gap_len  = max_gap_len;

	self->profile = dz_init_profile(&alloc, &conf);
	return(self);
}

#ifdef DZ_FULL_LENGTH_BONUS
static __dz_vectorize
dz_t *dz_init(
	int8_t const *score_matrix,
	uint16_t gap_open,
	uint16_t gap_extend,
	uint64_t max_gap_len,
	uint16_t full_length_bonus)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len,
		full_length_bonus
	));
}

#else
static __dz_vectorize
dz_t *dz_init(
	int8_t const *score_matrix,
	uint16_t gap_open,
	uint16_t gap_extend,
	uint64_t max_gap_len)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len, 0
	));
}

#endif

static __dz_vectorize
void dz_destroy(dz_t *self)
{
	debug("self(%p)", self);
	if(self == NULL) { return; }
	dz_arena_destroy(self->mem);
	return;
}
static __dz_vectorize
void dz_flush(dz_t *self)
{
	dz_arena_flush(self->mem);
	// dz_arena(self)->stack.top = (uint8_t *)(self->root + 1);
	return;
}



/* extend and scan wrapper */
/**
 * @fn dz_extend, dz_scan
 * @brief dz_extend for semi-global alignment (head gaps are penalized for the both directions),
 * and dz_scan for searching anchor (head gap is only penalized for the query side).
 *
 * NOTE: DZ_N_AS_UNMATCHING_BASE is not recommended when dz_scan is used
 */
static __dz_vectorize
dz_forefront_t const *dz_extend(
	dz_t *self,
	dz_query_t const *query,
	dz_forefront_t const **forefronts, size_t fcnt,
	char const *ref, int32_t rlen, uint32_t rid)
{
	/* compose ref object */
	dz_ref_t r;
	r.ptr = ref;
	r.len = rlen < 0 ? -rlen : rlen;	/* always positive */
	r.id  = rid;
	r.dir = rlen < 0 ? -1 : 1;			/* -1 for reverse */
	r.init_s = INT16_MIN;

	debug("extend, self(%p), fcnt(%zu), query(%p), ref(%p, %s), rlen(%d), rid(%u)", self, fcnt, query, ref, ref, rlen, rid);
	return(dz_cff(dz_extend_core(self->mem, self->profile, query, &r, dz_cpstate(forefronts), fcnt)));
}
static __dz_vectorize
dz_forefront_t const *dz_scan(
	dz_t *self,
	dz_query_t const *query,
	dz_forefront_t const **forefronts, size_t fcnt,
	char const *ref, int32_t rlen, uint32_t rid)
{
	dz_ref_t r;
	r.ptr = ref;
	r.len = rlen < 0 ? -rlen : rlen;
	r.id  = rid;
	r.dir = rlen < 0 ? -1 : 1;
	r.init_s = 0;

	debug("scan, self(%p), fcnt(%zu), query(%p), ref(%p, %s), rlen(%d), rid(%u)", self, fcnt, query, ref, ref, rlen, rid);
	return(dz_cff(dz_extend_core(self->mem, self->profile, query, &r, dz_cpstate(forefronts), fcnt)));
}



#if defined(UNITTEST) && UNITTEST != 0

/* unittest stuffs */
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

#elif defined(DZ_NUCL_4BIT)
static size_t const dz_unittest_query_length = 560;
static char const *dz_unittest_query =
	"\x1\x4\x2\x8\x8\x8\x8\x2\x1\x8\x8\x2\x8\x4\x1\x2\x8\x4\x2\x1\x1\x2\x4\x4\x4\x2\x1\x1\x8\x1\x8\x4\x8\x2\x8"
	"\x2\x8\x4\x8\x4\x8\x4\x4\x1\x8\x8\x1\x1\x1\x1\x1\x1\x1\x4\x1\x4\x8\x4\x8\x2\x8\x4\x1\x8\x1\x4\x2\x1\x4\x2"
	"\x8\x8\x2\x8\x4\x1\x1\x2\x8\x4\x4\x8\x8\x1\x2\x2\x8\x4\x2\x2\x4\x8\x4\x1\x4\x8\x1\x1\x1\x8\x8\x1\x1\x1\x1"
	"\x8\x8\x8\x8\x1\x8\x8\x4\x1\x2\x8\x8\x1\x4\x4\x8\x2\x1\x2\x8\x1\x1\x1\x8\x1\x2\x8\x8\x8\x1\x1\x2\x2\x1\x1"
	"\x8\x1\x8\x1\x4\x4\x2\x1\x8\x1\x4\x2\x4\x2\x1\x2\x1\x4\x1\x2\x1\x4\x1\x8\x1\x1\x1\x1\x1\x8\x8\x1\x2\x1\x4"
	"\x1\x4\x8\x1\x2\x1\x2\x1\x1\x2\x1\x8\x2\x2\x1\x8\x4\x1\x1\x1\x2\x4\x2\x1\x8\x8\x1\x4\x2\x1\x2\x2\x1\x2\x2"
	"\x1\x8\x8\x1\x2\x2\x1\x2\x2\x1\x2\x2\x1\x8\x2\x1\x2\x2\x1\x8\x8\x1\x2\x2\x1\x2\x1\x4\x4\x8\x1\x1\x2\x4\x4"
	"\x8\x4\x2\x4\x4\x4\x2\x8\x4\x1\x2\x4\x2\x4\x8\x1\x2\x1\x4\x4\x1\x1\x1\x2\x1\x2\x1\x4\x1\x1\x1\x1\x1\x1\x4"
	"\x2\x2\x2\x4\x2\x1\x2\x2\x8\x4\x1\x2\x1\x4\x8\x4\x2\x4\x4\x4\x2\x8\x8\x8\x8\x8\x8\x8\x8\x8\x2\x4\x1\x2\x2"
	"\x1\x1\x1\x4\x4\x8\x1\x1\x2\x4\x1\x4\x4\x8\x1\x1\x2\x1\x1\x2\x2\x1\x8\x4\x2\x4\x1\x4\x8\x4\x8\x8\x4\x1\x1"
	"\x4\x8\x8\x2\x4\x4\x2\x4\x4\x8\x1\x2\x1\x8\x2\x1\x4\x8\x4\x4\x2\x1\x1\x1\x8\x4\x2\x1\x4\x1\x1\x2\x4\x8\x8"
	"\x8\x8\x2\x8\x4\x2\x4\x8\x4\x8\x8\x4\x2\x2\x4\x1\x8\x1\x8\x8\x2\x8\x4\x4\x1\x1\x1\x4\x2\x1\x1\x8\x4\x2\x2"
	"\x1\x4\x4\x2\x1\x4\x4\x4\x4\x2\x1\x4\x4\x8\x4\x4\x2\x2\x1\x2\x2\x4\x8\x2\x2\x8\x2\x8\x2\x8\x4\x2\x2\x2\x2"
	"\x2\x4\x2\x2\x1\x1\x1\x1\x8\x2\x1\x2\x2\x1\x1\x2\x2\x1\x2\x2\x8\x4\x4\x8\x4\x4\x2\x4\x1\x8\x4\x1\x8\x8\x4"
	"\x1\x1\x1\x1\x1\x1\x2\x2\x1\x8\x8\x1\x4\x2\x4\x4\x2\x2\x1\x4\x4\x1\x8\x4\x2\x8\x8\x8\x1\x2\x2\x2\x1\x1\x8"
	"\x1\x8\x2\x1\x4\x2\x4\x1\x8\x4\x2\x2\x4\x1\x1\x2\x4\x8\x1\x8\x8\x8\x8\x8\x4\x2\x2\x4\x1\x1\x2\x8\x8\x8\x8";

static int8_t const dz_unittest_score_matrix[16 * 16] = {
	2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2, -3,
	-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 2
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
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);
	ut_assert(dz_root(dz) != NULL);
	dz_destroy(dz);
}

#endif	/* defined(UNITTEST) && UNITTEST != 0 */



/* query profile calculation */
typedef __m128i (*dz_query_conv_bulk_t)(int8_t const *score_matrix, __m128i v);
typedef uint8_t (*dz_query_conv_single_t)(int8_t const *score_matrix, uint8_t c);

typedef struct {
	dz_query_conv_bulk_t bulk;
	dz_query_conv_single_t single;
} dz_query_conv_t;

/**
 * @fn dz_pack_query, dz_pack_query_reverse
 * packed sequence object is allocated inside the context so freed by dz_flush
 */
#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT) || defined(DZ_NUCL_4BIT)

/* profile vector is calculated on-the-fly */
static __dz_vectorize
dz_query_t *dz_pack_query_alloc_mem(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	char const *query,
	size_t qlen)
{
	/* allocate mem */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	dz_query_t *q = (dz_query_t *)dz_arena_malloc(mem,
		  sizeof(dz_query_t)		/* header */
		+ sizeof(__m128i)			/* score matrix (4 x 4) */
		+ clen						/* body (sequence) */
		+ sizeof(__m128i)			/* margin at the tail */
	);

	*q = (dz_query_t){
		.blen = qlen == 0 ? 0 : (clen / DZ_L),
		.q = query,
		.bonus = { 0 }
	};

	/*
	 * tentative support of full-length bonus;
	 * precaclulated bonus score is saved at q->bonus[0] and q->bonus[1];
	 * [0] for non-tail vector and [1] for tail vector
	 */
	q->bonus[DZ_L + (qlen % DZ_L)] = profile->bonus;
	return(q);
}

#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)

/* { ASCII / 2bit } -> 2bit conversion */
static __dz_vectorize
__m128i dz_query_conv_bulk_forward(int8_t const *score_matrix, __m128i v)
{
	dz_unused(score_matrix);

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
	__m128i const cv = _mm_load_si128((__m128i const *)conv);	/* conversion table */
	return(_mm_shuffle_epi8(cv, v));
}

static __dz_vectorize
uint8_t dz_query_conv_single_forward(int8_t const *score_matrix, uint8_t c)
{
	dz_unused(score_matrix);

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		#ifdef DZ_NUCL_ASCII
			/* ['A' & 0x0f] = A, ['C' & 0x0f] = C, ['G' & 0x0f] = G, ['T' & 0x0f] = T, ['U' & 0x0f] = T, ['N' & 0x0f] = N */
			0, qA, 0, qC, qT, qU, 0, qG, 0, 0, 0, 0, 0, 0, qN, 0
		#else /* DZ_NUCL_2BIT */
			qA, qC, qG, qT, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
		#endif
	};
	return(conv[c & 0x0f]);
}


static __dz_vectorize
__m128i dz_query_conv_bulk_reverse(int8_t const *score_matrix, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		#ifdef DZ_NUCL_ASCII
			0, qT, 0, qG, qA, qA, 0, qC, 0, 0, 0, 0, 0, 0, qN, 0
		#else /* DZ_NUCL_2BIT */
			qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
		#endif
	};
	__m128i const cv = _mm_load_si128((__m128i const *)conv);	/* conversion table */
	return(_mm_shuffle_epi8(cv, v));
}

static __dz_vectorize
uint8_t dz_query_conv_single_reverse(int8_t const *score_matrix, uint8_t c)
{
	dz_unused(score_matrix);

	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		#ifdef DZ_NUCL_ASCII
			0, qT, 0, qG, qA, qA, 0, qC, 0, 0, 0, 0, 0, 0, qN, 0
		#else /* DZ_NUCL_2BIT */
			qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
		#endif
	};
	return(conv[c & 0x0f]);
}


#elif defined(DZ_NUCL_4BIT)
static __dz_vectorize
__m128i dz_query_conv_bulk_forward(int8_t const *score_matrix, __m128i v)
{
	dz_unused(score_matrix);

	// print_vector_8(v);
	return(v);
}

static __dz_vectorize
uint8_t dz_query_conv_single_forward(int8_t const *score_matrix, uint8_t c)
{
	dz_unused(score_matrix);

	return(c & 0x0f);
}

/* complement for reverse */
static __dz_vectorize
__m128i dz_query_conv_bulk_reverse(int8_t const *score_matrix, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	__m128i const cv = _mm_load_si128((__m128i const *)comp);
	__m128i const w = _mm_shuffle_epi8(cv, v);

	// print_vector_8(v);
	// print_vector_8(cv);
	// print_vector_8(w);
	return(w);
}

static __dz_vectorize
uint8_t dz_query_conv_single_reverse(int8_t const *score_matrix, uint8_t c)
{
	dz_unused(score_matrix);

	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};

	debug("c(%x, %x)", c, comp[c & 0x0f]);
	return(comp[c & 0x0f]);
}

#endif


static __dz_vectorize
size_t dz_pack_query_forward_core(
	dz_query_t *q,
	dz_profile_t const *profile,
	dz_query_conv_t qconv,
	char const *query,
	size_t qlen)
{
	dz_unused(profile);		/* unused for nucleotide */

	size_t const clen = dz_roundup(qlen + 1, DZ_L);		/* FIXME: dz_roundup(qlen + 1, sizeof(__m128i)) ? */
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));

	/* until the end of the query sequence */
	__m128i pv = _mm_set1_epi8((int8_t)qN);			/* initial vector with N's */
	for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);

		/* shuffle: ASCII -> 2bit conversion */
		__m128i const tv = qconv.bulk(NULL, qv);

		/* shift by one to make room for the top row */
		_mm_store_si128((__m128i *)&q->arr[i],
			_mm_alignr_epi8(tv, pv, 15)
		);

		/* save for the next iteration */
		pv = tv;
	}

	/* continue the same conversion on the remainder */
	q->arr[blen] = _mm_extract_epi8(pv, 15);
	for(size_t i = blen; i < qlen; i++) {
		q->arr[i + 1] = (int8_t)qconv.single(NULL, query[i]);
	}

	/* fill tail margin */
	for(size_t i = qlen; i < clen; i++) {
		q->arr[i + 1] = qX;
	}

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(clen);
}
static __dz_vectorize
size_t dz_pack_query_reverse_core(
	dz_query_t *q,
	dz_profile_t const *profile,
	dz_query_conv_t qconv,
	char const *query,
	size_t qlen)
{
	dz_unused(profile);

	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));

	/* for reversing vector (bswapdq) */
	static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
		15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	};
	__m128i const rv = _mm_load_si128((__m128i const *)rev);

	/* until the end of the query sequence */
	__m128i pv = _mm_set1_epi8((int8_t)qN);
	for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
		// print_vector_8(qv);
		__m128i const tv = qconv.bulk(NULL, _mm_shuffle_epi8(qv, rv));

		/* shift by one to make room for the top row */
		_mm_store_si128((__m128i *)&q->arr[i],
			_mm_alignr_epi8(tv, pv, 15)
		);

		/* save for the next iteration */
		pv = tv;
	}

	/* continue the same conversion on the remainings */
	q->arr[blen] = _mm_extract_epi8(pv, 15);
	for(size_t i = blen; i < qlen; i++) {
		q->arr[i + 1] = (int8_t)qconv.single(NULL, query[qlen - 1 - i]);
	}

	/* fill tail margin */
	for(size_t i = qlen; i < clen; i++) {
		q->arr[i + 1] = qX;
	}

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(clen);
}

#elif defined(DZ_PROTEIN)

/* allocate precalculated profile vector */
static __dz_vectorize
dz_query_t *dz_pack_query_alloc_mem(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	char const *query,
	size_t qlen)
{
	/* allocate mem */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	dz_query_t *q = (dz_query_t *)dz_arena_malloc(mem,
		  sizeof(dz_query_t)		/* header */
		+ DZ_MAT_SIZE * clen		/* profile vectors */
		+ sizeof(__m128i)			/* margin at the tail */
	);

	*q = (dz_query_t){
		.blen = qlen == 0 ? 0 : (clen / DZ_L),
		.q = query,
		.bonus = { 0 }
	};
	q->bonus[DZ_L + (qlen % DZ_L)] = profile->bonus;
	return(q);
}

static __dz_vectorize
__m128i dz_query_conv_bulk_forward(int8_t const *score_matrix, __m128i v)
{
	__m128i const lcv = _mm_loadu_si128((__m128i const *)&score_matrix[0]);		/* lower conversion table */
	__m128i const hcv = _mm_loadu_si128((__m128i const *)&score_matrix[16]);	/* upper conversion table */

	/* apply lower and upper conversion mask */
	__m128i const lv = _mm_shuffle_epi8(lcv, v);
	__m128i const hv = _mm_shuffle_epi8(hcv, v);

	/* blend lower and upper with bit [4] */
	__m128i const tv = _mm_blendv_epi8(lv, hv, _mm_slli_epi32(v, 3));
	return(tv);
}

static __dz_vectorize
uint8_t dz_query_conv_single_forward(int8_t const *score_matrix, uint8_t c)
{
	/* conv[(uint8_t)query[i] & 0x1f]; */
	return(score_matrix[c & 0x1f]);
}

/* we don't discriminate forward and reverse for protein sequence */
static __dz_vectorize
__m128i dz_query_conv_bulk_reverse(int8_t const *score_matrix, __m128i v)
{
	return(dz_query_conv_bulk_forward(score_matrix, v));
}

static __dz_vectorize
uint8_t dz_query_conv_single_reverse(int8_t const *score_matrix, uint8_t c)
{
	return(dz_query_conv_single_forward(score_matrix, c));
}

static __dz_vectorize
size_t dz_pack_query_forward_core(
	dz_query_t *q,
	dz_profile_t const *profile,
	dz_query_conv_t qconv,
	char const *query,
	size_t qlen)
{
	/* constants */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));

	/* calculate profile for each base (residue) */
	for(size_t j = 0; j < DZ_MAT_SIZE; j++) {
		int8_t *a = (int8_t *)&q->arr[j * clen];
		a[0] = 0;				/* as -DZ_SCORE_OFS */

		int8_t const *score_matrix = &profile->matrix[j * DZ_MAT_SIZE];
		__m128i pv = _mm_setzero_si128();	/* as -DZ_SCORE_OFS */

		debug("j(%zu)", j);
		// print_vector_8(_mm_loadu_si128((__m128i const *)score_matrix));

		/* until the end of the query sequence */
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
			__m128i const tv = qconv.bulk(score_matrix, qv);

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i],
				_mm_alignr_epi8(tv, pv, 15)
			);

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */
		a[blen] = _mm_extract_epi8(pv, 15);
		for(size_t i = blen; i < qlen; i++) {
			a[i + 1] = (int8_t)qconv.single(score_matrix, query[i]);
		}
		for(size_t i = qlen; i < clen; i++) {
			a[i + 1] = 0;		/* as -DZ_SCORE_OFS */
		}
	}
	debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}
static __dz_vectorize
size_t dz_pack_query_reverse_core(
	dz_query_t *q,
	dz_profile_t const *profile,
	dz_query_conv_t qconv,
	char const *query,
	size_t qlen)
{
	/* constants */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));

	static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
		15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	};
	__m128i const rv = _mm_load_si128((__m128i const *)rev);

	for(size_t j = 0; j < DZ_MAT_SIZE; j++) {
		int8_t *a = (int8_t *)&q->arr[j * clen];
		a[0] = 0;

		int8_t const *score_matrix = &profile->matrix[j * DZ_MAT_SIZE];
		__m128i pv = _mm_setzero_si128();

		/* toward the head of the query sequence */
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
			__m128i const tv = qconv.bulk(score_matrix, _mm_shuffle_epi8(qv, rv));		/* bswapdq before apply conversion */

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i],
				_mm_alignr_epi8(tv, pv, 15)
			);

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */
		a[blen] = _mm_extract_epi8(pv, 15);
		for(size_t i = blen; i < qlen; i++) {
			a[i + 1] = (int8_t)qconv.single(score_matrix, query[qlen - 1 - i]);
		}
		for(size_t i = qlen; i < dz_roundup(qlen + 1, DZ_L); i++) {
			a[i + 1] = 0;		/* as -DZ_SCORE_OFS */
		}
	}

	debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}
#endif


static __dz_vectorize
dz_query_t *dz_pack_query_core(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	char const *query,
	int64_t qlen)
{
	dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, query, qlen);

	if(qlen >= 0) {
		dz_query_conv_t qconv = {
			.single = dz_query_conv_single_forward,
			.bulk   = dz_query_conv_bulk_forward
		};
		dz_pack_query_forward_core(q, profile, qconv, query, (size_t)qlen);
	} else {
		dz_query_conv_t qconv = {
			.single = dz_query_conv_single_reverse,
			.bulk   = dz_query_conv_bulk_reverse
		};
		/* qlen is negative */
		dz_pack_query_reverse_core(q, profile, qconv, &query[qlen], (size_t)-qlen);
	}
	return(q);
}


/* wrappers for compatibility */
static __dz_vectorize
dz_query_t *dz_pack_query(dz_t *self, char const *query, int64_t qlen)
{
	return(dz_pack_query_core(self->mem, self->profile, query, qlen));
}

static __dz_vectorize
dz_query_t *dz_pack_query_reverse(dz_t *self, char const *query, size_t qlen)
{
	dz_query_t *q = dz_pack_query_alloc_mem(self->mem, self->profile, query, qlen);
	dz_query_conv_t qconv = {
		.single = dz_query_conv_single_reverse,
		.bulk   = dz_query_conv_bulk_reverse
	};
	dz_pack_query_reverse_core(q, self->profile, qconv, query, qlen);
	return(q);
}

static __dz_vectorize
dz_query_t *dz_pack_query_forward(dz_t *self, char const *query, size_t qlen)
{
	dz_query_t *q = dz_pack_query_alloc_mem(self->mem, self->profile, query, qlen);
	dz_query_conv_t qconv = {
		.single = dz_query_conv_single_forward,
		.bulk   = dz_query_conv_bulk_forward
	};
	dz_pack_query_forward_core(q, self->profile, qconv, query, qlen);
	return(q);
}



#if defined(UNITTEST) && UNITTEST != 0
unittest() {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_unused(q);
	dz_destroy(dz);
}
#endif	/* defined(UNITTEST) && UNITTEST != 0 */


#if defined(UNITTEST) && UNITTEST != 0

/* short, exact matching sequences */
unittest( "extend.base" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 3; trial++) {
		dz_query_t const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			:              10
		);
		dz_forefront_t const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront == *dz_root(dz));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "\x1", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 2, 5), "max(%u), cell(%u)", forefront->max, dz_ut_sel(2, 2, 2, 5));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 4, 9), "max(%u), cell(%u)", forefront->max, dz_ut_sel(4, 4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "\x1\x4\x1\x8\x8\x8\x8", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28), "max(%u), cell(%u)", forefront->max, dz_ut_sel(9, 9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "\x1\x4\x1\x8\x8\x8\x8\x2", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 34), "max(%u), cell(%u)", forefront->max, dz_ut_sel(11, 11, 11, 34));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "\x1\x4\x1\x8\x8\x8\x8\x2\x1", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 39), "max(%u), cell(%u)", forefront->max, dz_ut_sel(13, 13, 13, 39));

		dz_unused(forefront);
	}
	dz_destroy(dz);
}
#ifndef DZ_PROTEIN
unittest( "extend.base.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	char revcomp[dz_unittest_query_length];
	for(size_t i = 0; i < dz_unittest_query_length; i++) {
		switch(dz_unittest_query[dz_unittest_query_length - i - 1]) {
			case 'A': revcomp[i] = 'T'; break;
			case 'C': revcomp[i] = 'G'; break;
			case 'G': revcomp[i] = 'C'; break;
			case 'T': revcomp[i] = 'A'; break;
			case rA: revcomp[i] = rT; break;
			case rC: revcomp[i] = rG; break;
			case rG: revcomp[i] = rC; break;
			case rT: revcomp[i] = rA; break;
		}
	}
	for(size_t trial = 0; trial < 3; trial++) {
		size_t length = trial == 0 ? dz_unittest_query_length
					  : trial == 1 ? 3
					  :              10;
		dz_query_t const *q = dz_pack_query_reverse(dz, &revcomp[dz_unittest_query_length - length], length);
		dz_forefront_t const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront == *dz_root(dz));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "\x1", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 2, 5));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 5);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "\x1\x4\x1\x8\x8\x8\x8", "MASLVQT"), 7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"AAAATCT"[7], &"\x0\x0\x0\x0\x3\x1\x3"[7], &"\x1\x1\x1\x1\x8\x2\x8"[7], "MASLVQT"), -7, 6);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "\x1\x4\x1\x8\x8\x8\x8\x2", "MASLVQTG"), 8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 34));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"GAAAATCT"[8], &"\x2\x0\x0\x0\x0\x3\x1\x3"[8], &"\x4\x1\x1\x1\x1\x8\x2\x8"[8], "MASLVQTG"), -8, 7);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 11, 28));

		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "\x1\x4\x1\x8\x8\x8\x8\x2\x1", "MASLVQTGK"), 9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 39));
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"TGAAAATCT"[9], &"\x3\x2\x0\x0\x0\x0\x3\x1\x3"[9], &"\x8\x4\x1\x1\x1\x1\x8\x2\x8"[9], "MASLVQTGK"), -9, 8);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 13, 28));

		dz_unused(forefront);
	}
	dz_destroy(dz);
}
#endif

/* a small graph */
unittest( "extend.small" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 4; trial++) {
		dz_query_t const *q = dz_pack_query(dz, dz_unittest_query,
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			: trial == 2 ? 4
			:              7
		);
		dz_forefront_t const *forefronts[5] = { NULL };

		/*
		 * AG---TTTT------CTGA
		 *   \ /    \    /
		 *    C      CATT
		 */
		forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
		ut_assert(forefronts[0] != NULL && forefronts[0]->max == dz_ut_sel(4, 4, 4, 9));

		forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
		ut_assert(forefronts[1] != NULL && forefronts[1]->max == dz_ut_sel(6, 6, 6, 14));
		if(trial == 1) { continue; }

		forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
		if(trial == 2) {
			ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(8, 8, 8, 18));
			continue;
		}
		ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(14, 14, 14, 32));
		if(trial == 3) { continue; }

		forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATT", "\x1\x0\x3\x3", "\x2\x1\x8\x8", "CKAK"), 4, 4);
		ut_assert(forefronts[3] != NULL && forefronts[3]->max == dz_ut_sel(22, 22, 22, 43));

		forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
		ut_assert(forefronts[4] != NULL && forefronts[4]->max == dz_ut_sel(30, 30, 30, 61));
	}
	dz_destroy(dz);
}

#endif	/* defined(UNITTEST) && UNITTEST != 0 */



/* max pos search */

static __dz_vectorize
int64_t dz_calc_max_rpos_core(dz_state_t const *ff)
{
	/* max-scoring column */
	dz_cap_t const *pcap = ff->max.cap;

	/* cap->rrem overlaps tail->rsave.rlen */
	dz_tail_t const *tail = dz_restore_tail(ff);
	int32_t rpos = (pcap == dz_ccap(tail)) ? 0 : pcap->rstate.rem;
	return((int64_t)rpos);			/* signed expansion */
}

static __dz_vectorize
uint64_t dz_finalize_qpos(size_t p, uint64_t eq)
{
	/* tzcntq is faster but avoid using it b/c it requires relatively newer archs */
	uint64_t zcnt = (eq - 1) & (~eq & 0x5555);	/* subq, andnq, andq; chain length == 2 */
	zcnt += zcnt>>2; zcnt &= 0x3333;			/* shrq, addq, andq; chain length might be shorter */
	zcnt += zcnt>>4; zcnt &= 0x0f0f;
	zcnt += zcnt>>8; zcnt &= 0x00ff;
	debug("found, eq(%zx), zcnt(%zu), idx(%zu)", (size_t)eq, (size_t)zcnt, (size_t)(p * DZ_L + zcnt));
	return(p * DZ_L + zcnt);
}

static __dz_vectorize
uint64_t dz_calc_max_qpos_core(dz_state_t const *ff)
{
	/* load max score */
	dz_cap_t const *pcap = ff->max.cap;
	__m128i const maxv = _mm_set1_epi16(dz_add_ofs(ff->max.inc));

	/* load query pointer */
	dz_query_t const *query = dz_restore_tail(ff)->query;

	/* loac column pointer */
	dz_swgv_t const *col = dz_restore_column(pcap);
	for(size_t p = pcap->range.spos; p < pcap->range.epos; p++) {
		__m128i const v = _mm_load_si128(&col[p].s);
		__m128i const s = dz_query_add_bonus(query, v, p);
		print_vector(s);

		uint64_t eq = _mm_movemask_epi8(_mm_cmpeq_epi16(s, maxv));
		if(eq == 0) { continue; }
		
		return(dz_finalize_qpos(p, eq));
	}
	return(0);
}

static __dz_vectorize
uint64_t dz_calc_max_pos_core(dz_state_t const *ff)
{
	dz_cap_t const *pcap = ff->max.cap;

	if(pcap == NULL) {
		/* is root */
		debug("pcap == NULL");
		dz_tail_t const *tail = dz_restore_tail(ff);
		return(((uint64_t)tail->rsave.len)<<32);
	}

	debug("ff(%p), pcap(%p)", ff, pcap);
	return((dz_calc_max_rpos_core(ff)<<32) | dz_calc_max_qpos_core(ff));
}


static __dz_vectorize
int64_t dz_calc_max_rpos(dz_t *self, dz_forefront_t const *ff)
{
	dz_unused(self);
	return(dz_calc_max_rpos_core(dz_cstate(ff)));
}

static __dz_vectorize
uint64_t dz_calc_max_qpos(dz_t *self, dz_forefront_t const *ff)
{
	dz_unused(self);
	return(dz_calc_max_qpos_core(dz_cstate(ff)));
}

static __dz_vectorize
uint64_t dz_calc_max_pos(dz_t *self, dz_forefront_t const *ff)
{
	dz_unused(self);
	return(dz_calc_max_pos_core(dz_cstate(ff)));
}



/* traceback */

/* alignment path object */
typedef struct dz_path_span_s {
	uint32_t id;
	uint32_t offset;
} dz_path_span_t;

/* alignment object */
typedef struct dz_alignment_s {
	dz_path_span_t const *span;
	uint8_t const *path;
	uint32_t span_length, path_length, ref_length, query_length;
	int32_t rrem, score;
	uint32_t mismatch_count, match_count, ins_count, del_count;
} dz_alignment_t;

/* working buffer */
typedef struct {
	/* dst arrays */
	dz_alignment_t *aln;
	struct {
		dz_path_span_t *ptr, *base;
	} span;
	struct {
		uint8_t *ptr, *base;
	} path;

	/* lengths */
	size_t rlen;
	dz_query_t const *query;
	dz_profile_t const *profile;

	/* match / gap counters */
	uint32_t cnt[4];

	/* cell indices */
	dz_cap_t const *pcap, *cap;
	size_t idx;
	uint16_t score, ie, de;
	uint8_t rch;
} dz_trace_work_t;

static __dz_vectorize
size_t dz_trace_vector_idx(size_t idx)
{
	return(idx / DZ_L);
}

static __dz_vectorize
size_t dz_trace_cell_idx(size_t idx)
{
	return(idx & (DZ_L - 1));
}

static __dz_vectorize
uint16_t dz_trace_score(size_t layer, dz_cap_t const *cap, size_t idx)
{
	dz_swgv_t const *col = dz_restore_column(cap);

	/* vector */
	size_t const vidx = dz_trace_vector_idx(idx);
	uint16_t const *v = (uint16_t const *)&col[vidx];

	/* cell */
	size_t const cidx = dz_trace_cell_idx(idx);
	return(v[layer * DZ_L + cidx]);
}

static __dz_vectorize
void dz_trace_push_span(dz_trace_work_t *w, size_t id)
{
	/* finzlize current section */
	w->span.ptr->offset = w->path.base - w->path.ptr;
	debug("offset(%u)", w->span.ptr->offset);

	/* open new bin */
	w->span.ptr--;
	w->span.ptr->id = id;
	return;
}

static __dz_vectorize
uint64_t dz_trace_reload_section(dz_trace_work_t *w, size_t layer)
{
	/* just load pcap if internal bridge */
	if(dz_is_ilink(w->pcap)) {
		w->pcap = dz_rewind_ilink(w->pcap);
		return(0);
	}

	/* merging vector; load contents to find an edge */
	dz_head_t const *head = dz_chead(w->pcap);
	uint16_t const prev_score = dz_trace_score(layer, w->cap, w->idx) + head->adj;
	debug("head(%p), prev_score(%u)", head, prev_score);

	/* load incoming vectors */
	size_t const vidx = dz_trace_vector_idx(w->idx);
	size_t const fcnt = head->fcnt;
	dz_state_t const **farr = dz_cpstate(w->pcap) - fcnt;
	for(size_t i = 0; i < fcnt; i++) {
		dz_state_t const *ff = farr[i];
		debug("i(%zu), ff(%p)", i, ff);

		/* adj[i] = w.max - (ffs[i]->max - ffs[i]->inc); base = max - inc */
		if(!dz_inside(ff->range.spos, vidx, ff->range.epos)) { continue; }

		dz_tail_t const *tail = dz_restore_tail(ff);
		uint16_t const s = dz_trace_score(layer, dz_ccap(tail), w->idx) + (ff->max.score - ff->max.inc);
		debug("s(%u), prev_score(%u)", s, prev_score);
		if(prev_score == s) {
			w->pcap = dz_cap(tail);
			return(0);
		}
	}

	/* something is wrong */
	dz_trap();
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_unwind_h(dz_trace_work_t *w, size_t layer)
{
	debug("pcap(%p), cap(%p)", w->pcap, w->cap);
	while(1) {
		w->cap  = w->pcap;
		w->pcap = dz_unwind_cap(w->cap);
		debug("pcap(%p), cap(%p), is_head(%u), is_root(%u)", w->pcap, w->cap, (uint32_t)dz_is_head(w->pcap), (uint32_t)dz_is_root(w->pcap));

		if(dz_likely(!dz_is_head(w->pcap))) { break; }	/* escape if not head */
		if(dz_is_root(w->pcap)) { break; }				/* escape if root */

		/* pcap is valid merging node */
		dz_trace_reload_section(w, layer);

		/* push segment info */
		dz_trace_push_span(w, dz_ctail(w->pcap)->rsave.id);
		w->score = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx);
	}

	/* return the reference-side base */
	w->rlen++;
	w->rch = w->pcap->rstate.ch;

	debug("pcap(%p, %x), cap(%p, %x)", w->pcap, w->pcap->rstate.ch, w->cap, w->cap->rstate.ch);
	return(0);
}


static __dz_vectorize
uint64_t dz_trace_test_idx(dz_trace_work_t *w, size_t ofs)
{
	size_t const vidx = dz_trace_vector_idx(w->idx - ofs);
	return(!dz_inside(w->pcap->range.spos, vidx, w->pcap->range.epos));
}

static __dz_vectorize
uint64_t dz_trace_unwind_v(dz_trace_work_t *w, size_t layer)
{
	dz_unused(layer);

	w->idx--;
	return(0);
}


static __dz_vectorize
void dz_trace_push_op(dz_trace_work_t *w, uint64_t op, uint16_t next_score)
{
	*--w->path.ptr = DZ_CIGAR_INTL>>(op<<3);
	w->cnt[op]++;
	w->score = next_score;
	debug("score(%u), op(%c)", next_score, *w->path.ptr);
	return;
}


static __dz_vectorize
uint64_t dz_trace_eat_match(dz_trace_work_t *w) {
	if(dz_trace_test_idx(w, 1)) { return(0); }

	/* skip if score does not match */
	uint16_t const s = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx - 1);
	uint16_t const p = dz_pair_score(w->profile, w->query, w->rch, w->idx);
	if(w->score != (uint16_t)(s + p)) { return(0); }
	debug("match, rch(%x)", w->rch);

	/* determine match state */
	uint64_t const eq = dz_pair_eq(w->profile, w->query, w->rch, w->idx);
	dz_trace_push_op(w, DZ_S_MATRIX + eq, s);

	/* unwind both */
	dz_trace_unwind_v(w, DZ_S_MATRIX);
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_ins(dz_trace_work_t *w) {
	if(dz_trace_test_idx(w, 1)) { return(0); }

	/* skip if score does not match */
	uint16_t const f = dz_trace_score(DZ_F_MATRIX, w->cap, w->idx);
	if(dz_likely(w->score != f)) { return(0); }
	debug("ins, score(%u), f(%u)", w->score, f);

	do {
		uint16_t const x = dz_trace_score(DZ_F_MATRIX, w->cap, w->idx - 1);
		debug("ins, score(%u), x(%u), ie(%u)", w->score, x, w->ie);
		if(w->score != x - w->ie) { break; }

		dz_trace_push_op(w, DZ_F_MATRIX, x);
		dz_trace_unwind_v(w, DZ_F_MATRIX);
	} while(!dz_trace_test_idx(w, 1));

	/* eat last column */
	dz_trace_push_op(w, DZ_F_MATRIX, dz_trace_score(DZ_S_MATRIX, w->cap, w->idx - 1));
	dz_trace_unwind_v(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_del(dz_trace_work_t *w) {
	if(dz_trace_test_idx(w, 0)) { return(0); }

	/* skip if score does not match */
	uint16_t const e = dz_trace_score(DZ_E_MATRIX, w->cap, w->idx);
	if(dz_likely(w->score != e)) { return(0); }
	debug("del, score(%u), f(%u)", w->score, e);

	do {
		uint16_t const x = dz_trace_score(DZ_E_MATRIX, w->pcap, w->idx);
		if(w->score != x - w->de) { break; }
		debug("del, score(%u), x(%u), de(%u)", w->score, x, w->de);

		dz_trace_push_op(w, DZ_E_MATRIX, x);
		dz_trace_unwind_h(w, DZ_E_MATRIX);
	} while(!dz_trace_test_idx(w, 0));

	/* eat last row */
	dz_trace_push_op(w, DZ_E_MATRIX, dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx));
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
void dz_trace_init_aln(dz_alignment_t *aln, dz_state_t const *ff, size_t idx)
{
	aln->rrem  = dz_calc_max_rpos_core(ff);
	aln->score = ff->max.score;
	aln->query_length = idx;

	debug("ff(%p), idx(%zu)", ff, idx);
	debug("score(%d)", ff->max.score);
	return;
}

static __dz_vectorize
void dz_trace_allocate_aln(dz_trace_work_t *w, dz_arena_t *mem, dz_state_t const *ff, size_t idx)
{
	/* allocate aln object */
	size_t aln_size = (sizeof(dz_alignment_t)
		+ (ff->cnt.section + 6) * sizeof(dz_path_span_t)
		+ dz_roundup(ff->cnt.column + idx + 1, 8) * sizeof(uint8_t)			/* +1 for tail '\0' */
	);
	w->aln = (dz_alignment_t *)dz_arena_malloc(mem, aln_size);
	dz_trace_init_aln(w->aln, ff, idx);

	/* slice section and path */
	dz_path_span_t *span = (dz_path_span_t *)(w->aln + 1);
	w->span.ptr  = span + ff->cnt.section + 4;
	w->span.base = span + ff->cnt.section + 4;

	uint8_t *path = (uint8_t *)(w->span.base + 2);
	w->path.ptr  = path + ff->cnt.column + idx;
	w->path.base = path + ff->cnt.column + idx;

	/* push tail sentinels */
	dz_tail_t const *tail = dz_restore_tail(ff);
	dz_trace_push_span(w, tail->rsave.id);
	*w->path.ptr = '\0';		/* make sure readable as C string */
	return;
}

static __dz_vectorize
void dz_trace_init_work(dz_trace_work_t *w, dz_profile_t const *profile, dz_state_t const *ff, size_t idx)
{
	/* save lengths */
	w->rlen    = 0;		/* is unknown before traceback */
	w->query   = dz_restore_tail(ff)->query;
	w->profile = profile;
	w->ie = profile->iev[0];
	w->de = profile->dev[0];

	/* clear counters */
	memset(w->cnt, 0, 4 * sizeof(uint32_t));

	/* load max column pointers */
	w->pcap = ff->max.cap;
	w->cap  = NULL;
	w->idx  = idx;
	w->score = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx);
	return;
}

static __dz_vectorize
void dz_trace_finalize_path(dz_alignment_t *aln, uint8_t *path, uint8_t *base)
{
	aln->path = path;
	aln->path_length = base - path;
	return;
}

static __dz_vectorize
void dz_trace_finalize_span(dz_alignment_t *aln, dz_path_span_t *span, dz_path_span_t *base, size_t plen)
{
	span++;				/* adjust root span */
	size_t const slen = base - span;

	/* save pointers */
	aln->span = span;
	aln->span_length = slen;

	/* offset fixup */
	span[0].offset = 0;
	for(size_t i = 1; i < slen; i++) {
		debug("i(%zu), plen(%zu), offset(%u, %u)", i, plen, span[i].offset, (uint32_t)(plen - span[i].offset));
		span[i].offset = plen - span[i].offset;
	}
	span[slen].offset = plen;
	return;
}

static __dz_vectorize
dz_alignment_t *dz_trace_finalize(dz_trace_work_t *w)
{
	/* check buffer overrun */
	assert((uintptr_t)(w->aln + 1) <= (uintptr_t)w->span.ptr + 1);
	assert((uintptr_t)w->span.base <= (uintptr_t)w->path.ptr);

	dz_trace_finalize_path(w->aln, w->path.ptr, w->path.base);
	dz_trace_finalize_span(w->aln, w->span.ptr, w->span.base, w->aln->path_length);

	/* save length */
	w->aln->ref_length = w->rlen;

	/* save match / gap counts */
	w->aln->del_count = w->cnt[DZ_E_MATRIX];
	w->aln->ins_count = w->cnt[DZ_F_MATRIX];
	w->aln->mismatch_count = w->cnt[DZ_S_MATRIX];
	w->aln->match_count = w->cnt[DZ_S_MATRIX + 1];
	return(w->aln);
}

static __dz_vectorize
dz_alignment_t *dz_trace_core(dz_arena_t *mem, dz_profile_t const *profile, dz_state_t const *ff)
{
	/* no object is returned for root */
	if(ff->max.cap == NULL) { return(NULL); }

	/* detect pos */
	uint64_t const idx = dz_calc_max_qpos_core(ff);	/* vector index, cell index */
	debug("idx(%zu), root(%p)", (size_t)idx, profile->root);

	/* init working buffer */
	dz_trace_work_t w;
	dz_trace_init_work(&w, profile, ff, idx);
	dz_trace_allocate_aln(&w, mem, ff, idx);

	/* core loop */
	dz_trace_unwind_h(&w, DZ_S_MATRIX);
	while(!dz_is_root(w.pcap)) {
		if(dz_trace_eat_match(&w)) { continue; }
		if(dz_trace_eat_ins(&w)) { continue; }
		if(dz_trace_eat_del(&w)) { continue; }

		dz_trap();
	}

	/* done */
	return(dz_trace_finalize(&w));
}

static __dz_vectorize
dz_alignment_t *dz_trace(dz_t *self, dz_forefront_t const *ff)
{
	return(dz_trace_core(self->mem, self->profile, dz_cstate(ff)));
}

#if defined(UNITTEST) && UNITTEST != 0

/* short, exact matching sequences */
unittest( "trace" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefront = NULL;
	dz_alignment_t *aln = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, "A", 1, 1);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "GC", 2, 2);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "TTTT", 4, 3);
	aln = dz_trace(dz, forefront);

	(void)aln;
	dz_destroy(dz);
}

#endif	/* defined(UNITTEST) && UNITTEST != 0 */


#ifdef __cplusplus
};	/* extern "C" { */
#endif

/**
 * end of dozeu.h
 */
