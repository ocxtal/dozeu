// $(CC) -O3 -march=native -DMAIN -o dozeu dozeu.c

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
 *   SSSE3 (Core 2 / Bobcat or later)
 *   { gcc, clang, icc } x { Linux, Mac OS }
 * input sequence must be null-terminated array of {'A','C','G','T','a','c','g','t'}, otherwise undefined, not required to be terminated with '\0'.
 *
 * BLAST X-drop DP
 *   https://github.com/elucify/blast-docs/wiki/The-Developer's-Guide-to-BLAST (the most datailed document that describes the original X-drop DP algorithm)
 *
 * Myers bit vector (describes vertical (horizontal) tiling)
 *   Gene Myers, A fast bit-vector algorithm for approximate string matching based on dynamic programming, JACM (1999)
 */
#include <assert.h>
#include <stdint.h>
#include <stddef.h>
#include <x86intrin.h>
#include "log.h"

#define UNITTEST_ALIAS_MAIN			0
#define UNITTEST_UNIQUE_ID			3213
#include "unittest.h"

unittest_config( .name = "dozeu" );
unittest() { debug("hello"); }

enum dz_alphabet { A = 0x00, C = 0x01, G = 0x02, T = 0x03, N = 0x04 };		/* internal encoding */

#define dz_unused(x)				(void)(x)
#define dz_likely(x)				__builtin_expect(!!(x), 1)
#define dz_unlikely(x)				__builtin_expect(!!(x), 0)
#define dz_roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )
#define dz_rounddown(x, base)		( (x) & ~((base) - 1) )
#define dz_max2(x, y)				( (x) < (y) ? (y) : (x) )
#define dz_min2(x, y)				( (x) < (y) ? (x) : (y) )

#define dz_sa_cat_intl(x, y)		x##y
#define dz_sa_cat(x, y)				dz_sa_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_sa_cat(_st, __LINE__)[(expr) ? 1 : -1]

#define DZ_MEM_MARGIN_SIZE			( 256 )
#define DZ_MEM_ALIGN_SIZE			( 16 )
#define DZ_MEM_INIT_SIZE			( 16 * 1024 * 1024 )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_CELL_MARGIN				( 32 )
#define DZ_CELL_MARGINED_MIN		( DZ_CELL_MIN + DZ_CELL_MARGIN )
#define DZ_CELL_MARGINED_MAX		( DZ_CELL_MAX - DZ_CELL_MARGIN )

/* query */
struct dz_query_s { uint64_t blen, _pad; uint8_t arr[]; };		/* preconverted query sequence; blen = roundup(qlen, L) / L; array must have 16-byte-length margin at the tail */
dz_static_assert(sizeof(struct dz_query_s) % sizeof(__m128i) == 0);

/* DP matrix structures */
struct dz_swgv_s { __m128i e, f, s; };							/* followed by dz_cap_s */
struct dz_cap_s { uint32_t spos, epos, rch, n_tails; };			/* followed by dz_tail_s; spos and epos are shared to tail_s */
struct dz_tail_s { uint32_t spos, epos, max, inc; struct dz_query_s const *query; struct dz_cap_s const *cap; };
dz_static_assert(sizeof(struct dz_swgv_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_cap_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_tail_s) % sizeof(__m128i) == 0);
#define dz_swgv(_p)					( (struct dz_swgv_s *)(_p) )
#define dz_cswgv(_p)				( (struct dz_swgv_s const *)(_p) )
#define dz_cap(_p)					( (struct dz_cap_s *)(_p) )
#define dz_ccap(_p)					( (struct dz_cap_s const *)(_p) )
#define dz_tail(_p)					( (struct dz_tail_s *)(_p) )
#define dz_ctail(_p)				( (struct dz_tail_s const *)(_p) )

/* context (constants and working buffers) */
struct dz_mem_block_s { struct dz_mem_block_s *next; size_t size; };
struct dz_stack_s { struct dz_mem_block_s *curr; uint8_t *top, *end; uint64_t _pad[3]; };
struct dz_mem_s { struct dz_mem_block_s blk; struct dz_stack_s stack; };
#define dz_mem_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )

struct dz_s { int8_t matrix[16]; uint16_t giv[8], gev[8], xt, _pad[11]; struct dz_tail_s const *root; };
dz_static_assert(sizeof(struct dz_s) % sizeof(__m128i) == 0);
#define dz_mem(_self)				( (struct dz_mem_s *)(_self) - 1 )
#define dz_root(_self)				( (struct dz_tail_s const **)(&_self->root) )

#define print_vector(v) { \
	debug("%s (%d, %d, %d, %d, %d, %d, %d, %d)", #v, \
	(int16_t)_mm_extract_epi16(v, 7), \
	(int16_t)_mm_extract_epi16(v, 6), \
	(int16_t)_mm_extract_epi16(v, 5), \
	(int16_t)_mm_extract_epi16(v, 4), \
	(int16_t)_mm_extract_epi16(v, 3), \
	(int16_t)_mm_extract_epi16(v, 2), \
	(int16_t)_mm_extract_epi16(v, 1), \
	(int16_t)_mm_extract_epi16(v, 0)); \
}

/**
 * @fn dz_malloc, dz_free
 * @brief aligned and margined malloc and free
 */
static
void *dz_malloc(
	size_t size)
{
	void *ptr = NULL;

	/* roundup to align boundary, add margin at the head and tail */
	size = dz_roundup(size, DZ_MEM_ALIGN_SIZE);
	if(posix_memalign(&ptr, DZ_MEM_ALIGN_SIZE, size + 2 * DZ_MEM_MARGIN_SIZE) != 0) {
		debug("posix_memalign failed");
		trap(); return(NULL);
	}
	debug("posix_memalign(%p), size(%lu)", ptr, size);
	return(ptr + DZ_MEM_MARGIN_SIZE);
}
static
void dz_free(
	void *ptr)
{
	debug("free(%p)", ptr - DZ_MEM_MARGIN_SIZE);
	free(ptr - DZ_MEM_MARGIN_SIZE);
	return;
}

unittest() {
	uint64_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(uint64_t i = 0; i < sizeof(size) / sizeof(uint64_t); i++) {
		uint8_t *p = dz_malloc(size[i]);
		ut_assert(p != NULL);
		volatile uint8_t a = p[-32], b = p[0], c = p[size[i]], d = p[size[i] + 32];
		dz_unused(a); dz_unused(b); dz_unused(c); dz_unused(d);
		memset(p, 0, size[i]);
		dz_free(p);
	}
}

/**
 * @fn dz_mem_init, dz_mem_destroy, dz_mem_add_stack, dz_mem_malloc
 * @brief stack chain
 */
static
struct dz_mem_s *dz_mem_init(
	size_t size)
{
	struct dz_mem_s *mem = (struct dz_mem_s *)dz_malloc(sizeof(struct dz_mem_s) + dz_roundup(size, DZ_MEM_ALIGN_SIZE));
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}

	/* init stack pointers */
	mem->stack.curr = &mem->blk;
	mem->stack.top = (uint8_t *)(mem + 1);
	mem->stack.end = (uint8_t *)mem + DZ_MEM_INIT_SIZE;

	/* init mem object */
	mem->blk = (struct dz_mem_block_s){ .next = NULL, .size = DZ_MEM_INIT_SIZE };
	return(mem);
}
static
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
	return;
}
static
uint64_t dz_mem_add_stack(
	struct dz_mem_s *mem,
	size_t size)
{
	debug("add_stack, ptr(%p)", mem->stack.curr->next);
	if(mem->stack.curr->next == NULL) {
		/* current stack is the tail of the memory block chain, add new block */
		size = dz_max2(
			size + dz_roundup(sizeof(struct dz_mem_block_s), DZ_MEM_ALIGN_SIZE),
			2 * mem->stack.curr->size
		);
		struct dz_mem_block_s *blk = dz_malloc(size);
		debug("malloc called, blk(%p)", blk);
		if(blk == NULL) { return(1); }

		/* link new node to the tail of the current chain */
		mem->stack.curr->next = blk;

		/* initialize the new memory block */
		blk->next = NULL;
		blk->size = size;
	}
	/* follow a forward link and init stack pointers */
	mem->stack.curr = mem->stack.curr->next;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem->stack.curr + 1), DZ_MEM_ALIGN_SIZE);
	mem->stack.end = (uint8_t *)mem->stack.curr + mem->stack.curr->size;
	return(0);
}
static
void *dz_mem_malloc(
	struct dz_mem_s *mem,
	size_t size)
{
	if(dz_mem_stack_rem(mem) < 4096) { dz_mem_add_stack(mem, 0); }
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += size;
	return(ptr);
}

unittest() {
	uint64_t size[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000 };
	for(uint64_t i = 0; i < sizeof(size) / sizeof(uint64_t); i++) {
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
#define _begin_matrix(_query) ({ \
	; \
})
#define _end_matrix(_p, _wp) ({ \
	/* create tail object */ \
	struct dz_tail_s *tail = (struct dz_tail_s *)&((struct dz_swgv_s *)(_p))[(_wp)->epos]; \
	tail->spos = (_wp)->spos; tail->epos = (_wp)->epos; \
	tail->max = (_wp)->max + (_wp)->inc; \
	tail->inc = (_wp)->inc; \
	tail->query = (_wp)->query; \
	tail->cap = (_wp)->cap; \
	debug("[%u, %u), inc(%u), cap(%p), max(%u)", w.spos, w.epos, w.inc, w.cap, w.max); \
	/* write back stack pointer */ \
	dz_mem(self)->stack.top = (uint8_t *)(tail + 1); \
	/* return the tail pointer */ \
	(struct dz_tail_s const *)tail; \
})
#define _begin_column_head(_spos, _epos, _rch, _tails, _n_tails) ({ \
	size_t _sp = (_spos), _ep = (_epos), _len = (_ep) - (_sp), _nt = (_n_tails); \
	/* calculate sizes */ \
	size_t tail_arr_size = dz_roundup(sizeof(struct dz_tail_s *) * _nt, sizeof(__m128i)); \
	size_t est_column_size = 2 * _len * sizeof(struct dz_swgv_s); \
	size_t next_req = tail_arr_size + est_column_size + sizeof(struct dz_cap_s); \
	debug("est_column_size(%lu), next_req(%lu)", est_column_size, next_req); \
	/* allocate from heap */ \
	if(dz_mem_stack_rem(dz_mem(self)) < next_req) { dz_mem_add_stack(dz_mem(self), 0); } \
	/* push tail pointers */ \
	struct dz_tail_s const **tail_arr = (struct dz_tail_s const **)(dz_mem(self)->stack.top + tail_arr_size); \
	for(uint64_t i = 0; i < _nt; i++) { tail_arr[-((int64_t)_nt) + i] = (struct dz_tail_s const *)_tails[i]; } \
	/* push cap info */ \
	struct dz_cap_s *cap = (struct dz_cap_s *)tail_arr; \
	cap->spos = cap->epos = 0;		/* head marked as zero */ \
	cap->rch = (_rch);				/* rch for the first column */ \
	cap->n_tails = _nt;			/* record n_tails */ \
	/* return the head pointer */ \
	(struct dz_swgv_s *)(cap + 1) - (_sp); \
})
#define _begin_column(_spos, _epos, _rch) ({ \
	size_t _sp = (_spos), _ep = (_epos), len = (_ep) - (_sp); \
	/* calculate sizes */ \
	size_t est_column_size = 2 * len * sizeof(struct dz_swgv_s); \
	size_t next_req = est_column_size + sizeof(struct dz_cap_s); \
	debug("est_column_size(%lu), next_req(%lu)", est_column_size, next_req); \
	/* allocate from heap */ \
	if(dz_mem_stack_rem(dz_mem(self)) < next_req) { dz_mem_add_stack(dz_mem(self), 0); } \
	/* push cap info */ \
	struct dz_cap_s *cap = (struct dz_cap_s *)dz_mem(self)->stack.top; \
	cap->spos = _sp; cap->epos = _ep;/* [spos, epos) for the previous column */ \
	cap->rch = (_rch);				/* record the rch for the next column */ \
	cap->n_tails = 0;				/* n_tails == 0 for internal caps */ \
	/* return the head pointer */ \
	(struct dz_swgv_s *)(cap + 1) - (_sp); \
})
#define _end_column(_p, _epos) ({ \
	/* write back the stack pointer and return a cap */ \
	uint8_t *_t = (uint8_t *)&((struct dz_swgv_s *)(_p))[_epos]; \
	dz_mem(self)->stack.top = _t; \
	(struct dz_cap_s *)_t; \
})
#define _load_vector(_p) \
	__m128i e = _mm_load_si128((__m128i const *)(&((struct dz_swgv_s const *)(_p))->e)); print_vector(e); \
	__m128i s = _mm_load_si128((__m128i const *)(&((struct dz_swgv_s const *)(_p))->s)); print_vector(s);
#define _update_vector(_p) { \
	__m128i qv = _mm_loadl_epi64((__m128i const *)&qp[(_p) * L]); \
	__m128i sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)self->matrix), _mm_or_si128(rv, qv))); print_vector(sc); \
	__m128i te = _mm_subs_epi16(_mm_max_epi16(e, _mm_subs_epi16(s, giv)), gev1); \
	 print_vector(_mm_alignr_epi8(s, ps, 14));  \
	__m128i ts = _mm_max_epi16(te, _mm_adds_epi16(sc, _mm_alignr_epi8(s, ps, 14))); ps = s; \
	__m128i tf = _mm_max_epi16(_mm_subs_epi16(ts, giv), _mm_subs_epi16(_mm_alignr_epi8(minv, f, 14), gev1)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 14), gev1)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 12), gev2)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 8), gev4)); \
	ts = _mm_max_epi16(ts, tf); \
	maxv = _mm_max_epi16(maxv, ts); \
	/* print_vector(te); */ print_vector(ts); /* print_vector(tf); print_vector(maxv); */ \
	e = te; f = tf; s = ts; \
}
#define _store_vector(_p) { \
	_mm_store_si128((__m128i *)(&((struct dz_swgv_s *)(_p))->e), e); \
	_mm_store_si128((__m128i *)(&((struct dz_swgv_s *)(_p))->f), f); \
	_mm_store_si128((__m128i *)(&((struct dz_swgv_s *)(_p))->s), s); \
}
#define _hmax_vector(_v) ({ \
	__m128i _t = _mm_max_epi16(_v, _mm_bsrli_si128(_v, 8)); \
	_t = _mm_max_epi16(_t, _mm_bsrli_si128(_t, 4)); \
	_t = _mm_max_epi16(_t, _mm_bsrli_si128(_t, 2)); \
	(int64_t)((int16_t)(_mm_extract_epi16(_t, 0))); \
})
#define _test_xdrop(_s, _xtv) ({ \
	__m128i xtest = _mm_cmpgt_epi16(_s, _xtv); \
	/* print_vector(_s); print_vector(_xtv); */ \
	_mm_test_all_zeros(xtest, xtest); \
})

/* obsolete * /
#define _push_cap(_p) ({ \
	struct dz_cap_s *cap = (struct dz_cap_s *)&(_p)[w.epos]; \
	cap->spos = w.spos; cap->epos = w.epos; \
	_mm_storel_epi64((__m128i *)&cap->rch, rv); \
	dz_mem(self)->stack.top = (uint8_t *)(cap + 1); \
	cap; \
})
*/

/**
 * @fn dz_init, dz_destroy
 */
static
struct dz_s *dz_init(
	int8_t const *score_matrix,		/* match award in positive, mismatch penalty in negative. s(A,A) at [0], s(A,C) at [1], ... s(T,T) at [15] where s(ref_base, query_base) is a score function */
	uint16_t gap_open, uint16_t gap_extend,		/* gap penalties in positive */
	uint64_t max_gap_len)			/* as X-drop threshold */
{
	uint64_t const L = sizeof(__m128i) / sizeof(uint16_t), gi = gap_open, ge = gap_extend;
	struct dz_mem_s *mem = dz_mem_init(DZ_MEM_INIT_SIZE);
	if(mem == NULL) {
		debug("failed to malloc memory");
		return(NULL);
	}
	struct dz_s *self = dz_mem_malloc(mem, sizeof(struct dz_s));

	/* constants */
	__m128i const giv = _mm_set1_epi16(gi);
	__m128i const gev = _mm_set1_epi16(ge);
	_mm_store_si128((__m128i *)self->matrix, _mm_loadu_si128((__m128i const *)score_matrix));
	_mm_store_si128((__m128i *)self->giv, giv);
	_mm_store_si128((__m128i *)self->gev, gev);
	self->xt = gi + ge * max_gap_len;			/* X-drop threshold */

	/* create root cap */
	struct dz_cap_s *cap = dz_mem_malloc(mem, sizeof(struct dz_cap_s));
	_mm_store_si128((__m128i *)cap, _mm_setzero_si128());

	/* calc vector length */
	max_gap_len = dz_roundup(max_gap_len, L);
	struct dz_tail_s w = { .epos = max_gap_len / L, .query = NULL }, *a = &w;		/* query = NULL for the first (root) column */

	/* malloc the first column */
	struct dz_swgv_s *dp = _begin_column_head(0, max_gap_len / L, 0, &a, 0);

	/* fill the root (the leftmost) column; first init vectors */
	__m128i s = _mm_setr_epi16(0, -(gi+ge), -(gi+2*ge), -(gi+3*ge), -(gi+4*ge), -(gi+5*ge), -(gi+6*ge), -(gi+7*ge));
	__m128i const e = _mm_set1_epi16(DZ_CELL_MIN), xtv = _mm_set1_epi16(-self->xt);

	/* until the X-drop test fails on all the cells in a vector */
	for(uint64_t p = 0; p < max_gap_len / L; p++) {
		__m128i const f = s;
		if(dz_unlikely(_test_xdrop(s, xtv))) { debug("p(%lu)", p); w.epos = p; break; }
		_store_vector(&dp[p]);
		if(p == 0) { s = _mm_setr_epi16(-gi, -(gi+ge), -(gi+2*ge), -(gi+3*ge), -(gi+4*ge), -(gi+5*ge), -(gi+6*ge), -(gi+7*ge)); }
		s = _mm_subs_epi16(s, _mm_slli_epi16(gev, 3));
	}

	/* done; create and return a tail object */
	self->root = _end_matrix(dp, &w);
	return(self);
}
static
void dz_destroy(
	struct dz_s *self)
{
	if(self == NULL) { return; }
	dz_mem_destroy(dz_mem(self));
	return;
}

#define DZ_UNITTEST_SCORE_PARAMS	(int8_t const []){ 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2, -3, -3, -3, -3, 2 }, 5, 1, 20
static uint8_t const *dz_unittest_query = (uint8_t const *)
	"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
	"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
	"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
	"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
	"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
	"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
	"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
	"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);
	ut_assert(dz->root != NULL);
	dz_destroy(dz);
}

/**
 * @fn dz_pack_query
 */
static
struct dz_query_s *dz_pack_query(
	struct dz_s *self,
	uint8_t const *query, uint64_t qlen)
{
	uint64_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + dz_roundup(qlen, sizeof(__m128i)));
	q->blen = dz_roundup(qlen, L) / L;
	q->arr[0] = 0;

	/* ASCII to 2-bit conversion table (shifted by two bits to be used as the upper (row) shuffle index) */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		['A' & 0x0f] = A<<2, ['C' & 0x0f] = C<<2,
		['G' & 0x0f] = G<<2, ['T' & 0x0f] = T<<2,
		['U' & 0x0f] = T<<2, ['N' & 0x0f] = N<<2
	};
	__m128i pv = _mm_setzero_si128();
	__m128i const cv = _mm_load_si128((__m128i const *)conv), fv = _mm_set1_epi8(0x0f);			/* conversion table and mask */

	/* until the end of the query sequence */
	for(uint64_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
		__m128i tv = _mm_shuffle_epi8(cv, _mm_and_si128(qv, fv));
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
	}

	/* continue the same conversion on the remainings */
	for(uint64_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[query[i]];
	}
	return(q);
}

unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, strlen((char const *)dz_unittest_query));
	dz_unused(q);
	dz_destroy(dz);
}

/**
 * @fn dz_extend
 */
static
struct dz_tail_s const *dz_extend(
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_tail_s const **tails, uint64_t n_tails,
	uint8_t const *ref, int64_t rlen)
{
	uint64_t const L = sizeof(__m128i) / sizeof(uint16_t);
	if(n_tails == 0 || rlen == 0) { return(NULL); }

	/* load constants */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		['A' & 0x0f] = A, ['C' & 0x0f] = C,
		['G' & 0x0f] = G, ['T' & 0x0f] = T,
		['U' & 0x0f] = T, ['N' & 0x0f] = N
	};
	__m128i const minv = _mm_set1_epi16(DZ_CELL_MIN);
	__m128i const giv = _mm_load_si128((__m128i const *)self->giv);
	__m128i const gev1 = _mm_load_si128((__m128i const *)self->gev);
	__m128i const gev2 = _mm_add_epi16(gev1, gev1);
	__m128i const gev4 = _mm_slli_epi16(gev1, 2);
	__m128i const gev8 = _mm_slli_epi16(gev1, 3);
	uint8_t const *qp = query->arr;
	uint64_t const qblen = query->blen;

	/* first iterate over the incoming edge objects to get the current max */
	uint64_t adj[n_tails];								/* S[0, 0] = 0 */
	struct dz_tail_s w = { .spos = UINT32_MAX, .query = query };
	for(uint64_t i = 0; i < n_tails; i++) {
		/* update max and pos */
		w.max = dz_max2(w.max, tails[i]->max);
		w.spos = dz_min2(w.spos, tails[i]->spos);
		w.epos = dz_max2(w.epos, tails[i]->epos);
		debug("i(%lu), [%u, %u), inc(%u), max(%u)", i, tails[i]->spos, tails[i]->epos, tails[i]->inc, tails[i]->max);
	}
	debug("start extension [%u, %u), inc(%u), max(%u), stack(%p, %p)", w.spos, w.epos, w.inc, w.max, dz_mem(self)->stack.top, dz_mem(self)->stack.end);

	for(uint64_t i = 0; i < n_tails; i++) {
		adj[i] = w.max - (tails[i]->max - tails[i]->inc);/* base = max - inc */
		debug("i(%lu), adj(%lu)", i, adj[i]);
	}

	/* X-drop DP */
	ptrdiff_t dir = rlen < 0 ? -1 : 1;					/* direction */
	uint8_t const *rp = rlen < 0 ? &ref[rlen - 1] : ref;
	uint8_t rch = conv[*rp & 0x0f];						/* load the first base */
	struct dz_swgv_s *pdp = _begin_column_head(w.spos, w.epos, rch, tails, n_tails);	/* allocate memory for the first column */
	for(uint64_t p = w.spos; p < w.epos; p++) {
		/* memset(pdp, 0xff, sizeof(struct dz_swgv_s) * (w.epos - w.spos)); */
		__m128i const e = _mm_set1_epi16(INT16_MIN), f = e, s = e; _store_vector(&pdp[p]);
	}

	/* paste the last vectors */
	for(uint64_t i = 0; i < n_tails; i++) {
		struct dz_swgv_s const *tdp = (struct dz_swgv_s const *)tails[i] - tails[i]->epos;
		__m128i const adjv = _mm_set1_epi16(adj[i]);
		for(uint64_t p = tails[i]->spos; p < tails[i]->epos; p++) {
			/* adjust offset */
			__m128i e = _mm_subs_epi16(_mm_load_si128(&tdp[p].e), adjv);
			__m128i f = _mm_subs_epi16(_mm_load_si128(&tdp[p].f), adjv);
			__m128i s = _mm_subs_epi16(_mm_load_si128(&tdp[p].s), adjv);

			/* read-max-write */
			_mm_store_si128(&pdp[p].e, _mm_max_epi16(e, _mm_load_si128(&pdp[p].e)));
			_mm_store_si128(&pdp[p].f, _mm_max_epi16(f, _mm_load_si128(&pdp[p].f)));
			_mm_store_si128(&pdp[p].s, _mm_max_epi16(s, _mm_load_si128(&pdp[p].s)));
		}
	}

	/* once; extend the first vector */ {
		/* init vectors; rch (reference-side base) is previously loaded on rch */
		__m128i f = minv, ps = minv, maxv = _mm_setzero_si128();
		__m128i const xtv = _mm_set1_epi16(-self->xt), rv = _mm_set1_epi8(rch);	/* next offset == current max thus X-drop threshold is always -xt */

		/* until the bottommost vertically placed band... */
		for(uint64_t p = w.spos; p < w.epos; p++) {
			_load_vector(&pdp[p]); _update_vector(p);
			if(dz_unlikely(_test_xdrop(s, xtv))) {	/* tag _unlikely to move out of the core loop */
				/* drop either leading or trailing vector, skip the tail extension when the tail is clipped */
				if(p == w.spos) { w.spos++; } else { w.epos = p; goto _merge_tail; }
			}
			_store_vector(&pdp[p]);
		}

		/* if reached the tail of the query sequence, finish the extension */
		if(w.epos >= qblen) { goto _merge_tail; }

		/* tail extension; clip the column length if too long */
		__m128i e = minv, s = minv; _update_vector(w.epos);
		do {
			if(_test_xdrop(s, xtv)) { break; }
			_store_vector(&pdp[w.epos]); w.epos++;
			f = _mm_subs_epi16(f, gev8); s = _mm_subs_epi16(s, gev8);
		} while(w.epos < qblen);
	_merge_tail:;
		struct dz_cap_s *cap = _end_column(pdp, w.epos);
		uint64_t inc = _hmax_vector(maxv);
		if(inc > w.inc) { w.inc = inc; w.cap = cap; }					/* 2 x cmov */
		debug("update inc(%lu), max(%u, %u, %p), cap(%p)", w.inc, w.max, w.max + w.inc, w.cap, cap);
	}

	/* extension loop */
	debug("rp(%p), dir(%ld)", rp, dir);
	while((rp += dir, --rlen > 0)) {
		debug("rp(%p), [%u, %u), stack(%p, %p)", rp, w.spos, w.epos, dz_mem(self)->stack.top, dz_mem(self)->stack.end);
		/* load next base */
		rch = conv[*rp & 0x0f];

		/* open a new column (reallocate if needed) */
		struct dz_swgv_s *dp = _begin_column(w.spos, w.epos, rch);

		/* init vectors */
		__m128i f = minv, ps = minv, maxv = _mm_setzero_si128();
		__m128i const xtv = _mm_set1_epi16(w.inc - self->xt), rv = _mm_set1_epi8(rch);

		/* until the bottommost vertically placed band... */
		for(uint64_t p = w.spos; p < w.epos; p++) {
			_load_vector(&pdp[p]); _update_vector(p);
			if(dz_unlikely(_test_xdrop(s, xtv))) {
				if(p == w.spos) { w.spos++; } else { w.epos = p + 1; goto _loop_tail; }
			}
			_store_vector(&dp[p]);
		}

		/* if reached the tail of the query sequence */
		if(w.epos >= qblen) { goto _loop_tail; }

		/* if not, extend until X-drop test fails for all the cells in a vector */
		__m128i e = minv, s = minv; _update_vector(w.epos);
		do {
			if(_test_xdrop(s, xtv)) { break; }
			_store_vector(&dp[w.epos]); w.epos++;
			f = _mm_subs_epi16(f, gev8); s = _mm_subs_epi16(s, gev8);
		} while(w.epos < qblen);

		/* create cap object that contains [spos, epos) range (for use in the traceback routine) */
	_loop_tail:;
		struct dz_cap_s *cap = _end_column(dp, w.epos); pdp = dp;
		uint64_t inc = _hmax_vector(maxv);
		if(inc > w.inc) { w.inc = inc; w.cap = cap; }					/* update max; the actual score (absolute score accumulated from the origin) is expressed as max + inc; 2 x cmov */
		debug("update inc(%lu), max(%u, %u, %p), cap(%p)", w.inc, w.max, w.max + w.inc, w.cap, cap);
		/* FIXME: rescue overflow */
	}

	/* update max vector (pdp always points at the last vector) */
	return(_end_matrix(pdp, &w));
}

/* short, exact matching sequences */
unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s const *q = dz_pack_query(dz, dz_unittest_query, strlen((char const *)dz_unittest_query));
	struct dz_tail_s const *tail = NULL;

	/* nothing occurs */
	tail = dz_extend(dz, q, NULL, 0, NULL, 0);
	ut_assert(tail == NULL);
	tail = dz_extend(dz, q, NULL, 0, (uint8_t const *)"", 0);
	ut_assert(tail == NULL);
	tail = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"", 0);
	ut_assert(tail == NULL);

	/* extend */
	tail = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"A", 1);
	ut_assert(tail != NULL && tail->max == 2);

	tail = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AG", 2);
	ut_assert(tail != NULL && tail->max == 4);

	tail = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AGATTTT", 7);
	ut_assert(tail != NULL && tail->max == 9);

	dz_destroy(dz);
}

/* a small graph */
unittest( .name = "small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, strlen((char const *)dz_unittest_query));
	struct dz_tail_s const *tails[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	tails[0] = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AG", 2);
	ut_assert(tails[0] != NULL && tails[0]->max == 4);

	tails[1] = dz_extend(dz, q, &tails[0], 1, (uint8_t const *)"C", 1);
	ut_assert(tails[1] != NULL && tails[1]->max == 6);

	tails[2] = dz_extend(dz, q, &tails[0], 2, (uint8_t const *)"TTTT", 4);
	ut_assert(tails[2] != NULL && tails[2]->max == 14);

	tails[3] = dz_extend(dz, q, &tails[2], 1, (uint8_t const *)"CATT", 4);
	ut_assert(tails[3] != NULL && tails[3]->max == 22);

	tails[4] = dz_extend(dz, q, &tails[2], 2, (uint8_t const *)"CTGA", 4);
	ut_assert(tails[4] != NULL && tails[4]->max == 30);

	dz_destroy(dz);
}

/**
 * @fn dz_trace
 */
static
uint8_t *dz_trace(
	struct dz_s *self,
	struct dz_tail_s const *tail)
{
	uint64_t const L = sizeof(__m128i) / sizeof(uint16_t);
	#define _is_head(_cap)				( *((uint64_t *)(_cap)) == 0 )		/* cap->spos == 0 && cap->epos == 0 */
	#define _dp(_cap)					( (struct dz_swgv_s const *)(_cap) - (_cap)->epos )
	#define _score(_l, _cap, _idx)		( ((int16_t const *)(&_dp(_cap)[(_idx) / L]._l))[(_idx) & (L - 1)] )	

	/* load max column pointers */
	struct dz_cap_s const *pcap = tail->cap, *cap = NULL;
	struct dz_query_s const *query = tail->query;

	/* detect pos */
	uint64_t idx = 0;										/* vector index, cell index */
	__m128i const maxv = _mm_set1_epi16(tail->inc);
	for(uint64_t p = pcap->spos; p < pcap->epos; p++) {
		__m128i const s = _mm_load_si128(&_dp(pcap)[p].s); print_vector(s);
		uint64_t eq = _mm_movemask_epi8(_mm_cmpeq_epi16(s, maxv));
		if(eq != 0) {
			/* tzcntq is faster but avoid using it b/c it requires relatively newer archs */
			uint64_t zcnt = (eq - 1) & (~eq & 0x5555);		/* subq, andnq, andq; chain length == 2 */
			zcnt += zcnt>>2; zcnt &= 0x3333;				/* shrq, addq, andq; chain length might be shorter */
			zcnt += zcnt>>4; zcnt &= 0x0f0f;
			zcnt += zcnt>>8; zcnt &= 0x00ff;
			debug("found, eq(%lx), zcnt(%lu), idx(%lu)", eq, zcnt, p * L + zcnt);
			idx = p * L + zcnt; break;
		}
	}

	/* traceback loop */
	uint32_t score = tail->inc;
	while(1) {
		#define _prev_cap() ({ \
			/* bare the prevous column */ \
			cap = pcap; pcap = (struct dz_cap_s const *)((struct dz_swgv_s const *)cap - cap->epos + cap->spos) - 1; \
			debug("cap(%p), pcap(%p), [%u, %u)", cap, pcap, pcap->spos, pcap->epos); \
			/* load the reference-side base */ \
			uint8_t rch = pcap->rch; \
			/* reload pcap pointer if root */ \
			if(_is_head(pcap)) { \
				debug("head, n_tails(%u)", pcap->n_tails); if(pcap->n_tails == 0) { debug("root"); break; } \
				struct dz_tail_s const **tail_arr = (struct dz_tail_s const **)pcap; \
				debug("idx(%lu), n_tails(%u), tail_arr(%p)", idx, pcap->n_tails, tail_arr); \
				pcap = (struct dz_cap_s const *)(tail_arr[-((int64_t)pcap->n_tails)]); \
			} \
			/* return the reference-side base */ \
			rch; \
		})
		debug("loop, idx(%lu)", idx);
		uint8_t rch = _prev_cap();
		debug("diagonal, idx(%lu), rbase(%x, %c), qbase(%c), c(%d), e(%d), f(%d), s(%d), score(%d)",
			idx, rch, "ACGTNNNNNNNNNNNN"[rch & 0xff], "ANNNCNNNGNNNTNNN"[query->arr[idx]],
			score, _score(e, cap, idx), _score(f, cap, idx), _score(s, pcap, idx - 1), self->matrix[rch | query->arr[idx]]);
		idx--;
		if(score == _score(s, pcap, idx - 1) + self->matrix[rch | query->arr[idx]]) {
			debug("hit diagonal"); continue;
		}
		if(score == _score(f, cap, idx)) {
			do {
				_prev_cap();
			} while(score == _score(f, cap, idx - 1) - self->gev[0]);
		}
		if(score == _score(e, cap, idx)) {
			do {
				_prev_cap();
			} while(score == _score(e, pcap, idx) + self->gev[0]);
		}
	}
	return(NULL);
}

/* short, exact matching sequences */
unittest( .name = "trace" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, strlen((char const *)dz_unittest_query));
	struct dz_tail_s const *tail = NULL;
	uint8_t *cigar = NULL;

	tail = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"A", 1);
	cigar = dz_trace(dz, tail);

	tail = dz_extend(dz, q, &tail, 1, (uint8_t const *)"GC", 2);
	cigar = dz_trace(dz, tail);

	tail = dz_extend(dz, q, &tail, 1, (uint8_t const *)"TTTT", 4);
	cigar = dz_trace(dz, tail);

	dz_destroy(dz);
}






