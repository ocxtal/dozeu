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
struct dz_cap_s { uint32_t spos, epos, _pad[2]; };				/* followed by dz_tail_s */
struct dz_tail_s { uint32_t spos, epos; uint64_t max, inc; struct dz_cap_s const *cap; };
dz_static_assert(sizeof(struct dz_swgv_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_cap_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_tail_s) % sizeof(__m128i) == 0);

/* context (constants and working buffers) */
struct dz_mem_block_s { struct dz_mem_block_s *next; size_t size; };
struct dz_stack_s { struct dz_mem_block_s *curr; uint8_t *top, *end; uint64_t _pad[3]; };
struct dz_mem_s { struct dz_mem_block_s blk; struct dz_stack_s stack; };
#define dz_mem_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )

struct dz_s { int8_t matrix[16]; uint16_t giv[8], gev[8], xt, _pad[11]; struct dz_tail_s const *root; };
dz_static_assert(sizeof(struct dz_s) % sizeof(__m128i) == 0);
#define dz_mem(_self)				( (struct dz_mem_s *)(_self) - 1 )
#define dz_root(_self)				( (struct dz_tail_s const *const *)(&_self->root) )

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
#define _alloc_column(_len) ({ \
	size_t next_req = 2 * (size_t)(_len) * sizeof(struct dz_swgv_s) + sizeof(struct dz_cap_s); \
	if(dz_mem_stack_rem(dz_mem(self)) < next_req) { dz_mem_add_stack(dz_mem(self), 0); } \
	uintptr_t _p = (uintptr_t)dz_mem(self)->stack.top; \
	(struct dz_swgv_s *)dz_roundup(_p, sizeof(__m128i)); \
})
#define _load_vector(_p) \
	__m128i e = _mm_load_si128((__m128i const *)(&((struct dz_swgv_s const *)(_p))->e)); print_vector(e); \
	__m128i s = _mm_load_si128((__m128i const *)(&((struct dz_swgv_s const *)(_p))->s)); print_vector(s);
#define _update_vector(_p) { \
	__m128i rv = _mm_set1_epi8(conv[*rp & 0x0f]), qv = _mm_loadl_epi64((__m128i const *)&qp[(_p) * L]); \
	__m128i sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)self->matrix), _mm_or_si128(rv, qv))); print_vector(sc); \
	__m128i te = _mm_subs_epi16(_mm_max_epi16(e, _mm_subs_epi16(s, giv)), gev1); \
	 print_vector(_mm_alignr_epi8(s, ps, 14));  \
	__m128i ts = _mm_max_epi16(te, _mm_adds_epi16(sc, _mm_alignr_epi8(s, ps, 14))); ps = s; \
	__m128i tf = _mm_max_epi16(_mm_subs_epi16(ts, giv), _mm_subs_epi16(_mm_alignr_epi8(minv, f, 14), gev1)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 14), gev1)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 12), gev2)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 8), gev4)); \
	ts = _mm_max_epi16(ts, tf); \
	mv = _mm_max_epi16(mv, ts); \
	/* print_vector(te); */ print_vector(ts); /* print_vector(tf); print_vector(mv); */ \
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
#define _push_cap(_p) ({ \
	struct dz_cap_s *cap = (struct dz_cap_s *)&(_p)[w.epos]; \
	cap->spos = w.spos; cap->epos = w.epos; \
	dz_mem(self)->stack.top = (uint8_t *)(cap + 1); \
	cap; \
})

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

	/* calc vector length, malloc the first column */
	max_gap_len = dz_roundup(max_gap_len, L);
	struct dz_swgv_s *dp = _alloc_column(max_gap_len / L);

	/* cap object is placed at the tail of the last column */
	struct dz_tail_s w = { .epos = max_gap_len / L };

	/* fill the root (the leftmost) column */
	__m128i s = _mm_setr_epi16(0, -(gi+ge), -(gi+2*ge), -(gi+3*ge), -(gi+4*ge), -(gi+5*ge), -(gi+6*ge), -(gi+7*ge));
	__m128i const e = _mm_set1_epi16(DZ_CELL_MIN), xtv = _mm_set1_epi16(-self->xt);
	for(uint64_t p = 0; p < max_gap_len / L; p++) {
		__m128i const f = s;
		if(dz_unlikely(_test_xdrop(s, xtv))) { debug("p(%lu)", p); w.epos = p; break; }
		_store_vector(&dp[p]);
		if(p == 0) { s = _mm_setr_epi16(-gi, -(gi+ge), -(gi+2*ge), -(gi+3*ge), -(gi+4*ge), -(gi+5*ge), -(gi+6*ge), -(gi+7*ge)); }
		s = _mm_subs_epi16(s, _mm_slli_epi16(gev, 3));
	}

	/* done */
	struct dz_tail_s *tail = (struct dz_tail_s *)(w.cap = (struct dz_cap_s const *)&dp[w.epos]);
	self->root = tail; *tail = w;
	dz_mem(self)->stack.top = (uint8_t *)(tail + 1);
	debug("[%u, %u), inc(%lu), cap(%p), max(%lu)", w.spos, w.epos, w.inc, w.cap, w.max);
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

	/* ASCII to 2-bit conversion */
	static uint8_t const conv[16] __attribute__(( aligned(16) )) = {
		['A' & 0x0f] = A<<2, ['C' & 0x0f] = C<<2,
		['G' & 0x0f] = G<<2, ['T' & 0x0f] = T<<2,
		['U' & 0x0f] = T<<2, ['N' & 0x0f] = N<<2
	};
	__m128i pv = _mm_setzero_si128();
	__m128i const cv = _mm_load_si128((__m128i const *)conv), fv = _mm_set1_epi8(0x0f);			/* conversion table and mask */
	for(uint64_t i = 0; i < dz_rounddown(qlen, sizeof(__m128i)); i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
		__m128i tv = _mm_shuffle_epi8(cv, _mm_and_si128(qv, fv));
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;
	}
	for(uint64_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[query[i]];
	}
	return(q);
}

unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	uint8_t const *query = (uint8_t const *)
		"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
		"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
		"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
		"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
		"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
		"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
		"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
		"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

	struct dz_query_s *q = dz_pack_query(dz, query, strlen((char const *)query));
	dz_unused(q);
	dz_destroy(dz);
}

/**
 * @fn dz_extend
 */
static
struct dz_tail_s *dz_extend(
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_tail_s const *const *edges, uint64_t n_edges,
	uint8_t const *ref, int64_t rlen)
{
	uint64_t const L = sizeof(__m128i) / sizeof(uint16_t);
	if(n_edges == 0 || rlen == 0) { return(NULL); }

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
	uint64_t adj[n_edges];								/* S[0, 0] = 0 */
	struct dz_tail_s w = { .spos = UINT32_MAX };
	for(uint64_t i = 0; i < n_edges; i++) {
		/* update max and pos */
		w.max = dz_max2(w.max, edges[i]->max);
		w.spos = dz_min2(w.spos, edges[i]->spos);
		w.epos = dz_max2(w.epos, edges[i]->epos);
		debug("i(%lu), [%u, %u), inc(%lu), max(%lu)", i, edges[i]->spos, edges[i]->epos, edges[i]->inc, edges[i]->max);
	}
	debug("start extension [%u, %u), inc(%lu), max(%lu), stack(%p, %p)", w.spos, w.epos, w.inc, w.max, dz_mem(self)->stack.top, dz_mem(self)->stack.end);

	for(uint64_t i = 0; i < n_edges; i++) {
		adj[i] = w.max - (edges[i]->max - edges[i]->inc);		/* base = max - inc */
		debug("i(%lu), adj(%lu)", i, adj[i]);
	}

	/* X-drop DP */
	ptrdiff_t dir = rlen < 0 ? -1 : 1;					/* direction */
	uint8_t const *rp = rlen < 0 ? &ref[rlen - 1] : ref;
	struct dz_swgv_s *pdp = _alloc_column(w.epos - w.spos) - w.spos;	/* allocate memory for the first column */
	for(uint64_t p = w.spos; p < w.epos; p++) {
		/* memset(pdp, 0xff, sizeof(struct dz_swgv_s) * (w.epos - w.spos)); */
		__m128i const e = _mm_set1_epi16(INT16_MIN), f = e, s = e; _store_vector(&pdp[p]);
	}

	/* paste the last vectors */
	for(uint64_t i = 0; i < n_edges; i++) {
		struct dz_swgv_s const *tdp = (struct dz_swgv_s const *)edges[i] - edges[i]->epos;
		__m128i const adjv = _mm_set1_epi16(adj[i]);
		for(uint64_t p = edges[i]->spos; p < edges[i]->epos; p++) {
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

	/* once */ {
		/* extend the first vector (clip the column length if too long) */
		__m128i f = minv, ps = minv, mv = _mm_setzero_si128();
		__m128i const xtv = _mm_set1_epi16(-self->xt);					/* next offset == current max thus X-drop threshold is always -xt */
		for(uint64_t p = w.spos; p < w.epos; p++) {
			_load_vector(&pdp[p]); _update_vector(p);
			if(dz_unlikely(_test_xdrop(s, xtv))) {	/* tag _unlikely to move out of the core loop */
				/* drop either leading or trailing vector, skip the tail extension when the tail is clipped */
				if(p == w.spos) { w.spos++; } else { w.epos = p; goto _merge_tail; }
			}
			_store_vector(&pdp[p]);
		}
		/* tail extension */
		if(w.epos >= qblen) { goto _merge_tail; }
		__m128i e = minv, s = minv; _update_vector(w.epos);
		do {
			if(_test_xdrop(s, xtv)) { break; }
			_store_vector(&pdp[w.epos]); w.epos++;
			f = _mm_subs_epi16(f, gev8); s = _mm_subs_epi16(s, gev8);
		} while(w.epos < qblen);
	_merge_tail:;
		struct dz_cap_s *cap = _push_cap(pdp);
		uint64_t inc = _hmax_vector(mv);
		if(inc > w.inc) { w.inc = inc; w.cap = cap; }					/* 2 x cmov */
		debug("update inc(%lu), max(%lu, %lu, %p), cap(%p)", w.inc, w.max, w.max + w.inc, w.cap, cap);
	}

	/* extension loop */
	debug("rp(%p), dir(%ld)", rp, dir);
	while((rp += dir, --rlen > 0)) {
		debug("rp(%p), [%u, %u), stack(%p, %p)", rp, w.spos, w.epos, dz_mem(self)->stack.top, dz_mem(self)->stack.end);
		struct dz_swgv_s *dp = _alloc_column(w.epos - w.spos) - w.spos;
		__m128i f = minv, ps = minv, mv = _mm_setzero_si128();
		__m128i const xtv = _mm_set1_epi16(w.inc - self->xt);
		for(uint64_t p = w.spos; p < w.epos; p++) {
			_load_vector(&pdp[p]); _update_vector(p);
			if(dz_unlikely(_test_xdrop(s, xtv))) {
				if(p == w.spos) { w.spos++; } else { w.epos = p + 1; goto _loop_tail; }
			}
			_store_vector(&dp[p]);
		}
		if(w.epos >= qblen) { goto _loop_tail; }
		__m128i e = minv, s = minv; _update_vector(w.epos);
		do {
			if(_test_xdrop(s, xtv)) { break; }
			_store_vector(&dp[w.epos]); w.epos++;
			f = _mm_subs_epi16(f, gev8); s = _mm_subs_epi16(s, gev8);
		} while(w.epos < qblen);
	_loop_tail:;
		struct dz_cap_s *cap = _push_cap(dp); pdp = dp;
		uint64_t inc = _hmax_vector(mv);
		if(inc > w.inc) { w.inc = inc; w.cap = cap; }					/* 2 x cmov */
		debug("update inc(%lu), max(%lu, %lu, %p), cap(%p)", w.inc, w.max, w.max + w.inc, w.cap, cap);
		/* FIXME: rescue overflow */
	}

	/* update max vector (pdp always points at the last vector) */
	struct dz_tail_s *tail = (struct dz_tail_s *)((struct dz_cap_s *)dz_mem(self)->stack.top - 1);
	tail->max = w.max + w.inc;
	tail->inc = w.inc;
	tail->cap = w.cap;
	dz_mem(self)->stack.top = (uint8_t *)(tail + 1);
	debug("stack(%p, %p), tail(%p), max(%lu), inc(%lu), cap(%p)", dz_mem(self)->stack.top, dz_mem(self)->stack.end, tail, tail->max, tail->inc, tail->cap);
	return(tail);
}

/* short, exact matching sequences */
unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	uint8_t const *query = (uint8_t const *)
		"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
		"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
		"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
		"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
		"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
		"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
		"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
		"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

	struct dz_query_s *q = dz_pack_query(dz, query, strlen((char const *)query));
	struct dz_tail_s *cap = NULL;

	/* nothing occurs */
	cap = dz_extend(dz, q, NULL, 0, NULL, 0);
	ut_assert(cap == NULL);
	cap = dz_extend(dz, q, NULL, 0, (uint8_t const *)"", 0);
	ut_assert(cap == NULL);
	cap = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"", 0);
	ut_assert(cap == NULL);

	/* extend */
	cap = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"A", 1);
	ut_assert(cap != NULL && cap->max == 2);

	cap = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AG", 2);
	ut_assert(cap != NULL && cap->max == 4);

	cap = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AGATTTT", 7);
	ut_assert(cap != NULL && cap->max == 9);

	dz_destroy(dz);
}

/* a small graph */
unittest( .name = "small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	uint8_t const *query = (uint8_t const *)
		"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
		"TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
		"TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
		"ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
		"CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
		"GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
		"AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
		"AAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT";

	struct dz_query_s *q = dz_pack_query(dz, query, strlen((char const *)query));
	struct dz_tail_s const *cap[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	cap[0] = dz_extend(dz, q, dz_root(dz), 1, (uint8_t const *)"AG", 2);
	ut_assert(cap[0] != NULL && cap[0]->max == 4);

	cap[1] = dz_extend(dz, q, &cap[0], 1, (uint8_t const *)"C", 1);
	ut_assert(cap[1] != NULL && cap[1]->max == 6);

	cap[2] = dz_extend(dz, q, &cap[0], 2, (uint8_t const *)"TTTT", 4);
	ut_assert(cap[2] != NULL && cap[2]->max == 14);

	cap[3] = dz_extend(dz, q, &cap[2], 1, (uint8_t const *)"CATT", 4);
	ut_assert(cap[3] != NULL && cap[3]->max == 22);

	cap[4] = dz_extend(dz, q, &cap[2], 2, (uint8_t const *)"CTGA", 4);
	ut_assert(cap[4] != NULL && cap[4]->max == 30);

	dz_destroy(dz);
}

/**
 * @fn dz_trace
 */
static
uint8_t *dz_trace(
	struct dz_mem_s *self,
	struct dz_tail_s *edges)
{


	return(NULL);
}







