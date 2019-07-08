
// $(CC) -O3 -march=native -DMAIN -o dozeu dozeu.c
// #define DEBUG
// #define DZ_PRINT_VECTOR

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
#define print_vector_8(v) ;

#endif


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

// #include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <x86intrin.h>



/* default encoding: ascii */
#if !defined(DZ_WRAPPED_API) || (DZ_WRAPPED_API != 0)
#  if !defined(DZ_PROTEIN) && !defined(DZ_NUCL_ASCII) && !defined(DZ_NUCL_2BIT) && !defined(DZ_NUCL_4BIT) && !defined(DZ_REALIGNER)
#    define DZ_NUCL_ASCII
#  endif
#endif



/* scoring and traceback configurations */
#ifndef DZ_CIGAR_OP
#  define DZ_CIGAR_OP				0x03040201
#endif
#define DZ_CIGAR_INTL				( (((uint32_t)DZ_CIGAR_OP)>>16) | (((uint32_t)DZ_CIGAR_OP)<<16) )

#ifndef dz_cmp_max
#  define dz_cmp_max(x, y)			( (x) > (y) )
#endif

#ifndef DZ_OVERFLOW_MARGIN
#  define DZ_OVERFLOW_MARGIN		( 256 )
#endif



/* activate wrapped API (enabled by default) */
#ifndef DZ_WRAPPED_API
#  define DZ_WRAPPED_API			( 1 )
#endif



/* vectorization directive */
#ifndef __x86_64__
#  error "x86_64 is required"
#endif
#ifndef __SSE4_1__
#  warning "SSE4.1 is automatically enabled in dozeu.h, please check compatibility to the system."
#  define __dz_vectorize			__attribute__(( target( "sse4.1" ) )) inline
#else
#  define __dz_vectorize			/* follow the compiler options */
#endif



/* inlining directive (FIXME: add appropriate __force_inline flag) */
#define __dz_force_inline			inline



/* warning suppressor */
#ifdef __cplusplus
#  define dz_unused(x)
#else
#  define dz_unused(x)				(void)(x)
#endif



/* unittest and debug configuration */
#if (defined(DEBUG) || (defined(UNITTEST) && UNITTEST != 0)) && !defined(__cplusplus)
#  include "debug.h"
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



/*
 * scoring configuration
 */
typedef struct dz_score_conf_s {
	/*
	 * dimension is determined by macro at compile time;
	 * DZ_REF_MAT_SIZE for row length and DZ_QUERY_MAT_SIZE for #columns
	 */
	int8_t const *score_matrix;		/* negative for mismatch */

	/* gap penalties */
	uint16_t ins_open;				/* gap penalties in positive */
	uint16_t ins_extend;
	uint16_t del_open;				/* gap penalties in positive */
	uint16_t del_extend;

	/* end controls */
	uint16_t scan_ref;				/* do not penalize head gap on reference if scan_ref != 0 */
	uint16_t full_length_bonus;		/* end-to-end mapping bonus; only activated when compiled with -DDZ_FULL_LENGTH_BONUS */

	/* X-drop threshold */
	size_t max_ins_len;				/* as X-drop threshold */
	size_t max_del_len;				/* as X-drop threshold */
} dz_score_conf_t;

/* external (custom) allocator for profile constructor */
typedef void *(*dz_malloc_t)(void *ctx, size_t size);
typedef void (*dz_free_t)(void *ctx, void *ptr);

typedef struct dz_allocator_s {
	void *ctx;
	dz_malloc_t fp;
} dz_allocator_t;

typedef struct dz_destructor_s {
	void *ctx;
	dz_free_t fp;
} dz_destructor_t;


/*
 * type declarations
 */
typedef struct dz_arena_s dz_arena_t;
typedef struct dz_profile_s dz_profile_t;
typedef struct dz_forefront_s dz_forefront_t;
typedef struct dz_state_s dz_state_t;
typedef struct dz_alignment_s dz_alignment_t;

/* pointer casting */
#define dz_state(_p)				( (dz_state_t *)(_p) )
#define dz_cstate(_p)				( (dz_state_t const *)(_p) )



#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)

/* dz_t: wrapper context */
typedef struct dz_s {
	dz_arena_t *mem;
	dz_profile_t *profile;
} dz_t;

#endif



/*
 * interface abstraction:
 * abstract query packer
 *
 * score_matrix_row == &profile->matrix[ridx * DZ_QUERY_MAT_SIZE]
 * dir == 0 for forward, 1 for reverse
 */
typedef size_t (*dz_query_calc_dim_t)(size_t qlen);						/* ref-side size */
typedef __m128i (*dz_query_conv_t)(int8_t const *score_matrix_row, uint32_t dir, __m128i v);

typedef struct {
	uint32_t dir;
	uint8_t invalid;						/* invalid base */

	/* we do not provide opaque pointer because we suppose all the information we need is contained in profile object */
	dz_query_calc_dim_t calc_dim;
	dz_query_conv_t conv;
} dz_pack_query_t;

typedef struct dz_query_s dz_query_t;		/* forward declaration */

static size_t dz_query_row_len(dz_query_t const *query);				/* query->blen * DZ_L */
static uint8_t const *dz_query_seq(dz_query_t const *query);				/* query->q */
static uint8_t const *dz_query_packed_array(dz_query_t const *query);	/* query->arr */

/*
 * core pack function forward declaration
 */
static dz_query_t *dz_pack_query_core(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack);



/*
 * interface abstraction:
 * abstract sequence fetcher
 */
typedef struct {
	uint32_t is_term;
	uint32_t ch;
} dz_fill_fetch_t;

/* exact banding for realigner */
typedef struct {
	uint32_t sidx, eidx;
} dz_fill_bound_t;

typedef struct {
	uint32_t match;
	int32_t score;
} dz_trace_match_t;

typedef dz_fill_fetch_t (*dz_fill_fetch_next_t)(void *opaque, int8_t const *score_matrix, dz_query_t const *query);
typedef __m128i (*dz_fill_get_profile_t)(void *opaque, int8_t const *score_matrix, dz_query_t const *query, size_t qidx);
typedef dz_fill_bound_t (*dz_fill_get_bound_t)(void *opaque);
typedef dz_trace_match_t (*dz_trace_get_match_t)(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch);

typedef struct {
	void *opaque;
	dz_fill_fetch_next_t fetch_next;
	dz_fill_get_profile_t get_profile;
	dz_fill_get_bound_t get_bound;		/* can be NULL */
} dz_fetcher_t;

/*
 * core fill / trace function forward declarations
 */
static dz_state_t const *dz_extend_core(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_fetcher_t *fetcher, dz_state_t const **ff, size_t fcnt);
static dz_alignment_t const *dz_trace_core(dz_arena_t *mem, dz_profile_t const *profile, dz_state_t const *ff, dz_trace_get_match_t get_match);




/* implementations */

#if defined(DZ_NUCL_ASCII)
/*
 * matrix size is 4 x 4
 * query sequence is converted to 2bit representation in pack_query
 */
#define DZ_REF_MAT_SIZE			( 4 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

static
size_t dz_query_calc_dim_ascii(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_ascii(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	/*
	 * ASCII to 2-bit conversion table (shifted by two bits to be used as the upper (row) shuffle index)
	 * @ABC_DEFG_HIJK_LMNO
	 * PQRS_TUVW_XYZ
	 *
	 * ['A' & 0x0f] = A
	 * ['C' & 0x0f] = C
	 * ['G' & 0x0f] = G
	 * ['T' & 0x0f] = T
	 * ['U' & 0x0f] = T
	 * ['N' & 0x0f] = N
	 */
	static uint8_t const table[16 * 2] __attribute__(( aligned(16) )) = {
		0, qA, 0, qC, qT, qU, 0, qG, 0, 0, 0, 0, 0, 0, qN, 0,
		0, qT, 0, qG, qA, qA, 0, qC, 0, 0, 0, 0, 0, 0, qN, 0
	};
	__m128i const cv = _mm_load_si128((__m128i const *)&table[16 * dir]);
	// print_vector_8(v);
	// print_vector_8(cv);
	// print_vector_8(_mm_shuffle_epi8(cv, v));
	return(_mm_shuffle_epi8(cv, v));
}

static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_ascii,
		.conv     = dz_query_conv_ascii
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}


typedef struct {
	__m128i rv;
	uint8_t conv[16];
	uint8_t const *p, *t;
	int64_t inc;
} dz_fetcher_ascii_t;

static __dz_vectorize
void dz_fill_fetcher_init_ascii(dz_fetcher_ascii_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	debug("ref(%p), rlen(%zu), p(%p), t(%p), dir(%u), inc(%ld)", ref, rlen, self->p, self->t, dir, self->inc);

	static uint8_t const conv[16 * 2] __attribute__(( aligned(16) )) = {
		/* ['A' & 0x0f] = A, ['C' & 0x0f] = C, ['G' & 0x0f] = G, ['T' & 0x0f] = T, ['U' & 0x0f] = T, ['N' & 0x0f] = N */
		0, rA, 0, rC, rT, rU, 0, rG, 0, 0, 0, 0, 0, 0, rN, 0,
		0, rT, 0, rG, rA, rA, 0, rC, 0, 0, 0, 0, 0, 0, rN, 0
	};

	/* create conversion table (accessed by indirect %rbp) */
	__m128i const cv = _mm_load_si128((__m128i const *)&conv[16 * dir]);
	_mm_store_si128((__m128i *)self->conv, cv);
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fill_fetch_next_ascii(dz_fetcher_ascii_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = self->conv[c & 0x0f];
	self->rv = _mm_set1_epi8(e);
	debug("c(%c), e(%x)", c, e);
	print_vector_8(self->rv);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = e
	});
}

static __dz_vectorize
__m128i dz_fill_get_profile_ascii(dz_fetcher_ascii_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const mv = _mm_load_si128((__m128i const *)score_matrix);
	__m128i const sc = _mm_shuffle_epi8(mv, _mm_or_si128(self->rv, qv));

	print_vector_8(qv);
	// print_vector_8(mv);
	// print_vector_8(sc);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_trace_get_match_ascii(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[(packed[qidx] | ch) & 0x1f],
		.match = packed[qidx] == (ch<<2)
	});
}

#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_ascii_t w __attribute__(( aligned(16) ));
	dz_fill_fetcher_init_ascii(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)dz_fill_fetch_next_ascii,
		.get_profile = (dz_fill_get_profile_t)dz_fill_get_profile_ascii,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace(dz_t *self, dz_forefront_t const *forefront)
{
	/* just wrap */
	return(dz_trace_core(self->mem, self->profile, dz_cstate(forefront), (dz_trace_get_match_t)dz_trace_get_match_ascii));
}
#endif




#elif defined(DZ_NUCL_2BIT)
/*
 * 4 x 4 matrix for 2bit encoding
 */
#define DZ_REF_MAT_SIZE			( 4 )
#define DZ_QUERY_MAT_SIZE		( 4 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

enum dz_alphabet {
	rA = 0x00, rC = 0x01, rG = 0x02, rT = 0x03, rU = 0x03,
	qA = 0x00, qC = 0x04, qG = 0x08, qT = 0x0c, qU = 0x0c,
	#ifdef DZ_N_AS_UNMATCHING_BASE
		rN = 0x04, qN = 0x02, qX = 0x02
	#else
		rN = 0x90, qN = 0x90, qX = 0x02		/* pshufb instruction clears the column when the 7-th bit of the index is set */
	#endif
};

static
size_t dz_query_calc_dim_2bit(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_2bit(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const table[16 * 2] __attribute__(( aligned(16) )) = {
		qA, qC, qG, qT, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN,		/* shift up by 2bits */
		qT, qG, qC, qA, qN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qN
	};
	__m128i const cv = _mm_load_si128((__m128i const *)&table[16 * dir]);

	return(_mm_shuffle_epi8(cv, v));
}

static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_2bit,
		.conv     = dz_query_conv_2bit
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}


typedef struct {
	__m128i rv;
	uint8_t const *p, *t;
	int64_t inc;
} dz_fetcher_2bit_t;

static __dz_vectorize
void dz_fill_fetcher_init_2bit(dz_fetcher_2bit_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fill_fetch_next_2bit(dz_fetcher_2bit_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = ((self->inc>>1) ^ c) & 0x03;
	self->rv = _mm_set1_epi8(e);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = e
	});
}

static __dz_vectorize
__m128i dz_fill_get_profile_2bit(dz_fetcher_2bit_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const mv = _mm_load_si128((__m128i const *)score_matrix);
	__m128i const sc = _mm_shuffle_epi8(mv, _mm_or_si128(self->rv, qv));
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_trace_get_match_2bit(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[(packed[qidx] | ch) & 0x1f],
		.match = packed[qidx] == (ch<<2)
	});
}

#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_2bit_t w __attribute__(( aligned(16) ));
	dz_fill_fetcher_init_2bit(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)dz_fill_fetch_next_2bit,
		.get_profile = (dz_fill_get_profile_t)dz_fill_get_profile_2bit,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace(dz_t *self, dz_forefront_t const *forefront)
{
	/* just wrap */
	return(dz_trace_core(self->mem, self->profile, dz_cstate(forefront), (dz_trace_get_match_t)dz_trace_get_match_2bit));
}
#endif




#elif defined(DZ_NUCL_4BIT)
/*
 * matrix size is 16 x 16
 */
#define DZ_REF_MAT_SIZE			( 16 )
#define DZ_QUERY_MAT_SIZE		( 16 )
#define DZ_TRANSPOSE_MATRIX		( 1 )		/* do not transpose */

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

static
size_t dz_query_calc_dim_4bit(size_t qlen)
{
	dz_unused(qlen);

	return(DZ_REF_MAT_SIZE);
}

static __dz_vectorize
__m128i dz_query_conv_4bit(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	if(dir == 0) { return(v); }		/* do nothing for forward, just copy */

	static uint8_t const comp[16] __attribute__(( aligned(16) )) = {
		0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e,
		0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b, 0x07, 0x0f
	};
	__m128i const cv = _mm_load_si128((__m128i const *)comp);
	__m128i const w = _mm_shuffle_epi8(cv, v);
	return(w);			/* complemented */
}

static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_4bit,
		.conv     = dz_query_conv_4bit
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}


typedef struct {
	__m128i mv;
	uint8_t const *p, *t;
	int64_t inc;
	uint64_t conv;
} dz_fetcher_4bit_t;

static __dz_vectorize
void dz_fill_fetcher_init_4bit(dz_fetcher_4bit_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	self->conv = dir ? 0xf7b3d591e6a2c480 : 0xfedcba9876543210;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fill_fetch_next_4bit(dz_fetcher_4bit_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(query);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const c = *self->p;
	uint32_t const e = (self->conv>>(4 * c)) & 0x0f;
	self->mv = _mm_load_si128((__m128i const *)&score_matrix[e * DZ_QUERY_MAT_SIZE]);

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = e
	});
}

static __dz_vectorize
__m128i dz_fill_get_profile_4bit(dz_fetcher_4bit_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);

	__m128i const qv = _mm_loadl_epi64((__m128i const *)&packed[qidx]);
	__m128i const sc = _mm_shuffle_epi8(self->mv, qv);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_trace_get_match_4bit(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = score_matrix[ch * DZ_QUERY_MAT_SIZE + packed[qidx]],
		.match = packed[qidx] == ch
	});
}

#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_4bit_t w __attribute__(( aligned(16) ));
	dz_fill_fetcher_init_4bit(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)dz_fill_fetch_next_4bit,
		.get_profile = (dz_fill_get_profile_t)dz_fill_get_profile_4bit,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace(dz_t *self, dz_forefront_t const *forefront)
{
	/* just wrap */
	return(dz_trace_core(self->mem, self->profile, dz_cstate(forefront), (dz_trace_get_match_t)dz_trace_get_match_4bit));
}
#endif




#elif defined(DZ_PROTEIN)
/*
 * matrix size is 32 x 32
 */
#define DZ_REF_MAT_SIZE			( 32 )
#define DZ_QUERY_MAT_SIZE		( 32 )
#define DZ_TRANSPOSE_MATRIX		( 1 )		/* do not transpose */

static
size_t dz_query_calc_dim_protein(size_t qlen)
{
	dz_unused(qlen);

	return(DZ_REF_MAT_SIZE);
}

static __dz_vectorize
__m128i dz_query_conv_protein(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(dir);

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
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_protein,
		.conv     = dz_query_conv_protein
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}


typedef struct {
	uint8_t const *p, *t;
	int64_t inc;
	uint8_t const *parr;
} dz_fetcher_protein_t;

static __dz_vectorize
void dz_fill_fetcher_init_protein(dz_fetcher_protein_t *self, uint8_t const *ref, size_t rlen, uint32_t dir)
{
	self->p = &ref[dir ? rlen - 1ULL : 0];
	self->t = &ref[dir ? (size_t)-1LL : rlen];
	self->inc = dir ? -1LL : 1LL;
	self->parr = NULL;
	return;
}

static __dz_vectorize
dz_fill_fetch_t dz_fill_fetch_next_protein(dz_fetcher_protein_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);

	if(self->p == self->t) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch base and convert to 2bit */
	uint32_t const e = *self->p & 0x1f;
	uint8_t const *packed = dz_query_packed_array(query);
	self->parr = &packed[e * dz_query_row_len(query)];

	/* forward pointer */
	self->p += self->inc;
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = e
	});
}

static __dz_vectorize
__m128i dz_fill_get_profile_protein(dz_fetcher_protein_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);
	dz_unused(query);

	__m128i const sc = _mm_loadl_epi64((__m128i const *)&self->parr[qidx]);
	return(_mm_cvtepi8_epi16(sc));
}

static __dz_vectorize
dz_trace_match_t dz_trace_get_match_protein(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	dz_unused(score_matrix);

	uint8_t const *seq = dz_query_seq(query);
	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score =  packed[ch * dz_query_row_len(query) + qidx],
		.match = (seq[qidx - 1] & 0x1f) == ch
	});
}

#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_protein_t w __attribute__(( aligned(16) ));
	dz_fill_fetcher_init_protein(&w, (uint8_t const *)ref, rlen, dir);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)dz_fill_fetch_next_protein,
		.get_profile = (dz_fill_get_profile_t)dz_fill_get_profile_protein,
		.get_bound   = NULL
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace(dz_t *self, dz_forefront_t const *forefront)
{
	/* just wrap */
	return(dz_trace_core(self->mem, self->profile, dz_cstate(forefront), (dz_trace_get_match_t)dz_trace_get_match_protein));
}
#endif




#elif defined(DZ_REALIGNER)

#define DZ_SCALING_INDEX		( 6 )
#define DZ_TRANSPOSE_MATRIX		( 0 )		/* do not transpose */

/*
 * we expect the input sequences are 2bit encoded
 */
static
size_t dz_query_calc_dim_realigner(size_t qlen)
{
	dz_unused(qlen);

	return(1);
}

static __dz_vectorize
__m128i dz_query_conv_realigner(int8_t const *score_matrix, uint32_t dir, __m128i v)
{
	dz_unused(score_matrix);

	static uint8_t const mask[16] __attribute__(( aligned(16) )) = {
		0x00, 0x02, 0x04, 0x06, 0x08, 0x08, 0x08, 0x08,
		0x06, 0x04, 0x02, 0x00, 0x08, 0x08, 0x08, 0x08
	};

	/* take complement and double */
	__m128i const mv = _mm_load_si128(mask);
	__m128i const cv = _mm_set1_epi8(dir ? 0x08 : 0x00);
	__m128i const x  = _mm_shuffle_epi8(mv, _mm_or_si128(cv, v));
	return(x);
}

static __dz_vectorize
dz_query_t *dz_pack_query_intl(dz_t *self, char const *query, size_t qlen, uint32_t dir)
{
	dz_pack_query_t pack = {
		.dir      = dir,
		.calc_dim = dz_query_calc_dim_realigner,
		.conv     = dz_query_conv_realigner
	};
	return(dz_pack_query_core(self->mem, self->profile, (uint8_t const *)query, qlen, &pack));
}


typedef struct dz_realigner_node_s {
	uint32_t a, c, g, t, gap;
	uint32_t consensus_mask;

	/* range on query */
	dz_fill_bound_t bound;

	/* links */
	struct dz_realigner_node_s *next, *prev;
} dz_realigner_node_t;

typedef struct {
	__m128i sv;
	dz_realigner_node_t const *p;
	uint32_t dir;
} dz_fetcher_realigner_t;

static __dz_vectorize
void dz_fill_fetcher_init_realigner(dz_fetcher_realigner_t *self, uint32_t dir, dz_realigner_node_t const *ref)
{
	self->p   = ref;
	self->dir = dir;
	return;
}


static __dz_vectorize
__m128i dz_fill_calc_delta_c(dz_realigner_node_t const *node)
{
	static uint8_t const mask[16] __attribute__(( aligned(16) )) = {
		0x01, 0x01, 0x02, 0x02, 0x04, 0x04, 0x08, 0x08,
		0x10, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
	};
	__m128i const m = _mm_load_si128((__m128i const *)mask);

	__m128i const z = _mm_setzero_si128();
	__m128i const v = _mm_set1_epi8(node->consensus_mask);
	__m128i const w = _mm_cmpeq_epi8(_mm_and_si128(v, m), z);
	__m128i const x = _mm_sub_epi16(z, w);
	__m128i const y = _mm_slli_epi16(x, DZ_SCALING_INDEX);
	return(y);
}

static __dz_vectorize
__m128i dz_fill_calc_n(__m128i s, __m128i t)
{
	__m128i const x = _mm_add_epi32(s, _mm_shuffle_epi32(s, 0x4e));
	__m128i const y = _mm_add_epi32(x, _mm_shuffle_epi32(x, 0xb1));
	__m128i const z = _mm_add_epi32(y, _mm_shuffle_epi32(t, 0x00));
	return(z);
}

static __dz_vectorize
__m128i dz_fill_calc_delta_a(dz_realigner_node_t const *node)
{
	/* it's possible to use ymm here for better performance */
	__m128i const x1 = _mm_loadu_si128((__m128i const *)&node->a);		/* (t, g, c, a) */
	__m128i const x2 = _mm_loadu_si128((__m128i const *)&node->gap);	/* (-, -, -, gap) */
	__m128i const n = dz_fill_calc_n(x1, x2);							/* a + c + g + t + gap */

	__m128i const y1 = _mm_slli_epi32(_mm_sub_epi32(n, x1), DZ_SCALING_INDEX);	/* (..., a + g + t + gap, c + g + t + gap)<<k */
	__m128i const y2 = _mm_slli_epi32(_mm_sub_epi32(n, x2), DZ_SCALING_INDEX);

	__m128 const m = _mm_cvtepi32_ps(n);
	__m128 const z1 = _mm_div_ps(_mm_cvtepi32_ps(y1), m);
	__m128 const z2 = _mm_div_ps(_mm_cvtepi32_ps(y2), m);

	__m128i const w1 = _mm_cvtps_epi32(z1);
	__m128i const w2 = _mm_cvtps_epi32(z2);
	return(_mm_packs_epi32(w1, w2));		/* convert to 16bit */
}

static __dz_vectorize
dz_fill_fetch_t dz_fill_fetch_next_realigner(dz_fetcher_realigner_t *self, int8_t const *score_matrix, dz_query_t const *query)
{
	dz_unused(score_matrix);

	if(self->p == NULL) {
		return((dz_fill_fetch_t){
			.is_term = 1,
			.ch      = 0
		});
	}

	/* fetch profile */
	dz_realigner_node_t const *node = self->p;

	/* compute deltas; both scaled */
	__m128i const dc = dz_fill_calc_delta_c(node);
	__m128i const da = dz_fill_calc_delta_a(node);
	self->sv = _mm_add_epi16(dc, da);

	/* forward pointer */
	dz_realigner_node_t const *link = &node->next;
	self->p = link[self->dir];
	return((dz_fill_fetch_t){
		.is_term = 0,
		.ch      = node->consensus_mask
	});
}

static __dz_vectorize
__m128i dz_fill_get_profile_realigner(dz_fetcher_realigner_t *self, int8_t const *score_matrix, dz_query_t const *query, size_t qidx)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);
	__m128i const v = _mm_loadu_si128((__m128i const *)&packed[qidx]);

	/* convert to (2x + 1, 2x) pairs */
	__m128i const w = _mm_sub_epi8(v, _mm_cmpeq_epi8(self->sv, self->sv));	/* +1 */
	__m128i const x = _mm_unpacklo_epi8(v, w);

	/* shuffle 16bit array */
	__m128i const y = _mm_shuffle_epi8(self->sv, x);
	return(y);
}

static __dz_vectorize
dz_fill_bound_t dz_fill_get_bound_realigner(dz_fetcher_realigner_t *self)
{
	dz_realigner_node_t const *node = self->p;
	return(node->bound);
}

static __dz_vectorize
dz_trace_match_t dz_trace_get_match_realigner(int8_t const *score_matrix, dz_query_t const *query, size_t qidx, uint32_t ch)
{
	dz_unused(score_matrix);

	uint8_t const *packed = dz_query_packed_array(query);

	return((dz_trace_match_t){
		.score = packed[2 * qidx],
		.match = packed[2 * qidx] == 2 * ch
	});
}

#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)
static __dz_vectorize
dz_state_t const *dz_extend_intl(dz_arena_t *mem, dz_profile_t const *profile, dz_query_t const *query, dz_state_t const **ff, size_t fcnt, char const *ref, size_t rlen, uint32_t dir)
{
	dz_fetcher_realigner_t w __attribute__(( aligned(16) ));
	dz_fill_fetcher_init_realigner(&w, dir, (dz_realigner_node_t const *)ref);

	dz_fetcher_t fetcher = {
		.opaque      = (void *)&w,
		.fetch_next  = (dz_fill_fetch_next_t)dz_fill_fetch_next_realigner,
		.get_profile = (dz_fill_get_profile_t)dz_fill_get_profile_realigner,
		.get_bound   = (dz_fill_get_bound_t)dz_fill_get_bound_realigner
	};
	return(dz_extend_core(mem, profile, query, &fetcher, ff, fcnt));
}

static __dz_vectorize
dz_alignment_t const *dz_trace(dz_t *self, dz_forefront_t const *forefront)
{
	/* just wrap */
	return(dz_trace_core(self->mem, self->profile, dz_cstate(forefront), (dz_trace_get_match_t)dz_trace_get_match_realigner));
}
#endif

#endif




/*
 * core implementations:
 * utility macros
 */
#define dz_pp_cat_intl(x, y)		x##y
#define dz_pp_cat(x, y)				dz_pp_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_pp_cat(_st_, __LINE__)[(expr) ? 1 : -1]
#define dz_trap()					{ *((volatile uint8_t *)NULL); }



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



/* control */
#define dz_likely(x)				__builtin_expect(!!(x), 1)
#define dz_unlikely(x)				__builtin_expect(!!(x), 0)
#define dz_roundup(x, base)			( ((x) + (base) - 1) & ~((base) - 1) )
#define dz_rounddown(x, base)		( (x) & ~((base) - 1) )
#define dz_max2(x, y)				( (x) < (y) ? (y) : (x) )
#define dz_min2(x, y)				( (x) < (y) ? (x) : (y) )



/* memory and pointer handling */
#define dz_add_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) + (size_t)(_ofs)) )
#define dz_sub_ptr(_x, _ofs)		( (void *)((uint8_t *)(_x) - (size_t)(_ofs)) )

#define DZ_MEM_MARGIN_SIZE			( 256 )
#define DZ_MEM_ALIGN_SIZE			( 16 )
#define DZ_MEM_INIT_SIZE			( 16 * 1024 * 1024 )

#define dz_inside(x, y, z)			( ((uint64_t)(y) - (uint64_t)(x)) < ((uint64_t)(z) - (uint64_t)(x)) )
#define dz_loadu_u64(p)				({ uint8_t const *_p = (uint8_t const *)(p); *((uint64_t const *)_p); })
#define dz_storeu_u64(p, e)			{ uint8_t *_p = (uint8_t *)(p); *((uint64_t *)(_p)) = (e); }



/* score handling */
#define dz_add_ofs(_x)				( (uint16_t)((_x) + 0x8000) )
#define dz_rm_ofs(_x)				( (int32_t)((_x) - 0x8000) )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_CELL_MARGIN				( 32 )
#define DZ_CELL_MARGINED_MIN		( DZ_CELL_MIN + DZ_CELL_MARGIN )
#define DZ_CELL_MARGINED_MAX		( DZ_CELL_MAX - DZ_CELL_MARGIN )
#define DZ_L						( sizeof(__m128i) / sizeof(uint16_t) )

#define DZ_SCORE_OFS				( 64 )
#define DZ_HEAD_RCH					( 0x40 )
#define DZ_ROOT_RCH					( 0x80 )



/* vectorization supplementary */
#ifdef __SSE4_1__
#  define dz_is_all_zero(x)			( _mm_test_all_zeros((x), (x)) == 1 )
#else
// #  define dz_is_all_zero(x)			( _mm_movemask_epi8((x)) == 0 )
#endif




/*
 * memory allocation (arena):
 */
typedef struct dz_arena_block_s {
	struct dz_arena_block_s *next;
	size_t size;
} dz_arena_block_t;

typedef struct dz_stack_s {
	dz_arena_block_t *curr;
	uint8_t *top, *end;
	uint64_t _pad;
} dz_stack_t;

/* typedef */
struct dz_arena_s {
	dz_arena_block_t blk;
	dz_stack_t stack;
	uint64_t _pad[2];
} /* dz_arena_t */;
dz_static_assert(sizeof(dz_arena_t) == 64);
#define dz_arena_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )


/* aligned and margined malloc-free pair */
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

/* chained stack */
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
	mem->stack.end = (uint8_t *)mem->stack.curr + mem->stack.curr->size;		/* margin already added at the tail by dz_malloc */
	return(0);
}
static __dz_force_inline
void *dz_arena_malloc(dz_arena_t *mem, size_t size)
{
	size = dz_roundup(size, sizeof(__m128i));
	debug("size(%zu), stack(%p, %p, %zu)", size, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem));

	if(size > dz_arena_stack_rem(mem) || dz_arena_stack_rem(mem) < 4096) {
		dz_arena_add_stack(mem, size);
	}
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += size;

	debug("size(%zu), stack(%p, %p, %zu)", size, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem));
	return(ptr);
}


/*
 * save and restore malloc state; NOTE: unstable
 */
typedef struct dz_freeze_s {
	dz_stack_t stack;
} dz_freeze_t;
dz_static_assert(sizeof(dz_freeze_t) == 32);

static __dz_force_inline
dz_freeze_t *dz_arena_freeze(dz_arena_t *mem)
{
	dz_stack_t stack = mem->stack;

	dz_freeze_t *fz = dz_arena_malloc(mem, sizeof(dz_freeze_t));
	fz->stack = stack;
	return(fz);
}
static __dz_force_inline
void dz_arena_restore(dz_arena_t *mem, dz_freeze_t const *fz)
{
	mem->stack = fz->stack;
	return;
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



/*
 * score profile object
 */
/* typedef */
struct dz_profile_s {
	uint16_t iiv[8], iev[8];
	uint16_t div[8], dev[8];	/* deletion: open and extend */

	/* root column */
	struct dz_forefront_s const *root;

	/* constants */
	uint16_t xt, bonus, _pad, init;

	uint32_t _pad2;
	uint32_t max_ins_len, max_del_len;
	uint32_t size;				/* object size for copying */

	int8_t matrix[];
} /* dz_profile_t */;
dz_static_assert(sizeof(dz_profile_t) % sizeof(__m128i) == 0);





/*
 * query packing:
 *
 * blen = roundup(qlen, DZ_L) / DZ_L
 * array must have 16-byte-length margin at the tail
 */
/* typedef */
struct dz_query_s {
	size_t blen;		/* #vectors */
	uint8_t const *seq;	/* query string (encoded and aligned) */
	int16_t bonus[2 * DZ_L];
	uint8_t arr[];
} /* dz_query_t */;
dz_static_assert(sizeof(dz_query_t) % sizeof(__m128i) == 0);

/* full-length bonus */
static __dz_vectorize
__m128i dz_query_get_bonus(dz_query_t const *query, size_t qblk)
{
	#ifdef DZ_FULL_LENGTH_BONUS
		__m128i const *pbonus = &((__m128i const *)query->arr)[-2];
		__m128i const bv = _mm_load_si128(&pbonus[qblk == query->blen - 1]);

	#else
		dz_unused(query);
		dz_unused(qblk);

		__m128i const bv = _mm_setzero_si128();

	#endif
	return(bv);
}

static __dz_vectorize
size_t dz_query_row_len(dz_query_t const *query)
{
	return(query->blen * DZ_L);
}

static __dz_vectorize
uint8_t const *dz_query_seq(dz_query_t const *query)
{
	return(query->seq);
}

static __dz_vectorize
uint8_t const *dz_query_packed_array(dz_query_t const *query)
{
	return(query->arr);
}


/* profile vector is calculated on-the-fly */
static __dz_vectorize
dz_query_t *dz_pack_query_alloc_mem(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const dim  = pack->calc_dim(clen);

	/* allocate mem */
	dz_query_t *q = (dz_query_t *)dz_arena_malloc(mem,
		  sizeof(dz_query_t)	/* header */
		+ clen * dim			/* packed sequence (body) */
		+ sizeof(__m128i)		/* margin at the tail */
	);

	/* save length and seq pointers */
	*q = (dz_query_t){
		.blen  = qlen == 0 ? 0 : (clen / DZ_L),
		.seq   = query,
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

static __dz_vectorize
size_t dz_pack_query_forward_core(dz_query_t *q, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	/*
	 * determine lengths
	 *                         clen               blen
	 * qlen == k * DZ_L        (k + 1) * DZ_L      k      * DZ_L
	 * qlen == k * DZ_L - 1     k      * DZ_L     (k - 1) * DZ_L
	 * qlen == k * DZ_L - l     k      * DZ_L     (k - 1) * DZ_L
	 */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));
	size_t const dim = pack->calc_dim(clen);
	debug("clen(%zu), blen(%zu), dim(%zu)", clen, blen, dim);

	/* calculate profile for each possible reference-side base (residue) */
	for(size_t j = 0; j < dim; j++) {
		debug("j(%zu)", j);

		/* extract array pointers for this row */
		int8_t const *score_matrix = &profile->matrix[j * DZ_QUERY_MAT_SIZE];
		uint8_t *a = &q->arr[j * clen];

		/* until the end of the query sequence */
		__m128i pv = _mm_setzero_si128();	/* as -DZ_SCORE_OFS */
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
			// print_vector_8(qv);
			__m128i const tv = pack->conv(score_matrix, 0, qv);

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15));

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */ {
			char qbuf[16] = { 0 };
			for(size_t i = 0; i < qlen - blen; i++) {	/* qlen - blen < 16 */
				qbuf[i] = query[i + blen];
			}

			__m128i const qv = _mm_loadu_si128((__m128i const *)&qbuf[0]);
			__m128i const tv = pack->conv(score_matrix, 0, qv);
			_mm_storeu_si128((__m128i *)&a[blen], _mm_alignr_epi8(tv, pv, 15));
		}

		/* clear tail remainders */
		for(size_t i = qlen + 1; i < clen; i++) {
			a[i] = pack->invalid;		/* as -DZ_SCORE_OFS */
		}
	}
	// debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}

static __dz_vectorize
size_t dz_pack_query_reverse_core(dz_query_t *q, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	/* constants */
	size_t const clen = dz_roundup(qlen + 1, DZ_L);
	size_t const blen = dz_rounddown(qlen, sizeof(__m128i));
	size_t const dim = pack->calc_dim(clen);
	debug("clen(%zu), blen(%zu), dim(%zu)", clen, blen, dim);

	static uint8_t const rev[16] __attribute__(( aligned(16) )) = {
		15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
	};
	__m128i const rv = _mm_load_si128((__m128i const *)rev);


	/* calculate profile for each possible reference-side base (residue) */
	for(size_t j = 0; j < dim; j++) {
		debug("j(%zu)", j);

		/* extract array pointers for this row */
		int8_t const *score_matrix = &profile->matrix[j * DZ_QUERY_MAT_SIZE];
		uint8_t *a = &q->arr[j * clen];

		/* toward the head of the query sequence */
		__m128i pv = _mm_setzero_si128();
		for(size_t i = 0; i < blen; i += sizeof(__m128i)) {
			__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
			__m128i const tv = pack->conv(score_matrix, 1, _mm_shuffle_epi8(qv, rv));		/* bswapdq before apply conversion */

			/* shift by one to make room for the top row */
			_mm_storeu_si128((__m128i *)&a[i], _mm_alignr_epi8(tv, pv, 15));

			/* save for the next iteration */
			pv = tv;
		}

		/* continue the same conversion on the remainings */ {
			char qbuf[16] = { 0 };
			for(size_t i = 0; i < qlen - blen; i++) {
				qbuf[i] = query[qlen - blen - 1 - i];
				debug("i(%zu), c(%c, %u), c(%c, %u)", i, query[qlen - blen - 1 - i], query[qlen - blen - 1 - i], qbuf[i], qbuf[i]);
			}

			__m128i const qv = _mm_loadu_si128((__m128i const *)&qbuf[0]);
			__m128i const tv = pack->conv(score_matrix, 1, qv);		/* already reversed */
			_mm_storeu_si128((__m128i *)&a[blen], _mm_alignr_epi8(tv, pv, 15));
		}

		/* clear tail remainders */
		for(size_t i = qlen + 1; i < clen; i++) {
			a[i] = pack->invalid;		/* as -DZ_SCORE_OFS */
		}
	}
	// debug("qlen(%lu), q(%s)", qlen, query);
	return(clen);
}

/* direction-folded API */
static __dz_vectorize
dz_query_t *dz_pack_query_core(dz_arena_t *mem, dz_profile_t const *profile, uint8_t const *query, size_t qlen, dz_pack_query_t const *pack)
{
	dz_query_t *q = dz_pack_query_alloc_mem(mem, profile, query, qlen, pack);
	if(pack->dir) {
		dz_pack_query_reverse_core(q, profile, query, qlen, pack);
	} else {
		dz_pack_query_forward_core(q, profile, query, qlen, pack);
	}
	return(q);
}




/*
 * DP matrix internal objects
 */

/* Smith-Waterman-Gotoh layer indices */
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




/* coordinate and match information for traceback */
typedef struct {
	uint16_t ch;		/* fetched char (after conversion) */
	uint16_t adj;		/* difference of offset to the previous column */
	uint32_t idx;		/* forwarded length */
} dz_col_tracker_t;

/* placed just after every score vector to indicate the length */
typedef struct dz_range_s {
	uint32_t sblk, eblk;
} dz_range_t;
#define dz_range(_p)				( (dz_range_t *)(_p) )
#define dz_crange(_p)				( (dz_range_t const *)(_p) )

typedef struct {
	uint32_t column;	/* column count */
	uint32_t section;	/* section count */
} dz_ref_cnt_t;

typedef struct {
	int32_t score;		/* without offset */
	uint32_t inc;		/* with offset */
	struct dz_cap_s const *cap;
} dz_max_t;

/* as working buffer */
/* typedef */
struct dz_state_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;
	dz_max_t max;
} /* dz_state_t */;
#define dz_pstate(_p)				( (dz_state_t **)(_p) )
#define dz_cpstate(_p)				( (dz_state_t const **)(_p) )


/* struct dz_forefront_s is alias of dz_state_t */
/* typedef */
struct dz_forefront_s {
	dz_range_t range;
	dz_ref_cnt_t cnt;
	int32_t max, inc;
	struct dz_cap_s const *cap;
} /* dz_forefront_t */;
dz_static_assert(sizeof(dz_forefront_t) == sizeof(dz_state_t));
dz_static_assert(offsetof(dz_forefront_t, cnt)   == offsetof(dz_state_t, cnt));
dz_static_assert(offsetof(dz_forefront_t, range) == offsetof(dz_state_t, range));
dz_static_assert(offsetof(dz_forefront_t, max)   == offsetof(dz_state_t, max.score));

/* pointer casting */
#define dz_cff(_p)					( (dz_forefront_t const *)(_p) )
#define dz_cpff(_p)					( (dz_forefront_t const **)(_p) )

/* control */
#define dz_profile_root(_profile)	( dz_cpff(&(_profile)->root) )
#define dz_root(_self)				( dz_profile_root((_self)->profile) )

/* X-drop */
#define dz_is_terminated(_ff)		( dz_cff(_ff)->range.sblk >= dz_cff(_ff)->range.eblk )



/* merge incoming vectors; using SIMD max to avoid test-and-branch */
#if 1
typedef struct {
	__m128i range_cnt;	/* dz_range_t and dz_ref_cnt_t pair on xmm */
	__m128i max;		/* dz_max_t on xmm */
} dz_state_vec_t;

static __dz_vectorize
dz_state_vec_t dz_merge_state_intl(dz_state_t const **ff, size_t fcnt)
{
	/* for negating sblk to take min */
	static int32_t const neg[4] __attribute__(( aligned(16) )) = { -1, 0, 0, 0 };
	__m128i const negv = _mm_load_si128((__m128i const *)neg);

	/* fold max; two working registers can be concatenated into __m256i */
	__m128i range_cnt = _mm_setzero_si128();
	__m128i max       = _mm_setzero_si128();
	for(size_t i = 0; i < fcnt; i++) {
		/* load: sblk, eblk, ccnt, scnt */
		__m128i const tv  = _mm_loadu_si128((__m128i const *)&ff[i]->range.sblk);
		__m128i const rcv = _mm_xor_si128(tv, negv);	/* negate sblk */

		/* load: max.score */
		__m128i const mv  = _mm_loadu_si128((__m128i const *)&ff[i]->max.score);

		/* update max vectors; signed max for score */
		range_cnt = _mm_max_epu32(range_cnt, rcv);
		max       = _mm_max_epi32(max,       mv);

		debug("i(%zu), fcnt(%zu), ff(%p), range(%u, %u, %u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
			i, fcnt, ff[i], ff[i]->range.sblk, ff[i]->range.eblk, (uint32_t)_mm_extract_epi32(range_cnt, 0), (uint32_t)_mm_extract_epi32(range_cnt, 1),
			ff[i]->cnt.column, ff[i]->cnt.section, ff[i]->max.score, ff[i]->max.inc
		);
	}

	return((dz_state_vec_t){
		.range_cnt = _mm_xor_si128(range_cnt, negv),	/* flip sblk again before save */
		.max       = _mm_and_si128(max,       negv)		/* clear out inc and cap */
	});
}

static __dz_vectorize
dz_state_t dz_merge_state(dz_state_t const **ff, size_t fcnt)
{
	/* (cap = NULL, inc = dz_add_ofs(0), score) */
	static int32_t const init[4] __attribute__(( aligned(16) )) = { 0, dz_add_ofs(0), 0, 0 };
	__m128i const iv = _mm_load_si128((__m128i const *)init);

	/* merge */
	dz_state_vec_t merged = dz_merge_state_intl(ff, fcnt);

	/* save them all */
	dz_state_t state __attribute__(( aligned(16) ));
	_mm_store_si128((__m128i *)&state.range.sblk, merged.range_cnt);				/* (eblk, sblk, section, column) */
	_mm_store_si128((__m128i *)&state.max.score,  _mm_add_epi32(iv, merged.max));	/* (cap = NULL, inc = dz_add_ofs(0), score) */

	debug("merged state, range(%u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
		state.range.sblk, state.range.eblk, state.cnt.column, state.cnt.section, state.max.score, state.max.inc
	);
	return(state);
}
#else
static __dz_vectorize
dz_state_t dz_merge_state(dz_state_t const **ff, size_t fcnt)
{
	/* load constant */
	static int32_t const neg[4] __attribute__(( aligned(16) )) = { -1, 0, 0, 0 };	/* negate spos to take min */
	__m128i const nv = _mm_load_si128((__m128i const *)neg);

	/* working registers */
	__m128i mv0 = nv;
	__m128i mv1 = _mm_setzero_si128();

	for(size_t i = 0; i < fcnt; i++) {
		/* sblk, eblk, ccnt, scnt */
		__m128i const v0 = _mm_loadu_si128((__m128i const *)&ff[i]->range.sblk);
		mv0 = _mm_max_epu32(mv0, _mm_xor_si128(v0, nv));

		/* max */
		__m128i const v1 = _mm_loadu_si128((__m128i const *)&ff[i]->max.score);
		mv1 = _mm_max_epu32(mv1, v1);

		debug("i(%zu), fcnt(%zu), ff(%p), range(%u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
			i, fcnt, ff[i], ff[i]->range.sblk, ff[i]->range.eblk, ff[i]->cnt.column, ff[i]->cnt.section, ff[i]->max.score, ff[i]->max.inc
		);
	}

	/* save */
	dz_state_t state __attribute__(( aligned(16) ));
	_mm_store_si128((__m128i *)&state.range.sblk, _mm_xor_si128(mv0, nv));	/* (eblk, sblk, section, column) */
	_mm_store_si128((__m128i *)&state.max.score,  _mm_and_si128(mv1, nv));	/* (cap = NULL, inc = 0, score) */
	state.max.inc = dz_add_ofs(0);

	debug("merged state, range(%u, %u), ccnt(%u), scnt(%u), max(%d), inc(%d)",
		state.range.sblk, state.range.eblk, state.cnt.column, state.cnt.section, state.max.score, state.max.inc
	);
	return(state);
}
#endif

static __dz_vectorize
void dz_fixup_state(dz_state_t *state, dz_query_t const *query)
{
	state->range.eblk = dz_min2(state->range.eblk, query->blen);
	return;
}

static __dz_vectorize
void dz_finalize_state(dz_state_t *state, dz_col_tracker_t const *tracker)
{
	size_t const cols  = tracker->idx;
	state->cnt.column += cols;
	state->cnt.section++;

	uint32_t const inc = state->max.inc;		/* with offset; FIXME: fold tracker->adj ?? */
	state->max.score += dz_rm_ofs(inc);
	debug("idx(%zu), score(%d), inc(%u)", cols, state->max.score, inc);
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
	uint16_t xt, _pad[6], init;

	/* incoming matrices */
	struct {
		dz_state_t const **ptr;
		size_t cnt;
	} incoming;

	/* states */
	dz_state_t state;
	dz_col_tracker_t tracker;
} dz_work_t;

static __dz_vectorize
void dz_load_state(dz_work_t *w, dz_state_t const **ff, size_t fcnt)
{
	/* save pointer */
	w->incoming.ptr = ff;
	w->incoming.cnt = fcnt;

	/* init column state */
	w->tracker.ch  = DZ_HEAD_RCH;
	w->tracker.adj = 0;
	w->tracker.idx = 0;

	/* merge state (before merging vectors) */
	w->state = dz_merge_state(ff, fcnt);
	return;
}

static __dz_vectorize
void dz_load_score_vector(dz_work_t *w, dz_profile_t const *profile)
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
	w->xt   = profile->xt;
	w->init = profile->init;
	return;
}




/* DP matrix scaffolds */

/* placed at the head */
typedef struct dz_head_s {
	int32_t base;		/* base score; always positive */
	uint32_t zero;		/* always zero */

	uint16_t magic;		/* DZ_HEAD_RCH for merged head; internal link otherwise */
	uint16_t zero2;
	uint32_t fcnt;		/* #merged vectors (n_forefronts) */
} dz_head_t;
#define dz_head(_p)					( (dz_head_t *)(_p) )
#define dz_chead(_p)				( (dz_head_t const *)(_p) )

/* followed by dz_tail_t; sblk and eblk are shared to forefront_s */
typedef struct dz_cap_s {
	dz_range_t range;	/* column range */
	dz_col_tracker_t tracker;		/* ch, adj, idx */
} dz_cap_t;
dz_static_assert(sizeof(dz_cap_t) % sizeof(__m128i) == 0);
#define dz_cap(_p)					( (dz_cap_t *)(_p) )
#define dz_ccap(_p)					( (dz_cap_t const *)(_p) )

/* dz_head_t and dz_cap_t are compatibile for head two elements */
dz_static_assert(offsetof(dz_head_t, magic) == offsetof(dz_cap_t, tracker.ch));
dz_static_assert(offsetof(dz_head_t, zero)  == offsetof(dz_cap_t, range.eblk));

/* DP matrix object (exported) */
typedef struct dz_tail_s {
	dz_state_t state;
	dz_query_t const *query;
} dz_tail_t;
// dz_static_assert(sizeof(dz_tail_t) % sizeof(__m128i) == 0);
#define dz_tail(_p)			( (dz_tail_t *)(_p) )
#define dz_ctail(_p)		( (dz_tail_t const *)(_p) )





/* metadata handling */
#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)

static dz_tail_t const *dz_restore_tail(dz_state_t const *ff);

/* reference sequence metadata (will not be exposed outside) */
typedef struct dz_meta_s {
	uint32_t rid, dir;
} dz_meta_t;
#define dz_meta(_p)					( (dz_meta_t *)(_p) )
#define dz_cmeta(_p)				( (dz_meta_t const *)(_p) )
#define DZ_TAIL_MARGIN				( sizeof(dz_meta_t) )

static __dz_vectorize
void dz_save_meta(dz_state_t const *ff, dz_meta_t meta)
{
	dz_tail_t const *tail = dz_restore_tail(ff);
	dz_meta_t *p = dz_meta(tail + 1);
	*p = meta;
	return;
}

static __dz_vectorize
dz_meta_t dz_extract_meta(dz_tail_t const *tail)
{
	dz_meta_t const *p = dz_cmeta(tail + 1);
	dz_meta_t meta = *p;
	return(meta);
}

static __dz_vectorize
uint32_t dz_extract_id(dz_tail_t const *tail)
{
	dz_meta_t const *p = dz_cmeta(tail + 1);
	return(p->rid);
}

#else

#define DZ_TAIL_MARGIN				( 0 )
static __dz_vectorize
uint32_t dz_extract_id(dz_tail_t const *tail)
{
	dz_unused(tail);
	return(0);
}

#endif	/* DZ_WRAPPED_API */




/* scaffold handling */

static __dz_force_inline
size_t dz_calc_column_size(size_t sblk, size_t eblk)
{
	size_t const span = eblk - sblk;
	size_t const column_size = sizeof(dz_swgv_t) * (span + dz_max2(span, 8));
	size_t const size = (
		  sizeof(dz_cap_t)	/* header size */
		+ column_size		/* estimated column size */
		+ sizeof(dz_state_t)	/* tail cap */
		+ DZ_TAIL_MARGIN	/* metadata */
	);
	return(size);
}

static __dz_vectorize
dz_swgv_t *dz_init_root_head(dz_head_t *head)
{
	head->base  = 0;
	head->zero  = 0;		/* is_head */

	head->magic = DZ_ROOT_RCH | DZ_HEAD_RCH;
	head->zero2 = 0;
	head->fcnt  = 0;		/* is_root */
	return(dz_swgv(head + 1));
}

static __dz_vectorize
uint64_t dz_is_root(dz_cap_t const *cap)
{
	/* test for root tail */
	// return(cap->range.eblk == 0);
	return((dz_chead(cap)->magic & DZ_ROOT_RCH) != 0);
}

static __dz_vectorize
dz_cap_t *dz_slice_head(dz_work_t *w, dz_arena_t *mem)
{
	/* forefront pointer array size */
	size_t const fsize = dz_roundup(
		sizeof(dz_state_t *) * w->incoming.cnt,
		sizeof(__m128i)
	);

	/* column size (including head and internal caps) */
	size_t const csize = dz_calc_column_size(w->state.range.sblk, w->state.range.eblk);

	/* allocate mem */
	dz_head_t *head = dz_head(dz_add_ptr(dz_reserve_stack(mem, fsize + csize), fsize));
	debug("csize(%zu), fsize(%zu), stack(%p, %p, %zu), head(%p)", csize, fsize, mem->stack.top, mem->stack.end, dz_arena_stack_rem(mem), head);

	/* fill magic numbers */
	head->magic = DZ_HEAD_RCH;
	head->zero  = 0;		/* overlaps range.eblk */

	/* push head-cap info */
	head->fcnt  = w->incoming.cnt;
	head->base  = w->state.max.score;	/* save merging adjustment */
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
dz_swgv_t *dz_slice_column(dz_cap_t *prev_cap, size_t sblk)
{
	dz_swgv_t *col = dz_swgv(prev_cap + 1);
	return(col - sblk);
}

static __dz_vectorize
dz_cap_t *dz_cap_column(dz_swgv_t *col, size_t eblk)
{
	/* put cap at the tail of current column */
	return(dz_cap(&col[eblk]));
}

static __dz_vectorize
dz_swgv_t const *dz_restore_column(dz_cap_t const *cap)
{
	return(dz_cswgv(cap) - cap->range.eblk);
}

static __dz_vectorize
dz_cap_t const *dz_unwind_cap(dz_cap_t const *cap)
{
	dz_swgv_t const *col = dz_restore_column(cap);
	debug("cap(%p), range(%u, %u), col(%p), prev_cap(%p)", cap, cap->range.sblk, cap->range.eblk, col, dz_ccap(col + cap->range.sblk) - 1);
	return(dz_ccap(col + cap->range.sblk) - 1);
}


static __dz_vectorize
dz_state_t const *dz_slice_tail(dz_work_t *w, dz_swgv_t *col, dz_arena_t *mem, dz_query_t const *query)
{
	dz_tail_t *tail = dz_tail(dz_cap_column(col, w->state.range.eblk));
	debug("col(%p), range(%u, %u), tail(%p, %p)", col, w->state.range.sblk, w->state.range.eblk, tail, &tail->state);

	/* just copy */
	tail->state = w->state;
	tail->query = query;

	/* write back stack pointer */
	size_t const size = dz_roundup(
		sizeof(dz_tail_t) + DZ_TAIL_MARGIN,		/* make room for metadata */
		sizeof(__m128i)
	);
	dz_save_stack(mem, dz_add_ptr(tail, size));
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
	debug("tail(%p, %p), range(%u, %u), col(%p)", tail, ff, ff->range.sblk, ff->range.eblk, dz_cswgv(tail) - tail->state.range.eblk);
	return(dz_cswgv(tail) - tail->state.range.eblk);
}


/* internal link (ilink) */
static __dz_vectorize
dz_swgv_t *dz_slice_ilink(dz_cap_t const *prev_cap, uint8_t *ptr)
{
	/* save link */
	dz_state_t const **ff = ((dz_state_t const **)ptr) + 2;
	ff[-1] = dz_state(&prev_cap->range);

	/* save head info */
	dz_head_t *head = dz_head(ff);
	head->base  = 0;
	head->zero  = 0;

	/* overlaps with dz_col_tracker_t */
	head->magic = DZ_HEAD_RCH + 1;	/* to make the value different from head */
	head->zero2 = 0;
	head->fcnt  = 1;		/* #incoming vectors == 1 */
	return(dz_swgv(head + 1));
}

static __dz_vectorize
uint64_t dz_is_ilink(dz_cap_t const *cap)
{
	return(cap->tracker.ch != DZ_HEAD_RCH);
}

static __dz_vectorize
dz_cap_t const *dz_rewind_ilink(dz_cap_t const *cap)
{
	dz_state_t const **ff = dz_cpstate(cap);
	return(dz_ccap(dz_restore_tail(ff[-1])));
}

static __dz_vectorize
dz_swgv_t *dz_slice_cap_core(dz_work_t *w, dz_swgv_t *prev_col, dz_arena_t *mem)
{
	dz_cap_t *cap = dz_cap_column(prev_col, w->state.range.eblk);
	debug("prev_col(%p), cap(%p), range(%u, %u), tracker(%u, %u)", prev_col, cap, w->state.range.sblk, w->state.range.eblk, w->tracker.ch, w->tracker.idx);

	/* save column tracker and range for traceback */
	cap->range   = w->state.range;
	cap->tracker = w->tracker;			/* char saved; adj cleared; column count (idx) saved */
	dz_save_stack(mem, cap + 1);		/* save stack pointer */

	/* add stack if memory starved */
	size_t size = 2 * sizeof(dz_state_t *);
	size += dz_calc_column_size(w->state.range.sblk, w->state.range.eblk);

	/* allocate memory from stack */
	uint8_t *ptr = dz_reserve_stack(mem, size);
	debug("cap(%p, %p)", cap, dz_cap(ptr) - 1);
	if(dz_likely(dz_cap(ptr) == cap + 1)) { return(dz_swgv(cap + 1)); }

	/* new stack allocated */
	return(dz_slice_ilink(cap, ptr));
}

static __dz_vectorize
dz_swgv_t *dz_slice_cap(dz_work_t *w, dz_swgv_t *prev_col, dz_arena_t *mem)
{
	dz_swgv_t *col = dz_slice_cap_core(w, prev_col, mem);
	return(col - w->state.range.sblk);
}




/* scan */
static __dz_vectorize
uint64_t dz_is_scan(dz_profile_t const *profile)
{
	return(profile->init == (uint16_t)dz_add_ofs(0));
}

/* merge incoming vectors */
static __dz_vectorize
void dz_merge_init_col(dz_work_t *w, dz_swgv_t *col)
{
	__m128i const minv = _mm_set1_epi16(dz_add_ofs(INT16_MIN));
	dz_swgv_t const v = { minv, minv, minv };

	for(uint64_t p = w->state.range.sblk; p < w->state.range.eblk; p++) {
		dz_store_swgv(&col[p], v);
	}
	return;
}

static __dz_vectorize
uint16_t dz_merge_calc_adj(dz_work_t *w, dz_state_t const *ff)
{
	/* base score at the head of this segment */
	int32_t const base = ff->max.score - dz_rm_ofs(ff->max.inc);

	/* calc decrease amount from the current max */
	int32_t const dec  = w->state.max.score - base;		/* we are sure the value is always >= 0 */
	return(dec < UINT16_MAX ? dec : UINT16_MAX);		/* clip if overflowed */
}

static __dz_vectorize
void dz_merge_fold_col(dz_work_t *w, dz_swgv_t *col, dz_profile_t const *profile)
{
	/* squash max if scan */
	if(dz_is_scan(profile)) { w->state.max.score = 0; }

	/* for each incoming vectors */
	for(size_t i = 0; i < w->incoming.cnt; i++) {
		dz_state_t const *ff = w->incoming.ptr[i];

		uint16_t const adj = dz_is_scan(profile) ? 0 : dz_merge_calc_adj(w, ff);
		__m128i const adjv = _mm_set1_epi16(adj);
		debug("adj(%d)", adj);

		dz_swgv_t const *prev_col = dz_restore_tail_column(ff);
		for(uint64_t p = ff->range.sblk; p < ff->range.eblk; p++) {
			dz_swgv_t const v = dz_load_swgv(&prev_col[p]);
			dz_swgv_t const x = dz_load_swgv(&col[p]);

			/* adjust offset; unsigned saturated subtraction; clipped when underflowed */
			dz_swgv_t const u = dz_max_swgv(x, dz_subs_swgv(v, adjv));
			print_vector(u.s);
			dz_store_swgv(&col[p], u);
		}
	}
	return;
}

static __dz_vectorize
dz_swgv_t *dz_merge_column(dz_work_t *w, dz_cap_t *prev_cap, dz_profile_t const *profile)
{
	/* prev_cap is head cap */
	dz_swgv_t *col = dz_slice_column(prev_cap, w->state.range.sblk);

	dz_merge_init_col(w, col);
	dz_merge_fold_col(w, col, profile);
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
	/* S at i == 0 */
	__m128i isv = _mm_load_si128((__m128i const *)&(&w->init)[-7]);

	/* FIXME: offsetted */
	fw->f    = w->minv;
	fw->e    = w->minv;
	fw->s    = w->minv;
	fw->ps   = w->state.range.sblk == 0 ? isv : w->minv;	/* ofs(0) if scan, min otherwise */
	fw->maxv = w->minv;
	fw->xtv  = _mm_set1_epi16(w->state.max.inc - w->xt);	/* offset already included; next offset == current max thus X-drop threshold is always -xt */
	return;
}



static __dz_vectorize
void dz_fill_update_vector(dz_work_t *w, dz_fill_work_t *fw, __m128i sc, __m128i bv)
{
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

	/* done */
	fw->f = tf;
	fw->e = te;
	fw->s = us;

	print_vector(us); /* print_vector(bv); */
	return;
}

static __dz_vectorize
void dz_fill_update_maxv(dz_fill_work_t *fw)
{
	/* update max */
	#ifdef DZ_FULL_LENGTH_BONUS
		fw->maxv = _mm_max_epu16(fw->maxv, _mm_add_epi16(us, bv));
	#else
		dz_unused(bv);
		fw->maxv = _mm_max_epu16(fw->maxv, us);
	#endif
	print_vector(fw->maxv);
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
uint16_t dz_fill_fold_max(__m128i v)
{
	__m128i t = _mm_max_epu16(v, _mm_srli_si128(v, 8));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 4));
	t = _mm_max_epu16(t, _mm_srli_si128(t, 2));
	return(_mm_extract_epi16(t, 0));		/* with offset */
}

static __dz_vectorize
uint64_t dz_fill_test_xdrop(dz_fill_work_t *fw)
{
	__m128i xtest = _mm_subs_epu16(fw->s, fw->xtv);

	/* print_vector(_s); print_vector(_xtv); */
	return(dz_is_all_zero(xtest));
}

static __dz_vectorize
void dz_fill_remove_top(dz_work_t *w, dz_fill_work_t *fw)
{
	/* this function is called when the topmost vector is removed by X-drop criterion */

	/* first adjust pointer and index */
	w->state.range.sblk++;
	fw->col--;

	/* break dependency link from the removed cell */
	fw->f = w->minv;
	return;
}


static __dz_vectorize
uint64_t dz_fill_is_bottom(dz_work_t *w, dz_query_t const *query)
{
	/* eblk promoted to size_t */
	return(w->state.range.eblk >= query->blen);
}

static __dz_vectorize
uint64_t dz_fill_column_body(dz_work_t *w, dz_fill_work_t *fw, dz_fetcher_t *fetcher, dz_profile_t const *profile, dz_query_t const *query)
{
	/* return 0 if continue to tail-fill loop */
	size_t p = --w->state.range.sblk;	/* put margin at the head for simplicity of the outer loop */
	do {
		w->state.range.sblk++;

		/* clip head */
		fw->col--;
		fw->f = w->minv;

		/* load-update-save loop; fallthrough to return(is_bottom) */
		while(dz_unlikely(++p >= w->state.range.eblk)) {
			/* load previous DP matrix vectors */
			dz_fill_load_vector(fw, p);

			/* load score profile and update vectors */
			__m128i const sc = fetcher->get_profile(fetcher->opaque, profile->matrix, query, p * DZ_L);
			__m128i const bv = dz_query_get_bonus(query, p);
			dz_fill_update_vector(w, fw, sc, bv);

			/* save if at least one cell survived */
			if(dz_unlikely(dz_fill_test_xdrop(fw))) {
				goto _xdrop;			/* skip the return right after this loop */
			}

			/* save and update max vector */
			dz_fill_update_maxv(fw);
			dz_fill_store_vector(fw, p);

		}
		return(dz_fill_is_bottom(w, query));	/* extend further if query remains */

		/* break */
	_xdrop:;
	} while(p == w->state.range.sblk);	/* continue the loop when the xdrop is invoked at the top */

	/* write back tail when xdrop */
	w->state.range.eblk = p;

	debug("range(%u, %u), blen(%zu)", w->state.range.sblk, w->state.range.eblk, query->blen);
	return(1);							/* do not continue fill-in when xdrop */
}

static __dz_vectorize
void dz_fill_column_tail(dz_work_t *w, dz_fill_work_t *fw, dz_fetcher_t *fetcher, dz_profile_t const *profile, dz_query_t const *query)
{
	/* forefront extension; clip the column length if too long */
	fw->e = w->minv;
	fw->s = w->minv;

	/* update the last vector */
	__m128i const sc = fetcher->get_profile(fetcher->opaque, profile->matrix, query, w->state.range.eblk * DZ_L);
	__m128i const bv = dz_query_get_bonus(query, w->state.range.eblk);
	dz_fill_update_vector(w, fw, sc, bv);

	/* until all the cells drop */
	do {
		if(dz_fill_test_xdrop(fw)) { break; }
		dz_fill_store_vector(fw, w->state.range.eblk);

		fw->f = _mm_subs_epu16(fw->f, w->iev8);
		fw->s = _mm_subs_epu16(fw->s, w->iev8);
	} while(++w->state.range.eblk < query->blen);
	return;
}

static __dz_vectorize
void dz_fill_update_max(dz_work_t *w, dz_fill_work_t *fw)
{
	/* update max; the actual score (absolute score accumulated from the origin) is expressed as max + inc; 2 x cmov */
	uint32_t const inc = dz_fill_fold_max(fw->maxv);/* offset NOT removed */

	/* with offset */
	if(dz_cmp_max(inc, w->state.max.inc)) {
		w->state.max.inc = inc;
		w->state.max.cap = dz_cap_column(fw->col, w->state.range.eblk);
	}
	debug("inc(%d), base(%d), score(%d), cap(%p, %p)", dz_rm_ofs(inc), w->state.max.score, w->state.max.score + dz_rm_ofs(inc), dz_cap_column(fw->col, w->state.range.eblk), w->state.max.cap);
	return;
}

static __dz_vectorize
dz_swgv_t *dz_fill_column(dz_work_t *w, dz_swgv_t *prev_col, dz_arena_t *mem, dz_fetcher_t *fetcher, dz_profile_t const *profile, dz_query_t const *query)
{
	dz_fill_work_t fw;

	/* slice new cap */
	fw.prev_col = prev_col;
	fw.col = dz_slice_cap(w, prev_col, mem);

	/* init vectors */
	dz_fill_work_init(w, &fw);

	/* if reached the forefront of the query sequence, finish the extension */
	if(!dz_fill_column_body(w, &fw, fetcher, profile, query)) {
		dz_fill_column_tail(w, &fw, fetcher, profile, query);
	}

	/* create cap object that contains [sblk, eblk) range (for use in the traceback routine) */
	dz_fill_update_max(w, &fw);
	return(fw.col);
}

static __dz_vectorize
uint64_t dz_is_end(dz_work_t *w)
{
	return((w->state.range.eblk - w->state.range.sblk) == 0);
}

static __dz_vectorize
uint64_t dz_test_overflow(dz_work_t *w)
{
	return(w->state.max.inc > (UINT16_MAX - DZ_OVERFLOW_MARGIN));	/* with offset */
}

static __dz_vectorize
void dz_reset_base(dz_work_t *w, dz_swgv_t *col)
{
	uint16_t const adj = dz_rm_ofs(w->state.max.inc);	/* removing offset */
	__m128i const adjv = _mm_set1_epi16(adj);

	/* restore cap pointer; cap before col holds range for previous col */
	dz_cap_t *cap = dz_cap(col + w->state.range.sblk) - 1;
	cap->tracker.adj    = adj;		/* overwrite adj */
	w->state.max.score += adj;
	w->state.max.inc    = dz_add_ofs(0);
	debug("reset base, adj(%u), cap(%p)", adj, cap);

	/* subtract adjustment */
	for(uint64_t p = w->state.range.sblk; p < w->state.range.eblk; p++) {
		/* load-compensate-store */
		dz_swgv_t const v  = dz_load_swgv(&col[p]);
		dz_swgv_t const cv = dz_subs_swgv(v, adjv);
		dz_store_swgv(&col[p], cv);
	}
	return;
}


static __dz_vectorize
uint64_t dz_fetch_next(dz_work_t *w, dz_fetcher_t *fetcher, dz_profile_t const *profile, dz_query_t const *query)
{
	dz_fill_fetch_t next = fetcher->fetch_next(fetcher->opaque, profile->matrix, query);
	if(next.is_term) { return(1); }

	/* FIXME: update adj?? */
	w->tracker.ch = next.ch;
	w->tracker.idx++;
	return(0);
}

static __dz_vectorize
dz_state_t const *dz_extend_core(
	dz_arena_t *mem,
	dz_profile_t const *profile,
	dz_query_t const *query,
	dz_fetcher_t *fetcher,
	dz_state_t const **ff,
	size_t fcnt)
{
	if(fcnt == 0) { return(NULL); }		/* invalid */

	/* init working buffer */
	dz_work_t w __attribute__(( aligned(16) ));
	dz_load_state(&w, ff, fcnt);		/* iterate over the incoming edge objects to get the current max and range */
	dz_fixup_state(&w.state, query);
	// debug("(%u, %u)", w.tracker.ch, w.tracker.idx);

	/* load constants */
	dz_load_score_vector(&w, profile);

	/* save head info */
	dz_cap_t *cap = dz_slice_head(&w, mem);	/* range obtained in dz_load_state */
	dz_save_incoming(cap, ff, fcnt);

	/* merge incoming vectors */
	dz_swgv_t *col = dz_merge_column(&w, cap, profile);
	debug("merged, cap(%p), col(%p)", cap, col);

	/* until X-drop */
	while(!dz_is_end(&w)) {
		if(dz_fetch_next(&w, fetcher, profile, query)) { break; }

		/* forward one */
		col = dz_fill_column(&w, col, mem, fetcher, profile, query);

		/* place internal merging vector if score exceeded floor */
		if(dz_unlikely(dz_test_overflow(&w))) {
			dz_reset_base(&w, col);
		}
	}

	/* update max and range */
	dz_finalize_state(&w.state, &w.tracker);

	/* save them in tail object */
	return(dz_slice_tail(&w, col, mem, query));
}




/*
 * profile constructor
 */
static __dz_vectorize
dz_profile_t *dz_alloc_profile(dz_allocator_t *alloc)
{
	dz_profile_t *profile = (dz_profile_t *)alloc->fp(alloc->ctx,
		  sizeof(dz_profile_t)
		+ DZ_REF_MAT_SIZE * DZ_QUERY_MAT_SIZE
		+ 2 * sizeof(__m128i)
	);
	return(profile);
}

static __dz_vectorize
void dz_init_score_matrix(dz_profile_t *profile, int8_t const *score_matrix)
{
	/* constants */
	__m128i const tmat = _mm_loadu_si128((__m128i const *)score_matrix);
	__m128i const ofs = _mm_set1_epi8(DZ_SCORE_OFS);
	_mm_store_si128((__m128i *)&profile->matrix[0],
		_mm_add_epi8(tmat, ofs)
	);

	/* transpose */
	for(size_t i = 0; i < DZ_QUERY_MAT_SIZE; i++) {
		for(size_t j = 0; j < DZ_REF_MAT_SIZE; j++) {

			#if defined(DZ_TRANSPOSE_MATRIX) && (DZ_TRANSPOSE_MATRIX != 0)
				int8_t *p = &profile->matrix[j * DZ_QUERY_MAT_SIZE + i];
			#else
				int8_t *p = &profile->matrix[i * DZ_REF_MAT_SIZE + j];
			#endif

			/* copy */
			*p = score_matrix[i * DZ_REF_MAT_SIZE + j] + DZ_SCORE_OFS;
		}
	}

	/* padding invalid reference base */ {
		int8_t *p = &profile->matrix[DZ_REF_MAT_SIZE * DZ_QUERY_MAT_SIZE];
		_mm_store_si128((__m128i *)p, _mm_setzero_si128());
	}
	return;
}

static __dz_vectorize
void dz_init_gap_penalties(dz_profile_t *profile, dz_score_conf_t const *conf)
{
	/* insertion penalties */
	__m128i const ii = _mm_set1_epi16(conf->ins_open);
	__m128i const ie = _mm_set1_epi16(conf->ins_extend);
	_mm_store_si128((__m128i *)profile->iiv, ii);
	_mm_store_si128((__m128i *)profile->iev, ie);

	/* deletion penalties */
	__m128i const di = _mm_set1_epi16(conf->del_open);
	__m128i const de = _mm_set1_epi16(conf->del_extend);
	_mm_store_si128((__m128i *)profile->div, di);
	_mm_store_si128((__m128i *)profile->dev, de);

	/* X-drop threshold */
	uint16_t const gi = dz_max2(conf->ins_open,   conf->del_open);
	uint16_t const ge = dz_max2(conf->ins_extend, conf->del_extend);
	profile->xt = gi + ge * dz_max2(conf->max_ins_len, conf->max_del_len);	/* FIXME: X-drop threshold */
	profile->bonus = conf->full_length_bonus;
	profile->init  = dz_add_ofs(DZ_CELL_MIN);
	profile->max_ins_len = conf->max_ins_len;	/* save raw value */
	profile->max_del_len = conf->max_del_len;	/* save raw value */
	debug("gi(%u), ge(%u), xdrop_threshold(%u), full_length_bonus(%u), max_gap_len(%u, %u)", gi, ge, profile->xt, profile->bonus, profile->max_ins_len, profile->max_del_len);

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
dz_state_t const *dz_init_root_cap(dz_tail_t *tail, size_t eblk)
{
	/* initialize state */
	dz_state_t state;
	state.cnt.column  = 0;
	state.cnt.section = 0;
	state.range.sblk  = 0;
	state.range.eblk  = eblk;
	state.max.score   = 0;						/* without offset */
	state.max.inc     = dz_add_ofs(0);			/* with offset */
	state.max.cap     = NULL;
	tail->state = state;

	/* dummy */
	tail->query = NULL;
	return(&tail->state);
}

static __dz_vectorize
void dz_init_root(dz_profile_t *profile, dz_score_conf_t const *conf, dz_allocator_t *alloc)
{
	/* calc vector length; query = NULL for the first (root) column */
	size_t const max_ins_len = dz_roundup(conf->max_ins_len, DZ_L);
	size_t const blen = max_ins_len / DZ_L;

	/* fill root head */
	dz_head_t *head = dz_slice_root_head(alloc, blen);
	dz_swgv_t *col = dz_init_root_head(head);

	/* DP matrix */
	size_t eblk = dz_init_root_column(profile, conf, col, blen);

	/* tail object */
	dz_tail_t *tail = dz_tail(dz_cap_column(col, eblk));
	dz_state_t const *root = dz_init_root_cap(tail, eblk);
	debug("col(%p), range(%u, %u), tail(%p, %p)", col, root->range.sblk, root->range.eblk, tail, root);

	profile->root = dz_cff(root);
	return;
}

static __dz_vectorize
dz_profile_t *dz_init_profile(dz_allocator_t *alloc, dz_score_conf_t const *conf)
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

static __dz_vectorize
void dz_destroy_profile(dz_destructor_t *dtor, dz_profile_t *profile)
{
	dz_tail_t const *tail = dz_restore_tail(dz_cstate(profile->root));
	dz_cap_t const *head = dz_unwind_cap(dz_ccap(tail));

	dtor->fp(dtor->ctx, (void *)head);
	dtor->fp(dtor->ctx, profile);
	return;
}




/*
 * maximum-scoring position search
 */
typedef struct {
	uint32_t rpos, qpos;
} dz_max_pos_t;

static __dz_vectorize
int64_t dz_calc_max_rpos_core(dz_state_t const *ff)
{
	/* max-scoring column */
	dz_cap_t const *cap = ff->max.cap;
	if(cap == NULL) { return(0); }

	// dz_tail_t const *tail = dz_restore_tail(ff);
	dz_cap_t const *prev_cap = dz_unwind_cap(cap);
	uint32_t const rpos = (dz_is_head(prev_cap) ? 0 : prev_cap->tracker.idx);

	debug("tail(%p), prev_cap(%p, %u, %u, %u, %u), is_head(%lu), rpos(%u)",
		dz_restore_tail(ff),
		prev_cap, prev_cap->range.sblk, prev_cap->range.eblk, prev_cap->tracker.ch, prev_cap->tracker.idx, dz_is_head(prev_cap),
		rpos
	);
	return(rpos);
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
	dz_cap_t const *pcap = ff->max.cap;
	if(pcap == NULL) { return(0); }

	/* load max score */
	__m128i const maxv = _mm_set1_epi16(ff->max.inc);	/* offset already included in ff->max.inc */

	/* load query pointer */
	dz_query_t const *query = dz_restore_tail(ff)->query;

	/* loac column pointer */
	dz_swgv_t const *col = dz_restore_column(pcap);
	for(size_t p = pcap->range.sblk; p < pcap->range.eblk; p++) {
		__m128i const v = _mm_load_si128(&col[p].s);

		/* adjust full-length bonus if needed */
		#ifdef DZ_FULL_LENGTH_BONUS
			__m128i const b = dz_query_get_bonus(query, p);
			__m128i const s = _mm_add_epi16(v, b);
		#else
			dz_unused(query);
			__m128i const s = v;
		#endif

		print_vector(s);

		uint64_t eq = _mm_movemask_epi8(_mm_cmpeq_epi16(s, maxv));
		if(eq == 0) { continue; }
		
		return(dz_finalize_qpos(p, eq));
	}
	return(0);
}

static __dz_vectorize
dz_max_pos_t dz_calc_max_pos_core(dz_state_t const *ff)
{
	dz_cap_t const *pcap = ff->max.cap;
	if(pcap == NULL) {
		/* is root */
		debug("pcap == NULL");
		return((dz_max_pos_t){ 0, 0 });
	}

	debug("ff(%p), pcap(%p)", ff, pcap);
	return((dz_max_pos_t){
		.rpos = dz_calc_max_rpos_core(ff),
		.qpos = dz_calc_max_qpos_core(ff)
	});
}




/*
 * traceback
 */
/* alignment path object */
typedef struct dz_path_span_s {
	uint32_t id;
	uint32_t offset;
} dz_path_span_t;

/* alignment object */
/* typedef */
struct dz_alignment_s {
	int32_t score;
	uint32_t unused;
	uint32_t ref_length, query_length;

	/* path */
	dz_path_span_t const *span;
	uint8_t const *path;
	uint32_t span_length, path_length;

	/* gap counts */
	uint32_t mismatch_count, match_count, ins_count, del_count;
} /* dz_alignment_t */;

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

	/* match / gap counters */
	uint32_t cnt[4];

	dz_cap_t const *pcap, *ccap;
	size_t idx;
	uint16_t score, ie, de, rch;
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
uint16_t dz_trace_score_adj(dz_trace_work_t *w)
{
	/* w->score is of w->ccap; corresponding adj is saved one before */
	return(w->score + w->pcap->tracker.adj);
}
static __dz_vectorize
uint16_t dz_trace_score_raw(dz_trace_work_t *w)
{
	/* raw score for internal comparison (for vertical (F; ins) transition) */
	return(w->score);
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
int32_t dz_trace_merge_score(dz_state_t const *ff, size_t idx, size_t layer)
{
	/* extract tail */
	dz_tail_t const *tail = dz_restore_tail(ff);

	/* restore offsetted score; adj == 0 for head */
	int32_t const base = ff->max.score - dz_rm_ofs(ff->max.inc);
	int32_t const s = dz_trace_score(layer, dz_ccap(tail), idx) + base;
	return(s);
}

static __dz_vectorize
uint64_t dz_trace_rewind_merge(dz_trace_work_t *w, size_t layer)
{
	/* merging vector; load contents to find an edge */
	dz_head_t const *head = dz_chead(w->pcap);
	int32_t const prev_score = dz_trace_score(layer, w->ccap, w->idx) + head->base;
	debug("head(%p), prev_score(%u)", head, dz_rm_ofs(prev_score));

	/* load incoming vectors */
	size_t const vidx = dz_trace_vector_idx(w->idx);
	size_t const fcnt = head->fcnt;
	dz_state_t const **farr = dz_cpstate(w->pcap) - fcnt;
	for(size_t i = 0; i < fcnt; i++) {
		dz_state_t const *ff = farr[i];
		debug("i(%zu), ff(%p)", i, ff);

		/* skip if the cell is out of the vector */
		if(!dz_inside(ff->range.sblk, vidx, ff->range.eblk)) { continue; }

		/* adj[i] = w.max - (ffs[i]->max - ffs[i]->inc); base = max - inc */
		int32_t const score = dz_trace_merge_score(ff, w->idx, layer);
		debug("score(%u), prev_score(%u)", dz_rm_ofs(score), dz_rm_ofs(prev_score));
		if(prev_score == score) {
			w->pcap = dz_ccap(dz_restore_tail(ff));
			return(1);
		}
	}
	return(0);
}

static __dz_vectorize
uint64_t dz_trace_reload_section(dz_trace_work_t *w, size_t layer)
{
	/* just load pcap if internal bridge */
	if(dz_is_ilink(w->pcap)) {
		w->pcap = dz_rewind_ilink(w->pcap);
		return(0);
	}

	/* load merging vector */
	if(dz_trace_rewind_merge(w, layer)) {
		return(0);
	}

	/* something is wrong... */
	dz_trap();
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_unwind_h(dz_trace_work_t *w, size_t layer)
{
	debug("pcap(%p), cap(%p)", w->pcap, w->ccap);
	while(1) {
		w->ccap = w->pcap;
		w->pcap = dz_unwind_cap(w->ccap);
		debug("pcap(%p), cap(%p), is_head(%u), is_root(%u)", w->pcap, w->ccap, (uint32_t)dz_is_head(w->pcap), (uint32_t)dz_is_root(w->pcap));

		if(dz_likely(!dz_is_head(w->pcap))) { break; }	/* escape if not head */
		if(dz_is_root(w->pcap)) { break; }				/* escape if root */

		/* pcap is valid merging node */
		dz_trace_reload_section(w, layer);

		/* push segment info */
		uint32_t const id = dz_extract_id(dz_ctail(w->pcap));
		dz_trace_push_span(w, id);

		/* reload score */
		w->score = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx);
	}

	/* return the reference-side base */
	w->rlen++;
	w->rch = w->pcap->tracker.ch;

	debug("pcap(%p, %x), cap(%p, %x)", w->pcap, w->pcap->tracker.ch, w->ccap, w->ccap->tracker.ch);
	return(0);
}


static __dz_vectorize
uint64_t dz_trace_test_idx(dz_trace_work_t *w, dz_cap_t const *cap, size_t ofs)
{
	size_t const vidx = dz_trace_vector_idx(w->idx - ofs);
	return(!dz_inside(cap->range.sblk, vidx, cap->range.eblk));
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
	debug("score(%u, %d), op(%c), idx(%zu)", next_score, dz_rm_ofs(next_score), *w->path.ptr, w->idx);
	return;
}


static __dz_vectorize
uint64_t dz_trace_eat_match(dz_trace_work_t *w, dz_profile_t const *profile, dz_trace_get_match_t get_match) {
	if(dz_trace_test_idx(w, w->pcap, 1)) { debug("test match out of range, idx(%zu), range(%u, %u)", w->idx, w->pcap->range.sblk, w->pcap->range.eblk); return(0); }

	/* get diagonal score for this cell */
	uint16_t const s = dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx - 1);
	dz_trace_match_t m = get_match(profile->matrix, w->query, w->idx, w->rch);

	/* skip if score does not match */
	debug("test match, rch(%x), adj_score(%d, %d), match(%d, %u)", w->rch, dz_rm_ofs(dz_trace_score_adj(w)), dz_rm_ofs(s + m.score - DZ_SCORE_OFS), m.score, m.match);
	if(dz_trace_score_adj(w) != (uint16_t)(s + m.score - DZ_SCORE_OFS)) { return(0); }

	/* determine match state; compensate adj before saving */
	dz_trace_push_op(w, DZ_S_MATRIX + m.match, s);

	/* unwind both */
	dz_trace_unwind_v(w, DZ_S_MATRIX);
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_ins(dz_trace_work_t *w) {
	if(dz_trace_test_idx(w, w->ccap, 1)) { debug("test ins out of range, idx(%zu), range(%u, %u)", w->idx, w->ccap->range.sblk, w->ccap->range.eblk); return(0); }

	/* skip if score does not match */
	uint16_t const f = dz_trace_score(DZ_F_MATRIX, w->ccap, w->idx);
	if(dz_likely(dz_trace_score_raw(w) != f)) { debug("test ins score unmatch, idx(%zu), raw_score(%d), f(%u)", w->idx, dz_rm_ofs(dz_trace_score_raw(w)), f); return(0); }
	debug("ins, raw_score(%d), f(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), f, w->idx);

	while(!dz_trace_test_idx(w, w->ccap, 2)) {		/* do not move to F matrix when gap longer than 2 is not possible */
		uint16_t const x = dz_trace_score(DZ_F_MATRIX, w->ccap, w->idx - 1);
		debug("ins, raw_score(%d), x(%u), ie(%u), idx(%zu)", dz_rm_ofs(dz_trace_score_raw(w)), x, w->ie, w->idx);
		if(dz_trace_score_raw(w) != x - w->ie) { break; }

		dz_trace_push_op(w, DZ_F_MATRIX, x);
		dz_trace_unwind_v(w, DZ_F_MATRIX);
	}

	/* eat last column */
	dz_trace_push_op(w, DZ_F_MATRIX, dz_trace_score(DZ_S_MATRIX, w->ccap, w->idx - 1));
	dz_trace_unwind_v(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
uint64_t dz_trace_eat_del(dz_trace_work_t *w) {
	if(dz_trace_test_idx(w, w->pcap, 0)) { debug("test del out of range, idx(%zu), range(%u, %u)", w->idx, w->pcap->range.sblk, w->pcap->range.eblk); return(0); }

	/* skip if score does not match */
	uint16_t const e = dz_trace_score(DZ_E_MATRIX, w->ccap, w->idx);
	if(dz_likely(dz_trace_score_raw(w) != e)) { debug("test del score unmatch, idx(%zu), raw_score(%d), e(%u)", w->idx, dz_rm_ofs(dz_trace_score_raw(w)), e); return(0); }
	debug("del, raw_score(%d), f(%u)", dz_rm_ofs(dz_trace_score_raw(w)), e);

	do {
		uint16_t const x = dz_trace_score(DZ_E_MATRIX, w->pcap, w->idx);
		if(dz_trace_score_adj(w) != x - w->de) { break; }
		debug("del, adj_score(%d), x(%u), de(%u)", dz_rm_ofs(dz_trace_score_adj(w)), x, w->de);

		dz_trace_push_op(w, DZ_E_MATRIX, x);
		dz_trace_unwind_h(w, DZ_E_MATRIX);
	} while(!dz_trace_test_idx(w, w->pcap, 0));

	/* eat last row */
	dz_trace_push_op(w, DZ_E_MATRIX, dz_trace_score(DZ_S_MATRIX, w->pcap, w->idx));
	dz_trace_unwind_h(w, DZ_S_MATRIX);
	return(1);
}

static __dz_vectorize
void dz_trace_init_aln(dz_alignment_t *aln, dz_state_t const *ff, size_t idx)
{
	// aln->rrem  = dz_calc_max_rpos_core(ff);
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
		+ dz_roundup(ff->cnt.column + idx + 16, 8) * sizeof(uint8_t)	/* +1 for tail '\0' */
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
	uint32_t const id = dz_extract_id(tail);
	dz_trace_push_span(w, id);
	// *w->path.ptr = '\0';		/* make sure readable as C string */
	_mm_storeu_si128((__m128i *)w->path.ptr, _mm_setzero_si128());
	return;
}

static __dz_vectorize
void dz_trace_init_work(dz_trace_work_t *w, dz_profile_t const *profile, dz_state_t const *ff, size_t idx)
{
	/* save lengths */
	w->rlen    = 0;		/* is unknown before traceback */
	w->query   = dz_restore_tail(ff)->query;
	w->ie = profile->iev[0];
	w->de = profile->dev[0];

	/* clear counters */
	memset(w->cnt, 0, 4 * sizeof(uint32_t));

	/* load max column pointers */
	w->pcap = ff->max.cap;
	w->ccap = NULL;
	w->idx  = idx;

	/* compensate adj when save; make sure ccap == NULL */
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
dz_alignment_t const *dz_trace_finalize(dz_trace_work_t *w)
{
	/* check buffer overrun */
	debugblock({
		assert((uintptr_t)(w->aln + 1) <= (uintptr_t)w->span.ptr + 1, "(%p, %p)", w->aln + 1, w->span.ptr);
		assert((uintptr_t)w->span.base <= (uintptr_t)w->path.ptr, "(%p, %p)", w->span.base, w->path.ptr);
	});

	dz_trace_finalize_path(w->aln, w->path.ptr, w->path.base);
	dz_trace_finalize_span(w->aln, w->span.ptr, w->span.base, w->aln->path_length);

	/* save length */
	w->aln->ref_length = w->rlen - 1;

	/* save match / gap counts */
	w->aln->del_count = w->cnt[DZ_E_MATRIX];
	w->aln->ins_count = w->cnt[DZ_F_MATRIX];
	w->aln->mismatch_count = w->cnt[DZ_S_MATRIX];
	w->aln->match_count = w->cnt[DZ_S_MATRIX + 1];
	return(w->aln);
}

static __dz_vectorize
dz_alignment_t const *dz_trace_core(dz_arena_t *mem, dz_profile_t const *profile, dz_state_t const *ff, dz_trace_get_match_t get_match)
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
		if(dz_trace_eat_match(&w, profile, get_match)) { continue; }
		if(dz_trace_eat_ins(&w)) { continue; }
		if(dz_trace_eat_del(&w)) { continue; }

		dz_trap();
	}

	/* done */
	return(dz_trace_finalize(&w));
}




/*
 * wrapped APIs
 */
#if defined(DZ_WRAPPED_API) && (DZ_WRAPPED_API != 0)

static __dz_vectorize
dz_t *dz_initx(
	int8_t const *score_matrix,
	uint16_t ins_open,
	uint16_t ins_extend,
	uint16_t del_open,
	uint16_t del_extend,
	uint64_t max_ins_len,
	uint64_t max_del_len,
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
	conf.scan_ref     = 0;
	conf.full_length_bonus = full_length_bonus;
	conf.max_ins_len  = max_ins_len;
	conf.max_del_len  = max_del_len;

	self->profile = dz_init_profile(&alloc, &conf);
	return(self);
}

#ifdef DZ_FULL_LENGTH_BONUS
static __dz_vectorize
dz_t *dz_init(int8_t const *score_matrix, uint16_t gap_open, uint16_t gap_extend, uint64_t max_gap_len, uint16_t full_length_bonus)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len, max_gap_len,
		full_length_bonus
	));
}
#else		/* DZ_FULL_LENGTH_BONUS */
static __dz_vectorize
dz_t *dz_init(int8_t const *score_matrix, uint16_t gap_open, uint16_t gap_extend, uint64_t max_gap_len)
{
	return(dz_initx(score_matrix,
		gap_open, gap_extend,
		gap_open, gap_extend,
		max_gap_len, max_gap_len, 0
	));
}
#endif		/* DZ_FULL_LENGTH_BONUS */

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
	return;
}

/* packing query */
static __dz_vectorize
dz_query_t *dz_pack_query(dz_t *self, char const *query, int64_t qlen)
{
	uint32_t const dir = qlen < 0;
	char const *q  = dir ? &query[qlen] : query;
	size_t const l = dir ? -qlen : qlen;
	debug("q(%p), l(%zu), dir(%u)", q, l, dir);

	return(dz_pack_query_intl(self, q, l, dir));	/* internally calls dz_pack_query_core */
}

static __dz_vectorize
dz_query_t *dz_pack_query_forward(dz_t *self, char const *query, size_t qlen)
{
	return(dz_pack_query_intl(self, query, qlen, 0));
}

static __dz_vectorize
dz_query_t *dz_pack_query_reverse(dz_t *self, char const *query, size_t qlen)
{
	return(dz_pack_query_intl(self, query, qlen, 1));
}

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
	uint32_t const dir = rlen < 0 ? 1 : 0;		/* 1 for reverse */
	char const *p  = rlen < 0 ? &ref[rlen] : ref;
	size_t const l = rlen < 0 ? -rlen : rlen;	/* always positive */
	// debug("extend, self(%p), fcnt(%zu), query(%p), ref(%p), rlen(%d), rid(%u), p(%p), l(%zu), dir(%u)", self, fcnt, query, ref, rlen, rid, p, l, dir);

	dz_state_t const *tail = dz_extend_intl(self->mem, self->profile, query, dz_cpstate(forefronts), fcnt, p, l, dir);
	if(tail == NULL) { return(NULL); }

	/* save id and direction */
	dz_save_meta(tail, ((dz_meta_t){ .rid = rid, .dir = dir }));
	return(dz_cff(tail));
}

static __dz_vectorize
dz_forefront_t const *dz_scan(
	dz_t *self,
	dz_query_t const *query,
	dz_forefront_t const **forefronts, size_t fcnt,
	char const *ref, int32_t rlen, uint32_t rid)
{
	uint32_t const dir = rlen < 0 ? 1 : 0;		/* 1 for reverse */
	char const *p  = rlen < 0 ? &ref[rlen] : ref;
	size_t const l = rlen < 0 ? -rlen : rlen;	/* always positive */
	debug("scan, self(%p), fcnt(%zu), query(%p), ref(%p), rlen(%d), rid(%u), p(%p), l(%zu), dir(%u)", self, fcnt, query, ref, rlen, rid, p, l, dir);

	/* dup profile object */
	uint8_t buf[self->profile->size];
	dz_profile_t *pdup = (dz_profile_t *)buf;
	memcpy(pdup, self->profile, self->profile->size);
	pdup->init = dz_add_ofs(0);

	/* extend with modified profile */
	dz_state_t const *tail = dz_extend_intl(self->mem, pdup, query, dz_cpstate(forefronts), fcnt, p, l, dir);
	if(tail == NULL) { return(NULL); }

	/* save id and direction */
	dz_save_meta(tail, ((dz_meta_t){ .rid = rid, .dir = dir }));
	return(dz_cff(tail));
}

static __dz_vectorize
int64_t dz_calc_max_rpos(dz_t *self, dz_forefront_t const *ff)
{
	dz_unused(self);

	/* zero if invalid */
	if(ff == NULL) { return(0); }

	/* calculate tail rpos, convert to signed */
	int64_t rpos = dz_calc_max_rpos_core(dz_cstate(ff));

	/* fold direction in */
	dz_tail_t const *tail = dz_restore_tail(dz_cstate(ff));
	dz_meta_t meta = dz_extract_meta(tail);
	return(meta.dir ? -rpos : rpos);	/* negate if reverse */
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

	if(ff == NULL) { return(0); }
	dz_max_pos_t pos = dz_calc_max_pos_core(dz_cstate(ff));

	/* fold direction */
	dz_tail_t const *tail = dz_restore_tail(dz_cstate(ff));
	dz_meta_t meta = dz_extract_meta(tail);
	int64_t const rpos  = pos.rpos;
	int64_t const rdpos = meta.dir ? -rpos : rpos;
	return((uint64_t)(rdpos<<32) | (uint64_t)pos.qpos);
}




/*
 * unittests for wrapped APIs
 */
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

unittest() {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_unused(q);
	dz_destroy(dz);
}

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
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(0, 0, 0, 0), "max(%u), cell(%u)", forefront->max, dz_ut_sel(0, 0, 0, 0));

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
#if 0
		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, dz_root(dz), 1, "", 0, 3);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(0, 0, 0, 0), "max(%u), cell(%u)", forefront->max, dz_ut_sel(0, 0, 0, 0));

		/* extend */
		forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("A", "\x0", "\x1", "M"), 1, 4);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 2, 5));
#endif
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

unittest( "calc_max.small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefronts[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_qpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_pos(dz, forefronts[0]) == ((2ULL<<32) | 2ULL));

	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
	ut_assert(dz_calc_max_rpos(dz, forefronts[1]) == 1);
	ut_assert(dz_calc_max_qpos(dz, forefronts[1]) == 3);
	ut_assert(dz_calc_max_pos(dz, forefronts[1]) == ((1ULL<<32) | 3ULL));

	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
	ut_assert(dz_calc_max_rpos(dz, forefronts[2]) == 4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[2]) == 7);
	ut_assert(dz_calc_max_pos(dz, forefronts[2]) == ((4ULL<<32) | 7ULL));

	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATG", "\x1\x0\x3\x2", "\x2\x1\x8\x4", "CKAM"), 4, 4);
	ut_assert(dz_calc_max_rpos(dz, forefronts[3]) == 3);
	ut_assert(dz_calc_max_qpos(dz, forefronts[3]) == 10);
	ut_assert(dz_calc_max_pos(dz, forefronts[3]) == ((3ULL<<32) | 10ULL));

	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
	ut_assert(dz_calc_max_rpos(dz, forefronts[4]) == 4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[4]) == 15);
	ut_assert(dz_calc_max_pos(dz, forefronts[4]) == ((4ULL<<32) | 15ULL));

	dz_destroy(dz);
}

unittest( "calc_max.small.revcomp" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefronts[5] = { NULL };

	/*
	 * AG---TTTT------CTGA
	 *   \ /    \    /
	 *    C      CATT
	 */
	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"CT"[2], &"\x1\x3"[2], &"\x2\x8"[2], &"AM"[2]), -2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefronts[0]) == -2, "%ld", dz_calc_max_rpos(dz, forefronts[0]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[0]) == 2);
	ut_assert(dz_calc_max_pos(dz, forefronts[0]) == ((((uint64_t)-2LL)<<32) | 2LL));

	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel(&"G"[1], &"\x2"[1], &"\x4"[1], &"T"[1]), -1, 2);
	ut_assert(dz_calc_max_rpos(dz, forefronts[1]) == -1, "%ld", dz_calc_max_rpos(dz, forefronts[1]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[1]) == 3);
	ut_assert(dz_calc_max_pos(dz, forefronts[1]) == ((((uint64_t)-1LL)<<32) | 3LL));

	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel(&"AAAA"[4], &"\x0\x0\x0\x0"[4], &"\x1\x1\x1\x1"[4], &"TQVL"[4]), -4, 3);
	ut_assert(dz_calc_max_rpos(dz, forefronts[2]) == -4, "%ld", dz_calc_max_rpos(dz, forefronts[2]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[2]) == 7);
	ut_assert(dz_calc_max_pos(dz, forefronts[2]) == ((((uint64_t)-4LL)<<32) | 7LL));

	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel(&"CATG"[4], &"\x1\x0\x3\x2"[4], &"\x2\x1\x8\x4"[4], &"MAKC"[4]), -4, 4);
	ut_assert(dz_calc_max_rpos(dz, forefronts[3]) == -3, "%lu", dz_calc_max_rpos(dz, forefronts[3]));
	ut_assert(dz_calc_max_qpos(dz, forefronts[3]) == 10);
	ut_assert(dz_calc_max_pos(dz, forefronts[3]) == ((((uint64_t)-3LL)<<32) | 10LL));

	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel(&"TCAG"[4], &"\x3\x1\x0\x2"[4], &"\x8\x2\x1\x4"[4], &"LTLQ"[4]), -4, 5);
	ut_assert(dz_calc_max_rpos(dz, forefronts[4]) == -4);
	ut_assert(dz_calc_max_qpos(dz, forefronts[4]) == 15);
	ut_assert(dz_calc_max_pos(dz, forefronts[4]) == ((((uint64_t)-4LL)<<32) | 15LL));

	dz_destroy(dz);
}

unittest( "calc_max.null" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	struct dz_forefront_s const *forefront = NULL;

	forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("TT", "\x3\x3", "\x8\x8", "KC"), 2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_qpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_pos(dz, forefront) == 0);

	forefront = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"AA"[2], &"\x0\x0"[2], &"\x1\x1"[2], &"KC"[2]), -2, 1);
	ut_assert(dz_calc_max_rpos(dz, forefront) == 0, "%ld", dz_calc_max_rpos(dz, forefront));
	ut_assert(dz_calc_max_qpos(dz, forefront) == 0);
	ut_assert(dz_calc_max_pos(dz, forefront) == 0);

	dz_destroy(dz);
}

/* short, exact matching sequences */
unittest( "trace" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefronts[5] = { NULL };
	dz_alignment_t const *aln = NULL;

	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("AG", "\x0\x2", "\x1\x4", "MA"), 2, 1);
	aln = dz_trace(dz, forefronts[0]);

	ut_assert(aln->score == dz_ut_sel(4, 4, 4, 9));
	ut_assert(aln->ref_length   == 2);
	ut_assert(aln->query_length == 2);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 2);
	ut_assert(strcmp((char const *)aln->path, "==") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 2);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "\x2", "T"), 1, 2);
	aln = dz_trace(dz, forefronts[1]);

	ut_assert(aln->score == dz_ut_sel(6, 6, 6, 14));
	ut_assert(aln->ref_length   == 3);
	ut_assert(aln->query_length == 3);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 2);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 3);
	ut_assert(strcmp((char const *)aln->path, "===") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 3);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "\x8\x8\x8\x8", "LVQT"), 4, 3);
	aln = dz_trace(dz, forefronts[2]);

	ut_assert(aln->score == dz_ut_sel(14, 14, 14, 32));
	ut_assert(aln->ref_length   == 7);
	ut_assert(aln->query_length == 7);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 3);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 7);
	ut_assert(strcmp((char const *)aln->path, "=======") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 7);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATG", "\x1\x0\x3\x2", "\x2\x1\x8\x4", "GKAM"), 4, 4);
	aln = dz_trace(dz, forefronts[3]);

	ut_assert(aln->score == dz_ut_sel(20, 20, 20, 47));
	ut_assert(aln->ref_length   == 10);
	ut_assert(aln->query_length == 10);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 4);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 10);
	ut_assert(strcmp((char const *)aln->path, "==========") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 10);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "\x2\x8\x4\x1", "QLTL"), 4, 5);
	aln = dz_trace(dz, forefronts[4]);

	ut_assert(aln->score == dz_ut_sel(25, 25, 25, 64));
	ut_assert(aln->ref_length   == 15);
	ut_assert(aln->query_length == 15);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 5);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);
	ut_assert(aln->span[4].id     == 5);
	ut_assert(aln->span[4].offset == 11);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

unittest( "trace.revcomp" ) {
	dz_t *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	dz_query_t *q = dz_pack_query(dz, dz_unittest_query, dz_unittest_query_length);
	dz_forefront_t const *forefronts[5] = { NULL };
	dz_alignment_t const *aln = NULL;

	forefronts[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel(&"CT"[2], &"\x1\x3"[2], &"\x2\x8"[2], &"AM"[2]), -2, 1);
	aln = dz_trace(dz, forefronts[0]);

	ut_assert(aln->score == dz_ut_sel(4, 4, 4, 9));
	ut_assert(aln->ref_length   == 2);
	ut_assert(aln->query_length == 2);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 1);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 2);
	ut_assert(strcmp((char const *)aln->path, "==") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 2);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel(&"G"[1], &"\x2"[1], &"\x4"[1], &"T"[1]), -1, 2);
	aln = dz_trace(dz, forefronts[1]);

	ut_assert(aln->score == dz_ut_sel(6, 6, 6, 14));
	ut_assert(aln->ref_length   == 3);
	ut_assert(aln->query_length == 3);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 2);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 3);
	ut_assert(strcmp((char const *)aln->path, "===") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 3);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel(&"AAAA"[4], &"\x0\x0\x0\x0"[4], &"\x1\x1\x1\x1"[4], &"TQVL"[4]), -4, 3);
	aln = dz_trace(dz, forefronts[2]);

	ut_assert(aln->score == dz_ut_sel(14, 14, 14, 32));
	ut_assert(aln->ref_length   == 7);
	ut_assert(aln->query_length == 7);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 3);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 7);
	ut_assert(strcmp((char const *)aln->path, "=======") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 7);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel(&"CATG"[4], &"\x1\x0\x3\x2"[4], &"\x2\x1\x8\x4"[4], &"MAKG"[4]), -4, 4);
	aln = dz_trace(dz, forefronts[3]);

	ut_assert(aln->score == dz_ut_sel(20, 20, 20, 47));
	ut_assert(aln->ref_length   == 10);
	ut_assert(aln->query_length == 10);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 4);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 10);
	ut_assert(strcmp((char const *)aln->path, "==========") == 0);

	ut_assert(aln->mismatch_count == 0);
	ut_assert(aln->match_count    == 10);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);



	forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel(&"TCAG"[4], &"\x3\x1\x0\x2"[4], &"\x8\x2\x1\x4"[4], &"LTLQ"[4]), -4, 5);
	aln = dz_trace(dz, forefronts[4]);

	ut_assert(aln->score == dz_ut_sel(25, 25, 25, 64));
	ut_assert(aln->ref_length   == 15);
	ut_assert(aln->query_length == 15);

	ut_assert(aln->span != NULL);
	ut_assert(aln->span_length    == 5);
	ut_assert(aln->span[0].id     == 1);
	ut_assert(aln->span[0].offset == 0);
	ut_assert(aln->span[1].id     == 2);
	ut_assert(aln->span[1].offset == 2);
	ut_assert(aln->span[2].id     == 3);
	ut_assert(aln->span[2].offset == 3);
	ut_assert(aln->span[3].id     == 4);
	ut_assert(aln->span[3].offset == 7);
	ut_assert(aln->span[4].id     == 5);
	ut_assert(aln->span[4].offset == 11);

	ut_assert(aln->path != NULL);
	ut_assert(aln->path_length    == 15);
	ut_assert(strcmp((char const *)aln->path, "==========X====") == 0);

	ut_assert(aln->mismatch_count == 1);
	ut_assert(aln->match_count    == 14);
	ut_assert(aln->ins_count      == 0);
	ut_assert(aln->del_count      == 0);

	(void)aln;
	dz_destroy(dz);
}

/* gap counts */
unittest( "trace.gapcount" ) {
	struct dz_s *dz = dz_init(dz_unittest_score_matrix, 5, 1, 20, 10);
	ut_assert(dz != NULL);

	char const *query = dz_ut_sel(
		"ACACTTCTAGACTTTACCACTA",
		"\x0\x1\x0\x1\x3\x3\x1\x3\x0\x2\x0\x1\x3\x3\x3\x0\x1\x1\x0\x1\x3\x0",
		"\x1\x2\x1\x2\x8\x8\x2\x8\x1\x4\x1\x2\x8\x8\x8\x1\x2\x2\x1\x2\x8\x1",
		"MSALLILALVGAAVAPLEDDDK"
	);
	struct dz_query_s const *q = dz_pack_query_forward(dz, query, 22);


	struct dz_forefront_s const *ff[10] = { 0 };
	ff[0] = dz_extend(dz, q, dz_root(dz), 1, dz_ut_sel("ACAC", "\x0\x1\x0\x1", "\x1\x2\x1\x2", "MSAL"), 4, 0);

	ff[1] = dz_extend(dz, q, &ff[0], 1, dz_ut_sel("TTGT", "\x3\x3\x2\x3", "\x8\x8\x4\x8", "LILA"), 4, 1);
	ff[2] = dz_extend(dz, q, &ff[0], 1, dz_ut_sel("ATCC", "\x0\x3\x1\x1", "\x1\x8\x2\x2", "KKLG"), 4, 2);
	ff[3] = dz_extend(dz, q, &ff[1], 2, dz_ut_sel("AGAC", "\x0\x2\x0\x1", "\x1\x4\x1\x2", "LVGA"), 4, 3);

	ff[4] = dz_extend(dz, q, &ff[3], 1, dz_ut_sel("T", "\x3", "\x8", "A"), 1, 4);
	ff[5] = dz_extend(dz, q, &ff[3], 2, dz_ut_sel("TTCTA", "\x3\x3\x1\x3\x0", "\x8\x8\x2\x8\x1", "AVAFP"), 5, 5);

	ff[6] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("A", "\x0", "\x1", "A"), 1, 6);
	ff[7] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("C", "\x1", "\x2", "E"), 1, 6);
	ff[8] = dz_extend(dz, q, &ff[5], 1, dz_ut_sel("G", "\x2", "\x4", "M"), 1, 6);
	ff[9] = dz_extend(dz, q, &ff[6], 3, dz_ut_sel("CACGG", "\x1\x0\x1\x2\x2", "\x2\x1\x2\x4\x4", "EDDCL"), 5, 9);

	/* detect max */
	struct dz_forefront_s const *max = NULL;
	for(size_t i = 0; i < 10; i++) {
		if(max == NULL || ff[i]->max > max->max) { max = ff[i]; }
	}
	ut_assert(max == ff[9]);

	/* traceback */
	struct dz_alignment_s const *aln = dz_trace(dz, max);
	ut_assert(aln->ref_length     == dz_ut_sel(23, 23, 23, 23), "%u", aln->ref_length);
	ut_assert(aln->query_length   == dz_ut_sel(22, 22, 22, 22), "%u", aln->query_length);
	ut_assert(aln->score == dz_ut_sel(33, 33, 33, 88));
	ut_assert(strcmp((char const *)aln->path, dz_ut_sel("======X=======I======XX", "======X=======I======XX", "======X=======I======XX", "===============I=X===XX")) == 0, "%s", aln->path);

	ut_assert(aln->path_length    == dz_ut_sel(23, 23, 23, 23), "%u", aln->path_length);
	ut_assert(aln->span_length    == 6, "%u", aln->span_length);

	ut_assert(aln->span[0].id     == 0, "%u", aln->span[0].id);
	ut_assert(aln->span[0].offset == 0, "%u", aln->span[0].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[0].offset], dz_ut_sel("====", "====", "====", "===="), 4) == 0);

	ut_assert(aln->span[1].id     == 1, "%u", aln->span[1].id);
	ut_assert(aln->span[1].offset == 4, "%u", aln->span[1].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[1].offset], dz_ut_sel("==X=", "==X=", "==X=", "===="), 4) == 0);

	ut_assert(aln->span[2].id     == 3, "%u", aln->span[2].id);
	ut_assert(aln->span[2].offset == 8, "%u", aln->span[2].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[2].offset], dz_ut_sel("====", "====", "====", "===="), 4) == 0);

	ut_assert(aln->span[3].id     == 5, "%u", aln->span[3].id);
	ut_assert(aln->span[3].offset == 12, "%u", aln->span[3].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[3].offset], dz_ut_sel("==I==", "==I==", "==I==", "===I="), 5) == 0);

	ut_assert(aln->span[4].id     == 6, "%u", aln->span[4].id);
	ut_assert(aln->span[4].offset == 17, "%u", aln->span[4].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[4].offset], dz_ut_sel("=", "=", "=", "X"), 1) == 0);

	ut_assert(aln->span[5].id     == 9, "%u", aln->span[5].id);
	ut_assert(aln->span[5].offset == 18, "%u", aln->span[5].offset);
	ut_assert(strncmp((char const *)&aln->path[aln->span[5].offset], dz_ut_sel("===XX", "===XX", "===XX", "===XX"), 5) == 0);

	dz_destroy(dz);
}


static __dz_force_inline
uint8_t dz_unittest_random_char(void)
{
	#ifndef DZ_PROTEIN
		switch(rand() % 4) {
			case 0:  return(dz_ut_sel('A', 0, 1, 0));
			case 1:  return(dz_ut_sel('C', 1, 2, 0));
			case 2:  return(dz_ut_sel('G', 2, 4, 0));
			case 3:  return(dz_ut_sel('T', 3, 8, 0));
			default: return(dz_ut_sel('A', 0, 1, 0));
		}
	#else
		return(rand() % 20);
	#endif
}

typedef struct {
	char *seq;
	size_t len;
} dz_unittest_seq_t;

static __dz_force_inline
dz_unittest_seq_t dz_unittest_rand_seq(size_t len)
{
	char *seq = malloc(sizeof(char) * (len + 1));
	for(size_t i = 0; i < len; i++) {
		seq[i] = dz_unittest_random_char();
	}
	seq[len] = '\0';

	return((dz_unittest_seq_t){
		.seq = seq,
		.len = len
	});
}

static __dz_force_inline
dz_unittest_seq_t dz_unittest_mod_seq(char const *seq, size_t len)
{
	char *mod = malloc(sizeof(char) * (4 * len + 1024));
	size_t j = 0;
	for(size_t i = 0; i < len; i++) {
		switch(rand() % 20) {
			case 0: break;	/* del */
			case 1:			/* ins */
				mod[j++] = dz_unittest_random_char();
				mod[j++] = seq[i];
				break;
			// case 2:
				// mod[j++] = dz_unittest_random_char();
				// break;
			default:
				mod[j++] = seq[i];
				break;
		}
	}
	mod[j] = '\0';
	return((dz_unittest_seq_t){
		.seq = mod,
		.len = j
	});
}


unittest( "long" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	size_t const len = 16384;
	dz_unittest_seq_t query = dz_unittest_rand_seq(len);
	struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
	ut_assert(q != NULL);

	struct dz_forefront_s const *ff = dz_extend(dz, q, dz_root(dz), 1, query.seq, query.len, 0);
	ut_assert(ff != NULL);
	#ifndef DZ_PROTEIN
		/* FIXME: accumulate score for protein */
		ut_assert(ff->max == (int32_t)(dz_unittest_score_matrix[0] * query.len));
	#endif

	struct dz_alignment_s const *aln = dz_trace(dz, ff);
	ut_assert(aln != NULL);

	free(query.seq);
	dz_destroy(dz);
}

unittest( "long.gaps" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t i = 0; i < 256; i++) {
		size_t const len = 16384;
		dz_unittest_seq_t query = dz_unittest_rand_seq(len);
		struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
		ut_assert(q != NULL);

		dz_unittest_seq_t ref = dz_unittest_mod_seq(query.seq, query.len);
		struct dz_forefront_s const *ff = dz_extend(dz, q, dz_root(dz), 1, ref.seq, ref.len, 0);
		ut_assert(ff != NULL);
		/* FIXME? we don't exactly know how large the score is */

		uint64_t pos = dz_calc_max_pos(dz, ff);
		ut_assert(pos != 0);
		/* FIXME? we don't exactly know where the maxpos is */

		struct dz_alignment_s const *aln = dz_trace(dz, ff);
		ut_assert(aln != NULL);		/* we only expect the trace back is done without failure */

		free(query.seq);
		free(ref.seq);
	}
	dz_destroy(dz);
}

static __dz_force_inline
dz_unittest_seq_t dz_unittest_sample_seq(dz_unittest_seq_t *seq, size_t len, size_t cnt)
{
	size_t const rlen = 1 * 1024 * 1024;
	size_t const blen = cnt * len + 2 * rlen;
	char *buf = malloc(sizeof(char) * blen);

	/*
	 * compose graph as following:
	 *
	 *   B---D---F---H     X
	 *  / \ / \ / \ /       \
	 * A   X   X   X   ...   Z
	 *  \ / \ / \ / \       /
	 *   C---E---G---I     Y
	 *
	 * and build a seq by sampling either top or bottom nodes.
	 */
	size_t bpos = 0;
	memcpy(&buf[bpos], seq[0].seq, seq[0].len); bpos += seq[0].len;
	for(size_t i = 1; i < cnt - 1; i += 2) {
		size_t j = rand() & 0x01ULL;
		memcpy(&buf[bpos], seq[i + j].seq, seq[i + j].len); bpos += seq[i + j].len;
	}

	/* clip tail */
	size_t const tlen = seq[cnt - 1].len - dz_min2(seq[cnt - 1].len - 1, 8192);
	memcpy(&buf[bpos], seq[cnt - 1].seq, tlen); bpos += tlen;

	/* add random bases */
	dz_unittest_seq_t rseq = dz_unittest_rand_seq(rlen);
	memcpy(&buf[bpos], rseq.seq, rseq.len); bpos += rseq.len;
	free(rseq.seq);

	/* cap */
	buf[bpos] = '\0';

	#if 0
		/* put gaps and substitutions */
		dz_unittest_seq_t m = dz_unittest_mod_seq(buf, bpos);

		free(buf);
		return(m);
	#else
		return((dz_unittest_seq_t){
			.seq = buf,
			.len = bpos
		});
	#endif
}

unittest( "long.graph" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	size_t const cnt = 128;
	dz_unittest_seq_t ref[cnt];
	memset(ref, 0, cnt * sizeof(dz_unittest_seq_t));

	size_t const len = 4 * 1024 * 1024;
	for(size_t i = 0; i < cnt; i++) {
		ref[i] = dz_unittest_rand_seq((i % 3) == 0 ? 16 : len);
	}

	dz_unittest_seq_t query = dz_unittest_sample_seq(ref, len, cnt);
	struct dz_query_s const *q = dz_pack_query_forward(dz, query.seq, query.len);
	ut_assert(q != NULL);

	/* root */
	struct dz_forefront_s const *ff[cnt];
	ff[0] = dz_extend(dz, q, dz_root(dz), 1, ref[0].seq, ref[0].len, 0);
	ut_assert(ff[0] != NULL);
	ut_assert(ff[0]->max > 0, "%d", ff[0]->max);		/* score positive */

	/* extend over net graph */
	int32_t prev_max = ff[0]->max;
	for(size_t i = 1; i < cnt - 1; i += 2) {
		struct dz_forefront_s const *incoming[2];
		size_t icnt = 0;
		if(i == 1) {
			incoming[icnt++] = ff[0];
		} else {
			if(ff[i - 2] != NULL) { incoming[icnt++] = ff[i - 2]; }
			if(ff[i - 1] != NULL) { incoming[icnt++] = ff[i - 1]; }
		}

		debug("first");
		ff[i] = dz_extend(dz, q, incoming, icnt, ref[i].seq, ref[i].len, i);
		ut_assert(ff[i] != NULL);

		debug("second");
		ff[i + 1] = dz_extend(dz, q, incoming, icnt, ref[i + 1].seq, ref[i + 1].len, i + 1);
		ut_assert(ff[i + 1] != NULL);

		/* either is unterminated */
		ut_assert(!dz_is_terminated(ff[i]) || !dz_is_terminated(ff[i + 1]));

		/* score increasing */
		ut_assert(dz_max2(ff[i]->max, ff[i + 1]->max) > prev_max, "i(%zu), (%d, %d), %d", i, ff[i]->max, ff[i + 1]->max, prev_max);
		prev_max = dz_max2(ff[i]->max, ff[i + 1]->max);
	}

	/* tail */
	ff[cnt - 1] = dz_extend(dz, q,
		&ff[cnt - 3], 2,
		ref[cnt - 1].seq, ref[cnt - 1].len, cnt - 1
	);
	ut_assert(ff[cnt - 1] != NULL);
	ut_assert(ff[cnt - 1]->max > prev_max);
	ut_assert(ff[cnt - 1]->cap != NULL);
	/* we don't exactly know how large the score is */

	struct dz_alignment_s const *aln = dz_trace(dz, ff[cnt - 1]);
	ut_assert(aln != NULL);		/* we only expect the trace back is done without failure */


	free(query.seq);
	for(size_t i = 0; i < cnt; i++) {
		free(ref[i].seq);
	}
	dz_destroy(dz);
}



#endif	/* defined(UNITTEST) && UNITTEST != 0 */
#endif	/* DZ_WRAPPED_API */



#ifdef __cplusplus
};	/* extern "C" { */
#endif

/**
 * end of dozeu.h
 */
