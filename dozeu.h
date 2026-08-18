/*
 * Copyright (c) 2026, NVIDIA CORPORATION. All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 */

// $(CC) -O3 -march=native -DMAIN -o dozeu dozeu.c
//#ifdef DZ_QUAL_ADJ
//#define DEBUG
//#endif
//#define DZ_PRINT_VECTOR
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

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/x86/sse4.1.h"

// Only turn on extern "C" after including code that can use C++ features.
#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>

#ifndef DZ_INCLUDE_ONCE

#define DZ_NUM_QUAL_SCORES 64
#define DZ_QUAL_MATRIX_SIZE (16 * DZ_NUM_QUAL_SCORES)

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

#if defined(DZ_NUCL_2BIT) && defined(DZ_QUAL_ADJ)
#error "Base quality adjustment is only supported in ASCII configuration"
#endif

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
#ifdef DZ_FULL_LENGTH_BONUS
#define dz_end_bonus(_self, _query, _i)            ((_i) / 8 == (_query)->blen - 1 ? (_query)->bonus[8 + ((_i) & 7)] : 0)
#else
#define dz_end_bonus(_self, _query, _i)            0;
#endif

/* get the first index of the quals from a packed query (only valid if DZ_QUAL_ADJ is defined) */
#define dz_quals(_query)                           ( (uint8_t const *) &(_query)->arr[(_query)->arr_end] )
/* get the base quality adjusted matrix (only valid if DZ_QUAL_ADJ is defined) */
#define dz_qual_matrix(_self)                      ( (int8_t const *)((_self) + 1) )

#define dz_pair_score(_self, _q, _r, _i)	       ( (_self)->matrix[((_r) | (_q)->arr[(_i)]) & 0x1f] + dz_end_bonus(_self, _q, _i))

#define dz_qual_adj_pair_score(_self, _q, _r, _i)  ( dz_qual_matrix(_self)[(((uint32_t) dz_quals(_q)[(_i)]) << 5) | ((uint32_t)(_r)) | ((uint32_t)((_q)->arr[(_i)] & 0x1f))] + dz_end_bonus(_self, _q, _i) )

#define dz_pair_eq(_self, _q, _r, _i)		       ( (uint32_t)((_q)->arr[(_i)]) == ((uint32_t)(_r)<<2) )

#else // ! DZ_PROTEIN

/* protein */
#  ifndef DZ_MAT_SIZE
#    define DZ_MAT_SIZE				( 32 )
#  endif
#ifdef DZ_FULL_LENGTH_BONUS
#define dz_end_bonus(_self, _query, _i)            ((_i) / 8 == (_query)->blen - 1 ? (_query)->bonus[8 + ((_i) & 7)] : 0)
#else
#define dz_end_bonus(_self, _query, _i)            0;
#endif
#define dz_pair_score(_self, _q, _r, _i)	( (int8_t)((_q)->arr[(_r) * (_q)->blen * L + (_i)]) + dz_end_bonus(_self, _q, _i))
#define dz_pair_eq(_self, _q, _r, _i)		( (uint32_t)((_q)->q[(_i) - 1] & 0x1f) == (uint32_t)(_r) )
#endif


/* utils */
#define dz_pp_cat_intl(x, y)		x##y
#define dz_pp_cat(x, y)				dz_pp_cat_intl(x, y)
#define dz_static_assert(expr)		typedef char dz_pp_cat(_st_, __LINE__)[(expr) ? 1 : -1]
#define dz_trap()					{ *((volatile uint8_t *)NULL); }
#ifdef DZ_PROTEIN
	dz_static_assert(DZ_MAT_SIZE <= 32);
#endif // DZ_PROTEIN

#if (defined(DEBUG) || defined(UNITTEST))
#  include "log.h"
#else
#  define debug(...)                ;
#endif
                     
#ifdef UNITTEST
#  define UNITTEST_ALIAS_MAIN		0
#  define UNITTEST_UNIQUE_ID		3213
#  include "unittest.h"
unittest_config( "dozeu" );
unittest() { debug("hello"); }
#else
#  define unittest(...)				static void dz_pp_cat(dz_unused_, __LINE__)(void)
#  define ut_assert(...)			;
#  define trap()					;
#endif

#define dz_die(...) { \
	dz_die_impl(__VA_ARGS__, ""); \
}
#define dz_die_impl(fmt, ...) { \
	fprintf(stderr, "[%s: %s(%d)] " fmt "%s\n", __FILE__, __func__, __LINE__, __VA_ARGS__); \
    exit(1); \
}

#define dz_error                    dz_die

#if defined(DZ_NUCL_ASCII)
#  define DZ_UNITTEST_INDEX			0
#elif defined(DZ_NUCL_2BIT)
#  define DZ_UNITTEST_INDEX			1
#elif defined(DZ_PROTEIN)
#  define DZ_UNITTEST_INDEX			2
#endif
#define dz_ut_sel(a, b, c)			( (DZ_UNITTEST_INDEX == 0) ? (a) : ((DZ_UNITTEST_INDEX == 1) ? (b) : (c)) )

/* vectorize */
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

#define dz_is_all_zero(x)			( _mm_test_all_zeros((x), (x)) == 1 )

#define DZ_MEM_MARGIN_SIZE			( 256 )
#define DZ_MEM_ALIGN_SIZE			( 16 )
#define DZ_MEM_INIT_SIZE			( 16 * 1024 * 1024 )

#define DZ_CELL_MIN					( INT16_MIN )
#define DZ_CELL_MAX					( INT16_MAX )
#define DZ_CELL_MARGIN				( 32 )
#define DZ_CELL_MARGINED_MIN		( DZ_CELL_MIN + DZ_CELL_MARGIN )
#define DZ_CELL_MARGINED_MAX		( DZ_CELL_MAX - DZ_CELL_MARGIN )
                     
/* query; preconverted query sequence; blen = roundup(qlen, L) / L; array must have 16-byte-length margin at the tail */
struct dz_query_s {
    uint32_t blen;
    uint32_t arr_end;
    char const *q;
    int16_t bonus[16];
    uint8_t arr[];
};
dz_static_assert(sizeof(struct dz_query_s) % sizeof(__m128i) == 0);

/* node (reference) */
struct dz_node_s { int32_t id, len; uint8_t const *ptr; };
dz_static_assert(sizeof(struct dz_node_s) % sizeof(__m128i) == 0);

/* DP matrix structures */

/*
 * The Dozeu DP matrix for a node looks like:
 * - A head column, containing
 *		- An array of pointers to forefronts for preceeding nodes
 *		- A head cap (same size as a normal cap)
 *	- 0 or more internal columns, which are:
 *		- An array slice of score vector items
 *		- A cap, containing:
 *			- The range of the array slice
 *			- Other data
 *	- A final column, which is:
 *		- An array slice of score vector items
 *		- A forefront (a special final cap), containing
 *			- The range of the array slice
 *			- Other data
 *
 * So, each slice has the cap for the preceeding slice before it, and its own
 * range stored immediately after it. This is because we dynamically stop
 * slices early, and we put the final range actually used once we know it.
 *
 * We are able to break this structure in memory! If the next column's slice
 * and cap could be too big to fit right after the previous column in the
 * contiguous memory available, we start a new head column with the previous
 * column's cap as a forefront, and all the forefront data except the cap left
 * uninitialized. We call this an "internal bridge".
 */

/* Score vector item. Array of these will be followed by a dz_cap_s, which leads with a range. */
struct dz_swgv_s { __m128i e, f, s; };
/* placed just after every score vector (as part of a cap) to indicate the length */
struct dz_range_s { uint32_t spos, epos; };

/* A special kind of initial cap. Identical in size to a cap. */
struct dz_head_s {
	struct dz_range_s r;
	uint32_t rch, n_forefronts;
	struct dz_range_s fr;
	uint32_t _pad[2];
};

/* followed by dz_forefront_s; range spos and epos are "shared to" (???) forefront_s */
struct dz_cap_s {
	struct dz_range_s r;
	uint32_t rch; int32_t rrem;
	struct dz_range_s fr; // the range that should continue forward
	uint32_t _pad[2];
};
/*
 * A forefront. Appears at the end of a matrix.
 *
 * Is also a special kind of cap.
 */
struct dz_forefront_s {
	struct dz_range_s r; // the range that applies to the preceding column
	uint32_t rid;
	int32_t rlen;
    struct dz_range_s fr; // the range that should continue forward
	uint32_t rsum, rcnt; int32_t max, inc;
    uint8_t _pad[8];
	struct dz_query_s const *query;
	struct dz_cap_s const *mcap;
};
/*
 * When initializing forefronts, we want to be able to fill in the pad without saying 0s everywhere.
 */
#define DZ_FOREFRONT_PAD {0, 0, 0, 0, 0, 0, 0, 0}

/*
 * This holds the forefront for the alignment root, and the x-drop threshold,
 * computed from the max gap length.
 */
struct dz_alignment_init_s {
    struct dz_forefront_s const *root;
    uint16_t xt;
};
                     
dz_static_assert(sizeof(struct dz_swgv_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_cap_s) % sizeof(__m128i) == 0);
dz_static_assert(sizeof(struct dz_forefront_s) % sizeof(__m128i) == 0);
/* _unwind_cap steps back over a head cap by sizeof(struct dz_cap_s) */
dz_static_assert(sizeof(struct dz_head_s) == sizeof(struct dz_cap_s));
/* the traceback reads fr through head, cap and forefront pointers interchangeably */
dz_static_assert(offsetof(struct dz_head_s, fr) == offsetof(struct dz_cap_s, fr));
dz_static_assert(offsetof(struct dz_cap_s, fr) == offsetof(struct dz_forefront_s, fr));
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

/*
 * Header for a block of memory used in an arena.
 *
 * Exists at the beginning of the block.
 */
struct dz_mem_block_s { struct dz_mem_block_s *next; size_t size; };
/*
 * Bookkeeping information for the arena allocator.
 *
 * Exists at the beginning of the first block, just after its block header.
 */
struct dz_stack_s { struct dz_mem_block_s *curr; uint8_t *top, *end; uint64_t being_allocated; uint64_t reservation_remaining; uint8_t *dz_flush_top; };
/*
 * Represents a Dozeu memory allocation arena.
 *
 * Exists at the beginning of the arena's first block, and contains its block
 * header and also its stack bookkeeping information.
 */
struct dz_mem_s { struct dz_mem_block_s blk; struct dz_stack_s stack; };
/* Compute the number of bytes available to allocate in the current active block of an arena. */
#define dz_mem_stack_rem(_mem)		( (size_t)((_mem)->stack.end - (_mem)->stack.top) )

struct dz_s {
    int8_t matrix[32];
    uint16_t giv[8], gev[8], riv[8], rev[8]; // padding ensures memory alignment
    int8_t protein_matrix[];
};
dz_static_assert(sizeof(struct dz_s) % sizeof(__m128i) == 0);
#define dz_mem(_self)				( (struct dz_mem_s *)(_self) - 1 )
/* get the base quality adjusted matrix (only valid if DZ_QUAL_ADJ is defined) */
#define dz_qual_matrix(_self)       ( (int8_t const *)((_self) + 1) )

//#define dz_root(_self)				( (struct dz_forefront_s const **)(&_self->root) )
#define dz_is_terminated(_ff)		( dz_cff(_ff)->r.spos >= dz_cff(_ff)->r.epos )
#define dz_gt(_ff)					( dz_cff(_ff)->inc > 0 )
#define dz_geq(_ff)					( dz_cff(_ff)->mcap != NULL )

#ifdef DZ_PRINT_VECTOR
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
#define print_vector8(v) { \
debug("%s (%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d)", #v, \
    (int)(int8_t)_mm_extract_epi8(v, 15), \
    (int)(int8_t)_mm_extract_epi8(v, 14), \
    (int)(int8_t)_mm_extract_epi8(v, 13), \
    (int)(int8_t)_mm_extract_epi8(v, 12), \
    (int)(int8_t)_mm_extract_epi8(v, 11), \
    (int)(int8_t)_mm_extract_epi8(v, 10), \
    (int)(int8_t)_mm_extract_epi8(v, 9), \
    (int)(int8_t)_mm_extract_epi8(v, 8), \
    (int)(int8_t)_mm_extract_epi8(v, 7), \
    (int)(int8_t)_mm_extract_epi8(v, 6), \
    (int)(int8_t)_mm_extract_epi8(v, 5), \
    (int)(int8_t)_mm_extract_epi8(v, 4), \
    (int)(int8_t)_mm_extract_epi8(v, 3), \
    (int)(int8_t)_mm_extract_epi8(v, 2), \
    (int)(int8_t)_mm_extract_epi8(v, 1), \
    (int)(int8_t)_mm_extract_epi8(v, 0)); \
}
#else
#define print_vector(v) ;
#define print_vector8(v) ;
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
                     
#endif // DZ_INCLUDE_ONCE
                     
                     
#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)

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
                     
#endif

#ifndef DZ_INCLUDE_ONCE
                     
/**
 * @fn dz_mem_init, dz_mem_destroy, dz_mem_add_stack, dz_mem_malloc, dz_mem_flush
 * @brief stack chain
 */

/*
 * "Flush" a memory arena.
 *
 * Retain all blocks allocated from the system allcoator, but internally mark
 * all the space as free.
 */
static __dz_force_inline
void dz_mem_flush(
	struct dz_mem_s *mem)
{
	mem->stack.curr = &mem->blk;
	mem->stack.top = (uint8_t *)dz_roundup((uintptr_t)(mem + 1), DZ_MEM_ALIGN_SIZE);
    // TODO: i think this margin is redundant, since dz_malloc adds the margin past the end
    // of the requested size
    mem->stack.end = (uint8_t *)mem + mem->blk.size;// - DZ_MEM_MARGIN_SIZE;
	mem->stack.being_allocated = 0;
    mem->stack.reservation_remaining = 0;
	return;
}

/*
 * Create a new memory arena.
 *
 * The arena will initially have space for internal allocation structures, and
 * at least the given number of bytes of user data.
 *
 * If allocation fails, NULL is returned.
 */
static __dz_force_inline
struct dz_mem_s *dz_mem_init(
	size_t size)
{
	size = dz_max2(DZ_MEM_INIT_SIZE, sizeof(struct dz_mem_s) + dz_roundup(size, DZ_MEM_ALIGN_SIZE));
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

/*
 * Free a memory arena created with dz_mem_init().
 *
 * The arena must not be NULL.
 */
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

/*
 * Advance a memory arena to a next block that can allocate the given number of bytes.
 *
 * Always advances to a new block; assumes sufficient space is definitely not
 * available in the current block.
 *
 * Handles creating new blocks for the arena, or replacing old blocks that are
 * not large enough.
 *
 * Returns 0 on success, and 1 if memory could not be allocated.
 */
static __dz_force_inline
uint64_t dz_mem_add_stack(
	struct dz_mem_s *mem,
	size_t size)
{
    debug("add_stack, ptr(%p)", mem->stack.curr->next);
	if(mem->stack.curr->next == NULL
       || dz_unlikely(mem->stack.curr->next->size < size + dz_roundup(sizeof(struct dz_mem_block_s), DZ_MEM_ALIGN_SIZE))) {
        if (mem->stack.curr->next != NULL) {
            /* didn't allocate enough memory earlier, throw out all subsequent blocks so we can start again bigger */
            struct dz_mem_block_s *blk = mem->stack.curr->next;
            while (blk != NULL) {
                struct dz_mem_block_s *next = blk->next;
                dz_free(blk);
                blk = next;
            }
            mem->stack.curr->next = NULL;
        }
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

/*
 * Dump stack information to standard error.
 */
static __dz_force_inline
void dz_mem_log_stack(
    struct dz_mem_s *mem)
{
    fprintf(stderr, "Dozeu stack blocks:\n");
    struct dz_mem_block_s *blk = &(mem->blk);
    while (blk != NULL) {
        fprintf(stderr, "\t%p\t%lu\n", blk, blk->size);
        blk = blk->next;
    }
}

/*
 * Allocate memory from an arena.
 *
 * Returns NULL if memory could not be allocated.
 */
static __dz_force_inline
void *dz_mem_malloc(
	struct dz_mem_s *mem,
	size_t size)
{
	/* Make sure to maintain memory alignment.
	 *
	 * We need to precompute the aligned size instead of just bumping top by
	 * the aligned size, because we need to make sure the alignment padding
	 * memory is actually available in the space between top and end.
	 * Otherwise, top could pass end, and dz_mem_stack_rem() could overflow.
	 */
	size = dz_roundup(size, sizeof(__m128i));

    if(dz_mem_stack_rem(mem) < size) {
        if(dz_mem_add_stack(mem, size)) {
			/* Report a failed allocation. */
			dz_die("Allocation of %lu bytes failed", size);
            return NULL;
		}
    }
	void *ptr = (void *)mem->stack.top;
	mem->stack.top += size;
    return(ptr);
}

/*
 * Prepare a contiguous run of memory from an arena.
 *
 * Returns 0 if successful, and 1 on failure.
 */
static __dz_force_inline
uint64_t dz_mem_stream_reserve(
	struct dz_mem_s *mem,
	size_t size)
{
	debug("Reserve %lu bytes for stream", size);
    mem->stack.reservation_remaining = size;
	/* We cheat and just use a whole block as a run. */
	if(dz_likely(dz_mem_stack_rem(mem) < size)) {
		if (dz_mem_add_stack(mem, size)) {
            dz_error("Could not reserve %lu bytes", size);
            return 1;
        }
        return 0;
	}
	return 0;
}

/*
 * Get the amount of stream reservation remaining.
 */
static __dz_force_inline
uint64_t dz_mem_stream_remaining(
	struct dz_mem_s *mem)
{
	return mem->stack.reservation_remaining;
}

/*
 * Begin an allocation of up to the given number of bytes from the current
 * stream reservation.
 *
 * Returns a pointer to them, or NULL if they do not fit.
 */
static __dz_force_inline
void *dz_mem_stream_alloc_begin(
	struct dz_mem_s *mem,
	size_t size)
{
	if(dz_mem_stack_rem(mem) < size || mem->stack.being_allocated != 0) {
        if (dz_mem_stack_rem(mem) < size) {
		    dz_error("Asked to allocate %lu bytes, but only %lu bytes remain on stack", size, dz_mem_stack_rem(mem));
        } else if (mem->stack.being_allocated != 0) {
            dz_error("Asked to allocate %lu bytes, but another allocaton of %lu bytes is already in progress", size, mem->stack.being_allocated);
        }
        return NULL;
	}
    if (mem->stack.reservation_remaining < size) {
        dz_error("Asked to allocate %lu bytes, but only %lu bytes remain in reservation", size, mem->stack.reservation_remaining);
        return NULL;
    }
	void *ptr = (void *)mem->stack.top;
	debug("Begin allocating %lu bytes (%p) from %lu in reservation", size, ptr, mem->stack.reservation_remaining);
	/* Reserve this area above the stack */
	mem->stack.being_allocated = size;
    return(ptr);
}

/*
 * Get the start address of the active allocation.
 *
 * Returns NULL if no allocation is active.
 */
static __dz_force_inline
void *dz_mem_stream_alloc_current(
	struct dz_mem_s *mem)
{
	if(mem->stack.being_allocated == 0) {
		dz_error("Nothing currently being allocated");
        return NULL;
	}
	return (void *)mem->stack.top;
}


/*
 * End an active allocation, having used the given number of bytes.
 *
 * Returns 0 on success, or 1 if too many bytes were used.
 */
static __dz_force_inline
uint64_t dz_mem_stream_alloc_end(
	struct dz_mem_s *mem,
	size_t size)
{
	debug("Finish allocating %lu bytes (%p)", size, mem->stack.top);
    if(size > mem->stack.being_allocated) {
		/* Too much memory was used vs. what we set aside for the stream. */
		debug("That was too many");
		return 1;
	}
	mem->stack.top += size;
	mem->stack.being_allocated = 0;
    mem->stack.reservation_remaining -= size;
	return 0;
}

/*
 * Allocate the given number of bytes from the current stream reservation.
 *
 * Returns a pointer to them, or NULL if they do not fit.
 */
static __dz_force_inline
void *dz_mem_stream_alloc(
	struct dz_mem_s *mem,
	size_t size)
{
	/* TODO: Optimize instead of using the async primitives. */
	void *ptr = dz_mem_stream_alloc_begin(mem, size);
	if(dz_mem_stream_alloc_end(mem, size)) {
		debug("Stream alloc failed");
		ptr = NULL;
	}
	return ptr;
}

/*
 * Allocate the given number of items of the given size from the current
 * reservation, right-justified, in a block aligned to the given alignemnt.
 *
 * Returns a pointer to them, or NULL if they do not fit.
 */
static __dz_force_inline
void *dz_mem_stream_right_alloc(
	struct dz_mem_s *mem,
	size_t items,
	size_t size,
	size_t alignment)
{
	size_t data_size = size * items;
	size_t alloc_size = dz_roundup(data_size, alignment);

	void *ptr = dz_mem_stream_alloc(mem, alloc_size);
	if(ptr == NULL) {
		debug("Stream right alloc failed");
		return ptr;
	}
	void* right_aligned_ptr = (void*)((uint8_t*)ptr + (alloc_size - data_size));
	debug("Allocated %lu items of %lu bytes each padded to %lu: block %p and right-aligned array %p", items, size, alignment, ptr, right_aligned_ptr);
	return right_aligned_ptr;
}
                     
#endif // DZ_INCLUDE_ONCE

#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)
                     
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
                     
#endif

/**
 * vector update macros
 */

#define _calc_max_slice_size(_sp, _ep) ({ \
    size_t max_slice_size = (2 * ((_ep) - (_sp)) * sizeof(struct dz_swgv_s)); \
    debug("sp(%lu), ep(%lu), max_slice_size(%lu)", (size_t) (_sp), (size_t) (_ep), max_slice_size); \
    max_slice_size; \
})
/*
 * Determine how many bytes are needed to store forefront pointers, the head
 * cap, the column vector, and the cap at the end from _end_column, which may
 * need to be a (larger) forefront if this ends up being the last column.
 * _sp: start position in the column
 * _ep: end position in the column
 * _nt: number of forefronts that need pointers
 */
#define _calc_next_size(_sp, _ep, _nt) ({ \
	size_t forefront_arr_size = dz_roundup(sizeof(struct dz_forefront_s *) * (_nt), sizeof(__m128i)); \
	size_t est_column_size = _calc_max_slice_size(_sp, _ep); \
    size_t head_cap_size = sizeof(struct dz_head_s); \
    size_t end_cap_size = sizeof(struct dz_forefront_s); \
	size_t next_req = forefront_arr_size + head_cap_size + est_column_size + end_cap_size; \
	debug("sp(%lu), ep(%lu), nt(%lu), forefront_arr_size(%lu), head_cap_size(%lu), est_column_size(%lu), end_cap_size(%lu), next_req(%lu)", (size_t) (_sp), (size_t) (_ep), (size_t) (_nt), forefront_arr_size, head_cap_size, est_column_size, end_cap_size, next_req); \
	next_req; \
})
/*
 * Save an array of forefront pointers to the stack, followed by a head cap
 * that indicates the number of forefronts before it.
 *
 * Space must have been reserved in the arena.
 */
#define _init_cap(_adj, _rch, _forefronts, _n_forefronts) ({ \
	/* push forefront pointers */ \
	size_t forefront_arr_size = dz_roundup(sizeof(struct dz_forefront_s *) * (_n_forefronts), sizeof(__m128i)); \
    struct dz_forefront_s const **dst = (struct dz_forefront_s const **)(dz_mem_stream_right_alloc(dz_mem(self), (int64_t)(_n_forefronts), sizeof(struct dz_forefront_s *), sizeof(__m128i))); \
	struct dz_forefront_s const **src = (struct dz_forefront_s const **)(_forefronts); \
	for(size_t i = 0; i < (_n_forefronts); i++) { dst[i] = src[i]; } \
	/* push head-cap info */ \
	struct dz_head_s *_head = (struct dz_head_s *)(dz_mem_stream_alloc(dz_mem(self), sizeof(struct dz_head_s))); \
	(_head)->r.spos = (_adj);					/* save merging adjustment */ \
	(_head)->r.epos = 0;						/* head marked as zero */ \
	(_head)->rch = (_rch);						/* rch for the first column */ \
	(_head)->n_forefronts = (_n_forefronts);	/* record n_forefronts */ \
	(_head)->fr.spos = 0; (_head)->fr.epos = 0;	/* a head caps no column, so it carries no forward range */ \
	debug("create head cap(%p), n_forefronts(%lu)", _head, (uint64_t)(_n_forefronts)); \
	dz_cap(_head); \
})
/*
 * Reserve space for the forefront pointers, head cap, and column vector, in
 * that order.
 *
 * Return the slice array address, for which the _spos to _epos range
 * corresponds to the column vector.
 *
 * The returned pointer must be passed to _end_column() before _begin_column()
 * can be called. The final column's slice array address must be passed to
 * _end_matrix() before the arena can be used again for something else.
 */
#define _begin_column_head(_spos, _epos, _adj, _forefronts, _n_forefronts) ({ \
    debug("Beginning column head"); \
    /* calculate sizes */ \
	size_t next_req = _calc_next_size(_spos, _epos, _n_forefronts); \
	/* reserve stream of memory to push onto */  \
	dz_mem_stream_reserve(dz_mem(self), next_req); \
	/* push head-cap */ \
	struct dz_cap_s *cap = _init_cap(_adj, 0xff, _forefronts, _n_forefronts); \
	/* start array slice allocation */ \
	struct dz_swgv_s *slice_data = dz_swgv(dz_mem_stream_alloc_begin(dz_mem(self), _calc_max_slice_size(_spos, _epos))); \
	/* return array pointer */ \
	slice_data - (_spos); \
})

/*
 * Fill in a terminating cap for the most recent column passed to _end_column,
 * using the range information in the current active allocation.
 * Either immediately after it or in a new contiguous run of memory, reserve
 * space for a new column. If the new column is not immediately after the old
 * column, copy a (fake?) version of the old column's cap to be before it.
 *
 * Return the slice array address, for which the _spos to _epos range
 * corresponds to the column vector.
 *
 * The returned pointer must be passed to _end_column() before another column
 * can be begun, and the last column must be passed to _end_matrix() before the
 * arena can be used again.
 */
#define _begin_column(_w, _rch, _rlen) ({ \
    debug("Beginning column"); \
	/* push cap info */ \
	struct dz_cap_s *cap = dz_cap(dz_mem_stream_alloc_current(dz_mem(self))); \
	dz_mem_stream_alloc_end(dz_mem(self), sizeof(struct dz_cap_s)); \
	cap->rch = (_rch);							/* record rch for the next column */ \
	cap->rrem = (_rlen);						/* record rlen for use in traceback */ \
	/* calculate size of next column if it is here */ \
	size_t next_req_here = _calc_next_size((_w).fr.spos, (_w).fr.epos, 0); \
	if(dz_likely(dz_mem_stream_remaining(dz_mem(self)) < next_req_here)) { \
		/* Next column will not fit here. */ \
		debug("create internal bridge"); \
		size_t next_req_split = _calc_next_size((_w).fr.spos, (_w).fr.epos, 1); \
        /* We start a new head column with the previous column's cap pointed to as if it were a forefront. */ \
		dz_mem_stream_reserve(dz_mem(self), next_req_split); \
		cap = _init_cap(0, _rch, &cap, 1); \
	} \
	debug("create column(%p), [%u, %u), span(%u), rrem(%ld), max(%d), inc(%d)", cap, (_w).fr.spos, (_w).fr.epos, (_w).r.epos - (_w).r.spos, (_rlen), (_w).max, (_w).inc); \
	/* start array slice allocation */ \
	struct dz_swgv_s *slice_data = dz_swgv(dz_mem_stream_alloc_begin(dz_mem(self), _calc_max_slice_size((_w).fr.spos, (_w).fr.epos))); \
	/* return array pointer */ \
	slice_data - (_w).fr.spos; \
})
/*
 * Given the slice array address for a column begun with _begin_column() or
 * _begin_column_head(), finish the column.
 *
 * Note that we have actually used the _spos to _epos part of the column's
 * vector, and save range data after that vector recording the range its slice
 * covers.
 *
 * Returns the address of the filled-in range, which it stores in a current
 * allocation, to be made into a cap by _begin_column() or into a forefront by
 * _end_matrix().
 */
 // adding forward range arg
#define _end_column(_p, _spos, _epos, _fr) ({ \
    debug("Ending column"); \
	/* finish the slice data allocation with what was actually used */ \
	dz_mem_stream_alloc_end(dz_mem(self), ((_epos) - (_spos)) * sizeof(struct dz_swgv_s)); \
	/* immediately next in memory, allocate up to a forefront, as a range */ \
	struct dz_range_s *r = dz_range(dz_mem_stream_alloc_begin(dz_mem(self), sizeof(struct dz_forefront_s))); \
	debug("create range(%p), [%u, %u), forward [%u, %u)", r, (_spos), (_epos), (_fr).spos, (_fr).epos); \
	r->spos = (_spos); r->epos = (_epos); \
	/* the traceback needs the forward range as well as the retained one */ \
	dz_cap(r)->fr = (_fr); \
	/* Return it as a cap */ \
	(struct dz_cap_s *)r; \
})
/* 
 * Given the slice array address of the final column, which has been ended (so
 * its filled-in range is in an active forefront-sized allocation), turn that
 * filled-in range into a full forefront.
 *
 * After calling this, no allocation will be active and no memory adjacency
 * relationships will need to be maintained.
 */
#define _end_matrix(_p, _wp, _rrem) ({ \
    debug("Ending matrix"); \
	/* create forefront object */ \
	struct dz_forefront_s *forefront = dz_forefront(dz_mem_stream_alloc_current(dz_mem(self))); \
	dz_mem_stream_alloc_end(dz_mem(self), sizeof(struct dz_forefront_s)); \
	forefront->rid = (_wp)->rid; \
	forefront->rlen = (_wp)->rlen; \
    forefront->fr = (_wp)->fr; \
	forefront->rsum = (_wp)->rsum + ((_wp)->rlen < 0 ? -1 : 1) * ((_wp)->rlen - (_rrem)); \
	forefront->rcnt = (_wp)->rcnt + 1; \
	forefront->max = (_wp)->max + (_wp)->inc; \
	forefront->inc = (_wp)->inc; \
	forefront->query = (_wp)->query; \
	forefront->mcap = (_wp)->mcap; \
	debug("create forefront(%p), [%u, %u), [%u, %u), max(%d), inc(%d), rlen(%d), query(%p), rrem(%d)", forefront, (_wp)->r.spos, (_wp)->r.epos, (_wp)->fr.spos, (_wp)->fr.epos, (_wp)->max, (_wp)->inc, (int32_t)(_wp)->rlen, (_wp)->query, (int32_t)(_rrem)); \
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
#ifdef DZ_QUAL_ADJ
#define _init_rch(_query, _rt, _rrem) \
    _init_bonus(_query); \
    uint32_t rch = conv[_rt[-_rrem] & 0x0f]; \
    /* debug("rch(%c, %u, %x)", _rt[-_rrem], rch, rch); */ \
    uint8_t const *parr = (_query)->arr; \
    uint8_t const *qarr = dz_quals(_query); \
    __m128i const rv = _mm_set1_epi8(rch);

#define _calc_score_profile(_i) ({ \
    /* load qual and query index values */ \
    __m128i qxv = _mm_cvtepi8_epi16(_mm_loadl_epi64((__m128i const *)&qarr[(_i) * L])); \
    __m128i sxv = _mm_or_si128(rv, _mm_loadl_epi64((__m128i const *)&parr[(_i) * L])); /* note: left in 8-bit */\
    __m128i sc; \
    if (_mm_movemask_epi8(_mm_cmpeq_epi8(qxv, _mm_alignr_epi8(qxv, qxv, 2))) == 0xffff) { \
        /* we have uniform quality values, so we can pull from just one matrix */ \
        int8_t const *qual_mat = dz_qual_matrix(self) + 32 * (*((int16_t*)&qxv)); \
        sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)qual_mat), sxv)); \
    }\
    else { \
        /* combine the qual and sequence indexes */\
        __m128i xv = _mm_or_si128(_mm_slli_epi16(qxv, 5), \
                                  _mm_and_si128(_mm_cvtepi8_epi16(sxv), _mm_set1_epi16(0x1f)));\
        print_vector(xv);\
        /* get the scores from the quality matrices (unvectorized, unfortunately) */ \
        int8_t const *qual_mat = dz_qual_matrix(self); \
        sc = _mm_setr_epi16(qual_mat[_mm_extract_epi16(xv, 0)], \
                            qual_mat[_mm_extract_epi16(xv, 1)], \
                            qual_mat[_mm_extract_epi16(xv, 2)], \
                            qual_mat[_mm_extract_epi16(xv, 3)], \
                            qual_mat[_mm_extract_epi16(xv, 4)], \
                            qual_mat[_mm_extract_epi16(xv, 5)], \
                            qual_mat[_mm_extract_epi16(xv, 6)], \
                            qual_mat[_mm_extract_epi16(xv, 7)]);\
    } \
    print_vector(sc)\
    sc; \
})
#else
#define _init_rch(_query, _rt, _rrem) \
    _init_bonus(_query); \
    uint32_t rch = conv[_rt[-_rrem] & 0x0f]; \
    /* debug("rch(%c, %u, %x)", _rt[-_rrem], rch, rch); */ \
    uint8_t const *parr = (_query)->arr; \
    __m128i const rv = _mm_set1_epi8(rch);
                     
#define _calc_score_profile(_i) ({ \
	__m128i qv = _mm_loadl_epi64((__m128i const *)&parr[(_i) * L]); \
	__m128i sc = _mm_cvtepi8_epi16(_mm_shuffle_epi8(_mm_load_si128((__m128i const *)self->matrix), _mm_or_si128(rv, qv))); \
    print_vector(sc); \
     /* print_vector(_mm_cvtepi8_epi16(rv)); print_vector(_mm_cvtepi8_epi16(qv)); */ \
	sc; \
})
#endif
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
	__m128i sc = _add_bonus(_p, _calc_score_profile(_p)); \
	__m128i te = _mm_subs_epi16(_mm_max_epi16(e, _mm_subs_epi16(s, giv)), gev1); \
	/* print_vector(_mm_alignr_epi8(s, ps, 14)); print_vector(sc); */ \
	__m128i ts = _mm_max_epi16(te, _mm_adds_epi16(sc, _mm_alignr_epi8(s, ps, 14))); \
    ps = s; \
    __m128i tf = _mm_subs_epi16(_mm_max_epi16(_mm_subs_epi16(_mm_alignr_epi8(ts, pts, 14), giv), _mm_alignr_epi8(minv, f, 14)), gev1); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 14), gev1)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 12), gev2)); \
	tf = _mm_max_epi16(tf, _mm_subs_epi16(_mm_alignr_epi8(tf, minv, 8), gev4)); \
	ts = _mm_max_epi16(ts, tf); print_vector(ts); \
	maxv = _mm_max_epi16(maxv, ts); \
	/* print_vector(te); print_vector(_add_bonus(_p, ts)); print_vector(tf); print_vector(maxv);*/ \
	e = te; f = tf; s = ts; \
    pts = ts; \
}
#define _store_vector(_p) { \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->e), e); \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->f), f); \
	_mm_store_si128((__m128i *)(&dz_swgv(_p)->s), s); \
}
#define _hmax_vector(_v) ({ \
	__m128i _t = _mm_max_epi16(_v, _mm_srli_si128(_v, 8)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 4)); \
	_t = _mm_max_epi16(_t, _mm_srli_si128(_t, 2)); \
	((int16_t)(_mm_extract_epi16(_t, 0))); \
})
#define _test_xdrop(_s, _xtv) ({ \
	__m128i xtest = _mm_cmpgt_epi16(_s, _xtv); \
	/* print_vector(_s); print_vector(_xtv); */ \
	dz_is_all_zero(xtest); \
})

/**
 * @fn dz_init, dz_destroy
 */
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_s *dz_qual_adj_init_intl(
#else
struct dz_s *dz_init_intl(
#endif
	int8_t const *score_matrix,		/* match award in positive, mismatch penalty in negative. s(A,A) at [0], s(A,C) at [1], ... s(T,T) at [15] where s(ref_base, query_base) is a score function */
#ifdef DZ_QUAL_ADJ
    int8_t const *qual_adj_score_matrix, /* series of DZ_NUM_QUAL_SCORES adjusted score matrices for each phred base qual (no offset) */
#endif
	uint16_t gap_open,				/* gap penalties in positive */
	uint16_t gap_extend)
{
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

    // We need to work out the contiguous block size we need for the dz_s, its matrix data, and the terminal cap.
    // Since they need to be adjacent (TODO: do they?) we use the streaming allocator.
    size_t self_and_matrix_size = sizeof(struct dz_s);
    #if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
        #ifdef DZ_QUAL_ADJ
            self_and_matrix_size += 2 * DZ_QUAL_MATRIX_SIZE;
        #endif
    #else
            self_and_matrix_size += DZ_MAT_SIZE * DZ_MAT_SIZE + 2 * sizeof(__m128i);
    #endif
    size_t cap_size = _calc_next_size(0, 0, 0);
    dz_mem_stream_reserve(mem, self_and_matrix_size + cap_size);

    // TODO: what is the point of the latter 16 bytes all being zero?
	#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
        /* allocate with or without space afterwards for the quality adjusted matrices */
        struct dz_s *self = (struct dz_s *)dz_mem_stream_alloc(mem, self_and_matrix_size);

		/* constants */
		__m128i const tmat = _mm_loadu_si128((__m128i const *)score_matrix);
		_mm_store_si128((__m128i *)&self->matrix[0], tmat);
		_mm_store_si128((__m128i *)&self->matrix[16], _mm_setzero_si128());
    
        #ifdef DZ_QUAL_ADJ
            /* load the quality adjusted matrices immediately after the dz_s */
            __m128i *qual_adj_mats = (__m128i *) (self + 1);
            for(uint64_t i = 0; i < DZ_NUM_QUAL_SCORES; ++i) {
                _mm_store_si128(qual_adj_mats + (i * 2), _mm_loadu_si128(((__m128i const *)qual_adj_score_matrix) + i));
                _mm_store_si128(qual_adj_mats + (i * 2 + 1), _mm_setzero_si128());
            }
        #endif
	#else
		struct dz_s *self = (struct dz_s *)dz_mem_stream_alloc(mem, self_and_matrix_size);

		/* clear the first matrix field for protein */
		_mm_store_si128((__m128i *)&self->matrix[0], _mm_setzero_si128());
		_mm_store_si128((__m128i *)&self->matrix[16], _mm_setzero_si128());

		/* transpose */
		for(uint64_t i = 0; i < DZ_MAT_SIZE; i++) {
			for(uint64_t j = 0; j < DZ_MAT_SIZE; j++) {
				self->protein_matrix[j * DZ_MAT_SIZE + i] = score_matrix[i * DZ_MAT_SIZE + j];
			}
		}
	#endif

	__m128i const giv = _mm_set1_epi16(gap_open);
	__m128i const gev = _mm_set1_epi16(gap_extend);
	_mm_store_si128((__m128i *)self->giv, giv);
	_mm_store_si128((__m128i *)self->gev, gev);
	debug("gi(%u), ge(%u)", gap_open, gap_extend);

    /* initialize and memoize the vectors we need to compute the root */
	__m128i riv = _mm_setr_epi16(
        0,
        -(gap_open+gap_extend),
        -(gap_open+2*gap_extend),
        -(gap_open+3*gap_extend),
        -(gap_open+4*gap_extend),
        -(gap_open+5*gap_extend),
        -(gap_open+6*gap_extend),
        -(gap_open+7*gap_extend)
    );
    _mm_store_si128((__m128i *)self->riv, riv);
    __m128i rev = _mm_setr_epi16(
        -gap_open,
        -(gap_open+gap_extend),
        -(gap_open+2*gap_extend),
        -(gap_open+3*gap_extend),
        -(gap_open+4*gap_extend),
        -(gap_open+5*gap_extend),
        -(gap_open+6*gap_extend),
        -(gap_open+7*gap_extend)
    );
    _mm_store_si128((__m128i *)self->rev, rev);
    
    // precompute the end-cap, which will always live next to the dz_s in this mem block
    struct dz_cap_s *cap = _init_cap(0, 0xff, NULL, 0);

    // Save the stack top so we can tell if we flushed correctly.
    // This isn't quite sufficient to just reconstruct by restoring stored state; we also need curr to go to the right block.
    // TODO: Just store curr and recreate top/end from it and the known allocation size?
    dz_mem(self)->stack.dz_flush_top = dz_mem(self)->stack.top;
    
	return(self);
}

#ifdef DZ_FULL_LENGTH_BONUS
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_s *dz_qual_adj_init(
#else
struct dz_s *dz_init(
#endif
	int8_t const *score_matrix,
#ifdef DZ_QUAL_ADJ
    int8_t const *qual_adj_score_matrix, /* series of DZ_NUM_QUAL_SCORES adjusted score matrices for each phred base qual (no offset) */
#endif
	uint16_t gap_open,
	uint16_t gap_extend)
{
#ifdef DZ_QUAL_ADJ
    return(dz_qual_adj_init_intl(score_matrix, qual_adj_score_matrix, gap_open, gap_extend));
#else
	return(dz_init_intl(score_matrix, gap_open, gap_extend));
#endif
}

#else // DZ_FULL_LENGTH_BONUS
  
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_s *dz_qual_adj_init(
#else
struct dz_s *dz_init(
#endif
	int8_t const *score_matrix,
#ifdef DZ_QUAL_ADJ
    int8_t const *qual_adj_score_matrix, /* series of DZ_NUM_QUAL_SCORES adjusted score matrices for each phred base qual (no offset) */
#endif
	uint16_t gap_open,
	uint16_t gap_extend)//,
	//uint64_t max_gap_len)
{
#ifdef DZ_QUAL_ADJ
    return(dz_qual_adj_init_intl(score_matrix, qual_adj_score_matrix, gap_open, gap_extend, 0));
#else
    return(dz_init_intl(score_matrix, gap_open, gap_extend, 0));
#endif
}
#endif // DZ_FULL_LENGTH_BONUS

#ifndef DZ_INCLUDE_ONCE
                              
// TODO: qual: i think this doesn't change between qual/non-qual

static __dz_vectorize
//#ifdef DZ_QUAL_ADJ
//struct dz_alignment_init_s dz_qual_adj_align_init(
//#else
struct dz_alignment_init_s dz_align_init(
    //#endif
    struct dz_s *self,
    uint32_t max_gap_len)
{
    /* fill the root (the leftmost) column; first init vectors */
    
    /* calc vector length */
    size_t const L = sizeof(__m128i) / sizeof(uint16_t);
    size_t max_gap_vec_len = dz_roundup(max_gap_len, L) / L;
    
    struct dz_forefront_s w = { { 0, (uint32_t) max_gap_vec_len}, 0, 0, { 0, (uint32_t) max_gap_vec_len}, 0, 0, 0, 0, DZ_FOREFRONT_PAD, NULL, NULL };
    struct dz_forefront_s *a = &w;
    
    /* malloc the first column */
    struct dz_swgv_s *dp = _begin_column_head(0, max_gap_vec_len, 0, &a, 0);
    
    int16_t xt = self->giv[0] + self->gev[0] * max_gap_len;
    // calculate the x-drop threshold for this gap length
    __m128i xtv = _mm_set1_epi16(-xt);
        
    /* e, f, and s needed for _store_vector */
    __m128i const e = _mm_set1_epi16(DZ_CELL_MIN);
    /* until the X-drop test fails on all the cells in a vector */
    __m128i s = _mm_load_si128((__m128i const *)self->riv);
    for(size_t p = 0; p < max_gap_vec_len; p++) {
        __m128i const f = s;
        if (_test_xdrop(s, xtv)) {
            debug("p(%lu)", p);
            w.r.epos = p;
            break;
        }
        _store_vector(&dp[p]);
        if (p == 0) {
            s = _mm_load_si128((__m128i const *)self->rev);
        }
        s = _mm_subs_epi16(s, _mm_slli_epi16(_mm_load_si128((__m128i const *)self->gev), 3));
    }
    w.fr.epos = w.r.epos;
    
    /* done; create forefront object */
	// additional fr arg
    _end_column(dp, w.r.spos, w.r.epos, w.fr);
    
    /* package forefront and xt, return */
    struct dz_alignment_init_s aln_init;
    aln_init.root = _end_matrix(dp, &w, 0);
    aln_init.xt = xt;
    return(aln_init);
}
                              
#endif // DZ_INCLUDE_ONCE
                              
#ifndef DZ_INCLUDE_ONCE
/*
 * Destroy a dz_s object when you are done with it.
 */
static __dz_vectorize
void dz_destroy(
	struct dz_s *self)
{
	debug("self(%p)", self);
	if(self == NULL) { return; }
	dz_mem_destroy(dz_mem(self));
	return;
}

/*
 * Limit the bytes of memory retained by a dz_s object to under the given limit.
 * Limit must be greater than the memory required for the initialized, empty state.
 * Can only be called after dz_flush/dz_qual_adj_flush.
 */
static __dz_vectorize
void dz_trim(
	struct dz_s *self,
    size_t max_bytes)
{
    struct dz_mem_block_s *keep_block = &(dz_mem(self)->blk);
    size_t bytes_found = 0;
    while (keep_block != NULL) {
        bytes_found += keep_block->size;
        if (keep_block->next != NULL && bytes_found + keep_block->next->size > max_bytes) {
            /* Cut the free list past here. */
            struct dz_mem_block_s *drop_block = keep_block->next;
            while (drop_block != NULL) {
                struct dz_mem_block_s *next = drop_block->next;
                dz_free(drop_block);
                drop_block = next;
            }
            keep_block->next = NULL;
        }
        keep_block = keep_block->next;
    }

    return;
}
#endif // DZ_INCLUDE_ONCE
                              
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
void dz_qual_adj_flush(
#else
void dz_flush(
#endif
	struct dz_s *self)
{
    // point the mem back at the initial block
	dz_mem_flush(dz_mem(self));
    
    // Re-make the initial allocation, which will re-use any additional memory block it needed.
    // TODO: Share code with dz_init_intl/dz_qual_adj_init_intl which must match this to be correct.
    size_t self_and_matrix_size = sizeof(struct dz_s);
    #if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)
        #ifdef DZ_QUAL_ADJ
            self_and_matrix_size += 2 * DZ_QUAL_MATRIX_SIZE;
        #endif
    #else
            self_and_matrix_size += DZ_MAT_SIZE * DZ_MAT_SIZE + 2 * sizeof(__m128i);
    #endif
    size_t cap_size = _calc_next_size(0, 0, 0);
    dz_mem_stream_reserve(dz_mem(self), self_and_matrix_size + cap_size);


    // Re-allocate self and the matrices
    struct dz_s *new_self = (struct dz_s *)dz_mem_stream_alloc(dz_mem(self), self_and_matrix_size);
    if (new_self != self) {
        dz_error("Lost ourself when resetting memory arena! We should be at %p but actually re-allocated at %p!", self, new_self);
    }
    // Re-allocate the head cap (with no forefronts and no writing).
    // TODO: Must match the behavior of _init_cap
    if (dz_mem_stream_alloc(dz_mem(self), sizeof(struct dz_head_s)) == NULL) {
        dz_error("Could not re-allocate cap after self and matrix!");
    }

    // Now all the stack fields should be right, in particular stack.curr and
    // stack.end, which are now in the block we had to go to to fit all this,
    // if DZ_MEM_INIT_SIZE was too small.
    
    if(dz_mem(self)->stack.top != dz_mem(self)->stack.dz_flush_top) {
        // We didn't get the right stack top at the end of all that.
        dz_error("Could not recreate post-init state when flushing! Stack top should be %p but is actually %p!", dz_mem(self)->stack.dz_flush_top, dz_mem(self)->stack.top);
    }

	return;
}
                  
#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)
                       
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

#endif // DZ_PROTEIN
                       
#endif // !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)

#ifndef DZ_UNITTEST_SCORE_PARAMS
#ifdef DZ_FULL_LENGTH_BONUS
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1
#else
#define DZ_UNITTEST_SCORE_PARAMS	dz_unittest_score_matrix, 5, 1
#endif
#endif
                       
#ifndef DZ_UNITTEST_MAX_GAP_LEN
#define DZ_UNITTEST_MAX_GAP_LEN 20
#endif

#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)
                       
unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);
	//ut_assert(dz->root != NULL);
	dz_destroy(dz);
}
                  
                  
#endif

/**
 * @fn dz_pack_query, dz_pack_query_reverse
 * packed sequence object is allocated inside the context so freed by dz_flush
 */
#if defined(DZ_NUCL_ASCII) || defined(DZ_NUCL_2BIT)

static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_query_s *dz_qual_adj_pack_query_forward(
#else
struct dz_query_s *dz_pack_query_forward(
#endif
	struct dz_s *self,
    char const *query,
#ifdef DZ_QUAL_ADJ
    uint8_t const *qual,
#endif
#ifdef DZ_FULL_LENGTH_BONUS
    int8_t full_length_bonus,
#endif
    size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
    // TODO: qual: make sure I'm doing this right
    size_t pack_size = dz_roundup(qlen + 1, sizeof(__m128i));
    #ifdef DZ_QUAL_ADJ
        struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + 2 * pack_size);
    #else
        struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + pack_size);
    #endif
	*q = (struct dz_query_s){
		.blen = (uint32_t) ((qlen == 0 ? 0 : dz_roundup(qlen + 1, L) / L)),
        .arr_end = (uint32_t) pack_size,
		.q = query,
		.bonus = { 0 }
	};

    /*
	 * tentative support of full-length bonus; precaclulated bonus score is saved at q->bonus[0] and q->bonus[1]; [0] for non-tail vector and [1] for tail vector
	 */
    #ifdef DZ_FULL_LENGTH_BONUS
        q->bonus[L + (qlen % L)] = full_length_bonus;
    #endif
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
	for(size_t i = 0, end = dz_rounddown(qlen, sizeof(__m128i)); i < end; i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[i]);
		__m128i tv = _mm_shuffle_epi8(cv, _mm_and_si128(qv, fv));
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); /* shift by one to make room for a base */
        pv = tv;
	}

	/* continue the same conversion on the remainings */
	// _mm_store_si128((__m128i *)&q->arr[dz_rounddown(qlen, sizeof(__m128i))], _mm_srli_si128(pv, 15));
	q->arr[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
	for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[(uint8_t)query[i] & 0x0f];
	}
    // TODO: does he need the + 1 here? maybe that's why he is buffering the extra sizeof(__m128i) in the alloc?
	for(size_t i = qlen, end = dz_roundup(qlen + 1, sizeof(__m128i)); i < end; i++) {
		q->arr[i + 1] = qS;
	}
    
    // TODO: qual: probably could make this faster if we apply the <<5 once here and store in 16 bit values
    
    #ifdef DZ_QUAL_ADJ
        __m128i const qmax = _mm_set1_epi8(DZ_NUM_QUAL_SCORES - 1);
        __m128i qpv = _mm_setzero_si128();
    
        /* until the end of the query sequence */
        for(size_t i = 0, end = dz_rounddown(qlen, sizeof(__m128i)); i < end; i += sizeof(__m128i)) {
            __m128i const qtv = _mm_min_epu8(_mm_loadu_si128((__m128i const *)&qual[i]), qmax);
            _mm_store_si128((__m128i *)&q->arr[pack_size + i], _mm_alignr_epi8(qtv, qpv, 15)); /* shift by one to make room for a base */
            qpv = qtv;
        }
        
        /* continue the same conversion on the remainings */
        q->arr[pack_size + dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(qpv, 15);
        for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
            q->arr[pack_size + i + 1] = dz_min2(qual[i], DZ_NUM_QUAL_SCORES - 1);
        }
        for(size_t i = qlen, end = dz_roundup(qlen + 1, sizeof(__m128i)); i < end; i++) {
            q->arr[pack_size + i + 1] = 0;
        }
    #endif

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(q);
}
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_query_s *dz_qual_adj_pack_query_reverse(
#else
struct dz_query_s *dz_pack_query_reverse(
#endif
    struct dz_s *self,
    char const *query,
#ifdef DZ_QUAL_ADJ
    uint8_t const *qual,
#endif
#ifdef DZ_FULL_LENGTH_BONUS
    int8_t full_length_bonus,
#endif
    size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
    // TODO: qual: make sure I'm doing this right
    size_t pack_size = dz_roundup(qlen + 1, sizeof(__m128i));
    #ifdef DZ_QUAL_ADJ
        struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + 2 * pack_size);
    #else
        struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + pack_size);
    #endif
	*q = (struct dz_query_s){
		.blen = (uint32_t) (qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L)),
        .arr_end = (uint32_t) pack_size,
		.q = query,
		.bonus = { 0 }
	};
    #ifdef DZ_FULL_LENGTH_BONUS
        q->bonus[L + (qlen % L)] = full_length_bonus;
    #endif

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
	__m128i const cv = _mm_load_si128((__m128i const *)conv), rv = _mm_load_si128((__m128i const *)rev);

	/* until the end of the query sequence */
	for(size_t i = 0, end = dz_rounddown(qlen, sizeof(__m128i)); i < end; i += sizeof(__m128i)) {
		__m128i const qv = _mm_loadu_si128((__m128i const *)&query[qlen - 16 - i]);
		__m128i tv = _mm_shuffle_epi8(_mm_shuffle_epi8(cv, _mm_and_si128(qv, fv)), rv);
		_mm_store_si128((__m128i *)&q->arr[i], _mm_alignr_epi8(tv, pv, 15)); pv = tv;			/* shift by one to make room for a base */
	}

	/* continue the same conversion on the remainings */
	q->arr[dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(pv, 15);
	for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
		q->arr[i + 1] = conv[(uint8_t)query[qlen - 1 - i] & 0x0f];
	}
	for(size_t i = qlen, end = dz_roundup(qlen + 1, sizeof(__m128i)); i < end; i++) {
		q->arr[i + 1] = qS;
	}
    
    // TODO: qual: finish this rev
    #ifdef DZ_QUAL_ADJ
        
        __m128i const qmax = _mm_set1_epi8(DZ_NUM_QUAL_SCORES - 1);
        __m128i qpv = _mm_setzero_si128();
        
        /* until the end of the query sequence */
        for(size_t i = 0, end = dz_rounddown(qlen, sizeof(__m128i)); i < end; i += sizeof(__m128i)) {
            __m128i qtv = _mm_min_epu8(_mm_shuffle_epi8(_mm_loadu_si128((__m128i const *)&qual[qlen - 16 - i]), rv), qmax);
            _mm_store_si128((__m128i *)&q->arr[pack_size + i], _mm_alignr_epi8(qtv, qpv, 15));
            qpv = qtv;            /* shift by one to make room for a base */
        }
        
        /* continue the same conversion on the remainings */
        q->arr[pack_size + dz_rounddown(qlen, sizeof(__m128i))] = _mm_extract_epi8(qpv, 15);
        for(size_t i = dz_rounddown(qlen, sizeof(__m128i)); i < qlen; i++) {
            q->arr[pack_size + i + 1] = dz_min2(qual[qlen - 1 - i], DZ_NUM_QUAL_SCORES - 1);
        }
        for(size_t i = qlen, end = dz_roundup(qlen + 1, sizeof(__m128i)); i < end; i++) {
            q->arr[pack_size + i + 1] = 0;
        }
    #endif

	debug("qlen(%lu), q(%.*s)", qlen, (int)qlen, query);
	return(q);
}

#else /* DZ_PROTEIN */
static __dz_vectorize
struct dz_query_s *dz_pack_query_forward(
	struct dz_s *self,
	char const *query,
#ifdef DZ_FULL_LENGTH_BONUS
    int8_t full_length_bonus,
#endif
    size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + DZ_MAT_SIZE * dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
    
    #ifdef DZ_FULL_LENGTH_BONUS
        q->bonus[L + (qlen % L)] = full_length_bonus;
    #endif

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
			a[i + 1] = -1;
		}
	}
	debug("qlen(%lu), q(%s)", qlen, query);
	return(q);
}
static __dz_vectorize
struct dz_query_s *dz_pack_query_reverse(
	struct dz_s *self,
	char const *query,
#ifdef DZ_FULL_LENGTH_BONUS
    int8_t full_length_bonus,
#endif
    size_t qlen)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	struct dz_query_s *q = (struct dz_query_s *)dz_mem_malloc(dz_mem(self), sizeof(struct dz_query_s) + DZ_MAT_SIZE * dz_roundup(qlen + 1, L) + sizeof(__m128i));
	*q = (struct dz_query_s){
		.blen = qlen == 0 ? 0 : (dz_roundup(qlen + 1, L) / L),
		.q = query,
		.bonus = { 0 }
	};
    #ifdef DZ_FULL_LENGTH_BONUS
        q->bonus[L + (qlen % L)] = full_length_bonus;
    #endif

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
			a[i + 1] = -1;
		}
	}

	debug("qlen(%lu), q(%s)", qlen, query);
	return(q);
}
#endif
                                                  
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_query_s *dz_qual_adj_pack_query(
#else
struct dz_query_s *dz_pack_query(
#endif
	struct dz_s *self,
	char const *query,
#ifdef DZ_QUAL_ADJ
    const uint8_t *qual,
#endif
#ifdef DZ_FULL_LENGTH_BONUS
    int8_t full_length_bonus,
#endif
    int64_t qlen)
{
	if(qlen >= 0) {
#ifdef DZ_QUAL_ADJ
        return(dz_qual_adj_pack_query_forward(self, query, qual, full_length_bonus, (size_t)qlen));
#else
        return(dz_pack_query_forward(self, query, full_length_bonus, (size_t)qlen));
#endif
	} else {
#ifdef DZ_QUAL_ADJ
        return(dz_qual_adj_pack_query_reverse(self, &query[qlen], &qual[qlen], full_length_bonus, (size_t)-qlen));
#else
        return(dz_pack_query_reverse(self, &query[qlen], full_length_bonus, (size_t)-qlen));
#endif
	}
}

#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)

unittest() {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query,
#ifdef DZ_FULL_LENGTH_BONUS
                                         0,
#endif
                                         
                                         dz_unittest_query_length);
	dz_unused(q);
	dz_destroy(dz);
}
                                          
#endif

#define _merge_column(w, adj, forefronts, n_forefronts, query, init_s) ({ \
	for(size_t i = 0; i < n_forefronts; i++) { \
		/* update max and pos */ \
		w.fr.spos = dz_min2(w.fr.spos, forefronts[i]->fr.spos); \
		w.fr.epos = dz_max2(w.fr.epos, forefronts[i]->fr.epos); \
		w.rsum = dz_max2(w.rsum, forefronts[i]->rsum); \
		w.rcnt = dz_max2(w.rcnt, forefronts[i]->rcnt); \
		if(init_s != 0) { w.max = dz_max2(w.max, forefronts[i]->max); } \
		debug("i(%lu, %p), r [%u, %u), fr [%u, %u),  inc(%d), max(%d)", i, forefronts[i], forefronts[i]->r.spos, forefronts[i]->r.epos, forefronts[i]->fr.spos, forefronts[i]->fr.epos, forefronts[i]->inc, forefronts[i]->max); \
	} \
	w.fr.spos = w.fr.spos > INT32_MAX ? 0 : w.fr.spos; \
	w.fr.epos = dz_min2(w.fr.epos, query->blen);						/* make sure epos and p-iteration index is no larger than blen */ \
	debug("start extension [%u, %u), inc(%d), max(%d), stack(%p, %p), rlen(%d)", w.fr.spos, w.fr.epos, w.inc, w.max, dz_mem(self)->stack.top, dz_mem(self)->stack.end, (int32_t)rlen); \
	if(init_s != 0) { \
		for(size_t i = 0; i < n_forefronts; i++) { \
			adj[i] = w.max - (forefronts[i]->max - forefronts[i]->inc);	/* base = max - inc */ \
			debug("i(%lu), adj(%lu)", i, adj[i]); \
		} \
	} \
	/* once; merge incoming vectors */ \
	struct dz_swgv_s *cdp = _begin_column_head(w.fr.spos, w.fr.epos, w.max, forefronts, n_forefronts);	/* allocate memory for the first column */ \
	for(uint64_t p = w.fr.spos; p < w.fr.epos; p++) { \
		/* memset(cdp, 0xff, sizeof(struct dz_swgv_s) * (w.r.epos - w.r.spos)); */ \
		__m128i const e = _mm_set1_epi16(INT16_MIN), f = e, s = e; _store_vector(&cdp[p]); \
	} \
	/* paste the last vectors */ \
	for(size_t i = 0; i < n_forefronts; i++) { \
        debug("Add column %lu", i); \
		struct dz_swgv_s const *tdp = dz_cswgv(forefronts[i]) - forefronts[i]->r.epos; \
		__m128i const adjv = _mm_set1_epi16(init_s == 0 ? 0 : adj[i]); \
		for(uint64_t p = forefronts[i]->fr.spos; p < forefronts[i]->fr.epos; p++) { \
			/* adjust offset */ \
			__m128i e = _mm_subs_epi16(_mm_load_si128(&tdp[p].e), adjv); \
			__m128i f = _mm_subs_epi16(_mm_load_si128(&tdp[p].f), adjv); \
			__m128i s = _mm_subs_epi16(_mm_load_si128(&tdp[p].s), adjv); \
			/* read-max-write */ \
			_mm_store_si128(&cdp[p].e, _mm_max_epi16(e, _mm_load_si128(&cdp[p].e))); \
			_mm_store_si128(&cdp[p].f, _mm_max_epi16(f, _mm_load_si128(&cdp[p].f))); \
			_mm_store_si128(&cdp[p].s, _mm_max_epi16(s, _mm_load_si128(&cdp[p].s))); \
			print_vector(s); \
		} \
	} \
    w.r = w.fr;\
	_end_column(cdp, w.fr.spos, w.fr.epos, w.fr); \
	cdp; \
})
#define _fill_column(w, pdp, query, rt, rrem, xt, init_s) ({ \
	/* load the base */ \
	_init_rch(query, rt, rrem); \
	struct dz_swgv_s *cdp = _begin_column(w, rch, rrem); \
	/* init vectors */ \
	__m128i f = minv, ps = _mm_set1_epi16(init_s), maxv = _mm_set1_epi16(INT16_MIN), pts = _mm_set1_epi16(INT16_MIN); \
	__m128i const xtv = _mm_set1_epi16(w.inc - xt);	/* next offset == current max thus X-drop threshold is always -xt */ \
	/* until the bottommost vertically placed band... */ \
    uint64_t first_drop = w.fr.epos;\
    /* we want to remember the backward facing range even if we x-drop the forward facing range */ \
    w.r = w.fr; \
	for(uint64_t p = w.fr.spos; p < w.fr.epos; p++) { \
		_load_vector(&pdp[p]); \
        _update_vector(p); \
		if(dz_unlikely(_test_xdrop(s, xtv))) {	/* mark _unlikely to move out of the core loop */ \
			/* drop either leading or trailing vector, skip the forefront extension when the forefront is clipped */ \
            pts = _mm_set1_epi16(INT16_MIN); \
			if(p == w.fr.spos) { \
                w.fr.spos++; \
            } else { \
                first_drop = dz_min2(p, first_drop); \
            } \
        } else { \
            first_drop = w.fr.epos; /* reset the xdrop to indicate that any xdrops we found are not contiguous */ \
        } \
		_store_vector(&cdp[p]); \
	} \
    if (first_drop < w.fr.epos) { \
        /* we xdropped a contiguous block of vectors at the end of this column */\
        w.fr.epos = first_drop; \
        goto dz_pp_cat(_forefront_, __LINE__); \
    } \
    /* if reached the forefront of the query sequence, finish the extension */ \
	if(w.fr.epos < query->blen) { \
		/* forefront extension; clip the column length if too long */ \
		__m128i e = minv, s = minv; _update_vector(w.fr.epos); \
		do { \
			if(_test_xdrop(s, xtv)) { break; } \
            _store_vector(&cdp[w.fr.epos]); w.fr.epos++; w.r.epos++; \
			f = _mm_subs_epi16(f, gev8); s = _mm_subs_epi16(s, gev8); \
		} while(w.fr.epos < query->blen); \
	} \
dz_pp_cat(_forefront_, __LINE__):; \
	/* create cap object that contains [spos, epos) range (for use in the traceback routine), adding forward range since that is now used for traceback*/ \
	struct dz_cap_s *cap = _end_column(cdp, w.r.spos, w.r.epos, w.fr); \
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
#ifdef DZ_QUAL_ADJ
struct dz_forefront_s const *dz_qual_adj_extend_intl(
#else
struct dz_forefront_s const *dz_extend_intl(
#endif
    struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid, uint16_t xt,
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

	__m128i const minv = _mm_set1_epi16(DZ_CELL_MIN);
	__m128i const giv = _mm_load_si128((__m128i const *)self->giv);
	__m128i const gev1 = _mm_load_si128((__m128i const *)self->gev);
	__m128i const gev2 = _mm_add_epi16(gev1, gev1);
	__m128i const gev4 = _mm_slli_epi16(gev1, 2);
	__m128i const gev8 = _mm_slli_epi16(gev1, 3);

    struct dz_forefront_s w = { { UINT32_MAX, 0 }, 0, 0, { UINT32_MAX, 0 }, 0, 0, 0, 0, DZ_FOREFRONT_PAD, NULL, NULL };	/* uint32_t spos, epos, max, inc; uint8_t _pad[8]; struct dz_query_s const *query; struct dz_cap_s const *cap; */
	w.rlen = rlen;
	w.rid = rid;
	w.query = query;
	/* first iterate over the incoming edge objects to get the current max */
	uint64_t adj[n_forefronts];							/* keep variable length array out of statement expression to avoid a bug of icc */
	struct dz_swgv_s *pdp = _merge_column(w, adj, forefronts, n_forefronts, query, init_s);
	/* fetch the first base */
	int64_t rrem = rlen, dir = rlen < 0 ? 1 : -1;
	uint8_t const *rt = (uint8_t const *)&ref[rrem - (rlen < 0)];

	while(rrem != 0 && w.fr.spos < w.fr.epos) {
		pdp = _fill_column(w, pdp, query, rt, rrem, xt, init_s); rrem += dir;
		debug("rrem(%d), [%u, %u)", rrem, w.fr.spos, w.fr.epos);
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
#ifdef DZ_QUAL_ADJ
struct dz_forefront_s const *dz_qual_adj_extend(
#else
struct dz_forefront_s const *dz_extend(
#endif
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid, uint16_t xt)
{
    debug("entering extend (external)");
#ifdef DZ_QUAL_ADJ
    return(dz_qual_adj_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, xt, INT16_MIN));
#else
	return(dz_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, xt, INT16_MIN));
#endif
}
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_forefront_s const *dz_qual_adj_scan(
#else
struct dz_forefront_s const *dz_scan(
#endif
	struct dz_s *self,
	struct dz_query_s const *query,
	struct dz_forefront_s const **forefronts, size_t n_forefronts,
	char const *ref, int32_t rlen, uint32_t rid, uint16_t xt)
{
	debug("dz_scan called");
#ifdef DZ_QUAL_ADJ
    return(dz_qual_adj_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, xt, 0));
#else
	return(dz_extend_intl(self, query, forefronts, n_forefronts, ref, rlen, rid, xt, 0));
#endif
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

#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)
                                              
/* short, exact matching sequences */
unittest( "extend.base" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 3; trial++) {
		struct dz_query_s const *q = dz_pack_query(dz, dz_unittest_query,
#ifdef DZ_FULL_LENGTH_BONUS
                                                   0,
#endif
			  trial == 0 ? dz_unittest_query_length
			: trial == 1 ? 3
			:              10
		);
		struct dz_forefront_s const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1, 0);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2, 0);
		ut_assert(forefront == NULL);
        struct dz_alignment_init_s aln_init = dz_align_init(dz, DZ_UNITTEST_MAX_GAP_LEN);
		forefront = dz_extend(dz, q, &aln_init.root, 1, "", 0, 3, aln_init.xt);
		ut_assert(forefront == aln_init.root);

		/* extend */
		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("A", "\x0", "M"), 1, 4, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 5));

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 5, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "MASLVQT"), 7, 6, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "MASLVQTG"), 8, 7, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 34));

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "MASLVQTGK"), 9, 8, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 39));

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
		struct dz_query_s const *q = dz_pack_query_reverse(dz, &revcomp[dz_unittest_query_length - length],
#ifdef DZ_FULL_LENGTH_BONUS
                                                           0,
#endif
                                                           length);
		struct dz_forefront_s const *forefront = NULL;

		/* nothing occurs */
		forefront = dz_extend(dz, q, NULL, 0, NULL, 0, 1, 0);
		ut_assert(forefront == NULL);
		forefront = dz_extend(dz, q, NULL, 0, "", 0, 2, 0);
		ut_assert(forefront == NULL);
        struct dz_alignment_init_s aln_init = dz_align_init(dz, DZ_UNITTEST_MAX_GAP_LEN);
		forefront = dz_extend(dz, q, &aln_init.root, 1, "", 0, 3, aln_init.xt);
		ut_assert(forefront == aln_init.root);

		/* extend */
		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("A", "\x0", "M"), 1, 4, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(2, 2, 5));

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 5, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(4, 4, 9));
		if(trial == 1) { continue; }

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTT", "\x0\x2\x0\x3\x3\x3\x3", "MASLVQT"), 7, 6, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28));
		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel(&"AAAATCT"[7], &"\x0\x0\x0\x0\x3\x1\x3"[7], "MASLVQT"), -7, 6, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(9, 9, 28));
		if(trial == 2) { continue; }

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTTC", "\x0\x2\x0\x3\x3\x3\x3\x1", "MASLVQTG"), 8, 7, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 34));
		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel(&"GAAAATCT"[8], &"\x2\x0\x0\x0\x0\x3\x1\x3"[8], "MASLVQTG"), -8, 7, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(11, 11, 28));

		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AGATTTTCA", "\x0\x2\x0\x3\x3\x3\x3\x1\x0", "MASLVQTGK"), 9, 8, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 39));
		forefront = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel(&"TGAAAATCT"[9], &"\x3\x2\x0\x0\x0\x0\x3\x1\x3"[9], "MASLVQTGK"), -9, 8, aln_init.xt);
		ut_assert(forefront != NULL && forefront->max == dz_ut_sel(13, 13, 28));

		(void)forefront;
	}
	dz_destroy(dz);
}
                                              
#endif // !DZ_PROTEIN

/* a small graph */
unittest( "extend.small" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	for(size_t trial = 0; trial < 4; trial++) {
		struct dz_query_s const *q = dz_pack_query(dz, dz_unittest_query,
#ifdef DZ_FULL_LENGTH_BONUS
                                                   0,
#endif
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
        struct dz_alignment_init_s aln_init = dz_align_init(dz, DZ_UNITTEST_MAX_GAP_LEN);
		forefronts[0] = dz_extend(dz, q, &aln_init.root, 1, dz_ut_sel("AG", "\x0\x2", "MA"), 2, 1, aln_init.xt);
		ut_assert(forefronts[0] != NULL && forefronts[0]->max == dz_ut_sel(4, 4, 9));

		forefronts[1] = dz_extend(dz, q, &forefronts[0], 1, dz_ut_sel("C", "\x1", "T"), 1, 2, aln_init.xt);
		ut_assert(forefronts[1] != NULL && forefronts[1]->max == dz_ut_sel(6, 6, 14));
		if(trial == 1) { continue; }

		forefronts[2] = dz_extend(dz, q, &forefronts[0], 2, dz_ut_sel("TTTT", "\x3\x3\x3\x3", "LVQT"), 4, 3, aln_init.xt);
		if(trial == 2) {
			ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(8, 8, 18));
			continue;
		}
		ut_assert(forefronts[2] != NULL && forefronts[2]->max == dz_ut_sel(14, 14, 32));
		if(trial == 3) { continue; }

		forefronts[3] = dz_extend(dz, q, &forefronts[2], 1, dz_ut_sel("CATT", "\x1\x0\x3\x3", "CKAK"), 4, 4, aln_init.xt);
		ut_assert(forefronts[3] != NULL && forefronts[3]->max == dz_ut_sel(22, 22, 43));

		forefronts[4] = dz_extend(dz, q, &forefronts[2], 2, dz_ut_sel("CTGA", "\x1\x3\x2\x0", "QLTL"), 4, 5, aln_init.xt);
		ut_assert(forefronts[4] != NULL && forefronts[4]->max == dz_ut_sel(30, 30, 61));
	}
	dz_destroy(dz);
}
                                              
#endif // !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)

// TODO: qual: this can probably be shared with qual/non-qual
       
#ifndef DZ_INCLUDE_ONCE
                                              
/**
 * @fn dz_calc_max_rpos
 */
static __dz_vectorize
int64_t dz_calc_max_rpos(
	struct dz_forefront_s const *forefront)
{
	struct dz_cap_s const *pcap = forefront->mcap;
	int32_t rpos = (pcap == dz_ccap(forefront)) ? 0 : pcap->rrem;
	return((int64_t)rpos);									/* signed expansion */
}

/**
 * @fn dz_calc_max_qpos
 */
static __dz_vectorize
uint64_t dz_calc_max_qpos(struct dz_forefront_s const *forefront)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
	#define _dp(_cap)					( (struct dz_swgv_s const *)(_cap) - (_cap)->r.epos )

	_init_bonus(forefront->query);
	struct dz_cap_s const *pcap = forefront->mcap;
	__m128i const maxv = _mm_set1_epi16(forefront->inc);
	/* changing this loop as well, since only the forward range can hold the max, and the traceback can only walk within it */
	for(uint64_t p = pcap->fr.spos; p < pcap->fr.epos; p++) {
		__m128i const s = _mm_load_si128(&_dp(pcap)[p].s); print_vector(s);
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
uint64_t dz_calc_max_pos(struct dz_forefront_s const *forefront)
{
	struct dz_cap_s const *pcap = forefront->mcap;
	if(pcap == NULL) { debug("pcap == NULL, rlen(%d)", (int32_t)forefront->rlen); return(((uint64_t)forefront->rlen)<<32); }							/* is root */
	debug("forefront(%p), pcap(%p), rrem(%d), rlen(%d)", forefront, pcap, (int32_t)pcap->rrem, (int32_t)forefront->rlen);
	return((dz_calc_max_rpos(forefront)<<32) | dz_calc_max_qpos(forefront));
}
                                              
#endif // DZ_INCLUDE_ONCE
                                              
/* traceback macros */
#define _is_head(_cap)				( (_cap)->r.epos == 0 )
#define _dp(_cap)					( dz_cswgv(_cap) - (_cap)->r.epos )
#define _vector_idx(_idx)			( (_idx) / L )
#define _cell_idx(_idx)				( (_idx) & (L - 1) )
#define _s(_l, _cap, _idx)			( ((int16_t const *)(&_dp(_cap)[_vector_idx(_idx)]._l))[_cell_idx(_idx)] )
#define _unwind_cap(_c)				( dz_ccap(dz_cswgv(_c) - (_c)->r.epos + (_c)->r.spos) - 1 )
#define _push_span(_id)				{ span->offset = path_base - path; span--; span->id = (_id); }
#define _idx_asc(_i)				( -((int64_t)n_forefronts) + (_i) )
#define _idx_dsc(_i)				( -1LL - (_i) )
#define _load_prev_cap(_l, _score, _idx) ({ \
	while(1) { \
		cap = pcap; pcap = _unwind_cap(cap); \
		debug("idx(%lu), cap(%p), pcap(%p), [%u, %u)", idx, cap, pcap, pcap->r.spos, pcap->r.epos); \
		/* escape if not head */ \
		if(dz_likely(!_is_head(pcap))) { \
            break; \
        } \
		/* escape if root */ \
		debug("head, n_forefronts(%u)", dz_chead(pcap)->n_forefronts); \
		if(pcap->rrem == 0) { debug("reached root"); goto _trace_tail; } \
		/* just load pcap if internal bridge */ \
		struct dz_forefront_s const **forefront_arr = (struct dz_forefront_s const **)pcap; \
		if(pcap->rch != 0xff) { debug("internal bridge, pcap(%p)", dz_ccap(forefront_arr[-1])); pcap = dz_ccap(forefront_arr[-1]); break; } \
		/* merging vector; load contents to find an edge */ \
		int16_t prev_score = _s(_l, cap, idx) + pcap->r.spos;		/* adj */ \
		size_t n_forefronts = dz_chead(pcap)->n_forefronts; \
		debug("merging vector, %s(%d), n_forefronts(%lu)", #_l, prev_score, n_forefronts); \
		for(size_t i = 0; i < n_forefronts; i++) { \
			struct dz_forefront_s const *ff = forefront_arr[_idx(i)]; \
			pcap = dz_ccap(ff); \
			debug("i(%lu), %s(%d, %d), max(%d), inc(%d), [%u, %u)", i, #_l, prev_score, _s(_l, pcap, idx) + (ff->max - ff->inc), ff->max, ff->inc, ff->r.spos, ff->r.epos); \
			/* adj[i] = w.max - (ffs[i]->max - ffs[i]->inc); base = max - inc */ \
			/* This is also needed for correctly respecting the band boundaries during traceback. Else out of band columns may be merged incorrectly */ \
			if(!dz_inside(pcap->fr.spos, _vector_idx(idx), pcap->fr.epos)) { debug("out of range(%u), [%u, %u)", _vector_idx(idx), pcap->fr.spos, pcap->fr.epos); continue; } \
			if(prev_score == _s(_l, pcap, idx) + (ff->max - ff->inc)) { debug("found, i(%lu), id(%u), prev_score(%d), score(%d), max(%d), inc(%d)", i, dz_cff(pcap)->rid, prev_score, _s(_l, pcap, idx), ff->max, ff->inc); break; } \
		} \
		_push_span(dz_cff(pcap)->rid);								/* push segment info */ \
		_score = _s(_l, pcap, idx); \
	} \
	/* return the reference-side base */ \
	ref_length++; pcap->rch; \
})

/**
 * @fn dz_trace
 */
static __dz_vectorize
#ifdef DZ_QUAL_ADJ
struct dz_alignment_s *dz_qual_adj_trace(
#else
struct dz_alignment_s *dz_trace(
#endif
	struct dz_s *self,
	struct dz_forefront_s const *forefront)
{
	size_t const L = sizeof(__m128i) / sizeof(uint16_t);
    if(forefront->mcap == NULL) { debug("mcap is null"); return(NULL); }
	/* detect pos */
	uint64_t idx = dz_calc_max_qpos(forefront);							/* vector index, cell index */
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
	int32_t score = _s(s, pcap, idx), cnt[4] = { 0 };
    
    #ifdef DZ_QUAL_ADJ
    #    define _pair_score(_self, _q, _r, _i)   ( dz_qual_adj_pair_score((_self), (_q), (_r), (_i)) )
    #else
    #    define _pair_score(_self, _q, _r, _i)   ( dz_pair_score((_self), (_q), (_r), (_i)) )
    #endif

	/* traceback loop */
	#define _debug(_label) { \
		debug("test %s, idx(%lu), rbase(%d, %c), qbase(%c), c(%d), e(%d), f(%d), s(%d), score(%d, %d, %lu), diag(%d)", \
			#_label, idx, rch, '-', '-', /*"ACGTNNNNNNNNNNNN"[rch & 0xff], "ANNNCNNNGNNNTNNN"[query->arr[idx]],*/ \
			score, _s(e, cap, idx), _s(f, cap, idx), _s(s, pcap, idx - 1), \
			_pair_score(self, query, rch, idx), (uint64_t)dz_pair_eq(self, query, rch, idx), \
			dz_end_bonus(self, query, idx), (score == (_s(s, pcap, idx - 1) + _pair_score(self, query, rch, idx)))); \
	}
	#define _match(_idx) { \
        _debug(M); \
        if(dz_inside(pcap->fr.spos, _vector_idx(idx - 1), pcap->fr.epos) \
           && score == (_s(s, pcap, idx - 1) + _pair_score(self, query, rch, idx))) { \
            uint64_t eq = dz_pair_eq(self, query, rch, idx); \
            *--path = DZ_CIGAR_OP>>(eq<<3); cnt[eq]++; \
            score = _s(s, pcap, idx - 1); idx--; rch = _load_prev_cap(s, score, _idx); \
            continue; \
		} \
	}
	#define _ins(_idx) { \
		if(dz_inside(cap->fr.spos, _vector_idx(idx - 1), cap->fr.epos) && score == _s(f, cap, idx)) { \
			_debug(I); \
			while(_vector_idx(idx - 1) >= cap->fr.spos && score != _s(s, cap, idx - 1) - self->gev[0] - self->giv[0]) { \
				*--path = (DZ_CIGAR_OP>>16) & 0xff; cnt[2]++; score = _s(f, cap, idx - 1); idx--; _debug(I); \
			} \
			*--path = (DZ_CIGAR_OP>>16) & 0xff; cnt[2]++; score = _s(s, cap, idx - 1); idx--; \
			continue; \
		} \
	}
	#define _del(_idx) { \
		if(dz_inside(pcap->fr.spos, _vector_idx(idx), pcap->fr.epos) && score == _s(e, cap, idx)) { \
			_debug(D); \
			while(dz_inside(pcap->fr.spos, _vector_idx(idx), pcap->fr.epos) && score == _s(e, pcap, idx) - self->gev[0]) { \
                *--path = (DZ_CIGAR_OP>>24) & 0xff; cnt[3]++; score = _s(e, pcap, idx); _load_prev_cap(e, score, _idx); _debug(D); \
			} \
            *--path = (DZ_CIGAR_OP>>24) & 0xff; cnt[3]++; score = _s(s, pcap, idx); rch = _load_prev_cap(s, score, _idx); \
			continue; \
		} \
	}

	uint32_t rch = _load_prev_cap(s, score, _idx_asc);
	while(1) {
        _match(_idx_asc);
        _ins(_idx_asc);
        _del(_idx_asc);
        dz_trap();
    }
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
	aln->rrem = dz_calc_max_rpos(forefront);
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
    
    #undef _pair_score
    #undef _debug
    #undef _match
    #undef _ins
    #undef _del
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

#if !defined(DZ_UNITTESTS_INCLUDED) && !defined(DZ_QUAL_ADJ)
                                         
/* short, exact matching sequences */
unittest( "trace" ) {
	struct dz_s *dz = dz_init(DZ_UNITTEST_SCORE_PARAMS);
	ut_assert(dz != NULL);

	struct dz_query_s *q = dz_pack_query(dz, dz_unittest_query,
#ifdef DZ_FULL_LENGTH_BONUS
                                         0,
#endif
                                         dz_unittest_query_length);
	struct dz_forefront_s const *forefront = NULL;
	struct dz_alignment_s *aln = NULL;

    struct dz_alignment_init_s aln_init = dz_align_init(dz, DZ_UNITTEST_MAX_GAP_LEN);
	forefront = dz_extend(dz, q, &aln_init.root, 1, "A", 1, 1, aln_init.xt);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "GC", 2, 2, aln_init.xt);
	aln = dz_trace(dz, forefront);

	forefront = dz_extend(dz, q, &forefront, 1, "TTTT", 4, 3, aln_init.xt);
	aln = dz_trace(dz, forefront);

	(void)aln;
	dz_destroy(dz);
}

#endif

#ifdef __cplusplus
};	/* extern "C" { */
#endif

// make sure all of the non-duplicated functions are only included once if we also include
// include the quality adjusted versions
#ifndef DZ_INCLUDE_ONCE
#define DZ_INCLUDE_ONCE
#endif

// make sure the unittests only get included once and with non-qual-adjusted functions
#if !defined(DZ_QUAL_ADJ) && !defined(DZ_UNITTESTS_INCLUDED)
#define DZ_UNITTESTS_INCLUDED
#endif

#undef SIMDE_ENABLE_NATIVE_ALIASES

/**
 * end of dozeu.h
 */
