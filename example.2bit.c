/**
 * Aligning ACACTTCTAGACTTTACCACTA to the following graph;
 * sequences are encoded as follows: {'A','C','G','T'} -> {0x0,0x1,0x2,0x3}
 *                         6
 *                         A
 *  0    1    3      5    /7\9
 *  ACAC-TTGT-AGAC---TTCTA-C-CACGG
 *      \    /    \ /     \ /
 *       ATCC      T       G
 *       2         4       8
 */

/* make sure POSIX APIs are properly activated */
#if defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE		200112L
#endif

#if defined(__darwin__) && !defined(_BSD_SOURCE)
#  define _BSD_SOURCE
#endif

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define UNITTEST				1
#define DZ_NUCL_2BIT								/* use 2bit encoding for input sequences */
#define DZ_FULL_LENGTH_BONUS						/* use full-length bonus feature */
#define DZ_CIGAR_OP				0x44493d58			/* 'D', 'I', '=', 'X'; the default is 0x04030201 */
#include "dozeu.h"

int main(int argc, char *argv[])
{
	unittest_main(argc, argv);						/* just for testing */

	/* init score matrix and memory arena */
	int8_t const M = 2, X = -3, GI = 5, GE = 1;		/* match, mismatch, gap open, and gap extend; g(k) = GI + k + GE for k-length gap */
	int8_t const xdrop_threshold = 70, full_length_bonus = 10;
	int8_t const score_matrix[16] = {
	/*                 ref-side  */
	/*                A  C  G  T */
	/*               00 01 10 11 */  
	/*        A 00 */ M, X, X, X,
	/* query- C 01 */ X, M, X, X,
	/*  side  G 10 */ X, X, M, X,
	/*        T 11 */ X, X, X, M
	};
	struct dz_s *dz = dz_init(
		score_matrix,
		GI, GE,
		xdrop_threshold,
		full_length_bonus
	);

	/* pack query */
	char const *query = "\x0\x1\x0\x1\x3\x3\x1\x3\x0\x2\x0\x1\x3\x3\x3\x0\x1\x1\x0\x1\x3\x0";
	struct dz_query_s const *q = dz_pack_query_forward(
	    dz,
	    query,						/* char const *seq: query sequence */
	    22							/* length */
	);

	/* init node array */
	struct dz_forefront_s const *ff[10] = { 0 };

	/* fill root: node 0 */
	ff[0] = dz_extend(
	    dz, q,
	    dz_root(dz), 1,				/* struct dz_forefront_s const *ff[], uint64_t n_ffs; incoming forefront and degree; dz_root(dz) for root node */
	    "\x0\x1\x0\x1", 4, 0		/* reference-side sequence, its length, and node id */
	);

	/* branching paths: 1, 2 -> 3 */
	ff[1] = dz_extend(dz, q, &ff[0], 1, "\x3\x3\x2\x3", 4, 1);
	ff[2] = dz_extend(dz, q, &ff[0], 1, "\x0\x3\x1\x1", 4, 2);
	ff[3] = dz_extend(dz, q, &ff[1], 2, "\x0\x2\x0\x1", 4, 3);		/* "&ff[1], 2" indicates ff[1] and ff[2] are incoming nodes */

	/* insertion: -, 4 -> 5 */
	ff[4] = dz_extend(dz, q, &ff[3], 1, "\x3", 1, 4);
	ff[5] = dz_extend(dz, q, &ff[3], 2, "\x3\x3\x1\x3\x0", 5, 5);	/* "&ff[3], 2" indicates ff[3] and ff[4] are incoming nodes */

	/* SNVs: 6, 7, 8 -> 9 */
	ff[6] = dz_extend(dz, q, &ff[5], 1, "\x0", 1, 6);
	ff[7] = dz_extend(dz, q, &ff[5], 1, "\x1", 1, 7);
	ff[8] = dz_extend(dz, q, &ff[5], 1, "\x2", 1, 8);
	ff[9] = dz_extend(dz, q, &ff[6], 3, "\x1\x0\x1\x2\x2", 5, 9);

	/* detect max */
	struct dz_forefront_s const *max = NULL;
	for(size_t i = 0; i < 10; i++) {
		if(max == NULL || ff[i]->max > max->max) { max = ff[i]; }
	}

	/* traceback */
	struct dz_alignment_s const *aln = dz_trace(
		dz,
		max
	);

	printf("ref_length(%u), query_length(%u), score(%d), path(%s)\n", aln->ref_length, aln->query_length, aln->score, aln->path);
	for(size_t i = 0; i < aln->span_length; i++) {
		struct dz_path_span_s const *s = &aln->span[i];
		printf("node_id(%u), subpath_length(%u), subpath(%.*s)\n",
			s->id,
			s[1].offset - s[0].offset,
			s[1].offset - s[0].offset, &aln->path[s->offset]
		);
	}

	/* clean up */
	dz_destroy(dz);
	return(0);
}
