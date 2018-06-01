/**
 * Mapping ACACTTCTAGACTTTACCACTA to the following graph
 *                         6
 *                         A
 *  0    1    3      5    /7\9
 *  ACAC-TTGT-AGAC---TTCTA-C-CACTA
 *      \    /    \ /     \ /
 *       ATCC      T       G
 *       2         4       8
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define DZ_CIGAR_OP				0x44493d58		/* 'D', 'I', '=', 'X'; the default is 0x04030201 */
#include "dozeu.h"

int main(int argc, char *argv[])
{
	/* init score matrix and memory arena */
	int8_t const M = 2, X = -3, GI = 5, GE = 1;		/* match, mismatch, gap open, and gap extend; g(k) = GI + k + GE for k-length gap */
	int8_t const xdrop_threshold = 70;
	int8_t const score_matrix[16] = {
	/*              ref-side  */
	/*             A  C  G  T */
	/*        A */ M, X, X, X,
	/* query- C */ X, M, X, X,
	/*  side  G */ X, X, M, X,
	/*        T */ X, X, X, M
	};
	struct dz_s *dz = dz_init(
		score_matrix,
		GI, GE,
		xdrop_threshold
	);

	/* pack query */
	char const *query = "ACACTTCTAGACTTTACCACTA";
	struct dz_query_s const *q = dz_pack_query_forward(
	    dz,
	    query,						/* char const *seq: query sequence */
	    strlen(query)				/* length */
	);

	/* init node array */
	struct dz_forefront_s const *ff[10] = { 0 };

	/* fill root */
	ff[0] = dz_extend(
	    dz, q,						/* context and query */
	    dz_root(dz), 1,				/* struct dz_forefront_s const *ff[], uint64_t n_ffs; incoming forefront and degree; dz_root(dz) for root node */
	    "ACAC", strlen("ACAC"), 0	/* reference-side sequence, its length, and node id */
	);

	/* branching path */
	ff[1] = dz_extend(dz, q, &ff[0], 1, "TTGT", strlen("TTGT"), 1);
	ff[2] = dz_extend(dz, q, &ff[0], 1, "ATCC", strlen("ATCC"), 2);
	ff[3] = dz_extend(dz, q, &ff[1], 2, "AGAC", strlen("AGAC"), 3);	/* "&ff[1], 2" indicates ff[1] and ff[2] are incoming nodes */

	/* deletion */
	ff[4] = dz_extend(dz, q, &ff[3], 1, "T", strlen("T"), 4);
	ff[5] = dz_extend(dz, q, &ff[3], 2, "TTCTA", strlen("TTCTA"), 5);	/* "&ff[3], 2" indicates ff[3] and ff[4] are incoming nodes */

	/* SNVs */
	ff[6] = dz_extend(dz, q, &ff[5], 1, "A", strlen("A"), 6);
	ff[7] = dz_extend(dz, q, &ff[5], 1, "C", strlen("C"), 7);
	ff[8] = dz_extend(dz, q, &ff[5], 1, "G", strlen("G"), 8);
	ff[9] = dz_extend(dz, q, &ff[6], 3, "CACTA", strlen("CACTA"), 9);

	/* traceback */
	struct dz_alignment_s const *aln = dz_trace(
		dz,
		ff[9]
	);

	printf("ref_length(%u), query_length(%u), path(%s)\n", aln->ref_length, aln->query_length, aln->path);
	for(size_t i = 0; i < aln->span_length; i++) {
		struct dz_path_span_s const *s = &aln->span[i];
		printf("node_id(%u), subpath_length(%u), subpath(%.*s)\n",
			s->id,
			s[1].offset - s[0].offset,
			s[1].offset - s[0].offset, &aln->path[s->offset]
		);
	}
	return(0);
}
