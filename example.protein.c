/**
 * Aligning MSALLILALVGAAVAPLEDDDK to the following graph
 *                         6
 *                         A
 *  0    1    3      5    /7\9
 *  MSAL-LILA-LVGA---AVAFP-E-EDDCL
 *      \    /    \ /     \ /
 *       KKLG      A       M
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

#define DZ_PROTEIN
#define DZ_MAT_SIZE				32
#define DZ_FULL_LENGTH_BONUS						/* use full-length bonus feature */
#define DZ_CIGAR_OP				0x44493d58			/* 'D', 'I', '=', 'X'; the default is 0x04030201 */
#include "dozeu.h"

int main(int argc, char *argv[])
{
	/* init score matrix and memory arena */
	int8_t const GI = 5, GE = 1;				/* match, mismatch, gap open, and gap extend; g(k) = GI + k + GE for k-length gap */
	int8_t const xdrop_threshold = 70, full_length_bonus = 10;
	int8_t const raw_score_matrix[24 * 24] = {
	/*        BLOSUM62 obtained from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt                   */
	/*                                                     ref-side                                               */
	/*              A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X     */
	/*        A */  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,
	/*        R */ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,
	/*        N */ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,
	/*        D */ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,
	/*        C */  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
	/*        Q */ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,
	/*        E */ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
	/*        G */  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,
	/*        H */ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,
	/*        I */ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,
	/* query- L */ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,
	/*  side  K */ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,
	/*        M */ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,
	/*        F */ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,
	/*        P */ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
	/*        S */  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,
	/*        T */  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,
	/*        W */ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,
	/*        Y */ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,
	/*        V */  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,
	/*        B */ -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,
	/*        Z */ -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
	/*        X */  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,
	/*        * */ -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1
	};

	/* reorder BLOSUM62 score matrix so that `residue_in_ascii & 0x1f` to be index */
	char const *indextoascii = "ARNDCQEGHILKMFPSTWYVBZX ";
	int8_t score_matrix[DZ_MAT_SIZE * DZ_MAT_SIZE] = { 0 };
	for(size_t i = 0; i < 24; i++) {
		size_t ii = indextoascii[i] & 0x1f;
		for(size_t j = 0; j < 24; j++) {
			size_t jj = indextoascii[j] & 0x1f;
			score_matrix[ii * DZ_MAT_SIZE + jj] = raw_score_matrix[i * 24 + j];
		}
	}

	/* create context */
	struct dz_s *dz = dz_init(
		score_matrix,
		GI, GE,
		xdrop_threshold,
		full_length_bonus
	);

	/* pack query */
	char const *query = "MSALLILALVGAAVAPLEDDDK";
	struct dz_query_s const *q = dz_pack_query_forward(
	    dz,
	    query,						/* char const *seq: query sequence */
	    strlen(query)				/* length */
	);

	/* init node array */
	struct dz_forefront_s const *ff[10] = { 0 };

	/* fill root: node 0 */
	ff[0] = dz_extend(
	    dz, q,
	    dz_root(dz), 1,             /* struct dz_forefront_s const *ff[], uint64_t n_ffs; incoming forefront and degree; dz_root(dz) for root node */
	    "MSAL", strlen("MSAL"), 0   /* reference-side sequence, its length, and node id */
	);

	/* branching paths: 1, 2 -> 3 */
	ff[1] = dz_extend(dz, q, &ff[0], 1, "LILA", strlen("LILA"), 1);
	ff[2] = dz_extend(dz, q, &ff[0], 1, "KKLG", strlen("KKLG"), 2);
	ff[3] = dz_extend(dz, q, &ff[1], 2, "LVGA", strlen("LVGA"), 3);		/* "&ff[1], 2" indicates ff[1] and ff[2] are incoming nodes */

	/* insertion: -, 4 -> 5 */
	ff[4] = dz_extend(dz, q, &ff[3], 1, "A", strlen("A"), 4);
	ff[5] = dz_extend(dz, q, &ff[3], 2, "AVAFP", strlen("AVAFP"), 5);	/* "&ff[3], 2" indicates ff[3] and ff[4] are incoming nodes */

	/* SNVs: 6, 7, 8 -> 9 */
	ff[6] = dz_extend(dz, q, &ff[5], 1, "A", strlen("A"), 6);
	ff[7] = dz_extend(dz, q, &ff[5], 1, "E", strlen("E"), 7);
	ff[8] = dz_extend(dz, q, &ff[5], 1, "M", strlen("M"), 8);
	ff[9] = dz_extend(dz, q, &ff[6], 3, "EDDCL", strlen("EDDCL"), 9);

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
