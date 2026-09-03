# dozeu.h

SIMD-parallel BLAST X-drop alignment for sequence graphs.

## Update on 2026/9

The changes developed while embedding dozeu in [vg](https://github.com/vgteam/vg) have been merged back into this repository. They refresh the public API and include the fixes used by current vg.

The main changes are:

* X-drop parameters are selected for each alignment with `dz_align_init` and passed to `dz_extend` or `dz_scan`.
* Full-length bonus is selected when packing each query instead of when creating the context.
* Quality-adjusted scoring is available through the `dz_qual_adj_*` API.
* The DP and traceback fixes from vg are included, notably fixes for merged graph paths, X-dropped ranges, full-length bonus, and forward/reverse symmetry.
* The allocator now reserves each variable-sized DP record before writing it. `dz_flush` reuses an arena for another alignment and `dz_trim` releases retained blocks above a requested limit.
* SSE4.1 operations are provided through [SIMDe](https://github.com/simd-everywhere/simde), removing the explicit x86-64 requirement.

## Features

* BLAST X-drop extension with affine gap penalties
* Linear-to-graph alignment, originally developed for vg
* ASCII nucleotide, 2-bit nucleotide, and integer protein inputs
* Per-read full-length bonus and base-quality-adjusted scoring
* Header-only C implementation with C++ linkage support

## Use

Initialize a scoring context, initialize an alignment root, pack a query, and extend the root through graph nodes:

```c
int8_t const score_matrix[16] = {
     2, -3, -3, -3,
    -3,  2, -3, -3,
    -3, -3,  2, -3,
    -3, -3, -3,  2
};

struct dz_s *dz = dz_init(score_matrix, 5, 1);
struct dz_alignment_init_s init = dz_align_init(dz, 65);
struct dz_query_s const *query = dz_pack_query_forward(
    dz, "ACACTTCTAGACTTTACCACTA", 10, 22
);

struct dz_forefront_s const *ff = dz_extend(
    dz, query, &init.root, 1,
    "ACAC", 4, 0, init.xt
);

struct dz_alignment_s const *alignment = dz_trace(dz, ff);
dz_destroy(dz);
```

`DZ_FULL_LENGTH_BONUS` enables the bonus argument used above. Without it, the query packers omit that argument. See `example.c`, `example.2bit.c`, and `example.protein.c` for complete graph examples.

Define `DZ_QUAL_ADJ` before including `dozeu.h` to use the quality-adjusted API. Score matrices are supplied to `dz_qual_adj_init`, and qualities are supplied to the `dz_qual_adj_pack_query*` functions.

Build and run the examples with:

```
make
./example
./example.2bit
./example.protein
```

Call `dz_flush` after an alignment when the context will be reused. After flushing, `dz_trim` can reduce the memory retained by the arena. Quality-adjusted contexts use `dz_qual_adj_flush`.

## Implementation

The query is packed into SIMD-width blocks before alignment. For nucleotide input, the score profile is produced in the DP loop with a byte shuffle instead of loading a precalculated profile. Protein input uses a packed score profile derived from the supplied matrix.

Each graph node is extended from an array of preceding forefronts. Incoming DP vectors are merged cell by cell, then extended across the node's reference sequence. A negative reference length scans the reverse complement from the supplied end pointer.

X-drop is tested per vector. A vector is removed only when all its cells fall outside the threshold. Forefronts keep separate ranges for forward extension and traceback so that a dropped traceback range does not incorrectly change later extension.

DP cells use a biased unsigned 16-bit representation. This retains the original score range while making saturated SIMD arithmetic handle the lower boundary correctly.

The allocator stores variable-sized columns, caps, forefront arrays, and traceback results in linked arena blocks. Streaming reservations ensure that a record is contiguous before it is written. The blocks remain available across `dz_flush` calls unless removed by `dz_trim`.

The maximum-scoring graph leaf is not tracked globally. The caller compares the `max` field of candidate forefronts and passes the chosen forefront to `dz_trace`.

For the BLAST X-drop algorithm, see the [Gapped Alignment Phase](https://github.com/elucify/blast-docs/wiki/Gapped-Alignment-Phase). Graph extension follows a POA-like merge of predecessor DP columns.

![An illustration of SIMD-parallel X-drop alignment to a graph](https://raw.githubusercontent.com/ocxtal/dozeu/master/fig/xdrop.png)

## Limitations

* Scores and DP cells use 16-bit storage.
* The implementation is organized around 128-bit SSE4.1 operations; SIMDe supplies native or emulated implementations.
* The caller owns graph traversal and selection of the maximum-scoring leaf.

## Acknowledgements

The current API, allocator, quality-adjusted scoring, and many correctness fixes were contributed through vgteam's use and maintenance of dozeu. Contributors to these changes include Adam M. Novak, Jordan M. Eizenga, Evan Nemerson, Michael R. Crusoe, and Muaaz Gul Awan.

Erik Garrison described the internals of vg and provided insightful comments about alignment algorithms. Toshiaki Katayama and the DBCLS staff organized the third RDF Summit in Kyoto, where work on dozeu began.

## License

MIT

Copyright 2018-2026, Hajime Suzuki and contributors
