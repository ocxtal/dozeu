
# dozeu.h

All we need is faster implementation of BLAST alignment algorithm.

## Features

* BLAST X-drop DP algorithm ([semi-global "extension" aligmnent of affine-gap penalty](https://github.com/elucify/blast-docs/))
* vectorized with SSE4.1 (run on Core 2 (2007) or later)
* linear-to-graph alignment (first developed for [vg align](https://github.com/vgteam/vg) subcommand)

## Example

```C

```

## Issues and Limitations

* Only supports nucleotide sequences (can be removed but it will be slower, dispatcher would be fine)
* newer instruction support (AVX / AVX2 / AVX-512) is not planned (keep it as simple as possible)

## License

MIT

Copyright 2018, Hajime Suzuki
