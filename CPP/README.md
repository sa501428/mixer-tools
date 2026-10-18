# SLICe in C++

Native C++17 implementation of the active Mixer Tools Java `slice` pipeline,
located alongside its Java reference. It reads `.hic` files through straw C++;
Java and Python are not runtime dependencies. It does not call the neighboring
`slice-cpp-latest` project or use its newer regional/count model.

The source baseline is Mixer Tools commit
`863f1f42738c52f25a167eefc6cd31b1321abdc4`. See [PORTING.md](PORTING.md)
for the Java-to-C++ map, exact defaults, validation, and compatibility differences.
The numerical calling pipeline is tested against the actual Java code, including
multiple inputs, k-means, k-medians and ENCODE export. Diagnostic UMAP coordinates
use a native implementation and are **not numerically identical to tagbio UMAP**.

## Build

Requirements: CMake 3.16+, a C++17 compiler, libcurl, zlib, zstd, and the local
straw reader. The default straw path matches this workspace layout:

```text
hic-algos/subcompartments/mixer-tools/CPP
hic-infra/straw/C++
```

From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

For another layout, add `-DSTRAW_SOURCE_DIR=/absolute/path/to/straw/C++`.
The straw tree must provide `straw.h`, `straw.cpp`, `straw_v10.cpp`, `hbs.cpp`,
and `norm_vectors.cpp`. This project links a native reader library and does not
build straw's CLI. To build without the test dependencies, use
`-DBUILD_TESTING=OFF`. Python 3 is needed for integration tests; a JDK enables
the additional differential tests against the adjacent Java source and JARs.
No test downloads data.

Optional installation:

```sh
cmake --install build --prefix /your/install/prefix
```

## Run

```sh
# Java-compatible calling defaults: 100 kb, k-medians, automatic stage norms.
./build/slice sample.hic 2,10,3 output sample

# The optional slice subcommand is also accepted.
./build/slice slice -r 100000 -c 1,2,3 -z 8 --random-seeds 42 \
  sample.hic 2,10,3 output sample

# Each input contributes its own independently compressed feature block.
./build/slice -z 8 first.hic,second.hic 2,10,3 output pooled

# ENCODE BED with all tested clusterings (include K=5 for reliable selection).
./build/slice --encode-mode sample.hic 2,10,3 output sample

# Expose the existing Java k-means switch, or run both methods.
./build/slice --kmeans sample.hic 2,10,3 output means
./build/slice --both sample.hic 2,10,3 output both

# Run the dormant Java-style translocation helper and apply its rectangles.
./build/slice --has-translocation sample.hic 2,10,3 output screened
```

`K0,KF,nK` means minimum K, **exclusive** upper bound, and successful attempts
per K. For example `2,10,3` runs K=2 through K=9, three attempts each, choosing
the lowest Java WCSS score for each K. It does not select a single final K.

## Pipeline

1. Select autosomes and build a shared bad-bin mask across all inputs. Using
   cis `SCALE` (fallback `KR`) normalization on trans contacts, compute average
   `log(1 + contact)` coverage and fraction of nonzero trans entries. Exclude
   empty bins, log coverage z > 5 and nonzero coverage z < -3.
2. Independently for each input and chromosome, build cis O/E at
   `max(requested_resolution, 100000)`. Smooth the stored expected-distance
   vector with Java's 5 Mb rolling-median procedure; use pseudocount 1 in both
   numerator and denominator. Exclude shared bad bins, the first three
   diagonals, zeros and absolute O/E > 30.
3. Compute pairwise Pearson correlations using jointly observed coordinates,
   apply the bounded Fisher transform, and mask separations < 5 Mb. Retain
   candidate rows with Java's `< .7` missing-value test. Run 20-iteration
   k-means with `floor(chromosome_length / 50 Mb) + 7` initial centroids, keep
   groups larger than two rows, and recompute mean centroids. Assign bins to
   the most correlated centroid only if correlation > .2. Propagate those
   groups to the requested resolution. The Java singleton-filter switch is
   false, so singleton assignments remain.
4. Sum trans contacts into chromosome-specific compression groups with
   `INTER_SCALE` (fallback `INTER_KR`). Same-chromosome feature blocks are NaN.
   Keep uncompressed genomic rows, and concatenate the feature blocks and
   weights from all inputs. Inputs require matching genomes and chromosomes.
5. Convert zeros to NaN, divide compressed sums by their group sizes, and
   apply `log1p`. Remove values above global z=5, recompute statistics, remove
   values outside [-2,2], and invert with `expm1`. Drop rows with more than
   `ceil(.7 * number_of_columns)` missing values. Compute positive-value
   column means and population standard deviations, z-score columns, and
   multiply by the square root of their compression-group sizes.
6. For diagnostics, generate the Java 50-centroid Fisher-transformed cosine
   feature matrix, then a 2D native UMAP projection. It does not enter calling.
7. Run k-means if requested. Before k-medians, apply the **second** square-root
   column-weight scaling used by Java to accommodate its L1 distance. Use
   200 iterations per fit, the Java RNG seed hierarchy, farthest-first
   initialization, missing-aware distances and mean/median centers. Sort
   resulting clusters by decreasing size, breaking ties by first member.
8. Retain the lowest-score successful attempt for every K. Include the Java
   penalty for clusters smaller than five bins. Export WCSS, Java's AIC/BIC
   expressions, the disabled-silhouette value 0 and maximum off-diagonal
   centroid correlation. Merge adjacent same-label bins into BED9 intervals.
9. In ENCODE mode, select the highest tested K with maximum centroid
   correlation in (0,.8), subject to the K=5 lower bound when present. Export
   the selected calls with all alternative K assignments, merging only where
   every assignment agrees. These are numeric labels, not inferred A1/B1 names.

## Options

| Option | Default / behavior |
| --- | --- |
| `-r`, `--resolutions LIST` | 100000 BP; first value used, as in Java. `--resolution` is an alias. Resolutions finer than 100 kb must divide 100000 exactly. Both requested and compression resolutions must be stored in every input. |
| `-c`, `--chromosomes LIST` | All autosomes; accepts names with or without `chr`. At least two selected autosomes are required. |
| `-z`, `--threads N` | 1; zero means 1, negative means hardware concurrency. |
| `--random-seeds LIST` | 22871; signed 64-bit Java seeds, last value takes effect. |
| `--subsample N` | 100; accepted and stored like Java's `NUM_PER_CENTROID`, which is unused by the active pipeline. |
| `--encode-mode` | Disabled; enables the additional ENCODE BED. |
| `-v`, `--verbose` | Print normalization/seed diagnostics and save cis O/E, correlation and centroid matrices for inspection. |
| `-h`, `--help`; `-V`, `--version` | Help and native version. |
| `--kmeans` | C++ extension exposing Java's `USE_KMEANS`; runs only k-means. |
| `--both` | C++ extension; k-means then k-medians, preserving Java's scaling order. |
| `--no-filter-outliers` | C++ extension exposing `FILTER_OUTLIERS=false`. |
| `--no-umap` | C++ extension; skips cosine projection and UMAP, leaving calls unchanged. Useful for small maps or large runs. |
| `--has-translocation` | C++ extension connecting the otherwise dormant Java helpers; see below. |
| `-k`, `--normalization NAME` | C++ extension overriding stage normalizations. Java's active `Slice.java` ignores this global parser option. Omit for faithful stage defaults. |
| `--max-retries N` | C++ safety extension; default 100 retries per successful attempt instead of Java's unbounded retry loop. |

The global Mixer parser also accepts `-m/--matrix-window-width`, `-w/--window`,
`--seed`, `--compare`, `--type`, and `--corr`. Java SLICe does not consume these;
the C++ parser accepts them with a warning. Use `--random-seeds` to set the seed.

There is **no active Java CLI switch to combine cis and trans features**.
Cis establishes the compression groups; trans supplies the final clustering
features. The unused combined-normalization method and disabled correlation/GMM
experimental branch are mapped in PORTING.md; they are not presented as working
Java CLI options.

## Translocation extension

The default Java SLICe path does not instantiate `SimpleTranslocationFinder` or
apply `TranslocationSet`, so the default native run does neither. Opting into
`--has-translocation` runs the helpers' logic after coverage filtering:

- Work at the coarsest stored BP resolution using GW_SCALE, falling back to GW_KR.
- Map shared bad bins to the coarse grid; compute the mean and sample standard
  deviation of recorded valid `log1p(contact)` entries across trans pairs.
- Mark z > 6 pixels; dilate by one bin; bridge gaps up to ten bins; use Java's
  expanding-rectangle extraction and retain rectangles with area > 20 cells.
- Apply rectangles only to the affected chromosome pair. Keep their contacts
  out of compressed sums, and adjust compressed means for eligible partner bins
  using the union of overlapping rectangles. Cis information is retained.

This is a heuristic suspected-translocation mask, not structural variant calling.
It is an explicitly documented extension because the active Java command does
not connect detection to compression. It does not rerun the shared coverage mask.

## Outputs

| File | Contents |
| --- | --- |
| `PREFIX_K_kmedians_clusters.bed` | Collapsed BED9 calls for each K, labels starting at 1 and original Java palette. |
| `PREFIX_K_kmeans_clusters.bed` | Corresponding k-means calls when enabled. |
| `kmedians_cluster_post_stats.npy`, `kmeans_cluster_post_stats.npy` | Float64, shape `(6, number_of_K_values)`; rows are K, WCSS, AIC, BIC, silhouette=0, maximum centroid correlation. |
| `slice_matrix.npy` | Float32 cleaned weighted matrix **before** the second k-medians scaling. |
| `genome_indices.npy` | Int32 `(chromosome_index, start, end)` aligned with matrix rows. |
| `initial_split.bed`, `problems.bed` | Input 0 compression assignments and unassigned/bad bins. Additional inputs use `_1`, `_2`, etc. before `.bed`. |
| `corr_slice_matrix.npy` | Float32 native port of Java's centroid-cosine diagnostic features. |
| `umap_projection.npy` | Float32 native 2D diagnostic coordinates. |
| `umap_chrom.svg`, `umap_PREFIX_K_METHOD_clusters.svg` | Native diagnostic scatter plots. Java produces PNG plots. |
| `slice_subcompartment_clusters.bed` | ENCODE mode, BED9 plus alternative clustering columns. |
| `suspected_translocations_RES_INPUT.bedpe`, `all_suspected_translocations_RES_INPUT.bedpe` | Opt-in retained rectangles and individual z > 6 pixels for each input. |
| `cis_oe_CHROM.npy`, `cis_corr_CHROM.npy`, `cis_centroids_CHROM.npy` | Verbose-only diagnostics, overwritten by later inputs as in Java's diagnostic convention. |

Missing/filtered bins remain gaps in BED outputs. Call intervals are clipped to
chromosome length. All calling matrices remain dense as in Java; memory scales
with the largest cis matrix and the genome-wide compressed matrix. Native UMAP
uses exact neighbors and can be expensive for large genomes. Use `--no-umap`
when only calls are required. Degenerate/sparse inputs fail with an explanation;
a failed run may leave useful intermediate outputs.
