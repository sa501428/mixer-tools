# Java reference and compatibility audit

Reference: Mixer Tools commit `863f1f42738c52f25a167eefc6cd31b1321abdc4`.
All paths in the table are relative to the Mixer Tools repository root.
Original MIT attribution is retained in `LICENSE`.

## Source map

| Java reference | Native implementation |
| --- | --- |
| `src/mixer/algos/Slice.java` and `src/mixer/clt/CommandLineParserForMixer.java` | `src/main.cpp`, `include/slice/slice.hpp`, `src/reader.cpp`: positional syntax, active options, defaults and stage normalizations. |
| `utils/slice/kmeans/FullGenomeOEWithinClusters.java` | `src/pipeline.cpp`: seed hierarchy, input concatenation, K sweep, successful attempts, scaling order, exports. All remaining Java paths below are under `src/mixer/`. |
| `utils/slice/cleaning/BadIndexFinder.java`, `GWRegionStatistics.java`, `CRecordUtils.java` | `src/reader.cpp`: pooled coverage and shared masks. |
| `utils/slice/cleaning/IndexOrderer.java`, `IntraMatrixCleaner.java`, `NearDiagonalTrim.java`, `QuickCentroids.java` | `src/reader.cpp`, `src/numerics.cpp`: O/E, missing masks, correlations, compression clustering, centroid means, threshold assignments and resolution propagation. |
| java-straw `HiCFileTools`, `ExtractingOEDataUtils`, `ExpectedValueFunctionImpl`, `QuickMedian` | `src/reader.cpp`: pseudocount O/E and corrected expected vectors. `.hic` access uses the native straw reader. |
| `utils/slice/matrices/SliceMatrix.java`, `CompositeGenomeWideMatrix.java` | `src/reader.cpp`, `src/pipeline.cpp`: uncompressed rows, compressed columns, structural missing blocks, weights and independent per-input feature concatenation. |
| `utils/slice/cleaning/SliceMatrixCleaner.java`, `cleaning/utils/RowCleaner.java`, `DimensionCleaner.java`; `utils/common/LogTools.java`, `ParallelizedStatTools.java`, `ZScoreTools.java` | `src/numerics.cpp`: log-space filtering, sparsity threshold, positive-value column statistics, z-scores and weights. |
| `RobustConcurrentKMeans.1.2.jar`: `SmartInitialization`, `Worker`, `ProtoCluster`, `QuickMedian`, robust L1/L2 metrics | `src/numerics.cpp`: Java RNG, farthest-first initialization, float-rounded assignment distances, first-minimum tie breaking, missing-aware means and medians, removal of empty clusters. Bundled bytecode was inspected to resolve dependency details. |
| `utils/slice/kmeans/GenomeWideKmeansRunner.java`, `ClusterTools.java`, `KmeansResult.java`, `KmeansEvaluator.java` | `src/pipeline.cpp`: sorted cluster IDs, small-group score penalty, WCSS divided by K, AIC/BIC expressions, zero silhouette, maximum centroid correlation. |
| `utils/slice/EncodeExportUtils.java`, `structures/ENCODESubcompartmentInterval.java` | `src/pipeline.cpp`: selected-K logic, alternative IDs, all-label collapse and historical output fields. |
| `utils/slice/structures/SliceUtils.java`, `SubcompartmentInterval.java`, `SubcompartmentColors.java` | `src/reader.cpp`, `src/pipeline.cpp`: BED9, colors, contiguous interval merging and chromosome coordinates. |
| `utils/slice/cleaning/SimilarityMatrixTools.java`, `utils/similaritymeasures/RobustCosineSimilarity.java` | `src/projection.cpp`: 50 k-medians groups, groups > 5 rows, mean prototypes, Fisher-transformed cosine features. |
| `utils/umap/UmapProjection.java`; bundled `umap-java.jar` | `src/projection.cpp`: same diagnostic role and configured neighbors/min-distance; native exact-neighbor fuzzy graph and negative-sampling optimization, with a different initialization/RNG trajectory. SVG plotting replaces Java's PNG plotting. |
| `utils/slice/cleaning/SimpleTranslocationFinder.java`, `Connected2DComponent.java`, `TranslocationSet.java` | `src/translocations.cpp`: coarse GW-normalized log-tail screen, dilation, bridging, expanding rectangles. Pair-specific application to compression is a C++ extension. |

## Which Java branches actually run?

`Slice.java` defaults to `USE_KMEANS=false`, `USE_KMEDIANS=true`,
`FILTER_OUTLIERS=true`, `PROJECT_TO_UMAP=true`, `USE_WEIGHTED_MEAN=false`, and
`USE_INTER_CORR_CLUSTERING=false`. The native defaults preserve that selection.
`--kmeans`, `--both`, `--no-filter-outliers`, and `--no-umap` expose selected
source-level switches without requiring recompilation.

Several classes in `utils/slice` are not on that active call path:

- `SimpleTranslocationFinder` is not constructed by `Slice` or its active matrix
  builders. The native `--has-translocation` extension explicitly connects it.
- `CompositeGenomeWideMatrix.makeComboNormMatrix` is unused. It concatenates
  INTER- and GW-normalized **trans** matrices; it is not a cis/trans switch.
- `CorrMatrixClusterer`, the GMM implementation, imputation, and archived cleaners
  are not reachable from the Java CLI with its checked-in flags. They are not
  included as native calling modes. In particular, `--corr` is a global parser
  option ignored by Java `Slice.java`, not a way to enable that GMM branch.
- Distance-based `OutlierCleaner` calls are commented out in `DimensionCleaner`.
  The active outlier processing is the global log-space thresholding in
  `SliceMatrixCleaner`; the native code preserves that distinction.
- `NUM_PER_CENTROID` is assigned by `--subsample` but is not consumed by the
  active SLICe pipeline. Hard-coded 50 diagnostic centroids remain 50.

## Compatibility choices

The native implementation preserves the active calling algorithm rather than
claiming byte-identical behavior for all JVM library outputs. Known differences:

1. **UMAP** uses native exact Euclidean neighbors, smooth-kNN scaling, fuzzy union
   and negative-sampling optimization with min-distance .5, spread 1 and seed 0.
   It uses random initialization rather than the bundled library's embedding
   implementation and may produce different geometry. This is diagnostic only;
   neither its coordinates nor the SVG color/plot layout enter calling.
2. **Interval collapse:** Java's loop can emit an extra first-bin interval in
   addition to its merged interval. Native BED files contain each called span
   once. Differential tests compare genomic calls, not that redundant record.
3. **Failure handling:** missing pairs, missing resolutions/vectors, impossible
   K, empty matrices and failed ENCODE selection produce errors instead of
   JVM exits, unbounded retries, array-index exceptions or waiting indefinitely.
   Retries default to 100. Constant-feature columns become missing rather than
   allowing infinities into distance calculations.
4. **Chromosome boundaries:** native interval construction consistently excludes
   coverage-filtered empty terminal bins and clips call endpoints. It avoids
   Java's conflicting `ceil(length/resolution)` versus `floor(...)+1` allocation
   rules and associated boundary indexing hazards.
5. **Multiple-input diagnostics:** all inputs contribute separate feature blocks
   as in Java. Native exports preserve input-specific initial BED diagnostics;
   Java overwrites its initial/problematical-bin files for later inputs.
6. **Normalization override / translocation masking:** these are documented
   extensions, not enabled in the default Java-compatible calling path. The
   default normalization ordering is SCALE/KR, INTER_SCALE/INTER_KR,
   GW_SCALE/GW_KR. Advertised label precedence follows `NormalizationPicker`;
   actual vector availability is validated before reading contacts.
7. **Numerics/concurrency:** float intermediates, missing masks, Java Random and
   seed nesting are retained. Reductions use stable genomic order. Assignment
   and independent correlations run in parallel. Floating-point contraction is
   disabled for calling code. Java's parallel global reductions can still cause
   small cross-platform rounding differences; no whole-genome bitwise guarantee
   is implied by the synthetic differential tests.
8. **Memory:** cis correlations and compressed features are dense, like Java.
   UMAP neighbor construction is exact and quadratic in rows; it uses O(N*50)
   retained graph storage plus per-worker O(N) scratch. The translocation bitmap
   has a checked 100-million-cell limit. Large/sparse maps should be run at
   appropriate resolutions; `--no-umap` removes the optional projection cost.

## Validation

`ctest` runs three suites when a JDK is available:

- `numerics`: Java RNG reference, missing-aware distances/correlation, mean and
  median centers, and degenerate-correlation behavior.
- `integration`: an independently written synthetic `.hic` v9 fixture, both
  methods, K sweep, ENCODE, matrix types, all-input feature concatenation,
  SCALE/KR fallback and explicit override, selected chromosomes, exact terminal
  lengths, opt-in translocation detection, deterministic worker counts,
  projection independence, and malformed/missing input errors.
- `java_parity`: compiles a temporary copy of the adjacent Java reference with
  its bundled JARs. Disables only diagnostic UMAP for core comparisons, then
  exposes the existing `USE_KMEANS` switch for a second comparison. Checks
  compression assignments, cleaned matrix/mask, genomic row coordinates, both
  clustering methods, all six score rows, multi-input behavior and ENCODE output.
  A third reference run with diagnostic UMAP enabled checks the cosine-feature
  matrix. Native UMAP coordinates are deliberately not treated as JVM parity.

The synthetic reference has 360 rows across three autosomes. Calling matrices
and k-medians scores matched exactly on the initial reference case; automated
comparisons allow 2e-6 absolute/relative rounding tolerance. These tests establish
implementation agreement on controlled data, not biological validation on a
real genome. Address/undefined-behavior sanitizer builds can run the native tests
with the same CMake source tree. On this macOS 26.6.2 / Xcode 17 host, the
AddressSanitizer runtime deadlocked during shadow-memory initialization before
`main` (confirmed by a process stack sample); that run did not validate the
program. The separate undefined-behavior sanitizer build passed both the native
numerical and end-to-end integration suites.

Bundled dependency SHA-256 values used for the Java reference:

```text
RobustConcurrentKMeans.1.2.jar
54e7b031d5d1a014a629aa56a010416685836f1a83eaad80f145813c95e26a06
java-straw.2.21.00.jar
442037e8aa74bf5659acbe0ca2789002a66c1b7cc6dfe0c9b74a94af1de67f6c
umap-java.jar
6892f0d1e458e4f364bd3f9819f351f185f00283457ff3ae6a5f1d7426cbfee3
```
