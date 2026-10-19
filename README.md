# MIXER Tools

**Machine Intelligence for Xtracting Essential Regions**  
**Machine Intelligence for Xtracting Enriched Relationships**

MIXER Tools is a comprehensive Java-based command-line toolkit for analyzing Hi-C contact maps and identifying chromatin subcompartments, inter-chromosomal interactions, and long-range enrichments. The toolkit provides advanced algorithms for clustering, dimensionality reduction, network analysis, and statistical evaluation of Hi-C data.

## Overview

MIXER Tools processes Hi-C data files (`.hic` format) to extract meaningful biological insights about chromatin organization. It supports multiple normalization methods, parallel processing, and various analytical workflows for studying 3D genome architecture.

**Version:** 5.03.00

## Features

- **Subcompartment identification** via clustering algorithms (SLICE/DICE)
- **Long-range enrichment detection** (SHUFFLE)
- **CHIC scoring** (Calculated Heuristic for Inspecting Clusters) for cluster quality assessment
- **Dimensionality reduction** using UMAP
- **Network analysis** for chromatin interactions
- **Multi-threaded processing** for improved performance
- **Multiple normalization methods** (KR, VC, VC_SQRT, etc.)
- **Flexible resolution and chromosome selection**

## Tools

### SLICE (Subcompartment Landscape Identification via Clustering Enrichments)

Identifies chromatin subcompartments by clustering enrichment patterns in Hi-C contact maps. Uses k-means or k-medians clustering to partition the genome into distinct subcompartments based on interaction patterns.

**Usage:**
```bash
mixer_tools slice [-r resolution] [-k normalization] [-v] <file.hic> <K0,KF,nK> <outfolder> <prefix>
```

**Parameters:**
- `K0` - Minimum number of clusters
- `KF` - Maximum number of clusters
- `nK` - Number of times to rerun k-means
- `-r, --resolutions` - Resolution in base pairs (default: 100000)
- `-k, --normalization` - Normalization type (NONE, VC, VC_SQRT, KR, GW_KR, GW_VC, INTER_KR, INTER_VC)
- `-v, --verbose` - Verbose output mode

**Example:**
```bash
mixer_tools slice -r 100000 -k KR input.hic 2,10,5 output_dir prefix
```

### DICE (Differentiating Inter-Chromosomal Enrichments)

Similar to SLICE but focuses on inter-chromosomal interactions to differentiate enrichment patterns between chromosomes.

**Usage:**
```bash
mixer_tools dice [-r resolution] [-k normalization] [-v] <file.hic> <K0,KF,nK> <outfolder> <prefix>
```

### SHUFFLE (Scoring Heuristics Used For Finding Long-range Enrichments)

Evaluates long-range enrichments by comparing observed contact patterns with shuffled/null models. Provides statistical scoring to identify significant long-range interactions. The shuffle analysis generates **CHIC scores** (Calculated Heuristic for Inspecting Clusters) to quantify how well subcompartment clusters are preserved in the contact maps.

**Usage:**
```bash
mixer_tools [intra-]shuffle-umap [-r resolution] [-k normalization] [-w window] [-v] <file.hic> <subcompartment.bed(s)> <outfolder> <prefix>
```

**Parameters:**
- `-w, --window` - Window size for analysis
- `<subcompartment.bed(s)>` - BED file(s) containing subcompartment annotations (comma-separated for multiple files)
- `<prefix>` - Output file prefix (comma-separated for multiple prefixes)

**Example:**
```bash
mixer_tools intra-shuffle-umap -r 100000 -k KR -w 20 input.hic compartments.bed output_dir analysis
```

### CHIC (Calculated Heuristic for Inspecting Clusters)

CHIC score is a quantitative metric generated during shuffle analysis that evaluates the quality and preservation of subcompartment clusters in Hi-C contact maps. It compares baseline (observed) contact patterns with shuffled (null model) patterns to assess cluster integrity.

**How CHIC Score Works:**

1. **Baseline Calculation**: Computes statistics on observed contact patterns within defined subcompartment clusters
2. **Shuffled Calculation**: Computes the same statistics on shuffled/null model contact patterns
3. **Score Calculation**: CHIC score = Baseline / Shuffled (higher scores indicate better cluster preservation)

**CHIC Score Metrics:**

The CHIC score is calculated using two complementary metrics:

- **Variation Score**: Measures the variance of contact frequencies within cluster regions. Lower variance in shuffled data relative to baseline indicates better cluster preservation.
- **Relative Entropy Score**: Uses Kullback-Leibler (KL) divergence to measure information divergence between observed and expected contact distributions. Higher KL divergence in baseline vs shuffled indicates stronger cluster structure.

**Output:**

CHIC scores are saved in `aggregate_scores_<prefix>.txt` files, which include:
- Individual CHIC scores for each chromosome or map type
- Baseline and shuffled values for both metrics
- Final CHIC score ratios
- Geometric mean of CHIC scores across all analyzed regions

**Interpretation:**

- **CHIC Score > 1**: Indicates that baseline (observed) patterns show stronger cluster structure than shuffled patterns - clusters are well-preserved
- **CHIC Score ≈ 1**: Suggests that observed patterns are similar to shuffled patterns - weak or no cluster structure
- **CHIC Score < 1**: Rare, but would indicate shuffled patterns show stronger structure than observed - may indicate data quality issues

Higher CHIC scores indicate more robust and biologically meaningful subcompartment clusters.

### UMAP (Uniform Manifold Approximation and Projection)

Performs dimensionality reduction on Hi-C contact maps using UMAP to visualize and analyze high-dimensional interaction data.

**Usage:**
```bash
mixer_tools [intra-]umap [-r resolution] [-k normalization] [-w window] [-v] <file.hic> <subcompartment.bed(s)> <outfolder> <prefix>
```

### Network

Analyzes chromatin interaction networks, identifying connected components and loop anchors in Hi-C data.

**Usage:**
```bash
mixer_tools network [-r resolution] [-k normalization] [-v] <file.hic> <loop_anchors.bedpe> <outfolder>
```

### ABA (Aggregate Block Analysis)

Performs aggregate block analysis on specified genomic regions, computing average contact frequencies across multiple regions.

**Usage:**
```bash
mixer_tools aba [-w window] [-r resolution(s)] [-c chromosomes] [-k normalization] [--threads num_threads] <hicFile(s)> <BedFile> <SaveFolder>
```

**Parameters:**
- `-w, --window` - Window size (default: 20)
- `-r, --resolutions` - Resolution(s) in base pairs (comma-separated, default: 100)
- `-c, --chromosomes` - Chromosome(s) to analyze (comma-separated)
- `--threads` - Number of CPU threads to use

**Example:**
```bash
mixer_tools aba -w 20 -r 10000,25000,50000 -c 1,2,3 -k KR --threads 4 input.hic regions.bed output_dir
```

### Walk (Chromosome Walker)

Performs chromosome walking analysis to identify linked regions across the genome assembly.

**Usage:**
```bash
mixer_tools walk [-r res1,res2] [-k normalization] [-w minimum_size] <hicFile> <assembly_file> <output_directory>
```

### FineTune

Fine-tunes clustering parameters and optimization settings for subcompartment identification.

**Usage:**
```bash
mixer_tools finetune [options] <file.hic> <outfolder>
```

### Rename

Utility tool for renaming or reformatting output files and annotations.

**Usage:**
```bash
mixer_tools rename [options] <input> <output>
```

## Common Options

Most commands support the following global options:

- `-h, --help` - Print help message
- `-v, --verbose` - Enable verbose output mode
- `-V, --version` - Print version information
- `-k, --normalization` - Normalization type: `NONE`, `VC`, `VC_SQRT`, `KR`, `GW_KR`, `GW_VC`, `INTER_KR`, or `INTER_VC` (default: KR)
- `-r, --resolutions` - Resolution(s) in base pairs (comma-separated for multiple)
- `-c, --chromosomes` - Chromosome(s) to process (comma-separated, e.g., `1,2,X`)
- `-z, --threads` - Number of CPU threads (use `-1` for all available cores)
- `-m, --matrix-window-width` - Matrix window width parameter

## Requirements

- **Java** (JRE/JDK 8 or higher)
- **Hi-C data files** in `.hic` format
- **JavaStraw library** (included in dependencies)

## Building

The project uses Java and includes JAR dependencies in the `lib/` directory:

- `commons-math3-3.6.1.jar` - Apache Commons Math library
- `jargs.jar` - Command-line argument parsing
- `java-straw.2.21.00.jar` - Hi-C file reading (JavaStraw)
- `RobustConcurrentKMeans.1.2.jar` - Robust k-means clustering
- `umap-java.jar` - UMAP implementation

To compile:
```bash
javac -cp "lib/general/*:lib/umap/*" -d build src/mixer/**/*.java
```

To create a JAR:
```bash
jar cfm mixer-tools.jar src/META-INF/MANIFEST.MF -C build .
```

## Input File Formats

### Hi-C Files (`.hic`)
Standard Hi-C contact map format as generated by Juicer or similar tools. The toolkit uses JavaStraw to read these files.

### BED Files
Standard BED format for genomic regions:
```
chr1    1000    2000    region1    0    +
```

### BEDPE Files
Paired-end BED format for loop anchors or interactions:
```
chr1    1000    2000    chr1    5000    6000    interaction1    0    +    +
```

## Output

Output files vary by tool but typically include:
- Clustering assignments (BED format)
- Contact matrices
- Statistical scores
- Visualization data
- UMAP projections

## Examples

### Basic subcompartment identification:
```bash
mixer_tools slice -r 100000 -k KR -v sample.hic 2,10,5 output slice_output
```

### Analyze specific chromosomes:
```bash
mixer_tools slice -r 100000 -c 1,2,3 -k KR sample.hic 2,10,5 output slice_output
```

### Multi-threaded analysis:
```bash
mixer_tools aba -r 10000,25000 --threads 8 sample.hic regions.bed output_dir
```

### UMAP projection with subcompartment annotations:
```bash
mixer_tools intra-umap -r 100000 -k KR sample.hic compartments.bed output umap_projection
```

## Citation

If you use MIXER Tools in your research, please cite:

```
MIXER Tools: Machine Intelligence for Xtracting Essential Regions
Version 5.03.00
Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## Acknowledgments

- **Rice University**
- **Baylor College of Medicine**
- **Aiden Lab**

## Additional Resources

- [Code of Conduct](CODE_OF_CONDUCT.md)
- [Contributing Guidelines](CONTRIBUTING.md)

## Support

For issues, questions, or contributions, please open an issue on the GitHub repository.

---

**Note:** This toolkit is designed for research use. Ensure you have appropriate computational resources for large Hi-C datasets, as processing can be memory and CPU intensive.

## Native C++ SLICe

A native C++17 port of the active Java SLICe pipeline is in [CPP](CPP/README.md).
It includes independently compressed multi-input features, the default
k-medians path, optional k-means, ENCODE export, native diagnostic projection,
and an opt-in translocation mask. See [the porting audit](CPP/PORTING.md) for
source mappings, defaults, Java differential tests, and compatibility details.

## Browser-ready subcompartment tracks

Call tracks use BED9: chromosome, zero-based start, exclusive end, cluster name,
score `0`, strand `.`, thick start/end equal to the interval, and `itemRgb`.
Touching intervals on the same chromosome with the same final label are merged;
gaps, different labels, and chromosome boundaries remain separate. Per-bin
analytical tables keep their original granularity. The first 20 states use a
fixed high-contrast categorical palette, with additional colors beyond 20 rather
than repeating the palette. Explicit unassigned state `0` is neutral gray.
Numeric state IDs have the same colors across the tools; this does not imply
that independently fitted states are biologically equivalent. String state
names receive distinct colors in sorted name order within each track.
For UCSC custom tracks, enable `itemRgb="On"` in the track configuration.
Existing result directories are not migrated; these changes apply to new runs.

In ENCODE mode, `slice_subcompartment_clusters.bed` is the merged BED9 track
for the selected clustering. Alternate-clustering metadata is retained in
`slice_subcompartment_clusterings.tsv`, whose first nine columns describe the
interval and whose final three columns retain the existing ENCODE metadata.
Those custom columns are not standard BED12 block fields.

## Cis coverage exclusion

Every calling run applies an independent cis coverage mask before clustering or
row pooling. At the analysis resolution, read raw (`NONE`) cis contacts and count
positive off-diagonal pixels incident on each row; the self-diagonal is excluded.
For multiple inputs, sum each input's nonzero observations per row. Compute the
mean and population standard deviation separately for each chromosome, using
nonempty rows. Always exclude empty rows, and exclude nonempty rows with
`z = (row_nonzeros - mean) / sd < -2.5`. With zero variance, retain the nonempty
rows and still exclude empty rows. The cutoff is enabled by default, independently
of normalization-vector validity and the existing trans coverage filters.

Excluded bins cannot enter a row pool or acquire a positive final call through
pool expansion or smoothing. BED files that include an explicit unassigned state
may retain these intervals as label `0`. Missing required cis matrices cause the
run to fail rather than silently accepting the rows.

The per-bin `cis_coverage.tsv` report records `cis_nonzero`, `zscore`, the
chromosome's `mean_nonzero` and `sd_nonzero`, and `retained` (1 or 0). SLICE writes
`cis_coverage_RESOLUTIONbp.tsv` for each evaluated grid instead; its final call
mask is also enforced at the output resolution. These reports describe the cis
mask alone; trans coverage, normalization and model-support checks can exclude
additional bins.
