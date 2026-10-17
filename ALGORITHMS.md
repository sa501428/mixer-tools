# MIXER Tools: Technical Algorithm Documentation

**Version:** 5.03.00  
**Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab**

This document provides comprehensive mathematical formulations, algorithmic descriptions, and implementation details for the core algorithms implemented in MIXER Tools: SLICE, ABA, CHIC, and SHUFFLE. This documentation is sufficient for reimplementation and suitable for publication.

---

## Table of Contents

1. [SLICE: Subcompartment Landscape Identification via Clustering Enrichments](#1-slice)
2. [ABA: Aggregate Block Analysis](#2-aba)
3. [CHIC: Calculated Heuristic for Inspecting Clusters](#3-chic)
4. [SHUFFLE: Scoring Heuristics Used For Finding Long-range Enrichments](#4-shuffle)
5. [Common Components and Utilities](#5-common-components)
6. [Implementation Details](#6-implementation-details)
7. [Validation and Quality Control](#7-validation-and-quality-control)

---

## 1. SLICE: Subcompartment Landscape Identification via Clustering Enrichments

### 1.1 Overview and Motivation

SLICE identifies chromatin subcompartments by clustering genomic regions based on their Hi-C contact enrichment patterns. The algorithm addresses the challenge of partitioning the genome into biologically meaningful subcompartments that reflect distinct 3D chromatin organization patterns. SLICE uses k-means or k-medians clustering on a composite feature matrix derived from intra-chromosomal, inter-chromosomal, and genome-wide contact patterns.

**Key Innovation:** SLICE integrates multiple scales of interaction (intra-chromosomal, inter-chromosomal, genome-wide) into a unified feature space, enabling robust subcompartment identification that captures both local and global chromatin organization principles.

### 1.2 Mathematical Formulation

#### 1.2.1 Feature Matrix Construction

For each genomic bin $i$ at resolution $r$ (typically 100 kb), SLICE constructs a feature vector $\mathbf{x}_i \in \mathbb{R}^d$ from three normalized contact matrices:

1. **Intra-chromosomal matrix** ($M_{\text{intra}} \in \mathbb{R}^{n \times n}$): Normalized using SCALE or KR normalization
2. **Inter-chromosomal matrix** ($M_{\text{inter}} \in \mathbb{R}^{n \times m}$): Normalized using INTER_SCALE or INTER_KR normalization  
3. **Genome-wide matrix** ($M_{\text{gw}} \in \mathbb{R}^{n \times p}$): Normalized using GW_SCALE or GW_KR normalization

where $n$ is the number of genomic bins, $m$ is the number of inter-chromosomal bins, and $p$ is the number of genome-wide bins.

**Feature Vector Construction:**

The feature vector for bin $i$ is constructed as:
$$\mathbf{x}_i = [\mathbf{m}_{\text{intra},i}, \mathbf{m}_{\text{inter},i}, \mathbf{m}_{\text{gw},i}]$$

where:
- $\mathbf{m}_{\text{intra},i} \in \mathbb{R}^{n}$ is the $i$-th row (or column, for symmetric matrices) of $M_{\text{intra}}$ representing intra-chromosomal contacts
- $\mathbf{m}_{\text{inter},i} \in \mathbb{R}^{m}$ is the $i$-th row of $M_{\text{inter}}$ representing inter-chromosomal contacts
- $\mathbf{m}_{\text{gw},i} \in \mathbb{R}^{p}$ is the $i$-th row of $M_{\text{gw}}$ representing genome-wide contacts

**Feature Dimension:**
$$d = n + m + p$$

**Composite Feature Matrix:**
$$X = [\mathbf{x}_1^T, \mathbf{x}_2^T, \ldots, \mathbf{x}_n^T]^T \in \mathbb{R}^{n \times d}$$

#### 1.2.2 Matrix Preprocessing and Cleaning

Before clustering, the matrices undergo several preprocessing steps:

**1. Log Transformation with Pseudocount:**
$$M'_{ij} = \log(M_{ij} + 1)$$

This transformation:
- Stabilizes variance across different contact frequency ranges
- Reduces the influence of extreme outliers
- Makes the data more suitable for distance-based clustering

**2. Compression (Optional):**
For compression factor $c > 1$, the matrix is downsampled by aggregating $c \times c$ blocks:
$$M''_{ij} = \sum_{k=0}^{c-1} \sum_{l=0}^{c-1} M'_{ic+k, jc+l}$$

This reduces computational complexity and noise while preserving large-scale patterns.

**3. Sparsity Filtering:**

Bins with insufficient data are identified and removed. For each bin $i$, compute:
- Row sparsity: $\text{sparsity}_{\text{row}}(i) = \frac{\text{count}(\text{NaN in row } i)}{n}$
- Column sparsity: $\text{sparsity}_{\text{col}}(i) = \frac{\text{count}(\text{NaN in col } i)}{n}$

**Bad Index Detection:**

For inter-chromosomal matrices:
- Threshold: $\text{sparsity} > 0.25$ (25% NaN allowed)
- Bins exceeding threshold are marked as "bad indices"

For intra-chromosomal matrices:
- Threshold: $\text{sparsity} > 0.70$ (70% NaN allowed)
- More lenient due to expected sparsity in intra-chromosomal contacts

**Coverage-Based Filtering:**

For each bin $i$, compute coverage statistics:
- Sum of contacts: $s_i = \sum_j M_{ij}$
- Number of non-zero contacts: $nz_i = \sum_j \mathbf{1}_{M_{ij} > 0}$

Compute genome-wide statistics:
- Mean sum: $\bar{s} = \frac{1}{n} \sum_i s_i$
- Std sum: $\sigma_s = \sqrt{\frac{1}{n-1} \sum_i (s_i - \bar{s})^2}$
- Mean non-zero: $\overline{nz} = \frac{1}{n} \sum_i nz_i$
- Std non-zero: $\sigma_{nz} = \sqrt{\frac{1}{n-1} \sum_i (nz_i - \overline{nz})^2}$

**Z-score filtering:**
$$z_{\text{sum}}(i) = \frac{s_i - \bar{s}}{\sigma_s}$$
$$z_{\text{nz}}(i) = \frac{nz_i - \overline{nz}}{\sigma_{nz}}$$

Bins are marked as bad if:
- $z_{\text{sum}}(i) > 5$ (excessive coverage) OR
- $z_{\text{nz}}(i) < -3$ (insufficient non-zero contacts)

**4. Matrix Cleaning:**

After identifying bad indices:
- Set all entries in bad rows/columns to NaN: $M'_{ij} = \text{NaN}$ if $i \in \text{bad}$ or $j \in \text{bad}$
- Remove bad indices from feature matrix: $X' = X[\text{good indices}, :]$

**5. Near-Diagonal Trimming (Intra-chromosomal):**

For intra-chromosomal matrices, trim entries near the diagonal (within pixel distance $d_{\text{trim}} = 3$):
$$M'_{ij} = \text{NaN} \quad \text{if } |i - j| \leq d_{\text{trim}}$$

This removes expected high-frequency contacts that don't contribute to subcompartment signal.

**6. Extreme Value Filtering:**

For observed/expected (OE) matrices, filter extreme values:
$$M'_{ij} = \text{NaN} \quad \text{if } |M_{ij}| > \theta_{\text{OE}} = 30$$

#### 1.2.3 K-Means/K-Medians Clustering

SLICE supports two clustering algorithms:

**K-Means (L2 Distance):**

Objective function:
$$\min_{\mathbf{C}, \boldsymbol{\mu}} \sum_{k=1}^{K} \sum_{i \in C_k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|_2^2$$

where:
- $K$ is the number of clusters (ranging from $K_0$ to $K_F$)
- $C_k$ is the set of bin indices assigned to cluster $k$
- $\boldsymbol{\mu}_k \in \mathbb{R}^d$ is the centroid (mean) of cluster $k$
- $\|\cdot\|_2$ denotes the L2 (Euclidean) norm

**Centroid Update (K-Means):**
$$\boldsymbol{\mu}_k = \frac{1}{|C_k|} \sum_{i \in C_k} \mathbf{x}_i$$

**K-Medians (L1 Distance):**

Objective function:
$$\min_{\mathbf{C}, \boldsymbol{\mu}} \sum_{k=1}^{K} \sum_{i \in C_k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|_1$$

where $\|\cdot\|_1$ denotes the L1 (Manhattan) distance.

**Centroid Update (K-Medians):**
For each dimension $j$:
$$\boldsymbol{\mu}_{kj} = \text{median}\{\mathbf{x}_{ij} : i \in C_k\}$$

K-medians is more robust to outliers and is the default method in SLICE.

**Distance Calculation:**

For k-means:
$$d(\mathbf{x}_i, \boldsymbol{\mu}_k) = \sqrt{\sum_{j=1}^{d} (\mathbf{x}_{ij} - \boldsymbol{\mu}_{kj})^2}$$

For k-medians:
$$d(\mathbf{x}_i, \boldsymbol{\mu}_k) = \sum_{j=1}^{d} |\mathbf{x}_{ij} - \boldsymbol{\mu}_{kj}|$$

**Assignment Step:**
$$C_k = \{i : k = \arg\min_{k'} d(\mathbf{x}_i, \boldsymbol{\mu}_{k'})\}$$

#### 1.2.4 Algorithm Pseudocode

```
ALGORITHM: SLICE Subcompartment Identification

INPUT:
  - Hi-C datasets: D = {D_1, D_2, ..., D_m}
  - Resolution: r (e.g., 100000 bp)
  - Cluster range: [K_0, K_F]
  - Number of runs per K: n_K
  - Maximum iterations: I_max = 200
  - Random seed: s

OUTPUT:
  - Subcompartment assignments for each K ∈ [K_0, K_F]
  - Cluster evaluation metrics

BEGIN
  // Step 1: Initialize
  chromosomeHandler ← getChromosomeHandler(D_1)
  badIndexFinder ← new BadIndexFinder(chromosomes, r, norms)
  
  // Step 2: Identify bad indices
  FOR each dataset D_i in D:
    FOR each chromosome pair (chr1, chr2):
      zd ← getMatrixZoomData(D_i, chr1, chr2, r)
      updateCoverageStats(chr1, chr2, zd)
    END FOR
  END FOR
  determineBadIndicesFromGenomewideStats()
  
  // Step 3: Construct feature matrices
  sliceMatrix ← new SliceMatrix(chromosomeHandler, D_1, norms[0], r, badIndexFinder)
  
  FOR i = 2 to m:
    additionalData ← new SliceMatrix(chromosomeHandler, D_i, norms[i], r, badIndexFinder)
    sliceMatrix.appendDataAlongExistingRows(additionalData)
  END FOR
  
  // Step 4: Clean matrices
  sliceMatrix.cleanUpMatricesBySparsity()
  
  // Step 5: Optional UMAP projection
  IF PROJECT_TO_UMAP:
    projection ← new UmapProjection(sliceMatrix, true)
    projection.runUmapAndColorByChromosome(outputDirectory)
  END IF
  
  // Step 6: Clustering for each K
  FOR K = K_0 to K_F:
    bestResult ← null
    bestScore ← +∞
    
    // Multiple runs for robustness
    FOR run = 1 to n_K:
      seed ← generator.nextLong()
      
      // Initialize centroids randomly
      centroids ← randomInitialize(X, K, seed)
      
      // K-means/K-medians iteration
      FOR iter = 1 to I_max:
        // Assignment step
        assignments ← assignToNearestCentroid(X, centroids, useKMedians)
        
        // Update step
        newCentroids ← updateCentroids(X, assignments, useKMedians)
        
        // Check convergence
        IF converged(centroids, newCentroids):
          BREAK
        END IF
        
        centroids ← newCentroids
      END FOR
      
      // Evaluate clustering
      wcss ← computeWCSS(X, assignments, centroids, useKMedians)
      
      IF wcss < bestScore:
        bestScore ← wcss
        bestResult ← (assignments, centroids, wcss)
      END IF
    END FOR
    
    // Export results
    exportSubcompartments(bestResult, K, prefix)
    updateEvaluator(bestResult, K)
  END FOR
  
  // Step 7: Export evaluation metrics
  exportEvaluatorInfo(evaluator, useKMedians)
END
```

#### 1.2.5 Initialization Strategy

**Random Initialization:**

For each cluster $k$, randomly select a data point as initial centroid:
$$\boldsymbol{\mu}_k^{(0)} = \mathbf{x}_{i_k} \quad \text{where } i_k \sim \text{Uniform}\{1, 2, \ldots, n\}$$

**Multiple Runs:**

The algorithm performs $n_K$ independent runs (default: 3) for each $K$ value, each with a different random seed. The run with the lowest WCSS is selected as the final result.

#### 1.2.6 Convergence Criteria

**Convergence Check:**

The algorithm converges when:
$$\max_k \|\boldsymbol{\mu}_k^{(t)} - \boldsymbol{\mu}_k^{(t-1)}\| < \epsilon$$

where $\epsilon$ is a small threshold (typically $10^{-6}$).

**Maximum Iterations:**

If convergence is not reached within $I_{\max} = 200$ iterations, the algorithm terminates and uses the current clustering.

#### 1.2.7 Cluster Evaluation Metrics

**1. Within-Cluster Sum of Squares (WCSS):**

$$S_{\text{WCSS}} = \frac{1}{K} \sum_{k=1}^{K} \sum_{i \in C_k} d(\mathbf{x}_i, \boldsymbol{\mu}_k)^2$$

where $d$ is the distance metric (L1 for k-medians, L2 for k-means).

**2. Akaike Information Criterion (AIC):**

$$\text{AIC} = S_{\text{WCSS}} + 2 \cdot d \cdot K$$

where $d$ is the feature dimension and $K$ is the number of clusters.

**3. Bayesian Information Criterion (BIC):**

$$\text{BIC} = S_{\text{WCSS}} + \frac{1}{2} \cdot d \cdot K \cdot \log(n)$$

where $n$ is the number of data points.

**4. Silhouette Coefficient (Optional):**

For each data point $i$:
- $a(i)$: Average distance to other points in same cluster
- $b(i)$: Minimum average distance to points in other clusters

$$s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}$$

Overall silhouette:
$$S_{\text{silhouette}} = \frac{1}{n} \sum_{i=1}^{n} s(i)$$

**5. Worst Inter-Cluster Correlation:**

Compute correlation matrix between cluster centroids:
$$R_{kl} = \text{corr}(\boldsymbol{\mu}_k, \boldsymbol{\mu}_l)$$

Worst correlation:
$$R_{\text{worst}} = \max_{k \neq l} R_{kl}$$

Lower values indicate better cluster separation.

**6. Minimum Cluster Size:**

Clusters with fewer than $n_{\min} = 5$ members are penalized in evaluation.

#### 1.2.8 Output Format

For each value of $K \in [K_0, K_F]$, SLICE outputs:

1. **BED file**: `{prefix}_{K}_{kmeans|kmedians}_clusters.bed`
   - Format: `chr start end cluster_id score strand`
   - Each row represents a genomic interval assigned to a cluster

2. **Evaluation metrics**: `{kmeans|kmedians}_cluster_post_stats.npy`
   - Matrix with columns: [K, WCSS, AIC, BIC, Silhouette, WorstCorr]

3. **Optional UMAP projections**: `{prefix}_{K}_{kmeans|kmedians}_clusters_umap.png`

#### 1.2.9 Parameter Selection

**Recommended Parameters:**

- **Resolution**: $r = 100000$ bp (100 kb) for mammalian genomes
- **Cluster range**: $K_0 = 2$, $K_F = 10-15$ depending on genome complexity
- **Number of runs**: $n_K = 3-5$ for robustness
- **Maximum iterations**: $I_{\max} = 200$ (usually converges much earlier)

**Choosing Optimal K:**

1. **Elbow method**: Plot WCSS vs. K, select K at the "elbow"
2. **Information criteria**: Select K that minimizes AIC or BIC
3. **Biological validation**: Choose K that produces biologically interpretable subcompartments

---

## 2. ABA: Aggregate Block Analysis

### 2.1 Overview and Motivation

ABA computes aggregate contact frequency matrices by averaging Hi-C contact patterns across multiple genomic regions (e.g., loop anchors, peaks, binding sites). This enables identification of common interaction patterns shared across functionally similar genomic elements. ABA extracts local contact matrices centered on each region and computes three types of normalized aggregate matrices to highlight different aspects of the interaction pattern.

**Key Innovation:** ABA provides multiple normalization strategies (raw, diagonal-normalized, super-diagonal-normalized) that reveal different aspects of interaction patterns, enabling comprehensive analysis of chromatin architecture at specific genomic features.

### 2.2 Mathematical Formulation

#### 2.2.1 Local Matrix Extraction

For a genomic region $R$ defined by:
- Chromosome: $c_R$
- Start position: $s_R$
- End position: $e_R$
- Center position: $c_R = \frac{s_R + e_R}{2}$

And window size $w$ (default: 20 bins), ABA extracts a local contact matrix $L_R$ of size $L \times L$ where $L = 2w + 1$:

**Bin Coordinates:**

The center bin is:
$$b_{\text{center}} = \left\lfloor \frac{c_R}{r} \right\rfloor$$

where $r$ is the resolution.

**Extraction Window:**

$$b_{\text{start}} = b_{\text{center}} - w$$
$$b_{\text{end}} = b_{\text{center}} + w + 1$$

**Local Matrix:**

$$L_R[i,j] = M[b_{\text{start}} + i, b_{\text{start}} + j]$$

for $i, j \in [0, L-1]$, where $M$ is the normalized Hi-C contact matrix at resolution $r$.

**Edge Case Handling:**

If $b_{\text{start}} < 0$ or $b_{\text{end}} > \text{chromosome length}$:
- Pad with NaN values
- Or adjust window size to fit within chromosome boundaries

#### 2.2.2 Aggregate Matrix Construction

Given $N$ regions $\{R_1, R_2, \ldots, R_N\}$, the raw aggregate matrix is computed as:

$$A_{\text{raw}}[i,j] = \frac{1}{N} \sum_{n=1}^{N} L_{R_n}[i,j]$$

This represents the average contact frequency at offset $(i,j)$ across all regions.

**Accumulation Strategy:**

Instead of computing the mean directly, ABA accumulates sums:
$$S[i,j] = \sum_{n=1}^{N} L_{R_n}[i,j]$$

Then normalizes:
$$A_{\text{raw}}[i,j] = \frac{S[i,j]}{N}$$

This approach:
- Reduces floating-point errors
- Enables incremental computation
- Handles missing data gracefully

#### 2.2.3 Normalization Methods

ABA computes three normalized versions of the aggregate matrix:

**1. Raw Aggregate (ABA):**

$$A_{\text{raw}}[i,j] = \frac{1}{N} \sum_{n=1}^{N} L_{R_n}[i,j]$$

This preserves absolute contact frequencies.

**2. Diagonal-Normalized Aggregate:**

For each local matrix $L_R$, compute diagonal mean:
$$\bar{d}_{\text{diag}}(L_R) = \frac{1}{L} \sum_{k=0}^{L-1} L_R[k,k] \cdot \mathbf{1}_{L_R[k,k] > 0}$$

Normalize:
$$L_R^{\text{diag}}[i,j] = \frac{L_R[i,j]}{\max(1, \bar{d}_{\text{diag}}(L_R))}$$

Aggregate:
$$A_{\text{diag}}[i,j] = \frac{1}{N} \sum_{n=1}^{N} L_{R_n}^{\text{diag}}[i,j]$$

This normalization:
- Removes distance-dependent decay
- Highlights relative interaction patterns
- Useful for identifying loops and domains

**3. Super-Diagonal-Normalized Aggregate:**

For each local matrix $L_R$, compute super-diagonal mean:
$$\bar{d}_{\text{supdiag}}(L_R) = \frac{1}{L-1} \sum_{k=0}^{L-2} L_R[k,k+1] \cdot \mathbf{1}_{L_R[k,k+1] > 0}$$

Normalize:
$$L_R^{\text{supdiag}}[i,j] = \frac{L_R[i,j]}{\max(1, \bar{d}_{\text{supdiag}}(L_R))}$$

Aggregate:
$$A_{\text{supdiag}}[i,j] = \frac{1}{N} \sum_{n=1}^{N} L_{R_n}^{\text{supdiag}}[i,j]$$

This normalization:
- Normalizes by immediate neighbors
- Highlights longer-range interactions
- Useful for identifying specific interaction patterns

#### 2.2.4 Algorithm Pseudocode

```
ALGORITHM: Aggregate Block Analysis (ABA)

INPUT:
  - Hi-C dataset: D
  - BED file: regions = {R_1, R_2, ..., R_N}
  - Resolution: r
  - Window size: w (default: 20)
  - Normalization type: norm

OUTPUT:
  - Raw aggregate matrix: A_raw
  - Diagonal-normalized aggregate: A_diag
  - Super-diagonal-normalized aggregate: A_supdiag

BEGIN
  // Initialize aggregate matrices
  L ← 2 * w + 1
  A_raw ← zeros(L, L)
  A_diag ← zeros(L, L)
  A_supdiag ← zeros(L, L)
  
  N_processed ← 0
  
  // Process each chromosome
  FOR each chromosome chr in chromosomes:
    zd ← getMatrixZoomData(D, chr, chr, r)
    
    // Get regions on this chromosome
    regions_chr ← filter(regions, chromosome == chr)
    
    // Process each region
    FOR each region R in regions_chr:
      // Extract local matrix
      center_bin ← floor((R.start + R.end) / 2 / r)
      bin_start ← center_bin - w
      bin_end ← center_bin + w + 1
      
      // Handle edge cases
      IF bin_start < 0 OR bin_end > chromosome_length:
        CONTINUE  // Skip regions near boundaries
      END IF
      
      // Extract local matrix
      L_R ← extractLocalBoundedRegion(zd, bin_start, bin_end, bin_start, bin_end, L, L, norm)
      
      // Clean NaN values
      cleanUpNaNs(L_R)
      
      // Add to raw aggregate
      A_raw ← A_raw + L_R
      
      // Compute normalizations
      d_diag ← diagonalMean(L_R)
      d_supdiag ← supDiagonalMean(L_R)
      
      // Normalize and add
      L_R_diag ← L_R / max(1, d_diag)
      L_R_supdiag ← L_R / max(1, d_supdiag)
      
      A_diag ← A_diag + L_R_diag
      A_supdiag ← A_supdiag + L_R_supdiag
      
      N_processed ← N_processed + 1
    END FOR
  END FOR
  
  // Final normalization
  A_raw ← A_raw / N_processed
  A_diag ← A_diag / N_processed
  A_supdiag ← A_supdiag / N_processed
  
  // Export matrices
  saveMatrix(A_raw, "ABA.npy")
  saveMatrix(A_diag, "diagNormedABA.npy")
  saveMatrix(A_supdiag, "supDiagNormedABA.npy")
END
```

#### 2.2.5 Helper Functions

**Diagonal Mean:**
```
FUNCTION diagonalMean(matrix):
  sum ← 0
  count ← 0
  FOR i = 0 to L-1:
    IF matrix[i,i] > 0:
      sum ← sum + matrix[i,i]
      count ← count + 1
    END IF
  END FOR
  RETURN sum / max(1, count)
END FUNCTION
```

**Super-Diagonal Mean:**
```
FUNCTION supDiagonalMean(matrix):
  sum ← 0
  count ← 0
  FOR i = 0 to L-2:
    IF matrix[i,i+1] > 0:
      sum ← sum + matrix[i,i+1]
      count ← count + 1
    END IF
  END FOR
  RETURN sum / max(1, count)
END FUNCTION
```

#### 2.2.6 Genome-Wide Aggregation

For genome-wide analysis across multiple chromosomes, ABA maintains separate aggregate matrices per chromosome and then combines them:

$$A_{\text{gw}}[i,j] = \frac{1}{N_{\text{total}}} \sum_{\text{all chromosomes}} \sum_{\text{regions in chr}} L_R[i,j]$$

where $N_{\text{total}}$ is the total number of regions across all chromosomes.

#### 2.2.7 Parallelization

ABA processes chromosomes in parallel using multiple threads:

```
PARALLEL PROCESSING:
  Thread pool ← createThreadPool(num_threads)
  
  FOR each chromosome chr:
    submitTask(() -> {
      processChromosome(chr, regions, D, r, w, norm)
    })
  END FOR
  
  waitForAllTasks()
  
  combineResults()
END
```

#### 2.2.8 Output Format

ABA outputs three NumPy arrays for each resolution:

1. **ABA.npy**: Raw aggregate contact matrix (L × L)
2. **supDiagNormedABA.npy**: Super-diagonal-normalized aggregate matrix (L × L)
3. **diagNormedABA.npy**: Diagonal-normalized aggregate matrix (L × L)

**Matrix Interpretation:**

- **Center** $(w, w)$: The anchor point (region center)
- **Diagonal**: Contacts at genomic distance $d = |i-j| \cdot r$
- **Off-diagonal**: Contacts between different genomic distances

---

## 3. CHIC: Calculated Heuristic for Inspecting Clusters

### 3.1 Overview and Motivation

CHIC (Calculated Heuristic for Inspecting Clusters) is a quantitative metric that evaluates how well subcompartment clusters are preserved in Hi-C contact maps. It compares baseline (observed) contact patterns with shuffled (null model) patterns using two complementary metrics: variance and relative entropy (KL divergence). CHIC provides a rigorous statistical framework for assessing cluster quality and biological significance.

**Key Innovation:** CHIC uses cluster-preserving shuffling to create a null model that maintains cluster boundaries while destroying spatial relationships, enabling quantification of cluster preservation strength.

### 3.2 Mathematical Formulation

#### 3.2.1 Cluster Index Mapping

Given subcompartment clusters defined by BED intervals, CHIC maps genomic bins to cluster IDs.

**Cluster Definition:**

For cluster $c$ with genomic intervals $\{I_1, I_2, \ldots, I_{n_c}\}$, the set of bin indices is:

$$C_{\text{row}}[c] = \left\{i : \exists I_k \text{ s.t. } \left\lfloor \frac{\text{start}(I_k)}{r} \right\rfloor \leq i \leq \left\lfloor \frac{\text{end}(I_k)}{r} \right\rfloor \right\}$$

Similarly for columns:
$$C_{\text{col}}[c] = \left\{j : \exists I_k \text{ s.t. } \left\lfloor \frac{\text{start}(I_k)}{r} \right\rfloor \leq j \leq \left\lfloor \frac{\text{end}(I_k)}{r} \right\rfloor \right\}$$

**Cluster Boundaries:**

For cluster $c$, define boundaries:
- Row boundaries: $b_{\text{row}}[c] = [\min(C_{\text{row}}[c]), \max(C_{\text{row}}[c])]$
- Column boundaries: $b_{\text{col}}[c] = [\min(C_{\text{col}}[c]), \max(C_{\text{col}}[c])]$

#### 3.2.2 Aggregate Matrix Construction

For a contact matrix $M \in \mathbb{R}^{n \times m}$ and cluster assignments, CHIC constructs an aggregate matrix $A \in \mathbb{R}^{K \times K}$ where each element $A[rI, cI]$ represents average contacts between cluster $rI$ and cluster $cI$:

$$A[rI, cI] = \frac{1}{|C_{\text{row}}[rI]| \cdot |C_{\text{col}}[cI]|} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} M[i,j]$$

where $K$ is the number of clusters.

#### 3.2.3 Variance Score

The variance score measures the variability of contact frequencies within cluster regions.

**Mean Contact Frequency:**

For cluster pair $(rI, cI)$:
$$\mu_{rI,cI} = \frac{1}{N_{rI,cI}} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} M[i,j]$$

where $N_{rI,cI} = |C_{\text{row}}[rI]| \cdot |C_{\text{col}}[cI]|$ is the number of elements in the cluster pair region.

**Variance Score:**

$$S_{\text{var}} = \frac{1}{N_{\text{total}}} \sum_{rI=1}^{K} \sum_{cI=1}^{K} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} (M[i,j] - \mu_{rI,cI})^2$$

where $N_{\text{total}} = \sum_{rI,cI} N_{rI,cI}$ is the total number of matrix elements.

**Interpretation:**

- **Low variance**: Contacts are uniform within cluster pairs (strong cluster structure)
- **High variance**: Contacts vary widely within cluster pairs (weak cluster structure)

#### 3.2.4 Relative Entropy (KL Divergence) Score

The KL divergence score measures the information divergence between observed and expected contact distributions.

**Observed Probability Distribution:**

$$p_{i,j} = \frac{M[i,j]}{\sum_{i'=1}^{n} \sum_{j'=1}^{m} M[i',j']}$$

This represents the empirical distribution of contacts.

**Expected Probability Distribution (Cluster-Based):**

For cluster pair $(rI, cI)$, the expected probability is computed as the mean contact frequency for that cluster pair, normalized by the total sum:
$$q_{rI,cI} = \frac{\mu_{rI,cI}}{S_{\text{total}}} = \frac{\frac{1}{N_{rI,cI}} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} M[i,j]}{\sum_{i'=1}^{n} \sum_{j'=1}^{m} M[i',j']}$$

where $S_{\text{total}} = \sum_{i'=1}^{n} \sum_{j'=1}^{m} M[i',j']$ is the total sum of all matrix elements.

**Note:** In the implementation, $S_{\text{total}}$ is computed as the sum of all cluster pair sums (since each matrix element belongs to exactly one cluster pair), which is mathematically equivalent.

This represents the expected distribution if contacts were uniform within cluster pairs.

**KL Divergence (Forward):**

$$S_{\text{KL}} = \sum_{rI=1}^{K} \sum_{cI=1}^{K} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} p_{i,j} \log\left(\frac{p_{i,j}}{q_{rI,cI}}\right)$$

**KL Divergence (Reverse):**

$$S_{\text{KL,rev}} = \sum_{rI=1}^{K} \sum_{cI=1}^{K} \sum_{i \in C_{\text{row}}[rI]} \sum_{j \in C_{\text{col}}[cI]} q_{rI,cI} \log\left(\frac{q_{rI,cI}}{p_{i,j}}\right)$$

**Properties:**

- $S_{\text{KL}} \geq 0$ (Gibbs' inequality)
- $S_{\text{KL}} = 0$ if and only if $p_{i,j} = q_{rI,cI}$ for all $(i,j)$ in cluster pair $(rI,cI)$
- Higher values indicate greater divergence from cluster-uniform expectation

#### 3.2.5 Shuffled Matrix Generation

To create a null model, CHIC generates shuffled matrices by randomly permuting bin indices within each cluster.

**Cluster-Preserving Shuffle:**

For each cluster $c$ with bins $\{b_1, b_2, \ldots, b_{n_c}\}$:

1. Generate random permutation: $\pi_c: \{1, 2, \ldots, n_c\} \rightarrow \{1, 2, \ldots, n_c\}$
2. Shuffled bin mapping: $b'_i = b_{\pi_c(i)}$

**Shuffled Matrix Construction:**

For shuffled row indices $\mathbf{i}_{\text{shuffled}}$ and column indices $\mathbf{j}_{\text{shuffled}}$:

$$M_{\text{shuffled}}[i', j'] = M[\mathbf{i}_{\text{shuffled}}[i'], \mathbf{j}_{\text{shuffled}}[j']]$$

**Key Property:**

The shuffle preserves:
- Cluster boundaries (bins stay within their clusters)
- Cluster sizes (number of bins per cluster unchanged)

The shuffle destroys:
- Spatial relationships (genomic order is randomized)
- Distance-dependent contact patterns

#### 3.2.6 Multiple Shuffle Rounds

CHIC performs $R = 50$ independent shuffle rounds to obtain robust statistics.

**Shuffle Round $r$:**

1. Generate shuffled indices: $\mathbf{i}^{(r)}_{\text{shuffled}}, \mathbf{j}^{(r)}_{\text{shuffled}}$
2. Construct shuffled matrix: $M^{(r)}_{\text{shuffled}}$
3. Compute scores: $S_{\text{var}}^{(r)}, S_{\text{KL}}^{(r)}$

**Averaged Scores:**

$$\bar{S}_{\text{var,shuffled}} = \frac{1}{R} \sum_{r=1}^{R} S_{\text{var}}^{(r)}$$

$$\bar{S}_{\text{KL,shuffled}} = \frac{1}{R} \sum_{r=1}^{R} S_{\text{KL}}^{(r)}$$

#### 3.2.7 CHIC Score Calculation

The CHIC score is the ratio of baseline to shuffled scores:

**Variance-based CHIC Score:**

$$\text{CHIC}_{\text{var}} = \frac{S_{\text{var}}}{\bar{S}_{\text{var,shuffled}}}$$

**KL Divergence-based CHIC Score:**

$$\text{CHIC}_{\text{KL}} = \frac{S_{\text{KL}}}{\bar{S}_{\text{KL,shuffled}}}$$

**Interpretation:**

- **CHIC > 1**: Baseline shows stronger cluster structure than shuffled (good cluster preservation)
- **CHIC ≈ 1**: Baseline and shuffled are similar (weak cluster structure)
- **CHIC < 1**: Shuffled shows stronger structure (rare, may indicate data quality issues)

**Statistical Significance:**

Higher CHIC scores indicate:
- Stronger preservation of cluster boundaries in contact patterns
- More biologically meaningful subcompartment assignments
- Better separation between clusters

#### 3.2.8 Geometric Mean

For multiple chromosomes or map types, CHIC computes the geometric mean:

$$\text{CHIC}_{\text{GM}} = \left(\prod_{m=1}^{M} \text{CHIC}_m\right)^{1/M}$$

where $M$ is the number of chromosomes or map types analyzed.

**Properties:**

- Geometric mean is less sensitive to outliers than arithmetic mean
- Appropriate for ratio-based metrics
- Preserves multiplicative relationships

#### 3.2.9 Algorithm Pseudocode

```
ALGORITHM: CHIC Score Calculation

INPUT:
  - Contact matrix: M (n × m)
  - Subcompartment BED file: clusters
  - Resolution: r
  - Number of shuffle rounds: R = 50

OUTPUT:
  - CHIC variance score
  - CHIC KL divergence score

BEGIN
  // Step 1: Map bins to clusters
  clusterToRowIndices ← populateCluster(chromosomes, offsets, clusters, r)
  clusterToColIndices ← populateCluster(chromosomes, offsets, clusters, r)
  
  // Step 2: Compute baseline scores
  baseline_var ← computeVarianceScore(M, clusterToRowIndices, clusterToColIndices)
  baseline_kl ← computeKLDivergenceScore(M, clusterToRowIndices, clusterToColIndices)
  
  // Step 3: Generate shuffled matrices and compute scores
  shuffled_var_sum ← 0
  shuffled_kl_sum ← 0
  
  FOR r = 1 to R:
    // Generate shuffled indices (preserving cluster boundaries)
    shuffledRowIndices ← shuffleByCluster(clusterToRowIndices, random_seed[r])
    shuffledColIndices ← shuffleByCluster(clusterToColIndices, random_seed[r])
    
    // Construct shuffled matrix
    M_shuffled ← constructShuffledMatrix(M, shuffledRowIndices, shuffledColIndices)
    
    // Compute scores
    var_score ← computeVarianceScore(M_shuffled, clusterToRowIndices, clusterToColIndices)
    kl_score ← computeKLDivergenceScore(M_shuffled, clusterToRowIndices, clusterToColIndices)
    
    shuffled_var_sum ← shuffled_var_sum + var_score
    shuffled_kl_sum ← shuffled_kl_sum + kl_score
  END FOR
  
  // Step 4: Average shuffled scores
  avg_shuffled_var ← shuffled_var_sum / R
  avg_shuffled_kl ← shuffled_kl_sum / R
  
  // Step 5: Compute CHIC scores
  CHIC_var ← baseline_var / avg_shuffled_var
  CHIC_kl ← baseline_kl / avg_shuffled_kl
  
  RETURN (CHIC_var, CHIC_kl)
END
```

#### 3.2.10 Helper Functions

**Variance Score Computation:**
```
FUNCTION computeVarianceScore(M, rowClusters, colClusters):
  // Compute means for each cluster pair
  meanMap ← new Map()
  numRegionMap ← new Map()
  
  FOR each cluster pair (rI, cI):
    sum ← 0
    count ← 0
    
    FOR i in rowClusters[rI]:
      FOR j in colClusters[cI]:
        sum ← sum + M[i,j]
        count ← count + 1
      END FOR
    END FOR
    
    meanMap[(rI,cI)] ← sum / count
    numRegionMap[(rI,cI)] ← count
  END FOR
  
  // Compute variance
  total_variance ← 0
  total_elements ← 0
  
  FOR each cluster pair (rI, cI):
    mu ← meanMap[(rI,cI)]
    
    FOR i in rowClusters[rI]:
      FOR j in colClusters[cI]:
        total_variance ← total_variance + (M[i,j] - mu)^2
        total_elements ← total_elements + 1
      END FOR
    END FOR
  END FOR
  
  RETURN total_variance / total_elements
END FUNCTION
```

**KL Divergence Score Computation:**
```
FUNCTION computeKLDivergenceScore(M, rowClusters, colClusters):
  // Compute total sum
  total_sum ← sum(M)
  
  // Compute observed probabilities
  p[i,j] ← M[i,j] / total_sum
  
  // Compute expected probabilities (cluster-based)
  qMap ← new Map()
  
  FOR each cluster pair (rI, cI):
    cluster_sum ← 0
    cluster_count ← 0
    
    FOR i in rowClusters[rI]:
      FOR j in colClusters[cI]:
        cluster_sum ← cluster_sum + M[i,j]
        cluster_count ← cluster_count + 1
      END FOR
    END FOR
    
    cluster_mean ← cluster_sum / cluster_count
    qMap[(rI,cI)] ← cluster_mean / total_sum
  END FOR
  
  // Compute KL divergence
  kl_divergence ← 0
  
  FOR each cluster pair (rI, cI):
    q ← qMap[(rI,cI)]
    
    FOR i in rowClusters[rI]:
      FOR j in colClusters[cI]:
        p_val ← p[i,j]
        IF p_val > 0 AND q > 0:
          kl_divergence ← kl_divergence + p_val * log(p_val / q)
        END IF
      END FOR
    END FOR
  END FOR
  
  RETURN kl_divergence
END FUNCTION
```

---

## 4. SHUFFLE: Scoring Heuristics Used For Finding Long-range Enrichments

### 4.1 Overview and Motivation

SHUFFLE evaluates long-range enrichments by comparing observed Hi-C contact patterns with shuffled null models. It generates aggregate contact matrices from multiple shuffled iterations and computes statistical scores (CHIC scores) to assess the significance of cluster-preserved interactions. SHUFFLE provides a rigorous framework for identifying biologically significant long-range chromatin interactions.

**Key Innovation:** SHUFFLE uses cluster-preserving shuffling combined with compression and aggregation to create robust null models that enable statistical evaluation of interaction enrichments while maintaining computational efficiency.

### 4.2 Mathematical Formulation

#### 4.2.1 Shuffling Strategy

SHUFFLE employs cluster-preserving shuffling: bins are randomly permuted within their assigned clusters, maintaining cluster boundaries while destroying inter-cluster spatial relationships.

**Cluster-Preserved Shuffle:**

For cluster $c$ with bins $\{b_1, b_2, \ldots, b_{n_c}\}$:

1. Generate random permutation: $\pi_c: \{1, 2, \ldots, n_c\} \rightarrow \{1, 2, \ldots, n_c\}$
2. Shuffled bin mapping: $b'_i = b_{\pi_c(i)}$

**Compression Factor:**

To reduce computational complexity, bins are grouped into compression blocks of size $c$ (default: $c = 8$). Only bins divisible by $c$ are used:

$$n_{\text{used}} = \left\lfloor \frac{n_c}{c} \right\rfloor \cdot c$$

**Shuffled Index Construction:**

```
FUNCTION getShuffledByClusterIndices(clusterToIndices, generator):
  allIndices ← []
  boundaries ← [0]
  clusterIDs ← []
  
  FOR each cluster c in sorted(clusterToIndices.keys()):
    indexList ← clusterToIndices[c]
    shuffle(indexList, generator)  // Random permutation
    
    // Use only bins divisible by compression factor
    numToUse ← (length(indexList) / compressionFactor) * compressionFactor
    
    FOR i = 0 to numToUse-1:
      allIndices.append(indexList[i])
    END FOR
    
    count ← count + (numToUse / compressionFactor)
    boundaries.append(count)
    clusterIDs.append(c)
  END FOR
  
  RETURN ShuffledIndices(allIndices, boundaries, clusterIDs)
END FUNCTION

**Note:** The boundaries are stored in compressed bin units (divided by compression factor), not original bin units. This is because the boundaries are used to index into the compressed matrix.
```

#### 4.2.2 Shuffled Matrix Construction

For a contact matrix $M$ and shuffled bin indices $\mathbf{i}_{\text{shuffled}}$ and $\mathbf{j}_{\text{shuffled}}$, the shuffled matrix is:

$$M_{\text{shuffled}}[i', j'] = M[\mathbf{i}_{\text{shuffled}}[i'], \mathbf{j}_{\text{shuffled}}[j']]$$

#### 4.2.3 Compression

To reduce computational complexity, SHUFFLE applies compression factor $c$:

$$M_{\text{compressed}}[i, j] = \sum_{k=0}^{c-1} \sum_{l=0}^{c-1} M_{\text{shuffled}}[ic + k, jc + l]$$

This aggregates $c \times c$ blocks into single elements, reducing matrix size by factor $c^2$.

**Rationale:**

- Reduces memory requirements
- Speeds up computation
- Preserves large-scale patterns
- Reduces noise

#### 4.2.4 Log Transformation

After compression, apply log transformation with pseudocount:

$$M_{\text{log}}[i,j] = \log(M_{\text{compressed}}[i,j] + 1)$$

This:
- Stabilizes variance
- Reduces influence of extreme values
- Makes data more suitable for statistical analysis

#### 4.2.5 Aggregate Matrix Generation

SHUFFLE performs $R = 50$ independent shuffling rounds. For each round $r$:

1. Generate shuffled indices: $\mathbf{i}^{(r)}_{\text{shuffled}}, \mathbf{j}^{(r)}_{\text{shuffled}}$
2. Construct shuffled matrix: $M^{(r)}_{\text{shuffled}}$
3. Apply compression: $M^{(r)}_{\text{compressed}}$
4. Apply log transformation: $M^{(r)}_{\text{log}}[i,j] = \log(M^{(r)}_{\text{compressed}}[i,j] + 1)$

**Aggregate Matrix:**

$$A_{\text{shuffled}} = \frac{1}{R} \sum_{r=1}^{R} M^{(r)}_{\text{log}}$$

This represents the expected contact pattern under the null model (cluster-preserving shuffle).

#### 4.2.6 Baseline Matrix

**Important Note:** In the SHUFFLE implementation, both baseline and shuffled scores are computed on the same aggregate shuffled matrix. The distinction is:

- **Baseline score**: Computed on the aggregate matrix treating the entire matrix as one region (no cluster boundaries)
- **Shuffled score**: Computed on the aggregate matrix using cluster boundaries

This approach evaluates how cluster boundaries affect the scoring metrics on the shuffled (null model) data.

The baseline (observed) matrix processing:

1. Apply compression: $M_{\text{compressed}}[i,j] = \sum_{k=0}^{c-1} \sum_{l=0}^{c-1} M[ic+k, jc+l]$
2. Apply log transformation: $A_{\text{baseline}} = \log(M_{\text{compressed}} + 1)$

However, in practice, SHUFFLE computes scores on the aggregate shuffled matrix for both baseline and shuffled cases, using different boundary definitions.

#### 4.2.7 Inter-Chromosomal Analysis

For inter-chromosomal analysis, SHUFFLE uses three mapping strategies:

**1. Odds vs Evens:**

Compare odd-numbered chromosomes with even-numbered chromosomes:
- Row chromosomes: $\{1, 3, 5, \ldots\}$
- Column chromosomes: $\{2, 4, 6, \ldots\}$

**2. Skip by Twos:**

Compare chromosomes with indices $i$ and $i+2$:
- Row chromosomes: $\{1, 2, 3, \ldots\}$
- Column chromosomes: $\{3, 4, 5, \ldots\}$

**3. First Half vs Second Half:**

Compare first half of chromosomes with second half:
- Row chromosomes: $\{1, 2, \ldots, \lfloor n/2 \rfloor\}$
- Column chromosomes: $\{\lfloor n/2 \rfloor + 1, \ldots, n\}$

Each strategy generates a separate inter-chromosomal contact matrix for analysis.

#### 4.2.8 Scoring

SHUFFLE uses the same scoring methods as CHIC (variance and KL divergence) to compare baseline and shuffled aggregate matrices. The scores are computed using the `ScoreContainer` class.

**Baseline Scores:**

Compute on aggregate shuffled matrix, treating entire matrix as one region (no cluster boundaries):
- $S_{\text{var,baseline}} = \text{variance}(A_{\text{shuffled}}, \text{no boundaries})$
- $S_{\text{KL,baseline}} = \text{KLDivergence}(A_{\text{shuffled}}, \text{no boundaries})$

**Shuffled Scores:**

Compute on aggregate shuffled matrix, using cluster boundaries:
- $S_{\text{var,shuffled}} = \text{variance}(A_{\text{shuffled}}, \text{with cluster boundaries})$
- $S_{\text{KL,shuffled}} = \text{KLDivergence}(A_{\text{shuffled}}, \text{with cluster boundaries})$

**Note:** Both scores are computed on the same aggregate shuffled matrix. The difference is whether cluster boundaries are used in the scoring calculation. This evaluates how cluster structure affects the metrics on null model data.

**CHIC Scores:**

$$\text{CHIC}_{\text{var}} = \frac{S_{\text{var,baseline}}}{S_{\text{var,shuffled}}}$$

$$\text{CHIC}_{\text{KL}} = \frac{S_{\text{KL,baseline}}}{S_{\text{KL,shuffled}}}$$

#### 4.2.9 Algorithm Pseudocode

```
ALGORITHM: SHUFFLE Analysis

INPUT:
  - Hi-C dataset: D
  - Subcompartment BED file: clusters
  - Resolution: r
  - Compression factor: c = 8
  - Number of shuffle rounds: R = 50
  - Analysis type: {intra, inter}

OUTPUT:
  - Aggregate matrices (baseline and shuffled)
  - CHIC scores

BEGIN
  // Step 1: Map bins to clusters
  clusterToRowIndices ← populateCluster(chromosomes, offsets, clusters, r)
  clusterToColIndices ← populateCluster(chromosomes, offsets, clusters, r)
  
  // Step 2: Construct original matrix (used for shuffling)
  IF analysis_type == "intra":
    M_original ← constructIntraMatrix(D, r, norm, intraType)
  ELSE:
    M_original ← constructInterMatrix(D, r, norm, mapType)
  END IF
  
  // Step 3: Generate shuffled matrices
  aggregate_shuffled ← new AggregateMatrix()
  
  seeds ← generateSeeds(R, random_generator)
  
  PARALLEL FOR r = 1 to R:
    // Generate shuffled indices
    shuffledRowIndices ← getShuffledByClusterIndices(clusterToRowIndices, seeds[r])
    shuffledColIndices ← getShuffledByClusterIndices(clusterToColIndices, seeds[r])
    
    // Construct shuffled matrix
    M_shuffled ← constructShuffledMatrix(M_original, shuffledRowIndices, shuffledColIndices)
    
    // Compress and log transform
    M_shuffled_compressed ← compress(M_shuffled, c)
    M_shuffled_log ← log(M_shuffled_compressed + 1)
    
    // Add to aggregate
    aggregate_shuffled.add(M_shuffled_log)
  END PARALLEL FOR
  
  // Step 5: Average shuffled matrices
  aggregate_shuffled.scaleForNumberOfRounds(R)
  A_shuffled ← aggregate_shuffled.getMatrix()
  
  // Step 6: Compute scores
  // Baseline scores: compute on aggregate matrix without cluster boundaries
  baseline_scores ← computeScores(A_shuffled, no_boundaries)
  // Shuffled scores: compute on aggregate matrix with cluster boundaries  
  shuffled_scores ← computeScores(A_shuffled, clusterToRowIndices, clusterToColIndices)
  
  // Step 7: Compute CHIC scores
  CHIC_var ← baseline_scores.var / shuffled_scores.var
  CHIC_kl ← baseline_scores.kl / shuffled_scores.kl
  
  // Step 8: Save results
  saveMatrixPNG(A_baseline, "baseline")
  saveMatrixPNG(A_shuffled, "shuffled")
  saveScores(CHIC_var, CHIC_kl, "aggregate_scores.txt")
END
```

#### 4.2.10 Parallelization

SHUFFLE processes shuffle rounds in parallel:

```
PARALLEL PROCESSING:
  Thread pool ← createThreadPool(num_threads)
  task_queue ← queue(range(1, R+1))
  
  WHILE task_queue not empty:
    r ← task_queue.pop()
    submitTask(() -> {
      processShuffleRound(r, seeds[r], clusterToRowIndices, clusterToColIndices)
    })
  END WHILE
  
  waitForAllTasks()
  
  combineResults()
END
```

#### 4.2.11 Output Format

SHUFFLE generates:

1. **Aggregate matrices** (PNG images):
   - `baseline.png`: Baseline (observed) aggregate matrix
   - `shuffled.png`: Shuffled (null model) aggregate matrix

2. **CHIC scores** (text file): `aggregate_scores_{prefix}.txt`
   - Variance-based CHIC scores
   - KL divergence-based CHIC scores
   - Geometric means across chromosomes/map types

---

## 5. Common Components and Utilities

### 5.1 Normalization Methods

All algorithms support multiple normalization types:

**1. NONE:**
No normalization applied.

**2. VC (Vanilla Coverage):**
Row and column normalization:
$$M'_{ij} = \frac{M_{ij}}{\sqrt{\sum_k M_{ik} \cdot \sum_k M_{kj}}}$$

**3. VC_SQRT (Square Root Vanilla Coverage):**
$$M'_{ij} = \frac{M_{ij}}{\sqrt[4]{\sum_k M_{ik} \cdot \sum_k M_{kj}}}$$

**4. KR (Knight-Ruiz Matrix Balancing):**
Iterative balancing algorithm that finds vectors $\mathbf{u}, \mathbf{v}$ such that:
$$M'_{ij} = M_{ij} \cdot u_i \cdot v_j$$

where row and column sums are balanced.

**5. GW_KR (Genome-Wide KR):**
KR normalization applied genome-wide across all chromosomes.

**6. INTER_KR (Inter-Chromosomal KR):**
KR normalization applied to inter-chromosomal contacts only.

### 5.2 Matrix Operations

**Log Transformation:**
$$f(x) = \log(x + 1)$$

**NaN Handling:**
- Zero out: Replace NaN with 0
- Impute: Replace NaN with mean/median of row/column
- Remove: Filter out rows/columns with excessive NaN

**Symmetry Enforcement:**
For intra-chromosomal matrices:
$$M'_{ij} = M'_{ji} = \frac{M_{ij} + M_{ji}}{2}$$

### 5.3 Data Structures

**ShuffledIndices:**
```java
class ShuffledIndices {
  List<Integer> allIndices;      // Shuffled bin indices
  Integer[] boundaries;          // Cluster boundaries
  Integer[] ids;                 // Cluster IDs
}
```

**ScoreContainer:**
```java
class ScoreContainer {
  float[][] baselines;            // Baseline scores [numMaps][numScores]
  float[][] shuffled;             // Shuffled scores [numMaps][numScores]
  float[][] ratios;               // CHIC scores [numMaps][numScores]
}
```

---

## 6. Implementation Details

### 6.1 Computational Complexity

**SLICE:**
- Feature matrix construction: $O(n \cdot d)$
- K-means iteration: $O(n \cdot d \cdot K \cdot I)$
- Overall: $O(n \cdot d \cdot K \cdot I \cdot n_K)$
  - $n$: number of bins
  - $d$: feature dimension
  - $K$: number of clusters
  - $I$: iterations
  - $n_K$: number of runs per K

**ABA:**
- Per region: $O(w^2)$
- Overall: $O(N \cdot w^2)$
  - $N$: number of regions
  - $w$: window size

**CHIC:**
- Baseline score: $O(n^2)$
- Shuffled scores: $O(R \cdot n^2)$
- Overall: $O(R \cdot n^2)$
  - $R$: number of shuffle rounds
  - $n$: matrix size

**SHUFFLE:**
- Per shuffle round: $O(n^2 / c^2)$ (with compression)
- Overall: $O(R \cdot n^2 / c^2)$
  - $R$: number of shuffle rounds
  - $n$: matrix size
  - $c$: compression factor

### 6.2 Memory Requirements

**SLICE:**
- Feature matrix: $O(n \cdot d)$
- Cluster assignments: $O(n)$
- Total: $O(n \cdot d)$

**ABA:**
- Aggregate matrices: $O(w^2)$
- Per-region matrix: $O(w^2)$
- Total: $O(w^2)$

**CHIC/SHUFFLE:**
- Contact matrix: $O(n^2)$
- Shuffled matrices: $O(R \cdot n^2 / c^2)$ (with compression)
- Total: $O(n^2 + R \cdot n^2 / c^2)$

### 6.3 Parallelization Strategies

**SLICE:**
- Parallel processing across multiple K values
- Parallel processing across multiple runs
- Parallel distance calculations

**ABA:**
- Parallel processing across chromosomes
- Thread-safe accumulation of aggregate matrices

**SHUFFLE:**
- Parallel processing across shuffle rounds
- Thread-safe aggregation

### 6.4 Numerical Stability

**Floating-Point Precision:**
- Use double precision for accumulations
- Convert to float for storage
- Handle underflow/overflow gracefully

**Division by Zero:**
- Use `max(1, denominator)` to avoid division by zero
- Handle NaN propagation carefully

**Logarithm:**
- Use `log(x + 1)` to avoid log(0)
- Check for negative values before log transformation

---

## 7. Validation and Quality Control

### 7.1 Input Validation

**Hi-C File Validation:**
- Check file format (.hic)
- Verify resolution availability
- Check normalization availability
- Validate chromosome names

**BED File Validation:**
- Check format (chr, start, end)
- Validate chromosome names
- Check coordinate ranges
- Handle overlapping intervals

### 7.2 Output Validation

**SLICE:**
- Verify all bins are assigned to clusters
- Check cluster sizes (minimum: 5 bins)
- Validate BED file format
- Check evaluation metrics are finite

**ABA:**
- Verify matrix dimensions ($L \times L$)
- Check for NaN values
- Validate normalization factors
- Verify region counts

**CHIC/SHUFFLE:**
- Verify CHIC scores are positive
- Check for convergence (scores stabilize)
- Validate aggregate matrices
- Verify score consistency across runs

### 7.3 Quality Metrics

**SLICE:**
- WCSS should decrease with increasing K (up to optimal K)
- Silhouette coefficient: $s \in [-1, 1]$, higher is better
- Worst correlation: lower is better (better separation)

**CHIC:**
- CHIC scores should be $> 1$ for meaningful clusters
- Variance between runs should be small
- Geometric mean should be consistent across chromosomes

### 7.4 Error Handling

**Missing Data:**
- Skip regions with insufficient data
- Log warnings for skipped regions
- Continue processing remaining regions

**Convergence Issues:**
- Maximum iteration limit prevents infinite loops
- Log convergence status
- Use best result even if not fully converged

**Memory Issues:**
- Compression reduces memory requirements
- Process chromosomes sequentially if needed
- Clear intermediate results

---

## 8. References and Source Code

This documentation is based on MIXER Tools source code (version 5.03.00). Key source files:

**SLICE:**
- `src/mixer/algos/Slice.java`
- `src/mixer/utils/slice/kmeans/FullGenomeOEWithinClusters.java`
- `src/mixer/utils/slice/kmeans/GenomeWideKmeansRunner.java`
- `src/mixer/utils/slice/kmeans/KmeansEvaluator.java`
- `src/mixer/utils/slice/cleaning/BadIndexFinder.java`

**ABA:**
- `src/mixer/algos/ABA.java`
- `src/mixer/utils/aba/ABADataStack.java`
- `src/mixer/utils/bed/BedFile.java`

**CHIC:**
- `src/mixer/utils/shuffle/CHICTools.java`
- `src/mixer/utils/shuffle/scoring/VarianceScoring.java`
- `src/mixer/utils/shuffle/scoring/KLDivergenceScoring.java`
- `src/mixer/utils/matrix/ScoreContainer.java`

**SHUFFLE:**
- `src/mixer/utils/shuffle/ShuffleAction.java`
- `src/mixer/utils/matrix/AggregateMatrix.java`

---

**Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab**

**License:** MIT License
