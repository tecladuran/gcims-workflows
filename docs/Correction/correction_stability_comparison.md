Comparison of Raw and Corrected GC-IMS Data
================
Tecla Duran Fort
2025-04-21

## 1. Introduction

This report compares the **raw** and **corrected** GC-IMS Peak Table
using the same evaluation metrics presented in the [Stability
Analysis](https://github.com/tecladuran/gcims-workflows/blob/main/docs/stability_analysis.md),
namely **Relative Standard Deviation (RSD)** and **variance explained by
external factors**.

The correction applied is based on the orthogonalization procedure
described in the [EPO Correction
Report](https://github.com/tecladuran/gcims-workflows/blob/main/docs/epo_correction.md).
Rather than repeating theoretical explanations, this document focuses on
quantifying the improvement in signal stability and reduction of
unwanted variability after correction.

## 2. Apply Correction

``` r
df <- read.csv("../data/peak_table_var.csv")
intensities <- df %>% select(starts_with("Cluster"))

# Apply sequential correction
corr_time <- correction(intensities, df$elapsed_time)
intensities_time_corr <- corr_time$corrected

corr_batch <- correction(intensities_time_corr, df$batch)
intensities_final <- corr_batch$corrected
```

------------------------------------------------------------------------

## 3. Relative Standard Deviation (RSD)

<img src="Comparison_corrected_files/figure-latex/rsd-comparison-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

## 4. Cluster Stability

To evaluate how much the correction improves the technical robustness of
the dataset, we compute the number of clusters whose Relative Standard
Deviation (RSD) falls below the 20% threshold, both **before** and
**after** correction. This threshold is commonly used as an orientative
benchmark in metabolomics quality control.

<div class="figure" style="text-align: center">

<img src="Comparison_corrected_files/figure-latex/rsd-summary-table-1.png" alt="Proportion of clusters with RSD below 20%, before and after correction"  />
<p class="caption">
Proportion of clusters with RSD below 20%, before and after correction
</p>

</div>

The correction process increases the proportion of stable clusters from
**22.6%** to **71%**, confirming that the removal of variance associated
with acquisition order and batch improves overall signal reliability.

------------------------------------------------------------------------

## 5. Explained Variance by External Factors

<img src="Comparison_corrected_files/figure-latex/explained-variance-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

## 6. Principal Component Analysis (PCA)

We now perform a new PCA on the corrected data to explore whether the
dominant sources of variation are still aligned with external variables.
The PCA on raw data already showed strong trends related to
`elapsed_time` and `batch`, as shown in previous reports.

### 6.1 PCA of Raw Data

<div class="figure" style="text-align: center">

<img src="Comparison_corrected_files/figure-latex/pca-raw-plots-1.png" alt="PCA of raw data colored by elapsed time (left) and by batch (right)"  />
<p class="caption">
PCA of raw data colored by elapsed time (left) and by batch (right)
</p>

</div>

<div class="figure" style="text-align: center">

<img src="Comparison_corrected_files/figure-latex/pca-variance-raw-1.png" alt="Explained variance by the first six PCA components (raw data)"  />
<p class="caption">
Explained variance by the first six PCA components (raw data)
</p>

</div>

------------------------------------------------------------------------

### 6.2. PCA of Corrected Data

<div class="figure" style="text-align: center">

<img src="Comparison_corrected_files/figure-latex/pca-corrected-color-1.png" alt="PCA of corrected data colored by elapsed time (left) and by batch (right)"  />
<p class="caption">
PCA of corrected data colored by elapsed time (left) and by batch
(right)
</p>

</div>

<div class="figure" style="text-align: center">

<img src="Comparison_corrected_files/figure-latex/pca-corrected-variance-1.png" alt="Variance explained by the first six PCA components (corrected data)"  />
<p class="caption">
Variance explained by the first six PCA components (corrected data)
</p>

</div>

Compared to the PCA of the raw data, the corrected data shows a more
homogeneous distribution of variance across components, and no evident
separation or gradient is observed when coloring by elapsed time or
batch. This suggests that external influences no longer dominate the
variance structure after correction.

------------------------------------------------------------------------

## 7. Conclusion

The applied correction successfully reduces the influence of time and
batch effects. The PCA projection confirms that corrected data no longer
follows the acquisition order. This confirms the effectiveness of the
orthogonalization strategy when applied directly to the raw intensity
matrix.
