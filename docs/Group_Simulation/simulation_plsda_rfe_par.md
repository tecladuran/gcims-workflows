Evaluating the Impact of Orthogonal Correction on Simulated Group
Discrimination
================
Tecla Duran Fort
2025-11-13

- <a href="#1-introduction" id="toc-1-introduction">1. Introduction</a>
- <a href="#2-load-data-and-metadata" id="toc-2-load-data-and-metadata">2.
  Load Data and Metadata</a>
  - <a href="#21-assign-classes" id="toc-21-assign-classes">2.1. Assign
    Classes</a>
  - <a href="#visualizing-group-vs-technical-effects"
    id="toc-visualizing-group-vs-technical-effects">Visualizing Group vs
    Technical Effects</a>
- <a href="#3-function-definitions" id="toc-3-function-definitions">3.
  Function Definitions</a>
  - <a href="#31-simulate-controlled-shift-in-one-cluster"
    id="toc-31-simulate-controlled-shift-in-one-cluster">3.1 Simulate
    Controlled Shift in One Cluster</a>
  - <a href="#32-compute-auc-from-a-pls-da-model"
    id="toc-32-compute-auc-from-a-pls-da-model">3.2 Compute AUC from a
    PLS-DA Model</a>
  - <a href="#33-evaluate-auc-before-and-after-correction"
    id="toc-33-evaluate-auc-before-and-after-correction">3.3 Evaluate AUC
    Before and After Correction</a>
    - <a href="#procedure" id="toc-procedure">Procedure</a>
- <a href="#5-run-simulation-across-clusters"
  id="toc-5-run-simulation-across-clusters">5. Run Simulation Across
  Clusters</a>
  - <a href="#51-visualisation-of-auc-by-cluster"
    id="toc-51-visualisation-of-auc-by-cluster">5.1. Visualisation of AUC by
    Cluster</a>
- <a href="#6-robustness-check" id="toc-6-robustness-check">6. Robustness
  Check</a>
  - <a href="#61-repeated-random-assignments"
    id="toc-61-repeated-random-assignments">6.1. Repeated Random
    Assignments</a>
  - <a href="#62-global-validation-across-clusters-and-perturbation-levels"
    id="toc-62-global-validation-across-clusters-and-perturbation-levels">6.2.
    Global Validation Across Clusters and Perturbation Levels</a>
    - <a href="#run-simulations" id="toc-run-simulations">Run Simulations</a>
    - <a href="#summary-statistics-by-perturbation-level"
      id="toc-summary-statistics-by-perturbation-level">Summary Statistics by
      Perturbation Level</a>

# 1. Introduction

This document evaluates the impact of an orthogonal projection
correction applied in the sample space to GC-IMS peak table data. The
analysis is based on a dataset where a controlled perturbation is
introduced in a specific signal to simulate a potential biomarker shift.

We begin by applying this perturbation to a selected cluster and
splitting the dataset into training and test sets using stratified
sampling. A PLS-DA model is then trained on the raw data and evaluated
through ROC curves and AUC values. The same procedure is repeated after
correcting the data using a two-step orthogonal projection based on
elapsed time and batch.

To better understand the model behaviour, we visualize the latent space
projections (PLS-DA scores) before and after correction. Finally, this
process is extended to multiple clusters and perturbation levels to
assess the consistency of the correction effect, and summary statistics
are reported across all simulations.

# 2. Load Data and Metadata

The dataset used in this simulation corresponds to a peak table derived
from repeated measurements of a pooled urine sample, with external
variables `elapsed_time` and `batch` included.

To ensure reproducibility of the initial random group assignment and all
subsequent analyses, we explicitly fix the random seed at the beginning
of the workflow (`set.seed(1234)`). This allows the entire simulation
pipeline to be rerun with identical outputs.

## 2.1. Assign Classes

In this first step, artificial group labels are randomly assigned to the
samples. This creates a balanced division into “Control” and “Cancer”
groups, which do not correspond to any real condition but serve to
simulate a classification scenario.

This is the **initial random assignment**, used for the first round of
simulations and visualizations. In later sections, the group labels will
be reassigned multiple times with different random seeds to evaluate the
impact of random variation and ensure that results are not dependent on
a specific group configuration.

## Visualizing Group vs Technical Effects

Before running the simulation, we visualize the alignment between the
randomly assigned `Group` labels and the technical variables
`elapsed_time` and `batch`. All values are scaled to the same range for
direct comparison.

![](simulation_plsda_rfe_par_files/figure-gfm/plot-alignment-1.png)<!-- -->

# 3. Function Definitions

## 3.1 Simulate Controlled Shift in One Cluster

This function modifies the values of a selected cluster in the “Cancer”
group by adding a perturbation proportional to its mean intensity. The
goal is to simulate a consistent increase in signal for one variable,
mimicking the presence of a potential biomarker.

## 3.2 Compute AUC from a PLS-DA Model

This section describes the procedure used to evaluate the classification
performance of a PLS-DA model using a **nested cross-validation
framework** combined with **iterative feature elimination based on VIP
scores**. Although the code returns fold-wise predictions, the final AUC
is computed afterwards from these predictions.

**1. Stratified fold generation**

The function `make_folds()` creates *k* stratified folds by:

- identifying sample indices for each class,
- shuffling them,
- splitting them into *k* equally sized groups,
- and merging the groups across classes.

This ensures that every fold preserves the original class balance
between *Control* and *Cancer*, preventing bias during training and
evaluation.

**2. Outer cross-validation loop**

The main evaluation is performed through an outer *k*-fold
cross-validation:

- In each outer fold, one subset is held out as the **test set**.
- The remaining samples form the **training set**.
- Only the predictions obtained on these outer test sets are used for
  the final performance estimation (e.g., AUC).

This provides an unbiased assessment of the model’s generalisation
capability.

**3. Recursive feature elimination based on VIP**

Within each outer fold, variable selection is performed on the training
set through an iterative procedure:

**3.1 Inner cross-validation**

A second level of cross-validation (inner folds) is used to:

- select the optimal number of latent variables (LVs),
- and evaluate the performance of the current set of features.

For every inner fold:

1.  A PLS-DA model is trained with the current feature set.
2.  Predictions are computed for all possible numbers of components up
    to `ncomp_max`.
3.  Classification accuracy is computed for each number of components.

The best-performing number of components (based on mean inner accuracy)
is retained.

**3.2 VIP-based elimination**

A full PLS-DA model is then fitted using all inner-optimal parameters.
VIP scores are computed, and a small fraction of the least informative
variables (5% or at least one) is removed.

This process repeats until:

- no improvement in inner CV accuracy is observed, or
- only one feature remains.

The best feature subset and the corresponding optimal number of
components are stored for the outer fold.

**4. Final model training within each outer fold**

Once the optimal feature subset and number of components are identified:

1.  A final PLS-DA model is trained on the full training set of that
    outer fold.
2.  Predictions are generated for the outer test set.
3.  The predicted probabilities, predicted classes, and true labels are
    stored.

These results constitute the fold-wise outputs used to compute the final
AUC.

**5. AUC computation**

After all outer folds are processed, the test-set predictions from all
folds are concatenated. The **Area Under the ROC Curve (AUC)** is then
computed using:

- true binary labels, and
- predicted probabilities for the positive class (*Cancer*).

The AUC summarises the discriminative power of the model:

- **0.5** → random guessing
- **1.0** → perfect separation

## 3.3 Evaluate AUC Before and After Correction

This function compares the classification performance of a PLS-DA model
trained on:

1.  The original but **modified** dataset, where a controlled shift has
    been introduced in one cluster to simulate a biomarker-like effect
    in the “Cancer” group.
2.  The same dataset after applying orthogonal projection in the sample
    space to remove the influence of technical variables.

The correction method is based on a custom implementation of projection
onto the orthogonal complement of the space spanned by external
covariates. This projection is performed sequentially for:

- `elapsed_time`, to account for intra-batch variation;
- `batch`, to account for inter-batch variation.

This correction strategy is described in detail here:  
[Orthogonal Correction in the Sample
Space](https://github.com/tecladuran/gcims-workflows/blob/main/docs/Correction/orthogonal_correction.md)

### Procedure

- The matrix of cluster intensities (`X_raw`) and group labels (`y`) are
  extracted from the **synthetically altered** dataset.
- A PLS-DA model is trained on the uncorrected data using the function
  `compute_auc_plsda()` (see Section 3.2), and its AUC is calculated.
- The data is then corrected using orthogonal projections for the two
  technical covariates.
- A second PLS-DA model is trained on the corrected data, and its AUC is
  computed using the same split.

The function returns a named numeric vector with two values: - `Raw`:
AUC obtained using the modified but uncorrected data; - `Corrected`: AUC
obtained after correction.

This allows for direct evaluation of how removing technical variation
affects the detectability of a synthetic class difference.

# 5. Run Simulation Across Clusters

In this step, we systematically evaluate how classification performance
changes as the simulated “biomarker” effect increases, both before and
after correction.

For each cluster and for a range of perturbation levels (from 0 to 1),
the following procedure is applied:

1.  A synthetic effect is introduced in the “Cancer” group for the
    selected cluster using `simulate_shift()`. The intensity increase is
    proportional to the cluster’s mean signal.
2.  The resulting modified dataset is passed to `evaluate_auc()`, which
    computes the AUC of a PLS-DA classifier before and after correction.
3.  This is repeated for all clusters and perturbation levels, producing
    a complete matrix of AUC values across different conditions.

To ensure reproducibility, a fixed seed (`seed = 42`) is used during
each evaluation.  
This means that the random partitioning into training and test sets is
**identical every time** the same simulation is run. As a result,
differences in performance between the raw and corrected data can be
attributed **solely to the correction itself**, not to randomness in the
group split.

The results are stored as a long-format dataframe (`auc_df`), with
columns: - `Cluster`: the cluster where the perturbation was applied -
`Percent`: the perturbation level (from 0 to 1) - `Condition`: “Raw” or
“Corrected” - `AUC`: the classification performance at that setting

This data is later used for visualizing performance trends across
clusters.

## 5.1. Visualisation of AUC by Cluster

This section visualizes how classification performance evolves for each
cluster as the simulated perturbation increases. For every cluster, we
plot the AUC of the PLS-DA classifier as a function of the perturbation
proportion, separately for the raw and corrected datasets.

Each plot shows two lines:

- **Raw** (uncorrected data)
- **Corrected** (after orthogonal projection)

These curves allow us to see how sensitive each cluster is to the
simulated effect, and how much the correction improves detection for
each level of perturbation. Clusters where the corrected AUC rises more
steeply indicate a stronger benefit from the correction method.

Each plot is titled with the cluster ID (e.g., “Cluster16”) and uses a
consistent y-axis range from 0.3 to 1.0 for comparability across
clusters.

![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-1.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-2.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-3.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-4.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-5.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-6.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-7.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-8.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-9.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-10.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-11.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-12.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-13.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-14.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-15.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-16.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-17.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-18.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-19.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-20.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-21.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-22.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-23.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-24.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-25.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-26.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-27.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-28.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-29.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-30.png)<!-- -->![](simulation_plsda_rfe_par_files/figure-gfm/plot-results-31.png)<!-- -->

# 6. Robustness Check

## 6.1. Repeated Random Assignments

To evaluate the stability of the observed AUC improvements, we repeat
the group assignment multiple times with different random seeds. This
controls for potential biases introduced by a single, lucky alignment
between group labels and technical variables.

We compute the AUC with and without orthogonal correction across 200
repetitions for a single cluster (`Cluster48`) perturbation of 20% and
visualize the distribution.

![](simulation_plsda_rfe_par_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## 6.2. Global Validation Across Clusters and Perturbation Levels

To confirm that the improvement provided by othogonal correction is
consistent across clusters and perturbation levels, we repeat the
simulation for all clusters and multiple perturbation levels, with
repeated random group assignments.

### Run Simulations

### Summary Statistics by Perturbation Level

To assess whether the correction provides a statistically significant
improvement in classification performance, we perform a **one-sided
paired Wilcoxon signed-rank test** at each perturbation level.

Each perturbation level (e.g., 0.05, 0.10, …) includes results from
**multiple clusters**, and for each cluster the simulation is repeated
**10 times** using different random group assignments. As a result,
every perturbation level includes a large number of paired AUC values
comparing the raw and corrected data under the exact same conditions
(same cluster, perturbation level, and group assignment).

This non-parametric test compares the AUC values before and after
correction for each simulation run, under the null hypothesis that the
correction does **not** increase performance. The test is **paired**,
meaning that each corrected AUC is compared to its corresponding raw AUC
under the exact same conditions (same cluster, perturbation level, and
group assignment), thereby controlling for random variation.

The Wilcoxon test is well suited for AUC values, which are bounded and
may be non-normally distributed. However, due to the **large number of
simulations per condition**, the test has very high statistical power.
This means that **p-values are often extremely small**, even when the
improvement is modest.

Therefore, **what truly matters in this context is not the p-value
itself, but the size and consistency of the improvement**. For this
reason, we report: - `mean_gain`: the average AUC increase after
correction. - `prop_improved`: the proportion of simulations in which
the corrected AUC is higher than the raw AUC.

Together, these metrics help quantify the practical benefit of the
correction across a range of perturbation intensities.

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Summary of Improvement After Correction
</caption>
<thead>
<tr>
<th style="text-align:right;">
Perturbation
</th>
<th style="text-align:right;">
p_value
</th>
<th style="text-align:right;">
mean_gain
</th>
<th style="text-align:right;">
prop_improved
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
1.0000
</td>
<td style="text-align:right;">
-0.0270
</td>
<td style="text-align:right;">
0.3774
</td>
</tr>
<tr>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.9356
</td>
<td style="text-align:right;">
-0.0066
</td>
<td style="text-align:right;">
0.4774
</td>
</tr>
<tr>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.3722
</td>
<td style="text-align:right;">
0.0010
</td>
<td style="text-align:right;">
0.5032
</td>
</tr>
<tr>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.4307
</td>
<td style="text-align:right;">
-0.0003
</td>
<td style="text-align:right;">
0.4806
</td>
</tr>
<tr>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.8261
</td>
<td style="text-align:right;">
-0.0019
</td>
<td style="text-align:right;">
0.4613
</td>
</tr>
<tr>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.9830
</td>
<td style="text-align:right;">
-0.0019
</td>
<td style="text-align:right;">
0.4484
</td>
</tr>
<tr>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.9972
</td>
<td style="text-align:right;">
-0.0020
</td>
<td style="text-align:right;">
0.3645
</td>
</tr>
<tr>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
1.0000
</td>
<td style="text-align:right;">
-0.0033
</td>
<td style="text-align:right;">
0.2839
</td>
</tr>
<tr>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.9983
</td>
<td style="text-align:right;">
-0.0027
</td>
<td style="text-align:right;">
0.2774
</td>
</tr>
<tr>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.9992
</td>
<td style="text-align:right;">
-0.0028
</td>
<td style="text-align:right;">
0.2677
</td>
</tr>
<tr>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.9985
</td>
<td style="text-align:right;">
-0.0025
</td>
<td style="text-align:right;">
0.2419
</td>
</tr>
<tr>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.9822
</td>
<td style="text-align:right;">
-0.0015
</td>
<td style="text-align:right;">
0.2323
</td>
</tr>
<tr>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.9968
</td>
<td style="text-align:right;">
-0.0015
</td>
<td style="text-align:right;">
0.2226
</td>
</tr>
<tr>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.9682
</td>
<td style="text-align:right;">
-0.0008
</td>
<td style="text-align:right;">
0.2226
</td>
</tr>
<tr>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.9925
</td>
<td style="text-align:right;">
-0.0009
</td>
<td style="text-align:right;">
0.1903
</td>
</tr>
<tr>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.9153
</td>
<td style="text-align:right;">
-0.0004
</td>
<td style="text-align:right;">
0.1742
</td>
</tr>
<tr>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.9993
</td>
<td style="text-align:right;">
-0.0007
</td>
<td style="text-align:right;">
0.1548
</td>
</tr>
<tr>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.9997
</td>
<td style="text-align:right;">
-0.0006
</td>
<td style="text-align:right;">
0.1387
</td>
</tr>
<tr>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.9377
</td>
<td style="text-align:right;">
-0.0005
</td>
<td style="text-align:right;">
0.1645
</td>
</tr>
<tr>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.9601
</td>
<td style="text-align:right;">
-0.0003
</td>
<td style="text-align:right;">
0.1419
</td>
</tr>
<tr>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.9560
</td>
<td style="text-align:right;">
-0.0001
</td>
<td style="text-align:right;">
0.1387
</td>
</tr>
</tbody>
</table>
