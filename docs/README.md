# GC-IMS Workflows — Documentation Overview

---

## Contents

### [Preprocessing and Feature Extraction](Preprocessing/)

- **[Full Workflow: From Raw Data to Peak Table](Preprocessing/Full_workflow.md)**  
  Complete pipeline from raw GC-IMS signals to a processed and corrected peak table.  
  *Includes:* Filtering, alignment, peak detection, clustering, imputation, and baseline correction.

- **[Baseline Correction](Preprocessing/baseline_correction.md)**  
  Corrects overestimation in peak areas due to background signal.  
  *Method:* Cluster-specific residual volume estimation.


---

### [Stability](Stability/)
- **[Stability Analysis](Stability/stability_analysis.md)**  
  Quantifies variability using RSD.   
  *Result:* Only ~23% of clusters are stable (<20% RSD).

---

### [Analysis of Variability Sources ](Variability_Sources/)
- **[Linearity Report](Variability_Sources/linearity_report.md)**  
  Analyzes signal drift within and across batches using elapsed time, batch index, and storage time.  
  *Outcome:* Supports the use of batch as an ordinal temporal proxy.
- **[Variability Sources](Variability_Sources/variability_sources.md)**
  Quantifies variability in cluster intensities explained by elapsed time and batch effects using R² and PCA.
  *Outcome: Elapsed time and batch are major sources of signal variability (70%).*

---

### [Correction](Correction/)
- **[Orthogonal Correction Report](Correction/orthogonal_correction.md)**  
  Implements linear orthogonalization to remove variability due to elapsed time and batch.  
  *Validation:* PCA and explained variance analysis before/after correction.

- **[Comparison of Corrected vs Raw Data](Correction/correction_stability_comparison.md)**  
  Assesses RSD and PCA differences pre/post correction.  
  *Result:* Stable clusters increase from 22.6% to 71%.

---

### [Group Simulation](Group_Simulation/)
- **[Simulation PLS-DA Report](Group_Simulation/simulation_plsda_rfe_par.md)**  
  Simulates subtle biomarker shifts and compares classification performance before and after correction.  
  *Metric:* AUC (PLS-DA).  
  *Result:* Systematic improvement with correction.

---

## Author

**Tecla Duran Fort**  
Universitat de Barcelona
Signal and Information Processing for Sensing Systems – IBEC  
Barcelona, 2025
