# GC-IMS Workflows — Documentation Overview

---

## Contents

### [Preprocessing](Preprocessing/)

- **[Full Workflow: From Raw Data to Peak Table (Markdown)](Preprocessing/Full_workflow.md)**  
  Complete pipeline from raw GC-IMS signals to a processed and corrected peak table.  
  *Includes:* Filtering, alignment, peak detection, clustering, imputation, and baseline correction.

- **[Baseline Correction (PDF)](Preprocessing/baseline_correction.md)**  
  Corrects overestimation in peak areas due to background signal.  
  *Method:* Cluster-specific residual volume estimation.


---

### [Linearity](Linearity/)
- **[Linearity Report (Markdown)](Linearity/linearity_report.md)**  
  Analyzes signal drift within and across batches using elapsed time, batch index, and storage time.  
  *Outcome:* Supports the use of batch as an ordinal temporal proxy.

---

### [Stability](Stability/)
- **[Stability Analysis (Markdown)](Stability/stability_analysis.md)**  
  Quantifies variability using RSD and variance explained by external factors.  
  *Result:* Only ~23% of clusters are stable (<20% RSD).

---

### [Correction](Correction/)
- **[EPO Correction Report (Markdown)](Correction/linear_correction.md)**  
  Implements linear orthogonalization to remove variability due to elapsed time and batch.  
  *Validation:* PCA and explained variance analysis before/after correction.

- **[Comparison of Corrected vs Raw Data (Markdown)](Correction/correction_stability_comparison.md)**  
  Assesses RSD and PCA differences pre/post correction.  
  *Result:* Stable clusters increase from 22.6% to 71%.

---

### [Group Simulation](Group_Simulation/)
- **[Simulation PLS-DA Report (Markdown)](Group_Simulation/simulation_plsda.md)**  
  Simulates subtle biomarker shifts and compares classification performance before and after correction.  
  *Metric:* AUC (PLS-DA).  
  *Result:* Systematic improvement with correction.

---

## Author

**Tecla Duran Fort**  
Universitat de Barcelona
Signal and Information Processing for Sensing Systems – IBEC  
Barcelona, 2025
