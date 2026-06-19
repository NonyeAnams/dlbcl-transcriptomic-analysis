# Transcriptomic and Survival Analysis of Diffuse Large B-Cell Lymphoma (DLBCL) Using Python

An end-to-end reproducible bioinformatics pipeline built in Python to process microarray expression data, identify differentially expressed transcripts, perform pathway enrichment mapping, and model clinical time-to-event risk using public datasets.

**Dataset:** GEO GSE10846
**Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
**Cohort analyzed:** 200 R-CHOP treated DLBCL patients


## Abstract

**Background:** Diffuse Large B-Cell Lymphoma (DLBCL) exhibits profound molecular heterogeneity, marked by distinct cell-of-origin phenotypes including Activated B-Cell-like (ABC) and Germinal Center B-Cell-like (GCB) subtypes. While R-CHOP stands as the standard therapeutic intervention, clinical outcomes vary dramatically across these groups. 

**Methods:** An end-to-end reproducible pipeline was built in Python using public microarray expression profiles from GEO Series GSE10846 (N=200 R-CHOP treated patients). Transcripts were aggregated using mean probe intensities, and unsupervised global variance profiles were mapped using Principal Component Analysis (PCA). Supervised differential expression testing was conducted via parallel Ordinary Least Squares (OLS) modeling with Benjamini-Hochberg False Discovery Rate (FDR) adjustments. Biological programs were determined using Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA). Clinical outcome correlation was modeled using non-parametric Kaplan-Meier estimation and semi-parametric Multivariate Cox Proportional Hazards regression. 

**Results:** Supervised statistical modeling revealed 4,901 significantly dysregulated transcripts (FDR < 0.05), with canonical markers (*BCL6*, *MME*, *SERPINA9*) overexpressed in GCB and inflammatory/activation network markers (*BATF*, *IRF4*) elevated in ABC. GSEA mapped systemic pathobiological programs driving the aggressive ABC state to MYC targets, IL-6/JAK/STAT3 signaling, and Interferon-Gamma (IFN-γ) responses. Whole-cohort survival models implicated both MYC and IFN-γ pathway scores as univariate liabilities. However, subtype-stratified Kaplan-Meier analysis uncovered a molecular paradox: higher IFN-γ signaling confers a significant survival benefit within the GCB subtype ($p=0.0468$) but breaks down into non-significance within ABC ($p=0.4312$). Final multivariate Cox regression resolved this relationship, suggesting that the individual pathway signatures are heavily confounded by underlying disease epidemiology. Controlling for covariates isolated the binary cell-of-origin subtype as the dominant independent predictor, with the ABC entity conferring a striking 4.34-fold increase in risk of mortality ($p=0.0062$). 

**Conclusions:** This project establishes a robust, reproducible multi-stage computational framework linking raw transcriptomic feature metrics to macroscopic pathway architecture and definitive, stratified clinical hazard outcomes.


## Key Analytical Outputs

![Volcano Plot](../results/figures/volcano_plot_rchop.png)
![Leading-Edge Heatmap](../results/figures/ifn_heatmap.png)
![Risk Stratified KM Curve](../results/figures/KM_table_risk_stratified.png)


## Key Findings

- 4,901 genes were significantly dysregulated between ABC and GCB DLBCL subtypes (FDR < 0.05).
- Known GCB markers (BCL6, MME, SERPINA9) were highly expressed in GCB tumors.
- ABC tumors demonstrated enrichment of inflammatory and immune-related pathways.
- GSEA identified MYC Targets, IL6/JAK/STAT3 signaling, and IFN-γ response as dominant programs.
- Elevated IFN-γ activity was associated with improved survival within the GCB subtype.
- Multivariate Cox modeling identified molecular subtype as the strongest independent predictor of outcome (HR = 4.34).


## Analysis Workflow

| Dataset | Processing | Exploration | Statistics | Biology | Survival | Prognosis |
|----------|----------|----------|----------|----------|----------|----------|
| GSE10846 | Preprocessing | PCA | Differential Expression | ORA + GSEA | Kaplan-Meier & Cox | Risk Stratification |

  
## Pipeline Architecture & Core Visualizations

The analytical engine is structured into 5 sequential Jupyter Notebooks, progressing from raw data acquisition to clinical translation:

### 1. Data Cleaning & Alignment Verification (`notebooks/1_preprocessing.ipynb`)
* Automated download of GSE10846 and GPL570 platform configurations.
* Aggregation of 23,520 probe intensities into single gene-level mean profiles.
* Extraction and curation of clinical variables, subsetting to a homogeneous cohort (N=200 R-CHOP-treated samples).

### 2. Global Exploratory Data Analysis (`notebooks/2_pca_exploratory.ipynb`)
* Isolated the top 1,000 Highly Variable Genes (HVGs) based on cross-sample variance.
* Unsupervised PCA projections revealed that molecular subtype coordinates group along PC2, whereas clinical stage and survival status are stochastically distributed.

### 3. Supervised Differential Expression Testing (`notebooks/3_differential_expression.ipynb`)
* Fit 23,520 individual linear models via OLS to evaluate transcriptional shifts between subtypes.
* Controlled type-I errors using the Benjamini-Hochberg FDR procedure, discovering 4,901 significantly dysregulated genes.

### 4. Pathway Mapping & Functional Landscapes (`notebooks/4_pathway_enrichment.ipynb`)
* Conducted Over-Representation Analysis (ORA) via Enrichr to map highly altered transcripts ($|\log_2\text{FC}| > 1.0$).
* Executed continuous GSEA Prerank modeling against MSigDB Hallmark sets, revealing systemic upregulation of cell-cycle accelerators (MYC targets), pro-survival signaling (IL-6/JAK/STAT3), and microenvironmental inflammation (IFN-γ Response) in the ABC subtype.

### 5. Time-to-Event Survival Analysis (`notebooks/5_survival_analysis.ipynb`)
* Evaluated continuous and categorical variables using Kaplan-Meier survival curves and Cox proportional hazard models.
* Uncovered a critical biological inversion: high IFN-γ expression improves survival outcomes for GCB patients ($p=0.0468$) but has no impact on ABC outcomes ($p=0.4312$).
* Implemented multivariate modeling to reveal that MYC and interferon pathobiological signatures are epidemiologically confounded by the underlying cell-of-origin distribution, identifying the ABC subtype as the true independent driver of clinical risk ($\text{HR} = 4.34$, $p=0.0062$).
* Validated mathematical stability by running Schoenfeld residual diagnostics, confirming strict compliance with the proportional hazards assumption.


## Reproducibility Guide

To execute this pipeline locally, clone the repository and install the tracked dependencies inside a virtual environment:

```bash
git clone [https://github.com/nonyeanams/dlbcl-transcriptomic-analysis.git](https://github.com/nonyeanams/dlbcl-transcriptomic-analysis.git)
cd dlbcl-transcriptomic-analysis
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```


## Repository Structure

```
dlbcl-transcriptomic-analysis/
├── README.md              
├── requirements.txt         
├── data/
│   ├── raw/                
│   └── processed/          
├── notebooks/
│   ├── 1_data_preprocessing.ipynb
│   ├── 2_exploratory_pca.ipynb
│   ├── 3_differential_expression.ipynb
│   ├── 4_pathway_enrichment.ipynb
│   └── 5_survival_analysis.ipynb
└── results/
    ├── figures/            
    └── tables/            
```


## Limitations

- Analysis was performed on microarray rather than RNA-seq data.
- Findings were derived from a single public cohort and were not externally validated.
- Pathway scores were generated using simple gene-set averaging rather than more advanced methods such as GSVA or ssGSEA.
- The study is observational and cannot establish causal biological mechanisms.


## Future Work

Potential extensions include:

- Validation in independent DLBCL cohorts
- RNA-seq based replication using TCGA or GEO datasets
- GSVA/ssGSEA pathway scoring
- Development of multi-gene prognostic signatures
- Integration of mutation and transcriptomic profiles
- Machine learning models for risk prediction