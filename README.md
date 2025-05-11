# Machine learning for early detection of leukaemia from RNA-seq data using predictive and survival modelling.
This project applies machine learning techniques to RNA-seq gene expression data for early detection and prognosis of Acute Myeloid Leukaemia (AML). It includes differential gene expression analysis, classification of FLT3-ITD mutation status using multiple ML models, SHAP-based interpretability, functional enrichment, and survival modelling using TCGA-LAML data. ![image](https://github.com/user-attachments/assets/fc7e824f-d0ce-4b88-aa60-08a7342f8fca)

# Project Overview
Acute Myeloid Leukaemia (AML) is a genetically heterogeneous and aggressive blood cancer characterised by uncontrolled proliferation of myeloid precursor cells and poor patient outcomes. Despite advances in genomics, timely diagnosis and effective risk stratification remain major clinical challenges. Current diagnostic tools often overlook early molecular signals, particularly those associated with high-risk mutations such as FLT3-ITD.
This study aimed to address this gap by developing machine learning (ML) models for early AML classification and prognostic biomarker discovery using RNA sequencing data. Five algorithms, Random Forest (RF), XGBoost (XGB), Support Vector Machine (SVM), Least Absolute Shrinkage and Selection Operator (LASSO) and Neural Network (NN), were trained and evaluated, producing a consensus set of top-ranked genes. SHAP (SHapley Additive exPlanations) analysis revealed PDE4B, SOCS2, SPINK2, MIR155 and HOX Genes as key predictors contributing to AML classification.
Functional enrichment analyses (Gene Ontology, KEGG, Reactome, DISGENET) confirmed strong associations with haematopoiesis, angiogenesis and immune pathways, particularly involving the Renin-Angiotensin System and interleukin signalling. DISGENET findings also revealed subtype-specific links to AML-M1 and AML-M2. Kaplan-Meier and Cox regression survival analyses identified CIBAR1P1, SOCS2, MRC1 and LGALS3BP as significantly associated with prognosis. A LASSO-derived 8-gene signature stratified patients into high- and low-risk groups with strong predictive power (p = 2e-04; C-index = 0.72).
This integrative ML and transcriptomic approach successfully uncovered novel and known AML biomarkers, offering a robust framework to improve diagnosis, risk prediction and personalised treatment strategies in AML.

## Key Features
- RNA-seq preprocessing and DEG filtering
- Supervised ML model training:  
  - Random Forest  
  - XGBoost  
  - SVM  
  - LASSO  
  - Neural Networks
- SMOTE oversampling for class imbalance
- SHAP-based model interpretation
- Consensus gene selection across models
- Functional enrichment (GO, KEGG, Reactome, DISGENET)
- Survival analysis using TCGA-LAML:
  - Kaplan–Meier plots
  - Cox regression (univariate & multivariate)
  - LASSO-derived multigene prognostic signature

## Tools and Languages
- **R**: limma, survival, glmnet, ggplot2, VennDiagram  
- **Python**: scikit-learn, shap, pandas, matplotlib  
- **DAVID**: Functional enrichment  

## Results
Key results include:
- DEGs result (Volcano plot)
- Top 20 feature selections
- ROC and PR curves for model performance
- SHAP plots for gene importance
- Venn diagrams for cross-model gene overlap
- Functional enrichment dot plots
- Kaplan–Meier and forest plots for prognostic validation

## Data Access

This project uses publicly available gene expression and clinical data from **GEO** and **TCGA-LAML** via [https://www.ncbi.nlm.nih.gov/geo/] and [https://portal.gdc.cancer.gov/projects/TCGA-LAML] See `Rscript scripts/data_preprocessin.R` for instructions on how to download and format the data.

## How to Run
- **Run data collection/preprocessing in R**: Rscript scripts/data_preprocessin.R
- **Run differential expression in R**: Rscript scripts/differential_expression.R
- **Run machine learning models in R**: Rscript scripts/machine_learning_models.R
- **Run SHAP interpretability in Python**: python scripts/shap_interpretation.py
- **Run model comparism**: Rscript scripts/model comparism.R
- **Run enrichment analysis**: Rscript scripts/enrichment_analysis.R
- **Run survival analysis**: Rscript scripts/survival_analysis.R
- **Run lasso cox model**: Rscript scripts/Run lasso cox model.R

# Acknowledgements
This work was completed as part of the MSc Bioinformatics program at Teesside University. Special thanks to academic supervisors, Family and friends, GEO/TCGA data providers, open-source contributors.
