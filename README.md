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
This work was completed as part of the MSc Bioinformatics programme at Teesside University.
Special thanks to:
- Academic supervisors and programme leaders
- GEO and TCGA for open-access datasets
- The open-source bioinformatics community
- Family and friends for support throughout the project

## References
- Alshamleh, I., Kurrle, N., Makowka, P., Bhayadia, R., Kumar, R., Süsser, S., Seibert, M., Ludig, D., Wolf, S., Koschade, S.E. and Stoschek, K. (2023). PDP1 is a key metabolic gatekeeper and modulator of drug resistance in FLT3-ITD-positive acute myeloid leukemia. Leukemia, 37(12), pp.2367-2382.
- Bellantuono, L., Tommasi, R., Pantaleo, E., Verri, M., Amoroso, N., Crucitti, P., Di Gioacchino, M., Longo, F., Monaco, A., Naciu, A.M. and Palermo, A. (2023). An eXplainable Artificial Intelligence analysis of Raman spectra for thyroid cancer diagnosis. Scientific Reports, 13(1), p.16590.
- Borah, K., Das, H.S., Seth, S., Mallick, K., Rahaman, Z. and Mallik, S. (2024). A review on advancements in feature selection and feature extraction for high-dimensional NGS data analysis. Functional & Integrative Genomics, 24(5), p.139.
- Cheng, Y., Xu, S.M., Santucci, K., Lindner, G. and Janitz, M. (2024). Machine learning and related approaches in transcriptomics. Biochemical and Biophysical Research Communications, p.150225.
- Dakal, T.C., Xu, C. and Kumar, A. (2025). Advanced computational tools, artificial intelligence and machine-learning approaches in gut microbiota and biomarker identification. Frontiers in Medical Technology, 6, p.1434799.
- Ding, Z., Wu, J., Ye, Y., Zhong, Y., Yan, L. and Wang, P. (2025). A novel signature predicts prognosis in pancreatic cancer based on tumor membrane-associated genes. Heliyon, 11(4).
- Duan, Y. and Xu, X. (2023). A signature based on anoikis-related genes for the evaluation of prognosis, immunoinfiltration, mutation, and therapeutic response in ovarian cancer. Frontiers in Endocrinology, 14, p.1193622.
- Fan, Z. (2023). Elucidating complex biological interactions using computational techniques (Doctoral dissertation, University of Pittsburgh).
- Hemminki, K., Hemminki, J., Försti, A. and Sud, A. (2023). Survival in hematological malignancies in the Nordic countries through a half century with correlation to treatment. Leukemia, 37(4), pp.854-863.
- Kalyakulina, A., Yusipov, I., Kondakova, E., Bacalini, M.G., Franceschi, C., Vedunova, M. and Ivanchenko, M. (2023). Small immunological clocks identified by deep learning and gradient boosting. Frontiers in Immunology, 14, p.1177611.
- Liao, K.M., Liu, C.F., Chen, C.J., Feng, J.Y., Shu, C.C. and Ma, Y.S. (2023). Using an artificial intelligence approach to predict the adverse effects and prognosis of tuberculosis. Diagnostics, 13(6), p.1075.
- Lin, C., Lin, K., Li, P., Yuan, H., Lin, X., Dai, Y., Zhang, Y., Xie, Z., Liu, T. and Wei, C. (2024). A genomic instability-associated lncRNA signature for predicting prognosis and biomarkers in lung adenocarcinoma. Scientific Reports, 14(1), p.14460.
- Mak, T.K., Li, K., Zhao, Z., Wang, K., Zeng, L., He, Q., Lu, W., Chen, W., He, Y., Li, J. and Zhang, C. (2025). m6A demethylation of NNMT in CAFs promotes gastric cancer progression by enhancing macrophage M2 polarization. Cancer Letters, 611, p.217422.
- Martin, J. (2024). AlphaFold2 predicts whether proteins interact amidst confounding structural compatibility. Journal of Chemical Information and Modeling, 64(5), pp.1473-1480.
- Mei, Y., Liu, Y., Liu, W., Chen, M., Liu, X., Wang, S., Mou, J., Xing, H., Tang, K., Tian, Z. and Rao, Q. (2024). Identifying ADGRG1 as a specific marker for tumor-reactive T cells in acute myeloid leukemia. Experimental Hematology & Oncology, 13(1), p.92.
- Onoja, A. (2023). An integrated interpretable machine learning framework for high-dimensional multi-omics datasets.
- Raptis, S., Ilioudis, C. and Theodorou, K. (2024). From pixels to prognosis: unveiling radiomics models with SHAP and LIME for enhanced interpretability. Biomedical Physics & Engineering Express, 10(3), p.035016.
- Rizwan, S., Tabassum, F. and Islam, S. (2023). An Ensemble Method for Cancer Classification and Identification of Cancer-Specific Genes from Genomic Data (Doctoral dissertation, Department of Computer Science and Engineering (CSE), Islamic University of Technology (IUT), Board Bazar, Gazipur-1704, Bangladesh).
- Salmoiraghi, S. (2023). Molecular Genetics of Myeloid Malignancies Drives Modern Treatment Options.
- Shevyrev, D., Tereshchenko, V., Berezina, T.N. and Rybtsov, S. (2023). Hematopoietic stem cells and the immune system in development and aging. International Journal of Molecular Sciences, 24(6), p.5862.
- Silva, R., Riedel, C., Reboul, J., Guibert, B., Ruffle, F., Gallopin, M., Gilbert, N., Boureux, A. and Commes, T. (2024). Comparing machine learning models for predicting mutation status in Acute Myeloid Leukemia patients using RNA-seq data. bioRxiv, pp.2024-11.
- Sohal, I.S. and Kasinski, A.L. (2023). Emerging diversity in extracellular vesicles and their roles in cancer. Frontiers in oncology, 13, p.1167717.
- Suravajhala, P.N. and Bizzaro, J.W. (eds.) (2025). Next-generation sequencing: Standard operating procedures and applications. 1st ed. Boca Raton: CRC Press. Available at: https://doi.org/10.1201/9781003354062 (Accessed: 2 May 2025).
- Surblys, V., Kozłowski, E., Matijošius, J., Gołda, P., Laskowska, A. and Kilikevičius, A. (2024). Accelerometer-Based Pavement Classification for Vehicle Dynamics Analysis Using Neural Networks. Applied Sciences, 14(21), p.10027.
- Tsai, H.H.D., Oware, K.D., Wright, F.A., Chiu, W.A. and Rusyn, I. (2025). A workflow for human health hazard evaluation using transcriptomic data and Key Characteristics-based gene sets. Toxicological Sciences, p.kfaf036.
- Vieira, F.M.G. (2023). Unveiling Novel Glioma Biomarkers Through Multi-Omics Integration and Classification (Master's thesis, Universidade NOVA de Lisboa (Portugal)).
- Wang, T., Cui, S., Lyu, C., Wang, Z., Li, Z., Han, C., Liu, W., Wang, Y. and Xu, R. (2024). Molecular precision medicine: Multi-omics-based stratification model for acute myeloid leukemia. Heliyon, 10(17).
- Xie, J., Zhu, Y., Yang, Z., Yu, Z., Yang, M. and Wang, Q. (2025). An integrative analysis reveals cancer risk associated with artificial sweeteners. Journal of Translational Medicine, 23(1), p.32.
- Yin, C., Li, Y., Zhang, C., Zang, S., Wang, Z., Yan, X., Ma, T., Li, X. and Li, W. (2024). Sequential gene expression analysis of myelodysplastic syndrome transformation identifies HOXB3 and HOXB7 as the novel targets for mesenchymal cells in disease. BMC cancer, 24(1), p.111.
- Yu, J., Mei, J., Zuo, D., Zhang, M., Yu, S., Li, F., Wang, J., Bi, D., Ma, S., Wang, J. and Yin, Z.J. (2024). Inflammatory factor-mediated miR-155/SOCS1 signaling axis leads to Treg impairment in systemic lupus erythematosus. International Immunopharmacology, 141, p.113013.
- Zhang, M., Zhu, L., Liang, S., Mao, Z., Li, X., Yang, L., Yang, Y., Wang, K., Wang, P. and Chen, W. (2024). Pulmonary function test-related prognostic models in non-small cell lung cancer patients receiving neoadjuvant chemoimmunotherapy. Frontiers in Oncology, 14, p.1411436.
- Zhu, Z., Jiang, L. and Ding, X. (2023). Advancing breast cancer heterogeneity analysis: Insights from genomics, transcriptomics and proteomics at bulk and single-cell levels. Cancers, 15(16), p.4164.
