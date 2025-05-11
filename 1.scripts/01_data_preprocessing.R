# TITLE: Machine Learning (ML) models for early detection of leukaemia from RNA sequencing data.


# Acute Myeloid Leukaemia (AML) is an aggressive and heterogeneous blood cancer characterised by the uncontrolled proliferation of immature myeloid cells, leading to bone marrow failure and high mortality.
# One of the most clinically significant subtypes is defined by FLT3-ITD mutations, which are associated with rapid disease progression, poor prognosis, and resistance to treatment.
# While RNA sequencing (RNA-seq) has advanced our understanding of AML at the transcriptomic level, traditional diagnostic workflows often fail to capture the complexity of high-risk cases.

# ML offers a powerful solution by uncovering complex gene expression patterns that inform classification and prognosis.
# This study presents a comparative ML framework for classifying FLT3-ITD mutation status in AML using RNA-seq data. Multiple algorithms (RF, XGB, SVM, LASSO and NN) were trained through DEG analysis.
# SHAP-based interpretation was applied to identify consistent predictive biomarkers.

# To assess the clinical and biological relevance of these genes, functional enrichment analysis and survival modelling were conducted.
# This approach aims to improve early diagnosis and risk stratification in AML through interpretable and biologically grounded ML models.

# All analyses, including RNA-seq data preprocessing, DEG, ML, cross-model comparison, biologicaL and survival modelling, were conducted using R.
# with packages such as caret, glmnet, randomForest, xgboost and survival).
# Functional enrichment analysis was performed using the DAVID Bioinformatics Resources via its online interface and custom gene lists exported from R.

# Data Collection and Preprocessing

# Load All Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "biomaRt", "ggplot2", "dplyr"))

# Load libraries
library(GEOquery)
library(biomaRt)
library(limma)
library(ggplot2)
library(dplyr)
library(caret)
library(smotefamily)
library(ggrepel)

# Download and Load GEO Data
gse <- getGEO("GSE6891", GSEMatrix = TRUE, getGPL = FALSE)
expr_matrix <- exprs(gse[[1]])
pheno <- pData(gse[[1]])

# Check dimensions clearly
dim(expr_matrix)  

# Manual download method (alternative if above fails)
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6891/matrix/GSE6891_series_matrix.txt.gz",
              destfile = "GSE6891_series_matrix.txt.gz")
gse <- getGEO(filename = "GSE6891_series_matrix.txt.gz")
expr_matrix <- exprs(gse)
pheno <- pData(gse)

# Check dimensions clearly
dim(expr_matrix)  

# Read top lines to find where the table begins
lines <- readLines("GSE6891_series_matrix.txt.gz")

start_line <- grep("!series_matrix_table_begin", lines)
end_line <- grep("!series_matrix_table_end", lines)

# Extract only the expression matrix block
expr_data <- read.delim("GSE6891_series_matrix.txt.gz", 
                        skip = start_line, 
                        nrows = end_line - start_line - 1,
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# Set probe IDs as rownames
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]

# Convert to numeric matrix
expr_matrix <- as.matrix(sapply(expr_data, as.numeric))
rownames(expr_matrix) <- rownames(expr_data)

# Check dimensions
dim(expr_matrix)

set.seed(42)  # For reproducibility
# expr_mapped is your gene x sample matrix

# Select 5 random genes to test
sample_genes <- sample(rownames(expr_mapped), 5)

# Run Shapiro-Wilk test for each gene
shapiro_results <- sapply(sample_genes, function(gene) {
  shapiro.test(as.numeric(expr_mapped[gene, ]))$p.value
})

# Display p-values
shapiro_results

# For a single gene
gene <- sample_genes[1]
qqnorm(as.numeric(expr_mapped[gene, ]), main = paste("QQ Plot for", gene))
qqline(as.numeric(expr_mapped[gene, ]), col = "red")

# Probe-to-Gene Mapping

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

probe_ids <- rownames(expr_matrix)
chunks <- split(probe_ids, ceiling(seq_along(probe_ids)/1000))

mapped_list <- lapply(chunks, function(chunk) {
  getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
        filters = "affy_hg_u133_plus_2",
        values = chunk,
        mart = ensembl)
})

mapped_probes <- do.call(rbind, mapped_list) %>%
  distinct(affy_hg_u133_plus_2, .keep_all = TRUE)

# Map probes to genes and filter
expr_filtered <- expr_matrix[rownames(expr_matrix) %in% mapped_probes$affy_hg_u133_plus_2, ]
gene_symbols <- mapped_probes$hgnc_symbol[match(rownames(expr_filtered), mapped_probes$affy_hg_u133_plus_2)]
rownames(expr_filtered) <- gene_symbols

# Aggregate duplicated genes (average expression)
expr_mapped <- aggregate(expr_filtered, by = list(GeneSymbol = rownames(expr_filtered)), FUN = mean)
rownames(expr_mapped) <- expr_mapped$GeneSymbol
expr_mapped <- expr_mapped[, -1]

dim(expr_mapped)  # Should be Genes x Samples

# Phenotype data and Group processing

group_simplified <- trimws(pheno$`flt3 itd mutation:ch1`)
table(group_simplified)

valid_samples <- which(!is.na(group_simplified) & group_simplified %in% c("Pos", "Neg"))
group <- factor(group_simplified[valid_samples])

# Subset data correctly
pheno_sub <- pheno[valid_samples, ]
expr_sub <- expr_mapped[, valid_samples]

table(group)  # Check group distribution

# Ensure expression data is numeric

expr_sub <- as.matrix(expr_sub)  
mode(expr_sub) <- "numeric"  # force numeric

# Confirm numeric matrix
is.numeric(expr_sub)  # Should return TRUE
