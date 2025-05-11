# Survival Analysis Using the TCGA-LAML Cohort

# Survival analysis was performed on TCGA-LAML data to evaluate the prognostic value of consensus genes. Patients were stratified by median gene expression and Kaplan–Meier curves were used to compare survival.
# Univariate and multivariate Cox regression models were applied to assess individual and adjusted gene-level associations with overall survival, with results visualised using forest plots.


# 📦 Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(survival)
library(survminer)

# 📥 Step 1: Query & Download RNA-seq Data
query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# 📊 Step 2: Extract Clinical Data for Survival
clinical_df <- as.data.frame(colData(data))
clinical_surv <- data.frame(
  sample_id = rownames(clinical_df),
  time = clinical_df$days_to_death,
  status = ifelse(clinical_df$vital_status == "Alive", 0, 1)
)
clinical_surv <- clinical_surv[!is.na(clinical_surv$time), ]

# 🧬 Step 3: Extract Expression Matrix & Match Samples
expr_matrix <- assay(data, "unstranded")
shared_samples <- intersect(colnames(expr_matrix), clinical_surv$sample_id)
expr_matrix <- expr_matrix[, shared_samples]
clinical_surv <- clinical_surv[clinical_surv$sample_id %in% shared_samples, ]

# 🧠 Step 4: Map Ensembl IDs to Gene Symbols
rownames(expr_matrix) <- gsub("\\..*$", "", rownames(expr_matrix))
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(expr_matrix),
  mart = mart
)
gene_map <- gene_map[gene_map$hgnc_symbol != "", ]
expr_matrix <- expr_matrix[rownames(expr_matrix) %in% gene_map$ensembl_gene_id, ]
rownames(expr_matrix) <- gene_map$hgnc_symbol[match(rownames(expr_matrix), gene_map$ensembl_gene_id)]

# 🧬 Step 5: The 13 consistently selected genes across all models
common_genes <- Reduce(intersect, list(top20_genes_rf, top20_genes_xgb, svm_top20_genes, lasso_top20_genes, nn_top20_genes))
cat("Consistently selected genes across models:", common_genes, "\n")


# 🔗 Step 6: Merge Expression & Clinical Data
expr_top <- as.data.frame(t(expr_matrix[common_genes, ]))
expr_top$sample_id <- rownames(expr_top)
merged_data <- merge(clinical_surv, expr_top, by = "sample_id")


# Kaplan-Meier Survival for 12 Intersect Genes

# Optional: Save all plots to a PDF file

for (gene in setdiff(common_genes, "MIR155")) {
  
  gene_expr <- merged_data[[gene]]
  
  # Skip genes with missing or constant expression
  if (all(is.na(gene_expr)) || length(unique(gene_expr)) < 2) next
  
  # Divide into High/Low by median
  median_val <- median(gene_expr, na.rm = TRUE)
  merged_data$group <- ifelse(gene_expr >= median_val, "High", "Low")
  
  # Skip if no variation in group
  if (length(unique(merged_data$group)) < 2) next
  
  # Fit survival
  fit <- survfit(Surv(time, status) ~ group, data = merged_data)
  
  # Plot
  print(
    ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
               title = paste("Kaplan-Meier Survival: ", gene),
               legend.title = "Expression Level",
               legend.labs = c("High", "Low"))
  )
}

# Kaplan-Meier Survival Plots
library(gridExtra)

# Create an empty list to store KM plots
km_plot_list <- list()

# Loop through genes (excluding MIR155)
for (gene in setdiff(common_genes, "MIR155")) {
  
  gene_expr <- merged_data[[gene]]
  
  # Skip genes with NA or constant expression
  if (all(is.na(gene_expr)) || length(unique(gene_expr)) < 2) next
  
  # Create High/Low group based on median
  median_val <- median(gene_expr, na.rm = TRUE)
  merged_data$group <- ifelse(gene_expr >= median_val, "High", "Low")
  
  # Skip if grouping fails
  if (length(unique(merged_data$group)) < 2) next
  
  # Fit survival model
  fit <- survfit(Surv(time, status) ~ group, data = merged_data)
  
  # Create KM plot
  km_plot <- ggsurvplot(
    fit,
    data = merged_data,
    pval = TRUE,
    pval.coord = c(900, 0.4),      # Adjust X/Y for your data range
    pval.size = 4,
    legend.title = "Expression",
    legend.labs = c("High", "Low"),
    title = gene
  )
  
  # Adjust title formatting
  km_plot$plot <- km_plot$plot + theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  )
  
  # Store only the plot part
  km_plot_list[[gene]] <- km_plot$plot
}


# Display first 6 KM plots in a grid
grid.arrange(grobs = km_plot_list[1:6], ncol = 3)


# Get total number of plots
n <- length(km_plot_list)

# Display the last 6
grid.arrange(grobs = km_plot_list[(n-5):n], ncol = 3)

# Univariate Cox Regression for 12 Genes

# Create empty results list
cox_results <- data.frame(
  Gene = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through 12 genes (excluding MIR155)
for (gene in setdiff(common_genes, "MIR155")) {
  
  gene_expr <- merged_data[[gene]]
  
  # Skip genes with constant or NA values
  if (all(is.na(gene_expr)) || length(unique(gene_expr)) < 2) next
  
  # Cox model
  cox_model <- coxph(Surv(time, status) ~ gene_expr, data = merged_data)
  summary_model <- summary(cox_model)
  
  # Extract results
  HR <- summary_model$coefficients[,"exp(coef)"]
  CI <- summary_model$conf.int[,c("lower .95", "upper .95")]
  p <- summary_model$coefficients[,"Pr(>|z|)"]
  
  # Append to results
  cox_results <- rbind(cox_results, data.frame(
    Gene = gene,
    HR = HR,
    CI_lower = CI[1],
    CI_upper = CI[2],
    p_value = p
  ))
}

# View results sorted by p-value
cox_results_13 <- cox_results[order(cox_results$p_value), ]
print(cox_results_13)

# Optional: Save to CSV
write.csv(cox_results_13, "TCGA_LAML_Cox_Regression_TopGenes.csv13", row.names = FALSE)


# Assuming cox_results (univariate Cox regression) is ready
library(ggplot2)

# Order genes by HR or p-value (optional)
cox_results <- cox_results[order(cox_results$p_value), ]
cox_results$Gene <- factor(cox_results$Gene, levels = rev(cox_results$Gene))

ggplot(cox_results, aes(x = HR, y = Gene)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "darkblue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Univariate Cox Regression for 12 Consistently Selected Genes",
    x = "Hazard Ratio (HR)",
    y = "Genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )

# Load required package
library(forestplot)

# Select top N genes
top_n <- 12
top_cox <- head(cox_results, top_n)
top_cox$Gene


label_text <- cbind(
  c("Genes", as.character(top_cox$Gene)),  # <-- convert factor to character!
  c("Hazard Ratio (95% CI)", 
    sprintf("%.2f (%.2f–%.2f)", top_cox$HR, top_cox$CI_lower, top_cox$CI_upper)),
  c("p-value", 
    sprintf("%.4f", top_cox$p_value))
)


# Prepare data for plotting
hr_data <- data.frame(
  mean  = c(NA, top_cox$HR),
  lower = c(NA, top_cox$CI_lower),
  upper = c(NA, top_cox$CI_upper)
)

# Plot forest plot
forestplot(
  labeltext = label_text,
  mean  = hr_data$mean,
  lower = hr_data$lower,
  upper = hr_data$upper,
  zero = 1,                          # Reference line at HR = 1
  boxsize = 0.2,
  xlog = FALSE,
  col = fpColors(box = "royalblue", line = "darkblue", summary = "red"),
  title = "Univariate Forest Plot for 12 Consistently Selected Genes",
  txt_gp = fpTxtGp(                 # Font settings (no bolding)
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab  = gpar(cex = 1)
  )
)

# Multivariate Analysis (Adjusted for age and gender)

# Load required package
library(survival)

# Merge all data (clinical + expression)

colnames(clinical_df)[grepl("age|gender", colnames(clinical_df))]

# Recreate clinical_surv
clinical_surv <- data.frame(
  sample_id = rownames(clinical_df),
  time = clinical_df$days_to_death,
  status = ifelse(clinical_df$vital_status == "Alive", 0, 1)
)

# Filter to valid survival rows
clinical_surv <- clinical_surv[!is.na(clinical_surv$time), ]

# Check it worked
head(clinical_surv)
nrow(clinical_surv)

# Prepare clinical covariates with sample_id column
clinical_covars <- clinical_df[, c("age_at_diagnosis", "gender")]
clinical_covars$sample_id <- rownames(clinical_df)

# Merge: survival + clinical covariates + expression data
merged_data <- Reduce(function(x, y) merge(x, y, by = "sample_id"),
                      list(clinical_surv, clinical_covars, expr_top))

# Confirm it worked
nrow(merged_data)
head(merged_data)

# Multivariate Cox Regression for 12 Genes
# Create a results table
multi_results <- data.frame()

# Loop through 12 genes (excluding MIR155)
for (gene in setdiff(common_genes, "MIR155")) {
  
  gene_expr <- merged_data[[gene]]
  
  # Skip if missing/constant
  if (all(is.na(gene_expr)) || length(unique(gene_expr)) < 2) next
  
  # Fit multivariate Cox model
  model <- coxph(Surv(time, status) ~ gene_expr + age_at_diagnosis + gender, data = merged_data)
  s <- summary(model)
  
  # Extract gene-specific results (first row of coefficients)
  multi_results <- rbind(multi_results, data.frame(
    Gene = gene,
    HR = s$coefficients["gene_expr", "exp(coef)"],
    CI_lower = s$conf.int["gene_expr", "lower .95"],
    CI_upper = s$conf.int["gene_expr", "upper .95"],
    p_value = s$coefficients["gene_expr", "Pr(>|z|)"]
  ))
}

# Sort by p-value
multi_results <- multi_results[order(multi_results$p_value), ]
print(multi_results)

# Optional: Save
write.csv(multi_results, "Multivariate_Cox_Results_Adjusted.csv", row.names = FALSE)

# Multivariate Cox Regression for Consistently Selected Genes
# Load ggplot2
library(ggplot2)

# Assuming multi_results (multivariate results) is ready
# Ensure gene names are unique for plotting
multi_results$Gene <- make.unique(as.character(multi_results$Gene))

# Then factor by p-value ordering
multi_results$Gene <- factor(multi_results$Gene, levels = rev(multi_results$Gene))

ggplot(multi_results, aes(x = HR, y = Gene)) +
  geom_point(size = 3, color = "royalblue") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "darkblue") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Multivariate Cox Regression for 12 Consistently Selected Genes",
    subtitle = "Adjusted for Age at Diagnosis and Gender",
    x = "Hazard Ratio (HR)",
    y = "Genes"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0, size = 11)
  )

# Multivariate Forest Plot: Cox Regression for Consistently Selected Genes

# Create label text matrix
label_text <- cbind(
  c("Genes", as.character(multi_results$Gene)),  # Ensure gene names display
  c("Hazard Ratio (95% CI)",
    sprintf("%.4f (%.4f–%.4f)", multi_results$HR, multi_results$CI_lower, multi_results$CI_upper)),
  c("p-value", sprintf("%.4f", multi_results$p_value))
)

# Create numeric matrix for forestplot
hr_data <- structure(list(
  mean  = c(NA, multi_results$HR),
  lower = c(NA, multi_results$CI_lower),
  upper = c(NA, multi_results$CI_upper)
), class = "data.frame")

# Generate the forestplot
forestplot(
  labeltext = label_text,
  mean = hr_data$mean,
  lower = hr_data$lower,
  upper = hr_data$upper,
  zero = 1,
  boxsize = 0.2,
  xlog = FALSE,
  col = fpColors(box = "blue", line = "darkblue", summary = "blue"),
  title = "Multivariate Forest Plot for 12 Consistently Selected Genes",
  txt_gp = fpTxtGp(
    label = gpar(fontfamily = "", cex = 0.9),
    ticks = gpar(cex = 0.9),
    xlab = gpar(cex = 1)
  )
)
