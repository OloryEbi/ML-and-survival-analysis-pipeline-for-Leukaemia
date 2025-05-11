# LASSO-Cox Regression and Risk Stratification

# A LASSO-Cox model was used to identify prognostic genes and compute a risk score for each patient.
# Based on this score, patients were stratified into high- and low-risk groups, with survival outcomes compared using Kaplan–Meier analysis.


# LASSO COX REG.
# Load necessary libraries
library(glmnet)
library(survival)

# Step 1: Filter out samples with non-positive survival time
merged_data <- merged_data[merged_data$time > 0, ]

# Step 2: Create survival object
y <- Surv(time = merged_data$time, event = merged_data$status)

# Step 3: Define selected genes (if not already available)
# Replace this with your actual list
selected_genes <- c("CIBAR1P1", "ENPP2", "MIR155", "LGALS3BP", "HOXB3", 
                    "TRGC2", "HOXB6", "NKX2-3", "ADGRG1", "SNRPN", 
                    "SOCS2", "SPNS3", "C3orf80", "MRC1", "HOXB2", 
                    "PBX3", "IQCJ-SCHIP1", "LAPTM4B", "HOXA5", 
                    "PDE4B", "MMP2", "SPINK2", "HOXA9", "IGFBP2", 
                    "CTSG", "TRH", "CPA3")

missing_genes <- setdiff(selected_genes, colnames(merged_data))
cat("❌ Missing genes:\n", missing_genes, "\n")

selected_genes_present <- intersect(selected_genes, colnames(merged_data))


# Step 4: Subset expression matrix for selected genes
X <- as.matrix(merged_data[, selected_genes_present])


# Step 5: Fit LASSO Cox model
set.seed(123)
cv_lasso <- cv.glmnet(
  x = X,
  y = y,
  family = "cox",
  alpha = 1,
  nfolds = 10
)

# Step 6: Get best lambda and selected genes
best_lambda <- cv_lasso$lambda.min
lasso_coef <- coef(cv_lasso, s = best_lambda)

selected_genes_final <- rownames(lasso_coef)[as.vector(lasso_coef) != 0]
selected_genes_final <- selected_genes_final[selected_genes_final != "(Intercept)"]

# Print result
cat("✅ Final LASSO-selected genes:", selected_genes_final, "\n")


# CV Error vs. log(lambda)
plot(cv_lasso)
abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)
title("Cross-Validation Curve (LASSO Cox)", line = 2.5)


# Refit full glmnet model without CV for full coefficient path
fit_lasso <- glmnet(
  x = X,
  y = y,
  family = "cox",
  alpha = 1
)

# Plot coefficient paths vs log(lambda)
plot(fit_lasso, xvar = "lambda", label = TRUE)
abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)
title("LASSO Coefficient Paths", line = 2.5)

# Final Multivariate Cox Model (Post-LASSO)
# `selected_genes` should already be available from previous step
print(selected_genes)

# 1. Extract expression values of final selected genes
X_selected <- merged_data[, selected_genes_final]

# 2. Combine with survival data
cox_input <- cbind(
  time = merged_data$time,
  status = merged_data$status,
  X_selected
)

# 3. Fit Cox model
y_final <- Surv(cox_input$time, cox_input$status)
final_cox_model <- coxph(y_final ~ ., data = as.data.frame(X_selected))

# 4. View summary
summary(final_cox_model)

cox_summary <- summary(final_cox_model)
cox_df <- data.frame(
  Gene = rownames(cox_summary$coefficients),
  HR = exp(cox_summary$coefficients[, "coef"]),
  CI_lower = exp(cox_summary$conf.int[, "lower .95"]),
  CI_upper = exp(cox_summary$conf.int[, "upper .95"]),
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

# View or save
print(cox_df)
write.csv(cox_df, "Final_Multivariate_Cox_Model.csv", row.names = FALSE)

# Forest Plot of Final LASSO Cox Model
library(forestplot)
# Create label text
label_text <- cbind(
  c("Genes", cox_df$Gene),
  c("Hazard Ratio (95% CI)", sprintf("%.2f (%.2f–%.2f)", cox_df$HR, cox_df$CI_lower, cox_df$CI_upper)),
  c("p-value", sprintf("%.4f", cox_df$p_value))
)

# Numeric values for plot
hr_data <- data.frame(
  mean = c(NA, cox_df$HR),
  lower = c(NA, cox_df$CI_lower),
  upper = c(NA, cox_df$CI_upper)
)

# Plot
forestplot(
  labeltext = label_text,
  mean = hr_data$mean,
  lower = hr_data$lower,
  upper = hr_data$upper,
  zero = 1,
  boxsize = 0.2,
  xlog = FALSE,
  col = fpColors(box = "royalblue", line = "darkblue", summary = "blue"),
  title = "Final LASSO Cox Model: Forest Plot",
  txt_gp = fpTxtGp(label = gpar(cex = 0.9), ticks = gpar(cex = 0.8), xlab = gpar(cex = 1))
)

# Build Risk Score from LASSO Genes
# Use fitted model coefficients
lasso_risk_score <- predict(final_cox_model, newdata = as.data.frame(X_selected), type = "lp")

# Add to merged_data
merged_data$risk_score <- lasso_risk_score

# Split into high and low risk using median
merged_data$risk_group <- ifelse(merged_data$risk_score >= median(merged_data$risk_score), "High", "Low")
merged_data$risk_group <- factor(merged_data$risk_group, levels = c("Low", "High"))

fit_risk <- survfit(Surv(time, status) ~ risk_group, data = merged_data)

# 2. Load the library
library(survminer)

# 3. Now run your KM plot
ggsurvplot(
  fit_risk,
  data = merged_data,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("forestgreen", "red"),
  title = "Kaplan-Meier Survival Based on LASSO Risk Score",
  legend.title = "Risk Group",
  legend.labs = c("Low", "High"),
  risk.table.height = 0.25
)
