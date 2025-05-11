# Model Evaluation and Visualisation

# Model performance was assessed on the test set using metrics such as accuracy, sensitivity, specificity, PPV, NPV, balanced accuracy, and AUC.
# Visualisations included confusion matrices, ROC and PR curves, prediction probability histograms and correlation heatmaps for top-ranked genes.

# Cross-Referencing of Gene Signatures Across ML Models
# The top 20 genes from each ML model were compared to identify shared features.
# A five-way Venn diagram highlighted overlapping genes, defining a consensus predictive signature.
# Genes with low or absent expression were excluded from further analysis to ensure biological relevance.


# COMPARE MODEL RESULT
# Align Predictions and Outcomes

test_y_numeric <- ifelse(test_y == "Pos", 1, 0)

# Compare Confusion Matrices
library(patchwork)

# Define a function to plot confusion matrix
plot_confusion_matrix <- function(cm, title) {
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=title, x="Actual", y="Predicted") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10)
    )
}

# Create individual confusion matrix plots
cm_rf_plot <- plot_confusion_matrix(rf_cm_top20, "Random Forest")
cm_xgb_plot <- plot_confusion_matrix(xgb_cm_top20, "XGBoost")
cm_svm_plot <- plot_confusion_matrix(svm_cm_top20, "SVM")
cm_lasso_plot <- plot_confusion_matrix(lasso_top20_cm, "LASSO")
cm_nn_plot <- plot_confusion_matrix(nn_cm_top20, "Neural Network")

# Combine using patchwork
# Keep aspect ratio consistent for better alignment
confusion_plot <- (cm_rf_plot | cm_xgb_plot) / (cm_svm_plot | cm_lasso_plot | cm_nn_plot) +
  plot_annotation(
    title = "Confusion Matrix Comparison (Top 20 Genes)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  )

# Display combined plot
confusion_plot

# Summary dataframe for top 20 genes results (Bar graph)

# Extract Accuracy, Sensitivity, and Specificity directly from the confusion matrix for the top 20 genes.
summary_results_top20 <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "LASSO", "Neural Network"),
  Accuracy = c(rf_cm_top20$overall['Accuracy'],
               xgb_cm_top20$overall['Accuracy'],
               svm_cm_top20$overall['Accuracy'],
               lasso_top20_cm$overall['Accuracy'],
               nn_cm_top20$overall['Accuracy']),
  
  Sensitivity = c(rf_cm_top20$byClass['Sensitivity'],
                  xgb_cm_top20$byClass['Sensitivity'],
                  svm_cm_top20$byClass['Sensitivity'],
                  lasso_top20_cm$byClass['Sensitivity'],
                  nn_cm_top20$byClass['Sensitivity']),
  
  Specificity = c(rf_cm_top20$byClass['Specificity'],
                  xgb_cm_top20$byClass['Specificity'],
                  svm_cm_top20$byClass['Specificity'],
                  lasso_top20_cm$byClass['Specificity'],
                  nn_cm_top20$byClass['Specificity'])
)

# Print table clearly
print(summary_results_top20)

library(tidyr)

# Step 1: Reshape to long format
results_long_top20 <- pivot_longer(summary_results_top20, 
                                   cols = c(Accuracy, Sensitivity, Specificity),
                                   names_to = "Metric", 
                                   values_to = "Value")

# Step 2: Plot
ggplot(results_long_top20, aes(x = reorder(Model, Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Value, 4)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = c("Accuracy" = "skyblue", "Sensitivity" = "green", "Specificity" = "orange")) +
  ylim(0, 1) +
  labs(
    title = "Model Performance Comparison (Top 20 Genes)",
    x = "Machine Learning Models",
    y = "Performance Metric Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
  )

# Load required libraries
library(reshape2)

# Create a performance data frame
model_perf <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "LASSO", "Neural Network"),
  Accuracy = c(0.913, 0.8913, 0.9022, 0.8804, 0.8804),
  Balanced_Accuracy = c(0.8901, 0.8752, 0.8701, 0.8803, 0.8803),
  AUC_ROC = c(0.968, 0.968, 0.919, 0.949, 0.946),
  AUC_PR = c(0.895, 0.923, 0.832, 0.854, 0.863)
)

# Melt the data into long format
model_perf_long <- melt(model_perf, id.vars = "Model", 
                        variable.name = "Metric", value.name = "Score")

# Plot
ggplot(model_perf_long, aes(x = Model, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  coord_cartesian(ylim = c(0.8, 1)) +
  labs(
    title = "Model Performance Comparison (Top 20 Genes)",
    x = "Models",
    y = "Scores"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0, face = "bold", size = 11)  # Center and bold
  ) +
  scale_fill_brewer(palette = "Set1")

# Load required libraries
library(reshape2)

# Create performance data frame
model_perf <- data.frame(
  Model = c("Random Forest", "XGBoost", "SVM", "LASSO", "Neural Network"),
  Accuracy = c(0.913, 0.8913, 0.9022, 0.8804, 0.8804),
  Balanced_Accuracy = c(0.8901, 0.8752, 0.8701, 0.8803, 0.8803),
  AUC_ROC = c(0.968, 0.968, 0.919, 0.949, 0.946),
  AUC_PR = c(0.895, 0.923, 0.832, 0.854, 0.863)
)

# Melt the data into long format
model_perf_long <- melt(model_perf, id.vars = "Model", 
                        variable.name = "Metric", value.name = "Score")

# Plot with score labels
ggplot(model_perf_long, aes(x = Model, y = Score, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(Score, 3)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  coord_cartesian(ylim = c(0.8, 1.02)) +
  labs(
    title = "Model Performance Comparison (Top 20 Genes)",
    x = "Model",
    y = "Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # Center and bold the title
  ) +
  scale_fill_brewer(palette = "Set1")


# Overlay and Compare ROC Curves

library(pROC)
library(ggplot2)

# Extract AUC values
auc_rf <- round(auc(roc_rf_top20), 3)
auc_xgb <- round(auc(roc_xgb_top20), 3)
auc_svm <- round(auc(roc_svm_top20), 3)
auc_lasso <- round(auc(roc_lasso_top20), 3)
auc_nn <- round(auc(roc_nn_top20), 3)

# Create ROC data for each model
roc_combined_df <- rbind(
  data.frame(Specificity = 1 - roc_rf_top20$specificities, Sensitivity = roc_rf_top20$sensitivities, 
             Model = paste0("Random Forest (AUC = ", auc_rf, ")")),
  data.frame(Specificity = 1 - roc_xgb_top20$specificities, Sensitivity = roc_xgb_top20$sensitivities, 
             Model = paste0("XGBoost (AUC = ", auc_xgb, ")")),
  data.frame(Specificity = 1 - roc_svm_top20$specificities, Sensitivity = roc_svm_top20$sensitivities, 
             Model = paste0("SVM (AUC = ", auc_svm, ")")),
  data.frame(Specificity = 1 - roc_lasso_top20$specificities, Sensitivity = roc_lasso_top20$sensitivities, 
             Model = paste0("LASSO (AUC = ", auc_lasso, ")")),
  data.frame(Specificity = 1 - roc_nn_top20$specificities, Sensitivity = roc_nn_top20$sensitivities, 
             Model = paste0("Neural Network (AUC = ", auc_nn, ")"))
)

# Create a color palette (no need for matching names)
model_colors <- c("brown", "orange", "blue", "green", "red")

# Plot ROC curves
ggplot(roc_combined_df, aes(x = Specificity, y = Sensitivity, color = Model, linetype = Model)) +
  geom_line(size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = model_colors) +
  scale_linetype_manual(values = rep("solid", 5)) +
  labs(
    title = "ROC Curve Comparison (Top 20 Genes)",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.65, 0.25),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11)  # ← Center, bold, and size
  )


# Overlay and Compare Precision-Recall Curves

library(PRROC)

# Extract curves
pr_rf_curve    <- as.data.frame(pr_rf_top20$curve)
pr_xgb_curve   <- as.data.frame(pr_xgb_top20$curve)
pr_svm_curve   <- as.data.frame(pr_svm_top20$curve)
pr_lasso_curve <- as.data.frame(pr_lasso_top20$curve)
pr_nn_curve    <- as.data.frame(pr_nn_top20$curve)

# Extract AUC values
auc_pr_rf    <- round(pr_rf_top20$auc.integral, 3)
auc_pr_xgb   <- round(pr_xgb_top20$auc.integral, 3)
auc_pr_svm   <- round(pr_svm_top20$auc.integral, 3)
auc_pr_lasso <- round(pr_lasso_top20$auc.integral, 3)
auc_pr_nn    <- round(pr_nn_top20$auc.integral, 3)

# Combine into a single data frame
pr_combined_df <- rbind(
  data.frame(Recall = pr_rf_curve$V1, Precision = pr_rf_curve$V2, 
             Model = paste0("Random Forest (AUC = ", auc_pr_rf, ")")),
  data.frame(Recall = pr_xgb_curve$V1, Precision = pr_xgb_curve$V2, 
             Model = paste0("XGBoost (AUC = ", auc_pr_xgb, ")")),
  data.frame(Recall = pr_svm_curve$V1, Precision = pr_svm_curve$V2, 
             Model = paste0("SVM (AUC = ", auc_pr_svm, ")")),
  data.frame(Recall = pr_lasso_curve$V1, Precision = pr_lasso_curve$V2, 
             Model = paste0("LASSO (AUC = ", auc_pr_lasso, ")")),
  data.frame(Recall = pr_nn_curve$V1, Precision = pr_nn_curve$V2, 
             Model = paste0("Neural Network (AUC = ", auc_pr_nn, ")"))
)

# Get unique model labels for color mapping
model_labels <- unique(pr_combined_df$Model)

# Set custom colors
model_colors <- c("blue", "red", "green", "brown", "orange")
names(model_colors) <- model_labels

# Plot PR curves
ggplot(pr_combined_df, aes(x = Recall, y = Precision, color = Model, linetype = Model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = model_colors) +
  scale_linetype_manual(values = rep("solid", length(model_labels))) +
  labs(
    title = "Precision-Recall Curve Comparison (Top 20 Genes)",
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.65, 0.25),
    legend.title = element_blank(),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)  # ← Center, bold, and size
  )


# Overlay and Compare (Histogram) Prediction Probability Distribution

library(gridExtra)
library(grid)

# Create a dataframe for the top 20 genes predictions and true labels
predictions_df_top20 <- data.frame(
  True_Label = as.factor(test_y),
  RF = rf_probs_top20,
  XGB = xgb_probs_top20,
  SVM = svm_probs_top20,
  LASSO = as.vector(lasso_top20_probs),
  NN = as.vector(nn_probs_top20)
)

# Random Forest Prediction Distribution
rf_plot <- ggplot(predictions_df_top20, aes(x = RF, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "Random Forest") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11))

# XGBoost
xgb_plot <- ggplot(predictions_df_top20, aes(x = XGB, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "XGBoost") +  # <- fixed quote
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11))

# SVM
svm_plot <- ggplot(predictions_df_top20, aes(x = SVM, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "SVM") +  # <- fixed quote
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11))

# LASSO
lasso_plot <- ggplot(predictions_df_top20, aes(x = LASSO, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "LASSO") +  # <- fixed quote
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11))

# Neural Network
nn_plot <- ggplot(predictions_df_top20, aes(x = NN, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "Neural Network") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 11))

# Combine all Histogram plots together using grid.arrange()
grid.arrange(
  rf_plot, xgb_plot, svm_plot, lasso_plot, nn_plot,
  ncol = 2,
  top = textGrob(
    "Prediction Probability Distribution by Model (Top 20 Genes)",
    gp = gpar(fontsize = 12, fontface = "bold"),
    hjust = 0.5
  )
)


# Compare Correlation Heatmaps for all models

library(corrplot)
library(reshape2)
library(dplyr)

# Compute correlation matrices for each model's top 20 genes
rf_corr_matrix_top20 <- cor(train_x_top20_rf)
xgb_corr_matrix_top20 <- cor(train_x_top20_xgb)
svm_corr_matrix_top20 <- cor(train_x_top20_svm)
lasso_corr_matrix_top20 <- cor(train_x_top20)
nn_corr_matrix_top20 <- cor(train_x_top20)

# Convert correlation matrices to long format using melt()
rf_corr <- melt(rf_corr_matrix_top20) %>% mutate(Model = "Random Forest")
xgb_corr <- melt(xgb_corr_matrix_top20) %>% mutate(Model = "XGBoost")
svm_corr <- melt(svm_corr_matrix_top20) %>% mutate(Model = "SVM")
lasso_corr <- melt(lasso_corr_matrix_top20) %>% mutate(Model = "LASSO")
nn_corr <- melt(nn_corr_matrix_top20) %>% mutate(Model = "Neural Network")

# Combine all correlation matrices into one dataframe
combined_corr <- rbind(rf_corr, xgb_corr, svm_corr, lasso_corr, nn_corr)

# Average correlation values across models
combined_corr_avg <- combined_corr %>%
  group_by(Var1, Var2) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

# Plot the combined correlation matrix
ggplot(combined_corr_avg, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  
  # Define consistent color scale
  scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    mid = "white", 
    midpoint = 0, 
    limits = c(-1, 1)
  ) +
  
  # Title and labels
  labs(
    title = "Combined Correlation Heatmap for all models (Top 20 Genes)",
    x = "Genes",
    y = "Genes"
  ) +
  
  # Theme settings
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.title = element_text(size = 10)
  )


library(pROC)
library(patchwork)

# Function to plot individual ROC curve
plot_roc_curve <- function(roc_obj, model_name, auc_value, color) {
  df <- data.frame(
    Specificity = 1 - roc_obj$specificities,
    Sensitivity = roc_obj$sensitivities
  )
  
  ggplot(df, aes(x = Specificity, y = Sensitivity)) +
    geom_line(color = color, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(
      title = paste0(model_name, " (AUC = ", auc_value, ")"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text = element_text(size = 10)
    )
}


# Calculate AUCs
auc_rf <- round(auc(roc_rf_top20), 3)
auc_nn <- round(auc(roc_nn_top20), 3)
auc_xgb <- round(auc(roc_xgb_top20), 3)
auc_svm <- round(auc(roc_svm_top20), 3)
auc_lasso <- round(auc(roc_lasso_top20), 3)


# Generate plots
roc_rf_plot <- plot_roc_curve(roc_rf_top20, "Random Forest", auc_rf, "blue")
roc_nn_plot <- plot_roc_curve(roc_nn_top20, "Neural Network", auc_nn, "orange")
roc_xgb_plot <- plot_roc_curve(roc_xgb_top20, "XGBoost", auc_xgb, "red")
roc_svm_plot <- plot_roc_curve(roc_svm_top20, "SVM", auc_svm, "green")
roc_lasso_plot <- plot_roc_curve(roc_lasso_top20, "LASSO", auc_lasso, "brown")


# Combine using patchwork (same style as confusion matrices)
roc_combined_plot <- (roc_rf_plot | roc_nn_plot) / (roc_xgb_plot | roc_svm_plot | roc_lasso_plot) +
  plot_annotation(
    title = "ROC Curve Comparison (Top 20 Genes)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  )

# Display combined plot
roc_combined_plot

library(PRROC)
library(patchwork)

# Extract AUC values
auc_rf    <- round(pr_rf_top20$auc.integral, 3)
auc_nn    <- round(pr_nn_top20$auc.integral, 3)
auc_xgb   <- round(pr_xgb_top20$auc.integral, 3)
auc_svm   <- round(pr_svm_top20$auc.integral, 3)
auc_lasso <- round(pr_lasso_top20$auc.integral, 3)


# Function to plot individual PR curve
plot_pr_curve <- function(pr_curve_df, model_name, auc_value, color) {
  ggplot(pr_curve_df, aes(x = V1, y = V2)) +
    geom_line(color = color, size = 1.2) +
    labs(
      title = paste0(model_name, " (AUC = ", auc_value, ")"),
      x = "Recall",
      y = "Precision"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text = element_text(size = 10)
    )
}

# Generate plots with actual AUCs
pr_rf_plot    <- plot_pr_curve(pr_rf_curve, "Random Forest", auc_rf, "blue")
pr_nn_plot    <- plot_pr_curve(pr_nn_curve, "Neural Network", auc_nn, "orange")
pr_xgb_plot   <- plot_pr_curve(pr_xgb_curve, "XGBoost", auc_xgb, "red")
pr_svm_plot   <- plot_pr_curve(pr_svm_curve, "SVM", auc_svm, "green")
pr_lasso_plot <- plot_pr_curve(pr_lasso_curve, "LASSO", auc_lasso, "brown")


# Combine using patchwork
pr_combined_plot <- (pr_rf_plot | pr_nn_plot) / (pr_xgb_plot | pr_svm_plot | pr_lasso_plot) +
  plot_annotation(
    title = "Precision-Recall Curve Comparison (Top 20 Genes)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  )

# Display
pr_combined_plot

# Cross-referencing the Top Genes Across Models
# Load necessary libraries
library(VennDiagram)
library(grid)
library(dplyr)

# Example: Identify common genes across models
common_genes <- Reduce(intersect, list(
  top20_genes_rf,
  top20_genes_xgb,
  svm_top20_genes,
  lasso_top20_genes,
  nn_top20_genes
))

cat("Consistently selected genes across models:\n")
print(common_genes)

# Create Venn Diagram (but don't draw yet)
venn <- venn.diagram(
  x = list(
    RandomForest   = top20_genes_rf,
    SVM            = svm_top20_genes,
    XGBoost        = top20_genes_xgb,
    NeuralNetwork  = nn_top20_genes,
    LASSO          = lasso_top20_genes
  ),
  filename = NULL,
  fill = c("skyblue", "green", "red", "grey", "yellow"),
  alpha = 0.7,
  cex = 1.2,
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.dist = 0.05,
  lwd = 1.5,
  margin = 0.05,
  main = "Venn Diagram of Top 20 Genes Across Models",
  main.cex = 1,
  main.fontface = "bold",
  main.pos = c(0.5, 1.05)
)

# Clear the page and set custom viewport to scale + shift the diagram up
grid.newpage()
pushViewport(viewport(y = 0.52, height = 0.9))  # Shift up + shrink
grid.draw(venn)
popViewport()

# Add note below the diagram
grid.text(
  "Note: MIR155 was excluded from survival analysis due to zero expression across all samples.",
  x = 0.5, y = 0.05,
  gp = gpar(fontsize = 11, fontface = "italic", col = "black")
)

# Identify non-intersecting genes
non_common_rf <- setdiff(top20_genes_rf, common_genes)
non_common_xgb <- setdiff(top20_genes_xgb, common_genes)
non_common_svm <- setdiff(svm_top20_genes, common_genes)
non_common_lasso <- setdiff(lasso_top20_genes, common_genes)
non_common_nn <- setdiff(nn_top20_genes, common_genes)

# Create table with equal length using NA padding
max_length <- max(
  length(non_common_rf),
  length(non_common_xgb),
  length(non_common_svm),
  length(non_common_lasso),
  length(non_common_nn),
  length(common_genes)
)

gene_table <- data.frame(
  RF = c(non_common_rf, rep(NA, max_length - length(non_common_rf))),
  XGBoost = c(non_common_xgb, rep(NA, max_length - length(non_common_xgb))),
  SVM = c(non_common_svm, rep(NA, max_length - length(non_common_svm))),
  LASSO = c(non_common_lasso, rep(NA, max_length - length(non_common_lasso))),
  NN = c(non_common_nn, rep(NA, max_length - length(non_common_nn))),
  Intersection = c(common_genes, rep(NA, max_length - length(common_genes)))
)

# Display the table
print(gene_table)

# Export to CSV (Optional)
write.csv(gene_table, "Gene_Comparison_Table.csv", row.names = FALSE)

# Convert to long format using tidyr (Optional)
gene_table_long <- gene_table %>%
  pivot_longer(cols = everything(), names_to = "Model", values_to = "Gene") %>%
  drop_na()

# View long-format table
print(gene_table_long)

# Identify intersection of top genes across models
common_genes <- Reduce(intersect, list(top20_genes_rf, top20_genes_xgb, svm_top20_genes, lasso_top20_genes, nn_top20_genes))

# Print the list of common genes
print(common_genes)

# Display the list in a cleaner format
cat("Consistently selected genes across models:\n", paste(common_genes, collapse = ", "), "\n")

# Optional: Save to a text or CSV file
write.csv(data.frame(Gene = common_genes), "Common_Genes.csv", row.names = FALSE)

ml_genes <- data.frame(
  RF = c("HOXB3", "SPNS3", "SNRPN", "HOXB2", "HOXB6", "NKX2-3", "LAPTM4B", rep(NA, 6)),
  XGBoost = c("HOXB3", "SPNS3", "HOXB6", "IQCJ-SCHIP1", "PBX3", "TRH", "SNRPN", rep(NA, 6)),
  SVM = c("HOXB3", "HOXB6", "IQCJ.SCHIP1", "MMP2", "HOXA5", "C3orf80", "IGFBP2", rep(NA, 6)),
  LASSO = c("HOXB6", "IQCJ-SCHIP1", "SPNS3", "MMP2", "TRH", "IGFBP2", "HOXB3", rep(NA, 6)),
  NN = c("TRH", "SNRPN", "NKX2-3", "IGFBP2", "IQCJ-SCHIP1", "C3orf80", "LAPTM4B", rep(NA, 6)),
  Intersection = c("SOCS2", "MRC1", "CIBAR1P1", "ENPP2", "MIR155", "LGALS3BP", 
                   "TRGC2", "SPINK2", "CTSG", "PDE4B", "ADGRG1", "CPA3", "HOXA9"),
  stringsAsFactors = FALSE
)

# Get unique gene list (all models + intersection)
all_genes <- unique(na.omit(unlist(ml_genes)))

# Preview
print(all_genes)

# Save as text file for DAVID
write.table(all_genes, "top20_ML_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
