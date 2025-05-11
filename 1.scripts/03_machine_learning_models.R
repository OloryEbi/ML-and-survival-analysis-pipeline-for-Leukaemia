# Machine Learning Modelling Workflow

# Five supervised ML models (RF, XGB, SVM, LASSO and NN) were trained to classify FLT3-ITD status using RNA-seq data.
# The dataset was split into 80% training and 20% testing, with SMOTE applied to the training set to address class imbalance.

# Model Training and Feature Selection
# Each model was trained on the balanced data, top 20 predictive genes were selected and models were retrained using these features to enhance performance and interpretability.


set.seed(42)
library(randomForest)
library(caret)
library(pROC)
library(PRROC)
library(corrplot)

# Align test data features to SMOTE data
common_features <- intersect(colnames(train_smote)[-ncol(train_smote)], colnames(test_x))
train_x_smote <- train_smote[, common_features]
test_x_aligned <- test_x[, common_features]

# Validate alignment
if (!identical(colnames(train_x_smote), colnames(test_x_aligned))) {
  stop("Error: Training and test sets don't match after alignment!")
} else {
  cat("✅ Training and test sets aligned successfully.\n")
}

# Train Random Forest using SMOTE-balanced data
rf_model <- randomForest(train_x_smote, train_smote$Group, ntree = 500, importance = TRUE)

# Predict on the aligned test set
rf_predictions <- predict(rf_model, test_x_aligned)
rf_probs <- predict(rf_model, test_x_aligned, type = "prob")[, "Pos"]

# Feature Importance
varImpPlot(rf_model)

rf_importance <- importance(rf_model)
rf_importance_df <- data.frame(
  Gene = rownames(rf_importance),
  Importance = rf_importance[, 1]
)
rf_importance_df <- rf_importance_df[order(-rf_importance_df$Importance), ]

# Plot Top 20 Important Genes
ggplot(head(rf_importance_df, 20), aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Random Forest Top 20 Important Genes", x = "Genes", y = "Importance Score") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  )

# Subset top 20 genes
top20_genes_rf <- rf_importance_df$Gene[1:20]
train_x_top20_rf <- train_x_smote[, top20_genes_rf]
test_x_top20_rf  <- test_x_aligned[, top20_genes_rf]

# Retrain RF on top 20 genes
set.seed(42)
rf_model_top20 <- randomForest(
  train_x_top20_rf, 
  train_smote$Group, 
  ntree = 500,
  importance = TRUE
)

# Predict using the top 20 genes
rf_predictions_top20 <- predict(rf_model_top20, test_x_top20_rf)
rf_probs_top20 <- predict(rf_model_top20, test_x_top20_rf, type = "prob")[, "Pos"]


# Confusion Matrix
rf_cm_top20 <- confusionMatrix(rf_predictions_top20, test_y)
print(rf_cm_top20)

# Plot Confusion Matrix Function
plot_confusion_matrix <- function(pred, actual, model_name){
  cm <- confusionMatrix(pred, actual)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(data=cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=paste("Confusion Matrix -", model_name),
         x="Actual", y="Predicted") +
    theme_minimal()
}

# Plot Confusion Matrix
plot_confusion_matrix(rf_predictions_top20, test_y, "Random Forest (Top 20 Genes)")

# ROC Curve
roc_rf_top20 <- roc(as.numeric(test_y == "Pos"), rf_probs_top20)

# Plot ROC Curve
plot(roc_rf_top20, col = "blue", main = "ROC Curve - RF (Top 20 Genes)", legacy.axes = TRUE)
abline(a=0, b=1, lty=2, col="gray")

# Calculate and print AUC
auc_rf_top20 <- auc(roc_rf_top20)
cat("Random Forest (Top 20 Genes) AUC:", round(auc_rf_top20, 3), "\n")

# Enhanced ROC plot
roc_df_rf_top20 <- data.frame(
  Specificity = 1 - roc_rf_top20$specificities, 
  Sensitivity = roc_rf_top20$sensitivities
)

ggplot(roc_df_rf_top20, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("ROC Curve - RF (Top 20 Genes, AUC =", round(auc_rf_top20, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Precision-Recall Curve
pr_rf_top20 <- pr.curve(scores.class0 = rf_probs_top20, 
                        weights.class0 = as.numeric(test_y == "Pos"), 
                        curve = TRUE)

# Plot PR Curve
plot(pr_rf_top20, col = "blue", main = "RF Precision-Recall Curve (Top 20 Genes)")
cat("RF PR AUC (Top 20 Genes):", round(pr_rf_top20$auc.integral, 3), "\n")

# Probability Histogram
predictions_rf_top20_df <- data.frame(
  Probability = rf_probs_top20,
  True_Label = test_y
)

ggplot(predictions_rf_top20_df, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  theme_minimal() +
  labs(
    title = "RF Prediction Distribution (Top 20 Genes)",
    x = "Predicted Probability",
    y = "Count"
  )

# Correlation Heatmap
rf_corr_matrix_top20 <- cor(train_x_top20_rf)

corrplot(rf_corr_matrix_top20, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full",
         tl.col = "black", tl.cex = 0.8,
         title = "RF Correlation Matrix (Top 20 Genes)",
         mar = c(0,0,1,0))


library(xgboost)
library(caret)

# Identify Common Features between SMOTE dataset and test set
common_features <- intersect(colnames(train_smote)[-ncol(train_smote)], colnames(test_x))
train_x_smote <- train_smote[, common_features]
test_x_aligned <- test_x[, common_features]

# Ensure consistent column order
train_x_smote <- train_x_smote[, order(colnames(train_x_smote))]
test_x_aligned <- test_x_aligned[, order(colnames(test_x_aligned))]

# Verify alignment
if (!identical(colnames(train_x_smote), colnames(test_x_aligned))) {
  stop("❌ Columns mismatch between training and test sets!")
} else {
  cat("✅ Training and test sets aligned successfully.\n")
}

# Convert labels to numeric (0,1)
train_y_smote_numeric <- as.numeric(train_smote$Group) - 1
test_y_numeric <- as.numeric(test_y) - 1

# Prepare matrices for XGBoost
dtrain <- xgb.DMatrix(data = as.matrix(train_x_smote), label = train_y_smote_numeric)
dtest  <- xgb.DMatrix(data = as.matrix(test_x_aligned), label = test_y_numeric)

# Train XGBoost model
set.seed(42)
xgb_model <- xgboost(
  data = dtrain,
  nrounds = 100,
  objective = "binary:logistic",
  eval_metric = "auc",
  verbose = 0
)

# Predict probabilities on test data
xgb_probs <- predict(xgb_model, dtest)

# Predicted classes
xgb_pred_classes <- factor(ifelse(xgb_probs > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos"))

# Confusion Matrix
xgb_cm <- confusionMatrix(xgb_pred_classes, test_y)
print(xgb_cm)

# Extract Feature Importance from trained XGBoost model
xgb_importance <- xgb.importance(
  feature_names = colnames(train_x_smote),
  model = xgb_model
)

# View top features
head(xgb_importance, 10)

# Extract Top 20 Genes based on 'Gain'
top20_genes_xgb <- xgb_importance$Feature[1:20]

# Plot Top 20 Genes
ggplot(xgb_importance[1:20, ], aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "XGBoostTop 20 Important Genes",
    x = "Genes",
    y = "Importance Score"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))


# Retrain model on top 20 genes
train_x_top20_xgb <- train_x_smote[, top20_genes_xgb]
test_x_top20_xgb  <- test_x_aligned[, top20_genes_xgb]

# XGBoost model on top 20 features
set.seed(42)
dtrain_top20 <- xgb.DMatrix(data = as.matrix(train_x_top20_xgb), label = train_y_smote_numeric)
dtest_top20  <- xgb.DMatrix(data = as.matrix(test_x_top20_xgb), label = test_y_numeric)

xgb_model_top20 <- xgboost(
  data = dtrain_top20,
  nrounds = 100,
  objective = "binary:logistic",
  eval_metric = "auc",
  verbose = 0
)

# Predictions with top 20 genes
xgb_probs_top20 <- predict(xgb_model_top20, dtest_top20)
xgb_pred_classes_top20 <- factor(ifelse(xgb_probs_top20 > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos"))

# Confusion Matrix (Top 20 Genes)
xgb_cm_top20 <- confusionMatrix(xgb_pred_classes_top20, test_y)
print(xgb_cm_top20)

# Plot Confusion Matrix Function
plot_confusion_matrix <- function(pred, actual, model_name){
  cm <- confusionMatrix(pred, actual)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(data=cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=paste("Confusion Matrix -", model_name),
         x="Actual", y="Predicted") +
    theme_minimal()
}

# Plot confusion matrix (top 20)
plot_confusion_matrix(xgb_pred_classes_top20, test_y, "XGBoost (Top 20 Genes)")

# ROC Curve
roc_xgb_top20 <- roc(test_y_numeric, xgb_probs_top20)
plot(roc_xgb_top20, col = "red", main = "ROC Curve - XGBoost (Top 20 Genes)", legacy.axes = TRUE)
abline(a=0, b=1, lty=2, col="gray")

# Calculate and print AUC
auc_xgb_top20 <- auc(roc_xgb_top20)
cat("XGBoost (Top 20 Genes) AUC:", round(auc_xgb_top20, 3), "\n")

# ggplot ROC Curve
roc_df_xgb_top20 <- data.frame(
  Specificity = 1 - roc_xgb_top20$specificities, 
  Sensitivity = roc_xgb_top20$sensitivities
)

ggplot(roc_df_xgb_top20, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "red", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("ROC Curve - XGBoost (Top 20 Genes, AUC =", round(auc_xgb_top20, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Precision-Recall Curve
pr_xgb_top20 <- pr.curve(scores.class0 = xgb_probs_top20, 
                         weights.class0 = test_y_numeric, 
                         curve = TRUE)

plot(pr_xgb_top20, col = "red", main = "XGBoost PR Curve (Top 20 Genes)")
cat("XGBoost PR AUC (Top 20 Genes):", round(pr_xgb_top20$auc.integral, 3), "\n")

# Probability Histogram
predictions_xgb_top20_df <- data.frame(
  Probability = xgb_probs_top20,
  True_Label = test_y
)

ggplot(predictions_xgb_top20_df, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  theme_minimal() +
  labs(
    title = "XGBoost Prediction Probability Distribution (Top 20 Genes)",
    x = "Predicted Probability (Positive Class)",
    y = "Count"
  )

# Correlation Matrix Heatmap
xgb_corr_matrix_top20 <- cor(train_x_top20_xgb)

corrplot(xgb_corr_matrix_top20, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full",
         tl.col = "black", tl.cex = 0.8,
         title = "XGBoost Correlation Matrix (Top 20 Genes)",
         mar = c(0,0,1,0))


# For XGBoost SHAP in Python
write.csv(train_x_top20_xgb, "train_x_xgb.csv", row.names = FALSE)
write.csv(test_x_top20_xgb, "test_x_xgb.csv", row.names = FALSE)
write.csv(data.frame(Group = train_smote$Group), "train_y_rf.csv", row.names = FALSE)



library(e1071)
library(caret)

# Ensure common features between SMOTE-balanced dataset and test set
common_features <- intersect(colnames(train_smote)[-ncol(train_smote)], colnames(test_x))
train_x_smote <- train_smote[, common_features]
test_x_aligned <- test_x[, common_features]

# Ensure consistent column order
train_x_smote <- train_x_smote[, order(colnames(train_x_smote))]
test_x_aligned <- test_x_aligned[, order(colnames(test_x_aligned))]

# Verify alignment
if (!identical(colnames(train_x_smote), colnames(test_x_aligned))) {
  stop("❌ Columns mismatch!")
} else {
  cat("✅ Training and test sets aligned successfully.\n")
}

# Standardize column names
colnames(train_x_smote) <- make.names(colnames(train_x_smote))
colnames(test_x_aligned) <- make.names(colnames(test_x_aligned))

# Combine training data with SMOTE labels
train_data_smote <- data.frame(train_x_smote, Group = train_smote$Group)

# Hyperparameter tuning (Radial kernel)
set.seed(42)

library(e1071)

tuned_svm <- tune(
  e1071::svm,  # explicitly use svm from e1071
  Group ~ .,
  data = train_data_smote,
  kernel = "radial",
  ranges = list(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1)),
  probability = TRUE
)

best_svm_model <- tuned_svm$best.model
cat("\n✅ Best parameters from tuning:\n")
print(tuned_svm$best.parameters)

# Predict on test set
svm_predictions <- predict(best_svm_model, test_x_aligned, probability = TRUE)
svm_probs <- attr(svm_predictions, "probabilities")[, "Pos"]

# Confusion Matrix
svm_cm <- confusionMatrix(svm_predictions, test_y)
print(svm_cm)

# Train linear SVM for feature importance extraction
svm_linear <- svm(Group ~ ., data = train_data_smote, kernel = "linear", scale = TRUE)

# Calculate absolute weights (importance)
svm_weights <- t(svm_linear$coefs) %*% svm_linear$SV
svm_importance <- abs(svm_weights)

# Create importance data frame
svm_importance_df <- data.frame(
  Gene = colnames(train_x_smote),
  Importance = as.vector(svm_importance)
)

# Rank and select top 20
svm_importance_df <- svm_importance_df[order(-svm_importance_df$Importance), ]
svm_top20_genes <- svm_importance_df$Gene[1:20]

# Plot top 20 important genes
ggplot(head(svm_importance_df, 20), aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "green") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "SVM Top 20 Important Genes",
    x = "Genes",
    y = "Importance score"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Top 20 gene subset
svm_top20_genes <- svm_importance_df$Gene[1:20]

train_x_top20_svm <- train_x_smote[, svm_top20_genes]
test_x_top20_svm  <- test_x_aligned[, svm_top20_genes]

train_data_top20_svm <- data.frame(train_x_top20_svm, Group = train_smote$Group)

# Hyperparameter tuning on top 20 genes
set.seed(42)

tuned_svm_top20 <- tune(
  e1071::svm,  # ✅ this avoids ambiguity
  Group ~ .,
  data = train_data_top20_svm,
  kernel = "radial",
  ranges = list(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1)),
  probability = TRUE
)


best_svm_model_top20 <- tuned_svm_top20$best.model

# Predict on test set (top 20)
svm_predictions_top20 <- predict(best_svm_model_top20, test_x_top20_svm, probability = TRUE)
svm_probs_top20 <- attr(svm_predictions_top20, "probabilities")[, "Pos"]

# Confusion Matrix (top 20 genes)
svm_cm_top20 <- confusionMatrix(svm_predictions_top20, test_y)
print(svm_cm_top20)

# Confusion Matrix Plot Function
plot_confusion_matrix <- function(pred, actual, model_name){
  cm <- confusionMatrix(pred, actual)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(data=cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=paste("Confusion Matrix -", model_name),
         x="Actual", y="Predicted") +
    theme_minimal()
}

plot_confusion_matrix(svm_predictions_top20, test_y, "SVM (Top 20 Genes)")

# ROC Curve
roc_svm_top20 <- roc(as.numeric(test_y == "Pos"), svm_probs_top20)
plot(roc_svm_top20, col = "green", main = "SVM ROC Curve (Top 20 Genes)", legacy.axes = TRUE)
abline(a=0, b=1, lty=2, col="gray")

auc_svm_top20 <- auc(roc_svm_top20)
cat("SVM (Top 20 Genes) AUC:", round(auc_svm_top20, 3), "\n")

# Enhanced ROC plot
roc_df_svm_top20 <- data.frame(
  Specificity = 1 - roc_svm_top20$specificities, 
  Sensitivity = roc_svm_top20$sensitivities
)

ggplot(roc_df_svm_top20, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "green", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("ROC Curve - SVM (Top 20 Genes, AUC =", round(auc_svm_top20, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Precision-Recall Curve
pr_svm_top20 <- pr.curve(scores.class0 = svm_probs_top20, 
                         weights.class0 = as.numeric(test_y == "Pos"), 
                         curve = TRUE)

plot(pr_svm_top20, col = "green", main = "SVM Precision-Recall Curve (Top 20 Genes)")
cat("SVM PR AUC (Top 20 Genes):", round(pr_svm_top20$auc.integral, 3), "\n")

# Probability Histogram
predictions_df_top20 <- data.frame(
  Probability = svm_probs_top20,
  True_Label = test_y
)

ggplot(predictions_df_top20, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "SVM Prediction Probability Distribution (Top 20 Genes)",
       x = "Predicted Probability (Positive Class)",
       y = "Count") +
  theme_minimal()

# Correlation Matrix Heatmap
svm_corr_matrix_top20 <- cor(train_x_top20_svm)

corrplot(svm_corr_matrix_top20, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full",
         tl.col = "black", tl.cex = 0.8,
         title = "SVM Correlation Matrix (Top 20 Genes)",
         mar = c(0,0,1,0))


# Export training features
write.csv(train_x_smote, "train_x_svm.csv", row.names = FALSE)

# Export test features
write.csv(test_x_aligned, "test_x_svm.csv", row.names = FALSE)

# Convert to data frame
train_y_svm <- data.frame(Group = train_smote$Group)

# Export
write.csv(train_y_svm, "train_y_rf.csv", row.names = FALSE)  # You can reuse for other models


library(glmnet)
library(caret)

# Step 1: Align training and test features
common_features <- intersect(colnames(train_smote)[-ncol(train_smote)], colnames(test_x))
train_x_smote <- train_smote[, common_features]
test_x_aligned <- test_x[, common_features]

# Step 2: Sort columns alphabetically for consistency
train_x_smote <- train_x_smote[, order(colnames(train_x_smote))]
test_x_aligned <- test_x_aligned[, order(colnames(test_x_aligned))]

# Step 3: Check alignment
stopifnot(identical(colnames(train_x_smote), colnames(test_x_aligned)))
cat("✅ Features aligned successfully.\n")

# Step 4: Convert to matrices
train_x_lasso <- as.matrix(train_x_smote)
test_x_lasso  <- as.matrix(test_x_aligned)

# Step 5: Convert labels to numeric
train_y_lasso <- ifelse(train_smote$Group == "Pos", 1, 0)
test_y_lasso  <- ifelse(test_y == "Pos", 1, 0)

# Step 6: Train LASSO logistic regression with CV
set.seed(42)
lasso_model <- cv.glmnet(
  x = train_x_lasso,
  y = train_y_lasso,
  family = "binomial",
  alpha = 1
)
best_lambda <- lasso_model$lambda.min
cat("✅ Optimal lambda:", best_lambda, "\n")

# Step 7: Predictions & Evaluation
lasso_probs <- predict(lasso_model, s = best_lambda, newx = test_x_lasso, type = "response")
lasso_pred_classes <- factor(ifelse(lasso_probs > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos"))
lasso_cm <- confusionMatrix(lasso_pred_classes, test_y)
print(lasso_cm)

# Step 8: Feature importance (non-zero coefficients)
lasso_coef <- coef(lasso_model, s = best_lambda)
lasso_importance_df <- data.frame(
  Gene = rownames(lasso_coef),
  Importance = as.vector(lasso_coef)
)
lasso_importance_df <- subset(lasso_importance_df, Gene != "(Intercept)" & Importance != 0)
lasso_importance_df <- lasso_importance_df[order(-abs(lasso_importance_df$Importance)), ]

# Step 9: Plot Top 20 Genes
ggplot(head(lasso_importance_df, 20), aes(x = reorder(Gene, abs(Importance)), y = abs(Importance))) +
  geom_bar(stat = "identity", fill = "brown") +
  coord_flip() +
  theme_minimal() +
  labs(title = "LASSO Top 20 Important Genes",
       x = "Genes",
       y = "Importance Score") +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Step 10: Get Top 20 Genes
lasso_top20_genes <- lasso_importance_df$Gene[1:20]

# Step 11: Retrain on Top 20 Genes
train_x_top20 <- train_x_smote[, lasso_top20_genes]
test_x_top20  <- test_x_aligned[, lasso_top20_genes]

set.seed(42)
lasso_top20_model <- cv.glmnet(
  x = as.matrix(train_x_top20),
  y = train_y_lasso,
  family = "binomial",
  alpha = 1
)
best_lambda_top20 <- lasso_top20_model$lambda.min
cat("✅ Optimal lambda (top 20 genes):", best_lambda_top20, "\n")

# Step 12: Predict & Evaluate (Top 20)
lasso_top20_probs <- predict(lasso_top20_model, s = best_lambda_top20, newx = as.matrix(test_x_top20), type = "response")
lasso_top20_pred_classes <- factor(ifelse(lasso_top20_probs > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos"))


# Confusion Matrix (Top 20)
lasso_top20_cm <- confusionMatrix(lasso_top20_pred_classes, test_y)
print(lasso_top20_cm)

# Confusion Matrix Visualization
plot_confusion_matrix <- function(pred, actual, model_name){
  cm <- confusionMatrix(pred, actual)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(data=cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=paste("Confusion Matrix -", model_name),
         x="Actual", y="Predicted") +
    theme_minimal()
}

plot_confusion_matrix(lasso_top20_pred_classes, test_y, "LASSO (Top 20 Genes)")

# ROC Curve
roc_lasso_top20 <- roc(test_y_lasso, as.vector(lasso_top20_probs))
plot(roc_lasso_top20, col = "brown", main = "ROC Curve - LASSO (Top 20 Genes)", legacy.axes = TRUE)
abline(a=0, b=1, lty=2, col="grey")

auc_lasso_top20 <- auc(roc_lasso_top20)
cat("LASSO (Top 20 Genes) AUC:", round(auc_lasso_top20, 3), "\n")

# Enhanced ROC plot
roc_df_lasso_top20 <- data.frame(
  Specificity = 1 - roc_lasso_top20$specificities, 
  Sensitivity = roc_lasso_top20$sensitivities
)

ggplot(roc_df_lasso_top20, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "brown", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("ROC Curve - LASSO (Top 20 Genes, AUC =", round(auc_lasso_top20, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Precision-Recall Curve
pr_lasso_top20 <- pr.curve(
  scores.class0 = as.vector(lasso_top20_probs),
  weights.class0 = test_y_lasso,
  curve = TRUE
)

plot(pr_lasso_top20, col = "brown", main = "LASSO PR Curve (Top 20 Genes)")
cat("LASSO (Top 20 Genes) PR AUC:", round(pr_lasso_top20$auc.integral, 3), "\n")

# Probability Histogram
predictions_top20_df <- data.frame(
  Probability = as.vector(lasso_top20_probs),
  True_Label = factor(test_y, levels=c("Neg", "Pos"))
)

ggplot(predictions_top20_df, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "LASSO Prediction Distribution (Top 20 Genes)",
       x = "Predicted Probability (Positive Class)",
       y = "Count") +
  theme_minimal()

# Correlation Matrix Heatmap
lasso_top20_corr_matrix <- cor(train_x_top20)

corrplot(lasso_top20_corr_matrix, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full",
         tl.col = "black", tl.cex = 0.7,
         title = "LASSO Correlation Matrix (Top 20 Genes)",
         mar = c(0,0,1,0))


write.csv(train_x_lasso, "train_x_lasso.csv", row.names = FALSE)

write.csv(test_x_lasso, "test_x_lasso.csv", row.names = FALSE)

write.csv(data.frame(Group = train_y_lasso), "train_y_lasso.csv", row.names = FALSE)


library(nnet)
library(NeuralNetTools)
library(caret)

# Step 1: Align features
common_features <- intersect(colnames(train_smote)[-ncol(train_smote)], colnames(test_x))
train_x_nn <- train_smote[, common_features]
test_x_nn  <- test_x[, common_features]

# Step 2: Sort features alphabetically
train_x_nn <- train_x_nn[, order(colnames(train_x_nn))]
test_x_nn  <- test_x_nn[, order(colnames(test_x_nn))]

# Step 3: Ensure column alignment
stopifnot(identical(colnames(train_x_nn), colnames(test_x_nn)))
cat("✅ Neural Network: Features aligned successfully.\n")

# Step 4: Scale data
train_x_scaled <- scale(train_x_nn)
test_x_scaled  <- scale(test_x_nn, 
                        center = attr(train_x_scaled, "scaled:center"),
                        scale = attr(train_x_scaled, "scaled:scale"))

# Step 5: Convert group labels to binary
train_y_bin <- ifelse(train_smote$Group == "Pos", 1, 0)

# Step 6: Train neural network model (single output node)
set.seed(42)
nn_model <- nnet(
  x = train_x_scaled,
  y = train_y_bin,
  size = 10,
  decay = 0.01,
  maxit = 500,
  linout = FALSE,
  softmax = FALSE
)

# Step 7: Calculate feature importance using Garson's algorithm
importance_df <- garson(nn_model, bar_plot = FALSE)

# Step 8: Create importance data frame
nn_importance_df <- data.frame(
  Gene = rownames(importance_df),
  Importance = importance_df$rel_imp
)

# Step 9: Order by importance
nn_importance_df <- nn_importance_df[order(-nn_importance_df$Importance), ]

# Step 10: Plot Top 20 Important Genes
ggplot(head(nn_importance_df, 20), aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Neural Netwok Top 20 Important Genes",
    x = "Genes",
    y = "Importance Sore"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Step 11: Extract Top 20 Genes
nn_top20_genes <- nn_importance_df$Gene[1:20]

# Step 12: Subset training/test data with top 20 genes
train_x_top20 <- train_x_scaled[, nn_top20_genes]
test_x_top20  <- test_x_scaled[, nn_top20_genes]

# Step 13: Train neural network using top 20 genes (multiclass format)
set.seed(42)
nn_model_top20 <- nnet(
  x = train_x_top20,
  y = class.ind(train_smote$Group),  # Converts to 2-class matrix
  size = 10,
  decay = 0.01,
  maxit = 500,
  softmax = TRUE
)

# Step 14: Predict probabilities on test set
nn_probs_top20 <- predict(nn_model_top20, test_x_top20, type = "raw")[, "Pos"]
nn_pred_classes_top20 <- factor(ifelse(nn_probs_top20 > 0.5, "Pos", "Neg"), levels = c("Neg", "Pos"))

# Step 15: Evaluate performance
nn_cm_top20 <- confusionMatrix(nn_pred_classes_top20, test_y)
print(nn_cm_top20)


# Confusion Matrix Visualization
plot_confusion_matrix <- function(pred, actual, model_name){
  cm <- confusionMatrix(pred, actual)
  cm_df <- as.data.frame(cm$table)
  colnames(cm_df) <- c("Prediction", "Actual", "Freq")
  
  ggplot(data=cm_df, aes(x=Actual, y=Prediction, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), color="white", size=6) +
    scale_fill_gradient(low="blue", high="red") +
    labs(title=paste("Confusion Matrix -", model_name),
         x="Actual", y="Predicted") +
    theme_minimal()
}

plot_confusion_matrix(nn_pred_classes_top20, test_y, "Neural Network (Top 20 Genes)")

# ROC Curve
roc_nn_top20 <- roc(as.numeric(test_y == "Pos"), nn_probs_top20)
plot(roc_nn_top20, col = "orange", main = "ROC Curve - NN (Top 20 Genes)", legacy.axes = TRUE)
abline(a=0, b=1, lty=2, col="gray")

auc_nn_top20 <- auc(roc_nn_top20)
cat("NN (Top 20 Genes) AUC:", round(auc_nn_top20, 3), "\n")

roc_df_nn <- data.frame(
  Specificity = 1 - roc_nn_top20$specificities, 
  Sensitivity = roc_nn_top20$sensitivities
)

ggplot(roc_df_nn, aes(x = Specificity, y = Sensitivity)) +
  geom_line(color = "orange", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("ROC Curve - NN (Top 20 Genes, AUC =", round(auc_nn_top20, 3), ")"),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

# Precision-Recall Curve
pr_nn_top20 <- pr.curve(scores.class0 = nn_probs_top20, 
                        weights.class0 = as.numeric(test_y == "Pos"), 
                        curve = TRUE)

plot(pr_nn_top20, col = "orange", main = "NN Precision-Recall Curve (Top 20 Genes)")
cat("NN PR AUC (Top 20 Genes):", round(pr_nn_top20$auc.integral, 3), "\n")

# Probability Distribution
predictions_df_top20 <- data.frame(
  Probability = nn_probs_top20,
  True_Label = factor(test_y, levels = c("Neg", "Pos"))
)

ggplot(predictions_df_top20, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("red", "darkblue")) +
  labs(title = "NN Prediction Probability Distribution (Top 20 Genes)",
       x = "Predicted Probability (Positive Class)",
       y = "Count") +
  theme_minimal()

# Correlation Matrix Heatmap
nn_corr_matrix_top20 <- cor(train_x_top20)

corrplot(nn_corr_matrix_top20, method = "color", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         type = "full",
         tl.col = "black", tl.cex = 0.8,
         title = "NN Correlation Matrix (Top 20 Genes)",
         mar = c(0,0,1,0))


write.csv(train_x_scaled, "train_x_nn.csv", row.names = FALSE)

write.csv(test_x_scaled, "test_x_nn.csv", row.names = FALSE)

write.csv(data.frame(Group = train_y_binary), "train_y_nn.csv", row.names = FALSE)
