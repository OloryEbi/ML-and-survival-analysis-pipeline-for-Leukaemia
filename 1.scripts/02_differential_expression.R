# Differential Gene Expression Analysis

# Differential expression between FLT3-ITD positive and negative AML groups was assessed using the limma package in R.
# Genes with an absolute log₂ fold change > 1 and FDR-adjusted p-value < 0.05 were considered significant
# Results were visualised with volcano plots using ggplot2, with top genes annotated by fold change.


# Differential Expression Analysis (limma)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Linear modeling with limma
fit <- lmFit(expr_sub, design)
contrast_matrix <- makeContrasts(Pos - Neg, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract all DEGs
deg_results <- topTable(fit2, adjust = "fdr", number = Inf)

# Filter significant DEGs (Adj.P.Val < 0.05, logFC > |1|)
deg_filtered <- deg_results %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Quick summary
cat("Total significant DEGs:", nrow(deg_filtered), "\n")
head(deg_filtered)

# Volcano Plot
deg_results$Significance <- ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Up",
                                   ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Down", "Not"))

library(ggrepel)

top_labeled <- deg_results %>%
  dplyr::filter(Significance != "Not") %>%
  dplyr::top_n(10, wt = abs(logFC))

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(data = top_labeled, aes(label = rownames(top_labeled)), size = 3) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.title = element_blank())


top_labeled <- deg_results %>%
  dplyr::filter(Significance != "Not") %>%
  dplyr::top_n(10, wt = abs(logFC))

rownames(top_labeled)

# Prepare Expression Data for Machine Learning
selected_genes <- deg_filtered %>% rownames()
expr_df <- as.data.frame(t(expr_sub[selected_genes, ]))

# Attach group labels
expr_df$Group <- group

# Check dimensions and structure
dim(expr_df)
table(expr_df$Group)

# Train/Test Split

set.seed(42)
train_index <- createDataPartition(expr_df$Group, p = 0.8, list = FALSE)
train_data <- expr_df[train_index, ]
test_data <- expr_df[-train_index, ]

train_x <- train_data[, -ncol(train_data)]
train_y <- train_data$Group
test_x <- test_data[, -ncol(test_data)]
test_y <- test_data$Group

# Handle Class Imbalance (SMOTE)

train_y_numeric <- as.numeric(train_y) - 1
smote_result <- SMOTE(X = train_x, target = train_y_numeric, K = 5)
train_smote <- smote_result$data
colnames(train_smote)[ncol(train_smote)] <- "Group"
train_smote$Group <- factor(train_smote$Group, levels = c(0, 1), labels = levels(train_y))

# Verify SMOTE results
table(train_smote$Group)
