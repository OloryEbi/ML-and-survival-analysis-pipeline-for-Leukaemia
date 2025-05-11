# Biological Relevance and Functional Enrichment Analysis

# Functional enrichment was performed on the consensus gene set using DAVID, covering GO terms (BP, MF, CC), KEGG, Reactome, and DISGENET.
# Terms with FDR-adjusted p-values < 0.05 were considered significant. Results were visualised in R using dot plots and bar charts for clear interpretation.


# Load Required Libraries for Biological Relevance (DAVID Enrichment)

# Load the package
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# KEGG Pathway Enrichment
kegg_enrich <- enrichKEGG(
  gene = genes_entrez,
  organism = 'hsa',
  pvalueCutoff = 0.05
)
kegg_results <- as.data.frame(kegg_enrich)
cat("Number of significant KEGG pathways:", nrow(kegg_results), "\n")
head(kegg_results, 10)


dotplot(kegg_enrich, showCategory = 10, title = "Top KEGG Pathways") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10)  # Centered, bold, bigger
  )


library(ggplot2)

# Full GO:BP dataset from DAVID screenshot (16 terms)
go_bp_david <- data.frame(
  Term = c(
    "Embryonic skeletal system morphogenesis",
    "Anterior/posterior pattern specification",
    "Regulation of transcription by RNA polymerase II",
    "Positive regulation of transcription by RNA polymerase II",
    "Definitive hemopoiesis",
    "Mammary gland alveolus development",
    "Cellular response to lipopolysaccharide",
    "Leukocyte migration",
    "Thyroid gland development",
    "Angiogenesis",
    "Cell migration",
    "Cellular response to hormone stimulus",
    "Extracellular matrix disassembly",
    "Response to retinoic acid",
    "Response to estrogen",
    "Response to mechanical stimulus"
  ),
  Count = c(5, 5, 7, 6, 2, 2, 3, 2, 2, 3, 3, 3, 2, 2, 2, 2),
  PValue = c(2.1e-7, 5.3e-6, 1.6e-2, 1.6e-2, 2.1e-2, 2.2e-2, 3.5e-2, 3.9e-2,
             3.9e-2, 4.1e-2, 4.5e-2, 5.4e-2, 5.5e-2, 5.7e-2, 6.9e-2, 6.9e-2)
)

# Add formatted label for P-values
go_bp_david$logP <- -log10(go_bp_david$PValue)
go_bp_david$label <- format(go_bp_david$PValue, scientific = TRUE, digits = 2)

# Add leukemia-related flag
go_bp_david$Highlight <- ifelse(go_bp_david$Term %in% c(
  "Definitive hemopoiesis",
  "Leukocyte migration",
  "Angiogenesis",
  "Cell migration",
  "Extracellular matrix disassembly"
), "Leukemia-related", "Other")


ggplot(go_bp_david, aes(x = Count, y = reorder(Term, Count))) +
  geom_point(size = 5, color = "steelblue") +
  geom_text(aes(label = label), hjust = -0.2, size = 3.5) +
  labs(
    title = "GO: Biological Processes Enrichment",
    x = "Gene Count", y = "GO Term"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  xlim(0, max(go_bp_david$Count) + 2)


ggplot(go_bp_david, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "GO: Biological Processes Enrichment",
    x = "GO Term", y = expression(-log[10](P-value))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  ylim(0, max(go_bp_david$logP) + 1)


library(ggplot2)

# GO:MF DAVID data
go_mf_david <- data.frame(
  Term = c(
    "RNA pol II cis-regulatory DNA binding",
    "Transcription activator activity (RNA pol II)",
    "Transcription factor activity (RNA pol II)",
    "Sequence-specific dsDNA binding",
    "Scavenger receptor activity",
    "DNA-binding transcription factor activity"
  ),
  Count = c(7, 5, 7, 5, 2, 4),
  PValue = c(3.0e-3, 3.1e-3, 4.2e-3, 4.8e-3, 3.4e-2, 3.7e-2)
)

# Compute logP and formatted labels
go_mf_david$logP <- -log10(go_mf_david$PValue)
go_mf_david$label <- format(go_mf_david$PValue, scientific = TRUE, digits = 2)

# Add highlight label
go_mf_david$Highlight <- ifelse(go_mf_david$Term %in% c(
  "RNA pol II cis-regulatory DNA binding",
  "Transcription activator activity (RNA pol II)",
  "Transcription factor activity (RNA pol II)"
), "Leukemia-related", "Other")


ggplot(go_mf_david, aes(x = Count, y = reorder(Term, Count))) +
  geom_point(size = 5, color = "steelblue") +
  geom_text(aes(label = label), hjust = -0.2, size = 3.5) +
  labs(
    title = "GO: Molecular Function Enrichment",
    x = "Gene Count", y = "GO Term"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  xlim(0, max(go_mf_david$Count) + 2)


ggplot(go_mf_david, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "GO: Molecular Function Enrichment",
    x = "GO Term", y = expression(-log[10](P-value))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  ylim(0, max(go_mf_david$logP) + 1)


library(ggplot2)

# GO:CC DAVID data
go_cc_david <- data.frame(
  Term = c("Chromatin", 
           "Extracellular region", 
           "Secretory granule", 
           "Collagen-containing extracellular matrix", 
           "Extracellular space"),
  Count = c(7, 9, 3, 4, 6),
  PValue = c(1.9e-3, 3.3e-3, 1.1e-2, 1.1e-2, 8.6e-2)
)

# Add logP and p-value labels
go_cc_david$logP <- -log10(go_cc_david$PValue)
go_cc_david$label <- format(go_cc_david$PValue, scientific = TRUE, digits = 2)

go_cc_david$Highlight <- ifelse(go_cc_david$Term %in% c(
  "Chromatin", "Extracellular region"
), "Leukemia-related", "Other")


ggplot(go_cc_david, aes(x = Count, y = reorder(Term, Count))) +
  geom_point(size = 5, color = "steelblue") +
  geom_text(aes(label = label), hjust = -0.2, size = 3.5) +
  labs(
    title = "GO: Cellular Component Enrichment",
    x = "Gene Count", y = "GO Term"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  xlim(0, max(go_cc_david$Count) + 2)


ggplot(go_cc_david, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "GO: Cellular Component Enrichment",
    x = "GO Term", y = expression(-log[10](P-value))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  ylim(0, max(go_cc_david$logP) + 1)


# KEGG Pathway

gene_symbols <- rownames(deg_filtered)  # or your top gene list

gene_df <- bitr(gene_symbols,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Extract Entrez IDs
genes_entrez <- gene_df$ENTREZID

kegg_enrich <- enrichKEGG(
  gene = genes_entrez,
  organism = 'hsa',
  pvalueCutoff = 0.05
)

head(kegg_enrich)

dotplot(kegg_enrich, 
        showCategory = 10, 
        title = "KEGG Pathway Enrichment (Dotplot)",
        font.size = 11) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )


# Reactome enrichment data (manually extracted)
reactome_df <- data.frame(
  Term = c(
    "Regulation of IGF transport and uptake",
    "Metabolism of Angiotensinogen to Angiotensins",
    "Infection with Mycobacterium tuberculosis",
    "Activation of Matrix Metalloproteinases",
    "Bacterial Infection Pathways",
    "Signaling by Interleukins",
    "Peptide hormone metabolism"
  ),
  Count = c(3, 2, 2, 2, 2, 3, 2),
  PValue = c(7.6e-3, 1.8e-2, 2.9e-2, 3.5e-2, 7.4e-2, 8.6e-2, 9.1e-2)
)

# Add transformed p-values and formatted labels
reactome_df$logP <- -log10(reactome_df$PValue)
reactome_df$label <- format(reactome_df$PValue, scientific = TRUE, digits = 2)

# Add pathway labels
reactome_df$Highlight <- ifelse(reactome_df$Term %in% c(
  "Regulation of IGF transport and uptake",
  "Signaling by Interleukins",
  "Activation of Matrix Metalloproteinases"
), "Leukemia-related", "Other")


ggplot(reactome_df, aes(x = Count, y = reorder(Term, Count))) +
  geom_point(size = 5, color = "darkgreen") +
  geom_text(aes(label = label), hjust = -0.2, size = 3.5) +
  labs(
    title = "Reactome Pathway Enrichment",
    x = "Gene Count", y = "Reactome Pathway"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  xlim(0, max(reactome_df$Count) + 2)


ggplot(reactome_df, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "Reactome Pathway Enrichment",
    x = "Reactome Pathway", y = expression(-log[10](P-value))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  ylim(0, max(reactome_df$logP) + 1)


library(ggplot2)

# Full DISGENET data
disgenet_df <- data.frame(
  Disease = c(
    "Liver Cirrhosis", "Fibrosis, Liver", "Acute Myeloid Leukemia, M1",
    "Acute Myeloid Leukemia (AML-M2)", "Narcolepsy-Cataplexy Syndrome", "Narcolepsy",
    "Leukemia, Myelocytic, Acute", "Neoplasm Metastasis", "Diaphragmatic Hernia",
    "Alcoholic Intoxication, Chronic", "Insulin Resistance", "Insulin Sensitivity",
    "Dermatitis, Allergic Contact", "Alcohol Abuse"
  ),
  Count = c(3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2),
  PValue = c(
    1.1e-2, 1.1e-2, 1.6e-2, 1.7e-2, 2.4e-2, 2.6e-2, 2.7e-2,
    3.1e-2, 3.2e-2, 3.6e-2, 7.9e-2, 9.1e-2, 9.2e-2, 9.9e-2
  )
)

# Add -log10(p-value) and formatted labels
disgenet_df$logP <- -log10(disgenet_df$PValue)
disgenet_df$label <- format(disgenet_df$PValue, scientific = TRUE, digits = 2)

# Add a label column to flag leukemia terms
disgenet_df$LeukemiaFlag <- ifelse(grepl("Leukemia", disgenet_df$Disease, ignore.case = TRUE), "Leukemia-related", "Other")

ggplot(disgenet_df, aes(x = Count, y = reorder(Disease, Count))) +
  geom_point(size = 5, color = "firebrick") +
  geom_text(aes(label = label), hjust = -0.2, size = 3.5) +
  labs(
    title = "DISGENET Gene-Disease Associations",
    x = "Gene Count", y = "Disease"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  xlim(0, max(disgenet_df$Count) + 2)


ggplot(disgenet_df, aes(x = reorder(Disease, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(
    title = "DISGENET Gene-Disease Associations",
    x = "Disease", y = expression(-log[10](P-value))
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 12)
  ) +
  ylim(0, max(disgenet_df$logP) + 1)

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Select Top 3 Terms by P-value for Each GO Category
top3_bp <- go_bp_david %>% arrange(PValue) %>% slice_head(n = 3)
top3_mf <- go_mf_david %>% arrange(PValue) %>% slice_head(n = 3)
top3_cc <- go_cc_david %>% arrange(PValue) %>% slice_head(n = 3)

# 2. Calculate shared maximum for logP for axis alignment
max_logP_top3 <- max(top3_bp$logP, top3_mf$logP, top3_cc$logP)

# 3. Define a shared minimal theme
custom_theme <- theme_bw(base_size = 9) +
  theme(
    plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 2),
    plot.title = element_text(size = 10, face = "bold", hjust = 1),
    legend.position = "none"
  )

# 4. GO:Biological Process (BP) Plot
plot_bp_top3 <- ggplot(top3_bp, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "orange", width = 0.6) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.2) +
  coord_flip() +
  labs(title = "A. Top 3 GO: Biological Processes",
       x = "GO Term", y = expression(-log[10](P-value))) +
  custom_theme +
  ylim(0, max_logP_top3 + 0.5)

# 5. GO:Molecular Function (MF) Plot
plot_mf_top3 <- ggplot(top3_mf, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "darkgreen", width = 0.6) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.2) +
  coord_flip() +
  labs(title = "B. Top 3 GO: Molecular Functions",
       x = "GO Term", y = expression(-log[10](P-value))) +
  custom_theme +
  ylim(0, max_logP_top3 + 0.5)

# 6. GO:Cellular Component (CC) Plot
plot_cc_top3 <- ggplot(top3_cc, aes(x = reorder(Term, logP), y = logP)) +
  geom_bar(stat = "identity", fill = "firebrick", width = 0.6) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.2) +
  coord_flip() +
  labs(title = "C. Top 3 GO: Cellular Components",
       x = "GO Term", y = expression(-log[10](P-value))) +
  custom_theme +
  ylim(0, max_logP_top3 + 0.5)

# 7. Combine All Three Plots Vertically
combined_top3_plot <- plot_bp_top3 / plot_mf_top3 / plot_cc_top3 + 
  plot_layout(guides = "collect")

# 8. Display Combined Plot
combined_top3_plot
