########################## PART I ##############################
# --------------load necessary packages--------------------
library(biomaRt)
library(openxlsx)
library(readxl)
library(tidyverse)
library(writexl)

# --------------load datasets--------------------
expression_raw <- read_tsv ("GDC TCGA GBM gene expression dataset.tsv")
phenotype_raw <- read_tsv("GDC TCGA GBM phenotype dataset.tsv")

# --------------filter phenotype dataset--------------------
phenotype_filtered <- phenotype_raw %>%
  filter(
    !is.na(days_to_death.demographic),
    tumor_descriptor.samples == "Primary",
    tissue_type.samples == "Tumor",
    prior_treatment.diagnoses == "No"
  )

# --------------extract matched samples--------------------
# Get unique sample IDs from phenotype dataset - column 1
phenotype_samples <- unique(phenotype_filtered[[1]])

# Get unique sample IDs from gene expression dataset - row 1
expression_samples <- colnames(expression_raw)

# Find overlapping samples between both datasets
matched_columns <- intersect(phenotype_samples, expression_samples)

# Matched phenotype dataset
matched_phenotype <- phenotype_filtered[phenotype_filtered[[1]] %in% matched_columns, ]

# Matched gene expression dataset
gene_id_col <- colnames(expression_raw)[1]
matched_expression <- expression_raw[, c(gene_id_col, matched_columns), drop = FALSE]

# ----------------count no. of SS, IS, LS samples----------------------
# Count the sum of short survival (SS) samples
count <- sum(matched_phenotype$days_to_death.demographic < 365)
print(count) # 61

# Count the sum of intermediate survival (IS1) samples
count <- sum(matched_phenotype$days_to_death.demographic >= 365 & matched_phenotype$days_to_death.demographic <= 730)
print(count) # 43

# Count the sum of long survival (LS) samples
count <- sum(matched_phenotype$days_to_death.demographic > 730)
print(count) # 19

# ----------------assign survival groups---------------------
matched_phenotype$survival_group <- ifelse(matched_phenotype$days_to_death.demographic < 365, "SS",
                                           ifelse(matched_phenotype$days_to_death.demographic <= 730, "IS",
                                                  ifelse(matched_phenotype$days_to_death.demographic > 730, "LS", NA)))

# Remove rows with missing survival group
matched_phenotype <- matched_phenotype %>% filter(!is.na(survival_group))

# Convert survival group to a factor - Important for DESeq2 
matched_phenotype$survival_group <- factor(matched_phenotype$survival_group)

# --------------connect to Ensembl------------------
# Connect to Ensembl homo sapien genes database
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "asia")

# Define the requested attributes
attributes <- c("ensembl_gene_id", "hgnc_symbol")

# Extract the protein coding genes
protein_coding_ref <- getBM(attributes = attributes,
                            filters = "biotype",
                            values = "protein_coding",
                            mart = ensembl)

# To get the version of biomaRt I used
packageVersion("biomaRt")     # version = 2.66.2

# To get the version of ensembl I used 
# Release 116, June 2026, GRCh38.p14

# --------------filter gene expression dataset for protein coding genes ONLY------------------
# Normalize Ensembl IDs: Remove version numbers
cleaned_expr_data <- matched_expression %>%
  mutate(Ensembl_ID_base = str_remove(Ensembl_ID, "\\.\\d+$")) %>%
  mutate(across(where(is.numeric), ~ round(2^.x - 1)))

# Match gene expression data
gene_expression_data <- cleaned_expr_data %>%
  inner_join(protein_coding_ref, by = c("Ensembl_ID_base" = "ensembl_gene_id")) %>%
  
  # Handle Ensembl ID duplicates - keeping the one with the highest expression 
  group_by(Ensembl_ID_base) %>%
  slice_max(order_by = rowSums(across(where(is.numeric))), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Handle External Gene Name duplicates - keeping the one with the highest expression
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  group_by(hgnc_symbol) %>%
  slice_max(order_by = rowSums(across(where(is.numeric))), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Clean up
  dplyr::select(-Ensembl_ID_base)

# Result: 19411 protein coding genes

# --------------reorder columns--------------------
gene_expression_data <- gene_expression_data %>% dplyr::relocate(hgnc_symbol, .after = Ensembl_ID)

# --------------download final results--------------------
write_xlsx(gene_expression_data, "Final_Gene_Expression.xlsx") # 19411 pc genes
write_xlsx(matched_phenotype,    "Final_Phenotype.xlsx") # 123 samples 


########################## PART II ##############################
# --------------load necessary packages--------------------
library(airway)
library(apeglm)
library(ashr)
library(cowplot)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)
library(tibble)

# --------------load datasets--------------------
counts_data <- read_excel("Final_Gene_Expression.xlsx")
sample_info <- read_excel("Final_Phenotype.xlsx")

# --------------check the TCGA-06-0211 vial pair--------------------
a <- "TCGA-06-0211-01A"; b <- "TCGA-06-0211-01B"

count_mat <- as.matrix(counts_data[, !(names(counts_data) %in% c("Ensembl_ID", "hgnc_symbol"))])
rownames(count_mat) <- counts_data$hgnc_symbol
expr_log <- log2(count_mat + 1)

x <- expr_log[, a]; y <- expr_log[, b]

identical(x, y)                                   # FALSE
cor(x, y); cor(x, y, method = "spearman")         # x = 0.9795503 | y = 0.9827059
cc <- cor(expr_log); diag(cc) <- NA               # sample -01A is closest to -01B
sort(cc[a, ], decreasing = TRUE)[1:5]             

# MA plot
M <- x - y; A <- (x + y) / 2                      
plot(A, M, pch = ".", xlab = "Mean log2(count+1)", ylab = "01A - 01B")
abline(h = c(-1, 0, 1), col = c("grey", "red", "grey"))
colSums(count_mat[, c(a, b)])                     # A:46852894 | B:51319317 
colSums(count_mat[, c(a, b)] > 0)                 # A:17096 | B: 17242 

# --------------keep one sample per patient: drop -01A--------------------
counts_data <- counts_data %>% dplyr::select(-all_of(a))
sample_info <- sample_info %>% filter(sample != a)

# --------------build matrix and phenotype table--------------------
counts_matrix <- counts_data %>%
  dplyr::select(-Ensembl_ID) %>%
  column_to_rownames(var = "hgnc_symbol")

sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$sample

# Enforce the same sample order, then check
counts_matrix <- counts_matrix[, rownames(sample_info)]
stopifnot(ncol(counts_matrix) == 122,
          nrow(sample_info) == 122,
          !(a %in% colnames(counts_matrix)),
          all(colnames(counts_matrix) == rownames(sample_info)))

# --------------variance-stabilised expression--------------------
dds_122 <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_matrix)),
                                  colData   = sample_info,
                                  design    = ~ 1)
vst_mat <- assay(vst(dds_122, blind = TRUE))
dim(vst_mat)                                      # 19411 x 122

# --------------save for downstream scripts--------------------
saveRDS(counts_matrix, "counts_matrix_122.rds")
saveRDS(sample_info,   "sample_info_122.rds")
saveRDS(vst_mat,       "vst_mat_122.rds")

# --------------construct a DESeq2 data set--------------------
# Convert from charactor to factor 
sample_info$survival_group <- factor(sample_info$survival_group)

counts_matrix <- round(as.matrix(counts_matrix))
mode(counts_matrix) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData   = sample_info,
                              design    = ~ survival_group)

dds$survival_group <- relevel(dds$survival_group, ref = "LS")
levels(dds$survival_group)   

# --------------PCA quality control--------------------
vsd <- vst(dds, blind = TRUE)

pca_dat <- plotPCA(vsd, intgroup = "survival_group", returnData = TRUE)
percentVar <- round(100 * attr(pca_dat, "percentVar"))
percentVar

pca_dat$survival_group <- factor(pca_dat$survival_group, levels = c("SS","IS","LS"))

ggplot(pca_dat, aes(PC1, PC2, colour = survival_group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# ------------------PCA assessment--------------------
# one-way ANOVA
summary(aov(PC1 ~ survival_group, data = pca_dat))   # (F(2,119) = 0.046, p = 0.955)
summary(aov(PC2 ~ survival_group, data = pca_dat))   # (F(2,119) = 4.826, p = 0.00965**)

# Tukey test
TukeyHSD(aov(PC2 ~ survival_group, data = pca_dat)) # 95% CI 

# IS-SS: 6.5295 | p = 0.022 | 0.7847388 to 12.274306 
# LS-SS: 7.3880 | p = 0.057 | -0.1807797 to 14.956786
# LS-IS: 0.8585 | p = 0.964 | -7.0619671 to 8.778928

# So, already these numbers tell me that there is a difference between SS and LS, 
# but it mainly lies between the SS-IS margin then the IS-LS margin

# --------------run deseq2--------------------
dds <- DESeq(dds)

# --------------save deseq2 results--------------------
saveRDS(dds, "dds_gbm_survival_groups.rds")

# --------------load deseq2 results--------------------
dds <- readRDS("dds_gbm_survival_groups.rds")

# --------------deseq2 results--------------------
# 1. SS vs LS
res_SS_vs_LS <- results(dds, contrast = c("survival_group", "SS", "LS"), alpha = 0.05)
summary(res_SS_vs_LS)   # 208 at padj<0.05, pre-LFC filter (117 upregulated & 91 downregulated)

# 2. IS vs LS
res_IS_vs_LS <- results(dds, contrast = c("survival_group", "IS", "LS"), alpha = 0.05)
summary(res_IS_vs_LS)   # 54 at padj<0.05, pre-LFC filter (22 upregulated & 32 downregulated)

# 3. SS vs IS — IS as the reference
res_SS_vs_IS <- results(dds, contrast = c("survival_group", "SS", "IS"), alpha = 0.05)
summary(res_SS_vs_IS)   # 267 at padj<0.05, pre-LFC filter (225 upregulated & 42 downregulated)

# --------------post-hoc deseq2 results --------------------
sig_SS_vs_LS <- rownames(res_SS_vs_LS)[which(res_SS_vs_LS$padj < 0.05 & abs(res_SS_vs_LS$log2FoldChange) > 1)]
sig_IS_vs_LS <- rownames(res_IS_vs_LS)[which(res_IS_vs_LS$padj < 0.05 & abs(res_IS_vs_LS$log2FoldChange) > 1)]
sig_SS_vs_IS <- rownames(res_SS_vs_IS)[which(res_SS_vs_IS$padj < 0.05 & abs(res_SS_vs_IS$log2FoldChange) > 1)]

# --------------save post-hoc deseq results-----------------
sig <- function(r) rownames(r)[which(r$padj < 0.05 & abs(r$log2FoldChange) > 1)]
deg_union <- unique(c(sig(res_SS_vs_LS), sig(res_IS_vs_LS), sig(res_SS_vs_IS)))
length(deg_union)                         # 291 unique sDEGs

# -----save the 291 unique gene names-----
saveRDS(deg_union, "deg_union_122.rds") 

# Full DESeq2 results (adjusted p-value + log2FC) for downstream scripts
saveRDS(list(SS_vs_LS = res_SS_vs_LS,
             IS_vs_LS = res_IS_vs_LS,
             SS_vs_IS = res_SS_vs_IS),
        "res_122.rds")

# --------------number of sDEGs--------------------
length(sig_SS_vs_LS) # 160 DEGs
length(sig_IS_vs_LS) # 47  DEGs
length(sig_SS_vs_IS) # 144 DEGs 

all_sig <- unique(c(sig_SS_vs_LS, sig_IS_vs_LS, sig_SS_vs_IS))
length(all_sig) # 291 DEGs 

# -----------upregulated and downregulated sDEGs--------------
up_SS_vs_LS   <- rownames(res_SS_vs_LS)[which(res_SS_vs_LS$padj < 0.05 & res_SS_vs_LS$log2FoldChange > 1)]
down_SS_vs_LS <- rownames(res_SS_vs_LS)[which(res_SS_vs_LS$padj < 0.05 & res_SS_vs_LS$log2FoldChange < -1)]

up_IS_vs_LS   <- rownames(res_IS_vs_LS)[which(res_IS_vs_LS$padj < 0.05 & res_IS_vs_LS$log2FoldChange > 1)]
down_IS_vs_LS <- rownames(res_IS_vs_LS)[which(res_IS_vs_LS$padj < 0.05 & res_IS_vs_LS$log2FoldChange < -1)]

up_SS_vs_IS   <- rownames(res_SS_vs_IS)[which(res_SS_vs_IS$padj < 0.05 & res_SS_vs_IS$log2FoldChange > 1)]
down_SS_vs_IS <- rownames(res_SS_vs_IS)[which(res_SS_vs_IS$padj < 0.05 & res_SS_vs_IS$log2FoldChange < -1)]

# --------------number of upregulated and downregulated sDEGs  --------------------
length(up_SS_vs_LS);   length(down_SS_vs_LS)   # 113 up,  47 down
length(up_IS_vs_LS);   length(down_IS_vs_LS)   #  21 up,  26 down
length(up_SS_vs_IS);   length(down_SS_vs_IS)   # 130 up,  14 down

# --------------venn diagram--------------------
venn_list <- list("SS vs LS" = sig_SS_vs_LS,
                  "IS vs LS" = sig_IS_vs_LS,
                  "SS vs IS" = sig_SS_vs_IS)

venn_plot <- ggVennDiagram(venn_list,
                           label               = "both",   # count + percentage
                           label_percent_digit = 0,
                           label_alpha         = 0,
                           label_size          = 5,
                           set_size            = 5.5) +
  ggplot2::scale_fill_gradient(low = "white", high = "grey85", guide = "none")

# Bold set names and region labels
for (i in seq_along(venn_plot$layers)) {
  if (class(venn_plot$layers[[i]]$geom)[1] %in% c("GeomText", "GeomLabel")) {
    venn_plot$layers[[i]]$aes_params$fontface <- "bold"
  }
}

venn_plot

# --------------volcano plot--------------------
vol <- function(res, ttl, legend = "none") {
  EnhancedVolcano::EnhancedVolcano(res,
                                   lab            = rownames(res),
                                   x              = "log2FoldChange",
                                   y              = "padj",
                                   title          = ttl,
                                   subtitle       = NULL,
                                   caption        = NULL,
                                   xlab           = bquote(~log[2] ~ "fold change"),
                                   ylab           = bquote(~-log[10] ~ "adjusted" ~ italic(P)),
                                   pCutoff        = 0.05,
                                   FCcutoff       = 1,
                                   drawConnectors = TRUE,
                                   legendPosition = legend)
}

p1 <- vol(res_SS_vs_LS, "SS vs LS")
p2 <- vol(res_IS_vs_LS, "IS vs LS")
p3 <- vol(res_SS_vs_IS, "SS vs IS")

volcano_fig <- cowplot::plot_grid(p1, p2, p3, nrow = 1, labels = c("A", "B", "C"))
volcano_fig

ggplot2::ggsave("Figure_volcano.png", volcano_fig,
                width = 18, height = 7, dpi = 300, bg = "white")

# --------------volcano plot x venn diagram figure--------------------
# 2 x 2 figure 
fig2 <- cowplot::plot_grid(p1, p2, p3, venn_plot,
                           nrow = 2, ncol = 2,
                           labels = c("A", "B", "C", "D"),
                           label_size = 16)
fig2

ggplot2::ggsave("Figure2_DEGs.png", fig2,
                width = 14, height = 12, dpi = 300, bg = "white")

# --------------sDEGs info ----------------
deg_table <- function(res, contrast_name) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(abs(log2FoldChange) > 1,  padj < 0.05) %>%
    mutate(contrast  = contrast_name,
           direction = if_else(log2FoldChange > 0, "Up", "Down"))
}

deg_all <- bind_rows(
  deg_table(res_SS_vs_LS, "SS vs LS"),
  deg_table(res_IS_vs_LS, "IS vs LS"),
  deg_table(res_SS_vs_IS, "SS vs IS")
) %>%
  dplyr::select(gene, contrast, direction, baseMean,
                log2FoldChange, lfcSE, pvalue, padj) %>%
  arrange(contrast, direction, desc(abs(log2FoldChange)))

View(deg_all)

# --------------download sDEGs info ----------------
write.csv(deg_all, "results/DEGs_all_contrasts.csv", row.names = FALSE)

########################## PART III ##############################
# --------------load necessary packages--------------------
library(clusterProfiler)
library(dplyr)
library(enrichplot)  
library(ggh4x)
library(ggplot2)
library(msigdbr)
library(org.Hs.eg.db)
library(ReactomePA)
library(tidyr)
library(tidytext)

###################### GO: ORA #######################
# -------A single function that does GO for all three contrasts-----
run_go <- function(res) {
  genes_up   <- rownames(subset(res, padj < 0.05 & log2FoldChange >  1))
  genes_down <- rownames(subset(res, padj < 0.05 & log2FoldChange < -1))
  universe   <- rownames(subset(res, !is.na(padj)))
  
  go_one <- function(genes, ont) {
    if (length(genes) == 0) return(NULL)
    ego <- enrichGO(gene = genes, universe = universe,
                    OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                    ont = ont, pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                    minGSSize = 5, maxGSSize = 5000, readable = TRUE)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
    simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }
  
  list(
    n_up = length(genes_up), n_down = length(genes_down), n_universe = length(universe),
    up   = list(BP = go_one(genes_up, "BP"),   CC = go_one(genes_up, "CC"),   MF = go_one(genes_up, "MF")),
    down = list(BP = go_one(genes_down, "BP"), CC = go_one(genes_down, "CC"), MF = go_one(genes_down, "MF"))
  )
}

# -------run all three contrasts and save------
go_results <- list(
  SS_vs_LS = run_go(res_SS_vs_LS),
  IS_vs_LS = run_go(res_IS_vs_LS),
  SS_vs_IS = run_go(res_SS_vs_IS)
)

# To save GO results
saveRDS(go_results, "go_results.rds")

# To read GO results
go_results <- readRDS("go_results.rds")

# -------number of simplified GO terms per contrast, direction and ontology-------
term_counts <- do.call(rbind, lapply(names(go_results), function(ct) {
  r <- go_results[[ct]]
  do.call(rbind, lapply(c("up", "down"), function(dir) {
    data.frame(
      contrast  = ct,
      direction = dir,
      n_genes   = if (dir == "up") r$n_up else r$n_down,
      BP = if (is.null(r[[dir]]$BP)) 0 else nrow(as.data.frame(r[[dir]]$BP)),
      CC = if (is.null(r[[dir]]$CC)) 0 else nrow(as.data.frame(r[[dir]]$CC)),
      MF = if (is.null(r[[dir]]$MF)) 0 else nrow(as.data.frame(r[[dir]]$MF))
    )
  }))
}))

term_counts

# ------------plot: helper functions------------
prep_go <- function(obj, ont_label, n = 5) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df <- df %>% slice_min(p.adjust, n = n, with_ties = FALSE)
  gr <- strsplit(df$GeneRatio, "/")
  df$GeneRatio_num <- vapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  df$ontology <- ont_label
  df[, c("Description", "GeneRatio_num", "Count", "p.adjust", "ontology")]
}

go_dot <- function(x, ttl) {
  df <- bind_rows(prep_go(x$BP, "BP"), prep_go(x$CC, "CC"), prep_go(x$MF, "MF"))
  if (nrow(df) == 0) return(NULL)
  df$ontology <- factor(df$ontology, levels = c("BP", "CC", "MF"))
  ggplot(df, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num),
                 size = Count, colour = p.adjust)) +
    geom_point() +
    facet_grid2(ontology ~ ., scales = "free_y", space = "free_y",
                strip = strip_themed(
                  background_y = elem_list_rect(fill = c("#1B9E77", "#D95F02", "#7570B3")),
                  text_y = elem_list_text(colour = "white", face = "bold"))) +
    scale_colour_gradient(low = "red", high = "blue") +
    scale_y_discrete(labels = scales::label_wrap(35)) +   # wrap long GO terms
    labs(x = "Gene ratio", y = NULL, title = ttl, colour = "Adjusted P", size = "Gene count") +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.y  = element_text(face = "bold", size = 11),
          plot.title   = element_text(face = "bold", size = 14))
}

# ---------figure of the two contrasts with enough enriched terms-------
fig_go <- cowplot::plot_grid(
  go_dot(go_results$SS_vs_LS$up, "SS vs LS: upregulated in SS"),
  go_dot(go_results$SS_vs_IS$up, "SS vs IS: upregulated in SS"),
  ncol = 1, align = "v", axis = "lr",
  labels = c("A", "B"), label_size = 16)

fig_go
ggsave("Figure3A_GO.png", fig_go, width = 12, height = 14, dpi = 300, bg = "white")

# --------small results displayed in the console: methods + results---------
as.data.frame(go_results$IS_vs_LS$up$BP)[, c("Description", "Count", "p.adjust")]
as.data.frame(go_results$SS_vs_LS$down$CC)[, c("Description", "Count", "p.adjust")]

# ---- supplementary table: GO (all simplified significant terms) ----
grp <- list(
  SS_vs_LS = c("SS", "LS"),
  SS_vs_IS = c("SS", "IS"),
  IS_vs_LS = c("IS", "LS")
)

go_table <- bind_rows(lapply(names(go_results), function(ct) {
  bind_rows(lapply(c("up", "down"), function(dir) {
    bind_rows(lapply(c("BP", "CC", "MF"), function(ont) {
      o <- go_results[[ct]][[dir]][[ont]]
      if (is.null(o) || nrow(as.data.frame(o)) == 0) return(NULL)
      as.data.frame(o) %>%
        transmute(Comparison     = gsub("_", " ", ct),
                  `Higher in`    = if (dir == "up") grp[[ct]][1] else grp[[ct]][2],
                  Ontology       = ont,
                  `GO ID`        = ID,
                  Term           = Description,
                  `Gene ratio`   = GeneRatio,
                  `Background ratio` = BgRatio,
                  `Gene count`   = Count,
                  `P value`      = signif(pvalue, 3),
                  `Adjusted P`   = signif(p.adjust, 3),
                  `q value`      = signif(qvalue, 3),
                  Genes          = gsub("/", ", ", geneID))
    }))
  }))
}))

go_notes <- data.frame(Column = names(go_table), Description = c(
  "Survival group comparison (first group vs reference)",
  "Group in which the input genes were more highly expressed",
  "GO ontology: BP (biological process), CC (cellular component), MF (molecular function)",
  "Gene Ontology identifier",
  "GO term name",
  "Input genes annotated to the term / total input genes",
  "Background genes annotated to the term / total background genes",
  "Number of input genes annotated to the term",
  "Unadjusted P value (hypergeometric test)",
  "Benjamini-Hochberg adjusted P value",
  "q value",
  "Input genes annotated to the term"))

write_xlsx(list(GO_enrichment = go_table, Notes = go_notes),
           "Supplementary_Table_GO.xlsx")

##################### HALLMARK GSEA #######################
# ---------Hallmark gene sets----------
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# ----------function that runs Hallmark GSEA----------------
run_gsea <- function(res) {
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$stat), ]
  ranks  <- sort(setNames(res_df$stat, rownames(res_df)), decreasing = TRUE)
  
  set.seed(42)
  GSEA(geneList = ranks, TERM2GENE = hallmark,
       pvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 5000,
       eps = 0, seed = TRUE)
}

# ----------run all three contrasts once and save results--------------
if (!file.exists("gsea_results.rds")) {
  gsea_results <- list(
    SS_vs_LS = run_gsea(res_SS_vs_LS),
    IS_vs_LS = run_gsea(res_IS_vs_LS),
    SS_vs_IS = run_gsea(res_SS_vs_IS)
  )
  saveRDS(gsea_results, "gsea_results.rds")
}
gsea_results <- readRDS("gsea_results.rds")

# ----------number of significant gene sets per contrast and direction-----------
sapply(gsea_results, function(x) {
  d <- as.data.frame(x)
  c(positive_NES = sum(d$NES > 0), negative_NES = sum(d$NES < 0))
})

# ------------------bar plot-----------------------
gsea_bar <- function(x, ttl, pos_lab, neg_lab, n = 10) {
  df <- as.data.frame(x)
  if (nrow(df) == 0) return(NULL)
  df$label     <- gsub("_", " ", sub("^HALLMARK_", "", df$ID))
  df$direction <- factor(ifelse(df$NES > 0, pos_lab, neg_lab), levels = c(neg_lab, pos_lab))
  df <- df %>% group_by(direction) %>%
    slice_max(abs(NES), n = n, with_ties = FALSE) %>% ungroup()
  
  ggplot(df, aes(x = reorder(label, NES), y = NES, fill = direction)) +
    geom_col() + coord_flip() +
    geom_hline(yintercept = 0, colour = "grey40") +
    scale_fill_discrete(drop = FALSE) +
    labs(x = NULL, y = "Normalised enrichment score (NES)", title = ttl) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(face = "bold", size = 10),
          plot.title  = element_text(face = "bold", size = 14))
}

# ---------------publishable GSEA figure: A-C-----------------
fig_gsea <- plot_grid(
  gsea_bar(gsea_results$SS_vs_LS, "SS vs LS", "Higher in SS", "Higher in LS"),
  gsea_bar(gsea_results$IS_vs_LS, "IS vs LS", "Higher in IS", "Higher in LS"),
  gsea_bar(gsea_results$SS_vs_IS, "SS vs IS", "Higher in SS", "Higher in IS"),
  ncol = 1, align = "v", axis = "lr",
  labels = c("A", "B", "C"), label_size = 16)

fig_gsea
ggsave("Figure3B_GSEA.png", fig_gsea, width = 10, height = 16, dpi = 300, bg = "white")

# ----------------combining all enrichment diagrams: GO & GSEA------------
# Left column: GO (A, B)
go_col <- cowplot::plot_grid(
  go_dot(go_results$SS_vs_LS$up, "SS vs LS: higher in SS"),
  go_dot(go_results$SS_vs_IS$up, "SS vs IS: higher in SS"),
  ncol = 1, align = "v", axis = "lr",
  labels = c("A", "B"), label_size = 16)

# Right column: GSEA (C, D, E) — same contrast order as GO, heights matched to bar counts
gsea_col <- cowplot::plot_grid(
  gsea_bar(gsea_results$SS_vs_LS, "SS vs LS", "Higher in SS", "Higher in LS"),
  gsea_bar(gsea_results$SS_vs_IS, "SS vs IS", "Higher in SS", "Higher in IS"),
  gsea_bar(gsea_results$IS_vs_LS, "IS vs LS", "Higher in IS", "Higher in LS"),
  ncol = 1, align = "v", axis = "lr",
  labels = c("C", "D", "E"), label_size = 16,
  rel_heights = c(15, 18, 6) + 6)   # bars + room for title, axis and legend

# Side by side
fig_enrich <- cowplot::plot_grid(go_col, gsea_col, ncol = 2, rel_widths = c(1.2, 1))

fig_enrich
ggsave("Figure3_enrichment.png", fig_enrich, width = 16, height = 13, dpi = 300, bg = "white")

# ----------------saving the full list of enrichment terms (GO & GSEA)------------
go_results   <- readRDS("go_results.rds")
gsea_results <- readRDS("gsea_results.rds")

# ---- supplementary table: Hallmark GSEA (all significant gene sets) ----
gsea_table <- bind_rows(lapply(names(gsea_results), function(ct) {
  d <- as.data.frame(gsea_results[[ct]])
  if (nrow(d) == 0) return(NULL)
  d %>%
    transmute(Comparison    = gsub("_", " ", ct),
              `Gene set`    = gsub("_", " ", sub("^HALLMARK_", "", ID)),
              `Enriched in` = ifelse(NES > 0, grp[[ct]][1], grp[[ct]][2]),
              `Set size`    = setSize,
              NES           = round(NES, 2),
              `P value`     = signif(pvalue, 3),
              `Adjusted P`  = signif(p.adjust, 3),
              `Leading-edge genes` = gsub("/", ", ", core_enrichment)) %>%
    arrange(desc(NES))
}))

gsea_notes <- data.frame(Column = names(gsea_table), Description = c(
  "Survival group comparison (first group vs reference)",
  "MSigDB Hallmark gene set",
  "Group in which the gene set was enriched (positive NES = first group)",
  "Number of genes in the gene set present in the ranked list",
  "Normalised enrichment score",
  "Unadjusted P value",
  "Benjamini-Hochberg adjusted P value",
  "Core genes driving the enrichment signal"))

write_xlsx(list(Hallmark_GSEA = gsea_table, Notes = gsea_notes),
           "Supplementary_Table_GSEA.xlsx")

######################## KEGG ##########################
# --------------function that runs KEGG (up and down)---------------
run_kegg <- function(res) {
  genes_up   <- rownames(subset(res, padj < 0.05 & log2FoldChange >  1))
  genes_down <- rownames(subset(res, padj < 0.05 & log2FoldChange < -1))
  universe   <- rownames(subset(res, !is.na(padj)))
  
  to_entrez <- function(x) {
    if (length(x) == 0) return(character(0))
    suppressWarnings(bitr(x, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
  }
  up_e <- to_entrez(genes_up); down_e <- to_entrez(genes_down); uni_e <- to_entrez(universe)
  
  kegg_one <- function(genes) {
    if (length(genes) == 0) return(NULL)
    k <- enrichKEGG(gene = genes, universe = uni_e, organism = "hsa",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                    minGSSize = 5, maxGSSize = 5000)
    if (is.null(k) || nrow(as.data.frame(k)) == 0) return(NULL)
    setReadable(k, org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  list(
    mapped = c(up = length(genes_up),     up_entrez = length(up_e),
               down = length(genes_down), down_entrez = length(down_e),
               universe = length(universe), universe_entrez = length(uni_e)),
    up   = kegg_one(up_e),
    down = kegg_one(down_e)
  )
}

# --------------run all three once and save results--------------
if (!file.exists("kegg_results.rds")) {
  kegg_results <- list(
    SS_vs_LS = run_kegg(res_SS_vs_LS),
    IS_vs_LS = run_kegg(res_IS_vs_LS),
    SS_vs_IS = run_kegg(res_SS_vs_IS)
  )
  attr(kegg_results, "date_run") <- Sys.Date()   # KEGG access date for the methods
  saveRDS(kegg_results, "kegg_results.rds")
}
kegg_results <- readRDS("kegg_results.rds")
attr(kegg_results, "date_run")

# --------------entrez mapping counts and number of enriched pathways------------
sapply(kegg_results, function(r) r$mapped)
sapply(kegg_results, function(r) c(
  up   = if (is.null(r$up))   0 else nrow(as.data.frame(r$up)),
  down = if (is.null(r$down)) 0 else nrow(as.data.frame(r$down))))

# -------------dot plots of results (both directions)-------------
prep_kegg <- function(obj, dir_label, n = 10) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df <- df %>% slice_min(p.adjust, n = n, with_ties = FALSE)
  gr <- strsplit(df$GeneRatio, "/")
  df$GeneRatio_num <- vapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  df$Description <- sub(" - Homo sapiens \\(human\\)$", "", df$Description)
  df$direction <- dir_label
  df[, c("Description", "GeneRatio_num", "Count", "p.adjust", "direction")]
}

kegg_dot <- function(r, ttl, up_lab, down_lab) {
  df <- bind_rows(prep_kegg(r$up, up_lab), prep_kegg(r$down, down_lab))
  if (nrow(df) == 0) return(NULL)
  df$direction <- factor(df$direction, levels = c(up_lab, down_lab))
  ggplot(df, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num),
                 size = Count, colour = p.adjust)) +
    geom_point() +
    facet_grid2(direction ~ ., scales = "free_y", space = "free_y",
                strip = strip_themed(
                  background_y = elem_list_rect(fill = c("#D95F02", "#7570B3")),
                  text_y = elem_list_text(colour = "white", face = "bold"))) +
    scale_colour_gradient(low = "red", high = "blue") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
    labs(x = "Gene ratio", y = NULL, title = ttl,
         colour = "Adjusted P", size = "Gene count") +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.y  = element_text(face = "bold", size = 11),
          plot.title   = element_text(face = "bold", size = 14))
}

kegg_plots <- list(
  kegg_dot(kegg_results$SS_vs_LS, "SS vs LS", "Higher in SS", "Higher in LS"),
  kegg_dot(kegg_results$IS_vs_LS, "IS vs LS", "Higher in IS", "Higher in LS"),
  kegg_dot(kegg_results$SS_vs_IS, "SS vs IS", "Higher in SS", "Higher in IS")
)
kegg_plots <- Filter(Negate(is.null), kegg_plots)   # drop contrasts with no pathways

fig_kegg <- plot_grid(plotlist = kegg_plots, ncol = 1, align = "v", axis = "lr",
                      labels = LETTERS[seq_along(kegg_plots)], label_size = 12)
fig_kegg
ggsave("Figure3C_kegg.png", fig_kegg, width = 6.5, height = 9, dpi = 600, bg = "white")

# -------------------dot plot for the one result-----------------
df <- prep_kegg(kegg_results$SS_vs_IS$up, "Higher in SS")
df$panel <- "KEGG"   # label for the side strip, like BP/CC/MF in the GO plots

kegg_SSvIS <- ggplot(df, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num),
                             size = Count, colour = p.adjust)) +
  geom_point() +
  facet_grid2(panel ~ ., scales = "free_y", space = "free_y",
              strip = strip_themed(
                background_y = elem_list_rect(fill = "#1B9E77"),
                text_y = elem_list_text(colour = "white", face = "bold"))) +
  scale_colour_gradient(low = "red", high = "blue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = "Gene ratio", y = NULL, title = "SS vs IS: higher in SS",
       colour = "Adjusted P", size = "Gene count") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.y  = element_text(face = "bold", size = 11, lineheight = 0.9),
        plot.title   = element_text(face = "bold", size = 14))

kegg_SSvIS
ggsave("Supplementary_FigureS2_KEGG_SSvIS.png", kegg_SSvIS,
       width = 9, height = 6, dpi = 600, bg = "white")

# ---------------supplementary table: KEGG-------------------
kegg_table <- bind_rows(lapply(names(kegg_results), function(ct) {
  bind_rows(lapply(c("up", "down"), function(dir) {
    o <- kegg_results[[ct]][[dir]]
    if (is.null(o)) return(NULL)
    as.data.frame(o) %>%
      transmute(Comparison   = gsub("_", " ", ct),
                `Higher in`  = if (dir == "up") grp[[ct]][1] else grp[[ct]][2],
                `KEGG ID`    = ID,
                Pathway      = sub(" - Homo sapiens \\(human\\)$", "", Description),
                `Gene ratio` = GeneRatio,
                `Background ratio` = BgRatio,
                `Gene count` = Count,
                `P value`    = signif(pvalue, 3),
                `Adjusted P` = signif(p.adjust, 3),
                `q value`    = signif(qvalue, 3),
                Genes        = gsub("/", ", ", geneID))
  }))
}))
write_xlsx(list(KEGG_enrichment = kegg_table), "Supplementary_Table_KEGG.xlsx")

as.data.frame(kegg_results$SS_vs_IS$up)[, c("Description", "Count", "p.adjust")]
as.data.frame(kegg_results$SS_vs_LS$up)[, c("Description", "Count", "p.adjust", "geneID")]


########################## PART IV ##############################
# --------------load necessary packages--------------------
library(forestploter)
library(glmnet)
library(grid)
library(pheatmap)
library(STRINGdb)
library(survival)
library(survminer)
library(timeROC)

# -------------univariate cox regression----------------
# Load the 122-sample objects
pheno     <- readRDS("sample_info_122.rds")
vst_mat   <- readRDS("vst_mat_122.rds")
deg_union <- readRDS("deg_union_122.rds")

# Checks
stopifnot(nrow(pheno) == 122,
          ncol(vst_mat) == 122,
          all(colnames(vst_mat) == rownames(pheno)),
          all(deg_union %in% rownames(vst_mat)))
length(deg_union)   # 291 unique sDEGs                            

# ---- survival object ----
surv_time  <- pheno$days_to_death.demographic
surv_event <- rep(1, nrow(pheno))                # all patients deceased: no censoring
surv_obj   <- Surv(time = surv_time, event = surv_event)
sum(is.na(surv_time))                            # 0

# ---- expression of candidate genes, scaled so HR is per 1 SD ----
expr <- scale(t(vst_mat[deg_union, ]))
stopifnot(all(rownames(expr) == rownames(pheno)))

# ---- univariate Cox per gene ----
cox_res <- data.frame()
for (g in colnames(expr)) {
  s <- summary(coxph(surv_obj ~ expr[, g]))
  cox_res <- rbind(cox_res, data.frame(
    gene        = g,
    HR          = s$coefficients[1, "exp(coef)"],
    CI_lower    = s$conf.int[1, "lower .95"],
    CI_upper    = s$conf.int[1, "upper .95"],
    p           = s$coefficients[1, "Pr(>|z|)"],
    concordance = s$concordance[1]))
}

cox_res$padj <- p.adjust(cox_res$p, method = "BH")
cox_res <- cox_res[order(cox_res$p), ]

sum(cox_res$padj < 0.05)                              # 36 significant genes
range(cox_res$concordance[cox_res$padj < 0.05])       # range: 0.5298507 - 0.6168250

# Results:
# So, 36/291 genes are associated with survival and individually, 
# these genes are can modestly predict survival given the concordance index range. 
# Just above 0.5, meaning its slightly better than chance but still not as strong. 

# ---- proportional hazards check on significant genes ----
ph_p <- sapply(cox_res$gene[cox_res$padj < 0.05], function(g) {
  cox.zph(coxph(surv_obj ~ expr[, g]))$table[1, "p"]
})
names(ph_p)[ph_p < 0.05]                              # genes violating PH: PF4

# Check: PF4 violation
cox_res[cox_res$gene == "PF4", ]      # HR, padj, C-index
ph_p["PF4"]                           # Schoenfeld P
km[km$gene == "PF4", ]                # log-rank P and padj

# Results: 
# HR (1.341747), CI (1.101513 to 1.634375), P-value (0.003494947), Padj (0.03789058), C-index (0.5960651)
# Schoenfeld P (0.04563043 )
# Log-rank test P (0.07535186), Padj (0.1130278)

# ------describe the top genes-------
sig_cox <- cox_res[cox_res$padj < 0.05, ]
sig_cox$direction <- ifelse(sig_cox$HR > 1, "Risk (HR > 1)", "Protective (HR < 1)")
# 4 (HR < 1, protective) and 32 (HR > 1, risk)

table(sig_cox$direction)
head(sig_cox[, c("gene", "HR", "CI_lower", "CI_upper", "padj", "concordance")], 10)

# Univariate cox regression basically asks does survival change as the gene goes up?
# Log-rank sweek asks do patients with high versus low expression survive differently?

# --------log-rank sweep (robustness check on Cox-significant genes)-------
sig_genes <- cox_res$gene[cox_res$padj < 0.05]

km <- data.frame()
for (g in sig_genes) {
  grp   <- ifelse(expr[, g] > median(expr[, g]), "High", "Low")
  sdiff <- survdiff(surv_obj ~ grp)
  km <- rbind(km, data.frame(
    gene      = g,
    logrank_p = pchisq(sdiff$chisq, df = 1, lower.tail = FALSE)))
}

km$padj <- p.adjust(km$logrank_p, method = "BH")
km <- km[order(km$logrank_p), ]
sum(km$padj < 0.05)   

# 19 out of the 36 identified by univariate cox regression differ in survival
# based on whether the gene is expressed at a high orlow level. 

# -----visualise the 19 genes in the console----
km_sig <- merge(km[km$padj < 0.05, ],
                cox_res[, c("gene", "HR", "padj", "concordance")],
                by = "gene", suffixes = c("_logrank", "_cox"))
km_sig$direction <- ifelse(km_sig$HR > 1, "Risk (HR > 1)", "Protective (HR < 1)")
km_sig[order(km_sig$padj_logrank), ]

# Results: 3 protective genes and 16 risk genes 

# -----------forest plot: visualisation of univariate cox regression------
forest_df <- cox_res[cox_res$padj < 0.05 &
                       cox_res$gene %in% km$gene[km$padj < 0.05] &
                       cox_res$gene != "PF4", ]          # PF4 violates PH
forest_df <- forest_df[order(-forest_df$HR), ]           # risk genes on top
n <- nrow(forest_df)                                     # expect 19

stars <- as.character(cut(forest_df$padj, c(-Inf, 0.001, 0.01, 0.05, Inf),
                          labels = c("***", "**", "*", "")))
p_txt <- ifelse(forest_df$padj < 0.001, "<0.001",
                formatC(forest_df$padj, format = "f", digits = 3))

forest_df$y      <- n:1
forest_df$hr_txt <- sprintf("%.2f (%.2f\u2013%.2f)", forest_df$HR,
                            forest_df$CI_lower, forest_df$CI_upper)
forest_df$p_txt  <- trimws(paste(p_txt, stars))

# layout (all x positions are on the linear HR scale)
fx     <- c(0.6, 2.0)                     # forest range: widen if any CI is cut off
ticks  <- seq(0.6, 2.0, by = 0.2)         # even 0.2 increments
x_left <- -0.3
x_gene <- -0.25                           # gene column
x_hr   <- 2.2                             # HR (95% CI) column
x_p    <- 3.2                             # adjusted P column
x_max  <- 3.9
cap    <- 0.18                            # CI end-cap half-height
x_mid  <- mean(fx)

shade <- forest_df[forest_df$y %% 2 == 0, ]

fig_cox <- ggplot(forest_df) +
  # alternating row shading
  geom_rect(data = shade, aes(xmin = x_left, xmax = x_max, ymin = y - 0.5, ymax = y + 0.5),
            fill = "grey85") +
  # reference line at HR = 1
  geom_segment(aes(x = 1, xend = 1, y = 0.5, yend = n + 0.5),
               linetype = "dashed", colour = "grey30") +
  # CI lines + end caps
  geom_segment(aes(x = CI_lower, xend = CI_upper, y = y, yend = y), linewidth = 0.5) +
  geom_segment(aes(x = CI_lower, xend = CI_lower, y = y - cap, yend = y + cap), linewidth = 0.5) +
  geom_segment(aes(x = CI_upper, xend = CI_upper, y = y - cap, yend = y + cap), linewidth = 0.5) +
  # point estimate: black square
  geom_point(aes(x = HR, y = y), shape = 15, size = 2.6) +
  # text columns (gene names bold)
  geom_text(aes(x = x_gene, y = y, label = gene),   hjust = 0, size = 3.2, fontface = "bold") +
  geom_text(aes(x = x_hr,   y = y, label = hr_txt), hjust = 0, size = 3.2) +
  geom_text(aes(x = x_p,    y = y, label = p_txt),  hjust = 0, size = 3.2) +
  # headers
  annotate("text", x = x_gene, y = n + 1, label = "Gene",        hjust = 0, size = 3.4, fontface = "bold") +
  annotate("text", x = x_mid,  y = n + 2, label = "Hazard ratio",           size = 4.5, fontface = "bold") +
  annotate("text", x = x_hr,   y = n + 1, label = "HR (95% CI)", hjust = 0, size = 3.4, fontface = "bold") +
  annotate("text", x = x_p,    y = n + 1, label = "Adjusted P",  hjust = 0, size = 3.4, fontface = "bold") +
  # manual x axis under the forest
  annotate("segment", x = fx[1], xend = fx[2], y = 0.5, yend = 0.5) +
  annotate("segment", x = ticks, xend = ticks, y = 0.5, yend = 0.2) +
  annotate("text", x = ticks, y = -0.2, label = sprintf("%.1f", ticks), size = 3) +
  annotate("text", x = x_mid, y = -1.1,
           label = "Hazard ratio per 1 SD increase in expression", size = 3.2) +
  coord_cartesian(xlim = c(x_left, x_max), ylim = c(-1.5, n + 3), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

fig_cox

h <- 0.22 * n + 2
ggsave("Figure5_Cox_forest.png", fig_cox, width = 7.5, height = h, dpi = 600, bg = "white")

# --------------------supplementary table: univariate cox-----------------
uni_table <- data.frame(
  Gene                 = cox_res$gene,
  `Hazard ratio`       = round(cox_res$HR, 3),
  `95% CI lower`       = round(cox_res$CI_lower, 3),
  `95% CI upper`       = round(cox_res$CI_upper, 3),
  `P value`            = signif(cox_res$p, 3),
  `Adjusted P`         = signif(cox_res$padj, 3),
  `C-index`            = round(cox_res$concordance, 3),
  Direction            = ifelse(cox_res$HR > 1, "Risk", "Protective"),
  `PH test P`          = signif(ph_p[cox_res$gene], 3),
  `Log-rank P`         = signif(km$logrank_p[match(cox_res$gene, km$gene)], 3),
  `Log-rank adjusted P` = signif(km$padj[match(cox_res$gene, km$gene)], 3),
  check.names = FALSE)

uni_notes <- data.frame(
  Column = names(uni_table),
  Description = c(
    "Gene symbol (HGNC)",
    "Hazard ratio per 1 standard deviation increase in variance-stabilised expression",
    "Lower bound of the 95% confidence interval",
    "Upper bound of the 95% confidence interval",
    "Unadjusted P value (Wald test, univariate Cox regression)",
    "Benjamini-Hochberg adjusted P value across 291 sDEGs",
    "Concordance index of the single-gene model",
    "Risk: HR > 1 (higher expression, shorter survival); Protective: HR < 1",
    "Schoenfeld residual test of the proportional-hazards assumption (significant genes only); P < 0.05 indicates violation",
    "Log-rank P value comparing median-split high vs low expression (Cox-significant genes only)",
    "Benjamini-Hochberg adjusted log-rank P value across the 36 Cox-significant genes"))

write_xlsx(list(Univariate_Cox = uni_table, Notes = uni_notes),
           "Supplementary_Table_S4_Univariate_Cox.xlsx")


# ------LASSO Cox regression: feature selection---------------
# ---- Helper: genes and coefficients at a given lambda ----
get_panel <- function(cv, s) {
  cf <- as.matrix(coef(cv, s = s))
  p  <- data.frame(gene = rownames(cf)[cf[, 1] != 0], coef = cf[cf[, 1] != 0, 1])
  p[order(-abs(p$coef)), ]
}

# ---- Fit LASSO-Cox (fixed 10-fold split) ----
set.seed(42)
foldid   <- sample(rep(1:10, length.out = nrow(expr)))
cv_lasso <- cv.glmnet(expr, surv_obj, family = "cox", alpha = 1, foldid = foldid)
plot(cv_lasso)

cv_lasso$lambda.min; cv_lasso$lambda.1se

# Results: lambda.min (0.1513827) and lambda.1se (0.1910233)

# ---- Panels at both lambdas ----
panel_min <- get_panel(cv_lasso, "lambda.min")
panel_1se <- get_panel(cv_lasso, "lambda.1se")
c(lambda.min = nrow(panel_min), lambda.1se = nrow(panel_1se))
panel_min     # 17 genes
panel_1se     # 10 genes
all(panel_1se$gene %in% panel_min$gene)   # lambda.1se panel is a subset of lambda.min

# ---- Selection stability: 50 resampled partitions (both lambdas from the same fits) ----
stab_runs <- lapply(1:50, function(i) {
  set.seed(i)
  fid <- sample(rep(1:10, length.out = nrow(expr)))
  fit <- cv.glmnet(expr, surv_obj, family = "cox", alpha = 1, foldid = fid)
  list(min = get_panel(fit, "lambda.min")$gene,
       se  = get_panel(fit, "lambda.1se")$gene)
})                                               # takes a few minutes

stab_summary <- function(sel_list, panel) {
  sizes <- lengths(sel_list)
  n     <- as.integer(table(unlist(sel_list))[panel$gene]); n[is.na(n)] <- 0
  list(empty = sum(sizes == 0), sizes = summary(sizes),
       table = data.frame(gene = panel$gene, coef = round(panel$coef, 3),
                          direction    = ifelse(panel$coef > 0, "Risk", "Protective"),
                          pct_all      = round(100 * n / 50),
                          pct_nonempty = round(100 * n / max(sum(sizes > 0), 1))))
}

stab_min <- stab_summary(lapply(stab_runs, `[[`, "min"), panel_min)
stab_1se <- stab_summary(lapply(stab_runs, `[[`, "se"),  panel_1se)

c(empty_min = stab_min$empty, empty_1se = stab_1se$empty)
stab_min$sizes; stab_1se$sizes
stab_min$table
stab_1se$table

# Results: ADGRG6 (83%), FCN3 (3%), H2AC25 (97%), IGFBP6 (17%), LOXL1 (79%), NRXN3 (100%), 
# PGA3 (3%), PTPRN (97%), USP44 (17%), ZIC3 (100%) using lambda.1se percentages.  
# Lambda.1se is a subset of lambda.min. 
# So, removing those below 50%, ADGRG6, H2AC25, LOXL1,  NRXN3, PTPRN, ZIC3 are my candidate panel

# ----------------core panel: majority rule (>= 50% of non-empty lambda.1se repeats)---------
core <- stab_1se$table$gene[stab_1se$table$pct_nonempty >= 50]
core      # ADGRG6, H2AC25, LOXL1, NRXN3, PTPRN, ZIC3

# -------------save univariate and lasso-cox regression results-----------
saveRDS(list(cv_lasso = cv_lasso, foldid = foldid,
             panel_min = panel_min, panel_1se = panel_1se,
             stab_min = stab_min, stab_1se = stab_1se, stab_runs = stab_runs,
             core = core),
        "lasso_122.rds")
saveRDS(cox_res, "cox_res_122.rds")
saveRDS(core,    "core_panel_122.rds")

# --------------------supplementary table: lasso-cox-----------------
lambda_summary <- data.frame(
  Lambda                  = c("lambda.min", "lambda.1se"),
  Value                   = signif(c(cv_lasso$lambda.min, cv_lasso$lambda.1se), 4),
  `Genes selected`        = c(nrow(panel_min), nrow(panel_1se)),
  `Empty partitions (of 50)` = c(stab_min$empty, stab_1se$empty),
  check.names = FALSE)

lasso_table <- data.frame(
  Gene                          = panel_min$gene,
  `Coefficient (lambda.min)`    = round(panel_min$coef, 3),
  `Coefficient (lambda.1se)`    = round(panel_1se$coef[match(panel_min$gene, panel_1se$gene)], 3),
  Direction                     = ifelse(panel_min$coef > 0, "Risk", "Protective"),
  `Selected, lambda.min (% of 50)` =
    stab_min$table$pct_all[match(panel_min$gene, stab_min$table$gene)],
  `Selected, lambda.1se (% of non-empty)` =
    stab_1se$table$pct_nonempty[match(panel_min$gene, stab_1se$table$gene)],
  `Core panel`                  = ifelse(panel_min$gene %in% core, "Yes", "No"),
  check.names = FALSE)
lasso_table$`Selected, lambda.1se (% of non-empty)`[is.na(lasso_table$`Selected, lambda.1se (% of non-empty)`)] <- 0

lasso_notes <- data.frame(
  Column = c(names(lambda_summary), names(lasso_table)),
  Sheet  = c(rep("Lambda_summary", ncol(lambda_summary)), rep("LASSO_panel", ncol(lasso_table))),
  Description = c(
    "Penalty parameter choice from 10-fold cross-validation",
    "Value of the penalty parameter",
    "Number of genes with non-zero coefficients",
    "Number of 50 resampled 10-fold partitions in which no gene was selected",
    "Gene symbol (HGNC)",
    "LASSO-Cox coefficient (log hazard, per 1 SD) at lambda.min; shrunk, not an effect size",
    "LASSO-Cox coefficient at lambda.1se; blank if not selected",
    "Risk: coefficient > 0; Protective: coefficient < 0",
    "Percentage of 50 resampled partitions in which the gene was selected at lambda.min",
    "Percentage of non-empty partitions (29 of 50) in which the gene was selected at lambda.1se",
    "Core panel: selected in >= 50% of non-empty lambda.1se partitions"))

write_xlsx(list(Lambda_summary = lambda_summary,
                LASSO_panel    = lasso_table,
                Notes          = lasso_notes),
           "Supplementary_Table_S5_LASSO_Cox.xlsx")

# -------------------visualisation of results---------------------
# ---------- thesis table: core LASSO panel (6 genes) ----------
core_tab <- cox_res[match(core, cox_res$gene), ]          # univariate Cox rows for the core genes
core_tab$lasso_coef <- panel_1se$coef[match(core_tab$gene, panel_1se$gene)]   # LASSO coef at lambda.1se
core_tab$stability  <- stab_1se$table$pct_nonempty[match(core_tab$gene, stab_1se$table$gene)]
core_tab <- core_tab[order(-core_tab$lasso_coef), ]       # risk genes first, then protective

fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                            formatC(p, format = "f", digits = 3))

core_table <- data.frame(
  Gene                       = core_tab$gene,
  `LASSO coefficient`        = sprintf("%.3f", core_tab$lasso_coef),
  HR                         = sprintf("%.2f", core_tab$HR),
  `95% CI`                   = sprintf("%.2f\u2013%.2f", core_tab$CI_lower, core_tab$CI_upper),
  `P value`                  = fmt_p(core_tab$p),
  `Adjusted P`               = fmt_p(core_tab$padj),
  `Selected (%)`             = paste0(core_tab$stability, "%"),
  check.names = FALSE
)

core_table
write.csv(core_table, "Table_core_panel_6_genes.csv", row.names = FALSE)

# Copy to clipboard
write.table(core_table, "clipboard", sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------the rest of the results-----------------
lmin      <- log(cv_lasso$lambda.min)
l1se      <- log(cv_lasso$lambda.1se)
xlim_zoom <- c(-3.5, max(log(cv_lasso$lambda)))   # zoom to the region around the chosen lambdas

# - A. cross-validation curve -
# Why: shows why lambda.min and lambda.1se were chosen (both sit on the flat, low-error part of the curve)
cv_df   <- data.frame(loglambda = log(cv_lasso$lambda), cvm = cv_lasso$cvm,
                      lo = cv_lasso$cvlo, hi = cv_lasso$cvup)
cv_zoom <- cv_df[cv_df$loglambda > xlim_zoom[1], ]

pA <- ggplot(cv_df, aes(loglambda, cvm)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), colour = "grey70", width = 0) +
  geom_point(colour = "grey30", size = 1) +
  geom_vline(aes(xintercept = lmin, linetype = "lambda.min")) +
  geom_vline(aes(xintercept = l1se, linetype = "lambda.1se")) +
  scale_linetype_manual(values = c("lambda.min" = "dashed", "lambda.1se" = "dotted"),
                        breaks = c("lambda.min", "lambda.1se")) +
  coord_cartesian(xlim = xlim_zoom, ylim = range(cv_zoom$lo, cv_zoom$hi)) +
  labs(x = expression(log(lambda)), y = "Partial likelihood deviance", linetype = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# - B. coefficient paths -
# Why: shows which genes survive the strongest penalty (core genes are the last to shrink to zero)
beta  <- as.matrix(cv_lasso$glmnet.fit$beta)
paths <- as.data.frame(t(beta[panel_min$gene, ]))
paths$loglambda <- log(cv_lasso$glmnet.fit$lambda)
paths <- pivot_longer(paths, -loglambda, names_to = "gene", values_to = "coef")
paths$core <- factor(ifelse(paths$gene %in% core, "Core panel", "Other selected"),
                     levels = c("Core panel", "Other selected"))
paths <- paths[paths$loglambda > xlim_zoom[1], ]

pB <- ggplot(paths, aes(loglambda, coef, group = gene, colour = core)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = lmin, linetype = "dashed") +
  geom_vline(xintercept = l1se, linetype = "dotted") +
  scale_colour_manual(values = c("Core panel" = "black", "Other selected" = "grey70")) +
  coord_cartesian(xlim = xlim_zoom) +
  labs(x = expression(log(lambda)), y = "Coefficient", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# - C. selection stability at both lambdas -
# Why: shows why the six genes form the core panel (selected in >= 50% of non-empty lambda.1se partitions)
stab_df <- data.frame(
  gene      = panel_min$gene,
  direction = ifelse(panel_min$coef > 0, "Risk", "Protective"),
  `lambda.min (% of 50)` =
    stab_min$table$pct_all[match(panel_min$gene, stab_min$table$gene)],
  `lambda.1se (% of non-empty)` =
    stab_1se$table$pct_nonempty[match(panel_min$gene, stab_1se$table$gene)],
  check.names = FALSE)
stab_df[is.na(stab_df)] <- 0
stab_df$label <- ifelse(stab_df$gene %in% core, paste0(stab_df$gene, " *"), stab_df$gene)

stab_long <- pivot_longer(stab_df, c(`lambda.min (% of 50)`, `lambda.1se (% of non-empty)`),
                          names_to = "lambda", values_to = "pct")
stab_long$lambda <- factor(stab_long$lambda,
                           levels = c("lambda.min (% of 50)", "lambda.1se (% of non-empty)"))

gene_order <- stab_df$label[order(stab_df$`lambda.1se (% of non-empty)`,
                                  stab_df$`lambda.min (% of 50)`)]
stab_long$label <- factor(stab_long$label, levels = gene_order)   # core six at the top

pC <- ggplot(stab_long, aes(pct, label, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 50, linetype = "dashed") +
  facet_wrap(~ lambda) +
  scale_fill_manual(values = c("Risk" = "#F8766D", "Protective" = "#00BFC4")) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Selection frequency (%)", y = NULL, fill = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# - combine and save -
# Why: one supplementary figure showing how the core panel was selected
fig_lasso <- plot_grid(
  plot_grid(pA, pB, ncol = 2, labels = c("A", "B"), label_size = 12),
  pC, ncol = 1, rel_heights = c(1, 1.3), labels = c("", "C"), label_size = 12)

fig_lasso
ggsave("FigureS_LASSO.png", fig_lasso, width = 6.5, height = 8, dpi = 600, bg = "white")

# --------------expression of core genes by survival group---------------
# Why: shows whether the survival-associated genes also differ between SS, IS and LS
core_order <- c("NRXN3", "PTPRN", "LOXL1", "ZIC3", "H2AC25", "ADGRG6")  # risk genes, then protective

box_df <- as.data.frame(expr[, core_order])                 # z-scored vst expression
box_df$group <- factor(pheno$survival_group, levels = c("SS", "IS", "LS"))
box_long <- pivot_longer(box_df, all_of(core_order), names_to = "gene", values_to = "z")

# - Kruskal-Wallis per gene, BH-corrected across the 6 -
# Why: non-parametric test of a difference across the three groups
kw <- box_long |>
  group_by(gene) |>
  summarise(p = kruskal.test(z ~ group)$p.value, .groups = "drop") |>
  mutate(padj = p.adjust(p, method = "BH"),
         direction = ifelse(gene %in% c("NRXN3", "PTPRN", "LOXL1"), "Risk", "Protective"),
         lab = sprintf("%s (%s)\nKW adj. P = %.3g", gene, direction, padj))
kw

box_long$lab <- factor(kw$lab[match(box_long$gene, kw$gene)],
                       levels = kw$lab[match(core_order, kw$gene)])

n <- table(box_df$group)
levels(box_long$group) <- sprintf("%s\n(n = %d)", names(n), n)

# ----box plots----
fig_box <- ggplot(box_long, aes(group, z, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.35) +
  facet_wrap(~ lab, nrow = 2) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), guide = "none") +
  labs(x = "Survival group", y = "Expression (z-score)") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8))

fig_box
ggsave("Figure6_boxplots_core6.png", fig_box, width = 6.5, height = 5, dpi = 600, bg = "white")

# view of changes across the three survival groups
med <- aggregate(z ~ gene + group, data = box_long, FUN = median)
tidyr::pivot_wider(med, names_from = group, values_from = z)

# ------------Kaplan-Meier curves for the core genes (median split)--------------
# Why: shows what each gene's hazard ratio means in survival terms (high vs low expression)
km_order <- c("LOXL1", "NRXN3", "PTPRN",      # risk genes (top row)
              "ADGRG6", "H2AC25", "ZIC3")     # protective genes (bottom row)

# - log-rank P per gene, BH-corrected across the 6 -
# Why: tests whether high- and low-expression groups differ in survival
lr_p <- sapply(km_order, function(g) {
  grp   <- ifelse(expr[, g] > median(expr[, g]), "High", "Low")
  sdiff <- survdiff(surv_obj ~ grp)
  pchisq(sdiff$chisq, df = 1, lower.tail = FALSE)
})
lr_padj <- p.adjust(lr_p, method = "BH")

km_df <- data.frame(time = surv_time, status = surv_event)

# - one KM plot per gene -
km_plots <- lapply(km_order, function(g) {
  df <- km_df
  df$grp <- factor(ifelse(expr[, g] > median(expr[, g]), "High", "Low"),
                   levels = c("Low", "High"))
  fit <- surv_fit(Surv(time, status) ~ grp, data = df)
  hr  <- cox_res$HR[match(g, cox_res$gene)]
  
  ggsurvplot(fit, data = df,
             pval = sprintf("BH p = %.3g", lr_padj[g]),
             pval.size = 3, pval.coord = c(max(df$time) * 0.45, 0.9),
             conf.int = FALSE,
             palette = c("navy", "firebrick"),
             legend.labs = c("Low", "High"), legend.title = "",
             title = sprintf("%s (HR = %.2f)", g, hr),
             xlab = "Time (days)", ylab = "Survival probability",
             break.time.by = 500,
             risk.table = TRUE, risk.table.col = "strata",
             risk.table.y.text = FALSE,          # coloured bars instead of Low/High labels
             risk.table.fontsize = 2.6,
             risk.table.title = "No. at risk",
             surv.plot.height = 0.75,            # <--- ADDED: Controls proportion allocated to the plot vs table
             ggtheme = theme_classic(base_size = 9),
             tables.theme = theme_cleantable(base_size = 7))
})
  
# - arrange: arrange_ggsurvplots fills by column, so reorder to get rows -
# Row 1: LOXL1, NRXN3, PTPRN | Row 2: ADGRG6, H2AC25, ZIC3
km_grid <- arrange_ggsurvplots(km_plots[c(1, 4, 2, 5, 3, 6)],
                               ncol = 3, nrow = 2, print = FALSE)
km_grid
ggsave("Figure_core6_KM.png", km_grid, width = 10, height = 7, dpi = 300, bg = "white")

# - table of results -
data.frame(gene      = km_order,
           HR        = round(cox_res$HR[match(km_order, cox_res$gene)], 2),
           BH_p      = signif(lr_padj, 3),
           row.names = NULL)

# --------------combination of the two plots---------------------
table1_order <- c("NRXN3", "PTPRN", "LOXL1", "ZIC3", "H2AC25", "ADGRG6")
tags <- LETTERS[1:6]

# each gene = KM curve on top, number-at-risk table underneath
km_panels <- lapply(seq_along(table1_order), function(i) {
  p <- km_plots[[which(km_order == table1_order[i])]]
  
  top <- p$plot +
    labs(tag = tags[i]) +
    theme(legend.position = "none",
          plot.title   = element_text(size = 9, face = "bold"),
          axis.title.x = element_blank(),                 # the table carries the x title
          axis.title.y = element_text(size = 8),
          axis.text    = element_text(size = 7),
          plot.tag     = element_text(face = "bold", size = 12),
          plot.tag.position = "topleft")
  
  bottom <- p$table +
    theme(legend.position = "none",
          plot.title   = element_text(size = 7, face = "plain"),
          axis.title.x = element_text(size = 8),
          axis.text.x  = element_text(size = 7),
          axis.text.y  = element_blank(), axis.ticks.y = element_blank())
  
  top / bottom + plot_layout(heights = c(4, 1))
})

panel_a <- wrap_plots(km_panels, ncol = 3, nrow = 2)

panel_b <- fig_box +
  labs(tag = "G") +
  theme(strip.text = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 9),
        axis.text  = element_text(size = 8),
        plot.tag   = element_text(face = "bold", size = 12),
        plot.tag.position = "topleft")

# full A4 page: the boxplots get about half the height
combined_figure <- panel_a / panel_b + plot_layout(heights = c(1, 1.05))
combined_figure

ggsave("Figure5_km_and_boxplots.png", combined_figure,
       width = 8.3, height = 11.2, dpi = 600, bg = "white")

# -------------------risk score from the core panel---------------------
# Multivariable Cox model on the 6 core genes
# Why: gives each gene a weight, adjusted for the other five, to combine into one score
core <- c("LOXL1", "NRXN3", "PTPRN", "ADGRG6", "H2AC25", "ZIC3")
core_df <- as.data.frame(expr[, core])

cox_core <- coxph(surv_obj ~ ., data = core_df)
summary(cox_core)                                  

# - correlations between the six core genes -
# Why: genes that share information can lose significance after adjustment (e.g. PTPRN, ZIC3)
round(cor(core_df), 2)

# - proportional hazards check for the multivariable model -
# Why: the score assumes each gene's effect is constant over time
cox.zph(cox_core)

# - risk score = linear predictor (sum of coefficient x expression) -
# Why: one number per patient; higher = higher predicted risk
risk_score <- predict(cox_core, type = "lp")
summary(risk_score)

# - in-sample concordance (optimistic: model fitted on the same patients) -
concordance(cox_core)$concordance         # 0.7058345 (in-sample)

saveRDS(list(cox_core = cox_core, risk_score = risk_score, core = core),
        "risk_score_core6_122.rds")

# -------------------risk groups and survival---------------------
# - median split of the risk score -
# Why: divides patients into high- and low-risk groups of equal size
risk_grp <- factor(ifelse(risk_score > median(risk_score), "High risk", "Low risk"),
                   levels = c("Low risk", "High risk"))
table(risk_grp)

# - hazard ratio between the two risk groups -
# Why: summarises how much higher the risk is in the high-risk group
cox_grp <- coxph(surv_obj ~ risk_grp)
summary(cox_grp)$conf.int                    # HR (High vs Low) with 95% CI

# - median survival per risk group -
# Why: an intuitive summary of the separation, in days
surv_fit(Surv(surv_time, surv_event) ~ risk_grp,
         data = data.frame(surv_time, surv_event, risk_grp))

# - KM curve by risk group -
# Why: shows the separation in survival between high- and low-risk patients
km_risk_df <- data.frame(time = surv_time, status = surv_event, grp = risk_grp)
fit_risk   <- surv_fit(Surv(time, status) ~ grp, data = km_risk_df)
hr_grp     <- summary(cox_grp)$conf.int

km_risk <- ggsurvplot(fit_risk, data = km_risk_df,
                      pval = TRUE, pval.size = 3.2,
                      pval.coord = c(max(km_risk_df$time) * 0.45, 0.9),
                      conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.25,
                      palette = c("navy", "firebrick"),
                      legend.labs = c("Low risk", "High risk"), legend.title = "",
                      title = sprintf("Six-gene risk score (HR = %.2f, 95%% CI %.2f-%.2f)",
                                      hr_grp[1, "exp(coef)"], hr_grp[1, "lower .95"],
                                      hr_grp[1, "upper .95"]),
                      xlab = "Time (days)", ylab = "Survival probability",
                      ggtheme = theme_classic(base_size = 9),
                      tables.theme = theme_cleantable(base_size = 8))
km_risk

png("Figure_riskscore_KM.png", width = 6.5, height = 5.5, units = "in", res = 600)
print(km_risk)
dev.off()

# -------------------risk plot and heatmap---------------------
# - order patients by risk score -
# Why: the same left-to-right order is used in all three panels
ord <- order(risk_score)
rp_df <- data.frame(rank  = seq_along(ord),
                    score = risk_score[ord],
                    time  = surv_time[ord],
                    grp   = risk_grp[ord])
cut_rank <- sum(rp_df$grp == "Low risk") + 0.5          # position of the median cut-off

grp_cols <- c("Low risk" = "navy", "High risk" = "firebrick")

# - A. risk score -
# Why: shows the spread of the score and where the median splits patients
pA <- ggplot(rp_df, aes(rank, score, colour = grp)) +
  geom_point(size = 0.9) +
  geom_vline(xintercept = cut_rank, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = median(risk_score), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = grp_cols) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Risk score", colour = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top")

# - B. survival time -
# Why: shows that patients with higher risk scores tend to die sooner
pB <- ggplot(rp_df, aes(rank, time, colour = grp)) +
  geom_point(size = 0.9) +
  geom_vline(xintercept = cut_rank, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Survival (days)") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# - C. heatmap of the core genes -
# Why: shows how each gene's expression changes along the risk ranking
heat_order <- c("LOXL1", "NRXN3", "PTPRN", "ADGRG6", "H2AC25", "ZIC3")   # risk on top
heat_df <- as.data.frame(expr[ord, heat_order])
heat_df$rank <- seq_len(nrow(heat_df))
heat_long <- tidyr::pivot_longer(heat_df, -rank, names_to = "gene", values_to = "z")
heat_long$gene <- factor(heat_long$gene, levels = rev(heat_order))
heat_long$z <- pmax(pmin(heat_long$z, 2.5), -2.5)                       # cap extremes for colour

pC <- ggplot(heat_long, aes(rank, gene, fill = z)) +
  geom_tile() +
  geom_vline(xintercept = cut_rank, linetype = "dashed", colour = "grey20") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-2.5, 2.5), name = "Expression\n(z-score)") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Patients ranked by risk score", y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_blank(), legend.position = "right")

# - combine -
fig_risk <- cowplot::plot_grid(pA, pB, pC, ncol = 1, align = "v", axis = "lr",
                               rel_heights = c(1, 1, 1.1),
                               labels = c("A", "B", "C"), label_size = 12)
fig_risk
ggsave("Figure_risk_plot_heatmap.png", fig_risk, width = 6.5, height = 7, dpi = 600, bg = "white")

# ------schoenfeld residuals for the core model (supplementary)------
# Why: visual check that each gene's effect is constant over time
zph_plot <- ggcoxzph(cox.zph(cox_core), font.main = 9, font.x = 8, font.y = 8)
zph_plot

png("FigureS4_Schoenfeld.png", width = 10, height = 8, units = "in", res = 600)
print(zph_plot)
dev.off()

# -------------------out-of-fold performance---------------------
# Aim: To test the performance of the risk model

# Section A. Fixed-panel CV: oof_fixed
oof_fixed <- function(seed, k = 10) {
  set.seed(seed)
  fid <- sample(rep(1:k, length.out = nrow(expr)))
  oof <- numeric(nrow(expr))
  for (f in 1:k) {
    tr  <- fid != f
    fit <- coxph(surv_obj[tr] ~ ., data = core_df[tr, ])
    oof[!tr] <- predict(fit, newdata = core_df[!tr, ], type = "lp")
  }
  oof
}

# Section B. Full pipeline CV: oof_pipeline
oof_pipeline <- function(seed, k = 10) {
  set.seed(seed)
  fid <- sample(rep(1:k, length.out = nrow(expr)))
  oof <- numeric(nrow(expr))
  n_empty <- 0
  for (f in 1:k) {
    tr  <- fid != f
    cvf <- cv.glmnet(expr[tr, ], surv_obj[tr], family = "cox", alpha = 1)
    sel <- get_panel(cvf, "lambda.1se")$gene                
    if (length(sel) == 0) {                                
      oof[!tr] <- 0
      n_empty  <- n_empty + 1
      next
    }
    d_tr <- as.data.frame(expr[tr, sel, drop = FALSE])
    d_te <- as.data.frame(expr[!tr, sel, drop = FALSE])
    fit  <- coxph(surv_obj[tr] ~ ., data = d_tr)
    oof[!tr] <- predict(fit, newdata = d_te, type = "lp")
  }
  attr(oof, "n_empty") <- n_empty
  oof
}

# Section C. C-index: cidx()
cidx <- function(score) concordance(surv_obj ~ score, reverse = TRUE)$concordance

# repeat over 20 random splits
oofA <- sapply(1:20, oof_fixed)                         
oofB_list <- lapply(1:20, oof_pipeline)                      
oofB      <- sapply(oofB_list, as.numeric)                  
empty_B   <- sapply(oofB_list, attr, "n_empty")             

# summary of results
sum(empty_B)         # 133 out of 200 were empty
summary(empty_B)     # Min(5), 1stQ(6), Median(6.5), Mean(6.65), 3rdQ(7), Max(9)

cA <- apply(oofA, 2, cidx)
cB <- apply(oofB, 2, cidx)

# Section D. 
perf <- data.frame(
  Model  = c("In-sample (6 genes)", "Fixed-panel CV", "Full-pipeline CV"),
  Median = round(c(concordance(cox_core)$concordance, median(cA), median(cB)), 3),
  Min    = round(c(NA, min(cA), min(cB)), 3),
  Max    = round(c(NA, max(cA), max(cB)), 3))
perf

# Result: full-pipeline CV, median = 0.529 (0.486 - 0.591)
# Fixed panel CV, median = 0.688 (0.674 - 0.698)
# In-sample (6 genes) = 0.706 (No min and max)

# -------------------comparison with age---------------------
age     <- pheno$age_at_index.demographic
score_B <- rowMeans(oofB)                               

age_cmp <- c(age_alone      = cidx(age),
             score_alone    = cidx(score_B),
             age_plus_score = concordance(coxph(surv_obj ~ age + score_B))$concordance)
round(age_cmp, 3)

# Results: age_alone = 0.592, score_alone = 0.611, age + score = 0.619

cor(age, score_B, use = "complete.obs")    # 0.04053199

# -------------------time-dependent ROC (timeROC)---------------------
# Time-dependent ROC at 1, 2 and 3 years
roc <- timeROC(
  T = surv_time,
  delta = surv_event,
  marker = score_B,
  cause = 1,
  times = c(365, 730, 1095)
)

# Display AUCs as percentages
auc_percent <- round(roc$AUC * 100, 1)
names(auc_percent) <- c("1 year", "2 years", "3 years")
print(auc_percent)

# 1 year: 68.5
# 2 year: 69.4
# 3 year: 62.6

# Save the ROC curves
png("Figure_timeROC.png", width = 5, height = 5,
    units = "in", res = 600)

par(mar = c(4.5, 4.5, 1, 1))

plot(roc, time = 365, col = "#E41A1C",
     lwd = 2.5, title = FALSE)
plot(roc, time = 730, col = "#377EB8",
     lwd = 2.5, add = TRUE)
plot(roc, time = 1095, col = "#4DAF4A",
     lwd = 2.5, add = TRUE)

legend(
  "bottomright",
  bty = "n",
  lty = 1,
  lwd = 2.5,
  col = c("#E41A1C", "#377EB8", "#4DAF4A"),
  legend = sprintf(
    c("1 year (AUC %.1f%%)",
      "2 years (AUC %.1f%%)",
      "3 years (AUC %.1f%%)"),
    auc_percent
  )
)

dev.off()

# Number of patients surviving beyond each time point
sapply(c(365, 730, 1095), function(t) sum(surv_time > t))

# Results: 62 patients survived beyond 1 year, 19 patients survived beyond 2 year, 
# 7 patients survived beyond 3 year. 

# -------------------save---------------------
saveRDS(list(oofA = oofA, oofB = oofB, cA = cA, cB = cB, perf = perf,
             age_cmp = age_cmp, roc = roc),
        "oof_performance_122.rds")

# ------------combined diagram for all risk analyses------------
# ---- refit time-dependent ROC with all three time points ----
roc3 <- timeROC(T = surv_time, delta = surv_event, marker = score_B,
                cause = 1, times = c(365, 730, 1095))
round(roc3$AUC * 100, 1)          # check the three AUCs

# ---- D: risk-group KM + risk table ----
km_risk <- ggsurvplot(fit_risk, data = km_risk_df,
                      pval = TRUE, pval.size = 3.2,
                      pval.coord = c(max(km_risk_df$time) * 0.45, 0.9),
                      conf.int = FALSE,
                      palette = c("navy", "firebrick"),
                      legend.labs = c("Low risk", "High risk"), legend.title = "",
                      title = sprintf("Six-gene risk score (HR = %.2f, 95%% CI %.2f\u2013%.2f)",
                                      hr_grp[1, "exp(coef)"], hr_grp[1, "lower .95"],
                                      hr_grp[1, "upper .95"]),
                      xlab = "Time (days)", ylab = "Survival probability",
                      break.time.by = 500,
                      risk.table = TRUE, risk.table.fontsize = 3,
                      risk.table.y.text = FALSE,
                      risk.table.title = "No. at risk",
                      ggtheme = theme_classic(base_size = 9),
                      tables.theme = theme_cleantable(base_size = 8))

km_risk_panel <- (km_risk$plot + theme(plot.title = element_text(size = 9, face = "bold"),
                                       legend.position = "top")) /
  (km_risk$table + theme(legend.position = "none")) +
  plot_layout(heights = c(3.5, 1))

# ---- E: time-dependent ROC (1, 2 and 3 years) as a ggplot ----
yrs <- c(1, 2, 3)
roc_df <- do.call(rbind, lapply(1:3, function(i) {
  data.frame(FPR   = roc3$FP[, i],
             TPR   = roc3$TP[, i],
             curve = sprintf("%d year%s (AUC %.1f%%)", yrs[i],
                             ifelse(yrs[i] == 1, "", "s"), roc3$AUC[i] * 100))
}))
roc_df$curve <- factor(roc_df$curve, levels = unique(roc_df$curve))

roc_plot <- ggplot(roc_df, aes(FPR, TPR, colour = curve)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity", colour = NULL,
       title = "Time-dependent ROC") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        legend.position = c(0.68, 0.18))

# spacer so the ROC panel lines up with the KM curve, not the risk table
roc_panel <- roc_plot / plot_spacer() + plot_layout(heights = c(3.5, 1))

# ---- bottom row: D and E ----
bottom_panels <- cowplot::plot_grid(
  patchwork::patchworkGrob(km_risk_panel),
  patchwork::patchworkGrob(roc_panel),
  ncol = 2, rel_widths = c(1.1, 1),
  labels = c("D", "E"), label_size = 12)

# ---- full figure: A-C on top, D-E underneath ----
combined_survival_figure <- cowplot::plot_grid(
  fig_risk, bottom_panels,
  ncol = 1, rel_heights = c(1.25, 1))

combined_survival_figure

ggsave("Figure6_RiskScore_KM_TimeROC.png", combined_survival_figure,
       width = 8.3, height = 11.2, dpi = 600, bg = "white")

######################### PART V ##############################
# Aim: compare survival-associated sDEGs with Verhaak (2010) and Crisman (2016) subtype genes
# --------------load necessary packages--------------------
library(dplyr)
library(msigdbr)
library(pheatmap)
library(readxl)
library(survminer)
library(TCGAbiolinks)
library(tidyr)


# --------------thresholds--------------------
PADJ <- 0.05
LFC  <- 1

# --------------load the 122-sample objects--------------------
pheno   <- readRDS("sample_info_122.rds")
vst_mat <- readRDS("vst_mat_122.rds")
res     <- readRDS("res_122.rds")                   # list: SS_vs_LS, IS_vs_LS, SS_vs_IS
stopifnot(nrow(pheno) == 122, all(colnames(vst_mat) == rownames(pheno)))

surv_time <- pheno$days_to_death.demographic
surv_obj  <- Surv(surv_time, rep(1, nrow(pheno)))

# --------------DESeq2 results with gene symbols--------------------
res_list <- lapply(res, function(x) { x <- as.data.frame(x); x$symbol <- rownames(x); x })
universe <- unique(unlist(lapply(res_list, function(x) x$symbol[!is.na(x$padj)])))
length(universe)      # 18390 genes were tested (excluded genes with zero values, )                               

# --------------TCGA transcriptional subtype--------------------
# Why: needed by the marker heatmap and clustering sections below
sub_calls <- TCGAquery_subtype(tumor = "gbm")
sub <- sub_calls$Transcriptome.Subtype[match(pheno$submitter_id, sub_calls$patient)]
table(sub, useNA = "ifany")

# Results: 39 CL, 53 ME, 5 NE, 15 PN (10 NA) = 122 samples

# --------------build the sDEG list--------------------
# Why: one table of every significant gene in each contrast, with its direction
deg_df <- bind_rows(lapply(names(res_list), function(cn) {
  res_list[[cn]] %>%
    filter(!is.na(padj), padj < PADJ, abs(log2FoldChange) > LFC) %>%
    transmute(contrast = cn, symbol, log2FoldChange, padj,
              direction = ifelse(log2FoldChange > 0, "up", "down"))
}))

deg_df %>% dplyr::count(contrast, direction)
nrow(deg_df)                        # total sDEG entries across contrasts: 351 
length(unique(deg_df$symbol))       # unique sDEGs: 291

# --------------Verhaak subtype signatures (MSigDB C2:CGP)--------------------
# Why: published gene sets for the four GBM subtypes (classical, mesenchymal, neural, proneural)
cgp <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")

verhaak <- cgp %>%
  dplyr::filter(grepl("^VERHAAK_GLIOBLASTOMA", gs_name)) %>%
  dplyr::distinct(subtype = sub("VERHAAK_GLIOBLASTOMA_", "", gs_name), symbol = gene_symbol)

verhaak %>% dplyr::count(subtype)

# Results: 161 CL, 216 ME, 128 NE, 177 PN = 682

length(unique(verhaak$symbol))   # 682 unique Verhaak genes

# --------------overlap between sDEGs and Verhaak genes--------------------
# Why: which survival-associated genes are also subtype signature genes
deg_vs_verhaak <- deg_df %>%
  dplyr::left_join(verhaak, by = "symbol", relationship = "many-to-many") %>%
  dplyr::mutate(subtype = tidyr::replace_na(subtype, "-")) %>%
  dplyr::arrange(contrast, direction, subtype == "-", symbol)

overlap_hits <- deg_vs_verhaak %>%
  dplyr::filter(subtype != "-") %>%
  dplyr::select(contrast, direction, symbol, subtype, log2FoldChange, padj)
overlap_hits

overlap_summary <- deg_vs_verhaak %>%
  dplyr::count(contrast, direction, subtype, name = "n_genes") %>%
  tidyr::pivot_wider(names_from = subtype, values_from = n_genes, values_fill = 0)
overlap_summary

# --------------is the overlap more than expected by chance?--------------------
# Why: a few overlapping genes could happen at random; Fisher's test checks against the universe
sdeg_all <- unique(deg_df$symbol)
enrich <- verhaak %>%
  filter(symbol %in% universe) %>%
  group_by(subtype) %>%
  summarise(
    n_signature = n(),
    n_overlap   = sum(symbol %in% sdeg_all),
    expected    = round(n_signature * length(sdeg_all) / length(universe), 2),
    p = fisher.test(matrix(c(n_overlap,
                             length(sdeg_all) - n_overlap,
                             n_signature - n_overlap,
                             length(universe) - length(sdeg_all) - n_signature + n_overlap),
                           nrow = 2), alternative = "greater")$p.value,
    genes = paste(sort(symbol[symbol %in% sdeg_all]), collapse = ", "),
    .groups = "drop") %>%
  mutate(padj = p.adjust(p, method = "BH"))
enrich

# Result: 10/291 sDEGs overlap Verhaak genes (6 ME, 3 NE, 1 PN, 0 CL);
# no subtype enriched beyond chance (all adj. P > 0.5)

# all four Verhaak signatures together (one-sided Fisher's exact test)
v_in <- intersect(unique(verhaak$symbol), universe)        # signature genes that were tested (679)
o    <- sum(v_in %in% sdeg_all)                            # overlapping sDEGs (10)
e    <- length(v_in) * length(sdeg_all) / length(universe) # expected by chance (10.74)
p    <- fisher.test(matrix(c(o, length(sdeg_all) - o,
                             length(v_in) - o,
                             length(universe) - length(sdeg_all) - length(v_in) + o), nrow = 2),
                    alternative = "greater")$p.value       # 0.636

c(tested = length(v_in), overlap = o, expected = round(e, 2), P = round(p, 2))

# --------------------Verhaak et al., 2010: ClaNC centroids----------------------
# Why: gives each gene's average expression in each subtype, so we know if it is a high or low marker
clanc <- read_excel("Verhaak_clanc840_centroids.xls", sheet = 1, skip = 2)
names(clanc)                                         # check the column names
centroid_cols <- c("Classical", "Mesenchymal", "Neural", "Proneural")

# - assign each gene to the subtype where its centroid is most extreme -
# Why: a gene is a "high" marker if strongest positive, "low" if strongest negative
clanc_map <- clanc %>%
  dplyr::rename(symbol = `Gene Symbol`) %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(centroid_cols), as.numeric)) %>%
  dplyr::select(symbol, dplyr::all_of(centroid_cols)) %>%
  tidyr::pivot_longer(-symbol, names_to = "subtype", values_to = "centroid") %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(top    = subtype[which.max(centroid)], top_val = max(centroid),
                   bottom = subtype[which.min(centroid)], bot_val = min(centroid),
                   .groups = "drop") %>%
  dplyr::mutate(marker  = ifelse(abs(top_val) >= abs(bot_val), "high", "low"),
                subtype = tolower(ifelse(marker == "high", top, bottom))) %>%
  dplyr::select(symbol, subtype, marker)

dplyr::count(clanc_map, subtype, marker)

# - check that the 10 overlap genes are present under current names -
# Why: the file uses 2010 gene symbols; renamed genes would silently fail to match
setdiff(unique(overlap_hits$symbol), clanc_map$symbol)    # should be empty

# - sDEGs that are ClaNC markers, with direction in each contrast -
# Why: compare the direction in the survival contrasts with the subtype marker direction
overlap_clanc <- deg_df %>%
  dplyr::inner_join(clanc_map, by = "symbol") %>%
  dplyr::select(contrast, direction, symbol, subtype, marker, log2FoldChange, padj) %>%
  dplyr::arrange(subtype, symbol, contrast)
overlap_clanc

# -----------------------Crisman et al., 2016: 48-gene classifier--------------------
# Why: a second, independent subtype gene panel to compare against the sDEGs
crisman_48 <- c(
  "ACSL3", "ACYP2", "ADCY9", "ARPC1B", "BCAN", "BCL7A", "CASP5", "CDC25A",
  "CHD7", "CLASP2", "COL5A1", "CREB5", "CSPG5", "DCX", "DDX42", "EPHB4",
  "FEZF2", "FGFR3", "FNDC3B", "FZR1", "GSK3B", "GTF2F1", "IL1R1", "LAPTM5",
  "NDRG2", "NR2F6", "PAK3", "PIPOX", "PMP22", "PPM1D", "PRPSAP2", "PTPRC",
  "PURG", "RGS12", "SH3GL3", "SHC1", "SSH3", "TIMP1", "TMBIM1", "TMEM43",
  "TRIB2", "TRRAP", "UROS", "VAX2", "WASF1", "ZDHHC18", "ZEB2", "ZNF446")

sum(crisman_48 %in% universe)                         
deg_df %>% dplyr::filter(symbol %in% crisman_48)  

# Results: IL1R1 and TIMP1 are upregulated in SS vs IS, LFC ~ 1 and padj < 0.05

# Fisher's test
# - is the Crisman overlap more than expected by chance? (one-sided Fisher's exact test) -
# Why: same test as for the Verhaak signatures, so the two references are comparable
c_in <- intersect(crisman_48, universe)                     # Crisman genes that were tested
o_c  <- sum(c_in %in% sdeg_all)                            # overlapping sDEGs
e_c  <- length(c_in) * length(sdeg_all) / length(universe) # expected by chance
p_c  <- fisher.test(matrix(c(o_c, length(sdeg_all) - o_c,
                             length(c_in) - o_c,
                             length(universe) - length(sdeg_all) - length(c_in) + o_c), nrow = 2),
                    alternative = "greater")$p.value

c(tested = length(c_in), overlap = o_c, expected = round(e_c, 2), P = round(p_c, 2))

# Results: overlap (2), expected (0.76), P-value (0.18)

######################### PART VI ##############################
# Aim: To see whether I can propose a panel of genes to classify samples as CL, ME, NE, PN
library(caret)
library(HGNChelper)

# --------------load necessary datasets-------------
dds     <- readRDS("results/dds_glioblastoma_survival_groups.rds")
vst_mat <- readRDS("vst_mat_122.rds")

ncol(dds)                                          # 122 samples
stopifnot(setequal(colnames(vst_mat), colnames(dds)))
vst_mat <- vst_mat[, colnames(dds)]                

# ---------------subtype labels-----------------
# Why: attach each patient's TCGA transcriptome subtype to their sample
subtypes <- TCGAquery_subtype("gbm")

sub_df <- data.frame(sample  = colnames(dds),
                     patient = substr(colnames(dds), 1, 12)) %>%
  left_join(subtypes %>%
              dplyr::select(patient, Transcriptome.Subtype, Original.Subtype),
            by = "patient")

stopifnot(nrow(sub_df) == ncol(dds))               
table(sub_df$Transcriptome.Subtype, useNA = "ifany")
table(sub_df$Transcriptome.Subtype, sub_df$Original.Subtype, useNA = "ifany")
table(sub_df$Transcriptome.Subtype, colData(dds)$survival_group, useNA = "ifany")

colData(dds)$subtype <- sub_df$Transcriptome.Subtype[match(colnames(dds), sub_df$sample)]


# -------------keep labelled samples, three classes---------------
# Why: NE is dropped because it has only 5 samples and is thought to reflect
#      normal-brain contamination; unlabelled samples can't be used to train
keep <- !is.na(colData(dds)$subtype) & colData(dds)$subtype != "NE"
sum(keep)                                          # expect 107

labs4 <- droplevels(factor(colData(dds)$subtype[keep]))
table(labs4)                                       # 39 CL, 53 ME, 15 PN

# -------------candidate pool = all genes minus Verhaak's 840--------------
# Why: to test whether genes outside Verhaak's list can classify subtypes,
#      his genes must be excluded from the pool the new panel is chosen from

clanc    <- read_excel("Verhaak_clanc840_centroids.xls", skip = 2)
clanc840 <- clanc[[1]]
clanc840 <- unique(clanc840[!is.na(clanc840) & clanc840 != "Gene Symbol"])
length(clanc840)                                   # 840 genes

# ---------------update old symbols to current HGNC names---------------
# Why: Verhaak's list is from 2010, so some genes have since been renamed;
#      without this, renamed Verhaak genes would leak into the candidate pool
fixed   <- checkGeneSymbols(clanc840)
new_sym <- ifelse(is.na(fixed$Suggested.Symbol), fixed$x, fixed$Suggested.Symbol)
new_sym[new_sym == "HEXIM1 /// S1PR1"]      <- "S1PR1"
new_sym[new_sym == "HNRNPA3 /// HNRNPA3P1"] <- "HNRNPA3"
clanc840_updated <- unique(new_sym)

sum(clanc840_updated %in% rownames(vst_mat))       # 831 genes retained
setdiff(clanc840_updated, rownames(vst_mat))       

# ---------candidate pool----------
candidates <- setdiff(rownames(vst_mat), clanc840_updated)
length(candidates)                                 # 18580

expr4_cand <- vst_mat[candidates, keep]
dim(expr4_cand)                                    # 18580 x 107

# -------------one-vs-rest DE per subtype--------------
# Why: find non-Verhaak genes that are higher or lower in one subtype than in
#      the other two; this is used later to label which subtype each panel gene marks

stopifnot(all(candidates %in% rownames(dds)))      # dds and vst_mat use the same gene names

dds_sub <- dds[, keep]                             # new object; your survival dds is untouched
dds_sub$subtype <- droplevels(factor(dds_sub$subtype))
design(dds_sub) <- ~ 0 + subtype
dds_sub <- DESeq(dds_sub)

resultsNames(dds_sub)                              # expect subtypeCL, subtypeME, subtypePN
classes <- levels(dds_sub$subtype)
classes                                            # same order: CL, ME, PN

# ----------contrast: one subtype vs the mean of the other two-----------
# Why: weights of 1 for the subtype and -1/2 for each of the others compare
#      that subtype against the average of the rest

ovr <- lapply(seq_along(classes), function(i) {
  w <- rep(-1 / (length(classes) - 1), length(classes))
  w[i] <- 1
  r <- results(dds_sub, contrast = w, alpha = 0.05)
  r <- r[rownames(r) %in% candidates, ]            # keep non-Verhaak genes only
  r$padj <- p.adjust(r$pvalue, method = "BH")      # re-adjust within the candidate pool
  r
})
names(ovr) <- classes

# ---------number of candidate markers per subtype---------
sapply(ovr, function(r)
  c(up   = sum(r$padj < 0.05 & r$log2FoldChange >  1, na.rm = TRUE),
    down = sum(r$padj < 0.05 & r$log2FoldChange < -1, na.rm = TRUE)))

# Results: 
#      CL   ME   PN
# up   172 1469 1060
# down 928 439  951

# ---------sanity check: top 15 up-regulated genes per subtype--------
# Why: these should make biological sense (e.g. immune/myeloid genes for ME,
#      oligodendrocyte/neural genes for PN)

lapply(ovr, function(r) {
  r <- r[!is.na(r$padj) & r$padj < 0.05 & r$baseMean > 50 & r$log2FoldChange > 0, ]
  head(rownames(r[order(-r$stat), ]), 15)
})

# ---------- subtype versus survival group ----------
tab_sub <- table(Subtype = colData(dds_sub)$subtype, Group = colData(dds_sub)$survival_group)
tab_sub
set.seed(2026)
chisq.test(tab_sub, simulate.p.value = TRUE, B = 1e5)   # report chi-square (7.0477), df (NA) and this P (0.1323)

# ---------- panel-only versus Verhaak, paired within each repeat ----------
b <- do.call(rbind, lapply(runs, `[[`, "bal"))
d <- b[, "Panel"] - b[, "Verhaak"]
round(c(median = median(d), min = min(d), max = max(d)), 3)

# median (-0.031), min (-0.077), max (0.015)

# -------save results-----
# Why: DESeq on 107 samples is slow; reload this instead of rerunning
saveRDS(dds_sub, "results/dds_subtype_107.rds")
saveRDS(ovr, "subtype_ovr_107.rds")


# ----------------------elastic net: feature selection-------------------
# Why: 18 580 genes is too many for 107 samples; dropping low-expression genes
#      and keeping the 5,000 most variable uses no subtype labels, so no leakage

bm        <- rowMeans(counts(dds_sub, normalized = TRUE))
cand_expr <- candidates[bm[candidates] > 50]
length(cand_expr)         # 13 459 genes passed the expression filter

v       <- apply(expr4_cand[cand_expr, ], 1, var)
top5000 <- names(sort(v, decreasing = TRUE))[1:5000]

X <- t(expr4_cand[top5000, ])          # 107 samples x 5000 genes
y <- labs4
stopifnot(identical(rownames(X), colnames(dds_sub)))
dim(X)

# ---------helper functions--------
# Why: genes with a non-zero coefficient in any class (grouped = in or out
#      for all three classes together)

get_genes <- function(fit, s) {
  cf <- coef(fit, s = s)
  m  <- do.call(cbind, lapply(cf, as.matrix))
  setdiff(rownames(m)[rowSums(m != 0) > 0], "(Intercept)")
}

# Why: balanced accuracy = mean recall per class, so the 15 PN samples count
#      as much as the 53 ME samples
bal_acc <- function(truth, pred) {
  cm <- table(truth, factor(pred, levels = levels(truth)))
  mean(diag(cm) / rowSums(cm))
}

# Why: elastic net (alpha = 0.5) instead of LASSO because amplicon genes are
#      highly correlated; LASSO would keep one at random, elastic net keeps the group
fit_enet <- function(Xm, ym) {
  foldid <- createFolds(ym, k = 5, list = FALSE)   # stratified by subtype
  cv.glmnet(Xm, ym, family = "multinomial", type.multinomial = "grouped",
            alpha = 0.5, type.measure = "class", foldid = foldid)
}


# ---------------nested CV = honest accuracy estimate--------------
# Why: genes are selected inside each training fold and tested on the
#      held-out fold, so accuracy is not inflated by selection

set.seed(2026)
outer   <- createFolds(y, k = 5)
pred    <- character(length(y))
n_genes <- integer(length(outer))

for (f in seq_along(outer)) {
  te  <- outer[[f]]
  tr  <- setdiff(seq_along(y), te)
  cvf <- fit_enet(X[tr, ], y[tr])
  pred[te]   <- predict(cvf, X[te, ], s = "lambda.1se", type = "class")[, 1]
  n_genes[f] <- length(get_genes(cvf, "lambda.1se"))
}

table(truth = y, pred = factor(pred, levels = levels(y)))
# Results:
# prediction:
# truth CL ME PN
#    CL 33  6  0
#    ME  4 49  0
#    PN  3  1  11

bal_acc(y, pred)      # 0.8346718
n_genes               # 69  78 110  48  55                                       

# -----------final panel on all 107 samples----------
set.seed(2026)
cv_full <- fit_enet(X, y)
plot(cv_full)
panel <- get_genes(cv_full, "lambda.1se")
length(panel) # 107 genes 
panel

# ----------------stability across 50 resamples------------
# Why: panel size varied 48-110 across folds, so the single 107-gene panel is
#      not reliable; genes picked in most resamples are the robust core

set.seed(2026)
t0 <- Sys.time()
sel <- lapply(1:50, function(i) {
  idx <- createDataPartition(y, p = 0.8, list = FALSE)[, 1]   # stratified 80%
  get_genes(fit_enet(X[idx, ], y[idx]), "lambda.1se")
})
Sys.time() - t0                                    # run time: 3.56 mins

# Save straight away so you never have to rerun this
saveRDS(sel, "subtype_stability_sel.rds")

# -------------how often each gene was selected-----------
freq <- sort(table(unlist(sel)), decreasing = TRUE)
head(freq, 30)

summary(sapply(sel, length))                       # panel size across resamples

sum(freq >= 40)                                    # >= 80% of runs = 27
sum(freq >= 45)                                    # >= 90% of runs = 19
sum(freq == 50)                                    # every run =7 

# ----------how the 107-gene full-data panel holds up--------------
stab <- data.frame(gene = panel, freq = as.integer(freq[panel]))
stab$freq[is.na(stab$freq)] <- 0
stab <- stab[order(-stab$freq), ]
head(stab, 30)
table(cut(stab$freq, breaks = c(-1, 9, 24, 39, 50),
          labels = c("<10", "10-24", "25-39", "40-50")))

# < 10 (6), 10-24 (43), 25-39 (31), 40-50 (27)

saveRDS(list(sel = sel, freq = freq, stab = stab), "subtype_stability.rds")

# ---------- stable core panel ----------
core_sub <- names(freq)[freq >= 40]            # selected in >= 80% of the 50 resamples
length(core_sub)                               # 27

# ---------- Verhaak matrix (same samples as X) ----------
verhaak_genes <- intersect(clanc840_updated, rownames(vst_mat))     # 831
XV   <- t(vst_mat[verhaak_genes, keep])                             # 107 x 831
Xall <- cbind(XV, X)
stopifnot(identical(rownames(XV), rownames(X)), !any(duplicated(colnames(Xall))))

# ---------- one repeat of the swap experiment ----------
swap_once <- function(rep_id, r_cuts = c(0.5, 0.7), fixed_panel = NULL) {
  set.seed(1000 + rep_id)
  outer  <- createFolds(y, k = 5)
  models <- c("Verhaak", "Panel", paste0("Swap_r", r_cuts), paste0("Control_r", r_cuts))
  pred   <- setNames(lapply(models, function(m) character(length(y))), models)
  removed <- matrix(NA, length(outer), length(r_cuts), dimnames = list(NULL, paste0("r", r_cuts)))
  psize   <- integer(length(outer))
  
  for (f in seq_along(outer)) {
    te <- outer[[f]]; tr <- setdiff(seq_along(y), te)
    
    # 1. panel from non-Verhaak genes, chosen on the training samples only
    panel_f <- if (is.null(fixed_panel)) get_genes(fit_enet(X[tr, ], y[tr]), "lambda.1se") else fixed_panel
    stopifnot(length(panel_f) >= 2)
    psize[f] <- length(panel_f)
    
    # 2. highest |r| between each Verhaak gene and any panel gene (training samples only)
    maxr <- apply(abs(cor(XV[tr, ], X[tr, panel_f, drop = FALSE])), 1, max)
    
    # 3. feature sets
    sets <- list(Verhaak = colnames(XV), Panel = panel_f)
    for (j in seq_along(r_cuts)) {
      drop_g <- names(maxr)[maxr >= r_cuts[j]]
      removed[f, j] <- length(drop_g)
      sets[[paste0("Swap_r", r_cuts[j])]] <- c(setdiff(colnames(XV), drop_g), panel_f)
      rnd_v <- sample(colnames(XV), length(drop_g))                        # control
      rnd_n <- sample(setdiff(colnames(X), panel_f), length(panel_f))
      sets[[paste0("Control_r", r_cuts[j])]] <- c(setdiff(colnames(XV), rnd_v), rnd_n)
    }
    
    # 4. fit on training, predict the held-out fold
    for (m in models) {
      fit <- fit_enet(Xall[tr, sets[[m]], drop = FALSE], y[tr])
      pred[[m]][te] <- as.character(predict(fit, Xall[te, sets[[m]], drop = FALSE],
                                            s = "lambda.1se", type = "class")[, 1])
    }
  }
  list(bal = sapply(pred, function(p) bal_acc(y, p)),
       cm  = lapply(pred, function(p) table(truth = y, pred = factor(p, levels = levels(y)))),
       removed = removed, psize = psize)
}

# ---------- run: 10 repeats (roughly 10-20 minutes) ----------
runs       <- lapply(1:10, swap_once)                                   # panel chosen inside each fold
runs_fixed <- lapply(1:10, swap_once, fixed_panel = core_sub)           # 27-gene panel fixed (optimistic)
saveRDS(list(runs = runs, runs_fixed = runs_fixed, core_sub = core_sub), "subtype_swap_results.rds")

# ---------- what each core gene marks ----------
ovr_df <- lapply(ovr, as.data.frame)                    # CL, ME, PN one-vs-rest results
core_tab <- do.call(rbind, lapply(core_sub, function(g) {
  st <- sapply(ovr_df, function(r) r[g, "stat"])        # Wald statistic per contrast
  k  <- names(st)[which.max(abs(st))]                   # subtype where the gene differs most
  data.frame(gene = g, freq_pct = 2 * as.integer(freq[g]), subtype = k,
             direction = ifelse(st[k] > 0, "higher", "lower"),
             log2FC = round(ovr_df[[k]][g, "log2FoldChange"], 2),
             padj = signif(ovr_df[[k]][g, "padj"], 3), row.names = NULL)
}))
core_tab <- core_tab[order(core_tab$subtype, core_tab$direction, -core_tab$freq_pct), ]
core_tab

# ---------- summarise the swap experiment ----------
model_levels <- c("Verhaak", "Swap_r0.5", "Swap_r0.7", "Control_r0.5", "Control_r0.7", "Panel")
bal_df <- do.call(rbind, lapply(seq_along(runs), function(i)
  data.frame(rep = i, model = names(runs[[i]]$bal), bal = unname(runs[[i]]$bal))))
bal_df$model <- factor(bal_df$model, levels = model_levels)

bal_sum <- bal_df %>% group_by(model) %>%
  summarise(median = round(median(bal), 3), min = round(min(bal), 3), max = round(max(bal), 3),
            .groups = "drop")
bal_sum

# swapped model minus Verhaak model, paired within each repeat
bal_diff <- bal_df %>% pivot_wider(names_from = model, values_from = bal) %>%
  mutate(across(starts_with("Swap_"),    ~ . - Verhaak, .names = "diff_{.col}"),
         across(starts_with("Control_"), ~ . - Verhaak, .names = "diff_{.col}")) %>%
  dplyr::select(starts_with("diff_")) %>%
  summarise(across(everything(), list(median = median, min = min, max = max)))
round(t(bal_diff), 3)

rem <- do.call(rbind, lapply(runs, `[[`, "removed"))     # Verhaak genes removed per fold
apply(rem, 2, function(x) c(median = median(x), min = min(x), max = max(x)))
table(unlist(lapply(runs, `[[`, "psize")))               # panel sizes chosen inside folds

bal_fixed <- do.call(rbind, lapply(runs_fixed, `[[`, "bal"))
round(apply(bal_fixed, 2, function(x) c(median = median(x), min = min(x), max = max(x))), 3)

# ---------- Figure 7: confusion matrices + stable core panel ----------
# 1. Generate confusion matrix function
cm_plot <- function(model, title) {
  cm <- Reduce(`+`, lapply(runs, function(r) r$cm[[model]]))
  df <- as.data.frame(prop.table(cm, 1)); names(df)[3] <- "prop"
  ggplot(df, aes(pred, truth, fill = prop)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * prop)), size = 3.2) +
    scale_fill_gradient(low = "white", high = "#377EB8", limits = c(0, 1), guide = "none") +
    scale_y_discrete(limits = rev(levels(y))) +
    labs(x = "Predicted subtype", y = "TCGA subtype", title = title) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"), panel.grid = element_blank())
}

# 2. Build individual panels A, B, D, E
pA <- cm_plot("Verhaak", "Verhaak genes")
pB <- cm_plot("Panel",    "Non-Verhaak panel")

mod_cols <- c(Verhaak = "grey40", Swap_r0.5 = "#377EB8", Swap_r0.7 = "#377EB8",
              Control_r0.5 = "grey80", Control_r0.7 = "grey80", Panel = "#FF7F00")
mod_labs <- c(Verhaak = "Verhaak\ngenes", Swap_r0.5 = "Swap\n|r| \u2265 0.5", Swap_r0.7 = "Swap\n|r| \u2265 0.7",
              Control_r0.5 = "Random\nswap (0.5)", Control_r0.7 = "Random\nswap (0.7)", Panel = "Panel\nonly")

pD <- ggplot(bal_df, aes(model, bal, fill = model)) +
  geom_hline(yintercept = 1/3, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
  scale_fill_manual(values = mod_cols, guide = "none") +
  scale_x_discrete(labels = mod_labs) +
  coord_cartesian(ylim = c(0.75, 1)) +
  labs(x = NULL, y = "Balanced accuracy (nested CV)",
       title = "Performance after swapping correlated genes") +
  theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))

rem_df <- data.frame(cut = factor(rep(colnames(rem), each = nrow(rem)), levels = colnames(rem)),
                     n = as.vector(rem))
pE <- ggplot(rem_df, aes(cut, n)) +
  geom_boxplot(outlier.shape = NA, fill = "#377EB8", alpha = 0.6) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
  scale_x_discrete(labels = c(r0.5 = "|r| \u2265 0.5", r0.7 = "|r| \u2265 0.7")) +
  labs(x = NULL, y = "Verhaak genes removed", title = "Genes removed") +
  theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))

# 3. Combine into a 2x2 layout using patchwork (A | B over D | E)
row1 <- (pA | pB)
row2 <- (pD | pE) + plot_layout(widths = c(2.4, 1))

combined_fig <- (row1 / row2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 12))

# 4. Preview and save out to a high-res PNG
combined_fig
ggsave("Figure_ABDE_combined.png", combined_fig, width = 8, height = 8, dpi = 600, bg = "white")

# ---------- Figure 8: heatmap of the core panel ----------
surv_grp <- colData(dds_sub)$survival_group
genes_hm <- as.character(core_tab$gene)
hm       <- pmax(pmin(t(scale(X[, genes_hm])), 2.5), -2.5)           # genes x samples, z-scores
ord      <- order(y)
ann_col  <- data.frame(Subtype = y, Survival = surv_grp, row.names = rownames(X))
ann_row  <- data.frame(`Marker of` = core_tab$subtype, check.names = FALSE, row.names = genes_hm)
ann_cols <- list(Subtype = sub_cols, `Marker of` = sub_cols,
                 Survival = c(SS = "#E41A1C", IS = "#377EB8", LS = "#4DAF4A"))

pheatmap(hm[, ord], cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = ann_col[ord, , drop = FALSE], annotation_row = ann_row,
         annotation_colors = ann_cols, show_colnames = FALSE, fontsize = 8, fontsize_row = 7,
         gaps_col = cumsum(table(y)[-nlevels(y)]),
         color = colorRampPalette(c("navy", "white", "firebrick"))(100),
         filename = "Figure8_subtype_heatmap.png", width = 7.5, height = 6.5)

# ---------- supplementary tables ----------
write_xlsx(list(Core_panel_27 = core_tab, Swap_summary = as.data.frame(bal_sum),
                Verhaak_overlap = as.data.frame(overlap_hits), Enrichment = as.data.frame(enrich)),
           "Supplementary_Table_S6_subtype_panel.xlsx")

# -----------supplementary table: overlap with established genes----------
fisher_up <- function(set) {
  set <- intersect(set, universe); o <- sum(set %in% sdeg_all)
  data.frame(`Genes tested` = length(set), `Overlap with sDEGs` = o,
             `Expected by chance` = round(length(set) * length(sdeg_all) / length(universe), 2),
             P = signif(fisher.test(matrix(c(o, length(sdeg_all) - o, length(set) - o,
                                             length(universe) - length(sdeg_all) - length(set) + o), nrow = 2),
                                    alternative = "greater")$p.value, 3), check.names = FALSE)
}
tests <- rbind(
  cbind(Set = enrich$subtype, do.call(rbind, lapply(enrich$subtype, function(s)
    fisher_up(verhaak$symbol[verhaak$subtype == s])))),
  cbind(Set = "All Verhaak signatures", fisher_up(unique(verhaak$symbol))),
  cbind(Set = "Crisman 48 genes",       fisher_up(crisman_48)))
tests$`Adjusted P (BH, four Verhaak subtypes)` <- c(signif(enrich$padj, 3), NA, NA)

s6_overlap <- overlap_clanc %>%
  dplyr::transmute(Gene = symbol, Contrast = contrast, Direction = direction,
                   `log2 fold change` = round(log2FoldChange, 3), `Adjusted P` = signif(padj, 3),
                   `ClaNC subtype` = subtype, `ClaNC marker` = marker)

s6_crisman <- deg_df %>% dplyr::filter(symbol %in% crisman_48) %>%
  dplyr::transmute(Gene = symbol, Contrast = contrast, Direction = direction,
                   `log2 fold change` = round(log2FoldChange, 3), `Adjusted P` = signif(padj, 3))

tab_sub <- table(Subtype = colData(dds_sub)$subtype, Group = colData(dds_sub)$survival_group)
s6_sub  <- cbind(Subtype = rownames(tab_sub), as.data.frame.matrix(tab_sub)[, c("SS", "IS", "LS")])
rownames(s6_sub) <- NULL
chisq.test(tab_sub, simulate.p.value = TRUE, B = 1e5)         # report this P

s6_notes <- data.frame(
  Sheet = c("Verhaak_ClaNC_overlap", "Enrichment_tests", "Crisman_overlap", "Subtype_by_survival"),
  Description = c(
    "sDEGs that are Verhaak ClaNC genes, with contrast, direction (positive = higher in the first group), and whether the gene is a high or low marker of its subtype",
    "One-sided Fisher's exact test of the overlap between the 291 sDEGs and each signature, against the 18,390 tested genes; BH adjustment across the four Verhaak subtypes only",
    "sDEGs that are in the Crisman et al. (2016) 48-gene classifier",
    "Number of samples of each TCGA subtype in each survival group (107 labelled samples, NE excluded)"))

write_xlsx(list(Verhaak_ClaNC_overlap = s6_overlap, Enrichment_tests = tests,
                Crisman_overlap = s6_crisman, Subtype_by_survival = s6_sub, Notes = s6_notes),
           "Supplementary_Table_S6_subtype_overlap.xlsx")

mark_gene <- function(g) {
  st <- sapply(ovr_df, function(r) r[g, "stat"]); k <- names(st)[which.max(abs(st))]
  data.frame(Gene = g, `Marks subtype` = k, `Higher or lower` = ifelse(st[k] > 0, "higher", "lower"),
             `log2 fold change` = round(ovr_df[[k]][g, "log2FoldChange"], 2),
             `Adjusted P` = signif(ovr_df[[k]][g, "padj"], 3), check.names = FALSE, row.names = NULL)
}
core_genes <- as.character(core_tab$gene)
s7_core <- cbind(do.call(rbind, lapply(core_genes, mark_gene)),
                 `Selected in resamples (%)` = 2 * as.integer(freq[core_genes]))
s7_full <- cbind(do.call(rbind, lapply(stab$gene, mark_gene)),
                 `Selected in resamples (%)` = 2 * stab$freq,
                 `In stable core` = stab$gene %in% core_genes)

model_names <- c(Verhaak = "Verhaak genes", Panel = "Panel only",
                 Swap_r0.5 = "Swap (|r| >= 0.5)", Swap_r0.7 = "Swap (|r| >= 0.7)",
                 Removed_r0.5 = "Removal only (|r| >= 0.5)", Removed_r0.7 = "Removal only (|r| >= 0.7)",
                 Control_r0.5 = "Random swap (matched to 0.5)", Control_r0.7 = "Random swap (matched to 0.7)")
perf_table <- function(run_list, label) {
  b <- do.call(rbind, lapply(run_list, `[[`, "bal"))      # repeats x models
  d <- b - b[, "Verhaak"]                                 # paired difference from the Verhaak model
  data.frame(Analysis = label,
             Model = unname(ifelse(colnames(b) %in% names(model_names), model_names[colnames(b)], colnames(b))),
             Median = round(apply(b, 2, median), 3), Min = round(apply(b, 2, min), 3), Max = round(apply(b, 2, max), 3),
             `Diff from Verhaak, median` = round(apply(d, 2, median), 3),
             `Diff from Verhaak, min`    = round(apply(d, 2, min), 3),
             `Diff from Verhaak, max`    = round(apply(d, 2, max), 3),
             check.names = FALSE, row.names = NULL)
}
s7_perf <- rbind(perf_table(runs, "Panel selected inside each fold"),
                 perf_table(runs_fixed, "Core panel fixed (selected on all samples; optimistic)"))

rem <- do.call(rbind, lapply(runs, `[[`, "removed"))
s7_removed <- data.frame(Cutoff = sub("r", "|r| >= ", colnames(rem)),
                         Median = apply(rem, 2, median), Min = apply(rem, 2, min), Max = apply(rem, 2, max),
                         `Percent of Verhaak genes (median)` = round(100 * apply(rem, 2, median) / ncol(XV), 1),
                         check.names = FALSE, row.names = NULL)

s7_first <- data.frame(`TCGA subtype` = c("CL", "ME", "PN"),
                       `Predicted CL` = c(33, 4, 3), `Predicted ME` = c(6, 49, 1), `Predicted PN` = c(0, 0, 11),
                       check.names = FALSE)     # first nested CV run, seed 2026: balanced accuracy 0.835

s7_notes <- data.frame(
  Sheet = c("Core_panel_27", "Full_data_panel_107", "Performance", "Genes_removed", "First_nested_CV"),
  Description = c(
    "Genes selected in at least 80% of 50 stratified 80% resamples; subtype where the one-versus-rest Wald statistic was largest, and whether the gene is higher or lower there",
    "Genes with non-zero coefficients in the elastic net fitted on all 107 samples, with their selection frequency across the 50 resamples",
    "Balanced accuracy (nested cross-validation, 10 repeats of 5 folds): median, minimum and maximum, and the paired difference from the Verhaak model within each repeat",
    "Verhaak genes removed per fold because of absolute correlation with a panel gene (10 repeats x 5 folds = 50 values)",
    "Confusion counts of the first nested cross-validation run (seed 2026)"))

write_xlsx(list(Core_panel_27 = s7_core, Full_data_panel_107 = s7_full, Performance = s7_perf,
                Genes_removed = s7_removed, First_nested_CV = s7_first, Notes = s7_notes),
           "Supplementary_Table_S7_subtype_classifier.xlsx")

# ----------supplementary figure-------------
# Define subtype colors
sub_cols <- c(CL = "#E41A1C", ME = "#377EB8", PN = "#4DAF4A")

# Order genes by frequency percentage for clean display
core_tab$gene <- factor(core_tab$gene, levels = rev(core_tab$gene[order(core_tab$freq_pct, core_tab$gene)]))

# Create the standalone plot
fig_supp <- ggplot(core_tab, aes(freq_pct, gene, fill = subtype)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = direction), hjust = -0.1, size = 3) +
  scale_fill_manual(values = sub_cols, name = "Marker of") +
  scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 20)) +
  labs(
    x = "Selected in resamples (%)", 
    y = NULL, 
    title = "Stable core panel"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"), 
    axis.text.y = element_text(size = 9),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.direction = "horizontal"
  )

# Display plot in R
fig_supp

# Save as a high-resolution supplementary figure PNG
ggsave("FigureS_stable_core_panel.png", fig_supp, width = 6.5, height = 7.5, dpi = 600, bg = "white")

