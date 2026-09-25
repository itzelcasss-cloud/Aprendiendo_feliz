#Resultados de Tesis de Montserrat Ramirez 

library(limma)
library(GEOquery)
library(tidyverse)
BiocManager::install("edgeR")
library(edgeR)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library(apeglm)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")
library(enrichplot)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#GSE188915 - Diferencias de tejido endometrial con y sin tratamiento hormonal
gse_188915 <- getGEO("GSE188915", GSEMatrix = TRUE)
pheno_gse_188915 <- pData(gse_188915[[1]])
conteo_gse_188915 <- exprs((gse_188915[[1]]))

#esta incompleto el conteo
getGEOSuppFiles("GSE188915")
list.files("GSE188915")
conteo_gse_188915 <- read.csv(
  "GSE188915/GSE188915_Organoid_TMM_Normalised_MTremoved_CPM_Counts_ENSEMBL.csv.gz",
  row.names = 1,
  check.names = FALSE
)
dim(conteo_gse_188915)
summary(as.vector(as.matrix(conteo_gse_188915)))
dim(conteo_gse_188915)
colnames(conteo_gse_188915)
pheno_gse_188915$title
data.frame(
  sample = pheno_gse_188915$title,
  tratamiento = pheno_gse_188915$`hormonal medication:ch1`
)
colnames(conteo_gse_188915)[
  colnames(conteo_gse_188915) == "Sanple_4"
] <- "Sample_4"
conteo_gse_188915 <- conteo_gse_188915[, c(
  "Sample_1",
  "Sample_2",
  "Sample_3",
  "Sample_4",
  "Sample_5",
  "Sample_6",
  "Sample_7"
)]
data.frame(
  sample = pheno_gse_188915$title,
  treatment = pheno_gse_188915$`hormonal medication:ch1`
)

raw_gse188915 <- read.csv(
  "GSE188915/GSE188915_Organoid_Raw_Gene_Counts_ENSEMBL.csv.gz",
  row.names = 1,
  check.names = FALSE
)
colnames(raw_gse188915)[
  colnames(raw_gse188915) == "Sanple_4"
] <- "Sample_4"
raw_gse188915 <- raw_gse188915[, c(
  "Sample_1",
  "Sample_2",
  "Sample_3",
  "Sample_4",
  "Sample_5",
  "Sample_6",
  "Sample_7"
)]
colnames(raw_gse188915)
grupo_gse188915 <- factor(
  c(
    "None",
    "None",
    "OCP",
    "Mirena",
    "OCP",
    "None",
    "None"
  ),
  levels = c("None", "Mirena", "OCP")
)
table(grupo_gse188915)
design_gse188915 <- model.matrix(
  ~ 0 + grupo_gse188915
)
colnames(design_gse188915) <- c(
  "None",
  "Mirena",
  "OCP"
)
ncol(raw_gse188915) == nrow(design_gse188915)

#listo
view(pheno_gse_188915)


md_pheno_gse_188915 <- pheno_gse_188915 |>
  select( c(geo_accession,
            'hormonal medication:ch1',
            source_name_ch1,
            title))

md_pheno_gse_188915$title_2 <- colnames(conteo_gse_188915)
view(md_pheno_gse_188915)

md_pheno_gse_188915 <- 
  md_pheno_gse_188915 |> 
  rownames_to_column("id") |> 
  column_to_rownames("title_2")

md_pheno_gse_188915 <- 
  md_pheno_gse_188915 |> 
  dplyr::rename("hormone" = `hormonal medication:ch1`)

md_pheno_gse_188915 <- 
  md_pheno_gse_188915 |> 
  mutate(hormone = factor(hormone))
md_pheno_gse_188915$hormone <- relevel(md_pheno_gse_188915$hormone, "None")

#para quitar los "NA" y hacer más pequeño el conteo
rowSums(raw_gse188915 > 30) > (7*0.40) #tienen más de 30 en el 40% de los id

md_pheno_gse_188915 <- md_pheno_gse_188915 |> 
  mutate(treatment = factor(ifelse(hormone== "None", yes = "None", no = "hormonal")))

count_less <- raw_gse188915[rowSums(raw_gse188915 > 30) > (7*0.25),]
dim(count_less)
md_pheno_gse_188915$treatment <- relevel(md_pheno_gse_188915$treatment, "None") 

#DESEQ2
dds <- DESeqDataSetFromMatrix(countData = raw_gse188915,
                              colData = md_pheno_gse_188915,
                              design= ~ treatment)
dds <- DESeq(dds)
resultsNames(dds) 
res <- results(dds, name="treatment_hormonal_vs_None")
reslfc <- lfcShrink(dds, coef="treatment_hormonal_vs_None", type="apeglm")


#graficar 
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01,
  FCcutoff = 1,
  title = "Mirena",
  xlab = "log2 Fold Change",
  ylab = "-log10 adjusted p-value"
)

library(ggplot2)

volcano <- as.data.frame(res)
volcano$gene <- rownames(volcano)

# Quitar genes sin p-value
volcano <- volcano[!is.na(volcano$pvalue), ]

# Clasificación EXACTA según la tesis
volcano$category <- "no change"

volcano$category[
  !is.na(volcano$padj) &
    volcano$padj < 0.1 &
    volcano$log2FoldChange > 1.5
] <- "Up regulated"

volcano$category[
  !is.na(volcano$padj) &
    volcano$padj < 0.1 &
    volcano$log2FoldChange < -1.5
] <- "Down regulated"

# Eje Y
volcano$minuslog10p <- -log10(volcano$pvalue)

library(ggrepel)
volcano$gene_symbol <- sub(".*\\|", "", volcano$gene)
genes_label <- c(
  "ROBO1",
  "LRRN1",
  "MAP1B",
  "CCL2",
  "GDF15",
  "LDHA",
  "SCGB2A1",
  "SGIP1",
  "CAPN6",
  "IGFBP1"
)

ggplot(
  volcano,
  aes(
    x = log2FoldChange,
    y = minuslog10p,
    color = category
  )
) +
  geom_point(size = 1.5, alpha = 0.8) +
  
  geom_text_repel(
    data = subset(volcano, gene_symbol %in% genes_label),
    aes(label = gene_symbol),
    color = "black",
    size = 3,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Down regulated" = "#F8766D",
      "no change" = "#00BA38",
      "Up regulated" = "#619CFF"
    )
  ) +
  
  labs(
    x = "log2FoldChange",
    y = "-log(pvalue)",
    color = "category"
  ) +
  
  theme_gray()

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[order(res_df$pvalue), ]

head(
  res_df[, c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )],
  20
)


colData(dds)

table(dds$treatment)

colnames(counts(dds))

dim(counts(dds))

design(dds)

resultsNames(dds)

# Los genes más significativos de TU análisis
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[order(res_df$pvalue), ]

head(
  res_df[, c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
  20
)
# ¿Cuántos genes significativos tienes?
sum(res$padj < 0.1, na.rm = TRUE)

# Con los criterios escritos en la Figura 6
sum(
  res$padj < 0.1 &
    abs(res$log2FoldChange) > 1.5,
  na.rm = TRUE
)

sum(
  res$padj < 0.1 &
    res$log2FoldChange > 1.5,
  na.rm = TRUE
)

sum(
  res$padj < 0.1 &
    res$log2FoldChange < -1.5,
  na.rm = TRUE
)
