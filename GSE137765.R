#Resultados de Tesis de Montserrat Ramirez 

library(limma)
library(GEOquery)
library(tidyverse)

#GSE137765 - Diferencias en la expresión génica en sujetos bajo ningún tratamiento vs. tratamiento oral con Levonorgestrel
gse_137765 <- getGEO("GSE137765", GSEMatrix = TRUE)
pheno_gse_137765 <- pData(gse_137765[[1]])
conteo_gse_137765 <- exprs((gse_137765[[1]]))

#TEJIDO CERVICAL
md_pheno_gse137765 <- pheno_gse_137765 |> 
  select( c(geo_accession,
            `contraceptive used:ch1`,
             `age:ch1`,
            characteristics_ch1.3)) |> 
  rename(contraceptive = `contraceptive used:ch1`, 
         age = `age:ch1`) |> 
  filter(characteristics_ch1.3 == "tissue: cervical biopsy",
         !(contraceptive == "LNG-IUD"| contraceptive =="cu-IUD"))


md_pheno_gse137765 <- md_pheno_gse137765 |> 
  mutate(contraceptive = factor(contraceptive))

md_pheno_gse137765$contraceptive <- relevel(md_pheno_gse137765$contraceptive, ref = "control")

design_gse137765 <- model.matrix(~ 0 + md_pheno_gse137765$contraceptive)
design_gse137765

conteo_gse_137765 <- conteo_gse_137765[, md_pheno_gse137765$geo_accession]
all(colnames(conteo_gse_137765) == md_pheno_gse137765$geo_accession)

colnames(design_gse137765) <- make.names(colnames(design_gse137765))

fit_gse137765 <- lmFit(conteo_gse_137765, design_gse137765)


contraste_gse13776 <- makeContrasts(
  md_pheno_gse137765.contraceptivecontrol - md_pheno_gse137765.contraceptiveCOC,
  levels = design_gse137765)

fit2_gse137765 <- contrasts.fit(fit_gse137765, contraste_gse13776)
fit2_gse137765 <- eBayes(fit2_gse137765)

#resultados 
resultados_gse137765 <- topTable(
  fit2_gse137765,
  coef = 1,
  number = Inf,
  adjust.method = "BH"
)

#genes diferencilmente expresados 
deg_gse137765 <- resultados_gse137765[
  resultados_gse137765$adj.P.Val <0.05 &
    abs(resultados_gse137765$logFC) >=2.5,]

resultados_gse137765$significativo <- "No significativo"

resultados_gse137765$significativo[
  resultados_gse137765$adj.P.Val < 0.05 &
    resultados_gse137765$logFC >= 1
] <- "Sobreexpresado"

resultados_gse137765$significativo[
  resultados_gse137765$adj.P.Val < 0.05 &
    resultados_gse137765$logFC <= -1
] <- "Subexpresado"

ggplot(resultados_gse137765, aes(
  x = logFC,
  y = -log10(adj.P.Val),
  color = significativo
)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Clasificación"
  ) +
  theme_minimal()

#TEJIDO ENDOMETRIAL
conteo_gse_137765_endometrial <- exprs((gse_137765[[1]]))

md_pheno_gse137765_endometrial <- pheno_gse_137765 |> 
  select( c(geo_accession,
            `contraceptive used:ch1`,
            `age:ch1`,
            characteristics_ch1.3)) |> 
  rename(contraceptive = `contraceptive used:ch1`, 
         age = `age:ch1`) |> 
  filter(characteristics_ch1.3 == "tissue: endometrial biopsy",
         !(contraceptive == "LNG-IUD"| contraceptive =="cu-IUD"))


md_pheno_gse137765_endometrial <- md_pheno_gse137765_endometrial |> 
  mutate(contraceptive = factor(contraceptive))

md_pheno_gse137765_endometrial$contraceptive <- relevel(md_pheno_gse137765_endometrial$contraceptive, ref = "control")

design_gse137765_endometrial <- model.matrix(~ 0 + md_pheno_gse137765_endometrial$contraceptive)

conteo_gse_137765_endometrial <- conteo_gse_137765_endometrial[, md_pheno_gse137765_endometrial$geo_accession]
all(colnames(conteo_gse_137765_endometrial) == md_pheno_gse137765_endometrial$geo_accession)

colnames(design_gse137765_endometrial) <- make.names(colnames(design_gse137765_endometrial))

fit_gse137765_endometrial <- lmFit(conteo_gse_137765_endometrial, design_gse137765_endometrial)


contraste_gse13776_endometrial <- makeContrasts(
  md_pheno_gse137765_endometrial.contraceptivecontrol - md_pheno_gse137765_endometrial.contraceptiveCOC,
  levels = design_gse137765_endometrial)

fit2_gse137765_endometrial <- contrasts.fit(fit_gse137765_endometrial, contraste_gse13776_endometrial)
fit2_gse137765_endometrial <- eBayes(fit2_gse137765_endometrial)

#resultados 
resultados_gse137765_endometrial <- topTable(
  fit2_gse137765_endometrial,
  coef = 1,
  number = Inf,
  adjust.method = "BH"
)

#genes diferencilmente expresados 
deg_gse137765_endometrial <- resultados_gse137765_endometrial[
  resultados_gse137765_endometrial$adj.P.Val <0.05 &
    abs(resultados_gse137765_endometrial$logFC) >=2.5,]

resultados_gse137765_endometrial$significativo <- "No significativo"

resultados_gse137765_endometrial$significativo[
  resultados_gse137765_endometrial$adj.P.Val < 0.05 &
    resultados_gse137765_endometrial$logFC >= 1
] <- "Sobreexpresado"

resultados_gse137765_endometrial$significativo[
  resultados_gse137765_endometrial$adj.P.Val < 0.05 &
    resultados_gse137765_endometrial$logFC <= -1
] <- "Subexpresado"

ggplot(resultados_gse137765_endometrial, aes(
  x = logFC,
  y = -log10(adj.P.Val),
  color = significativo
)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Clasificación"
  ) +
  theme_minimal()
