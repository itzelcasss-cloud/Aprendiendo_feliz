library(limma)

File1 <- getGEO("GSE137765", GSEMatrix = TRUE)
pheno_gse13776 <- pData(File1[[1]])
conteo <- exprs((File1[[1]]))

metadata_ed_gse13776 <- pheno_gse13776 |> 
  select(geo_accession, `tissue:ch1`) |> 
  mutate(`tissue:ch1`= factor(`tissue:ch1`))

metadata_ed_gse13776$`tissue:ch1` <- relevel(metadata_ed_gse13776$`tissue:ch1`, ref = "cervical biopsy")

design <- model.matrix(~ 0 + metadata_ed_gse13776$`tissue:ch1`)
design

#asegurarse de que todos los nombres esten en el mismo orden
conteo <- conteo[, metadata_ed_gse13776$geo_accession]
all(colnames(conteo) == metadata_ed_gse13776$geo_accession)


fit <- lmFit(conteo, design)
fit

#asegurarse que los nombres de las columnas sean válidos
colnames(design) <- make.names(colnames(design))

contraste <- makeContrasts(
  metadata_ed_gse13776..tissue.ch1.cervical.biopsy - metadata_ed_gse13776..tissue.ch1.endometrial.biopsy,
  levels = design)

fit2 <- contrasts.fit(fit, contraste)
fit2 <- eBayes(fit2)

#resultados
resultados <- topTable(
  fit2, 
  coef = 1,
  number = Inf,
  adjust.method = "BH")

#seleccionar genes diferencialmente expresados
DEG <- resultados[
  resultados$adj.P.Val <0.001 &
    abs(resultados$logFC) >=2.5,]
DEG

library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)

#para saber el nombre de los genes 
mapIds(
  hugene10sttranscriptcluster.db,
  keys = rownames(DEG),
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

#agregar una columna que tenga el nombre del gen
DEG$symbol <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = rownames(DEG),
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

#ver los nombres de los genes SI HAY repertorio de GEO
annot <- getGEO("GPL6244", AnnotGPL = TRUE)
gpl <- Table(annot) |> 
  select("Gene symbol", ID)

#juntar los nombres de los genes (del repertorio que nos da GEO) con los que encontramos que estan diferenciados
gpl[gpl$ID %in% rownames(DEG),]

gpl <- gpl[gpl$ID %in% rownames(DEG),]

rownames(gpl) <- gpl$ID #para que el nombre de la fila no sea un numero, que sea el ID

merge(DEG,gpl,by=0)
