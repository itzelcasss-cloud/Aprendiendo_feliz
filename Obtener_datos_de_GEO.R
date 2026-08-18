library(GEOquery)
library(tidyverse)

File1 <- getGEO("GSE137765", GSEMatrix = TRUE)

pheno_gse13776 <- pData(File1[[1]])

conteo <- exprs((File1[[1]]))
view(conteo)

#edad de la paciente que tenga el conteo mas alto del gen 7892501

pheno_gse13776 <- pheno_gse13776 |> 
  rename(age = `agee`)
  

pheno_gse13776[,colnames(pheno_gse13776) == "age"]

conteo[rownames(conteo) == "7892501",]
max(conteo[rownames(conteo) == "7892501",])

pheno_gse13776 |> 
  filter(geo_accession %in% "GSM4086937") |> 
  select(age)

pheno_gse13776$age[pheno_gse13776$geo_accession == colnames(conteo)[max(conteo[rownames(conteo) == "7892501",])]]

#tipo de biopsia del paciente con el menor número de conteos del gen en la posicion 1032

rownames(conteo)[1032]
conteo[1032,]
min(conteo[1032,])
min_1032 <- colnames(conteo)[(min(conteo[1032,]))]

pheno_gse13776$`tissue:ch1`[pheno_gse13776$geo_accession %in% "GSM4086934"]

pheno_gse13776["GSM4086934", "tissue:ch1"]

GSE188915 <- getGEO("GSE188915", GSEMatrix = TRUE)
pheno_gse188915 <- pData(GSE188915[[1]])
conteo_gse188915 <- exprs ((GSE188915[[1]]))


