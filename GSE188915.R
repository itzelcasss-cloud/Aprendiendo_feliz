library(limma)
library(GEOquery)
library(tidyverse)
install.packages("affy")

#GSE188915 - Diferencias en la expresión génica en sujetos bajo ningún tratamiento vs. tratamiento con Levonorgestrel
gse_188915 <- getGEO("GSE188915", GSEMatrix = TRUE)
pheno_gse_188915 <- pData(gse_188915[[1]])

#porque no se pudo hacer con exprs - se tuvo que descargar a la computadora 
conteo_gse_188915 <- getGEOSuppFiles( "GSE188915" , fetch_files = TRUE, baseDir = "Downloads/", filter_regex = "^GSE188915_Organoid_Raw_Gene_Counts_ENSEMBL.csv.gz$") 
x <- list.files("Downloads/GSE188915/", full.names = TRUE, pattern = "\\.gz$")
conteo_gse_188915 <- read.delim(x, header = TRUE, sep = ",")

md_pheno_gse_188915 <- pheno_gse_188915 |> 
  select( c(geo_accession,
            source_name_ch1,
            characteristics_ch1,
            `age:ch1`))

md_pheno_gse_188915 <- md_pheno_gse_188915 |> 
  rename(source_name = source_name_ch1)
md_pheno_gse_188915 <- md_pheno_gse_188915 |> 
  rename(characteristics = characteristics_ch1)
md_pheno_gse_188915 <- md_pheno_gse_188915 |> 
  rename(age = `age:ch1`)

md_pheno_gse_188915 <- md_pheno_gse_188915 |> 
  mutate(characteristics= factor(characteristics))

md_pheno_gse_188915$characteristics <- relevel(md_pheno_gse_188915$characteristics, ref = "hormonal medication: None")

design <- model.matrix(~ 0 + md_pheno_gse_188915$characteristics)
design
