#EJERCICIOS

library(tidyverse)
names <- c("GSM531326", "GSM531329", "GSM531287", "GSM531289",
           "GSM531325", "GSM531290", "GSM531303", "GSM531307", "GSM531335",
           "GSM531330")

# Quedate con los pacientes que aparecen en names
metadata_gse21257|> 
  dplyr::filter(geo_accession %in% names)

metadata_gse21257[names, ]

# Quedate con los pacientes que aparecen en names y que son masculinos
metadata_gse21257 |> 
  filter(geo_accession %in% names, gender == "M")

metadata_gse21257[metadata_gse21257$geo_accession %in% names & metadata_gse21257$gender == "M", ]

# Quedate con pacientes que aparecen en names o que sean masculinos
metadata_gse21257 |> 
  filter(geo_accession %in% names | gender == "M")

metadata_gse21257[metadata_gse21257$geo_accession %in% names | metadata_gse21257$gender == "M", ]

# Obten los nombres de los pacientes que tienen huvos 2
metadata_gse21257 |> 
  rownames( group_by(huvos == 2) )

metadata_gse21257$geo_accession[metadata_gse21257$huvos == 2]

# Agrega una columna que sean puros numeros 1 y llamala column_1
metadata_gse21257 |> 
  mutate(column_1 = 1) 

metadata_gse21257$column_1 <- 1
metadata_gse21257
  
# Genera un histograma con la distribucion de edades
metadata_gse21257 |> 
  ggplot(aes(x=age_yr))+
  geom_histogram()

hist(metadata_gse21257$age_yr)  


# Obten la media, mediana, quantiles, maximo, minimo de la edad, extra si sacas desviacion estandar
summary(metadata_gse21257$age_yr)
min(metadata_gse21257$age_yr)
quantile(metadata_gse21257$age_yr)  
median(metadata_gse21257$age_yr)
mean(metadata_gse21257$age_yr)
max(metadata_gse21257$age_yr)
sd(metadata_gse21257$age_yr)

metadata_gse21257 |> 
  summarise(min = min(age_yr),
            q = quantile(age_yr, 0.75),  
            median = median(age_yr),
            mean = mean(age_yr),
            max = max(age_yr),
            sd = sd(age_yr),
            sum = sum(age_yr))

# Obten el subtipo histologico del paciente con mayor edad
metadata_gse21257[metadata_gse21257$age_yr == (max(metadata_gse21257$age_yr)), "hist_sub"]

metadata_gse21257 |> 
  filter(age_yr == max(age_yr)) |> 
  select(hist_sub)
