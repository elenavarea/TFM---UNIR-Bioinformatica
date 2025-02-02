

rm(list=ls())
gc()

library("SNPannot")

setwd("directorio")

load("tabla_snps.RData")


tabla_snps$CHR_POS <- sprintf("%s:%d", tabla_snps$CHR, tabla_snps$POS)

# Lista de SNPs
snp_list <- tabla_snps$CHR_POS

# Creo una lista para guardar los resultados
results <- list()

# Bucle para iterar sobre cada SNP
for (i in seq_along(snp_list)) {
  
  snp <- snp_list[i]
  
  print(i)
  
  # Ejecutar la función dbSNP_info ## r2=0.80
  res<- dbSNP_info(dat = snp, type = "pos", p = T, build = 37, r2 = 0.80,
                 pop = "EUR")
  
  # Guardo el resultado
  results[[snp]] <- res
  
}
rm(res)
# Mostrar los resultados
results

# lo paso a vector para luego pasarlo a dataframe
 term_vector <- sapply(results, function(x) if (!is.null(x$term)) x$term else NA)
 rsID_vector <- sapply(results, function(x) if (!is.null(x$rsID)) x$rsID else NA)
 gene_vector <- sapply(results, function(x) if (!is.null(x$gene)) x$gene else NA)
 alleles_vector <- sapply(results, function(x) if (!is.null(x$alleles)) x$alleles else NA)
 GRCh37_vector <- sapply(results, function(x) if (!is.null(x$GRCh37)) x$GRCh37 else NA)
 GRCh38_vector <- sapply(results, function(x) if (!is.null(x$GRCh38)) x$GRCh38 else NA)
 rs_ld_vector <- sapply(results, function(x) if (!is.null(x$rs_ld)) x$rs_ld else NA)

 results_df <- data.frame(
   term = term_vector,
   rsID = rsID_vector,
   gene = gene_vector,
   allelles = alleles_vector,
   GRCh37 = GRCh37_vector,
   GRCh38 = GRCh38_vector,
   rs_ld = rs_ld_vector,
   stringsAsFactors = FALSE
 )
 

 tabla_snps$log_OR <- log(tabla_snps$OR)
 tabla_snps$abs_log_OR <- abs(tabla_snps$log_OR)


library(dplyr)

# Dentro de los bloques de SNPs que están en DL, quedarme con los que tienen mayor beta absoluto
# Agrupar por cromosoma y seleccionar el SNP con el mayor beta 
  group_by(CHR) %>% # Agrupar por cromosoma
  slice_max(abs_log_OR, n = 1) %>%   # Seleccionar el SNP con el mayor beta absoluto
  ungroup()  
  
  save(snps_selected,file="directorio/snps_selected.RData")

#Extraer información del gen KDM4C  
genes_info <- snp_gene(dat = "KDM4C", type = "hgnc", p = T)
openxlsx::write.xlsx(genes_info,file="directorio/genes.info.xlsx")


