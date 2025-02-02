
rm(list=ls())
gc()

require(lme4)


# Extracción de dosis de los SNPs seleccionados

archivos<-list.files()

# cromosoma 4
print(paste("Estoy cargardo los datos del cromosoma", archivos[17]))

load(archivos[17])

df_snp1 <- data.frame(
  newid = df_ds_list$"4.114783874.A.G"$newid, 
  "4.114783874.A.G" = df_ds_list$"4.114783874.A.G"$snp
)

rm(df_ds_list)

# cromosoma 5
print(paste("Estoy cargardo los datos del cromosoma", archivos[18]))

load(archivos[18])

df_snp2 <- data.frame(
  newid = df_ds_list$"5.31067563.A.G"$newid, 
  "5.31067563.A.G" = df_ds_list$"5.31067563.A.G"$snp
)

rm(df_ds_list)

# cromosoma 9
print(paste("Estoy cargardo los datos del cromosoma", archivos[22]))

load(archivos[22])

df_snp3 <- data.frame(
  newid = df_ds_list$"9.7147976.G.A"$newid, 
  "9.7147976.G.A" = df_ds_list$"9.7147976.G.A"$snp
)

rm(df_ds_list)

# cromosoma 19
print(paste("Estoy cargardo los datos del cromosoma", archivos[11]))

load(archivos[11])

df_snp4 <- data.frame(
  newid = df_ds_list$"19.15955743.A.T"$newid, 
  "19.15955743.A.T" = df_ds_list$"19.15955743.A.T"$snp
)

rm(df_ds_list)

# Combinar los data frames
merged_df <- merge(df_snp1, df_snp2, by = "newid", all = TRUE)  
merged_df <- merge(merged_df, df_snp3, by = "newid", all = TRUE)  
merged_df <- merge(merged_df, df_snp4, by = "newid", all = TRUE) 

ps_as_tabla <- merged_df

rm(merged_df)


load("directorio/snps_selected.RData")

# Calcular el peso de cada SNP (log(OR))
snps_selected$log_OR <- log(snps_selected$OR)

# Cálculo del PGSas de SNPs seleccionados: multiplico las dosis(0 ,1 ,2) por los pesos(log_OR) y los sumo
ps_as_tabla$ps_as <- (ps_as_tabla$"X4.114783874.A.G" * snps_selected$log_OR[1]) + (ps_as_tabla$"X5.31067563.A.G" * snps_selected$log_OR[2]) + (ps_as_tabla$"X9.7147976.G.A" * snps_selected$log_OR[3]) + (ps_as_tabla$"X19.15955743.A.T" * snps_selected$log_OR[4])



##########################
# SACO LA MEDIANA DE AS EN UÑAS
#########################
load("directorio/arsenico.RData")

datos<-merge(ps_as_tabla,as_mcc[,c("newid","as_corr","tipocancer", "edad", "sexo", "ccaa")],all.x=T,all.y=FALSE)
datos<-datos[!is.na(datos$"as_corr"),]
mediana_as_log<-median(log(datos$"as_corr"))
datos$"outcome"<-ifelse(log(datos$"as_corr")>mediana_as_log,1,0)

mediana_as<-median((datos$"as_corr")/1000) # está en ppb (microg/Kg) y lo paso a microg/g dividiendo entre 1000

ps_as_tabla<-datos

rm(datos)

####################################################
## MODELIZACIÓN NIVELES DE AS CON BETAS GENERALES
####################################################

df<-ps_as_tabla

table(df$sexo)
prop.table(table(df$sexo))*100

mean(df$edad)
t.test(df$edad)

data.frame(table(df$ccaa))

data.frame(prop.table(table(df$ccaa))*100)

median(df$as_corr)/1000

hist(df$ps_as,main="",ylab="Frecuencia",xlab="PGS")


# ps_as en continua
sd(df$ps_as) 


df$pgs_sd<-df$ps_as/sd(df$ps_as) 

model1<-(glmer(outcome ~ pgs_sd + edad + sexo + (1 | ccaa), 
data = df, family = binomial))

ci_model1<-confint(model1,method="Wald")

tabla1_res<-data.frame(OR=exp(summary(model1)$"coefficients"[2,1]), L_CI=exp(ci_model1[3,1]) , 
                       U_CI=exp(ci_model1[3,2]), 
                       P=summary(model1)$"coefficients"[2,4])

# ps_as en cuartiles

cuartiles <- quantile(ps_as_tabla$ps_as)

													 
ps_as_tabla$ps_as_q <- cut(ps_as_tabla$ps_as, 
													                breaks = quantile(ps_as_tabla$ps_as), 
													                       include.lowest=T ,    
													                  right = TRUE)

model2<-(glmer(outcome ~ ps_as_q + edad + sexo + (1 | ccaa), data = ps_as_tabla, family = binomial))

ci_model2<-confint(model2,method="Wald")

tabla2_res<-data.frame(OR=exp(summary(model2)$"coefficients"[2:4,1]), L_CI=exp(ci_model2[3:5,1]) , 
                       U_CI=exp(ci_model2[3:5,2]), 
                       P=summary(model2)$"coefficients"[2:4,4])


tabla_unida<-rbind(tabla1_res,tabla2_res)

openxlsx::write.xlsx(tabla_unida,file="tabla_unida_res_PGS_as.xlsx")


######################################################
# MODELIZACION CANCER DE MAMA CON BETAS GENERALES
######################################################


rm(list=ls())
gc()

setwd("directorio")


load("datos.RData")


# Cálculo del PGSas para muestras caso-control de cáncer de mama: multiplico las dosis(0 ,1 ,2) por los pesos(log_OR)
ds_bc$cc_ps_as <- (ds_bc$"4:114783874:A:G" * snps_selected$log_OR[1]) + (ds_bc$"5:31067563:A:G" * snps_selected$log_OR[2]) + (ds_bc$"9:7147976:G:A" * snps_selected$log_OR[3]) + (ds_bc$"19:15955743:A:T" * snps_selected$log_OR[4])


length(intersect(bc_mcc$"newid",ds_bc$"newid"))
bc_mcc <- merge(bc_mcc, ds_bc[, c("newid", "cc_ps_as")], by = "newid")
bc_mcc <- merge(bc_mcc, casom_as[, c("newid", "as_corr")], by = "newid")
dim(bc_mcc)

df<-bc_mcc


table(df$casom,exclude=NULL)

################################
#######     CON PGS
################################

# Exportar los gráficos a un archivo PNG
png("res/hist_cc.png", width = 800, height = 600)

par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.1)

# Crear el primer histograma
hist(df$cc_ps_as[df$casom %in% "Caso_mama"], 
     main = "Casos", 
     ylab = "Frecuencia", 
     xlab = "PGS", 
     xlim = c(-1.5, 1.8))

# Crear el segundo histograma
hist(df$cc_ps_as[df$casom %in% "Control"], 
     main = "Controles", 
     ylab = "Frecuencia", 
     xlab = "PGS", 
     xlim = c(-1.5, 1.8))

dev.off()


# cc_ps_as sd

sd(df$cc_ps_as) 

df$pgs_sd<-df$cc_ps_as/sd(df$cc_ps_as) 

model1<-(glmer(casom ~ pgs_sd + edad + (1 | ccaa), 
data = df, family = binomial))

ci_model1<-confint(model1,method="Wald")

tabla1_res<-data.frame(OR=exp(summary(model1)$"coefficients"[2,1]), L_CI=exp(ci_model1[3,1]) , 
                       U_CI=exp(ci_model1[3,2]), 
                       P=summary(model1)$"coefficients"[2,4])


# ps_as en cuartiles

cuartiles <- quantile(df$cc_ps_as)
											 
df$ps_casom_q <- cut(df$cc_ps_as, breaks = quantile(df$cc_ps_as), 
include.lowest=T ,    
right = TRUE)

table(df$ps_casom_q,exclude=NULL)

model2<-(glmer(casom ~ ps_casom_q + edad + (1 | ccaa), data = df, family = binomial))

ci_model2<-confint(model2,method="Wald")

tabla2_res<-data.frame(OR=exp(summary(model2)$"coefficients"[2:4,1]), L_CI=exp(ci_model2[3:5,1]) , 
											                        U_CI=exp(ci_model2[3:5,2]), 
											                        P=summary(model2)$"coefficients"[2:4,4])


tabla_unida<-rbind(tabla1_res,tabla2_res)

openxlsx::write.xlsx(tabla_unida,file="tabla_unida_res_PGS_casom.xlsx")

#############################
###   CON AS_CORR: As en uñas
#############################
df$as_corr <- df$as_corr/1000 # lo paso a microg/g
# Exportar los gráficos a un archivo PNG
png("res/hist_as_corr.png", width = 800, height = 600)

par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.1)

# Crear el primer histograma
hist(df$as_corr[df$casom %in% "Caso_mama"], 
     main = "Casos", 
     ylab = "Frecuencia", 
     xlab = "As uñas (µg/g)")

# Crear el segundo histograma
hist(df$as_corr[df$casom %in% "Control"], 
     main = "Controles", 
     ylab = "Frecuencia", 
     xlab = "As uñas (µg/g)")

dev.off()

# as_corr sd

sd(df$as_corr) 


df$as_corr_sd<-df$as_corr/sd(df$as_corr) 

model1<-(glmer(casom ~ as_corr_sd + edad + (1 | ccaa), 
               data = df, family = binomial))

ci_model1<-confint(model1,method="Wald")

tabla1_res<-data.frame(OR=exp(summary(model1)$"coefficients"[2,1]), L_CI=exp(ci_model1[3,1]) , 
                       U_CI=exp(ci_model1[3,2]), 
                       P=summary(model1)$"coefficients"[2,4])


# ps_as en cuartiles

cuartiles <- quantile(df$as_corr)

df$as_corr_q <- cut(df$as_corr, breaks = quantile(df$as_corr), 
                     include.lowest=T ,    
                     right = TRUE)

table(df$as_corr_q,exclude=NULL)

model2<-(glmer(casom ~ as_corr_q + edad + (1 | ccaa), data = df, family = binomial))

ci_model2<-confint(model2,method="Wald")

tabla2_res<-data.frame(OR=exp(summary(model2)$"coefficients"[2:4,1]), L_CI=exp(ci_model2[3:5,1]) , 
                       U_CI=exp(ci_model2[3:5,2]), 
                       P=summary(model2)$"coefficients"[2:4,4])


tabla_unida<-rbind(tabla1_res,tabla2_res)

openxlsx::write.xlsx(tabla_unida,file="tabla_unida_res_as_corr_casom.xlsx")





