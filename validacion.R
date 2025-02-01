#####################################
# VALIDACION LEAVE ONE OUT
#####################################


indice_snps<-2:5
datos_list<-list()

for(i in 1:dim(df)[1]){
  
  print(i)
  snps<-c(indice_snps)
  datos_train<-df[-i,]
  datos_test<-df[i,]
  
  # snps vs genes
  
  tabla<-data.frame(SNP=names(datos_train)[indice_snps])
  tabla$indice<-seq(1,dim(tabla)[1],1)
  
  # Estimadores de asociacion para cada SNP en modelos genicos (datos_train)
  
  res<-list()
  
  for(j in 1:length(snps)){
    
    datos_train$"snp"<-datos_train[,snps[j]]
    
    model_train<-glmer(outcome ~ snp + edad + sexo + (1 | ccaa), data = datos_train, family = binomial)
    
    res[[j]]<-summary(model_train)$"coefficients"[2,1]
    
  }
  
  res_final<-plyr::ldply(res)
  
  tabla<-cbind(tabla,res_final)
  names(tabla)[3]<-"coef"
  
  datos_w<-datos_test 
  for(p in 1:length(snps)){
    
    peso<-tabla$coef[p]
    datos_w[,snps[p]]<-datos_w[,snps[p]]*peso
    
  }	
  
  datos_w$ps_loo<-apply(datos_w[,indice_snps],1,function(x)sum(x))
  
  datos_list[[i]]<-datos_w
  
  
}

datos_finales<-plyr::ldply(datos_list)

datos_finales$ps_loo_check<-apply(datos_finales[,indice_snps],1,function(x)sum(x))

sum(datos_finales$ps_loo_check==datos_finales$ps_loo)

plot(datos_finales$ps_as,datos_finales$ps_loo)


save(datos_finales,file="datos_con_ps_loo.RData")


########################################################
# Analisis de los datos de validacion
########################################################

rm(list=ls())
gc()

load("datos_con_ps_loo.RData")

df<-datos_finales

# ps_as en continua
sd(df$ps_loo) # 0.4475101

require(lme4)

df$pgs_loo_sd<-df$ps_loo/sd(df$ps_loo) 

model1<-(glmer(outcome ~ pgs_loo_sd + edad + sexo + (1 | ccaa), 
               data = df, family = binomial))

ci_model1<-confint(model1,method="Wald")

tabla1_res<-data.frame(OR=exp(summary(model1)$"coefficients"[2,1]), L_CI=exp(ci_model1[3,1]) , 
                       U_CI=exp(ci_model1[3,2]), 
                       P=summary(model1)$"coefficients"[2,4])

# ps_as en cuartiles

cuartiles <- quantile(df$ps_loo)


df$ps_loo_q <- cut(df$ps_loo,breaks =  quantile(df$ps_loo), 
                   include.lowest=T ,    
                   right = TRUE)

model2<-(glmer(outcome ~ ps_loo_q + edad + sexo + (1 | ccaa), data = df, family = binomial))

ci_model2<-confint(model2,method="Wald")

tabla2_res<-data.frame(OR=exp(summary(model2)$"coefficients"[2:4,1]), L_CI=exp(ci_model2[3:5,1]) , 
                       U_CI=exp(ci_model2[3:5,2]), 
                       P=summary(model2)$"coefficients"[2:4,4])


tabla_unida<-rbind(tabla1_res,tabla2_res)

openxlsx::write.xlsx(tabla_unida,file="tabla_unida_res_PGS_loo_as.xlsx")

hist(df$ps_loo,main="",ylab="Frecuencia",xlab="PGS As 1-out")