
rm(list=ls())
gc()

require(lme4)


setwd("directorio")

archivos<-list.files()


for(i in 1:length(archivos)){
  
	
	print(paste("Estoy cargardo los datos del cromosoma", archivos[i]))
	
	load(archivos[i])
	
	
	# Pega los datos de As en unas a cada base de datos (snp) del cromosoma en estudio
	
	print("Estoy pegando los datos de arsenico a las bases de datos")
	
	inicio<-Sys.time()
	for(j in 1:length(df_ds_list)){
		
		print(j)
		
		load("direcorio/datos_arsenico.RData")
		
		datos<-merge(df_ds_list[[j]],as_mcc[,c("newid","as_corr","tipocancer")],all.x=T,all.y=FALSE)
		datos<-datos[!is.na(datos$"as_corr"),]
		mediana_as<-median(log(datos$"as_corr"))
		datos$"outcome"<-ifelse(log(datos$"as_corr")>mediana_as,1,0)
	
		df_ds_list[[j]]<-datos
		
		rm(datos)
		
	}
	final<-Sys.time()
	
	print(final-inicio) 
	
	rm(inicio,final)
	
	
	func_analisis<-function(x){
		
		library(lme4)
		
		print(x$"nameSNP"[1])
		print("estoy con el modelo")
		model<-try(glmer(outcome ~ snp + edad + sexo + (1 | ccaa), data = x, family = binomial))
		
		if(class(model)=="try-error"){
			
			
			ci_model<-NA
			
			tabla_res<-data.frame(SNP=x$"nameSNP"[1],
			OR=NA, L_CI=NA , 
			U_CI=NA, 
			P=NA)
			

		}
		
		if(class(model)!="try-error"){
			
		
			print("estoy con el confint")

			ci_model<-confint(model,method="Wald")
	
			tabla_res<-data.frame(SNP=x$"nameSNP"[1],
			OR=exp(summary(model)$"coefficients"[2,1]), L_CI=exp(ci_model[3,1]) , 
			U_CI=exp(ci_model[3,2]), 
			P=summary(model)$"coefficients"[2,4])
		
		}
		
		return(tabla_res)
		
		rm(model,ci_model,tabla_res)
		
	}
	
	
	
	print("Estoy realizando los analisis")
	
	#######################
	# FORMATO PARALLEL
	#######################
	

	 library("parallel")

	 numWorkers <- detectCores()-1

	 cl <- makeCluster(numWorkers)

	 inicio<-Sys.time()
	
	 print("Empieza la funcion de analisis")

	 res <- parLapply(cl,df_ds_list,func_analisis)

	 print("Termina la funcion de analisis")
	
	 final<-Sys.time()
	
	 print(final-inicio)  

	 rm(inicio,final)

	 stopCluster(cl)
	
	 save(res,file=paste0("directorio/tabla_res_",archivos[i]))
	
	 print(paste("Se ha guardado el cromosoma", archivos[i]))
	
	 rm(res)
	
}

###########################################
########### UNION RESULTADOS ##############
###########################################

rm(list=ls())
gc()

require(lme4)


archivos<-list.files()

load(archivos[1])

tabla<-plyr::ldply(res)

info<-plyr::ldply(strsplit(tabla$SNP,split="[.]"))

tabla$"CHR"<-info$"V1"
tabla$"POS"<-info$"V2"
tabla$"REF"<-info$"V3"
tabla$"ALT"<-info$"V4"
tabla$".id"<-NULL

rm(info)

# Unimos las tablas de resultados

for(i in 2:length(archivos)){
	
	print(i)
	
	load(archivos[i])

	tabla_unir<-plyr::ldply(res)

	info<-plyr::ldply(strsplit(tabla_unir$SNP,split="[.]"))

	tabla_unir$"CHR"<-info$"V1"
	tabla_unir$"POS"<-info$"V2"
	tabla_unir$"REF"<-info$"V3"
	tabla_unir$"ALT"<-info$"V4"
	tabla_unir$".id"<-NULL
	
	tabla<-rbind(tabla,tabla_unir)
	
	rm(tabla_unir,info)
	
}

save(tabla,file="directorio/res_total.RData")


table(tabla$"CHR")

tabla$CHR<-as.numeric(tabla$CHR)

tabla$POS<-as.numeric(tabla$POS)

require("qqman")

# reviso los valores
range(tabla$P, na.rm = TRUE)
tabla <- tabla[!is.na(tabla$P) & is.finite(tabla$P), ]

range(tabla$OR, na.rm = TRUE)
tabla <- tabla[!is.na(tabla$OR) & is.finite(tabla$OR), ]
tabla <- tabla[tabla$OR != 0, ]
range(tabla$OR)




