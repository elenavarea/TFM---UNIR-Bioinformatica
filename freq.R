
rm(list=ls())
gc()


setwd("directorio")

archivos<-list.files()


for(i in 1:length(archivos)){
  
	
	print(paste("Estoy cargardo los datos del cromosoma", archivos[i]))
	
	load(archivos[i])
	
	
	func_maf<-function(x){
		
	  gt<-cut(x$snp,breaks=c(0,0.5,1.5,2),right=F,include.lowest=T)
		
		check<-as.numeric(sort(table(gt)))
		
		sum_allele<-c(c(check[1]*2)+c(check[2]*2)+c(check[3]*2))
		
		maf<-c(c(check[1]*2)+c(check[2]))/sum_allele*100
		
		return(maf)
		
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

	 res <- parLapply(cl,df_ds_list,func_maf)

	 print("Termina la funcion de analisis")
	
	 final<-Sys.time()
	
	 print(final-inicio)  

	 rm(inicio,final)

	 stopCluster(cl)
	
	 save(res,file=paste0("directorio/tabla_maf_",archivos[i]))
	
	 print(paste("Se ha guardado el cromosoma", archivos[i]))
	
	 rm(res)
	
	
}


#########################################################
##### UNION DE FRECUENCIAS A TABLA DE RESULTADOS ########
#########################################################
rm(list=ls())
gc()

archivos<-list.files()

load(archivos[1])

tabla_freq<-plyr::ldply(res)

colnames(tabla_freq)[1] <- "SNP"
colnames(tabla_freq)[2] <- "FREQ_ALT"

# Unimos las frecuencias a la tabla de resultados

for(i in 2:length(archivos)){
  
  print(i)
  
  load(archivos[i])
  
  tabla_unir<-plyr::ldply(res)
  
  colnames(tabla_unir)[1] <- "SNP"
  colnames(tabla_unir)[2] <- "FREQ_ALT"
  
  tabla_freq<-rbind(tabla_freq,tabla_unir)
  
  rm(tabla_unir)
  
}


load("directorio/res_total.RData")

resultado <- merge(tabla, tabla_freq, by = "SNP")

# Me quedo con los que tienen una frecuencia alélica mayor al 5%

resultado_filtro <- resultado[resultado$FREQ_ALT>=5,] 

save(resultado_filtro,file="directorio/res_freq_total.RData")


rm(tabla)

tabla <- resultado_filtro

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

# Saco el Manhattan plot
manhattan(tabla, chr="CHR", bp="POS", snp="SNP", p="P")

# Filtro para quedarme con los SNPs con p valor <= 1e-05
tabla_snps<-tabla[tabla$P<=1e-05,]

save(tabla_snps, file =  "directorio/tabla_snps.RData")


