
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
