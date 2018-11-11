source("functions.R")
source("load.R")




######PRIMERO, QUITAR MICROARREGLOS FEOS######


#Generamos boxplot
#boxplot(Ecoli_data[,2:51],las=2)
#boxplot(Ecoli_data[,52:101],las=2)
#boxplot(Ecoli_data[,102:151],las=2)
#boxplot(Ecoli_data[,152:201],las=2)
#boxplot(Ecoli_data[,202:214],las=2)
#Están especialmente feas las muestras 4, 211, 212. (únicas con cuartil_3 < 5)

#Generamos dendrograma para ver qué quitar
#data_dendr <- Ecoli_data[,-1] %>% `rownames<-`(1:nrow(Ecoli_data[,-1])) %>% `colnames<-`(1:ncol(Ecoli_data[,-1]))
#dd= dist2(log2(data_dendr))
#diag(dd)=0
#set.seed(1032)
#dd.row<-as.dendrogram(hclust(as.dist(dd)))
#row.ord<-order.dendrogram(dd.row) #organizacion
#legend=list(top=list(fun=dendrogramGrob,args=list(x=dd.row,side="top")))
#lp=levelplot(dd[row.ord,rev(row.ord)],xlab="", ylab="",legend=legend,  scales=list(x=list(rot=90)))
#lp
#No se ve mucho, pero por ahí hay un par muy poco correlacionadas con las demás. Miremos cuáles.
#which(apply(dd,MARGIN = 1,FUN = median)>0.22)
#La 210 también es relativamente diferente de las demás.

#Por lo anterior, quitamos las muestras 4, 210, 211, y 212
datos <- Ecoli_data %>% dplyr::select(-c(5,211,212,213))







####### AHORA, RESUMEN SONDAS GENES #####

#Nombres de filas y columnas
rownames(datos) <- datos[,1]
colnames(correspondencia)[1] <- colnames(datos)[1] <- "ProbeSet.id"
#Correspondencia (Resumen)
datos <- datos %>% full_join(correspondencia, by = "ProbeSet.id") %>% dplyr::select(-ProbeSet.id) %>%
  group_by(Gene.Symbol) %>% summarize_all(max) %>% drop_na(Gene.Symbol) %>% as.data.frame()
#Nombres de filas y columnas
rownames(datos) <- datos$Gene.Symbol
datos <- datos[,-1]







####### AHORA, MIRAMOS SI QUITAR GENES CON POCA VARIABILIDAD  #####

sds <- apply(datos,MARGIN = 1,FUN = sd)
#boxplot(sds)
means <- apply(datos, MARGIN = 1,FUN = mean)
#boxplot(means)
cvs <- sds/means
boxplot(cvs)

#ggplot(data=cvs,aes(x=.)) + geom_density(alpha=0.01)

check_pv<-which(cvs<quantile(cvs,0.001))
check_data_pv <- t(datos[check_pv,]) %>% as.data.frame() %>% 
  `colnames<-`(paste0("Pv",1:length(check_pv))) %>% mutate(Muestra = 1:dim(datos)[2]) %>%
  reshape2::melt(id.vars = "Muestra")
#ggplot(data=check_data_pv,aes(x=value,fill=variable,color=variable)) + geom_density(alpha =0.1) +
#  facet_grid(rows=vars(variable))

check_mv<-which(cvs<=quantile(cvs,0.152) & cvs>quantile(cvs,0.151))
check_data_mv <- t(datos[check_mv,]) %>% as.data.frame() %>% 
  `colnames<-`(paste0("Mv",1:length(check_mv))) %>% mutate(Muestra = 1:dim(datos)[2]) %>%
  reshape2::melt(id.vars = "Muestra")
#ggplot(data=check_data_mv,aes(x=value,fill=variable,color=variable)) + geom_density(alpha =0.1) +
#  facet_grid(rows=vars(variable))

check_lv<-which(cvs<=quantile(cvs,0.691) & cvs>quantile(cvs,0.69))
check_data_lv <- t(datos[check_lv,]) %>% as.data.frame() %>% 
  `colnames<-`(paste0("Lv",1:length(check_lv))) %>% mutate(Muestra = 1:dim(datos)[2]) %>%
  reshape2::melt(id.vars = "Muestra")
#ggplot(data=check_data_lv,aes(x=value,fill=variable,color=variable)) + geom_density(alpha =0.1) +
#  facet_grid(rows=vars(variable))


check_data_all<- rbind.data.frame(check_data_pv,check_data_mv,check_data_lv)
ggplot(data=check_data_all,aes(x=factor(0),y=value,fill=variable,color=variable)) +
  geom_violin(alpha =0.1) +
  facet_grid(cols=vars(variable))

#INTENTEMOS NO QUITANDO GENES... QUEDÉMONOS CON LOS 4267








######### ESTANDARIZACIÓN VSN ##########


#meanSdPlot(as.matrix(datos),ranks = T)
#Hay un poco de sobredispersión
datos1 <- justvsn(as.matrix(datos))
#meanSdPlot(datos1,ranks = T)
#Listo el pollo
datos1 %<>% as.data.frame()







###### EXPORTAR DATOS ########

datos <- datos1
saveobj("datos")

#cols <- c( "ALIAS")
#ensids <- keys(ecoli2.db)
#correspondencia <- AnnotationDbi::select(ecoli2.db, keys=ensids, columns=cols)


#sum(is.na(Ecoli_data))
