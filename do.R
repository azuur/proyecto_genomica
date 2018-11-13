source("functions.R")

readobj("datos")



CORR_TIME  <- system.time(CORR_NET  <- cor(t(datos)) %>% `diag<-`(0)); CORR_TIME
#MIM_P_NET <- build.mim(datos,"pearson") 
#MIM_TIME   <- system.time(MIM_E_NET <- build.mim(t(datos),estimator = "mi.mm",disc = "equalwidth"))
MIM_TIME   <- system.time(MIM_E_NET <- fastGeneMI::get.mim.MM(t(datos),discretisation = "equalwidth")); MIM_TIME
ARACN_TIME <- system.time(ARACN_NET <- minet::aracne(MIM_E_NET)); ARACN_TIME
#ARACN_TIME1 <- system.time(ARACN_NET1 <- fastGeneMI::(MIM_E_NET))
CLR_TIME   <- system.time(CLR_NET   <- clr(MIM_E_NET)); CLR_TIME
MRNET_TIME <- system.time(MRNET_NET <- mrnet(MIM_E_NET)); MRNET_TIME

ls() %>% savelist()


readobj("MIM_E_NET")

GENIE_TIME <- system.time(GENIE_NET <- GENIE3(as.matrix(datos),verbose = T,nCores=3))
ls() %>% savelist()

#EL LAMBDA ESTÁ FEO... EN EL SENTIDO DE SER MUY DISTINTO AL VALOR DE 1 ORIGINAL
NARRO_TIME <- system.time(NARRO_NET <- narromi(t(datos),MIM = MIM_E_NET,verbose = T,lambda=0.3,parallel = T))
ls() %>% savelist()

TIGRE_TIME <- system.time(TIGRE_NET <- tigress(t(datos),verbose = T)); TIGRE_TIME
ls() %>% savelist()






