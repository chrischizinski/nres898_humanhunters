#####################################
#Read data                        ###
#####################################
#Spatial data set
datS<-read.table(file="datS.txt",sep="\t",header=T)
names(datS)

#Adjacency matrix for spatial INLA model
adj_red<-read.table(file="data/adj.red.txt",sep="\t")
rownames(adj_red)<-as.character(seq(1,424,1))

#Make log-density data for PCA
sd<-cbind(datS[,c("Region","fylkenr","idnr","densMoose16","densRed16","densRoe16","densreinhunters","densmoosehunters"         
                   ,"densredhunters" ,"densroehunters","densHuman17" ,"kdensrein17")])   

sd$logdenshumans<-log(sd$densHuman17)
sd$logdensMoose<-log(sd$densMoose16+0.001)
sd$logdensRed<-log(sd$densRed16+0.001)
sd$logdensRoe<-log(sd$densRoe16+0.01)
sd$logdensRein<-log(sd$kdensrein17+0.00001)
sd$logreinshunterD<-log(sd$densreinhunter+0.01)
sd$logmoosehunterD<-log(sd$densmoosehunter+0.01)
sd$logredhunterD<-log(sd$densredhunter+0.01)
sd$logroehunterD<-log(sd$densroehunter+0.01)

sd1<-sd[,c("Region","fylkenr","idnr","logdenshumans","logdensMoose","logdensRed","logdensRoe","logdensRein","logreinshunterD",
           "logmoosehunterD","logredhunterD","logroehunterD")]


#############################################################################
### Spatial Analysis                                                      ###
#############################################################################

######################
#1. PCA analysis
######################
library(ggplot2)
library(ggfortify)
library(ggrepel)

#Remove northern counties 
sd2<-sd1[!(sd1$fylkenr %in% c("18","19","20")),]

#help(princomp)
#Use data with log-densities
xx<-sd2[,-c(1:3)]
names(xx)

#Make shorter names
names(xx)<-c("Inhabitants","Moose","Red_deer","Roe_deer","Reindeer","Reindeer_hunters","Moose_hunters","Reddeer_hunters","Roedeer_hunters")
pc.cr <- princomp(xx, cor = TRUE)
bplot<-autoplot(pc.cr, data = sd2, colour = 'Region',size=1.5,loadings=T,
                loadings.colour = grey(0.5),
                loadings.label = TRUE, loadings.label.size = 4.0,loadings.label.colour='darkblue',
                loadings.label.repel=T)+ theme_bw()+
  theme(legend.justification=c(0.0,0.0), legend.position=c(0.0,0.0),legend.title=element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(values=c("pink","snow3","lightsalmon","lightblue"))

bplot$layers[[2]]$aes_params$size <- 1.0 #Line thickness
bplot$layers[[2]]$geom_params$arrow$length <- unit(6, units = "points")
bplot

#PCA_Loadings
round(pc.cr$loadings[,1:3],4)


###########################
#2. Pairwise correlations
###########################

#Pairwise correlations

library(boot)

dxy <- data.frame(cbind(V1,V2))
Brep = 10000
n=length(V1)

pearson <- function(d,i=c(1:n)){
  d2 <- d[i,]
  return(cor(d2$V1,d2$V2))
}

spearman <- function(d,i=c(1:n)){
  d2 <- d[i,]
  return(cor(d2$V1,d2$V2,method="spearman"))
}

#Pairwise correlation with bootstrap confidence interval (results given in supplementary tables)

#Example: correlation between log hunter density (red deer) and log red deer density in region West
sdat1<-datS[datS$Region=="West",]
V1 <- log(sdat1$densredhunter+1)
V2 <- log(sdat1$densRed16+0.01)

Brep = 10000
n=length(V1)

dxy <- data.frame(cbind(V1,V2))

bootcorr <- boot(data=dxy,statistic=pearson,R=Brep)
#bootcorr  
bb<-boot.ci(bootcorr,conf=.95)
cc1<-c(bb[2]$t0,bb[5]$basic[1,c(4,5)])
round(cc1,2) #Pearson correlation coefficient with lower and upper limit of bootstrap CI

bootcorrS <- boot(data=dxy,statistic=spearman,R=Brep)
#bootcorrS  
sbb<-boot.ci(bootcorrS,conf=.95)
sc1<-c(sbb[2]$t0,sbb[5]$basic[1,c(4,5)])
round(sc1,2) #Spearman correlation coefficient with lower and upper limit of bootstrap CI

###################################################
#3. Spatial regression models in INLA (BYM model)
###################################################

require(sp)
require(rgdal)
require(INLA)
# 
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

sdat<-datS

adj_red[-1,-1] -> adj_red_2
rownames(adj_red_2)<-as.character(seq(1,424,1))
colnames(adj_red_2)<-as.character(seq(1,424,1))
adj1<-as.matrix(adj_red_2, "dgTMatrix")
#unique(sdat$idnr) #idnr refers to the corresponding row/column number in adj1 matrix

g = inla.read.graph(adj1)
summary(g)

dim(adj_red)

#Data set without north
sdatR<-sdat[sdat$Region!="North",]
sdatR$Region<-factor(sdatR$Region,levels=c("West","South","East"))
table(sdatR$Region)

#All Norway
table(sdat$Region)

#Results for Table 1

#Offset for density
sdat$E=sdat$totareal
sdatR$E=sdatR$totareal

#A. Red deer hunters
fbym1.st <- redhunters17 ~  1 +
  scale(log(densHuman17)) +
  scale(log(densRed16+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbym1 = inla(fbym1.st, 
             data = sdatR, 
             control.inla = list(tolerance = 1e-9), 
             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
             family = "poisson",
             offset=log(E))

summary(mbym1)

modR1<-rbind(summary(mbym1)$fixed[,1:6],
             summary(mbym1)$hyperpar)
modR1

#B. Moose hunters
fbymm1.st <- moosehunters17 ~  1 +
  scale(log(densHuman17)) +
  scale(log(densMoose16+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbymm1 = inla(fbymm1.st, 
              data = sdat, 
              control.inla = list(tolerance = 1e-9), 
              control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
              family = "poisson",
              offset=log(E))
summary(mbymm1)

modM1<-rbind(summary(mbymm1)$fixed[,1:6],
             summary(mbymm1)$hyperpar)
modM1

#C. Reindeer hunters
fbymrr1.st <- reinhunters17 ~  1 +
  scale(log(densHuman17)) +
  scale(log(kdensrein17+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbymrr1 = inla(fbymrr1.st, 
               data = sdatR, 
               control.inla = list(tolerance = 1e-9), 
               control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
               family = "poisson",
               offset=log(E))

summary(mbymrr1)
modRr1<-rbind(summary(mbymrr1)$fixed[,1:6],
             summary(mbymrr1)$hyperpar)
modRr1


#D. Roe deer hunters
fbymroe1.st <- roehunters ~  1 +
  scale(log(densHuman17)) +
  scale(log(densRoe16+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbymroe1 = inla(fbymroe1.st, 
               data = sdatR, 
               control.inla = list(tolerance = 1e-9), 
               control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
               family = "poisson",
               offset=log(E))

summary(mbymroe1)
modRo1<-rbind(summary(mbymroe1)$fixed[,1:6],
             summary(mbymroe1)$hyperpar)
modRo1

#Table 1
modR1
modM1
modRr1
modRo1


#Results for Table 2

#Offset for incidence
sdatR$E1<-sdatR$humans2017
sdat$E1<-sdat$humans2017

#A. Red deer hunter incidence
fbym.st <- redhunters17 ~  1 +
  scale(log(densRed16+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbym1i = inla(fbym.st, 
             data = sdatR, 
             control.inla = list(tolerance = 1e-9), 
             control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
             family = "poisson",
             offset=log(E1))
summary(mbym1i)


modRei<-rbind(summary(mbym1i)$fixed[,1:6],
               summary(mbym1i)$hyperpar)
modRei


#B. Moose hunters
#Alternative model formulation with bym (bym2 will not run)
fbymm.st <- moosehunters17 ~  1 +
  scale(log(densMoose16+0.01))*Region+
  f(idnr, model = "bym",graph=g)  

mbymm1i = inla(fbymm.st, 
               data = sdat, 
               control.inla = list(tolerance = 1e-9), 
               control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
               family = "poisson",
               offset=log(E1))
summary(mbymm1i)

#Vary by region and only spatial effect
fbm.st <- moosehunters17 ~  1 +
  scale(log(densMoose16+0.01))*Region+
  f(idnr, model = "besag",graph=g)  

mbm1i = inla(fbm.st, 
               data = sdat, 
               control.inla = list(tolerance = 1e-9), 
               control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
               family = "poisson",
               offset=log(E1))
summary(mbm1i)

mbymm1i$dic$dic
mbm1i$dic$dic
#select the region model 

modM1i<-rbind(summary(mbm1i)$fixed[,1:6],
              summary(mbm1i)$hyperpar)
modM1i


#C. Reindeer hunters
fbymrr.st <- reinhunters17 ~  1 +
  scale(log(kdensrein17+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbymrri = inla(fbymrr.st, 
               data = sdatR, 
               control.inla = list(tolerance = 1e-9), 
               control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
               family = "poisson",
               offset=log(E1))

summary(mbymrri)
modRri<-rbind(summary(mbymrri)$fixed[,1:6],
              summary(mbymrri)$hyperpar)
modRri


#D. Roe deer hunters
fbymroe.st <- roehunters17 ~  1 +
  scale(log(densRoe16+0.01))+
  f(idnr, model = "bym2",graph=g)  

mbymroe1i = inla(fbymroe.st, 
                data = sdatR, 
                control.inla = list(tolerance = 1e-9), 
                control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                family = "poisson",
                offset=log(E1))

summary(mbymroe1i)

fbymroe2.st <- roehunters ~  1 +
  scale(log(densRoe16+0.01))*Region +
  f(idnr, model = "bym2",graph=g)  

mbymroe2i = inla(fbymroe2.st, 
                 data = sdatR, 
                 control.inla = list(tolerance = 1e-9), 
                 control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                 family = "poisson",
                 offset=log(E1))

summary(mbymroe2i)

mbymroe1i$dic$dic
mbymroe2i$dic$dic #Better model


modRo2<-rbind(summary(mbymroe2i)$fixed[,1:6],
              summary(mbymroe2i)$hyperpar)
modRo2

#Table 2
modRei
modM1i
modRri
modRo2



