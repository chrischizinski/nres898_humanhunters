library(nlme)
library(MuMIn) #providing function for AICc (AIC corrected for small sample sizes)

tdat<-read.table(file="tdat.txt",header=T,sep="\t")
#table(tdat$Region,tdat$County)

###########
#West #####
###########
ddW<-tdat[tdat$Region=="West",]
unique(ddW$County)
vf=varIdent(form=~1|County)

#Delta AIC for table 3
fitW<-gls(R_redhunters~(R_redharvest.West_tl1)+relevel(as.factor(County),ref="14"),method="ML",dat=ddW)
summary(fitw)
fitW1<-gls(R_redhunters~(R_redharvest.West_tl1)+relevel(as.factor(County),ref="14"),weights=vf,method="ML",dat=ddW)
fitW_w<-gls(R_redhunters~relevel(as.factor(County),ref="14"),method="ML",dat=ddW)
fitW_f<-gls(R_redhunters~(R_redharvest.West_tl1),method="ML",dat=ddW)
fitWc<-gls(R_redhunters~(R_redharvest.West_tl1)+relevel(as.factor(County),ref="14"),
           correlation=corAR1(form=~1|as.factor(County)),method="ML",dat=ddW)

round(AICc(fitW,fitW1,fitW_w,fitW_f,fitWc)-AICc(fitW),1)

#Model estimates for table 3 by REML
tfitW<-gls(R_redhunters~(R_redharvest.West_tl1)+relevel(as.factor(County),ref="14"),method="REML",dat=ddW)
summary(tfitW)


############
#East ######
############
ddE<-tdat[tdat$Region=="East",]
unique(ddE$County)
vf=varIdent(form=~1|County)

#Delta AIC for table 3
fitE1<-gls(R_redhunters ~ R_redharvest.West_tl1+R_redharvest.East_tl1+R_moosehunters +as.factor(County),weights=vf,method="ML",dat=ddE)
summary(fit14)
fitE1_Mh<-gls(R_redhunters ~ R_redharvest.West_tl1 +R_redharvest.Region_tl1+as.factor(County),weights=vf,method="ML",dat=ddE)
fitE1_W<-gls(R_redhunters ~ R_redharvest.Region_tl1+R_moosehunters+as.factor(County),weights=vf,method="ML",dat=ddE)
fitE1_M<-gls(R_redhunters ~ R_redharvest.West_tl1+R_moosehunters +as.factor(County),weights=vf,method="ML",dat=ddE)
fitE1_f<-gls(R_redhunters ~ R_redharvest.West_tl1 +R_redharvest.Region_tl1+R_moosehunters,weights=vf,method="ML",dat=ddE)
fitE1c<-gls(R_redhunters ~ R_redharvest.West_tl1 +R_redharvest.Region_tl1+R_moosehunters+as.factor(County),weights=vf,
            correlation=corAR1(form=~1|as.factor(County) ),method="ML",dat=ddE)

round(AICc(fitE1,fitE1_W,fitE1_M,fitE1_Mh,fitE1_f,fitE1c)-AICc(fitE1),1)

#Model estimates for table 3 by REML
tfitE1<-gls(R_redhunters ~ R_redharvest.West_tl1+R_redharvest.Region_tl1+R_moosehunters+as.factor(County),weights=vf,method="ML",dat=ddE)
summary(tfitE1)


############
#South #####
############
ddS<-tdat[tdat$Region=="South",]
unique(ddS$County)
vf=varIdent(form=~1|County)

#Delta AIC for table 3
fitS1<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+I(year-2010)+I((year-2010)^2)+as.factor(County),method="ML",weights=vf,dat=ddS)
summary(fitS1)
fitS1_W<-gls(R_redhunters~(R_redharvest.Region_tl1)+I(year-2010)+I((year-2010)^2)+as.factor(County),method="ML",weights=vf,dat=ddS)
fitS1_R<-gls(R_redhunters~(R_redharvest.West_tl1)+I(year-2010)+I((year-2010)^2)+as.factor(County),method="ML",weights=vf,dat=ddS)
fitS1_Y2<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+I(year-2010)+as.factor(County),method="ML",weights=vf,dat=ddS)
fitS1_Y<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+as.factor(County),method="ML",weights=vf,dat=ddS)
fitS1_B<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+I(year-2010)+I((year-2010)^2),method="ML",weights=vf,dat=ddS)
fitS1c<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+I(year-2010)+I((year-2010)^2)+as.factor(County),
            method="ML",correlation=corAR1(form=~1|as.factor(County)),weights=vf,dat=ddS)

round(AICc(fitS1,fitS1_W,fitS1_R,fitS1_Y,fitS1_Y2,fitS1_B,fitS1c)-AICc(fitS1),1)

#Model estimates for table 3 by REML
tfitS1<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+I(year-2010)+I((year-2010)^2)+as.factor(County),method="REML",weights=vf,dat=ddS)
summary(tfitS1)


############
#North######
############
ddN<-tdat[tdat$Region=="Mid",]
unique(ddN$County)
vf1=varIdent(form=~1|County)

#Delta AIC for table 3
fitN1<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+factor(County),weights=vf1,method="ML",dat=ddN)
summary(fitN1)
fitN1_W<-gls(R_redhunters~(R_redharvest.Region_tl1)+factor(County),method="ML",weights=vf1,dat=ddN)
fitN1_R<-gls(R_redhunters~(R_redharvest.West_tl1)+factor(County),method="ML",weights=vf1,dat=ddN)
fitN1_f<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1),method="ML",weights=vf1,dat=ddN)
fitN1c<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+factor(County),weights=vf1,
            correlation=corAR1(form=~1|as.factor(County)),method="ML",dat=ddN)

round(AICc(fitN1,fitN1_W,fitN1_R,fitN1_f,fitN1c)-AICc(fitN1),1)

#Model estimates for table 3 by REML
tfitN1<-gls(R_redhunters~(R_redharvest.West_tl1)+ (R_redharvest.Region_tl1)+factor(County),weights=vf1,method="REML",dat=ddN)
summary(tfitN1)

###########
