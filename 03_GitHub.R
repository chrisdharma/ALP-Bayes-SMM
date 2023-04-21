load("cchsSNall.Rdata")

###Bayesian analysis
###Note the SN models were run with 10000 iterations, whereas the CCHS were run with 1000 iterations

library(dplyr)
library(R2jags)
library(fastDummies)
library(tidyverse)

set.seed(134)

cchsSNall <- cchsSNall %>% mutate(employ = case_when(employ == 1 ~ 1,
                                                     employ %in% c(2,3,4) ~ 2,
                                                     employ %in% c(5,6) ~ 3))

cchsSNall <- dummy_cols(cchsSNall, select_columns = "province", remove_first_dummy = TRUE)

cchsSNall2 <- cchsSNall %>% 
  dplyr::select(-c(HIV_test12m,STI_test12m,lastcondom_use,insurance,marijuana,year))

cchsSNall_msm <- cchsSNall2 %>% 
  filter(msm == 1) 

alldat_msm <- cchsSNall_msm %>% drop_na(age_grp,income,education,orient,msm,illicit_drug,have_pcp,smoke,employ,province,rural,race,existing_mh,immigration,alcohol)
alldat <- cchsSNall2 %>% drop_na(age_grp,income,education,orient,msm,illicit_drug,have_pcp,smoke,employ,province,rural,race,existing_mh,immigration,alcohol)

SN <- alldat %>% dplyr::filter(survey == 1)
cchsall <- alldat %>% dplyr::filter(survey == 0)

NumSN <- nrow(SN)
ageSN <- SN$age_grp
incomeSN <- SN$income
smokeSN <- SN$smoke
eduSN <- SN$education
DiscSN <- SN$disclosure
drugSN <- SN$illicit_drug
pcpSN <- SN$have_pcp
employSN <- SN$employ
provinceSN <- SN$province
raceSN <- SN$race
immigSN <- SN$immigration + 1
alcoholSN <- SN$alcohol
orientSN <- SN$orient
existingSN <- SN$existing_mh
provCentralSN <- SN$province_Central
provPrairiesSN <- SN$province_Prairies
provWestSN <- SN$province_WestCoast
provNorthSN <- SN$province_Northern
DeprSN <- SN$depress_binary
ruralSN <- SN$rural

catimmig<-length(unique(immigSN))
catemploy<-length(unique(employSN))
catedu<-length(unique(eduSN))

model_codeSN<-"
model{
  for (iSN in 1:NumSN) {
    DiscSN[iSN] ~ dbern(pdiscSN[iSN])
    logit(pdiscSN[iSN]) <- DI0 + DIage*ageSN[iSN] + DIexisting*existingSN[iSN] + DIWest*provWestSN[iSN] + DINorth*provNorthSN[iSN] + DIPrairies*provPrairiesSN[iSN] + DICentral*provCentralSN[iSN] + DIincome*incomeSN[iSN] +
    DIalcohol*alcoholSN[iSN] + DIrace*raceSN[iSN] + DIpcp*pcpSN[iSN] + DIsmoke*smokeSN[iSN] + DIdrug*drugSN[iSN] + DIrural*ruralSN[iSN] + DIedu[eduSN[iSN]] + DIimmig[immigSN[iSN]] + DIemploy[employSN[iSN]]
  }
  for(ic in 1:catimmig){
    DIimmig.temp[ic] ~ dnorm(0,0.001)
   }
   DIimmig <- DIimmig.temp - mean(DIimmig.temp)
   for(ib in 1:catemploy){
     DIemploy.temp[ib] ~ dnorm(0,0.001)
   }
   DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  for(id in 1:catedu){
     DIedu.temp[id] ~ dnorm(0,0.001)
   }
   DIedu <- DIedu.temp - mean(DIedu.temp)
#Priors
DI0 ~ dnorm(0,.001)
DIage~dnorm(0, 0.001)
DIexisting~dnorm(0, 0.001)
DIWest~dnorm(0, 0.001)
DINorth~dnorm(0, 0.001)
DIPrairies~dnorm(0, 0.001)
DICentral~dnorm(0, 0.001)
DIincome~dnorm(0, 0.001)
DIalcohol~dnorm(0, 0.001)
DIrace~dnorm(0, 0.001)
DIpcp~dnorm(0, 0.001)
DIsmoke~dnorm(0, 0.001)
DIdrug~dnorm(0, 0.001)
DIrural~dnorm(0, 0.001)
}
"

##The best initial value will be the ones from the regression model itself
DI.mod<-glm(DiscSN ~ ageSN + existingSN + provWestSN + provNorthSN + provPrairiesSN + provCentralSN + incomeSN +
              alcoholSN + raceSN + pcpSN + smokeSN + raceSN + drugSN + ruralSN + as.factor(eduSN) + as.factor(immigSN) + as.factor(employSN),family=binomial)

init.vals <- function(){
  list(
    DI0=rnorm(1,mean=summary(DI.mod)$coefficients[1,1],sd=summary(DI.mod)$coefficients[1,2]),
    DIage=rnorm(1,mean=summary(DI.mod)$coefficients[2,1],sd=summary(DI.mod)$coefficients[2,2]),
    DIexisting=rnorm(1,mean=summary(DI.mod)$coefficients[3,1],sd=summary(DI.mod)$coefficients[3,2]),
    DIWest=rnorm(1,mean=summary(DI.mod)$coefficients[4,1],sd=summary(DI.mod)$coefficients[4,2]),
    DINorth=rnorm(1,mean=summary(DI.mod)$coefficients[5,1],sd=summary(DI.mod)$coefficients[5,2]),
    DIPrairies=rnorm(1,mean=summary(DI.mod)$coefficients[6,1],sd=summary(DI.mod)$coefficients[6,2]),
    DICentral=rnorm(1,mean=summary(DI.mod)$coefficients[7,1],sd=summary(DI.mod)$coefficients[7,2]),
    DIincome=rnorm(1,mean=summary(DI.mod)$coefficients[8,1],sd=summary(DI.mod)$coefficients[8,2]),
    DIalcohol=rnorm(1,mean=summary(DI.mod)$coefficients[9,1],sd=summary(DI.mod)$coefficients[9,2]),
    DIrace=rnorm(1,mean=summary(DI.mod)$coefficients[10,1],sd=summary(DI.mod)$coefficients[10,2]),
    DIpcp=rnorm(1,mean=summary(DI.mod)$coefficients[11,1],sd=summary(DI.mod)$coefficients[11,2]),
    DIsmoke=rnorm(1,mean=summary(DI.mod)$coefficients[12,1],sd=summary(DI.mod)$coefficients[12,2]),
    DIdrug=rnorm(1,mean=summary(DI.mod)$coefficients[13,1],sd=summary(DI.mod)$coefficients[13,2]),
    DIrural=rnorm(1,mean=summary(DI.mod)$coefficients[14,1],sd=summary(DI.mod)$coefficients[14,2]))
}
init.vals()

jagsSN.dat=list("NumSN","DiscSN","ageSN","existingSN","drugSN","immigSN","employSN","ruralSN","eduSN",
                "pcpSN","smokeSN","raceSN","alcoholSN","incomeSN","provCentralSN","provPrairiesSN","provWestSN","provNorthSN",
                "catedu","catimmig","catemploy")

paramsSN=c("DI0","DIage","DIexisting","DIWest","DINorth","DIPrairies","DICentral","DIincome","DIalcohol","DIrace",
           "DIpcp","DIsmoke","DIdrug","DIrural","DIedu","DIemploy","DIimmig")

start.time <- Sys.time()
jagsSN=jags(data=jagsSN.dat, inits=init.vals,parameters.to.save=paramsSN, model.file=textConnection(model_codeSN),
            n.chains=3, n.iter=10000, n.burnin=100, n.thin=1,jags.seed=123)
print(jagsSN)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

summary=jagsSN$BUGSoutput$summary
DI0_m<- summary["DI0",]["mean"]
DI0_p<- 1/(summary["DI0",]["sd"])^2
DIage_m<- summary["DIage",]["mean"]
DIage_p<- 1/(summary["DIage",]["sd"])^2
DIexisting_m <- summary["DIexisting",]["mean"]
DIexisting_p <- 1/(summary["DIexisting",]["sd"])^2
DIWest_m <-summary["DIWest",]["mean"]
DIWest_p <- 1/(summary["DIWest",]["sd"])^2
DINorth_m <-summary["DINorth",]["mean"]
DINorth_p <- 1/(summary["DINorth",]["sd"])^2
DIPrairies_m <-summary["DIPrairies",]["mean"]
DIPrairies_p <- 1/(summary["DIPrairies",]["sd"])^2
DICentral_m <-summary["DICentral",]["mean"]
DICentral_p <- 1/(summary["DICentral",]["sd"])^2
DIincome_m <-summary["DIincome",]["mean"]
DIincome_p <- 1/(summary["DIincome",]["sd"])^2
DIalcohol_m <-summary["DIalcohol",]["mean"]
DIalcohol_p <- 1/(summary["DIalcohol",]["sd"])^2
DIrace_m <-summary["DIrace",]["mean"]
DIrace_p <- 1/(summary["DIrace",]["sd"])^2
DIpcp_m <-summary["DIpcp",]["mean"]
DIpcp_p <- 1/(summary["DIpcp",]["sd"])^2
DIsmoke_m <-summary["DIsmoke",]["mean"]
DIsmoke_p <- 1/(summary["DIsmoke",]["sd"])^2
DIdrug_m <-summary["DIdrug",]["mean"]
DIdrug_p <- 1/(summary["DIdrug",]["sd"])^2
DIrural_m <-summary["DIrural",]["mean"]
DIrural_p <- 1/(summary["DIrural",]["sd"])^2
DIedu_m1 <-summary["DIedu[1]",]["mean"]
DIedu_p1 <- 1/(summary["DIedu[1]",]["sd"])^2
DIedu_m2 <-summary["DIedu[2]",]["mean"]
DIedu_p2 <- 1/(summary["DIedu[2]",]["sd"])^2
DIedu_m3 <-summary["DIedu[3]",]["mean"]
DIedu_p3 <- 1/(summary["DIedu[3]",]["sd"])^2
DIemploy_m1 <-summary["DIemploy[1]",]["mean"]
DIemploy_p1 <- 1/(summary["DIemploy[1]",]["sd"])^2
DIemploy_m2 <-summary["DIemploy[2]",]["mean"]
DIemploy_p2 <- 1/(summary["DIemploy[2]",]["sd"])^2
DIemploy_m3 <-summary["DIemploy[3]",]["mean"]
DIemploy_p3 <- 1/(summary["DIemploy[3]",]["sd"])^2
DIimmig_m1 <-summary["DIimmig[1]",]["mean"]
DIimmig_p1 <- 1/(summary["DIimmig[1]",]["sd"])^2
DIimmig_m2 <-summary["DIimmig[2]",]["mean"]
DIimmig_p2 <- 1/(summary["DIimmig[2]",]["sd"])^2
DIimmig_m3 <-summary["DIimmig[3]",]["mean"]
DIimmig_p3 <- 1/(summary["DIimmig[3]",]["sd"])^2
# Time difference of 17.42 mins

cchsall_dep <- cchsall %>% drop_na(depress_binary)
cchsall_consult <- cchsall %>% drop_na(consult_mh)

NumCC <- nrow(cchsall_dep)
ageCC <- cchsall_dep$age_grp
incomeCC <- cchsall_dep$income
smokeCC <- cchsall_dep$smoke
eduCC <- cchsall_dep$education
drugCC <- cchsall_dep$illicit_drug
pcpCC <- cchsall_dep$have_pcp
employCC <- cchsall_dep$employ
provinceCC <- cchsall_dep$province
raceCC <- cchsall_dep$race
immigCC <- cchsall_dep$immigration + 1
alcoholCC <- cchsall_dep$alcohol
orientCC <- cchsall_dep$orient
existingCC <- cchsall_dep$existing_mh
provCentralCC <- cchsall_dep$province_Central
provPrairiesCC <- cchsall_dep$province_Prairies
provWestCC <- cchsall_dep$province_WestCoast
provNorthCC <- cchsall_dep$province_Northern
DeprCC <- cchsall_dep$depress_binary
consultCC <- cchsall_dep$consult_mh
ruralCC <- cchsall_dep$rural
wtsCC <- cchsall_dep$WTS_M_rescaled
MSMrep <- cchsall_dep$msm

model_codeCCHS_depress<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  DeprCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-DeprCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(DeprCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","DeprCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself
r.mod<-glm(cchsall_dep$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_depress),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

jagsCCHS_wts.mcmc<-as.mcmc(jagsCCHS_wts) 
summary(jagsCCHS_wts.mcmc)

library(coda)
geweke.diag(jagsCCHS_wts.mcmc)
# Z score of most less than |2|
geweke.plot(jagsCCHS_wts.mcmc)

gelman.diag(jagsCCHS_wts.mcmc)
gelman.plot(jagsCCHS_wts.mcmc)

if(T){
  SArray= jagsCCHS_wts$BUGSoutput$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){ 
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    xnul=locator(1)    
    acf( SArray[,1,nn], main=nn)  #note: this is only for 1st chain
    xnul=locator(1)
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
    xnul=locator(1)
  }
} 

NumCC <- nrow(cchsall_consult)
ageCC <- cchsall_consult$age_grp
incomeCC <- cchsall_consult$income
smokeCC <- cchsall_consult$smoke
eduCC <- cchsall_consult$education
drugCC <- cchsall_consult$illicit_drug
pcpCC <- cchsall_consult$have_pcp
employCC <- cchsall_consult$employ
provinceCC <- cchsall_consult$province
raceCC <- cchsall_consult$race
immigCC <- cchsall_consult$immigration + 1
alcoholCC <- cchsall_consult$alcohol
orientCC <- cchsall_consult$orient
existingCC <- cchsall_consult$existing_mh
provCentralCC <- cchsall_consult$province_Central
provPrairiesCC <- cchsall_consult$province_Prairies
provWestCC <- cchsall_consult$province_WestCoast
provNorthCC <- cchsall_consult$province_Northern
consultCC <- cchsall_consult$consult_mh
ruralCC <- cchsall_consult$rural
wtsCC <- cchsall_consult$WTS_M_rescaled
MSMrep <- cchsall_consult$msm

model_codeCCHS_consult<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  consultCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-consultCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(consultCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","consultCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_consult$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_consult),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

jagsCCHS_wts.mcmc<-as.mcmc(jagsCCHS_wts) 
summary(jagsCCHS_wts.mcmc)

library(coda)
geweke.diag(jagsCCHS_wts.mcmc)
geweke.plot(jagsCCHS_wts.mcmc)

gelman.diag(jagsCCHS_wts.mcmc)
gelman.plot(jagsCCHS_wts.mcmc)

if(T){
  SArray= jagsCCHS_wts$BUGSoutput$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){ 
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    xnul=locator(1)    
    acf( SArray[,1,nn], main=nn)  #note: this is only for 1st chain
    xnul=locator(1)
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
    xnul=locator(1)
  }
} 

####Next, rerun everything for gays only
SN_gay <- SN %>% dplyr::filter(orient == 2)
cchsall_gay_het <- cchsall %>% dplyr::filter(orient %in% c(1,2))

NumSN <- nrow(SN_gay)
ageSN <- SN_gay$age_grp
incomeSN <- SN_gay$income
smokeSN <- SN_gay$smoke
eduSN <- SN_gay$education
DiscSN <- SN_gay$disclosure
drugSN <- SN_gay$illicit_drug
pcpSN <- SN_gay$have_pcp
employSN <- SN_gay$employ
provinceSN <- SN_gay$province
raceSN <- SN_gay$race
immigSN <- SN_gay$immigration + 1
alcoholSN <- SN_gay$alcohol
orientSN <- SN_gay$orient
existingSN <- SN_gay$existing_mh
provCentralSN <- SN_gay$province_Central
provPrairiesSN <- SN_gay$province_Prairies
provWestSN <- SN_gay$province_WestCoast
provNorthSN <- SN_gay$province_Northern
DeprSN <- SN_gay$depress_binary
ruralSN <- SN_gay$rural

catimmig<-length(unique(immigSN))
catemploy<-length(unique(employSN))
catedu<-length(unique(eduSN))

model_codeSN<-"
model{
  for (iSN in 1:NumSN) {
    DiscSN[iSN] ~ dbern(pdiscSN[iSN])
    logit(pdiscSN[iSN]) <- DI0 + DIage*ageSN[iSN] + DIexisting*existingSN[iSN] + DIWest*provWestSN[iSN] + DINorth*provNorthSN[iSN] + DIPrairies*provPrairiesSN[iSN] + DICentral*provCentralSN[iSN] + DIincome*incomeSN[iSN] +
    DIalcohol*alcoholSN[iSN] + DIrace*raceSN[iSN] + DIpcp*pcpSN[iSN] + DIsmoke*smokeSN[iSN] + DIdrug*drugSN[iSN] + DIrural*ruralSN[iSN] + DIedu[eduSN[iSN]] + DIimmig[immigSN[iSN]] + DIemploy[employSN[iSN]]
  }
  for(ic in 1:catimmig){
    DIimmig.temp[ic] ~ dnorm(0,0.001)
   }
   DIimmig <- DIimmig.temp - mean(DIimmig.temp)
   for(ib in 1:catemploy){
     DIemploy.temp[ib] ~ dnorm(0,0.001)
   }
   DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  for(id in 1:catedu){
     DIedu.temp[id] ~ dnorm(0,0.001)
   }
   DIedu <- DIedu.temp - mean(DIedu.temp)
#Priors
DI0 ~ dnorm(0,.001)
DIage~dnorm(0, 0.001)
DIexisting~dnorm(0, 0.001)
DIWest~dnorm(0, 0.001)
DINorth~dnorm(0, 0.001)
DIPrairies~dnorm(0, 0.001)
DICentral~dnorm(0, 0.001)
DIincome~dnorm(0, 0.001)
DIalcohol~dnorm(0, 0.001)
DIrace~dnorm(0, 0.001)
DIpcp~dnorm(0, 0.001)
DIsmoke~dnorm(0, 0.001)
DIdrug~dnorm(0, 0.001)
DIrural~dnorm(0, 0.001)
}
"

##The best initial value will be the ones from the regression model itself
DI.mod<-glm(DiscSN ~ ageSN + existingSN + provWestSN + provNorthSN + provPrairiesSN + provCentralSN + incomeSN +
              alcoholSN + raceSN + pcpSN + smokeSN + raceSN + drugSN + ruralSN + as.factor(eduSN) + as.factor(immigSN) + as.factor(employSN),family=binomial)

init.vals <- function(){
  list(
    DI0=rnorm(1,mean=summary(DI.mod)$coefficients[1,1],sd=summary(DI.mod)$coefficients[1,2]),
    DIage=rnorm(1,mean=summary(DI.mod)$coefficients[2,1],sd=summary(DI.mod)$coefficients[2,2]),
    DIexisting=rnorm(1,mean=summary(DI.mod)$coefficients[3,1],sd=summary(DI.mod)$coefficients[3,2]),
    DIWest=rnorm(1,mean=summary(DI.mod)$coefficients[4,1],sd=summary(DI.mod)$coefficients[4,2]),
    DINorth=rnorm(1,mean=summary(DI.mod)$coefficients[5,1],sd=summary(DI.mod)$coefficients[5,2]),
    DIPrairies=rnorm(1,mean=summary(DI.mod)$coefficients[6,1],sd=summary(DI.mod)$coefficients[6,2]),
    DICentral=rnorm(1,mean=summary(DI.mod)$coefficients[7,1],sd=summary(DI.mod)$coefficients[7,2]),
    DIincome=rnorm(1,mean=summary(DI.mod)$coefficients[8,1],sd=summary(DI.mod)$coefficients[8,2]),
    DIalcohol=rnorm(1,mean=summary(DI.mod)$coefficients[9,1],sd=summary(DI.mod)$coefficients[9,2]),
    DIrace=rnorm(1,mean=summary(DI.mod)$coefficients[10,1],sd=summary(DI.mod)$coefficients[10,2]),
    DIpcp=rnorm(1,mean=summary(DI.mod)$coefficients[11,1],sd=summary(DI.mod)$coefficients[11,2]),
    DIsmoke=rnorm(1,mean=summary(DI.mod)$coefficients[12,1],sd=summary(DI.mod)$coefficients[12,2]),
    DIdrug=rnorm(1,mean=summary(DI.mod)$coefficients[13,1],sd=summary(DI.mod)$coefficients[13,2]),
    DIrural=rnorm(1,mean=summary(DI.mod)$coefficients[14,1],sd=summary(DI.mod)$coefficients[14,2]))
}
init.vals()

jagsSN.dat=list("NumSN","DiscSN","ageSN","existingSN","drugSN","immigSN","employSN","ruralSN","eduSN",
                "pcpSN","smokeSN","raceSN","alcoholSN","incomeSN","provCentralSN","provPrairiesSN","provWestSN","provNorthSN",
                "catedu","catimmig","catemploy")

paramsSN=c("DI0","DIage","DIexisting","DIWest","DINorth","DIPrairies","DICentral","DIincome","DIalcohol","DIrace",
           "DIpcp","DIsmoke","DIdrug","DIrural","DIedu","DIemploy","DIimmig")

start.time <- Sys.time()
jagsSN=jags(data=jagsSN.dat, inits=init.vals,parameters.to.save=paramsSN, model.file=textConnection(model_codeSN),
            n.chains=3, n.iter=10000, n.burnin=100, n.thin=1,jags.seed=123)
print(jagsSN)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

summary=jagsSN$BUGSoutput$summary
DI0_m<- summary["DI0",]["mean"]
DI0_p<- 1/(summary["DI0",]["sd"])^2
DIage_m<- summary["DIage",]["mean"]
DIage_p<- 1/(summary["DIage",]["sd"])^2
DIexisting_m <- summary["DIexisting",]["mean"]
DIexisting_p <- 1/(summary["DIexisting",]["sd"])^2
DIWest_m <-summary["DIWest",]["mean"]
DIWest_p <- 1/(summary["DIWest",]["sd"])^2
DINorth_m <-summary["DINorth",]["mean"]
DINorth_p <- 1/(summary["DINorth",]["sd"])^2
DIPrairies_m <-summary["DIPrairies",]["mean"]
DIPrairies_p <- 1/(summary["DIPrairies",]["sd"])^2
DICentral_m <-summary["DICentral",]["mean"]
DICentral_p <- 1/(summary["DICentral",]["sd"])^2
DIincome_m <-summary["DIincome",]["mean"]
DIincome_p <- 1/(summary["DIincome",]["sd"])^2
DIalcohol_m <-summary["DIalcohol",]["mean"]
DIalcohol_p <- 1/(summary["DIalcohol",]["sd"])^2
DIrace_m <-summary["DIrace",]["mean"]
DIrace_p <- 1/(summary["DIrace",]["sd"])^2
DIpcp_m <-summary["DIpcp",]["mean"]
DIpcp_p <- 1/(summary["DIpcp",]["sd"])^2
DIsmoke_m <-summary["DIsmoke",]["mean"]
DIsmoke_p <- 1/(summary["DIsmoke",]["sd"])^2
DIdrug_m <-summary["DIdrug",]["mean"]
DIdrug_p <- 1/(summary["DIdrug",]["sd"])^2
DIrural_m <-summary["DIrural",]["mean"]
DIrural_p <- 1/(summary["DIrural",]["sd"])^2
DIedu_m1 <-summary["DIedu[1]",]["mean"]
DIedu_p1 <- 1/(summary["DIedu[1]",]["sd"])^2
DIedu_m2 <-summary["DIedu[2]",]["mean"]
DIedu_p2 <- 1/(summary["DIedu[2]",]["sd"])^2
DIedu_m3 <-summary["DIedu[3]",]["mean"]
DIedu_p3 <- 1/(summary["DIedu[3]",]["sd"])^2
DIemploy_m1 <-summary["DIemploy[1]",]["mean"]
DIemploy_p1 <- 1/(summary["DIemploy[1]",]["sd"])^2
DIemploy_m2 <-summary["DIemploy[2]",]["mean"]
DIemploy_p2 <- 1/(summary["DIemploy[2]",]["sd"])^2
DIemploy_m3 <-summary["DIemploy[3]",]["mean"]
DIemploy_p3 <- 1/(summary["DIemploy[3]",]["sd"])^2
DIimmig_m1 <-summary["DIimmig[1]",]["mean"]
DIimmig_p1 <- 1/(summary["DIimmig[1]",]["sd"])^2
DIimmig_m2 <-summary["DIimmig[2]",]["mean"]
DIimmig_p2 <- 1/(summary["DIimmig[2]",]["sd"])^2
DIimmig_m3 <-summary["DIimmig[3]",]["mean"]
DIimmig_p3 <- 1/(summary["DIimmig[3]",]["sd"])^2

#Dep and consult gays in CCHS
cchsall_dep <- cchsall_gay_het %>% drop_na(depress_binary)
cchsall_consult <- cchsall_gay_het %>% drop_na(consult_mh)

NumCC <- nrow(cchsall_dep)
ageCC <- cchsall_dep$age_grp
incomeCC <- cchsall_dep$income
smokeCC <- cchsall_dep$smoke
eduCC <- cchsall_dep$education
drugCC <- cchsall_dep$illicit_drug
pcpCC <- cchsall_dep$have_pcp
employCC <- cchsall_dep$employ
provinceCC <- cchsall_dep$province
raceCC <- cchsall_dep$race
immigCC <- cchsall_dep$immigration + 1
alcoholCC <- cchsall_dep$alcohol
orientCC <- cchsall_dep$orient
existingCC <- cchsall_dep$existing_mh
provCentralCC <- cchsall_dep$province_Central
provPrairiesCC <- cchsall_dep$province_Prairies
provWestCC <- cchsall_dep$province_WestCoast
provNorthCC <- cchsall_dep$province_Northern
DeprCC <- cchsall_dep$depress_binary
consultCC <- cchsall_dep$consult_mh
ruralCC <- cchsall_dep$rural
wtsCC <- cchsall_dep$WTS_M_rescaled
MSMrep <- cchsall_dep$msm


model_codeCCHS_depress<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  DeprCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-DeprCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(DeprCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","DeprCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_dep$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_depress),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

if(T){
  SArray= jagsCCHS_wts$BUGSoutput$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){ 
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    xnul=locator(1)    
    acf( SArray[,1,nn], main=nn)  #note: this is only for 1st chain
    xnul=locator(1)
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
    xnul=locator(1)
  }
} 

NumCC <- nrow(cchsall_consult)
ageCC <- cchsall_consult$age_grp
incomeCC <- cchsall_consult$income
smokeCC <- cchsall_consult$smoke
eduCC <- cchsall_consult$education
drugCC <- cchsall_consult$illicit_drug
pcpCC <- cchsall_consult$have_pcp
employCC <- cchsall_consult$employ
provinceCC <- cchsall_consult$province
raceCC <- cchsall_consult$race
immigCC <- cchsall_consult$immigration + 1
alcoholCC <- cchsall_consult$alcohol
orientCC <- cchsall_consult$orient
existingCC <- cchsall_consult$existing_mh
provCentralCC <- cchsall_consult$province_Central
provPrairiesCC <- cchsall_consult$province_Prairies
provWestCC <- cchsall_consult$province_WestCoast
provNorthCC <- cchsall_consult$province_Northern
consultCC <- cchsall_consult$consult_mh
ruralCC <- cchsall_consult$rural
wtsCC <- cchsall_consult$WTS_M_rescaled
MSMrep <- cchsall_consult$msm

model_codeCCHS_consult<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  consultCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-consultCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(consultCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","consultCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_consult$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_consult),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

####Next, rerun everything for bis only
SN_bi <- SN %>% dplyr::filter(orient == 3)
cchsall_bi_het <- cchsall %>% dplyr::filter(orient %in% c(1,3))

NumSN <- nrow(SN_bi)
ageSN <- SN_bi$age_grp
incomeSN <- SN_bi$income
smokeSN <- SN_bi$smoke
eduSN <- SN_bi$education
DiscSN <- SN_bi$disclosure
drugSN <- SN_bi$illicit_drug
pcpSN <- SN_bi$have_pcp
employSN <- SN_bi$employ
provinceSN <- SN_bi$province
raceSN <- SN_bi$race
immigSN <- SN_bi$immigration + 1
alcoholSN <- SN_bi$alcohol
orientSN <- SN_bi$orient
existingSN <- SN_bi$existing_mh
provCentralSN <- SN_bi$province_Central
provPrairiesSN <- SN_bi$province_Prairies
provWestSN <- SN_bi$province_WestCoast
provNorthSN <- SN_bi$province_Northern
DeprSN <- SN_bi$depress_binary
ruralSN <- SN_bi$rural

catimmig<-length(unique(immigSN))
catemploy<-length(unique(employSN))
catedu<-length(unique(eduSN))

model_codeSN<-"
model{
  for (iSN in 1:NumSN) {
    DiscSN[iSN] ~ dbern(pdiscSN[iSN])
    logit(pdiscSN[iSN]) <- DI0 + DIage*ageSN[iSN] + DIexisting*existingSN[iSN] + DIWest*provWestSN[iSN] + DINorth*provNorthSN[iSN] + DIPrairies*provPrairiesSN[iSN] + DICentral*provCentralSN[iSN] + DIincome*incomeSN[iSN] +
    DIalcohol*alcoholSN[iSN] + DIrace*raceSN[iSN] + DIpcp*pcpSN[iSN] + DIsmoke*smokeSN[iSN] + DIdrug*drugSN[iSN] + DIrural*ruralSN[iSN] + DIedu[eduSN[iSN]] + DIimmig[immigSN[iSN]] + DIemploy[employSN[iSN]]
  }
  for(ic in 1:catimmig){
    DIimmig.temp[ic] ~ dnorm(0,0.001)
   }
   DIimmig <- DIimmig.temp - mean(DIimmig.temp)
   for(ib in 1:catemploy){
     DIemploy.temp[ib] ~ dnorm(0,0.001)
   }
   DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  for(id in 1:catedu){
     DIedu.temp[id] ~ dnorm(0,0.001)
   }
   DIedu <- DIedu.temp - mean(DIedu.temp)
#Priors
DI0 ~ dnorm(0,.001)
DIage~dnorm(0, 0.001)
DIexisting~dnorm(0, 0.001)
DIWest~dnorm(0, 0.001)
DINorth~dnorm(0, 0.001)
DIPrairies~dnorm(0, 0.001)
DICentral~dnorm(0, 0.001)
DIincome~dnorm(0, 0.001)
DIalcohol~dnorm(0, 0.001)
DIrace~dnorm(0, 0.001)
DIpcp~dnorm(0, 0.001)
DIsmoke~dnorm(0, 0.001)
DIdrug~dnorm(0, 0.001)
DIrural~dnorm(0, 0.001)
}
"

##The best initial value will be the ones from the regression model itself
DI.mod<-glm(DiscSN ~ ageSN + existingSN + provWestSN + provNorthSN + provPrairiesSN + provCentralSN + incomeSN +
              alcoholSN + raceSN + pcpSN + smokeSN + raceSN + drugSN + ruralSN + as.factor(eduSN) + as.factor(immigSN) + as.factor(employSN),family=binomial)

init.vals <- function(){
  list(
    DI0=rnorm(1,mean=summary(DI.mod)$coefficients[1,1],sd=summary(DI.mod)$coefficients[1,2]),
    DIage=rnorm(1,mean=summary(DI.mod)$coefficients[2,1],sd=summary(DI.mod)$coefficients[2,2]),
    DIexisting=rnorm(1,mean=summary(DI.mod)$coefficients[3,1],sd=summary(DI.mod)$coefficients[3,2]),
    DIWest=rnorm(1,mean=summary(DI.mod)$coefficients[4,1],sd=summary(DI.mod)$coefficients[4,2]),
    DINorth=rnorm(1,mean=summary(DI.mod)$coefficients[5,1],sd=summary(DI.mod)$coefficients[5,2]),
    DIPrairies=rnorm(1,mean=summary(DI.mod)$coefficients[6,1],sd=summary(DI.mod)$coefficients[6,2]),
    DICentral=rnorm(1,mean=summary(DI.mod)$coefficients[7,1],sd=summary(DI.mod)$coefficients[7,2]),
    DIincome=rnorm(1,mean=summary(DI.mod)$coefficients[8,1],sd=summary(DI.mod)$coefficients[8,2]),
    DIalcohol=rnorm(1,mean=summary(DI.mod)$coefficients[9,1],sd=summary(DI.mod)$coefficients[9,2]),
    DIrace=rnorm(1,mean=summary(DI.mod)$coefficients[10,1],sd=summary(DI.mod)$coefficients[10,2]),
    DIpcp=rnorm(1,mean=summary(DI.mod)$coefficients[11,1],sd=summary(DI.mod)$coefficients[11,2]),
    DIsmoke=rnorm(1,mean=summary(DI.mod)$coefficients[12,1],sd=summary(DI.mod)$coefficients[12,2]),
    DIdrug=rnorm(1,mean=summary(DI.mod)$coefficients[13,1],sd=summary(DI.mod)$coefficients[13,2]),
    DIrural=rnorm(1,mean=summary(DI.mod)$coefficients[14,1],sd=summary(DI.mod)$coefficients[14,2]))
}
init.vals()

jagsSN.dat=list("NumSN","DiscSN","ageSN","existingSN","drugSN","immigSN","employSN","ruralSN","eduSN",
                "pcpSN","smokeSN","raceSN","alcoholSN","incomeSN","provCentralSN","provPrairiesSN","provWestSN","provNorthSN",
                "catedu","catimmig","catemploy")

paramsSN=c("DI0","DIage","DIexisting","DIWest","DINorth","DIPrairies","DICentral","DIincome","DIalcohol","DIrace",
           "DIpcp","DIsmoke","DIdrug","DIrural","DIedu","DIemploy","DIimmig")

start.time <- Sys.time()
jagsSN=jags(data=jagsSN.dat, inits=init.vals,parameters.to.save=paramsSN, model.file=textConnection(model_codeSN),
            n.chains=3, n.iter=10000, n.burnin=100, n.thin=1,jags.seed=123)
print(jagsSN)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

summary=jagsSN$BUGSoutput$summary
DI0_m<- summary["DI0",]["mean"]
DI0_p<- 1/(summary["DI0",]["sd"])^2
DIage_m<- summary["DIage",]["mean"]
DIage_p<- 1/(summary["DIage",]["sd"])^2
DIexisting_m <- summary["DIexisting",]["mean"]
DIexisting_p <- 1/(summary["DIexisting",]["sd"])^2
DIWest_m <-summary["DIWest",]["mean"]
DIWest_p <- 1/(summary["DIWest",]["sd"])^2
DINorth_m <-summary["DINorth",]["mean"]
DINorth_p <- 1/(summary["DINorth",]["sd"])^2
DIPrairies_m <-summary["DIPrairies",]["mean"]
DIPrairies_p <- 1/(summary["DIPrairies",]["sd"])^2
DICentral_m <-summary["DICentral",]["mean"]
DICentral_p <- 1/(summary["DICentral",]["sd"])^2
DIincome_m <-summary["DIincome",]["mean"]
DIincome_p <- 1/(summary["DIincome",]["sd"])^2
DIalcohol_m <-summary["DIalcohol",]["mean"]
DIalcohol_p <- 1/(summary["DIalcohol",]["sd"])^2
DIrace_m <-summary["DIrace",]["mean"]
DIrace_p <- 1/(summary["DIrace",]["sd"])^2
DIpcp_m <-summary["DIpcp",]["mean"]
DIpcp_p <- 1/(summary["DIpcp",]["sd"])^2
DIsmoke_m <-summary["DIsmoke",]["mean"]
DIsmoke_p <- 1/(summary["DIsmoke",]["sd"])^2
DIdrug_m <-summary["DIdrug",]["mean"]
DIdrug_p <- 1/(summary["DIdrug",]["sd"])^2
DIrural_m <-summary["DIrural",]["mean"]
DIrural_p <- 1/(summary["DIrural",]["sd"])^2
DIedu_m1 <-summary["DIedu[1]",]["mean"]
DIedu_p1 <- 1/(summary["DIedu[1]",]["sd"])^2
DIedu_m2 <-summary["DIedu[2]",]["mean"]
DIedu_p2 <- 1/(summary["DIedu[2]",]["sd"])^2
DIedu_m3 <-summary["DIedu[3]",]["mean"]
DIedu_p3 <- 1/(summary["DIedu[3]",]["sd"])^2
DIemploy_m1 <-summary["DIemploy[1]",]["mean"]
DIemploy_p1 <- 1/(summary["DIemploy[1]",]["sd"])^2
DIemploy_m2 <-summary["DIemploy[2]",]["mean"]
DIemploy_p2 <- 1/(summary["DIemploy[2]",]["sd"])^2
DIemploy_m3 <-summary["DIemploy[3]",]["mean"]
DIemploy_p3 <- 1/(summary["DIemploy[3]",]["sd"])^2
DIimmig_m1 <-summary["DIimmig[1]",]["mean"]
DIimmig_p1 <- 1/(summary["DIimmig[1]",]["sd"])^2
DIimmig_m2 <-summary["DIimmig[2]",]["mean"]
DIimmig_p2 <- 1/(summary["DIimmig[2]",]["sd"])^2
DIimmig_m3 <-summary["DIimmig[3]",]["mean"]
DIimmig_p3 <- 1/(summary["DIimmig[3]",]["sd"])^2

cchsall_dep <- cchsall_bi_het %>% drop_na(depress_binary)
cchsall_consult <- cchsall_bi_het %>% drop_na(consult_mh)

NumCC <- nrow(cchsall_dep)
ageCC <- cchsall_dep$age_grp
incomeCC <- cchsall_dep$income
smokeCC <- cchsall_dep$smoke
eduCC <- cchsall_dep$education
drugCC <- cchsall_dep$illicit_drug
pcpCC <- cchsall_dep$have_pcp
employCC <- cchsall_dep$employ
provinceCC <- cchsall_dep$province
raceCC <- cchsall_dep$race
immigCC <- cchsall_dep$immigration + 1
alcoholCC <- cchsall_dep$alcohol
orientCC <- cchsall_dep$orient
existingCC <- cchsall_dep$existing_mh
provCentralCC <- cchsall_dep$province_Central
provPrairiesCC <- cchsall_dep$province_Prairies
provWestCC <- cchsall_dep$province_WestCoast
provNorthCC <- cchsall_dep$province_Northern
DeprCC <- cchsall_dep$depress_binary
consultCC <- cchsall_dep$consult_mh
ruralCC <- cchsall_dep$rural
wtsCC <- cchsall_dep$WTS_M_rescaled
MSMrep <- cchsall_dep$msm


model_codeCCHS_depress<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  DeprCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-DeprCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(DeprCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","DeprCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_dep$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_depress),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

if(T){
  SArray= jagsCCHS_wts$BUGSoutput$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){ 
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    xnul=locator(1)    
    acf( SArray[,1,nn], main=nn)  #note: this is only for 1st chain
    xnul=locator(1)
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
    xnul=locator(1)
  }
} 

NumCC <- nrow(cchsall_consult)
ageCC <- cchsall_consult$age_grp
incomeCC <- cchsall_consult$income
smokeCC <- cchsall_consult$smoke
eduCC <- cchsall_consult$education
drugCC <- cchsall_consult$illicit_drug
pcpCC <- cchsall_consult$have_pcp
employCC <- cchsall_consult$employ
provinceCC <- cchsall_consult$province
raceCC <- cchsall_consult$race
immigCC <- cchsall_consult$immigration + 1
alcoholCC <- cchsall_consult$alcohol
orientCC <- cchsall_consult$orient
existingCC <- cchsall_consult$existing_mh
provCentralCC <- cchsall_consult$province_Central
provPrairiesCC <- cchsall_consult$province_Prairies
provWestCC <- cchsall_consult$province_WestCoast
provNorthCC <- cchsall_consult$province_Northern
consultCC <- cchsall_consult$consult_mh
ruralCC <- cchsall_consult$rural
wtsCC <- cchsall_consult$WTS_M_rescaled
MSMrep <- cchsall_consult$msm

model_codeCCHS_consult<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  consultCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-consultCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(consultCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","consultCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_consult$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_consult),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

####Next, rerun everything for other MSM only, that has the lowest number of disclosures

SN_other <- SN %>% dplyr::filter(orient == 1)
cchsall_het <- cchsall %>% dplyr::filter(orient == 1)

NumSN <- nrow(SN_other)
ageSN <- SN_other$age_grp
incomeSN <- SN_other$income
smokeSN <- SN_other$smoke
eduSN <- SN_other$education
DiscSN <- SN_other$disclosure
drugSN <- SN_other$illicit_drug
pcpSN <- SN_other$have_pcp
employSN <- SN_other$employ
provinceSN <- SN_other$province
raceSN <- SN_other$race
immigSN <- SN_other$immigration + 1
alcoholSN <- SN_other$alcohol
orientSN <- SN_other$orient
existingSN <- SN_other$existing_mh
provCentralSN <- SN_other$province_Central
provPrairiesSN <- SN_other$province_Prairies
provWestSN <- SN_other$province_WestCoast
provNorthSN <- SN_other$province_Northern
DeprSN <- SN_other$depress_binary
ruralSN <- SN_other$rural

catimmig<-length(unique(immigSN))
catemploy<-length(unique(employSN))
catedu<-length(unique(eduSN))

model_codeSN<-"
model{
  for (iSN in 1:NumSN) {
    DiscSN[iSN] ~ dbern(pdiscSN[iSN])
    logit(pdiscSN[iSN]) <- DI0 + DIage*ageSN[iSN] + DIexisting*existingSN[iSN] + DIWest*provWestSN[iSN] + DINorth*provNorthSN[iSN] + DIPrairies*provPrairiesSN[iSN] + DICentral*provCentralSN[iSN] + DIincome*incomeSN[iSN] +
    DIalcohol*alcoholSN[iSN] + DIrace*raceSN[iSN] + DIpcp*pcpSN[iSN] + DIsmoke*smokeSN[iSN] + DIdrug*drugSN[iSN] + DIrural*ruralSN[iSN] + DIedu[eduSN[iSN]] + DIimmig[immigSN[iSN]] + DIemploy[employSN[iSN]]
  }
  for(ic in 1:catimmig){
    DIimmig.temp[ic] ~ dnorm(0,0.001)
   }
   DIimmig <- DIimmig.temp - mean(DIimmig.temp)
   for(ib in 1:catemploy){
     DIemploy.temp[ib] ~ dnorm(0,0.001)
   }
   DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  for(id in 1:catedu){
     DIedu.temp[id] ~ dnorm(0,0.001)
   }
   DIedu <- DIedu.temp - mean(DIedu.temp)
#Priors
DI0 ~ dnorm(0,.001)
DIage~dnorm(0, 0.001)
DIexisting~dnorm(0, 0.001)
DIWest~dnorm(0, 0.001)
DINorth~dnorm(0, 0.001)
DIPrairies~dnorm(0, 0.001)
DICentral~dnorm(0, 0.001)
DIincome~dnorm(0, 0.001)
DIalcohol~dnorm(0, 0.001)
DIrace~dnorm(0, 0.001)
DIpcp~dnorm(0, 0.001)
DIsmoke~dnorm(0, 0.001)
DIdrug~dnorm(0, 0.001)
DIrural~dnorm(0, 0.001)
}
"

##The best initial value will be the ones from the regression model itself
DI.mod<-glm(DiscSN ~ ageSN + existingSN + provWestSN + provNorthSN + provPrairiesSN + provCentralSN + incomeSN +
              alcoholSN + raceSN + pcpSN + smokeSN + raceSN + drugSN + ruralSN + as.factor(eduSN) + as.factor(immigSN) + as.factor(employSN),family=binomial)

init.vals <- function(){
  list(
    DI0=rnorm(1,mean=summary(DI.mod)$coefficients[1,1],sd=summary(DI.mod)$coefficients[1,2]),
    DIage=rnorm(1,mean=summary(DI.mod)$coefficients[2,1],sd=summary(DI.mod)$coefficients[2,2]),
    DIexisting=rnorm(1,mean=summary(DI.mod)$coefficients[3,1],sd=summary(DI.mod)$coefficients[3,2]),
    DIWest=rnorm(1,mean=summary(DI.mod)$coefficients[4,1],sd=summary(DI.mod)$coefficients[4,2]),
    DINorth=rnorm(1,mean=summary(DI.mod)$coefficients[5,1],sd=summary(DI.mod)$coefficients[5,2]),
    DIPrairies=rnorm(1,mean=summary(DI.mod)$coefficients[6,1],sd=summary(DI.mod)$coefficients[6,2]),
    DICentral=rnorm(1,mean=summary(DI.mod)$coefficients[7,1],sd=summary(DI.mod)$coefficients[7,2]),
    DIincome=rnorm(1,mean=summary(DI.mod)$coefficients[8,1],sd=summary(DI.mod)$coefficients[8,2]),
    DIalcohol=rnorm(1,mean=summary(DI.mod)$coefficients[9,1],sd=summary(DI.mod)$coefficients[9,2]),
    DIrace=rnorm(1,mean=summary(DI.mod)$coefficients[10,1],sd=summary(DI.mod)$coefficients[10,2]),
    DIpcp=rnorm(1,mean=summary(DI.mod)$coefficients[11,1],sd=summary(DI.mod)$coefficients[11,2]),
    DIsmoke=rnorm(1,mean=summary(DI.mod)$coefficients[12,1],sd=summary(DI.mod)$coefficients[12,2]),
    DIdrug=rnorm(1,mean=summary(DI.mod)$coefficients[13,1],sd=summary(DI.mod)$coefficients[13,2]),
    DIrural=rnorm(1,mean=summary(DI.mod)$coefficients[14,1],sd=summary(DI.mod)$coefficients[14,2]))
}
init.vals()

jagsSN.dat=list("NumSN","DiscSN","ageSN","existingSN","drugSN","immigSN","employSN","ruralSN","eduSN",
                "pcpSN","smokeSN","raceSN","alcoholSN","incomeSN","provCentralSN","provPrairiesSN","provWestSN","provNorthSN",
                "catedu","catimmig","catemploy")

paramsSN=c("DI0","DIage","DIexisting","DIWest","DINorth","DIPrairies","DICentral","DIincome","DIalcohol","DIrace",
           "DIpcp","DIsmoke","DIdrug","DIrural","DIedu","DIemploy","DIimmig")

start.time <- Sys.time()
jagsSN=jags(data=jagsSN.dat, inits=init.vals,parameters.to.save=paramsSN, model.file=textConnection(model_codeSN),
            n.chains=3, n.iter=10000, n.burnin=100, n.thin=1,jags.seed=123)
print(jagsSN)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

summary=jagsSN$BUGSoutput$summary
DI0_m<- summary["DI0",]["mean"]
DI0_p<- 1/(summary["DI0",]["sd"])^2
DIage_m<- summary["DIage",]["mean"]
DIage_p<- 1/(summary["DIage",]["sd"])^2
DIexisting_m <- summary["DIexisting",]["mean"]
DIexisting_p <- 1/(summary["DIexisting",]["sd"])^2
DIWest_m <-summary["DIWest",]["mean"]
DIWest_p <- 1/(summary["DIWest",]["sd"])^2
DINorth_m <-summary["DINorth",]["mean"]
DINorth_p <- 1/(summary["DINorth",]["sd"])^2
DIPrairies_m <-summary["DIPrairies",]["mean"]
DIPrairies_p <- 1/(summary["DIPrairies",]["sd"])^2
DICentral_m <-summary["DICentral",]["mean"]
DICentral_p <- 1/(summary["DICentral",]["sd"])^2
DIincome_m <-summary["DIincome",]["mean"]
DIincome_p <- 1/(summary["DIincome",]["sd"])^2
DIalcohol_m <-summary["DIalcohol",]["mean"]
DIalcohol_p <- 1/(summary["DIalcohol",]["sd"])^2
DIrace_m <-summary["DIrace",]["mean"]
DIrace_p <- 1/(summary["DIrace",]["sd"])^2
DIpcp_m <-summary["DIpcp",]["mean"]
DIpcp_p <- 1/(summary["DIpcp",]["sd"])^2
DIsmoke_m <-summary["DIsmoke",]["mean"]
DIsmoke_p <- 1/(summary["DIsmoke",]["sd"])^2
DIdrug_m <-summary["DIdrug",]["mean"]
DIdrug_p <- 1/(summary["DIdrug",]["sd"])^2
DIrural_m <-summary["DIrural",]["mean"]
DIrural_p <- 1/(summary["DIrural",]["sd"])^2
DIedu_m1 <-summary["DIedu[1]",]["mean"]
DIedu_p1 <- 1/(summary["DIedu[1]",]["sd"])^2
DIedu_m2 <-summary["DIedu[2]",]["mean"]
DIedu_p2 <- 1/(summary["DIedu[2]",]["sd"])^2
DIedu_m3 <-summary["DIedu[3]",]["mean"]
DIedu_p3 <- 1/(summary["DIedu[3]",]["sd"])^2
DIemploy_m1 <-summary["DIemploy[1]",]["mean"]
DIemploy_p1 <- 1/(summary["DIemploy[1]",]["sd"])^2
DIemploy_m2 <-summary["DIemploy[2]",]["mean"]
DIemploy_p2 <- 1/(summary["DIemploy[2]",]["sd"])^2
DIemploy_m3 <-summary["DIemploy[3]",]["mean"]
DIemploy_p3 <- 1/(summary["DIemploy[3]",]["sd"])^2
DIimmig_m1 <-summary["DIimmig[1]",]["mean"]
DIimmig_p1 <- 1/(summary["DIimmig[1]",]["sd"])^2
DIimmig_m2 <-summary["DIimmig[2]",]["mean"]
DIimmig_p2 <- 1/(summary["DIimmig[2]",]["sd"])^2
DIimmig_m3 <-summary["DIimmig[3]",]["mean"]
DIimmig_p3 <- 1/(summary["DIimmig[3]",]["sd"])^2

#Dep and consult other MSMs in CCHS
cchsall_dep <- cchsall_het %>% drop_na(depress_binary)
cchsall_consult <- cchsall_het %>% drop_na(consult_mh)

NumCC <- nrow(cchsall_dep)
ageCC <- cchsall_dep$age_grp
incomeCC <- cchsall_dep$income
smokeCC <- cchsall_dep$smoke
eduCC <- cchsall_dep$education
drugCC <- cchsall_dep$illicit_drug
pcpCC <- cchsall_dep$have_pcp
employCC <- cchsall_dep$employ
provinceCC <- cchsall_dep$province
raceCC <- cchsall_dep$race
immigCC <- cchsall_dep$immigration + 1
alcoholCC <- cchsall_dep$alcohol
orientCC <- cchsall_dep$orient
existingCC <- cchsall_dep$existing_mh
provCentralCC <- cchsall_dep$province_Central
provPrairiesCC <- cchsall_dep$province_Prairies
provWestCC <- cchsall_dep$province_WestCoast
provNorthCC <- cchsall_dep$province_Northern
DeprCC <- cchsall_dep$depress_binary
consultCC <- cchsall_dep$consult_mh
ruralCC <- cchsall_dep$rural
wtsCC <- cchsall_dep$WTS_M_rescaled
MSMrep <- cchsall_dep$msm

model_codeCCHS_depress<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  DeprCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-DeprCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(DeprCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","DeprCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("SumMSMCC","sumtotdep","ratedepMSM")

r.mod<-glm(cchsall_dep$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_depress),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100
# mean 
# 3.709096 
# > jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
# 2.5% 
# 3.030303 
# > jagsCCsum_wts["ratedepMSM",]["97.5%"]*100
# 97.5% 
# 7.844423 

if(T){
  SArray= jagsCCHS_wts$BUGSoutput$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){ 
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    xnul=locator(1)    
    acf( SArray[,1,nn], main=nn)  #note: this is only for 1st chain
    xnul=locator(1)
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
    xnul=locator(1)
  }
} 

NumCC <- nrow(cchsall_consult)
ageCC <- cchsall_consult$age_grp
incomeCC <- cchsall_consult$income
smokeCC <- cchsall_consult$smoke
eduCC <- cchsall_consult$education
drugCC <- cchsall_consult$illicit_drug
pcpCC <- cchsall_consult$have_pcp
employCC <- cchsall_consult$employ
provinceCC <- cchsall_consult$province
raceCC <- cchsall_consult$race
immigCC <- cchsall_consult$immigration + 1
alcoholCC <- cchsall_consult$alcohol
orientCC <- cchsall_consult$orient
existingCC <- cchsall_consult$existing_mh
provCentralCC <- cchsall_consult$province_Central
provPrairiesCC <- cchsall_consult$province_Prairies
provWestCC <- cchsall_consult$province_WestCoast
provNorthCC <- cchsall_consult$province_Northern
consultCC <- cchsall_consult$consult_mh
ruralCC <- cchsall_consult$rural
wtsCC <- cchsall_consult$WTS_M_rescaled
MSMrep <- cchsall_consult$msm

model_codeCCHS_consult<-"
model{
  for (iCC in 1:NumCC) {
    DiscCC[iCC] ~ dbern(pDisc[iCC])
    logit(pDisc[iCC]) <- r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]
    MSMrep[iCC] ~ dpois(pMSM[iCC] * wtsCC[iCC])
    pMSM[iCC] <- exp(-(DI0 + DIage*ageCC[iCC] + DIexisting*existingCC[iCC] + DIWest*provWestCC[iCC] + DINorth*provNorthCC[iCC] + DIPrairies*provPrairiesCC[iCC] + DICentral*provCentralCC[iCC] + DIincome*incomeCC[iCC] +
    DIalcohol*alcoholCC[iCC] + DIrace*raceCC[iCC] + DIpcp*pcpCC[iCC] + DIsmoke*smokeCC[iCC] + DIdrug*drugCC[iCC] + DIrural*ruralCC[iCC] + DIedu[eduCC[iCC]] + DIimmig[immigCC[iCC]] + DIemploy[employCC[iCC]]) + 
    (r0 + r.age*ageCC[iCC] + r.existing*existingCC[iCC] + r.West*provWestCC[iCC] + r.North*provNorthCC[iCC] + r.Prairies*provPrairiesCC[iCC] + r.Central*provCentralCC[iCC] + r.income*incomeCC[iCC] +
    r.alcohol*alcoholCC[iCC] + r.race*raceCC[iCC] + r.pcp*pcpCC[iCC] + r.smoke*smokeCC[iCC] + r.drug*drugCC[iCC] + r.rural*ruralCC[iCC] + r.edu[eduCC[iCC]] + r.immig[immigCC[iCC]] + r.employ[employCC[iCC]]))
    MSMCC[iCC] <- DiscCC[iCC] + (1-DiscCC[iCC])*MSMrep[iCC]
	  consultCC[iCC] ~ dpois(pDepCC[iCC] * wtsCC[iCC])
	logit(pDepCC[iCC]) <- Dp0 
	totdepMSM[iCC]<-consultCC[iCC]*MSMCC[iCC]
  }
  DIimmig.temp[1] ~ dnorm(DIimmig_m1,DIimmig_p1)
  DIimmig.temp[2] ~ dnorm(DIimmig_m2,DIimmig_p2)
  DIimmig.temp[3] ~ dnorm(DIimmig_m3,DIimmig_p3)
  for(ig in 1:catimmig){
    r.immig.temp[ig] ~ dnorm(0,0.001)
  }
  DIimmig <- DIimmig.temp - mean(DIimmig.temp)
  r.immig <- r.immig.temp - mean(r.immig.temp)
  DIemploy.temp[1] ~ dnorm(DIemploy_m1,DIemploy_p1)
  DIemploy.temp[2] ~ dnorm(DIemploy_m2,DIemploy_p2)
  DIemploy.temp[3] ~ dnorm(DIemploy_m3,DIemploy_p3)
  for(ie in 1:catemploy){
    r.employ.temp[ie] ~ dnorm(0,0.001)
  }
  DIemploy <- DIemploy.temp - mean(DIemploy.temp)
  r.employ <- r.employ.temp - mean(r.employ.temp)
  DIedu.temp[1] ~ dnorm(DIedu_m1,DIedu_p1)
  DIedu.temp[2] ~ dnorm(DIedu_m2,DIedu_p2)
  DIedu.temp[3] ~ dnorm(DIedu_m3,DIedu_p3)
  for(ie in 1:catedu){
    r.edu.temp[ie] ~ dnorm(0,0.001)
  }
  DIedu <- DIedu.temp - mean(DIedu.temp)
  r.edu <- r.edu.temp - mean(r.edu.temp)
SumMSMCC <- sum(MSMCC[])
sumtotdep<-sum(consultCC[])
ratedepMSM<-sum(totdepMSM[])/sum(MSMCC[])
#Priors
DI0 ~ dnorm(DI0_m,DI0_p)
DIage ~ dnorm(DIage_m,DIage_p)
DIexisting ~ dnorm(DIexisting_m,DIexisting_p)
DIWest ~ dnorm(DIWest_m,DIWest_p)
DINorth ~ dnorm(DINorth_m,DINorth_p)
DIPrairies ~ dnorm(DIPrairies_m,DIPrairies_p)
DICentral ~ dnorm(DICentral_m,DICentral_p)
DIincome ~ dnorm(DIincome_m,DIincome_p)
DIalcohol ~ dnorm(DIalcohol_m,DIalcohol_p)
DIrace ~ dnorm(DIrace_m,DIrace_p)
DIpcp ~ dnorm(DIpcp_m,DIpcp_p)
DIsmoke ~ dnorm(DIsmoke_m,DIsmoke_p)
DIdrug ~ dnorm(DIdrug_m,DIdrug_p)
DIrural ~ dnorm(DIrural_m,DIrural_p)
r0 ~ dnorm(0, 0.001)
r.age ~ dnorm(0,.001)
r.existing ~ dnorm(0,.001)
r.West ~ dnorm(0,.001)
r.North ~ dnorm(0,.001)
r.Prairies ~ dnorm(0,.001)
r.Central ~ dnorm(0,.001)
r.income ~ dnorm(0,.001)
r.alcohol ~ dnorm(0,.001)
r.race ~ dnorm(0,.001)
r.pcp ~ dnorm(0,.001)
r.smoke ~ dnorm(0,.001)
r.drug ~ dnorm(0,.001)
r.rural ~ dnorm(0,.001)
Dp0 ~ dnorm(0,.001)
DpMSM ~ dnorm(0,.001)
}
"

jagsCC.dat=list("NumCC","ageCC","existingCC","drugCC","immigCC","employCC","ruralCC","eduCC","consultCC","MSMrep",
                "pcpCC","smokeCC","raceCC","alcoholCC","incomeCC","provCentralCC","provPrairiesCC","provWestCC","provNorthCC",
                "catedu","catimmig","catemploy","wtsCC","DI0_m","DI0_p","DIage_m","DIage_p","DIexisting_m","DIexisting_p",
                "DIWest_m","DIWest_p","DINorth_m","DINorth_p","DIPrairies_m","DIPrairies_p",
                "DICentral_m","DICentral_p","DIincome_m","DIincome_p","DIalcohol_m","DIalcohol_p",
                "DIrace_m","DIrace_p","DIpcp_m","DIpcp_p","DIsmoke_m","DIsmoke_p",
                "DIdrug_m","DIdrug_p","DIrural_m","DIrural_p","DIedu_m1","DIedu_p1",
                "DIedu_m2","DIedu_p2","DIedu_m3","DIedu_p3","DIemploy_m1","DIemploy_p1",
                "DIemploy_m2","DIemploy_p2","DIemploy_m3","DIemploy_p3","DIimmig_m1","DIimmig_p1",
                "DIimmig_m2","DIimmig_p2","DIimmig_m3","DIimmig_p3")

paramsCCHS=c("ratedepMSM")

##The best initial value will be the ones from the regression model itself

r.mod<-glm(cchsall_consult$disclosure ~ ageCC + existingCC + provWestCC + provNorthCC + provPrairiesCC + provCentralCC + incomeCC +
             alcoholCC+ raceCC + pcpCC + smokeCC + raceCC + drugCC + ruralCC + as.factor(eduCC) + as.factor(immigCC) + as.factor(employCC),family=binomial)

init.vals <- function(){
  list(
    r0=rnorm(1,mean=summary(r.mod)$coefficients[1,1],sd=summary(r.mod)$coefficients[1,2]),
    r.age=rnorm(1,mean=summary(r.mod)$coefficients[2,1],sd=summary(r.mod)$coefficients[2,2]),
    r.existing=rnorm(1,mean=summary(r.mod)$coefficients[3,1],sd=summary(r.mod)$coefficients[3,2]),
    r.West=rnorm(1,mean=summary(r.mod)$coefficients[4,1],sd=summary(r.mod)$coefficients[4,2]),
    r.North=rnorm(1,mean=summary(r.mod)$coefficients[5,1],sd=summary(r.mod)$coefficients[5,2]),
    r.Prairies=rnorm(1,mean=summary(r.mod)$coefficients[6,1],sd=summary(r.mod)$coefficients[6,2]),
    r.Central=rnorm(1,mean=summary(r.mod)$coefficients[7,1],sd=summary(r.mod)$coefficients[7,2]),
    r.income=rnorm(1,mean=summary(r.mod)$coefficients[8,1],sd=summary(r.mod)$coefficients[8,2]),
    r.alcohol=rnorm(1,mean=summary(r.mod)$coefficients[9,1],sd=summary(r.mod)$coefficients[9,2]),
    r.race=rnorm(1,mean=summary(r.mod)$coefficients[10,1],sd=summary(r.mod)$coefficients[10,2]),
    r.pcp=rnorm(1,mean=summary(r.mod)$coefficients[11,1],sd=summary(r.mod)$coefficients[11,2]),
    r.smoke=rnorm(1,mean=summary(r.mod)$coefficients[12,1],sd=summary(r.mod)$coefficients[12,2]),
    r.drug=rnorm(1,mean=summary(r.mod)$coefficients[13,1],sd=summary(r.mod)$coefficients[13,2]),
    r.rural=rnorm(1,mean=summary(r.mod)$coefficients[14,1],sd=summary(r.mod)$coefficients[14,2]))
}
init.vals()

start.time <- Sys.time()
jagsCCHS_wts=jags(data=jagsCC.dat,inits=init.vals, parameters.to.save=paramsCCHS, model.file=textConnection(model_codeCCHS_consult),
                  n.chains=3, n.iter=1000, n.burnin=200, n.thin=1,jags.seed=123)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

jagsCCsum_wts<-jagsCCHS_wts$BUGSoutput$summary
print(jagsCCsum_wts)
jagsCCsum_wts["ratedepMSM",]["mean"]*100
jagsCCsum_wts["ratedepMSM",]["2.5%"]*100
jagsCCsum_wts["ratedepMSM",]["97.5%"]*100

jagsCCHS_wts.mcmc<-as.mcmc(jagsCCHS_wts) 
summary(jagsCCHS_wts.mcmc)

library(coda)
geweke.diag(jagsCCHS_wts.mcmc)
geweke.plot(jagsCCHS_wts.mcmc)

gelman.diag(jagsCCHS_wts.mcmc)
gelman.plot(jagsCCHS_wts.mcmc)



