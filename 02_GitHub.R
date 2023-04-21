library(haven)
library(dplyr)
library(naniar)
library(survey) 
library(tableone)
library(tidyverse)
library(boot)
library(rsample)

load("cchsSNall.Rdata")

cchsSNall <- cchsSNall %>% dplyr::mutate(employ = case_when(employ == 1 ~ 1,
                                                     employ %in% c(2,3,4) ~ 2,
                                                    employ %in% c(5,6) ~ 3))

names<-c("education","orient","immigration","existing_mh","consult_mh","alcohol","marijuana","smoke","insurance","illicit_drug","rural",
         "have_pcp","employ","province","race","depress_binary","HIV_test12m","STI_test12m","lastcondom_use","survey","disclosure")
cchsSNall[,names]<-lapply(cchsSNall[,names],factor)

svyds.ALP.CC.all  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = cchsSNall) ### 

svyCreateTableOne(vars=c("consult_mh"), factorVars = c("consult_mh"),strata="survey",data=svyds.ALP.CC.all)


cchsSNall2 <- cchsSNall %>% 
  dplyr::select(-c(HIV_test12m,STI_test12m,lastcondom_use,insurance,marijuana))

alldat <- cchsSNall2 %>% drop_na(age_grp,income,education,orient,msm,illicit_drug,have_pcp,smoke,employ,province,rural,race,existing_mh,immigration,alcohol)

bsw2015 <- read_sas("bsw2015.sas7bdat")
bsw2017 <- read_sas("bsw2017.sas7bdat")

bsw2015 <- bsw2015 %>% dplyr::select(-c(fwgt))
bsw2017 <- bsw2017 %>% dplyr::select(-c(fwgt))

#create bootstraps
B = 1000

cchs2015 <- alldat %>% dplyr::filter(year == 2015)
cchs2017 <- alldat %>% dplyr::filter(year == 2017)

sndat_all_bsw <- alldat %>% dplyr::filter(survey == 1)
  
for (i in 1:B) {
  sndat_all_bsw[[as.symbol(paste0('BSW', i))]] <- 1
}

cchs2015_bsw <- merge(x = cchs2015, y = bsw2015, by = "ADM_RNO", all.x = TRUE)
cchs2017_bsw <- merge(x = cchs2017, y = bsw2017, by = "ADM_RNO", all.x = TRUE)
cchsbsw <- rbind(cchs2015_bsw,cchs2017_bsw)
nyear <- 2

cchsbsw <- cchsbsw %>% dplyr::mutate_at(.,vars(matches("BSW")), funs(. / nyear))

alldat_gay_cchs <- cchsbsw %>% dplyr::filter(orient==2)
alldat_gay_sn <- sndat_all_bsw %>% dplyr::filter(orient==2)

alldat_gay <- alldat %>% dplyr::filter(orient == 2)

svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat_gay) ### 
Formula_fit = as.formula("survey ~ age_grp + age_grp^2 + income + smoke + education + illicit_drug + employ + smoke + alcohol + province + age_grp*income + age_grp*employ + rural + race + existing_mh + immigration") 
lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)    
beta.w = summary(lgtreg.w)$coeff[,1]   
alldat_gay$ps.w = lgtreg.w$fitted.values 

alldat_gay <- alldat_gay %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M_rescaled,ALP_weight)) 

svyds.ALP_gay  = svydesign(ids=~1, weight = ~ ALP_weight, data = alldat_gay) ### 

cchs_gay <- alldat_gay %>% filter(survey == 0)
SN_gay <- alldat_gay %>% filter(survey == 1)

tabw<-svyCreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=svyds.ALP_gay,strata="survey",test=FALSE,smd = TRUE)
tab2<-print(tabw,smd=TRUE)

svyds.ALP.SN.w  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN_gay) ### 
svyds.ALP.CC.w  = svydesign(ids=~1, weight = ~ ALP_weight, data = cchs_gay) ### 

point.est.gay.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w)
point.est.gay.cchs<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.CC.w)

#unadjusted SN
CreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=SN_gay)

nfeature = 2
prev.list<-vector("list",nfeature)

bt_samples_cchs <- bootstraps(alldat_gay_cchs, times = B)
bt_samples_sn <- bootstraps(alldat_gay_sn, times = B)

start.time <- Sys.time()
for (i in 1:B) {
  x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
  x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
  x.all<-rbind(x.cchs,x.sn)
  variable = paste0("BSW",i)
  svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = x.all) ### 
  lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)   
  beta.w = summary(lgtreg.w)$coeff[,1]   
  x.all$ps.w = lgtreg.w$fitted.values 
  SN_gay <- x.all %>% filter(survey == 1)
  SN_gay <- SN_gay %>%
    mutate(ALP_weight = (1-ps.w) /ps.w)
  svyds.ALP.SN.w  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN_gay)
  svyprev<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w)
  prev.list[[1]]<-c(prev.list[[1]],svyprev$CatTable$Overall$consult_mh$percent[2])
  prev.list[[2]]<-c(prev.list[[2]],svyprev$CatTable$Overall$depress_binary$percent[2])
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

prev.ci<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))
                  
prev.ci$variable <- c("consult mh","depression")
prev.ci$estimate <-c(point.est.gay.adj$CatTable$Overall$consult_mh$percent[2],point.est.gay.adj$CatTable$Overall$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci$lower[k]<-quantile(prev.list[[k]],.025)
  prev.ci$upper[k]<-quantile(prev.list[[k]],.975)
}
prev.ci

CreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=SN_gay)

depSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(consult_mh) %>% select(consult_mh)))
consultSN <- consultSN - 1
depSN <- depSN - 1

fc <- function(d, i){
  d2 <- d[i]
  return(mean(d2))
}

b1<-boot(data=consultSN,fc,R=B)
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN,fc,R=B)
boot.ci(b2, conf = 0.95,type = "perc")

#Try to generalize to all MSM from CCHS

svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat) 

Formula_fit = as.formula("survey ~ age_grp + age_grp^2 + income + income^2 + age_grp*education + age_grp*employ + age_grp*income + income*province + smoke + education + illicit_drug + employ + smoke + alcohol + province + rural + race + existing_mh + immigration") 

lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)    
beta.w = summary(lgtreg.w)$coeff[,1]   
alldat$ps.w = lgtreg.w$fitted.values 

alldat <- alldat %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M_rescaled,ALP_weight)) 

svyds.ALP_all  = svydesign(ids=~1, weight = ~ ALP_weight, data = alldat) ### 

tabw<-svyCreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=svyds.ALP_all,strata="survey",test=FALSE,smd = TRUE)
tab3<-print(tabw,smd=TRUE)

###Try all MSM only
alldat_cchs <- cchsbsw %>% dplyr::filter(msm==1)
alldat_sn <- sndat_all_bsw 
alldat_msm <- alldat %>% dplyr::filter(msm==1)
svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat_msm) 

Formula_fit = as.formula("survey ~ age_grp + age_grp^2 + income + income^2 + age_grp*education + age_grp*employ + age_grp*income + income*province + smoke + education + illicit_drug + employ + smoke + alcohol + province + rural + race + existing_mh + immigration") 

lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)    
beta.w = summary(lgtreg.w)$coeff[,1]   
alldat_msm$ps.w = lgtreg.w$fitted.values 

alldat_msm <- alldat_msm %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M_rescaled,ALP_weight)) 

##Check if the weights are correctly applied
svyds.ALP_all.msm  = svydesign(ids=~1, weight = ~ ALP_weight, data = alldat_msm) ### 

tabw<-svyCreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","age_grp","education","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=svyds.ALP_all.msm,strata="survey",test=FALSE,smd = TRUE)
tab4<-print(tabw,smd=TRUE)

#preselect the three categories, bi, other, and all, which are recalibrated on all MSMs from CCHS
allsn_bi <- alldat_msm %>% dplyr::filter(orient == 3 & survey == 1)
allsn_other <- alldat_msm %>% dplyr::filter(orient == 1 & survey == 1)
allsn <- alldat_msm %>% dplyr::filter(survey == 1)

#Check consult_mh, SN weighted, CCHS weighted, SN unweighted
svyds.ALP.SN.w.bi  = svydesign(ids=~1, weight = ~ ALP_weight, data = allsn_bi) ### 
svyds.ALP.SN.w.other  = svydesign(ids=~1, weight = ~ ALP_weight, data = allsn_other) ### 
svyds.ALP.SN.w.msm  = svydesign(ids=~1, weight = ~ ALP_weight, data = allsn) ### 

point.est.bi.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.bi)
point.est.other.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.other)
point.est.allmsm<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.msm)

#unadjusted SN
CreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=allsn_bi)
CreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=allsn_other)
CreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=allsn)

nfeature = 2
prev.list.bi<-vector("list",nfeature)
prev.list.other<-vector("list",nfeature)
prev.list.all.msm<-vector("list",nfeature)

bt_samples_cchs <- bootstraps(alldat_cchs, times = B)
bt_samples_sn <- bootstraps(alldat_sn, times = B)

start.time <- Sys.time()
for (i in 1:B) {
  x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
  x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
  x.all<-rbind(x.cchs,x.sn)
  variable = paste0("BSW",i)
  svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = x.all) ### 
  lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)   
  beta.w = summary(lgtreg.w)$coeff[,1]   
  x.all$ps.w = lgtreg.w$fitted.values 
  SN.all <- x.all %>% dplyr::filter(survey == 1)
  SN.all <- SN.all %>%
    mutate(ALP_weight = (1-ps.w) /ps.w)
  SN.bi <- SN.all %>% dplyr::filter(orient == 3)
  SN.other <- SN.all %>% dplyr::filter(orient == 1)
  svyds.ALP.SN.w.bi  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN.bi)
  svyds.ALP.SN.w.other  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN.other)
  svyds.ALP.SN.w.msm  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN.all) ### 
  svyprev.bi<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.bi)
  svyprev.other<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.other)
  svyprev.msm<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.SN.w.msm)
  prev.list.bi[[1]]<-c(prev.list.bi[[1]],svyprev.bi$CatTable$Overall$consult_mh$percent[2])
  prev.list.bi[[2]]<-c(prev.list.bi[[2]],svyprev.bi$CatTable$Overall$depress_binary$percent[2])
  prev.list.other[[1]]<-c(prev.list.other[[1]],svyprev.other$CatTable$Overall$consult_mh$percent[2])
  prev.list.other[[2]]<-c(prev.list.other[[2]],svyprev.other$CatTable$Overall$depress_binary$percent[2])
  prev.list.all.msm[[1]]<-c(prev.list.all.msm[[1]],svyprev.msm$CatTable$Overall$consult_mh$percent[2])
  prev.list.all.msm[[2]]<-c(prev.list.all.msm[[2]],svyprev.msm$CatTable$Overall$depress_binary$percent[2])
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

prev.ci.bi<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))
                  
prev.ci.bi$variable <- c("consult mh","depression")
prev.ci.bi$estimate <-c(point.est.bi.adj$CatTable$Overall$consult_mh$percent[2],point.est.bi.adj$CatTable$Overall$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.bi$lower[k]<-quantile(prev.list.bi[[k]],.025)
  prev.ci.bi$upper[k]<-quantile(prev.list.bi[[k]],.975)
}
prev.ci.bi

prev.ci.other<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.other$variable <- c("consult mh","depression")
prev.ci.other$estimate <-c(point.est.other.adj$CatTable$Overall$consult_mh$percent[2],point.est.other.adj$CatTable$Overall$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.other$lower[k]<-quantile(prev.list.other[[k]],.025)
  prev.ci.other$upper[k]<-quantile(prev.list.other[[k]],.975)
}
prev.ci.other

prev.ci.allmsm<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.allmsm$variable <- c("consult mh","depression")
prev.ci.allmsm$estimate <-c(point.est.allmsm$CatTable$Overall$consult_mh$percent[2],point.est.allmsm$CatTable$Overall$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.allmsm$lower[k]<-quantile(prev.list.all.msm[[k]],.025)
  prev.ci.allmsm$upper[k]<-quantile(prev.list.all.msm[[k]],.975)
}
prev.ci.allmsm

#Also bootstrap unadjusted estimates for bi and other MSM from SN

depSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(consult_mh) %>% select(consult_mh)))
consultSN.bi <- consultSN.bi - 1
depSN.bi <- depSN.bi - 1

fc <- function(d, i){
  d2 <- d[i]
  return(mean(d2))
}

b1<-boot(data=consultSN.bi,fc,R=B)
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.bi,fc,R=B)
boot.ci(b2, conf = 0.95,type = "perc")

depSN.other <- as.numeric(unlist(allsn_other %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.other <- as.numeric(unlist(allsn_other %>% drop_na(consult_mh) %>% select(consult_mh)))
consultSN.other <- consultSN.other - 1
depSN.other <- depSN.other - 1

b1<-boot(data=consultSN.other,fc,R=B)
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.other,fc,R=B)
boot.ci(b2, conf = 0.95,type = "perc")

depSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(consult_mh) %>% select(consult_mh)))
consultSN.allmsm <- consultSN.allmsm - 1
depSN.allmsm <- depSN.allmsm - 1

b1<-boot(data=consultSN.allmsm,fc,R=B)
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.allmsm,fc,R=B)
boot.ci(b2, conf = 0.95,type = "perc")

###CCHS bootstrap weights
svyds.ALP.CC.w  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat_cchs) ### 
point.est.cchs.all<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), data=svyds.ALP.CC.w)
point.est.cchs.orient<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), strata="orient",data=svyds.ALP.CC.w)

nfeature = 2 
bt_samples_cchs <- bootstraps(alldat_cchs, times = B)
prev.list.cchs.gay<-vector("list",nfeature)
prev.list.cchs.all<-vector("list",nfeature)
prev.list.cchs.bi<-vector("list",nfeature)
prev.list.cchs.other<-vector("list",nfeature)

start.time <- Sys.time()
for (i in 1:B) {
  all.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
  variable = paste0("BSW",i)
  svyds.ALP.cc.all  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = all.cchs)
  cchs.boot.orient<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"), strata="orient",data=svyds.ALP.cc.all)
  cchs.boot.all<-svyCreateTableOne(vars=c("consult_mh","depress_binary"), factorVars = c("consult_mh","depress_binary"),data=svyds.ALP.cc.all)
  prev.list.cchs.all[[1]]<-c(prev.list.cchs.all[[1]],cchs.boot.all$CatTable$Overall$consult_mh$percent[2])
  prev.list.cchs.all[[2]]<-c(prev.list.cchs.all[[2]],cchs.boot.all$CatTable$Overall$depress_binary$percent[2])
  prev.list.cchs.other[[1]]<-c(prev.list.cchs.other[[1]],cchs.boot.orient$CatTable$"1"$consult_mh$percent[2])
  prev.list.cchs.other[[2]]<-c(prev.list.cchs.other[[2]],cchs.boot.orient$CatTable$"1"$depress_binary$percent[2])
  prev.list.cchs.gay[[1]]<-c(prev.list.cchs.gay[[1]],cchs.boot.orient$CatTable$"2"$consult_mh$percent[2])
  prev.list.cchs.gay[[2]]<-c(prev.list.cchs.gay[[2]],cchs.boot.orient$CatTable$"2"$depress_binary$percent[2])
  prev.list.cchs.bi[[1]]<-c(prev.list.cchs.bi[[1]],cchs.boot.orient$CatTable$"3"$consult_mh$percent[2])
  prev.list.cchs.bi[[2]]<-c(prev.list.cchs.bi[[2]],cchs.boot.orient$CatTable$"3"$depress_binary$percent[2])
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

prev.ci.allmsm.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.allmsm.cchs$variable <- c("consult mh","depression")
prev.ci.allmsm.cchs$estimate <-c(point.est.cchs.all$CatTable$Overall$consult_mh$percent[2],point.est.cchs.all$CatTable$Overall$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.allmsm.cchs$lower[k]<-quantile(prev.list.cchs.all[[k]],.025)
  prev.ci.allmsm.cchs$upper[k]<-quantile(prev.list.cchs.all[[k]],.975)
}
prev.ci.allmsm.cchs

prev.ci.gay.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.gay.cchs$variable <- c("consult mh","depression")
prev.ci.gay.cchs$estimate <-c(point.est.cchs.orient$CatTable$"2"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"2"$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.gay.cchs$lower[k]<-quantile(prev.list.cchs.gay[[k]],.025)
  prev.ci.gay.cchs$upper[k]<-quantile(prev.list.cchs.gay[[k]],.975)
}
prev.ci.gay.cchs

prev.ci.bi.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.bi.cchs$variable <- c("consult mh","depression")
prev.ci.bi.cchs$estimate <-c(point.est.cchs.orient$CatTable$"3"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"3"$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.bi.cchs$lower[k]<-quantile(prev.list.cchs.bi[[k]],.025)
  prev.ci.bi.cchs$upper[k]<-quantile(prev.list.cchs.bi[[k]],.975)
}
prev.ci.bi.cchs

prev.ci.other.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.other.cchs$variable <- c("consult mh","depression")
prev.ci.other.cchs$estimate <-c(point.est.cchs.orient$CatTable$"1"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"1"$depress_binary$percent[2])

for (k in 1:nfeature) {
  prev.ci.other.cchs$lower[k]<-quantile(prev.list.cchs.other[[k]],.025)
  prev.ci.other.cchs$upper[k]<-quantile(prev.list.cchs.other[[k]],.975)
}
prev.ci.other.cchs

tab5<-CreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","age_grp","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=alldat_sn)
tab5<-print(tab5)
tab6<-CreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","age_grp","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=SN_gay)
tab6<-print(tab6)


