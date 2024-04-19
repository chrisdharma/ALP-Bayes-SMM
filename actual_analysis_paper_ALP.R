library(haven)
library(dplyr)
library(naniar)
library(survey) 
library(tableone)
library(tidyverse)
library(boot)
library(rsample)

load("cchsSNall.Rdata")

###Need to match the bootstraps
setwd("/Users/christofferdharma/Documents/CCHS/AllDat")

#Create data
svyds.ALP.CC.all  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = cchsSNall) ### 

svyCreateTableOne(vars=c("consult_mh"), factorVars = c("consult_mh"),strata="survey",data=svyds.ALP.CC.all)

cchsSNall2 <- cchsSNall %>% 
  dplyr::select(-c(HIV_test12m,STI_test12m,lastcondom_use,insurance,marijuana))

alldat <- cchsSNall2 %>% drop_na(age_grp,income,education,orient,msm,illicit_drug,have_pcp,smoke,employ,province,rural,race,existing_mh,immigration,alcohol)

#Load the bootstrap from CCHS
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

#Get the point estimates first
#Check if in the other program CCHS 2017 also very low prevalences
alldat_gay_cchs <- cchsbsw %>% dplyr::filter(orient==2)
alldat_gay_sn <- sndat_all_bsw %>% dplyr::filter(orient==2)

alldat_gay <- alldat %>% dplyr::filter(orient == 2)

svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat_gay) ### 
Formula_fit = as.formula("survey ~ age_grp + age_grp^2 + income + smoke + education + illicit_drug + employ + smoke + alcohol + province + age_grp*income + age_grp*employ + rural + race + existing_mh + immigration") 
#have_pcp is removed but this is the best performing one
lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)    
beta.w = summary(lgtreg.w)$coeff[,1]   
alldat_gay$ps.w = lgtreg.w$fitted.values 

alldat_gay <- alldat_gay %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M_rescaled,ALP_weight)) 

##Check if the weights are correctly applied
svyds.ALP_gay  = svydesign(ids=~1, weight = ~ ALP_weight, data = alldat_gay) ### 

cchs_gay <- alldat_gay %>% filter(survey == 0)
SN_gay <- alldat_gay %>% filter(survey == 1)

tabw<-svyCreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=svyds.ALP_gay,strata="survey",test=FALSE,smd = TRUE)
tab2<-print(tabw,smd=TRUE)
#very well balanced, so no changes necesarry 

#Check consult_mh, SN weighted, CCHS weighted, SN unweighted
svyds.ALP.SN.w  = svydesign(ids=~1, weight = ~ ALP_weight, data = SN_gay) ### 
svyds.ALP.CC.w  = svydesign(ids=~1, weight = ~ ALP_weight, data = cchs_gay) ### 

point.est.gay.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w)

point.est.gay.cchs<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh"), data=svyds.ALP.CC.w)

#unadjusted SN
CreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","srmh","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","srmh"), data=SN_gay)

nfeature = 7
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
  svyprev<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w)
  prev.list[[1]]<-c(prev.list[[1]],svyprev$CatTable$Overall$consult_mh$percent[2])
  prev.list[[2]]<-c(prev.list[[2]],svyprev$CatTable$Overall$depress_binary$percent[2])
  prev.list[[3]]<-c(prev.list[[3]],svyprev$CatTable$Overall$anxiety_binary$percent[2])
  prev.list[[4]]<-c(prev.list[[4]],svyprev$CatTable$Overall$poor_srmh$percent[2])
  prev.list[[5]]<-c(prev.list[[5]],svyprev$CatTable$Overall$loneliness_ucla_bin$percent[2])
  prev.list[[6]]<-c(prev.list[[6]],svyprev$ContTable$Overall[,"mean"])
  prev.list[[7]]<-c(prev.list[[7]],svyprev$ContTable$Overall[,"median"])
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
# > time.taken
# Time difference of 49.22 mins

prev.ci<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))
                  
prev.ci$variable <- c("consult mh","depression","anxiety","Poor SRMH","Loneliness binary","Loneliness Mean","Loneliness Median")
prev.ci$estimate <-c(point.est.gay.adj$CatTable$Overall$consult_mh$percent[2],point.est.gay.adj$CatTable$Overall$depress_binary$percent[2],
                     point.est.gay.adj$CatTable$Overall$anxiety_binary$percent[2],point.est.gay.adj$CatTable$Overall$poor_srmh$percent[2],
                     point.est.gay.adj$CatTable$Overall$loneliness_ucla_bin$percent[2],point.est.gay.adj$ContTable$Overall[,"mean"],point.est.gay.adj$ContTable$Overall[,"median"])

for (k in 1:nfeature) {
  prev.ci$lower[k]<-quantile(prev.list[[k]],.025)
  prev.ci$upper[k]<-quantile(prev.list[[k]],.975)
}
prev.ci

#             variable     lower  estimate    upper
# 1        consult mh 23.231537 30.583100 39.97747
# 2        depression 10.100880 14.775540 25.25119
# 3           anxiety 10.577163 16.869158 22.29837
# 4         Poor SRMH 12.784813 18.166529 26.93676
# 5 Loneliness binary 42.474008 48.411116 59.32547
# 6   Loneliness Mean  5.130182  5.385797  5.97049
# 7 Loneliness Median  5.000000  5.000000  6.00000

CreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=SN_gay)
#Calculate CI for unadjusted estimate from SN, for consistenty, bootstrap all estimates

depSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(consult_mh) %>% select(consult_mh)))
anxietySN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(anxiety_binary) %>% select(anxiety_binary)))
srmhSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(poor_srmh) %>% select(poor_srmh)))
lonelinessbinSN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(loneliness_ucla_bin) %>% select(loneliness_ucla_bin)))
lonelinessUCLASN <- as.numeric(unlist(alldat_gay_sn %>% drop_na(loneliness_ucla) %>% select(loneliness_ucla)))

# anxietySN and srmhSN are already 0, 1, so no need for -1

consultSN <- consultSN - 1
depSN <- depSN - 1

fc <- function(d, i){
  d2 <- d[i]
  return(mean(d2))
}

b1<-boot(data=consultSN,fc,R=B)
mean(consultSN)*100
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN,fc,R=B)
mean(depSN)*100
boot.ci(b2, conf = 0.95,type = "perc")

b3<-boot(data=anxietySN,fc,R=B)
mean(anxietySN)*100
boot.ci(b3, conf = 0.95,type = "perc")

b4<-boot(data=srmhSN,fc,R=B)
mean(srmhSN)*100
boot.ci(b4, conf = 0.95,type = "perc")

b5<-boot(data=lonelinessbinSN,fc,R=B)
mean(lonelinessbinSN)*100
boot.ci(b5, conf = 0.95,type = "perc")

b6<-boot(data=lonelinessUCLASN,fc,R=B)
mean(lonelinessUCLASN)
boot.ci(b6, conf = 0.95,type = "perc")

#Try to generalize to all MSM from CCHS

svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat) 

Formula_fit = as.formula("survey ~ age_grp + age_grp^2 + income + income^2 + age_grp*education + age_grp*employ + age_grp*income + income*province + smoke + education + illicit_drug + employ + smoke + alcohol + province + rural + race + existing_mh + immigration") 

lgtreg.w = svyglm(Formula_fit, family= quasibinomial, design = svyds.wei)    
beta.w = summary(lgtreg.w)$coeff[,1]   
alldat$ps.w = lgtreg.w$fitted.values 
# ps.w.c = alldat$ps.w[which(alldat$survey==1)] 
# ALP_wei = as.vector((1-ps.w.c)/ps.w.c)              

alldat <- alldat %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M_rescaled,ALP_weight)) 

##Check if the weights are correctly applied
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

point.est.bi.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.bi)
point.est.other.adj<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.other)
point.est.allmsm<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.msm)

#unadjusted SN
CreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=allsn_bi)
CreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=allsn_other)
CreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla","loneliness_ucla_bin"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=allsn)

nfeature = 7
prev.list.bi<-vector("list",nfeature)
prev.list.other<-vector("list",nfeature)
prev.list.all.msm<-vector("list",nfeature)

#So now this is the bootstrap for all MSM, for gays it was different
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

  svyprev.bi<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.bi)
  svyprev.other<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.other)
  svyprev.msm<-svyCreateTableOne(vars=c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin","loneliness_ucla"), factorVars = c("consult_mh","depress_binary","anxiety_binary","poor_srmh","loneliness_ucla_bin"), data=svyds.ALP.SN.w.msm)
 
  prev.list.bi[[1]]<-c(prev.list.bi[[1]],svyprev.bi$CatTable$Overall$consult_mh$percent[2])
  prev.list.bi[[2]]<-c(prev.list.bi[[2]],svyprev.bi$CatTable$Overall$depress_binary$percent[2])
  prev.list.bi[[3]]<-c(prev.list.bi[[3]],svyprev.bi$CatTable$Overall$anxiety_binary$percent[2])
  prev.list.bi[[4]]<-c(prev.list.bi[[4]],svyprev.bi$CatTable$Overall$poor_srmh$percent[2])
  prev.list.bi[[5]]<-c(prev.list.bi[[5]],svyprev.bi$CatTable$Overall$loneliness_ucla_bin$percent[2])
  prev.list.bi[[6]]<-c(prev.list.bi[[6]],svyprev.bi$ContTable$Overall[,"mean"])
  prev.list.bi[[7]]<-c(prev.list.bi[[7]],svyprev.bi$ContTable$Overall[,"median"])
  
  prev.list.other[[1]]<-c(prev.list.other[[1]],svyprev.other$CatTable$Overall$consult_mh$percent[2])
  prev.list.other[[2]]<-c(prev.list.other[[2]],svyprev.other$CatTable$Overall$depress_binary$percent[2])
  prev.list.other[[3]]<-c(prev.list.other[[3]],svyprev.other$CatTable$Overall$anxiety_binary$percent[2])
  prev.list.other[[4]]<-c(prev.list.other[[4]],svyprev.other$CatTable$Overall$poor_srmh$percent[2])
  prev.list.other[[5]]<-c(prev.list.other[[5]],svyprev.other$CatTable$Overall$loneliness_ucla_bin$percent[2])
  prev.list.other[[6]]<-c(prev.list.other[[6]],svyprev.other$ContTable$Overall[,"mean"])
  prev.list.other[[7]]<-c(prev.list.other[[7]],svyprev.other$ContTable$Overall[,"median"])
  
  prev.list.all.msm[[1]]<-c(prev.list.all.msm[[1]],svyprev.msm$CatTable$Overall$consult_mh$percent[2])
  prev.list.all.msm[[2]]<-c(prev.list.all.msm[[2]],svyprev.msm$CatTable$Overall$depress_binary$percent[2])
  prev.list.all.msm[[3]]<-c(prev.list.all.msm[[3]],svyprev.msm$CatTable$Overall$anxiety_binary$percent[2])
  prev.list.all.msm[[4]]<-c(prev.list.all.msm[[4]],svyprev.msm$CatTable$Overall$poor_srmh$percent[2])
  prev.list.all.msm[[5]]<-c(prev.list.all.msm[[5]],svyprev.msm$CatTable$Overall$loneliness_ucla_bin$percent[2])
  prev.list.all.msm[[6]]<-c(prev.list.all.msm[[6]],svyprev.msm$ContTable$Overall[,"mean"])
  prev.list.all.msm[[7]]<-c(prev.list.all.msm[[7]],svyprev.msm$ContTable$Overall[,"median"])
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
# > time.taken
# Time difference of 3.95 hours

prev.ci.bi<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))
    
prev.ci.bi$variable <- c("consult mh","depression","anxiety","Poor SRMH","Loneliness binary","Loneliness Mean","Loneliness Median")
prev.ci.bi$estimate <-c(point.est.bi.adj$CatTable$Overall$consult_mh$percent[2],point.est.bi.adj$CatTable$Overall$depress_binary$percent[2],
                        point.est.bi.adj$CatTable$Overall$anxiety_binary$percent[2],point.est.bi.adj$CatTable$Overall$poor_srmh$percent[2],
                        point.est.bi.adj$CatTable$Overall$loneliness_ucla_bin$percent[2],point.est.bi.adj$ContTable$Overall[,"mean"],point.est.bi.adj$ContTable$Overall[,"median"])

for (k in 1:nfeature) {
  prev.ci.bi$lower[k]<-quantile(prev.list.bi[[k]],.025)
  prev.ci.bi$upper[k]<-quantile(prev.list.bi[[k]],.975)
}
prev.ci.bi
# > prev.ci.bi
#             variable     lower  estimate     upper
# 1        consult mh 21.897325 33.785843 40.543372
# 2        depression 11.190990 18.098241 24.854441
# 3           anxiety 18.599272 27.830158 52.354655
# 4         Poor SRMH 21.492357 31.897641 53.292034
# 5 Loneliness binary 32.495815 47.850582 55.590103
# 6   Loneliness Mean  5.045932  5.399121  5.667292
# 7 Loneliness Median  5.000000  5.000000  6.000000

prev.ci.other<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.other$variable <- c("consult mh","depression","anxiety","Poor SRMH","Loneliness binary","Loneliness Mean","Loneliness Median")
prev.ci.other$estimate <-c(point.est.other.adj$CatTable$Overall$consult_mh$percent[2],point.est.other.adj$CatTable$Overall$depress_binary$percent[2],
                           point.est.other.adj$CatTable$Overall$anxiety_binary$percent[2],point.est.other.adj$CatTable$Overall$poor_srmh$percent[2],
                           point.est.other.adj$CatTable$Overall$loneliness_ucla_bin$percent[2],point.est.other.adj$ContTable$Overall[,"mean"],point.est.other.adj$ContTable$Overall[,"median"])


for (k in 1:nfeature) {
  prev.ci.other$lower[k]<-quantile(prev.list.other[[k]],.025)
  prev.ci.other$upper[k]<-quantile(prev.list.other[[k]],.975)
}
prev.ci.other

# > prev.ci.other
#           variable     lower  estimate     upper
# 1        consult mh 11.686983 20.335945 29.391743
# 2        depression  6.588438 18.389619 36.667756
# 3           anxiety  9.896555 21.124713 39.415299
# 4         Poor SRMH 13.309137 26.888196 44.042402
# 5 Loneliness binary 30.973294 45.391992 60.262320
# 6   Loneliness Mean  4.775667  5.279517  5.850465
# 7 Loneliness Median  4.000000  5.000000  6.000000

prev.ci.allmsm<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.allmsm$variable <- c("consult mh","depression","anxiety","Poor SRMH","Loneliness binary","Loneliness Mean","Loneliness Median")
prev.ci.allmsm$estimate <-c(point.est.allmsm$CatTable$Overall$consult_mh$percent[2],point.est.allmsm$CatTable$Overall$depress_binary$percent[2],
                            point.est.allmsm$CatTable$Overall$anxiety_binary$percent[2],point.est.allmsm$CatTable$Overall$poor_srmh$percent[2],
                            point.est.allmsm$CatTable$Overall$loneliness_ucla_bin$percent[2],point.est.allmsm$ContTable$Overall[,"mean"],point.est.allmsm$ContTable$Overall[,"median"])

for (k in 1:nfeature) {
  prev.ci.allmsm$lower[k]<-quantile(prev.list.all.msm[[k]],.025)
  prev.ci.allmsm$upper[k]<-quantile(prev.list.all.msm[[k]],.975)
}
prev.ci.allmsm

# > prev.ci.allmsm
#             variable     lower  estimate     upper
# 1        consult mh 26.189296 32.127218 42.172454
# 2        depression 11.604426 15.959957 23.979131
# 3           anxiety 13.679656 18.828822 27.946946
# 4         Poor SRMH 14.963663 20.768159 30.400540
# 5 Loneliness binary 42.119989 48.449003 55.242907
# 6   Loneliness Mean  5.168873  5.396254  5.703777
# 7 Loneliness Median  5.000000  5.000000  6.000000

#Also bootstrap unadjusted estimates for bi and other MSM from SN

fc <- function(d, i){
  d2 <- d[i]
  return(mean(d2))
}

depSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(consult_mh) %>% select(consult_mh)))
anxietySN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(anxiety_binary) %>% select(anxiety_binary)))
srmhSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(poor_srmh) %>% select(poor_srmh)))
lonelinessbinSN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(loneliness_ucla_bin) %>% select(loneliness_ucla_bin)))
lonelinessUCLASN.bi <- as.numeric(unlist(allsn_bi %>% drop_na(loneliness_ucla) %>% select(loneliness_ucla)))

consultSN.bi <- consultSN.bi - 1
depSN.bi <- depSN.bi - 1

b1<-boot(data=consultSN.bi,fc,R=B)
mean(consultSN.bi)*100
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.bi,fc,R=B)
mean(depSN.bi)*100
boot.ci(b2, conf = 0.95,type = "perc")

b3<-boot(data=anxietySN.bi,fc,R=B)
mean(anxietySN.bi)*100
boot.ci(b3, conf = 0.95,type = "perc")

b4<-boot(data=srmhSN.bi,fc,R=B)
mean(srmhSN.bi)*100
boot.ci(b4, conf = 0.95,type = "perc")

b5<-boot(data=lonelinessbinSN.bi,fc,R=B)
mean(lonelinessbinSN.bi)*100
boot.ci(b5, conf = 0.95,type = "perc")

b6<-boot(data=lonelinessUCLASN.bi,fc,R=B)
mean(lonelinessUCLASN.bi)
boot.ci(b6, conf = 0.95,type = "perc")

depSN.other <- as.numeric(unlist(allsn_other %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.other <- as.numeric(unlist(allsn_other %>% drop_na(consult_mh) %>% select(consult_mh)))
anxietySN.other <- as.numeric(unlist(allsn_other %>% drop_na(anxiety_binary) %>% select(anxiety_binary)))
srmhSN.other <- as.numeric(unlist(allsn_other %>% drop_na(poor_srmh) %>% select(poor_srmh)))
lonelinessbinSN.other <- as.numeric(unlist(allsn_other %>% drop_na(loneliness_ucla_bin) %>% select(loneliness_ucla_bin)))
lonelinessUCLASN.other <- as.numeric(unlist(allsn_other %>% drop_na(loneliness_ucla) %>% select(loneliness_ucla)))

consultSN.other <- consultSN.other - 1
depSN.other <- depSN.other - 1

b1<-boot(data=consultSN.other,fc,R=B)
mean(consultSN.other)*100
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.other,fc,R=B)
mean(depSN.other)*100
boot.ci(b2, conf = 0.95,type = "perc")

b3<-boot(data=anxietySN.other,fc,R=B)
mean(anxietySN.other)*100
boot.ci(b3, conf = 0.95,type = "perc")

b4<-boot(data=srmhSN.other,fc,R=B)
mean(srmhSN.other)*100
boot.ci(b4, conf = 0.95,type = "perc")

b5<-boot(data=lonelinessbinSN.other,fc,R=B)
mean(lonelinessbinSN.other)*100
boot.ci(b5, conf = 0.95,type = "perc")

b6<-boot(data=lonelinessUCLASN.other,fc,R=B)
mean(lonelinessUCLASN.other)
boot.ci(b6, conf = 0.95,type = "perc")

depSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(depress_binary) %>% select(depress_binary)))
consultSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(consult_mh) %>% select(consult_mh)))
anxietySN.allmsm <- as.numeric(unlist(allsn %>% drop_na(anxiety_binary) %>% select(anxiety_binary)))
srmhSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(poor_srmh) %>% select(poor_srmh)))
lonelinessbinSN.allmsm <- as.numeric(unlist(allsn %>% drop_na(loneliness_ucla_bin) %>% select(loneliness_ucla_bin)))
lonelinessUCLASN.allmsm <- as.numeric(unlist(allsn %>% drop_na(loneliness_ucla) %>% select(loneliness_ucla)))

consultSN.allmsm <- consultSN.allmsm - 1
depSN.allmsm <- depSN.allmsm - 1

b1<-boot(data=consultSN.allmsm,fc,R=B)
mean(consultSN.allmsm)*100
boot.ci(b1, conf = 0.95,type = "perc")

b2<-boot(data=depSN.allmsm,fc,R=B)
mean(depSN.allmsm)*100
boot.ci(b2, conf = 0.95,type = "perc")

b3<-boot(data=anxietySN.allmsm,fc,R=B)
mean(anxietySN.allmsm)*100
boot.ci(b3, conf = 0.95,type = "perc")

b4<-boot(data=srmhSN.allmsm,fc,R=B)
mean(srmhSN.allmsm)*100
boot.ci(b4, conf = 0.95,type = "perc")

b5<-boot(data=lonelinessbinSN.allmsm,fc,R=B)
mean(lonelinessbinSN.allmsm)*100
boot.ci(b5, conf = 0.95,type = "perc")

b6<-boot(data=lonelinessUCLASN.allmsm,fc,R=B)
mean(lonelinessUCLASN.allmsm)
boot.ci(b6, conf = 0.95,type = "perc")

###CCHS bootstrap weights
svyds.ALP.CC.w  = svydesign(ids=~1, weight = ~ WTS_M_rescaled, data = alldat_cchs) ### 
point.est.cchs.all<-svyCreateTableOne(vars=c("consult_mh","depress_binary","poor_srmh"), factorVars = c("consult_mh","depress_binary","poor_srmh"), data=svyds.ALP.CC.w)
point.est.cchs.orient<-svyCreateTableOne(vars=c("consult_mh","depress_binary","poor_srmh"), factorVars = c("consult_mh","depress_binary","poor_srmh"), strata="orient",data=svyds.ALP.CC.w)

nfeature = 3
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
  cchs.boot.orient<-svyCreateTableOne(vars=c("consult_mh","depress_binary","poor_srmh"), factorVars = c("consult_mh","depress_binary","poor_srmh"), strata="orient",data=svyds.ALP.cc.all)
  cchs.boot.all<-svyCreateTableOne(vars=c("consult_mh","depress_binary","poor_srmh"), factorVars = c("consult_mh","depress_binary","poor_srmh"),data=svyds.ALP.cc.all)
  prev.list.cchs.all[[1]]<-c(prev.list.cchs.all[[1]],cchs.boot.all$CatTable$Overall$consult_mh$percent[2])
  prev.list.cchs.all[[2]]<-c(prev.list.cchs.all[[2]],cchs.boot.all$CatTable$Overall$depress_binary$percent[2])
  prev.list.cchs.all[[3]]<-c(prev.list.cchs.all[[3]],cchs.boot.all$CatTable$Overall$poor_srmh$percent[2])
  
  prev.list.cchs.other[[1]]<-c(prev.list.cchs.other[[1]],cchs.boot.orient$CatTable$"1"$consult_mh$percent[2])
  prev.list.cchs.other[[2]]<-c(prev.list.cchs.other[[2]],cchs.boot.orient$CatTable$"1"$depress_binary$percent[2])
  prev.list.cchs.other[[3]]<-c(prev.list.cchs.other[[3]],cchs.boot.orient$CatTable$"1"$poor_srmh$percent[2])
  
  prev.list.cchs.gay[[1]]<-c(prev.list.cchs.gay[[1]],cchs.boot.orient$CatTable$"2"$consult_mh$percent[2])
  prev.list.cchs.gay[[2]]<-c(prev.list.cchs.gay[[2]],cchs.boot.orient$CatTable$"2"$depress_binary$percent[2])
  prev.list.cchs.gay[[3]]<-c(prev.list.cchs.gay[[3]],cchs.boot.orient$CatTable$"2"$poor_srmh$percent[2])
  
  prev.list.cchs.bi[[1]]<-c(prev.list.cchs.bi[[1]],cchs.boot.orient$CatTable$"3"$consult_mh$percent[2])
  prev.list.cchs.bi[[2]]<-c(prev.list.cchs.bi[[2]],cchs.boot.orient$CatTable$"3"$depress_binary$percent[2])
  prev.list.cchs.bi[[3]]<-c(prev.list.cchs.bi[[3]],cchs.boot.orient$CatTable$"3"$poor_srmh$percent[2])
  
}
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

# > time.taken
# Time difference of 1.5 hours

#This one is for the point estimates for CCHS

prev.ci.allmsm.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.allmsm.cchs$variable <- c("consult mh","depression","poor_srmh")
prev.ci.allmsm.cchs$estimate <-c(point.est.cchs.all$CatTable$Overall$consult_mh$percent[2],point.est.cchs.all$CatTable$Overall$depress_binary$percent[2],
                                 point.est.cchs.all$CatTable$Overall$poor_srmh$percent[2])

for (k in 1:nfeature) {
  prev.ci.allmsm.cchs$lower[k]<-quantile(prev.list.cchs.all[[k]],.025)
  prev.ci.allmsm.cchs$upper[k]<-quantile(prev.list.cchs.all[[k]],.975)
}
prev.ci.allmsm.cchs
# > prev.ci.allmsm.cchs
# variable     lower  estimate    upper
# 1 consult mh 15.094595 21.593910 29.25672
# 2 depression  3.727183  9.124652 18.03423
# 3  poor_srmh  6.975984 12.058501 19.40676

prev.ci.gay.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.gay.cchs$variable <- c("consult mh","depression","poor_srmh")
prev.ci.gay.cchs$estimate <-c(point.est.cchs.orient$CatTable$"2"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"2"$depress_binary$percent[2],
                              point.est.cchs.orient$CatTable$"2"$poor_srmh$percent[2])

for (k in 1:nfeature) {
  prev.ci.gay.cchs$lower[k]<-quantile(prev.list.cchs.gay[[k]],.025)
  prev.ci.gay.cchs$upper[k]<-quantile(prev.list.cchs.gay[[k]],.975)
}
prev.ci.gay.cchs
# > prev.ci.gay.cchs
# variable     lower  estimate    upper
# 1 consult mh 12.951667 21.336778 32.47429
# 2 depression  2.568024  9.609042 23.03188
# 3  poor_srmh  5.001418 10.825413 18.89602

prev.ci.bi.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.bi.cchs$variable <- c("consult mh","depression","poor_srmh")
prev.ci.bi.cchs$estimate <-c(point.est.cchs.orient$CatTable$"3"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"3"$depress_binary$percent[2],
                             point.est.cchs.orient$CatTable$"3"$poor_srmh$percent[2])

for (k in 1:nfeature) {
  prev.ci.bi.cchs$lower[k]<-quantile(prev.list.cchs.bi[[k]],.025)
  prev.ci.bi.cchs$upper[k]<-quantile(prev.list.cchs.bi[[k]],.975)
}
prev.ci.bi.cchs
# > prev.ci.bi.cchs
# variable     lower estimate    upper
# 1 consult mh 12.631851 22.37722 33.98534
# 2 depression  2.042538 10.77256 26.80901
# 3  poor_srmh  6.491465 14.13642 23.62245

prev.ci.other.cchs<-setNames(data.frame(matrix(ncol = 4, nrow = nfeature)), c("variable","lower","estimate","upper"))

prev.ci.other.cchs$variable <- c("consult mh","depression","poor_srmh")
prev.ci.other.cchs$estimate <-c(point.est.cchs.orient$CatTable$"1"$consult_mh$percent[2],point.est.cchs.orient$CatTable$"1"$depress_binary$percent[2],
                                point.est.cchs.orient$CatTable$"1"$poor_srmh$percent[2])

for (k in 1:nfeature) {
  prev.ci.other.cchs$lower[k]<-quantile(prev.list.cchs.other[[k]],.025)
  prev.ci.other.cchs$upper[k]<-quantile(prev.list.cchs.other[[k]],.975)
}
prev.ci.other.cchs
# > prev.ci.other.cchs
# variable   lower   estimate     upper
# 1 consult mh 3.30538 20.7698766 47.692808
# 2 depression 0.00000  0.6701766  7.136573
# 3  poor_srmh 0.00000 13.7803552 52.935725

tab5<-CreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","age_grp","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=alldat_sn)
tab5<-print(tab5)
tab6<-CreateTableOne(vars=c("immigration","age_grp","education","income","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), factorVars = c("immigration","education","age_grp","illicit_drug","have_pcp","smoke","employ","province","race","existing_mh","rural","alcohol"), data=SN_gay)
tab6<-print(tab6)

