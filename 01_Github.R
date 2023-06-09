###Program: 01_explore_PUMF

##Create dataset from PUMF and Sex Now original data

library(haven)
library(dplyr)
library(naniar)

cchs2015<-read.csv("cchs-82M0013-E-2015-2016-Annual-component_F1.csv")

nyear <- 2
cchsdat2015 <- cchs2015 %>% dplyr::filter(DHH_SEX == 1 & SDC_035 %in% c(1,2,3) & DHHGAGE >= 3) %>%
  dplyr::select(c(SDCDVIMM,SDCDGRES,EHG2DVR3,INCDVRCA,SDC_035,DHHGAGE,SXB_060,DEPDVSEV,GENDVHDI,GENDVMHI,GEO_PRV,GEODGHR4,SMKDVSTY,SDCDGCGT,MAC_005,LBFDVWSS,LBFDVPFT,LBFG10,
                  DRMDVLAY,PHC_020,INS_005,ALC_020,ALC_010,SXB_005,SXB_110,SXB_115,SXB_080,SXB_130,CCC_195,CCC_200,CMH_005,DRGDVYCM,WTS_M,SXBDVTST,ADM_RNO)) %>% 
  dplyr::mutate(dplyr::select(.,c(EHG2DVR3,SDC_035,SXB_060,DEPDVSEV,SDCDGRES,CCC_195,ALC_010,CCC_200,CMH_005,SDCDGCGT,INS_005,PHC_020,DRMDVLAY,SXB_080,GENDVMHI,GENDVHDI,DRGDVYCM,SXBDVTST)) %>% 
           replace_with_na_all(condition = ~.x %in% c(6,7,8,9))) %>%
  dplyr::mutate(select(.,c(INCDVRCA,ALC_020,SMKDVSTY)) %>% 
           replace_with_na_all(condition = ~.x %in% c(96,97,98,99))) %>%
  dplyr::rename(orient = SDC_035, age_grp = DHHGAGE, consult_mh = CMH_005, grh = GENDVHDI, srmh = GENDVMHI, income = INCDVRCA,
                insurance = INS_005, marijuana = DRGDVYCM, have_pcp = PHC_020, illicit_drug = DRMDVLAY, education = EHG2DVR3) %>%
  dplyr::mutate(ALC_020 = ifelse(ALC_010 == 2,1,ALC_020),
                SXB_110 = ifelse(SXB_005 == 2,0,SXB_110),
                SXB_130 = ifelse(SXB_005 == 2,0,SXB_130),
                SXB_080 = ifelse(SXB_005 == 2,0,SXB_080),
                SXBDVTST = ifelse(SXB_005 == 2,0,SXBDVTST)) %>%
  dplyr::mutate(immigration = case_when(SDCDVIMM==2 ~ 0,
                                        SDCDGRES==1 ~ 1,
                                        SDCDGRES==2 ~ 2),
                GEODGHR4 = paste(substr(GEODGHR4,1,2),substr(GEODGHR4,4,5),sep=""),
                rural = ifelse(GEODGHR4 %in% c("1204","1302","2403","2413","3527","3537","3544","3546","3568","4601","5913","5933","5941",
                                               "2407","2414","2415","2416","3530","3536","3551","3560","3565","3566","4605","4704","4706",
                                               "4710","4831","4832","4833","4834","5921","2406","3595","5932","3553","3570","5922","5923","5931"),0,1),
                grh = grh + 1,
                srmh = srmh + 1,
         existing_mh = case_when(CCC_195 == 1 | CCC_200 == 1 ~ 1,
                                 CCC_195 == 2 | CCC_200 == 2 ~ 0),
         depress_binary = case_when(DEPDVSEV >= 3 ~ 1,
                                    DEPDVSEV < 3 ~ 0),
         consult_mh = ifelse(consult_mh == 2,0,consult_mh),
         #race == 0 white 1 = non white in CCHS 1 = white 2 = non white
         race = case_when(SDCDGCGT == 2 ~ 1,
                          SDCDGCGT == 1 ~ 0),
         msm = case_when(orient %in% c(2,3) | SXB_060 == 1 ~ 1,
                         orient == 1 ~ 0),
         alcohol = case_when(ALC_020 == 1 ~ 0,
                             ALC_020 > 1 ~ 1),
         disclosure = case_when(msm == 1 ~ 1,
                                msm == 0 ~ 0),
         employ = case_when(MAC_005 == 4 ~ 4,
                            MAC_005 == 7 ~ 6,
                            MAC_005 %in% c(3,6,9) | LBFDVWSS == 3 ~ 5, 
                            LBFG10 == 2 ~ 3,
                            LBFDVPFT == 1 ~ 1,
                            LBFDVPFT == 2 ~ 2),
         province = case_when(GEO_PRV == 59 ~ "WestCoast",
                              GEO_PRV == 48 ~ "Prairies",
                              GEO_PRV == 46 ~ "Prairies",
                              GEO_PRV == 47 ~ "Prairies",
                              GEO_PRV == 35 ~ "Central",
                              GEO_PRV == 24 ~ "Central",
                              GEO_PRV == 13 ~ "Atlantic",
                              GEO_PRV == 11 ~ "Atlantic",
                              GEO_PRV == 12 ~ "Atlantic",
                              GEO_PRV == 10 ~ "Atlantic",
                              GEO_PRV == 61 ~ "Northern",
                              GEO_PRV == 62 ~ "Northern",
                              GEO_PRV == 60 ~ "Northern"),
         smoke = case_when(SMKDVSTY %in% c(1,2) ~ 1,
                           SMKDVSTY %in% c(3,4,5,6) ~ 0),
         HIV_test12m = case_when(SXB_115 %in% c(1,2,3) ~ 1,
                                 SXB_110 %in% c(0,2) | SXB_115 == 4 ~ 0),
         STI_test12m = case_when(SXBDVTST %in% c(1,2,3) ~ 1,
                                 SXB_130 %in% c(0,2) | SXBDVTST == 4 ~ 0),
         insurance = ifelse(insurance == 2,0,insurance),
         illicit_drug = ifelse(illicit_drug == 2,0,illicit_drug),
         marijuana = ifelse(marijuana == 2,0,marijuana),
         have_pcp = ifelse(have_pcp == 2,0,have_pcp),
         lastcondom_use = ifelse(SXB_080 == 2,0,SXB_080),
         WTS_M_rescaled = WTS_M / nyear) %>%
  dplyr::mutate(survey = 0,year=2015) %>%
  dplyr::select(c(education,income,orient,age_grp,msm,grh,srmh,illicit_drug,have_pcp,insurance,consult_mh,marijuana,smoke,employ,province,rural,ADM_RNO,
                  race,depress_binary,existing_mh,immigration,alcohol,HIV_test12m,STI_test12m,lastcondom_use,survey,disclosure,year,WTS_M,WTS_M_rescaled))

cchs2017<-read.csv("cchs-82M0013-E-2017-2018-Annual-component_F1.csv")

cchsdat2017 <- cchs2017 %>% dplyr::filter(DHH_SEX == 1 & SDC_035 %in% c(1,2,3) & DHHGAGE >= 3) %>%
  dplyr::select(c(SDCDVIMM,SDCDGRES,EHG2DVR3,INCDVRCA,SDC_035,DHHGAGE,DEPDVSEV,GENDVHDI,GENDVMHI,GEO_PRV,GEODGHR4,SMKDVSTY,SDCDGCGT,SDCDGRES,
                  DRMDVLAY,PHC_020,INS_005,ALC_020,ALC_010,CCC_195,CCC_200,CMH_005,DRGDVYCM,MACG005,LBFDVWSS,LBFG10,LBFDVPFT,WTS_M,ADM_RNO)) %>%
  dplyr::mutate(dplyr::select(.,c(EHG2DVR3,SDC_035,DEPDVSEV,SDCDGRES,CCC_195,ALC_010,CCC_200,CMH_005,SDCDGCGT,INS_005,PHC_020,DRMDVLAY,GENDVMHI,GENDVHDI,DRGDVYCM)) %>% 
                  replace_with_na_all(condition = ~.x %in% c(6,7,8,9))) %>%
  dplyr::mutate(select(.,c(INCDVRCA,ALC_020,SMKDVSTY)) %>% 
                  replace_with_na_all(condition = ~.x %in% c(96,97,98,99))) %>%
  dplyr::rename(orient = SDC_035, age_grp = DHHGAGE, consult_mh = CMH_005, grh = GENDVHDI, srmh = GENDVMHI, income = INCDVRCA,
                insurance = INS_005, marijuana = DRGDVYCM, have_pcp = PHC_020, illicit_drug = DRMDVLAY, education = EHG2DVR3) %>%
  dplyr::mutate(ALC_020 = ifelse(ALC_010 == 2,1,ALC_020)) %>%
  dplyr::mutate(immigration = case_when(SDCDVIMM==2 ~ 0,
                                        SDCDGRES==1 ~ 1,
                                        SDCDGRES==2 ~ 2),
                grh = grh + 1,
                srmh = srmh + 1,
                GEODGHR4 = paste(substr(GEODGHR4,1,2),substr(GEODGHR4,4,5),sep=""),
                rural = ifelse(GEODGHR4 %in% c("1204","1302","2403","2413","3527","3537","3544","3546","3568","4601","5913","5933","5941",
                                               "2407","2414","2415","2416","3530","3536","3551","3560","3565","3566","4605","4704","4706",
                                               "4710","4831","4832","4833","4834","5921","2406","3595","5932","3553","3570","5922","5923","5931"),0,1),
                province = case_when(GEO_PRV == 59 ~ "WestCoast",
                                     GEO_PRV == 48 ~ "Prairies",
                                     GEO_PRV == 46 ~ "Prairies",
                                     GEO_PRV == 47 ~ "Prairies",
                                     GEO_PRV == 35 ~ "Central",
                                     GEO_PRV == 24 ~ "Central",
                                     GEO_PRV == 13 ~ "Atlantic",
                                     GEO_PRV == 11 ~ "Atlantic",
                                     GEO_PRV == 12 ~ "Atlantic",
                                     GEO_PRV == 10 ~ "Atlantic",
                                     GEO_PRV == 61 ~ "Northern",
                                     GEO_PRV == 62 ~ "Northern",
                                     GEO_PRV == 60 ~ "Northern"),
                existing_mh = case_when(CCC_195 == 1 | CCC_200 == 1 ~ 1,
                                        CCC_195 == 2 | CCC_200 == 2 ~ 0),
                #Moderate or higher is the depression binary
                depress_binary = case_when(DEPDVSEV >= 3 ~ 1,
                                           DEPDVSEV < 3 ~ 0),
                consult_mh = ifelse(consult_mh == 2,0,consult_mh),
                #race == 0 white 1 = non white in CCHS 1 = white 2 = non white
                race = case_when(SDCDGCGT == 2 ~ 1,
                                 SDCDGCGT == 1 ~ 0),
                msm = case_when(orient %in% c(2,3) ~ 1,
                                orient == 1 ~ 0),
                alcohol = case_when(ALC_020 == 1 ~ 0,
                                    ALC_020 > 1 ~ 1),
                smoke = case_when(SMKDVSTY %in% c(1,2) ~ 1,
                                  SMKDVSTY %in% c(3,4,5,6) ~ 0),
                disclosure = case_when(msm == 1 ~ 1,
                                       msm == 0 ~ 0),
                employ = case_when(MACG005 == 3 ~ 4,
                                   MACG005 == 4 ~ 6,
                                   MACG005 %in% c(2,5) | LBFDVWSS == 3 ~ 5, 
                                   LBFG10 == 2 ~ 3,
                                   LBFDVPFT == 1 ~ 1,
                                   LBFDVPFT == 2 ~ 2),
                illicit_drug = ifelse(illicit_drug == 2,0,illicit_drug),
                marijuana = ifelse(marijuana == 2,0,marijuana),
                insurance = ifelse(insurance == 2,0,insurance),
                have_pcp = ifelse(have_pcp == 2,0,have_pcp),
                WTS_M_rescaled = WTS_M / nyear) %>%
  dplyr::mutate(survey = 0,year=2017,HIV_test12m = NA, STI_test12m = NA,lastcondom_use = NA) %>%
  dplyr::select(c(education,income,orient,age_grp,msm,grh,srmh,illicit_drug,have_pcp,insurance,consult_mh,marijuana,smoke,employ,province,rural,ADM_RNO,
                  race,depress_binary,existing_mh,immigration,alcohol,HIV_test12m,STI_test12m,lastcondom_use,survey,disclosure,year,WTS_M,WTS_M_rescaled))

cchsall <- rbind(cchsdat2015,cchsdat2017)

setwd("/Users/christofferdharma/Documents/SN2019/")
sn2019_text <- read.csv("Data_Sex Now 2019_text_incl.csv")
sn2019 <- read_sas("cleaned_data_quant.sas7bdat")

sn2019_text2 <- sn2019_text[,c("p0_1_anon_code","p47_199_last_male_anal_sex_condom","p47_199_last_male_anal_sex_no_condom")]

sn2019c <- merge(x = sn2019, y = sn2019_text2, by = "p0_1_anon_code", all.x = TRUE)
sn2019_clean <- sn2019c %>% filter(age >= 18) %>%
  mutate(age_grp = case_when(age >= 18 & age <= 19 ~ 3,
                             age >= 20 & age <= 24 ~ 4,
                             age >= 25 & age <= 29 ~ 5,
                             age >= 30 & age <= 34 ~ 6,
                             age >= 35 & age <= 39 ~ 7,
                             age >= 40 & age <= 44 ~ 8,
                             age >= 45 & age <= 49 ~ 9,
                             age >= 50 & age <= 54 ~ 10,
                             age >= 55 & age <= 59 ~ 11,
                             age >= 60 & age <= 64 ~ 12,
                             age >= 65 & age <= 69 ~ 13,
                             age >= 70 & age <= 74 ~ 14,
                             age >= 75 & age <= 79 ~ 15,
                             age >= 80 ~ 16),
       residence4 = case_when(p12_64_rurality == '1.Rural area (<1,000 people)' ~ 1,
                                p12_64_rurality == '2.Small city/town (1,000-29,999 people)' ~ 2,
                                p12_64_rurality == '3.Medium city/town (30,000-99,999 peopl' ~ 3,
                                p12_64_rurality == '4.Large urban centre (100,000+ people)' ~ 4),
       rural = case_when(residence4 == 4 ~ 0,
                         residence4 %in% c(1,2,3) ~ 1),
         orient = case_when(p7_42_sex_orientation_heterofle == 'Yes' ~ 1,
                            p7_42_sex_orientation_straight == 'Yes' ~ 1,
                            p7_42_sex_orientation_bi == 'Yes' ~ 3,
                            p7_42_sex_orientation_pansexual == 'Yes' ~ 3,
                            p7_42_sex_orientation_queer == 'Yes' ~ 3,
                            p7_42_sex_orientation_gay == 'Yes' ~ 2),
         immigration = case_when(immigration == 'Canadian born' ~ 0,
                                        immigration_year == '1. <= 2007' ~ 2,
                                        immigration_year %in% c('2. 2008-2013','3. 2014-2018') ~ 1),
         disclosure = case_when(p7_44_survey_reveal_sex_orienta == 'Likely' ~ 1,
                                p7_44_survey_reveal_sex_orienta == 'Very likely' ~ 1,
                                p7_44_survey_reveal_sex_orienta == 'Very unlikely' ~ 0,
                                p7_44_survey_reveal_sex_orienta == 'Unlikely' ~ 0),
         existing_mh = case_when(p6_37_disabilitiy_emotions == 'Always' ~ 1,
                                 p6_37_disabilitiy_emotions == 'Often' ~ 1,
                                 p6_37_disabilitiy_emotions == 'No (Never)' ~ 0,
                                 p6_37_disabilitiy_emotions == 'Sometimes' ~ 0),
         consult_mh = case_when(p37_145_service_counsellor == 'Yes' ~ 1,
                                p37_145_service_elder == 'Yes' ~ 1,
                                p37_145_service_knowledge == 'Yes' ~ 1,
                                p37_145_service_peer == 'Yes' ~ 1,
                                p37_145_service_psychiatrist == 'Yes' ~ 1,
                                p37_145_service_psychologist == 'Yes' ~ 1,
                                p37_145_service_sexologist == 'Yes' ~ 1,
                                p37_145_service_socialworker == 'Yes' ~ 1,
                                p37_145_service_none == 'Yes' ~ 0),
         alcohol = case_when(p40_155_substance_alcohol == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_alcohol %in% c('No','No drug use') ~ 0),
         benzo = case_when(p40_155_substance_benzo == 'Yes' ~ 1,
                           p39_151_substance_any == "No" ~ 0,
                           p40_155_substance_benzo %in% c('No','No drug use') ~ 0),
         marijuana = case_when(p40_155_substance_cannabis == 'Yes' ~ 1,
                              p39_151_substance_any == "No" ~ 0,
                              p40_155_substance_cannabis %in% c('No','No drug use') ~ 0),
         cocaine = case_when(p40_155_substance_cocaine == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_cocaine %in% c('No','No drug use') ~ 0),
         crack = case_when(p40_155_substance_crack == 'Yes' ~ 1,
                           p39_151_substance_any == "No" ~ 0,
                           p40_155_substance_crack %in% c('No','No drug use') ~ 0),
         crystal = case_when(p40_155_substance_crystal == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_crystal %in% c('No','No drug use') ~ 0),
         ecstasy = case_when(p40_155_substance_ecstasy == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_ecstasy %in% c('No','No drug use') ~ 0),
         fentanyl = case_when(p40_155_substance_fentanyl == 'Yes' ~ 1,
                              p39_151_substance_any == "No" ~ 0,
                              p40_155_substance_fentanyl %in% c('No','No drug use') ~ 0),
         GHB = case_when(p40_155_substance_GHB == 'Yes' ~ 1,
                         p39_151_substance_any == "No" ~ 0,
                         p40_155_substance_GHB %in% c('No','No drug use') ~ 0),
         heroin = case_when(p40_155_substance_heroin == 'Yes' ~ 1,
                            p39_151_substance_any == "No" ~ 0,
                            p40_155_substance_heroin %in% c('No','No drug use') ~ 0),
         ketamine = case_when(p40_155_substance_ketamine == 'Yes' ~ 1,
                              p39_151_substance_any == "No" ~ 0,
                              p40_155_substance_ketamine %in% c('No','No drug use') ~ 0),
         psychedelic = case_when(p40_155_substance_psych == 'Yes' ~ 1,
                                 p39_151_substance_any == "No" ~ 0,
                                 p40_155_substance_psych %in% c('No','No drug use') ~ 0),
         steroid = case_when(p40_155_substance_steroids == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_steroids %in% c('No','No drug use') ~ 0),
         smoke = case_when(p40_155_substance_tobacco == 'Yes' ~ 1,
                           p39_151_substance_any == "No" ~ 0,
                           p40_155_substance_tobacco %in% c('No','No drug use') ~ 0),
         poppers = case_when(p40_155_substance_poppers == 'Yes' ~ 1,
                             p39_151_substance_any == "No" ~ 0,
                             p40_155_substance_poppers %in% c('No','No drug use') ~ 0),
         insurance = case_when(p42_177_insurance_counselling == 'Yes' ~ 1,
                               p42_177_insurance_none == 'Yes' ~ 0,
                               p42_177_insurance_counselling == 'No' & p42_177_insurance_other == 'No' & p42_177_insurance_prescriptions == 'No' & p42_177_insurance_vaccinations == 'No' ~ 0,
                               p42_177_insurance_other == 'Yes' ~ 1,
                               p42_177_insurance_prescriptions == 'Yes' ~ 1,
                               p42_177_insurance_vaccinations == 'Yes' ~ 1),
         have_pcp = case_when(p42_179_regular_practitioner == 'Yes' ~ 1,
                              p42_179_regular_practitioner == 'No' ~ 0),
         grh = case_when(p42_170_general_health == 'Excellent' ~ 5,
                         p42_170_general_health == 'Very good' ~ 4,
                         p42_170_general_health == 'Good' ~ 3,
                         p42_170_general_health == 'Fair' ~ 2,
                         p42_170_general_health == 'Poor' ~ 1),
         # https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1110019301   
         #Income deciles are roughly similar
         income = case_when(p12_66_income %in% c('$10,000-$19,999','<$10,000') ~ 1,
                                  p12_66_income == '$20,000-$29,999' ~ 2,
                                  p12_66_income == '$30,000-$39,999' ~ 3,
                                  p12_66_income == '$40,000-$49,999' ~ 4,
                                  p12_66_income == '$50,000-$59,999' ~ 5,
                                  p12_66_income == '$60,000-$69,999' ~ 6,
                                  p12_66_income == '$70,000-$79,999' ~ 7,
                                  p12_66_income == '$80,000-$89,999' ~ 8,
                                  p12_66_income == '$90,000-$99,999' ~ 9,
                                  p12_66_income == '$100,000 or more' ~ 10),
         education = case_when(p12_68_education == 'Did not finish high school' ~ 1,
                               p12_68_education == 'High school, or equivalent' ~ 2,
                               p12_68_education %in% c('Post-secondary school (e.g. certificate, diploma)','Bachelors degree','Above a bachelors degree (e.g., masters, doctorate)') ~ 3),
         employ = case_when(p12_67_employment_fulltime == 'Yes' ~ 1,
                                p12_67_employment_parttime == 'Yes' ~ 2,
                                p12_67_employment_self == 'Yes' ~ 3,
                                p12_67_employment_student == 'Yes' ~ 4,
                                p12_67_employment_unemployed == 'Yes' ~ 5,
                                p12_67_employment_retired == 'Yes' ~ 6),
         ethnicitycchs = case_when(ethnicityall1 == 1 & morethan1ethc == 0 ~ "4. Black",
                                   ethnicityall2 == 1 & morethan1ethc == 0  ~ "7, 9 Arab / West Asian",
                                   ethnicityall3 == 1 & morethan1ethc == 0 ~ "3, 10, 11. East Asian",
                                   ethnicityall4 == 1 & morethan1ethc == 0 ~ "4. Black",
                                   ethnicityall5 == 1 & morethan1ethc == 0 ~ "4. Black",
                                   ethnicityall6 == 1 & morethan1ethc == 0 ~ "Indigenous",
                                   ethnicityall7 == 1 & morethan1ethc == 0 ~ "6. Latin",
                                   ethnicityall8 == 1 & morethan1ethc == 0 ~ "2. South Asian",
                                   ethnicityall9 == 1 & morethan1ethc == 0 ~ "5. South East Asian",
                                   ethnicityall10 == 1 & morethan1ethc == 0 ~ "1. White",
                                   ethnicityall11 == 1 & morethan1ethc == 0 ~ "12. Other",
                                   morethan1ethc == 1 ~ "Multiracial"),
         race = ifelse(ethnicitycchs == "1. White",0,1),
         province = case_when(p2_14_eligible_province == "British Columbia" ~ "WestCoast",
                              p2_14_eligible_province == "Alberta" ~ "Prairies",
                              p2_14_eligible_province == "Manitoba" ~ "Prairies",
                              p2_14_eligible_province == "Saskatchewan" ~ "Prairies",
                              p2_14_eligible_province == "Ontario" ~ "Central",
                              p2_14_eligible_province == "Quebec" ~ "Central",
                              p2_14_eligible_province == "New Brunswick" ~ "Atlantic",
                              p2_14_eligible_province == "Prince Edward Island" ~ "Atlantic",
                              p2_14_eligible_province == "Nova Scotia" ~ "Atlantic",
                              p2_14_eligible_province == "Newfoundland & Labrador" ~ "Atlantic",
                              p2_14_eligible_province == "Northwest Territories" ~ "Northern",
                              p2_14_eligible_province == "Nunavut" ~ "Northern",
                              p2_14_eligible_province == "Yukon" ~ "Northern"),
         STI_test12m = case_when(p15_75_last_STI_test %in% c("4-6 months ago","7-12 months ago","In the past 3 months") ~ 1,
                                 p15_75_last_STI_test %in% c("Longer than a year ago","Never") ~ 0),
         HIV_test12m = case_when(p31_114_HIV_last_test_date %in% c("4-6 months ago","7-12 months ago","In the past 3 months") ~ 1,
                                 p31_114_HIV_last_test_date %in% c("Longer than a year ago","Never","I have never tested for") ~ 0),
         lastcondom_use = case_when(p47_199_last_male_anal_sex_condom == "Yes, I used a condom" ~ 1,
                                    p36_137_HIRI_bottom_no_condom_P == "Yes (score of 10)" ~ 1,
                                    p36_139_HIRI_top_no_condom_P6M == "5+ times (score of 6)" ~ 1,
                                    p47_199_last_male_anal_sex_condom %in% c("8888: Has not had any male sex partners in past six months","8888: Never had sex with a man") ~ 2,
                                    p47_199_last_male_anal_sex_no_condom == "No, I did not use a condom" ~ 0,
                                    p47_203_prevention_condoms == "Yes" ~ 1),
         illicit_drug = case_when(cocaine == 1 | crystal == 1 | ecstasy == 1 | heroin == 1 | steroid == 1 | psychedelic == 1 | ketamine == 1 | crack == 1 | marijuana == 1 ~ 1,
                                   cocaine == 0 & crystal == 0 & ecstasy == 0 & heroin == 0 & steroid == 0 & psychedelic == 0 & ketamine == 0 & crack == 0 & marijuana == 0 ~ 0),
         survey = 1, year = 2019, msm = 1, WTS_M = 1, WTS_M_rescaled = 1, ADM_RNO= floor(runif(n=nrow(sn2019_clean),min=10000,max=99999))) %>%
  rename(srmh = mental_health_general,
         depress_binary = depression_binary) %>%
  select(c(education,income,orient,age_grp,msm,grh,srmh,illicit_drug,have_pcp,insurance,consult_mh,marijuana,smoke,employ,province,rural,ADM_RNO,
           race,depress_binary,existing_mh,immigration,alcohol,HIV_test12m,STI_test12m,lastcondom_use,survey,disclosure,year,WTS_M,WTS_M_rescaled))

cchsSNall<-rbind(cchsall,sn2019_clean)
table(sn2019_clean$orient)

save(cchsSNall,file="cchsSNall.Rdata")
