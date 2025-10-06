#Data analysis for PVOD project
library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(data.table)
library(ggplotify)
setwd('C:/mywd')

#load in the data
age_adjusted_gdf_v2 <- read.csv2("~/myfilepath/age_adjusted_gdf_v2.csv")
samples_for_Akiko_2nd_shipment <- read_excel("samples_for_Akiko_2nd_shipment.xlsx")
Updated_list_akiko <- read_excel("Updated_list_akiko.xlsx")

v1 <- age_adjusted_gdf_v2 %>% select(cohort_id, diagnosis, Age)
v2 <- samples_for_Akiko_2nd_shipment %>% select(SubjectID, diagnosis_verified, `age at sampling`)
names(v2) <- c('cohort_id', 'diagnosis', 'Age')

ages <- rbind(v1,v2)
ages <- ages %>% group_by(cohort_id) %>% filter(Age == max(Age))

names(ages) <- c('cohort_id', 'diagnosis', 'age_at_sampling')

list_akiko <- Updated_list_akiko %>% select(`Patient ID`, `PH Group`, `GDF15 (pg/ml)`, `Ext Subj`, Diagnosis)

#fix weird values in columns
list_akiko$`Ext Subj` <- ifelse(is.na(list_akiko$`Ext Subj`), list_akiko$`Patient ID`, list_akiko$`Ext Subj`)

#list_akiko <- list_akiko %>% filter(str_starts(`Ext Subj`, "PAH"))
list_akiko2 <- merge(ages, list_akiko, by.x='cohort_id', by.y='Ext Subj')
list_akiko2 <- unique(list_akiko2)

#find duplicates
duplicates <- list_akiko2 %>% group_by(cohort_id) %>% filter(n() > 1) %>% ungroup()
#no duplicates

missed <- list_akiko %>% filter(!`Ext Subj` %in% list_akiko2$cohort_id)
#3 people missed
add <- Updated_list_akiko %>% filter(`Patient ID` %in% missed$`Ext Subj`)
add <- add[,2:9]

list_akiko2 <- list_akiko2 %>% select(cohort_id, age_at_sampling, diagnosis, `GDF15 (pg/ml)`)
add <- add %>% select(`Patient ID`, Age, Diagnosis, `GDF15 (pg/ml)`)
names(add) <- c('cohort_id', 'age_at_sampling', 'diagnosis', 'GDF15 (pg/ml)')

list_definitive <- rbind(list_akiko2, add)

#load in clinical data
all_oc_ids <- readRDS("~/myfilepath/all_oc_ids.rds")
clin_dat_at_diagnosis_clean <- readRDS("~/myfilepath/clin_dat_at_diagnosis_clean.rds")
mutations_PAH_cohort_easy <- readRDS("~/myfilepath/mutations_PAH_cohort_easy.rds")

ids <- all_oc_ids %>% select(id, id_cohort, id_wgs)
ids <- merge(ids, mutations_PAH_cohort_easy[,c(1,5)], by.x='id_wgs', by.y='WGS.ID', all=T)
ids <- ids %>% filter(!is.na(id))
ids$BMPR2_mutation <- ifelse(is.na(ids$id_wgs), 'no sequencing', ids$BMPR2_mutation)
ids$BMPR2_mutation <- ifelse(is.na(ids$BMPR2_mutation), 'no mutation', ids$BMPR2_mutation)

#get clinical data ready
clin <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex, hb_pvr_calc, cbt_card_ntprobnp_ngpl, hb_pap_m, hb_cardiac_index_value_1, hb_cardiac_output_value_1, ep_1_distance_meters, ep_1_type_6mwt)
clin$MWD6 <- ifelse(clin$ep_1_type_6mwt == 'corridor', clin$ep_1_distance_meters, NA)
clin$ISWD <- ifelse(clin$ep_1_type_6mwt == 'shuttle', clin$ep_1_distance_meters, NA)
clin <- clin %>% select(-ep_1_type_6mwt, -ep_1_distance_meters)

clin <- merge(clin, ids, by='id')

list_definitive <- merge(list_definitive, clin, by.x='cohort_id', by.y='id_cohort', all=T)
list_definitive <- list_definitive %>% filter(!is.na(`GDF15 (pg/ml)`))

bl_df <- Updated_list_akiko %>% select(`Patient ID`, `Ext Subj`, Age, Diagnosis, Sex, `GDF15 (pg/ml)`, PVR, `mPAP (mmHg)`, `Cardiac index`, `NT-proBNP (ng/l)`, `6MWD distance`, `ISWD distance`, `BMPR2 mutation`)
bl_df$`Ext Subj` <- ifelse(is.na(bl_df$`Ext Subj`), bl_df$`Patient ID`, bl_df$`Ext Subj`)

ucsf_patients <- list_definitive %>% filter(is.na(id))
bl_df <- bl_df %>% filter(`Ext Subj` %in% ucsf_patients$cohort_id)

names(bl_df) <- c('id', 'cohort_id', 'age_at_sampling', 'diagnosis', 'sex', 'GDF15 (pg/ml)', 'hb_pvr_calc', 'hb_pap_m', 'hb_cardiac_index_value_1', 'cbt_card_ntprobnp_ngpl', 'MWD6', 'ISWD', 'BMPR2_mutation', 'hb_cardiac_output_value_1', 'id_wgs', 'age_diagnosis')
bl_df$hb_cardiac_output_value_1 <- NA
bl_df$age_diagnosis <- NA
bl_df$id_wgs <- NA

list_definitive2 <- rbind(list_definitive %>% filter(!is.na(id)), bl_df)
#now find duplicate

duplicates <- list_definitive2 %>% group_by(cohort_id) %>% filter(n() > 1) %>% ungroup()

list_definitive2 <- unique(list_definitive2)
#now works! 

#now format the data 
list_definitive2$sex <- ifelse(list_definitive2$sex %in% c('female', 'F'), 'female', 'male')

#some PVRs <10, this must be in WU
list_definitive2$hb_pvr_calc <- ifelse(list_definitive2$hb_pvr_calc <10, list_definitive2$hb_pvr_calc*80, list_definitive2$hb_pvr_calc)

#remove cols not relevant
list_definitive2 <- list_definitive2 %>% select(-id_wgs)
df_akiko <- list_definitive2
write_rds(df_akiko, file='GDF_file_Akiko_July25.rds')

#====================================================================================
#residualise the data
hist(df_akiko$`GDF15 (pg/ml)`)
hist(df_akiko$age_at_sampling)
#both not normal, although age is a bit more normal

df_akiko$log_age <- log(df_akiko$age_at_sampling +1)
df_akiko$log_GDF <- log(df_akiko$`GDF15 (pg/ml)`)

#do LOESS regression, as histograms still not normal!
library(mgcv)

# GAM: flexible, nonparametric model
rownames(df_akiko) <- df_akiko$cohort_id
model <- gam(log_GDF ~ s(log_age), data = df_akiko)
df_akiko$age_adjusted_GDF15 <- residuals(model)

#also run survival analysis
df_akiko$time_from_diagnosis_to_sampling <- df_akiko$age_at_sampling - df_akiko$age_diagnosis

#================================================================================================
#run CoxPH
clin_dat_at_diagnosis_clean <- readRDS("~/myfilepath/clin_dat_at_diagnosis_clean.rds")

#prepare mortality df
mort_df <-  clin_dat_at_diagnosis_clean %>% select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'sub_cause', 'sub_date', 'diagnosis_date', 'hb_pvr_calc', 'cbt_haem_wbc_x10e9pl', 'cbt_haem_hb_gpl', 'hb_rap_m', 'ec_tricuspid_apse', 'cbt_inflammation_crp_mgpl', 'egfr_mdrd')
mort_df <- unique(mort_df)
#make sure sub_date is complete and format the data
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$sub_date2 <- fifelse(is.na(mort_df$sub_date), ymd("2022-07-01"), mort_df$sub_date)
mort_df$surv_time <- (mort_df$sub_date2 - mort_df$diagnosis_date)
mort_df[c('surv_time', 'day')] <- str_split_fixed(mort_df$surv_time, ' ', 2)
mort_df$surv_time <- as.numeric(mort_df$surv_time)
mort_df$surv_time <- mort_df$surv_time/365.2422
mort_df$event <- ifelse(mort_df$sub_cause == 'death', 1, 0)
mort_df$event <- ifelse(is.na(mort_df$event), 0, mort_df$event)
mort_df$diagnosis_verified <- as.factor(mort_df$diagnosis_verified)
mort_df$sex <- as.factor(mort_df$sex)
mort_df <- mort_df %>% filter(!is.na(surv_time))
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)

#take survival up to 25 years
mort_df$event <- ifelse(mort_df$surv_time >=30, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >=30, 30, mort_df$surv_time)

mdf <- merge(mort_df[,c(1,17,19)], df_akiko, by='id')
#149 retained
mdf$surv_time_from_samping <- mdf$surv_time - mdf$time_from_diagnosis_to_sampling
#1 tiny error due to rounding error (last blood from last census, so negative adjusted surv time, make this just over 0)
mdf$surv_time_from_samping <- ifelse(mdf$surv_time_from_samping <0, 0.01, mdf$surv_time_from_samping)

mdf$GDF15 <- mdf$`GDF15 (pg/ml)`
mdf$GDF15_scale <- scale(mdf$GDF15)

#only have group 1 PAH & adult patients
mdf <- mdf %>% filter(diagnosis %in% c('IPAH', 'HPAH', 'PVOD') & age_diagnosis >=16)

#add UCSF patients with data
mdf$id<- as.character(mdf$id)
mdf[146,c(1,2,3,5,6,7, 9)] <- c('UCSF1', #patient data entered here (but not shown due to privacy concerns)
mdf[147,c(1,2,3,5,6,7, 9)] <- c('UCSF2', #patient data entered here (but not shown due to privacy concerns)
mdf[148,c(1,2,3,5,6,7, 9)] <- c('UCSF3', #patient data entered here (but not shown due to privacy concerns)

mdf[149,c(1,2,3,5,6,7, 9)] <- c('UCSF4', #patient data entered here (but not shown due to privacy concerns)
mdf[150,c(1,2,3,5,6,7, 9)] <- c('UCSF5', #patient data entered here (but not shown due to privacy concerns)
mdf[151,c(1,2,3,5,6,7, 9)] <- c('UCSF6', #patient data entered here (but not shown due to privacy concerns)
mdf[152,c(1,2,3,5,6,7, 9)] <- c('UCSF7', #patient data entered here (but not shown due to privacy concerns)

mdf$surv_time<- as.numeric(mdf$surv_time)
mdf$event <- as.numeric(mdf$event)
mdf$age_at_sampling <- as.numeric(mdf$age_at_sampling)

mdf$GDF15 <- mdf$`GDF15 (pg/ml)`
mdf$GDF15 <- as.numeric(mdf$GDF15)
mdf$GDF15_scale <- scale(mdf$GDF15)

sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sobj2 <- Surv(mdf$surv_time_from_samping, mdf$event, type='right')

cox <- coxph(sobj ~ GDF15_scale + sex + age_at_sampling, data=mdf)

ggforest(cox, data = mdf)
grid::grid.text("Scaled GDF15 survival association (from diagnosis)", 
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 10, fontface = "bold"))

#now adjust for disease subtype
cox <- coxph(sobj ~ GDF15_scale + sex + age_at_sampling + diagnosis, data=mdf)

ggforest(cox, data = mdf)
grid::grid.text("Scaled GDF15 survival association (from diagnosis)", 
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 10, fontface = "bold"))

#now construct for survival from age at sampling
cox <- coxph(sobj2 ~ GDF15_scale + sex + age_at_sampling, data=mdf)

ggforest(cox, data = mdf)
grid::grid.text("Scaled GDF15 survival association (from sampling)", 
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 10, fontface = "bold"))

#now adjust for disease subtype
cox <- coxph(sobj2 ~ GDF15_scale + sex + age_at_sampling + diagnosis, data=mdf)

ggforest(cox, data = mdf)
grid::grid.text("Scaled GDF15 survival association (from sampling)", 
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 10, fontface = "bold"))

#add in adjustment for disease duration
#now construct for survival from age at sampling
cox <- coxph(sobj2 ~ GDF15_scale + sex + age_at_sampling + time_from_diagnosis_to_sampling + diagnosis, data=mdf)

ggforest(cox, data = mdf)
grid::grid.text("Scaled GDF15 survival association (from sampling)", 
                x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 10, fontface = "bold"))


#=============================================================================================
#now analyse the longitudinal data
Updated_list_akiko <- read_excel("Updated_list_akiko.xlsx", 
                                 +     sheet = "Longitudinal samples")
longitudinal <- Updated_list_akiko
ages2 <- rbind(v1,v2)
ages2 <- ages2 %>% filter(cohort_id %in% longitudinal$cohort_id)
max_age <- ages2 %>% group_by(cohort_id) %>% filter(Age == max(Age))
max_age$max_age <- max_age$Age
min_age <- ages2 %>% group_by(cohort_id) %>% filter(Age == min(Age))
min_age$min_age <- min_age$Age

akiko_age <- merge(min_age, max_age, by='cohort_id')
longitudinal <- merge(longitudinal, akiko_age[,c(1,4,7)])

longitudinal <- longitudinal %>% pivot_longer(3:4, names_to = 'timepoint', values_to = 'GDF15')
longitudinal$age <- ifelse(longitudinal$timepoint == 'T1', longitudinal$min_age, longitudinal$max_age) 

df_diag <- df_akiko %>% select(cohort_id, age_diagnosis)
longitudinal <- merge(longitudinal, df_diag, by='cohort_id')
longitudinal$time_from_diagnosis <- longitudinal$age - longitudinal$age_diagnosis

hist(longitudinal$GDF15)
hist(longitudinal$time_from_diagnosis)

longitudinal$log_time_from_diagnois <- log(longitudinal$time_from_diagnosis + 1)
longitudinal$log_gdf <- log(longitudinal$GDF15 + 1)

#do histogram
hist(longitudinal$log_gdf)
hist(longitudinal$log_time_from_diagnois)

#not looking great but slightly better

model1 <- lm(log_gdf ~ log_time_from_diagnois + Subtype, data=longitudinal)
summary(model1)

write.csv2(df_akiko, file='baseline_file_Akiko.csv')

#==============================================================
#get the data
GDF_file_Akiko_July25 <- readRDS("~/myfilepath/GDF_file_Akiko_July25.rds")

pvod <- GDF_file_Akiko_July25  %>% filter(diagnosis == 'PVOD') %>% select(id)
#get meds
therapy_type_per_patient <- readRDS("~/myfilepath/therapy_type_per_patient.rds")
pvod_therapy <- merge(pvod, therapy_type_per_patient, by='id')
summary(as.factor(pvod_therapy$medication_group))