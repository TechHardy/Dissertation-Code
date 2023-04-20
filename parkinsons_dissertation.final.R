################################################################################################################################################
##
## Luke Pilling -- 2022.11.07

## analysis of baseline phenotypes with dementia and Parkinson's disease in UK Biobank 
## syntax to load and merge data files, and example syntax 

################################################################################################################################################
##
## load useful packages for data loading and manipulation
library(tidyverse)
library(lubridate)
library(vroom)
library(haven)
library(broom) ## useful for getting better summary output from model fit
library(lukesRlib)
library(gt)
library(ggplot2)


################################################################################################################################################
##
## load UK Biobank data -- rather than duplicate lots of files, we merge in what we need as we need it

## load baseline phenotypes (age, sex, assessment center, date of birth, date of attending assessment center, etc.)
ukb = read_dta("I:/Projects/BioBank/Luke/GRS/500K/UKB_baseline_phenos.2022-10-25.dta")

## load mortality data (includes date of death)
ukb_t = read_dta("I:/Projects/BioBank/14631_ageing-well/Death data/ukb14631_death_20230412.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")

## load dementia diagnosis
##    when merging just keep three phenotypes for now -- in future can use more -- variables described below
ukb_t = vroom("I:/Projects/BioBank/Luke/GWAS/phenos/dementia/ukb14631.dementia_hes2021gp2017.20220721.txt")
ukb = left_join(ukb, ukb_t %>% select(n_eid_14631,dementia_hes2021gp2017,dementia_hes2021gp2017_df,dementia_ad_clean_hes2021gp2017,dementia_ad_clean_hes2021gp2017_df,dementia_vasc_clean_hes2021gp2017,dementia_vasc_clean_hes2021gp2017_df), by="n_eid_14631")

## load Parkinson's disease diagnosis (includes binary yes/no variables, as well as date)
##    when merging just keep two phenotypes for now -- in future can use more
ukb_t = vroom("I:/Projects/BioBank/Luke/GWAS/phenos/movement_disorders/ukb14631.5outcomes_hes2021gp2017.20220707.txt")
ukb = left_join(ukb, ukb_t %>% select(n_eid_14631,parkinsons_hes2021gp2017,parkinsons_hes2021gp2017_df,essential_tremor_hes2021gp2017,essential_tremor_hes2021gp2017_df), by="n_eid_14631")

## load APOE genotype
ukb_t = read_dta("I:/Projects/_UKB_bulk/Genetics/SNPs/APOE_e4/APOE.haplotype.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")

## load iron serum biomarker polygenic scores from HUNT paper
ukb_t = read_dta("I:/Projects/_UKB_bulk/Genetics/GRS/v3/SNP_list_211021_ironHUNT/SNP_list_211021_ironHUNT.txt.451K.grs_list.txt.merged.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")

## load genetic covariates
ukb_t = read_dta("I:/Projects/BioBank/14631_ageing-well/gwas_vars/500K/UKB_500K_GWAS_VARS.dta")
ukb = left_join(ukb, ukb_t[,c(1,10,12:56)], by="n_eid_14631")

## load C282Y (rs1800562) genotype
ukb_t = read_dta("I:/Projects/_UKB_bulk/Genetics/UKB_WES_470K/HFE/ukb14631_hfe_composite_genotypes.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")



## load relevant MRI phenos from the GWAS
ukb_t = read_dta("I:/Projects/BioBank/Luke/GWAS/phenos/mri/ukb14631.mri39k.grey_matter_1.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")
ukb_t = read_dta("I:/Projects/BioBank/Luke/GWAS/phenos/mri/ukb14631.mri39k.t2star.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")

## load MRI covariates
ukb_t = read_dta("I:/Projects/BioBank/14631_ageing-well/data downloaded 2021 October - update all fields/ukb48975.assessment_info.dta")
ukb = left_join(ukb, ukb_t, by="n_eid_14631")



## load HES 2021
load("I:/Projects/BioBank/14631_ageing-well/HES up to 2021 Sept/ukb14631_HES_20220401.RDat")
ukb = left_join(ukb, ukb_hes2021, by="n_eid_14631")
ukb_t = NULL



################################################################################################################################################

## Hardy's part

# Filtered NA UKB for TSAT holes and c282y

natsatukb <- drop_na(ukb, tsat21hunt_grs_20)

natsatmale <- natsatukb %>% select(sex) %>% sum()

nac282yukb <- drop_na(ukb, rs1800562_g_a)

nac282ymale <- nac282yukb %>% select(sex) %>% sum()



# Counts the number of male participants that have incident Parkinson's

maleswpar <- ukb %>% filter(sex == 1 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) %>% count(parkinsons_hes2021gp2017==1)

# Counts the number of female participants that have incident Parkinson's

femaleswpar <- ukb %>% filter(sex == 0 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) %>% count(parkinsons_hes2021gp2017==1)



# Counts the number of male participants that have c282y

maleswhfe <- ukb %>% filter(rs1800562_g_a == 2 & sex == 1) %>% count()

# Counts the number of male participants that have no hfe

maleswnohfe <- nac282yukb %>% filter(rs1800562_g_a != 2 & sex == 1) %>% count()

# Counts the number of male participants that have incident Parkinson's and no hfe

maleswparnohfe <- nac282yukb %>% filter(rs1800562_g_a != 2 & sex == 1 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) %>% count()

# Counts the number of male participants that have incident Parkinson's and hfe

maleswparhfe <- ukb %>% filter(rs1800562_g_a == 2 & sex == 1 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) %>% count()



# Calculates the mean age of participants with incident parkinsons by sex

agefpar <- ukb %>% filter(sex == 0 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date) %>% select(age) %>% lapply(mean, na.rm = TRUE)

agempar <- ukb %>% filter(sex == 1 & parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date) %>% select(age) %>% lapply(mean, na.rm = TRUE)

agepar <- ukb %>% filter(parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date) %>% select(age) %>% lapply(mean, na.rm = TRUE)



# List of characteristics to characterise:
# Smoking status ratio of degenerated

# Association between sex

sexpar_lm <- glm(parkinsons_hes2021gp2017 ~ age   + sex, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)))

# Association between BMI

bmipar_lm <- glm(parkinsons_hes2021gp2017 ~ age + sex + bmi, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)))

# Association between age

agepar_lm <- glm(parkinsons_hes2021gp2017 ~ bmi + sex + age, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)))

# Association between Parkinsons and dementia - not used

pardem_log <- glm(parkinsons_hes2021gp2017 ~ age + sex + dementia_hes2021gp2017, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))



# Association of c282y with liver disorders

hfecirr <- glm(fibr_cirrhosis_hes2021v2 ~ age + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

hfelivdis <- glm(liver_disease_hes2021v2 ~ age + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

hfet2d <- glm(t2d_hes2021v2 ~ age + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

hfechd <- glm(chd_hes2021v2 ~ age + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))


hfecirrmen <- glm(fibr_cirrhosis_hes2021v2 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)  & sex==1), family = binomial(link="logit"))

hfelivdismen <- glm(liver_disease_hes2021v2 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)  & sex==1), family = binomial(link="logit"))

hfet2dmen <- glm(t2d_hes2021v2 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)  & sex==1), family = binomial(link="logit"))

hfechdmen <- parchd <- glm(chd_hes2021v2 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)  & sex==1), family = binomial(link="logit"))



# Association with coronary heart disease, type 2 diabetes, liver diseases and liver fibrosis

parlivfib <- glm(parkinsons_hes2021gp2017 ~ age + sex + fibr_cirrhosis_hes2021v2, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

parlivdis <- glm(parkinsons_hes2021gp2017 ~ age + sex + liver_disease_hes2021v2, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

part2d <- glm(parkinsons_hes2021gp2017 ~ age + sex + t2d_hes2021v2, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))

parchd <- glm(parkinsons_hes2021gp2017 ~ age + sex + chd_hes2021v2, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family = binomial(link="logit"))



## Analysis of Transferrin Saturation polygenic score (tsat21)

# Association between Parkinsons and transferrin saturation

partrans <- glm(parkinsons_hes2021gp2017 ~ age + sex + tsat21hunt_grs_20 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))

partransm <- glm(parkinsons_hes2021gp2017 ~ age + sex + tsat21hunt_grs_20 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(sex == 1 & (parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date))), family=binomial(link="logit"))

partransf <- glm(parkinsons_hes2021gp2017 ~ age + sex + tsat21hunt_grs_20 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(sex == 0 & (parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date))), family=binomial(link="logit"))

# Association between Parkinsons and transferrin saturation nohfe

partransnohfe <- glm(parkinsons_hes2021gp2017 ~ age + sex + tsat21hunt_nohfe_grs_19 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



## Analysis of Ferritin polygenic score (ferritin21)

# Association between Parkinsons and ferritin

parfer <- glm(parkinsons_hes2021gp2017 ~ age + sex   + ferritin21hunt_grs_65 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))

# Association between Parkinsons and ferritin nohfe

parfernohfe <- glm(parkinsons_hes2021gp2017 ~ age + sex   + ferritin21hunt_nohfe_grs_64 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



## Analysis of iron polygenic score (iron21)

# Association between Parkinsons and iron

parir <- glm(parkinsons_hes2021gp2017 ~ age + sex   + iron21hunt_grs_21 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))

# Association between Parkinsons and iron nohfe

parirnohfe <- glm(parkinsons_hes2021gp2017 ~ age + sex   + iron21hunt_nohfe_grs_20 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



## Analysis of iron binding capacity polygenic score (tibc21)

# Association between Parkinsons and iron binding capacity

partib <- glm(parkinsons_hes2021gp2017 ~ age + sex   + tibc21hunt_grs_42 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))

# Association between Parkinsons and iron binding capacity nohfe

partibnohfe <- glm(parkinsons_hes2021gp2017 ~ age + sex   + tibc21hunt_nohfe_grs_41 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



##  Analysis of HFE C282Y

# Parkinson's, males only
parhfem <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parhfef <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parhfe <- glm(parkinsons_hes2021gp2017 ~ age + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



##  Analysis of HFE H63D

# Parkinson's, males only
parh63Dm <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(rs1799945_c_g) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parh63Df <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(rs1799945_c_g) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parh63D <- glm(parkinsons_hes2021gp2017 ~ age + sex + as.factor(rs1799945_c_g) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



##  Analysis of HFE C282Y + H63D

# Parkinson's, males only
parhfehdm <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(c282y_h63d_genotypes) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parhfehdf <- glm(parkinsons_hes2021gp2017 ~ age + as.factor(c282y_h63d_genotypes) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parhfehd <- glm(parkinsons_hes2021gp2017 ~ age + sex + as.factor(c282y_h63d_genotypes) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = ukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



# Association of HFE c282Y controlled for liver disease

# Parkinson's, males only
parhfemlivdev <- glm(parkinsons_hes2021gp2017 ~ age + liver_disease_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parhfeflivdev <- glm(parkinsons_hes2021gp2017 ~ age + liver_disease_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parhfelivdev <- glm(parkinsons_hes2021gp2017 ~ age + liver_disease_hes2021v2 + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



# Association of HFE c282y controlled for liver cirrhosis

# Parkinson's, males only
parhfemcirr <- glm(parkinsons_hes2021gp2017 ~ age + fibr_cirrhosis_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parhfefcirr <- glm(parkinsons_hes2021gp2017 ~ age + fibr_cirrhosis_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parhfecirr <- glm(parkinsons_hes2021gp2017 ~ age + fibr_cirrhosis_hes2021v2 + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



# Association of HFE c282y controlled for t2d

# Parkinson's, males only
parhfemt2d <- glm(parkinsons_hes2021gp2017 ~ age + t2d_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==1), family=binomial(link="logit"))

# Parkinson's, females only
parhfeft2d <- glm(parkinsons_hes2021gp2017 ~ age + t2d_hes2021v2 + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date) & sex==0), family=binomial(link="logit"))

# Parkinson's
parhfet2d <- glm(parkinsons_hes2021gp2017 ~ age + t2d_hes2021v2 + sex + as.factor(rs1800562_g_a) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = nac282yukb %>% filter(parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date)), family=binomial(link="logit"))



# Association of pd tsat when c282y == 2

partransc282y <- glm(parkinsons_hes2021gp2017 ~ tsat21hunt_grs_20 + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(rs1800562_g_a == 2 & (parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date))), family=binomial(link="logit"))

partransc282ym <- glm(parkinsons_hes2021gp2017 ~ tsat21hunt_grs_20 + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(sex == 1 & rs1800562_g_a == 2 & (parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date))), family=binomial(link="logit"))

partransc282yf <- glm(parkinsons_hes2021gp2017 ~ tsat21hunt_grs_20 + age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = natsatukb %>% filter(sex == 0 & rs1800562_g_a == 2 & (parkinsons_hes2021gp2017 == 0 | (parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df > assessment_date))), family=binomial(link="logit"))



##Graphing

#Age histogram

png("ageplot.png")
ageplot <- ggplot(ukb, aes(x=age)) + geom_histogram(color="black", fill="white") + geom_vline(aes(xintercept=mean(age)),color="blue", linetype="dashed", linewidth=1)
ageplot
dev.off()

#Sex bar chart

png("sexplot.png")
malecount <- sum(ukb$sex == 1)
femalecount <- sum(ukb$sex == 0)
tempdf <- data.frame(sex=c("Male", "Female"), number=c(malecount, femalecount))
sexplot <- ggplot(tempdf, aes(x=sex,y=number)) +  geom_bar(stat="identity") + scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n=28))
sexplot
dev.off()

#Tsat Dist histogram/density plot

png("tsatdist.png")
tsatplot <- ggplot(natsatukb, aes(x=tsat21hunt_grs_20)) + geom_histogram(color="black", fill="white") + geom_vline(aes(xintercept=mean(tsat21hunt_grs_20)),color="blue", linetype="dashed", linewidth=1)
tsatplot
dev.off()

#Age with parkinsons histogram

png("paragedist.png")
ageparplot <- ggplot(ukb %>% filter(parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date), aes(x=age)) + geom_histogram(color="black", fill="white") + geom_vline(aes(xintercept=mean(age)),color="blue", linetype="dashed", linewidth=1)
ageparplot
dev.off()

#Tsat dist in c282y

png("tsatc282ydist.png")
tsatplotc282y <- ggplot(natsatukb %>% filter(rs1800562_g_a == 2), aes(x=tsat21hunt_grs_20)) + geom_histogram(color="black", fill="white") + geom_vline(aes(xintercept=mean(tsat21hunt_grs_20)),color="blue", linetype="dashed", linewidth=1)
tsatplotc282y
dev.off()

#Percentage of male normal pop with par vs percentage of c282y with par bar chart

png("EffectSizec282y.png")
normalpar <- nac282yukb %>% filter(parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date & sex == 1 & rs1800562_g_a != 2) %>% count()
normalpop <- nac282yukb %>% filter(sex == 1 & rs1800562_g_a != 2) %>% count()
normalper <- (normalpar/normalpop) * 100
c282ypar <- nac282yukb %>% filter(parkinsons_hes2021gp2017 == 1 & parkinsons_hes2021gp2017_df>assessment_date & sex == 1 & rs1800562_g_a == 2) %>% count()
c282ypop <- nac282yukb %>% filter(sex == 1 & rs1800562_g_a == 2) %>% count()
c282yper <- (c282ypar/c282ypop) * 100
tempdf <- data.frame(Phenotype=c("Normal male", "Homozygous p.c282y male"), Percentage=unlist(c(normalper, c282yper)))
c282yparplot <- ggplot(tempdf, aes(x=Phenotype,y=Percentage)) +  geom_bar(stat="identity") + scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n=20))
c282yparplot
dev.off()
