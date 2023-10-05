rm(list=ls())


################################################
## path & library setting
################################################
dir_working="/your/working/directory"
dir_genotyping="/where/you/saved/your/genetics/information"
dir_output="/where/you/will/save/your/output"

if (file.exists(dir_output)) {
} else {
  dir.create(dir_output)
}

setwd(dir_working)


library(tibble)
library(stringr)
library(dplyr)
library(tidyverse)
library(readxl)
library(coxme)
library(nlme)
library(kinship2)
library(car)
library(glmnet)


################################################
## data setting
################################################
snp_info = read_excel(paste0(dir_genotyping, "/your_snp_data_file.xlsx"), sheet = "pruned_snps_info") 
snp_info = snp_info[!duplicated(snp_info$`Variant ID`), ]

snp_data = read_excel(paste0(dir_genotyping, "/your_snp_data_file.xlsx"), sheet = "pruned_snps_data") 
snp_data = snp_data[, -c(1:6)]

brain_cov_data = read_excel(paste0(dir_genotyping, "/your_snp_data_file.xlsx"), sheet = "setting") 
brain_data = brain_cov_data[,c(5:6)]
cov_data = brain_cov_data[, c(2:4, 11:20)]


#-------------------------------------------------
# divide snp_info by neurotransmitter system
#-------------------------------------------------
snp_ACH = snp_info[grep("acetylcholine", snp_info$`NT system`), ]
snp_NE = snp_info[grep("epinephrine", snp_info$`NT system`), ]
snp_5HT = snp_info[grep("serotonin", snp_info$`NT system`), ]
snp_DA = snp_info[grep("dopamine", snp_info$`NT system`), ]

library(tibble)
snp_ACH = add_column(snp_ACH, data.frame(paste0(snp_ACH$`Variant ID`, '_', snp_ACH$Minor)), .after = "Variant ID")
snp_NE = add_column(snp_NE, data.frame(paste0(snp_NE$`Variant ID`, '_', snp_NE$Minor)), .after = "Variant ID")
snp_5HT = add_column(snp_5HT, data.frame(paste0(snp_5HT$`Variant ID`, '_', snp_5HT$Minor)), .after = "Variant ID")
snp_DA = add_column(snp_DA, data.frame(paste0(snp_DA$`Variant ID`, '_', snp_DA$Minor)), .after = "Variant ID")

snp_ACH = snp_ACH %>% rename(VariantID = 3)
snp_NE = snp_NE %>% rename(VariantID = 3)
snp_5HT = snp_5HT %>% rename(VariantID = 3)
snp_DA = snp_DA %>% rename(VariantID = 3)


#-------------------------------------------------
# divide snp_data by neurotransmitter system
#-------------------------------------------------
snp_data_ACH = snp_data[, intersect(names(snp_data),snp_ACH$VariantID)]
snp_data_NE = snp_data[, intersect(names(snp_data),snp_NE$VariantID)]
snp_data_5HT = snp_data[, intersect(colnames(snp_data),snp_5HT$VariantID)]
snp_data_DA = snp_data[, intersect(names(snp_data),snp_DA$VariantID)]



#-------------------------------------------------
# Kinship information
#-------------------------------------------------
demo_info = read_excel(paste(dir_genotyping, 'your_demographic_file.xlsx', sep = "/"), sheet = "info")
demo_info = demo_info[c("id","sex", "momid", "dadid","famid")]
demo_info_select = demo_info[1:674,]
demo_info_fixed = with(demo_info_select,fixParents(id = id, dadid = dadid, momid = momid, sex = sex))

gped = with(demo_info, pedigree(id = id, dadid = dadid, momid = momid, sex = sex)) 
kmat = kinship(gped) 

reltwins = read_excel(paste(dir_genotyping, 'your_demographic_file.xlsx', sep = "/"), sheet = "reltwins")
reltwins_2 = reltwins[c("id1", "id2", "code")]
reltwins_2 = as.matrix(reltwins_2)
colnames(reltwins_2)=c("id1","id2","code")

pedAll = with(demo_info, pedigree(id=id, dadid=dadid, momid=momid, sex=sex, relation=reltwins_2))
kmat2 = kinship(pedAll)


data = data.frame(demo_info_select, brain_data, cov_data, snp_data)



################################################
## model building 
## (based on stability selection results)
################################################
DV_list = names(brain_data) #c("FO_norm", "DV_norm")
NT_list = c("ACH", "DA", "NE", "5HT")


select_covs_list = list()

select_covs_list[["TP_norm_5HT"]] = c("1", "(1|id)", "n_ps3") #, "n_ps10"
select_covs_list[["FO_norm_ACH"]] = c("1", "(1|id)", "sex") #, "age"
select_covs_list[["FO_norm_DA"]] = c("1", "(1|id)", "sex") #, "age"


select_snps_list = list()

select_snps_list[["TP_norm_5HT"]] = c("rs2770295_C", "rs2070036_T") #, "rs10881831_G", "rs56015034_T"
select_snps_list[["FO_norm_ACH"]] = c("rs1803274_C") #, "rs6067_A"
select_snps_list[["FO_norm_DA"]] = c("rs4646318_G") #, "rs6356_C", "rs56087035_T", "rs73856267_C"

snps_counts = list()

for (model_name in names(select_snps_list)) {
  select_snps_names = select_snps_list[[model_name]]
  snps_counts[[model_name]] = sum(grepl("^rs", select_snps_names))
  }



mdl_full_formulas = list()

for (model_name in names(select_snps_list)){
  DV_norm = gsub(pattern = paste(NT_list, collapse = "|"), "", model_name) # Remove the patterns from NT_list in model_name using gsub
  DV_norm = sub("_$", "", DV_norm) # Remove any trailing underscore at the end of DV_norm

  # full model formula
  formula_str = paste(DV_norm, " ~ ", paste(select_covs_list[[model_name]], collapse = " + "), " + ", paste(select_snps_list[[model_name]], collapse = " + "), sep = "")
  mdl_full_formulas[[model_name]] = as.formula(formula_str)
  
  # full model fitting
  mdl = lmekin(mdl_full_formulas[[model_name]], data = data, varlist = kmat * 2)
  assign(paste0("mdl_full_", model_name), mdl)
  
  rm(list = c("DV_norm", "formula_str", "mdl"))
}


mdl_cov_formulas = list()

for (model_name in names(select_covs_list)){
  DV_norm = gsub(pattern = paste(NT_list, collapse = "|"), "", model_name) # Remove the patterns from NT_list in model_name using gsub
  DV_norm = sub("_$", "", DV_norm) # Remove any trailing underscore at the end of DV_norm

  # base model formula
  formula_str = paste(DV_norm, " ~ ", paste(select_covs_list[[model_name]], collapse = " + "), sep = "")
  mdl_cov_formulas[[model_name]] = as.formula(formula_str)
  
  # base model fitting
  mdl = lmekin(mdl_cov_formulas[[model_name]], data = data, varlist = kmat * 2)
  assign(paste0("mdl_cov_", model_name), mdl)
  
  rm(list = c("DV_norm", "formula_str", "mdl"))
}


################################################
## model testing
################################################
chsq = list()

for (model_name in names(select_snps_list)){
  DV_norm = gsub(pattern = paste(NT_list, collapse = "|"), "", model_name) # Remove the patterns from NT_list in model_name using gsub
  DV_norm = sub("_$", "", DV_norm) # Remove any trailing underscore at the end of DV_norm
  
  chsq[[model_name]] = 2 * (get(paste0("mdl_full_", model_name))$loglik - get(paste0("mdl_cov_", model_name))$loglik)
  
  rm(DV_norm)
}


pval = list()

for (model_name in names(select_snps_list)){
  DV_norm = gsub(pattern = paste(NT_list, collapse = "|"), "", model_name) # Remove the patterns from NT_list in model_name using gsub
  DV_norm = sub("_$", "", DV_norm) # Remove any trailing underscore at the end of DV_norm
  
  pval[[model_name]] = pchisq(chsq[[model_name]], df = snps_counts[[model_name]], lower.tail = FALSE)
  
  rm(DV_norm)
}








