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
library(ggplot2)
require(stabs)
require(lars)
require(glmnet)



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
demo_info = demo_info[c("id","sex", "momid", "dadid","famid_strata")]
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


data_full = data.frame(demo_info_select, brain_data, cov_data, snp_data)



################################################
## model building
################################################
DV_list = names(brain_data) #c("FO_norm", "DV_norm")
NT_list = c("ACH", "DA", "NE", "5HT")


snps_formulas = list()
for (DV_norm in DV_list) {
  for (NT in NT_list){
    snp_data_temp = get(paste0('snp_data_', NT)) # Load relevant SNP dataset
    
    #---------------------------------------------------------
    # remove VariantIDs (SNPs) with a single factor. 
    # A single factor predictor will not allow building a model
    #---------------------------------------------------------
    checkdata = data.frame(sapply(lapply(snp_data_temp, unique), length))
    checkdata = cbind(rownames(checkdata), data.frame(checkdata, row.names=NULL))
    colnames(checkdata) = c("VariantID","FactorNum")
    
    checkdata_idx1 = checkdata[grep("1", checkdata$FactorNum),]
    checkdata_idx1 = as.list(checkdata_idx1$VariantID)
    
    if (length(checkdata_idx1)==0){    } else {
      snp_data_temp = snp_data_temp[, -which(names(snp_data_temp) %in% checkdata_idx)]}
    
    snps_formulas[paste0(DV_norm, "_", NT)] = paste(names(snp_data_temp), collapse = " + ")
    
    rm(list=c("snp_data_temp", "checkdata", "checkdata_idx1"))
  }
}


# Initialize an empty list for model formulas
mdl_formulas = list()
mdl_basic_formula = "age + sex + FD + n_ps1 + n_ps2 + n_ps3 + n_ps4 + n_ps5 + n_ps6 + n_ps7 + n_ps8 + n_ps9 + n_ps10"

# Generate the model formulas
for (model_name in names(snps_formulas)){
  DV_norm = substr(model_name, 1, 7)
  formula_str = paste(DV_norm, "~", mdl_basic_formula, "+", snps_formulas[[model_name]], sep="")
  mdl_formulas[[model_name]] = formula_str
}


################################################
## stability selection
################################################
for (model_name in names(mdl_formulas)){
  DV_norm = substr(model_name, 1, 7)
  NT = substr(model_name, 9, nchar(model_name))
  
  snp_data_temp = get(paste0('snp_data_', NT))
  checkdata = data.frame(sapply(lapply(snp_data_temp, unique), length))
  checkdata = cbind(rownames(checkdata), data.frame(checkdata, row.names=NULL))
  colnames(checkdata) = c("VariantID","FactorNum")
  checkdata_idx1 = checkdata[grep("1", checkdata$FactorNum),]
  checkdata_idx1 = as.list(checkdata_idx1$VariantID)
  
  if (length(checkdata_idx1)==0){    } else {
    snp_data_temp = snp_data_temp[, -which(names(snp_data_temp) %in% checkdata_idx)]}
  
  
  X       = as.matrix(cbind(cov_data, snp_data_temp))
  Y       = brain_data[[DV_norm]]
  data    = data.frame(demo_info_select, brain_data[DV_norm], X) 
  
  
  stab_ss.glmlasso = stabsel(X, Y, 
                          fitfun = glmnet.lasso,
                          args.fitfun = list(type = "anticonservative"),
                          cutoff = 0.90, 
                          PFER = 1, 
                          B = 50,
                          folds = subsample(weights = rep(1, nrow(X)), B = 50, strata = as.factor(data$famid_strata)), #a weight matrix that represents the subsamples.
                          sampling.type = "SS",
                          assumption = "r-concave", 
                          papply = mclapply, 
                          mc.preschedule = FALSE,
                          verbose = TRUE, 
                          eval = TRUE)

 
  # plot(stab.glmlasso)
  
  rm(list = c("DV_norm", "NT", "snp_data_temp", "checkdata", "checkdata_idx1", "X", "Y", "data", "stab_ss.glmlasso")) #, "stab_ss.glmlasso_unimodal"
}
