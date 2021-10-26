# Age-specific effects of childhood BMI on multiple sclerosis risk: a Mendelian Randomisation study

## Requirements
- [MS GWAS data](https://imsgc.net/?page_id=31)
- [Childhood BMI data](https://www.medrxiv.org/node/346664.external-links.html)
- [Adult BMI GWAS data](https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz)
- [Birth Weight GWAS data](https://egg-consortium.org/BW5/Fetal_BW_European_meta.NG2019.txt.gz)
- TwoSampleMR

## Parts
- MR of childhood BMI per epoch --> MS risk
- Checking pleiotropy, heterogeneity, and MR-PRESSO
- MR of childhood BMI at each time point --> MS risk
- Sensitivity analysis without SNPs influencing adult BMI
- MR of birthweight --> MS (from larger GWAS)


Load the packages we'll need
````R
library(TwoSampleMR)
library(dplyr)
library(readr)
library(ggplot2)
library(readr)
library(gridExtra)
library(Kendall)
library(ggstance)
library(stringr)
````

Load the childhood BMI GWAS and format
````R
setwd("/data/Wolfson-UKBB-Dobson/cbmi_mr_ms")
cbmi = read_csv("top_hits_all_time_points.csv")

# restrict to pval < 5e-8
cbmi = cbmi %>% filter(P_BOLT_LMM<5e-8)
# format

cbmi = cbmi %>% rename("SNP"=snp,"effect_allele"=tested_allele,"other_allele"=reference_allele,"beta"=BETA,"se"=SE,"pval"=P_BOLT_LMM,"eaf"=tested_allele_frequency)

# exclude SE MHC
cbmi = cbmi %>% filter(!(chromosome==6 & position > 25000000 & position < 35000000))
# NB no SNPs excluded

# split into 4 cohorts
# within each cohort order by pval then keep the lowest one (i.e. where one SNP is associated with more than one time point keep the strongest association)
cbmi_birth_6w = cbmi %>% filter(age %in% c("Birth","6w")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE) %>% mutate(Phenotype="cbmi_birth_6w")
cbmi_3mo_1.5y = cbmi %>% filter(age %in% c("3m","6m","8m","1y","1.5y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE) %>% mutate(Phenotype="cbmi_3mo_1.5y")
cbmi_2y_5y = cbmi %>% filter(age %in% c("2y","3y","5y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE) %>% mutate(Phenotype="cbmi_2y_5y")
cbmi_7y_8y = cbmi %>% filter(age %in% c("7y","8y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE) %>% mutate(Phenotype="cbmi_7y_8y")

# format
cbmi_birth_6w = format_data(cbmi_birth_6w,type="exposure")
cbmi_3mo_1.5y = format_data(cbmi_3mo_1.5y,type="exposure")
cbmi_2y_5y = format_data(cbmi_2y_5y,type="exposure")
cbmi_7y_8y = format_data(cbmi_7y_8y,type="exposure")

````

Load the MS GWAS and format
````R
# read in mschip
ms_chip = read_table2("discovery_metav3.0.meta")
ms_chip = mutate(ms_chip, 'variant' = paste0(CHR,':',BP))

# calculate z score from pval
ms_dat = mutate(ms_chip,'z_score' = qnorm(1-(P/2)))
rm(ms_chip)
# beta from or
ms_dat = mutate(ms_dat,'logodds' = log(OR))
# se as beta/z
ms_dat = mutate(ms_dat, 'se' = sqrt((logodds/z_score)^2))
ms_dat = rename(ms_dat,'effect_allele' = A1, 'other_allele' = A2, 'beta' = logodds,'pval'=P)
ms_dat = select(ms_dat,CHR,BP,SNP,variant,effect_allele,other_allele,beta,se,pval)
ms_dat = ms_dat %>% filter(!is.na(beta))
ms_dat = format_data(ms_dat,type="outcome")
````

Harmonise datasets
````R
# write exposure SNPs (pre-harmonisation) to file
exposures = bind_rows(cbmi_birth_6w,
cbmi_3mo_1.5y,
cbmi_2y_5y,
cbmi_7y_8y)

get_f_stat = function(x){
# steiger and moe
mean_f2 = mean(x$beta.exposure^2/x$se.exposure^2)
return(mean_f2)
}
exposure_datasets = list(cbmi_birth_6w,cbmi_3mo_1.5y,cbmi_2y_5y,cbmi_7y_8y)
f2 = lapply(exposure_datasets,get_f_stat)


# harmonise
combo_cbmi_birth_6w = harmonise_data(cbmi_birth_6w,ms_dat)
combo_cbmi_3mo_1.5y = harmonise_data(cbmi_3mo_1.5y,ms_dat)
combo_cbmi_2y_5y = harmonise_data(cbmi_2y_5y,ms_dat)
combo_cbmi_7y_8y = harmonise_data(cbmi_7y_8y,ms_dat)

# clump
combo_cbmi_birth_6w = clump_data(combo_cbmi_birth_6w)
combo_cbmi_3mo_1.5y = clump_data(combo_cbmi_3mo_1.5y)
combo_cbmi_2y_5y = clump_data(combo_cbmi_2y_5y)
combo_cbmi_7y_8y = clump_data(combo_cbmi_7y_8y)


write_csv(combo_cbmi_birth_6w,"combo_cbmi_birth_6w.csv")
write_csv(combo_cbmi_3mo_1.5y,"combo_cbmi_3mo_1.5y.csv")
write_csv(combo_cbmi_2y_5y,"combo_cbmi_2y_5y.csv")
write_csv(combo_cbmi_7y_8y,"combo_cbmi_7y_8y.csv")

# do basic MR
res_epoch1 = mr(combo_cbmi_birth_6w)
res_epoch2 = mr(combo_cbmi_3mo_1.5y)
res_epoch3 = mr(combo_cbmi_2y_5y)
res_epoch4 = mr(combo_cbmi_7y_8y)

# plot results
res_epoch1$epoch = "Birth to 6 weeks"
res_epoch2$epoch = "3 months to 1.5 years"
res_epoch3$epoch = "2 years to 5 years"
res_epoch4$epoch = "7 years to 8 years"

mr_results = bind_rows(res_epoch1,res_epoch2,res_epoch3,res_epoch4)

# write results
write_csv(mr_results,"mr_results_primary_analysis.csv")

# make nice plots
mr_results$epoch = factor(mr_results$epoch,levels=c("Birth to 6 weeks","3 months to 1.5 years","2 years to 5 years","7 years to 8 years"))

p=ggplot(mr_results,aes(b,epoch,col=method))+geom_point()+geom_errorbarh(mapping=aes(xmin=b-1.96*se,xmax=b+1.96*se,height=0.1))+geom_vline(xintercept=0,alpha=0.1)+theme_bw()+facet_wrap(~method,ncol=1)+
labs(x="Beta (MR causal estimate)",y="Childhood BMI epoch",col="MR method")

png("cbmi_mr_fig1.png",res=300,units="in",height=8,width=8)
p
dev.off()

mr_results %>% filter(method=="Inverse variance weighted") %>% mutate(or = exp(b),lowerci = exp(b-1.96*se),upperci=exp(b+1.96*se))

p1=mr_scatter_plot(res_epoch1[3,],dat=combo_cbmi_birth_6w)$v8JIZj.ndQjSx + theme_bw() + labs(x="SNP effect on standardised BMI",y="SNP effect on MS (log OR)")+annotate("label",label="Birth - 6 weeks",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")

p2=mr_scatter_plot(res_epoch2[3,],dat=combo_cbmi_3mo_1.5y)$`76IAZ6.ndQjSx`+ theme_bw() + labs(x="SNP effect on standardised BMI",y="SNP effect on MS (log OR)")+annotate("label",label="3 months - 1.5 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")

p3=mr_scatter_plot(res_epoch3[3,],dat=combo_cbmi_2y_5y)$aSEBgq.ndQjSx+ theme_bw() + labs(x="SNP effect on standardised BMI",y="SNP effect on MS (log OR)")+annotate("label",label="2 years - 5 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")

p4=mr_scatter_plot(res_epoch4[3,],dat=combo_cbmi_7y_8y)$`37Ncgz.ndQjSx`+ theme_bw() + labs(x="SNP effect on standardised BMI",y="SNP effect on MS (log OR)")+annotate("label",label="7 years - 8 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")

png("cbmi_mr_scatterplots.png",res=300,units="in",height=8,width=8)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

# make forest plots
p1=mr_forest_plot(mr_singlesnp(combo_cbmi_birth_6w))$v8JIZj.ndQjSx+theme_bw()+annotate("label",label="Birth - 6 weeks",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")
p2=mr_forest_plot(mr_singlesnp(combo_cbmi_3mo_1.5y))$`76IAZ6.ndQjSx`+theme_bw()+annotate("label",label="3 months - 1.5 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")
p3=mr_forest_plot(mr_singlesnp(combo_cbmi_2y_5y))$aSEBgq.ndQjSx+theme_bw()+annotate("label",label="2 years - 5 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")
p4=mr_forest_plot(mr_singlesnp(combo_cbmi_7y_8y))$`37Ncgz.ndQjSx`+theme_bw()+annotate("label",label="7 years - 8 years",x = -Inf, y = Inf, hjust = 0, vjust = 1)+theme(legend.position="none")

png("cbmi_mr_forestplots.png",res=300,units="in",height=8,width=8)
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

````

Steiger filtering
````R

steiger_fx = function(x_dat) {
  # steiger and moe
  x_dat$samplesize.exposure=28681
  x_dat$samplesize.outcome=(14802+26703)
  x_dat$prevalence.outcome=14802/(14802+26703)
  x_dat$units.outcome="log odds"
  x_dat$units.exposure="z score"
  x_dat$r.outcome =get_r_from_lor(lor=x_dat$beta.outcome,prevalence=x_dat$prevalence.outcome,ncase = 14802,ncontrol=26703,af =x_dat$eaf.exposure)
  x_dat$r.exposure = get_r_from_pn(p=x_dat$pval.exposure,n=x_dat$samplesize.exposure)
  x_dat$rsq.outcome = x_dat$r.outcome^2
  x_dat$rsq.exposure = x_dat$r.exposure^2

  excluded_snps = steiger_filtering(x_dat) %>% filter(steiger_dir==FALSE)
  print(paste0("Number of excluded snps: ",nrow(excluded_snps)))
}

steiger_fx(combo_cbmi_birth_6w)
steiger_fx(combo_cbmi_3mo_1.5y)
steiger_fx(combo_cbmi_2y_5y)
steiger_fx(combo_cbmi_7y_8y)
````
Check pleiotropy & heterogeneity:
````R

pleiotropy=bind_rows(mr_pleiotropy_test(combo_cbmi_birth_6w),
mr_pleiotropy_test(combo_cbmi_3mo_1.5y),
mr_pleiotropy_test(combo_cbmi_2y_5y),
mr_pleiotropy_test(combo_cbmi_7y_8y))
write_csv(pleiotropy,"pleiotropy.csv")

heterogeneity=bind_rows(mr_heterogeneity(combo_cbmi_birth_6w),
mr_heterogeneity(combo_cbmi_3mo_1.5y),
mr_heterogeneity(combo_cbmi_2y_5y),
mr_heterogeneity(combo_cbmi_7y_8y))
write_csv(heterogeneity,"heterogeneity.csv")
````
Run MR-PRESSO:
````R
# mr presso
presso_fx = function(x_dat){
set.seed(1)
presso_res = run_mr_presso(x_dat,NbDistribution=1000)
presso_res
}

presso_fx(combo_cbmi_birth_6w)
presso_fx(combo_cbmi_3mo_1.5y)
presso_fx(combo_cbmi_2y_5y)
presso_fx(combo_cbmi_7y_8y)
````

Leave-one-out
````R
loo_snps = bind_rows(mr_leaveoneout(combo_cbmi_birth_6w),
mr_leaveoneout(combo_cbmi_3mo_1.5y) ,
mr_leaveoneout(combo_cbmi_2y_5y) ,
mr_leaveoneout(combo_cbmi_7y_8y) )

write_csv(loo_snps,"loo.csv")
````

Multivariable MR
````R
# read in Yengo et al GWAS
adult_bmi = read_tsv("Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")

adult_bmi_dat = adult_bmi %>%
rename("effect_allele"=Tested_Allele,
"other_allele"=Other_Allele,
"eaf"=Freq_Tested_Allele_in_HRS,
"beta"=BETA,
"se"=SE,
"pval"=P)

adult_bmi_dat$Phenotype = "a_bmi"
adult_bmi_dat = adult_bmi_dat %>% filter(SNP %in% exposures$SNP)

adult_bmi_dat = format_data(adult_bmi_dat,type="outcome")

# function to do MV MR
epochs = list(cbmi_birth_6w,cbmi_3mo_1.5y,cbmi_2y_5y,cbmi_7y_8y)

mv_mr_results = lapply(epochs,function(x){
message("Filtering ",x$exposure[1])

# combine cbmi and abmi, and harmonise
combo = harmonise_data(x,adult_bmi_dat)

# split harmonised data
abmi =  combo %>% filter(mr_keep==TRUE) %>% select(SNP,contains("outcome"),-samplesize.outcome)
cbmi =  combo %>% filter(mr_keep==TRUE) %>% select(SNP,contains("exposure"))

# change cols
colnames(abmi) = str_replace(colnames(abmi),"outcome","exposure")
abmi$id.exposure = "abmi"
cbmi$id.exposure = "cbmi"

combined_exposure = rbind(cbmi,abmi)
combined_exposure = clump_data(combined_exposure)

if(nrow(combined_exposure)>2){
  # harmonise
  harmonised_combined_exposure = mv_harmonise_data(combined_exposure,ms_dat)

  f_stat_cbmi = mean(harmonised_combined_exposure$exposure_beta[,2]^2/harmonised_combined_exposure$exposure_se[,2]^2)
  f_stat_abmi = mean(harmonised_combined_exposure$exposure_beta[,1]^2/harmonised_combined_exposure$exposure_se[,1]^2)

  message("Mean cbmi F stat for ",x$exposure[1]," is ",f_stat_cbmi)
  message("Mean abmi F stat for ",x$exposure[1]," is ",f_stat_abmi)

  # do MV MR
  mv_res1 = mv_multiple(harmonised_combined_exposure)
  mv_res2 = mv_residual(harmonised_combined_exposure)

  df = data.frame(cbind(harmonised_combined_exposure$exposure_beta,harmonised_combined_exposure$exposure_pval,harmonised_combined_exposure$outcome_beta,harmonised_combined_exposure$outcome_pval))
  df$SNP = rownames(df)
  colnames(df) = c("Beta_adult_bmi","Beta_childhood_bmi","Pval_adult_bmi","Pval_childhood_bmi","Beta_MS","Pval_MS","SNP")
  rownames(df)=NULL
  df = df %>% select(7,1,2,5,3,4,6) %>% mutate(concordant = ifelse((Beta_adult_bmi <0 & Beta_childhood_bmi<0 )|(Beta_adult_bmi>0 & Beta_childhood_bmi>0),"Yes","No")) %>% mutate(GWAS_sig_cbmi = ifelse(Pval_childhood_bmi<5e-8,"Yes","No")) %>% mutate(GWAS_sig_abmi = ifelse(Pval_adult_bmi<5e-8,"Yes","No")) %>% mutate(GWAS_suggestive_abmi = ifelse(Pval_adult_bmi<1e-5,"Yes","No"))

  df$exposure_epoch = as.character(x$exposure[1])
  # combine results
  mv_res1$result$method = "Multiple"
  mv_res2$result$method = "Residual"
  result_df = bind_rows(mv_res1$result,mv_res2$result)
  return(list(df,result_df))
} else {
  message("Insufficient numbers of SNPs. Skipping.")
}
})


overall_res_df = bind_rows(mv_mr_results[[1]][[2]],mv_mr_results[[2]][[2]],mv_mr_results[[4]][[2]]) %>% filter(exposure!="a_bmi") %>% mutate(OR = exp(b)) %>% mutate(lower_ci = exp(b-1.96*se)) %>% mutate(upper_ci = exp(b + 1.96*se))

overall_res_df = overall_res_df %>% select(-id.exposure,-id.outcome,-outcome)

univariable_ivw_results = mr_results %>% filter(method=="Inverse variance weighted") %>% mutate(OR = exp(b),lower_ci = exp(b-1.96*se),upper_ci=exp(b+1.96*se))
univariable_ivw_results$exposure = c("cbmi_birth_6w","cbmi_3mo_1.5y","cbmi_2y_5y","cbmi_7y_8y")
univariable_ivw_results = univariable_ivw_results %>% select(-id.exposure,-id.outcome,-outcome,-epoch)

overall_res_df = overall_res_df %>% bind_rows(univariable_ivw_results)

overall_res_df$epoch = factor(overall_res_df$exposure)
overall_res_df$epoch = recode(overall_res_df$exposure,
  "cbmi_birth_6w"="Birth to 6 weeks",
  "cbmi_3mo_1.5y"="3 months to 1.5 years",
  "cbmi_2y_5y"="2 years to 5 years",
  "cbmi_7y_8y"="7 years to 8 years")

overall_res_df$epoch = factor(overall_res_df$epoch,levels=c("Birth to 6 weeks","3 months to 1.5 years","2 years to 5 years","7 years to 8 years"))

p=ggplot(overall_res_df %>% filter(method %in% c("Inverse variance weighted","Residual")),aes(b,epoch,col=method))+
geom_point(position=position_dodgev(height=0.2))+
ggplot2::geom_errorbarh(mapping=aes(xmin=b-1.96*se,xmax=b+1.96*se,height=0.1),position=position_dodgev(height=0.2))+
geom_vline(xintercept=0,alpha=0.1)+
theme_bw()+
labs(x="Beta (MR causal estimate)",y="Childhood BMI epoch",col="MR method")

png("multivariable_mr.png",res=300,units="in",height=8,width=8)
p
dev.off()

mv_snp_df = bind_rows(mv_mr_results[[1]][[1]],mv_mr_results[[2]][[1]],mv_mr_results[[4]][[1]])
write_csv(mv_snp_df,"mv_snps.csv")

write_csv(overall_res_df,"mv_res_df.csv")

````

Restrict to SNPs not associated with adult BMI
````R
# read in Yengo et al GWAS
adult_bmi = read_tsv("Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt")

adult_bmi_dat = adult_bmi %>%
rename("effect_allele"=Tested_Allele,
"other_allele"=Other_Allele,
"eaf"=Freq_Tested_Allele_in_HRS,
"beta"=BETA,
"se"=SE,
"pval"=P)

adult_bmi_dat = format_data(adult_bmi_dat,type="exposure")

adult_bmi_sigsnps = adult_bmi %>%  filter(P<0.05)

# remove Yengo SNPs associated with adult BMI at p<0.05

abmi_removed_mr_res = lapply(epochs,function(x){
x = x %>% filter(SNP %in% adult_bmi_dat$SNP) %>% filter(!(SNP %in% adult_bmi_sigsnps$SNP))
if(nrow(x)>1){
  combo_x = harmonise_data(x,ms_dat)
  write_csv(path=paste0("abmi_removed_combo_",as.character(x$exposure[1]),".csv"),combo_x)
  res = mr(combo_x)
  print(mr_pleiotropy_test(combo_x))
  res$epoch = as.character(x$exposure[1])
  return(res)
}
})

mr_results = do.call("rbind",abmi_removed_mr_res)
mr_results = mr_results %>% mutate(OR = exp(b)) %>% mutate(lower_ci = exp(b-1.96*se)) %>% mutate(upper_ci = exp(b + 1.96*se))

# write results
write_csv(mr_results,"mr_results_abmi_removed.csv")

````

Repeat at individual time points
````R
cbmi$Phenotype = cbmi$age
cbmi_birth = cbmi %>% filter(age=="Birth")
cbmi_6w = cbmi %>% filter(age=="6w")
cbmi_3m = cbmi %>% filter(age=="3m")
cbmi_6m = cbmi %>% filter(age=="6m")
cbmi_8m = cbmi %>% filter(age=="8m")
cbmi_1y = cbmi %>% filter(age=="1y")
cbmi_1.5y = cbmi %>% filter(age=="1.5y")
cbmi_2y = cbmi %>% filter(age=="2y")
cbmi_3y = cbmi %>% filter(age=="3y")
cbmi_5y = cbmi %>% filter(age=="5y")
cbmi_7y = cbmi %>% filter(age=="7y")
cbmi_8y = cbmi %>% filter(age=="8y")

cbmi_birth = format_data(cbmi_birth,type="exposure")
cbmi_6w = format_data(cbmi_6w,type="exposure")
cbmi_3m = format_data(cbmi_3m,type="exposure")
cbmi_6m = format_data(cbmi_6m,type="exposure")
cbmi_8m = format_data(cbmi_8m,type="exposure")
cbmi_1y = format_data(cbmi_1y,type="exposure")
cbmi_1.5y = format_data(cbmi_1.5y,type="exposure")
cbmi_2y = format_data(cbmi_2y,type="exposure")
cbmi_3y = format_data(cbmi_3y,type="exposure")
cbmi_5y = format_data(cbmi_5y,type="exposure")
cbmi_7y = format_data(cbmi_7y,type="exposure")
cbmi_8y = format_data(cbmi_8y,type="exposure")

# harmonise with ms data
cbmi_birth = harmonise_data(cbmi_birth,ms_dat)
cbmi_6w = harmonise_data(cbmi_6w,ms_dat)
cbmi_3m = harmonise_data(cbmi_3m,ms_dat)
cbmi_6m = harmonise_data(cbmi_6m,ms_dat)
cbmi_8m = harmonise_data(cbmi_8m,ms_dat)
cbmi_1y = harmonise_data(cbmi_1y,ms_dat)
cbmi_1.5y = harmonise_data(cbmi_1.5y,ms_dat)
cbmi_2y = harmonise_data(cbmi_2y,ms_dat)
cbmi_3y = harmonise_data(cbmi_3y,ms_dat)
cbmi_5y = harmonise_data(cbmi_5y,ms_dat)
cbmi_7y = harmonise_data(cbmi_7y,ms_dat)
cbmi_8y = harmonise_data(cbmi_8y,ms_dat)

# clump
cbmi_birth =  clump_data(cbmi_birth )
cbmi_6w =  clump_data(cbmi_6w )
cbmi_3m =  clump_data(cbmi_3m )
cbmi_6m =  clump_data(cbmi_6m )
cbmi_8m =  clump_data(cbmi_8m )
cbmi_1y =  clump_data(cbmi_1y )
cbmi_1.5y =  clump_data(cbmi_1.5y )
cbmi_2y =  clump_data(cbmi_2y )
cbmi_3y =  clump_data(cbmi_3y )
cbmi_5y =  clump_data(cbmi_5y )
cbmi_7y =  clump_data(cbmi_7y )
cbmi_8y =  clump_data(cbmi_8y )

# run mr
res_cbmi_birth = mr(cbmi_birth )
res_cbmi_6w = mr(cbmi_6w )
res_cbmi_3m = mr(cbmi_3m )
res_cbmi_6m = mr(cbmi_6m )
res_cbmi_8m = mr(cbmi_8m )
res_cbmi_1.5y = mr(cbmi_1.5y )
res_cbmi_1y = mr(cbmi_1y )
res_cbmi_2y = mr(cbmi_2y )
res_cbmi_3y = mr(cbmi_3y )
res_cbmi_5y = mr(cbmi_5y )
res_cbmi_7y = mr(cbmi_7y )
res_cbmi_8y = mr(cbmi_8y )

individual_time_points = bind_rows(
  res_cbmi_birth,
  res_cbmi_6w,
  res_cbmi_3m,
  res_cbmi_6m,
  res_cbmi_8m,
  res_cbmi_1.5y,
  res_cbmi_1y,
  res_cbmi_2y,
  res_cbmi_3y,
  res_cbmi_5y,
  res_cbmi_7y,
  res_cbmi_8y)


# make nice plots
individual_time_points$exposure = factor(individual_time_points$exposure,levels=c("Birth","6w","3m","6m","8m","1y","1.5y","2y","3y","5y","7y","8y"))

library(ggstance)
p=ggplot(individual_time_points,aes(b,exposure,col=method))+geom_point(position=position_dodgev(height=0.5))+ggplot2::geom_errorbarh(mapping=aes(xmin=b-1.96*se,xmax=b+1.96*se,height=0.1),position=position_dodgev(height=0.5))+geom_vline(xintercept=0,alpha=0.1)+theme_bw()+
labs(x="Beta (MR causal estimate)",y="Childhood BMI epoch",col="MR method")

png("cbmi_mr_indiv_time_points.png",res=300,units="in",height=8,width=8)
p
dev.off()

write_csv(individual_time_points,"individual_time_points.csv")

betas = individual_time_points %>% filter(method=="Inverse variance weighted" | method == "Wald ratio")
betas = betas$b
MannKendall(betas)


# pleiotropy
pleio_cbmi_birth = mr_pleiotropy_test(cbmi_birth )
pleio_cbmi_6w = mr_pleiotropy_test(cbmi_6w )
pleio_cbmi_3m = mr_pleiotropy_test(cbmi_3m )
pleio_cbmi_6m = mr_pleiotropy_test(cbmi_6m )
pleio_cbmi_8m = mr_pleiotropy_test(cbmi_8m )
pleio_cbmi_1.5y = mr_pleiotropy_test(cbmi_1.5y )
pleio_cbmi_1y = mr_pleiotropy_test(cbmi_1y )
pleio_cbmi_2y = mr_pleiotropy_test(cbmi_2y )
pleio_cbmi_3y = mr_pleiotropy_test(cbmi_3y )
pleio_cbmi_5y = mr_pleiotropy_test(cbmi_5y )
pleio_cbmi_7y = mr_pleiotropy_test(cbmi_7y )
pleio_cbmi_8y = mr_pleiotropy_test(cbmi_8y )


pleio_timepoints = bind_rows(
  pleio_cbmi_birth,
  pleio_cbmi_6w,
  pleio_cbmi_3m,
  pleio_cbmi_6m,
  pleio_cbmi_8m,
  pleio_cbmi_1.5y,
  pleio_cbmi_1y,
  pleio_cbmi_2y,
  pleio_cbmi_3y,
  pleio_cbmi_5y,
  pleio_cbmi_7y,
  pleio_cbmi_8y)

  # heterogeneity
  hetero_cbmi_birth = mr_heterogeneity(cbmi_birth )
  hetero_cbmi_6w = mr_heterogeneity(cbmi_6w )
  hetero_cbmi_3m = mr_heterogeneity(cbmi_3m )
  hetero_cbmi_6m = mr_heterogeneity(cbmi_6m )
  hetero_cbmi_8m = mr_heterogeneity(cbmi_8m )
  hetero_cbmi_1.5y = mr_heterogeneity(cbmi_1.5y )
  hetero_cbmi_1y = mr_heterogeneity(cbmi_1y )
  hetero_cbmi_2y = mr_heterogeneity(cbmi_2y )
  hetero_cbmi_3y = mr_heterogeneity(cbmi_3y )
  hetero_cbmi_5y = mr_heterogeneity(cbmi_5y )
  hetero_cbmi_7y = mr_heterogeneity(cbmi_7y )
  hetero_cbmi_8y = mr_heterogeneity(cbmi_8y )


  hetero_timepoints = bind_rows(
    hetero_cbmi_birth,
    hetero_cbmi_6w,
    hetero_cbmi_3m,
    hetero_cbmi_6m,
    hetero_cbmi_8m,
    hetero_cbmi_1.5y,
    hetero_cbmi_1y,
    hetero_cbmi_2y,
    hetero_cbmi_3y,
    hetero_cbmi_5y,
    hetero_cbmi_7y,
    hetero_cbmi_8y)


write_csv(pleio_timepoints,"pleio_indiv_time_points.csv")
write_csv(hetero_timepoints,"hetero_indiv_time_points.csv")

````

# Negative control (birth weight, Warrington 2019)
````R
bw = read_table2("Fetal_BW_European_meta.NG2019.txt",col_types = cols_only(
  SNP = col_character(),
  ea = col_character(),
  nea = col_character(),
  eaf = col_double(),
  beta = col_double(),
  se = col_double(),
  p = col_double(),
  n = col_double(),
  rsid = col_character()
))

bw = bw %>%
select(-SNP) %>%
rename(
  "SNP" = rsid,
  "effect_allele"= ea,
  "other_allele"=nea,
  "pval"=p)

bw = bw %>% filter(SNP %in% ms_dat$SNP)
bw = bw %>% filter(pval < 5e-8)
bw = clump_data(bw)
bw = format_data(bw,type="exposure")
combo_dat = harmonise_data(bw,ms_dat)
res = mr(combo_dat)
mr_pleiotropy_test(combo_dat)
````
