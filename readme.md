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
````

Here are the results of the main IVW analysis:

Epoch   |   Beta    | SE    | P       | N SNPs
------- | -------   | ----- | ------- | --------
Birth to 6 weeks | -0.2089121 | 0.24341350 | 0.390748941 | 7
3 months to 1.5 years | 0.1766631 | 0.05558243 | 0.001480913 | 23
2 years to 5 years | 0.2670678 | 0.12266639 | 0.029466543 | 4
7 years to 8 years | 0.2890067 | 0.10974622 | 0.008453254 | 4

````R

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
Remove adult BMI SNPs
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

adult_bmi_sigsnps = adult_bmi %>%  filter(P<1e-5)

# remove Yengo SNPs associated with adult BMI at p<1e-5

# show excluded snps
cbmi_birth_6w %>% filter((SNP %in% adult_bmi_sigsnps$SNP))
cbmi_3mo_1.5y %>% filter((SNP %in% adult_bmi_sigsnps$SNP))
cbmi_2y_5y %>% filter((SNP %in% adult_bmi_sigsnps$SNP))
cbmi_7y_8y %>% filter((SNP %in% adult_bmi_sigsnps$SNP))

# exclude snps
cbmi_birth_6w = cbmi_birth_6w %>% filter(!(SNP %in% adult_bmi_sigsnps$SNP))
cbmi_3mo_1.5y = cbmi_3mo_1.5y %>% filter(!(SNP %in% adult_bmi_sigsnps$SNP))
cbmi_2y_5y = cbmi_2y_5y %>% filter(!(SNP %in% adult_bmi_sigsnps$SNP))
cbmi_7y_8y = cbmi_7y_8y %>% filter(!(SNP %in% adult_bmi_sigsnps$SNP))

# format
combo_cbmi_birth_6w = harmonise_data(cbmi_birth_6w,ms_dat)
combo_cbmi_3mo_1.5y = harmonise_data(cbmi_3mo_1.5y,ms_dat)
combo_cbmi_2y_5y = harmonise_data(cbmi_2y_5y,ms_dat)
combo_cbmi_7y_8y = harmonise_data(cbmi_7y_8y,ms_dat)

write_csv(combo_cbmi_birth_6w,"abmi_removed_combo_cbmi_birth_6w.csv")
write_csv(combo_cbmi_3mo_1.5y,"abmi_removed_combo_cbmi_3mo_1.5y.csv")
write_csv(combo_cbmi_2y_5y,"abmi_removed_combo_cbmi_2y_5y.csv")
write_csv(combo_cbmi_7y_8y,"abmi_removed_combo_cbmi_7y_8y.csv")

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
write_csv(mr_results,"mr_results_abmi_removed.csv")

# make nice plots
mr_results$epoch = factor(mr_results$epoch,levels=c("Birth to 6 weeks","3 months to 1.5 years","2 years to 5 years","7 years to 8 years"))

p=ggplot(mr_results,aes(b,epoch,col=method))+geom_point()+ggplot2::geom_errorbarh(mapping=aes(xmin=b-1.96*se,xmax=b+1.96*se,height=0.1))+geom_vline(xintercept=0,alpha=0.1)+theme_bw()+facet_wrap(~method,ncol=1)+
labs(x="Beta (MR causal estimate)",y="Childhood BMI epoch",col="MR method")

png("cbmi_mr_abmi_removed.png",res=300,units="in",height=8,width=8)
p
dev.off()
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
