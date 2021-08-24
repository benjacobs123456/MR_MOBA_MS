# Age-specific effects of childhood BMI on multiple sclerosis risk: a Mendelian Randomisation study

Requirements
- [MS GWAS data](https://imsgc.net/?page_id=31)
- [Childhood BMI data](https://www.medrxiv.org/node/346664.external-links.html)
- TwoSampleMR

Load the packages we'll need
````R
library(TwoSampleMR)
library(dplyr)
library(readr)
library(ggplot2)
library(readr)
````

Load the childhood BMI GWAS and format
````R
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
cbmi_birth_6w = cbmi %>% filter(age %in% c("Birth","6w")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE)
cbmi_3mo_1.5y = cbmi %>% filter(age %in% c("3m","6m","8m","1y","1.5y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE)
cbmi_2y_5y = cbmi %>% filter(age %in% c("2y","3y","5y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE)
cbmi_7y_8y = cbmi %>% filter(age %in% c("7y","8y")) %>% arrange(pval) %>% distinct(SNP,.keep_all=TRUE)

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
# format
combo_cbmi_birth_6w = harmonise_data(cbmi_birth_6w,ms_dat)
combo_cbmi_3mo_1.5y = harmonise_data(cbmi_3mo_1.5y,ms_dat)
combo_cbmi_2y_5y = harmonise_data(cbmi_2y_5y,ms_dat)
combo_cbmi_7y_8y = harmonise_data(cbmi_7y_8y,ms_dat)

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

mr_results$epoch = factor(mr_results$epoch,levels=c("Birth to 6 weeks","3 months to 1.5 years","2 years to 5 years","7 years to 8 years"))

p=ggplot(mr_results,aes(b,epoch,col=method))+geom_point(mapping=aes(size=-log10(pval)))+geom_errorbarh(mapping=aes(xmin=b-1.96*se,xmax=b+1.96*se,height=0.1))+geom_vline(xintercept=0,alpha=0.1)+theme_bw()+facet_wrap(~method,ncol=1)+
labs(x="Beta (MR causal estimate)",y="Childhood BMI epoch",col="MR method")

png("cbmi_mr_fig1.png",res=300,units="in",height=8,width=8)
p
dev.off()

mr_results %>% filter(method=="Inverse variance weighted")


````

Here are the results of the main IVW analysis:

Epoch   |   Beta    | SE    | P       | N SNPs
------- | -------   | ----- | ------- | --------
Birth to 6 weeks | -0.2089121 | 0.24341350 | 0.390748941 | 7
3 months to 1.5 years | 0.1766631 | 0.05558243 | 0.001480913 | 23
2 years to 5 years | 0.2670678 | 0.12266639 | 0.029466543 | 4
7 years to 8 years | 0.2890067 | 0.10974622 | 0.008453254 | 4

````R

mr_scatter_plot(res_epoch1,dat=combo_cbmi_birth_6w)
mr_scatter_plot(res_epoch2,dat=combo_cbmi_birth_6w)
mr_scatter_plot(res_epoch3,dat=combo_cbmi_birth_6w)
mr_scatter_plot(res_epoch4,dat=combo_cbmi_birth_6w)
````

Steiger filtering
````R

# steiger and moe
combo_cbmi_birth_6w$samplesize.exposure=28681
combo_cbmi_birth_6w$prevalence.outcome=14802/(14802+26703)
combo_cbmi_birth_6w$units.outcome="log odds"
combo_cbmi_birth_6w$units.exposure="z score"
combo_cbmi_birth_6w$r_out =get_r_from_lor(lor=combo_cbmi_birth_6w$beta.outcome,prevalence=combo_cbmi_birth_6w$prevalence.outcome,ncase = 14802,ncontrol=26703,af =combo_cbmi_birth_6w$eaf.exposure)
combo_cbmi_birth_6w$r_exp = get_r_from_pn(p=combo_cbmi_birth_6w$pval.exposure,n=combo_cbmi_birth_6w$samplesize.exposure)

mr_steiger(
  p_exp = combo_cbmi_birth_6w$pval.exposure,
  p_out = combo_cbmi_birth_6w$pval.outcome,
  r_exp = combo_cbmi_birth_6w$r_exp,
  r_out = combo_cbmi_birth_6w$r_out,
  n_exp = combo_cbmi_birth_6w$samplesize.exposure,
  n_out = (14802+26703)
  )

# Load the downloaded RData object. This loads the rf object
load("rf.rdata?dl=0")

# Obtain estimates from all methods, and generate data metrics
combo_dat$units.outcome="log odds"
combo_dat$units.exposure="z score"
combo_dat$prevalence.outcome = 37688/(37688+981372)
combo_dat$ncase.outcome=37688
combo_dat$ncontrol.outcome=981372

combo_dat$id.exposure="Lymphocyte count"
combo_dat$id.outcome="PD"
combo_dat$r.outcome = get_r_from_lor(combo_dat$beta.outcome,combo_dat$eaf.outcome,combo_dat$ncase.outcome,combo_dat$ncontrol.outcome,combo_dat$prevalence.outcome)
combo_dat$r.exposure = get_r_from_pn(combo_dat$pval.exposure,combo_dat$samplesize.exposure)
combo_dat = combo_dat %>% filter(mr_keep=="TRUE")

steiger = steiger_filtering(combo_dat) %>% filter(steiger_dir == "TRUE") %>% select(SNP)
combo_dat = combo_dat %>% filter(SNP %in% steiger$SNP)

````
