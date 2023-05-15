COlOC-SuSiE.R

```{r}
library(coloc)
library(Rfast)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)


### EUR COlOC-SuSiE

#1Load full BOLT summary statistics file (SAIGE for asthma(sex+age))
#2Extract only snps that pass QC from full summary stats file
#3Use only snps that pass QC to generate LD matrix using LDstore

df.df.b1=fread("~/Desktop/chr17folder/eur_gwas/df.df.b1")#BAS# (full chr17 file, 8,056 SNPs, study region after QC 3381snps)
df.df.b2=fread("~/Desktop/chr17folder/eur_gwas/df.df.b2")#BAS% 
df.e1=fread("~/Desktop/chr17folder/eur_gwas/df.e1")#EOS#
df.e2=fread("~/Desktop/chr17folder/eur_gwas/df.e2")#EOS%
df.l1=fread("~/Desktop/chr17folder/eur_gwas/df.l1")#LYM#
df.l2=fread("~/Desktop/chr17folder/eur_gwas/df.l2")#LYM%
df.m1=fread("~/Desktop/chr17folder/eur_gwas/df.m1")#MON#
df.m2=fread("~/Desktop/chr17folder/eur_gwas/df.m2")#MON%
df.n1=fread("~/Desktop/chr17folder/eur_gwas/df.n1")#NEU
df.n2=fread("~/Desktop/chr17folder/eur_gwas/df.n2")#NEU%
df.w1=fread("~/Desktop/chr17folder/eur_gwas/df.w1")#WBC#
df.ast4=fread("~/Desktop/chr17folder/eur_gwas/df.ast4")#18-40yrs
df.ast5=fread("~/Desktop/chr17folder/eur_gwas/df.ast5")#>40yrs
df.ast6=fread("~/Desktop/chr17folder/eur_gwas/df.ast6")#<18yrs
df.ast8=fread("~/Desktop/chr17folder/eur_gwas/df.ast8")#<40yrs
df.ast_sa=fread("~/Desktop/chr17folder/eur_gwas/df.age_sex")#asthma_s+a


###Load LD matrix for study region(37,281,157 - 38,876,409bp, hg19) generated using LDstoreV2
# use same matrix for COlOC-SuSiE, hyprcoloc and SuSiE
chr17.matrix=fread("~/Desktop/chr17folder/eur_gwas/asthma6.ld") dim(3381 3381)
chr17.matrix=as.matrix(chr17.matrix) ###convert to R matrix




#WBC trait vs WBC trait Colocalization test
##create a COlOC-SuSiE file for each phenotype
### Load the snp.list (the same snp.list used to generate the LDmatrix)

####BAS#
EUR.snplist=fread('~/Desktop/chr17folder/eur_gwas/asthma6.z') ### snplist (3381snps) with rsid,chr, BP, a1,a2, maf 
EUR.snplist=EUR.snplist[,c(1:6)]
EUR.snplist$BETA=df.b1[match(EUR.snplist$rsid, df.b1$SNP), "BETA"]
EUR.snplist$SE=df.b1[match(EUR.snplist$rsid, df.b1$SNP), "SE"]
EUR.snplist$varbeta=EUR.snplist$SE^2 #variance of beta
EUR.snplist$type="quant"  #### use "cc" for binary trait
EUR.snplist$N="408963"


## add the rsids to the matrix or the script wont run 
dimnames(chr17.matrix) <- list(EUR.snplist$rsid, EUR.snplist$rsid)
 

p1=EUR.snplist$rsid
p2=EUR.snplist$BETA
p3=EUR.snplist$varbeta
p4=EUR.snplist$position
p5=EUR.snplist$type
p6=EUR.snplist$maf
p7=as.numeric(EUR.snplist$N)
df1=list(p1,p2,p3,p4,p5,p6,p7,chr17.matrix)
names(df1) <- c("snp", "beta", "varbeta", "position","type","MAF","N", "LD") ##BAS#



###BAS%
EUR.snplist=fread('~/Desktop/chr17folder/eur_gwas/asthma6.z')
EUR.snplist=EUR.snplist[,c(1:6)]
EUR.snplist$BETA=df.b2[match(EUR.snplist$rsid, df.b2$SNP), "BETA"]
EUR.snplist$SE=df.b2[match(EUR.snplist$rsid, df.b2$SNP), "SE"]
EUR.snplist$varbeta=EUR.snplist$SE^2
EUR.snplist$type="quant" #### use "cc" for binary trait
EUR.snplist$N="408963"
p1=EUR.snplist$rsid
p2=EUR.snplist$BETA
p3=EUR.snplist$varbeta
p4=EUR.snplist$position
p5=EUR.snplist$type
p6=EUR.snplist$maf
p7=as.numeric(EUR.snplist$N)
df2=list(p1,p2,p3,p4,p5,p6,p7,chr17.matrix)
names(df2) <- c("snp", "beta", "varbeta", "position","type","MAF","N", "LD") ##BAS%

###EOS#
EUR.snplist=fread('~/Desktop/chr17folder/eur_gwas/asthma6.z')
EUR.snplist=EUR.snplist[,c(1:6)]
EUR.snplist$BETA=df.e1[match(EUR.snplist$rsid, df.e1$SNP), "BETA"]
EUR.snplist$SE=df.e1[match(EUR.snplist$rsid, df.e1$SNP), "SE"]
EUR.snplist$varbeta=EUR.snplist$SE^2
EUR.snplist$type="quant" #### use "cc" for binary trait
EUR.snplist$N="408963"
p1=EUR.snplist$rsid
p2=EUR.snplist$BETA
p3=EUR.snplist$varbeta
p4=EUR.snplist$position
p5=EUR.snplist$type
p6=EUR.snplist$maf
p7=as.numeric(EUR.snplist$N)
df3=list(p1,p2,p3,p4,p5,p6,p7,chr17.matrix)
names(df3) <- c("snp", "beta", "varbeta", "position","type","MAF","N", "LD") ##EOS#





#Run susie for each phenotype (first 3 WBC traits as an example )
s1.eur=runsusie(df1)# BAS#
s2.eur=runsusie(df2)# BAS%
s3.eur=runsusie(df3)# EOS#
#s4.eur=runsusie(df4)# EOS%
#s5.eur=runsusie(df5)# LYM#
#s6.eur=runsusie(df6)# LYM%
#s7.eur=runsusie(df7)# MON#
#s8.eur=runsusie(df8)# MON%
#s10.eur=runsusie(df10)# NEU%
#s11.eur=runsusie(df11)# WBC#
#s12.eur=runsusie(df12)# ast4
#s13.eur=runsusie(df13)# ast5#
#s14.eur=runsusie(df14)# ast6%
#s15.eur=runsusie(df15)# ast8#
#s16.eur=runsusie(df16)# ast_s+a



##### run COLOC-SuSiE in pairs for each trait-pair e.g BAS# and BAS%; BAS# and EOS#
if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res1=coloc.susie(s1.eur,s2.eur)
  print(susie.res1$summary)
}

View(susie.res1$summary)
## add label for each run 
susie.res1$summary$trait-pair="BAS#+BAS%"


if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res2=coloc.susie(s1.eur,s3.eur)
  print(susie.res2$summary)
}

View(susie.res2$summary)
## add label for each run 
susie.res2$summary$trait-pair="BAS#+EOS#"



###### bind all files 
all_susie.EUR=rbind(susie.res1$summary,susie.res2$summary) ##for BAS#-BAS% and BAS#-EOS# only


## add run number 

all_susie.EUR %>%
  all_susie.EUR=mutate(all_susie.EUR, 
    run = case_when(
      `trait-pair` == "BAS#$+BAS%" ~ 'EUR:1',
      `trait-pair` == "BAS#+EOS#" ~ 'EUR:2',))


  fwrite(all_suisie.EUR, "~/Desktop/chr17folder/eur_gwas/coloc_susie_runs/test1.coloc.susie.txt", sep="\t", row.names=F, quote=F)

```


