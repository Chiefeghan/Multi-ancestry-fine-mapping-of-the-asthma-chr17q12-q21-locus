hyprcoloc runs for EUR 

```{r}


library(coloc)
library(Rfast)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(hyprcoloc)



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


#### Select rsid and BETA from each file
betas=cbind(df.b1[,c(1)], df.b1[,c(11)],df.b2[,c(11)],
	df.e1[,c(11)],df.e2[,c(11)],df.l1[,c(11)],
	df.l2[,c(11)],df.m1[,c(11)],df.m2[,c(11)],df.n1[,c(11)],
	df.n2[,c(11)],df.w1[,c(11)],df.ast4[,c(8)], 
	df.ast5[,c(8)],df.ast6[,c(8)],df.ast8[,c(8)], df.ast_sa[,c(5)])

## Select rsid and SE from each file
se=cbind(df.b1[,c(1)], df.b1[,c(12)],df.b2[,c(12)],
	df.e1[,c(12)],df.e2[,c(12)],df.l1[,c(12)],
	df.l2[,c(12)],df.m1[,c(12)],df.m2[,c(12)],
	df.n1[,c(12)],df.n2[,c(12)],df.w1[,c(12)],
	df.ast4[,c(9)], df.ast5[,c(9)],df.ast6[,c(9)],
	df.ast8[,c(9)], df.ast_sa[,c(6)] )

#rename
colnames(betas)=c("SNP", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14", "T15", "T16") ###(8056 snps)
colnames(se)=c("SNP", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14","T15", "T16"). #(8056 snps)

#track which phenotypes are in which column
#T1=BAS#
#T2=BAS%
#T3=EOS#
#T4=EOS%
#T5=LYM#
#T6=LYM%
#T7=MON#
#T8=MON%
#T9=NEU#
#T10=NEU%
#T11=WBC#
#T12=asthma4
#T13=asthma5
#T14=asthma6
#T15=asthma8
#T16=asthma_sa

## extract the QC'd snps only (3381)
EUR.snplist=fread('~/Desktop/chr17folder/eur_gwas/asthma6.z') ### snplist (3381snps) with rsid,chr, BP, a1,a2, maf 

betas=betas[betas$SNP %in% EUR.snplist$SNP, ]. ##3381 snps
se=se[se$SNP %in% EUR.snplist$SNP,] ##3381 snps
nrow(se)
nrow(betas)


##load the LD matrix obatined using LDstore for the 3381 snps
chr17.matrix=fread("~/Desktop/chr17folder/eur_gwas/asthma6.ld") dim(3381 3381)
chr17.matrix=as.matrix(chr17.matrix) ###convert to R matrix



#include term to account for sample overlap
sample.overlap = matrix(1, dim(betas)[2], dim(betas)[2]) 
View(sample.overlap)

####sample.overlap matrix is longer than number of phenotypes considered by 1 because of rsid column 
##exclude the rsid column in the matrix to make it a 16X16 matrix
sample.overlap=sample.overlap[-17, -17] ##n=16

traits <- paste0("T", 1:dim(betas[,-c(17)])[2])
rsid <- betas$SNP
#check rsid length
length(rsid)



################## add phenotypic correlation matrix using cpassoc Zscore mat.

df.b1=betas[df.b1$SNP %in% EUR.snplist$SNP, ] ##3381
df.b2=betas[df.b2$SNP %in% EUR.snplist$SNP, ] ##3381
df.e1=betas[df.e1$SNP %in% EUR.snplist$SNP, ] ##3381
df.e2=betas[df.e2$SNP %in% EUR.snplist$SNP, ] ##3381
df.l1=betas[df.l1$SNP %in% EUR.snplist$SNP, ] ##3381
df.l2=betas[df.l2$SNP %in% EUR.snplist$SNP, ] ##3381
df.m1=betas[df.m1$SNP %in% EUR.snplist$SNP, ] ##3381
df.n1=betas[df.n1$SNP %in% EUR.snplist$SNP, ] ##3381
df.n2=betas[df.n2$SNP %in% EUR.snplist$SNP, ] ##3381
df.w1=betas[df.w1$SNP %in% EUR.snplist$SNP, ] ##3381
df.ast4=betas[df.ast4$SNP %in% EUR.snplist$SNP, ] ##3381
df.ast5=betas[df.ast5$SNP %in% EUR.snplist$SNP, ] ##3381
df.ast6=betas[df.ast6$SNP %in% EUR.snplist$SNP, ] ##3381
df.ast8=betas[df.ast8$SNP %in% EUR.snplist$SNP, ] ##3381
df.ast_sa=betas[df.ast_sa$SNP %in% EUR.snplist$SNP, ] ##3381


## snplist
WBC_matrix=cbind(df.b1$BETA/df.b1$SE,df.b2$BETA/df.b2$SE,
	df.e1$BETA/df.e1$SE,df.e2$BETA/df.e2$SE,
	df.l1$BETA/df.l1$SE, df.l2$BETA/df.l2$SE, df.m1$BETA/df.m1$SE,
	df.m2$BETA/df.m2$SE,
	df.n1$BETA/df.n1$SE, df.n2$BETA/df.n2$SE,
	df.w1$BETA/df.w1$SE, df.ast4$BETA/df.ast$SE,
	df.ast5$BETA/df.ast5$SE, df.ast6$BETA/df.ast6$SE,
	df.ast8$BETA/df.ast8$SE, df.ast_sa$BETA/df.ast_sa$SE )

WBC_matrix=as.data.frame(WBC_matrix)
colnames(WBC_matrix)=c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", 
	"Z7","Z8","Z9","Z10","Z11","Z12", "Z13", "Z14", "Z15", "Z16")


Pheno.mat<- wbc.2[,.(Z1, Z2, Z3, Z4, Z5, Z6, Z7,Z8,Z9,Z10,Z11,Z12, Z13, Z14, Z15, Z16)]

##generate matrix 
trait.cor <- cor(Pheno.mat)
rownames(trait.cor)<-NULL
colnames(trait.cor)<-NULL
################################################################


binary.traits=c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1) ##0 is quantitative and 1 is binary trait (asthma), ordered according to phenotypes in file i.e T1...T16
betas=column_to_rownames(betas, var = "SNP")
se=column_to_rownames(se, var = "SNP")
betas=as.matrix(betas)
se=as.matrix(se)


#summary of results across all phenotypes
ptm1 = proc.time();
res <- hyprcoloc(betas, se, trait.names=traits, snp.id=rsid, trait.cor = trait.cor, 
	binary.outcomes = binary.traits,
	ld.matrix = chr17.matrix, 
	sample.overlap = sample.overlap, 
	uniform.priors = FALSE);
time.corr = proc.time() - ptm1;
res



##check snpscores
ptm1 = proc.time(); 
res2 <- hyprcoloc(betas, se, trait.names=traits, snp.id=rsid, trait.cor = trait.cor, 
	binary.outcomes = binary.traits,
	ld.matrix = chr17.matrix, 
	sample.overlap = sample.overlap, 
	uniform.priors = FALSE,snpscores = TRUE);
time.corr = proc.time() - ptm2;
res2



##Inspect credible sets 
cred.sets(res1, value = 0.95)
cred.sets(res2, value = 0.95)


```







