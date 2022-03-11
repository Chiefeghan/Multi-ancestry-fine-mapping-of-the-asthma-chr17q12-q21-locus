#! /usr/bin/env Rscript
---
title: Figure-Manuscript- upset plot/overlap of 95% CS in EUR
author: "Chief Ben-Eghan"
date: "2/2/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# R Markdown
Load the required packages
```{r}
library(dplyr)
library(data.table)
library(ggplot2)
library(ggalluvial)
library(alluvial)
library(ggpubr)
library(grid)
library(reshape2)
library(grid)
library(devtools)
source_gist("524eade46135f6348140")
library(ggpmisc)
library(cowplot)
library(susieR)
library(UpSetR)
```


```{r}

##load 95% CS in EUR

data_wide4=fread("~/Desktop/data_wide4")


pdf("~/Desktop/upsetplot1.pdf", 12,10)
upset(data_wide4,sets =c("<18yrs", "<40yrs", "asthma(age+sex)","BASO#","BASO%","EOS#","EOS%", "LYM#","LYM%","MON#","MON%", "NEU#","NEU%","WBC#"),keep.order = TRUE,order.by = "freq",
  query.legend = "bottom", nsets = 15, number.angles = 0, point.size = 2.5, line.size = 1,nintersects=50,text.scale = c(2, 1.5, 1.8, 1.5, 1.5, 1.3),
  mainbar.y.label = "Overlapping variants in 95% CS", sets.x.label ="# of SNPs per Phenotype" ,mb.ratio = c(0.55, 0.45),
  queries = list(
  list(
    query = intersects,
    params = list("<18yrs", "asthma(age+sex)", "MON%", "<40yrs"), 
    color = "#56B4E9", 
    active = T,
    query.name = "rs4795399"),
  list(query=intersects, params=list("asthma(age+sex)", "EOS#", "EOS%", "LYM#", "LYM%", "<40yrs"),
       color="#E69F00",
       active=T, 
       query.name="rs112401631", size=4)))

dev.off()

```
