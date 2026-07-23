## Project: Tara Pacific sequencing - Clipperton atoll time-series
## Script purpose: Explore sequencing depth and community composition of original data
## Date: 23/07/2026
## Authors: James O'Brien

# Load packages
library(data.table)  # fread()
library(tidyverse)   # pipes, dplyr, ggplot2, tibble, forcats
library(rstudioapi)  # getActiveDocumentContext()
library(ggforce)     # geom_sina()
library(vegan)       # rrarefy(), vegdist()
library(ape)         # pcoa()

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
asv_abtab<-fread("../data/processed/asv_abtab_original.tsv",sep="\t",header=TRUE,data.table=FALSE) %>%
  column_to_rownames("id") %>% as.matrix()

metadat<-fread("../data/processed/metadat_original.tsv",header=TRUE,data.table=FALSE)

# Compute sequencing depth
seq_depth_df<-apply(asv_abtab,1,sum) %>% as.data.frame() %>% #sum all ASVs per sample
  rename("."="seq_depth") %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id") #append calculated sequencing depth to metadat

# Plot sequencing depth
ggplot(data=seq_depth_df,aes(x=fct_reorder(date,seq_depth,median),y=seq_depth,fill=date,colour=date)) +
  geom_violin(draw_quantiles=0.5,scale="width",alpha=0.5,show.legend=FALSE) +
  geom_sina(scale="width") +
  scale_y_log10() +
  ylab("Sequencing depth") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.position="none")
ggsave(filename="../data/processed/02_seq-depth.png",width=10,height=7,dpi=500)

# Rarefy at 1000 reads
rr_cutoff<-1000
asv_abtab_rr<-asv_abtab[which(apply(asv_abtab,1,sum)>=rr_cutoff),] #only keep samples >=1000 reads
asv_abtab_rr<-rrarefy(asv_abtab_rr,sample=rr_cutoff) #subsample to 1000 reads
asv_abtab_rr<-asv_abtab_rr[,which(apply(asv_abtab_rr,2,sum)>0)] #remove all-0 (non-present when downsampling) ASVs

# Compare communities (PCoA)
pcoa_res<-pcoa(D=vegdist(sqrt(asv_abtab_rr))) #sqrt-transform to lessen impact of abundant ASVs
pcoa_df<-pcoa_res$vectors %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(metadat,by="id")

# Adding Sample_type to discriminate samples from controls
pcoa_df <- pcoa_df %>%
  mutate(Sample_type = case_when(
      grepl("PCR", id, ignore.case = TRUE) ~ "PCR_negative_control",
      grepl("pos", id, ignore.case = TRUE) ~ "PCR_positive_control",
      grepl("ExB", id, ignore.case = TRUE) ~ "Extraction_blank",
      grepl("RB",  id, ignore.case = TRUE) ~ "Robot_blank",
      TRUE ~ "Sample"))

# Plot community distances (PCoA)
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=Sample_type)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.5) +
  theme_bw()
ggsave(filename="../data/processed/02_PCoA_by-sample-type.png",width=14,height=10,dpi=500)

# Plot community distances (PCoA) by sample type
ggplot(data=pcoa_df,aes(x=Axis.1,y=Axis.2,colour=Sample_type)) +
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2) +
  geom_point(size=4,alpha=0.5) +
  facet_wrap(~Sample_type, ncol=5) +
  guides(colour="none") +
  theme_bw()
ggsave(filename="../data/processed/02_PCoA_by-sample-type-faceted.png",width=14,height=10,dpi=500)