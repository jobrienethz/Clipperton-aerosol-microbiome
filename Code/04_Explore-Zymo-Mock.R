# Project: Tara Pacific sequencing - Clipperton atoll time-series
## Script purpose: Explore composition of and contamination in Zymo mock community samples with tidied metadata
## Date: 23/07/2026
## Authors: Fabienne Wiederkehr, James O'Brien

# Load packages
library(data.table)  # fread()
library(tidyverse)   # dplyr, tidyr, tibble, ggplot2, pipes, etc.
library(vegan)       # decostand(), vegdist()
library(Biostrings)  # readDNAStringSet()
library(rstudioapi)  # getActiveDocumentContext()

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
asv_dat<-fread("../data/processed/asv_dat_original.tsv",sep="\t",header=TRUE,data.table=FALSE)
metadat<-fread("../data/processed/metadat_original.tsv",header=TRUE,data.table=FALSE) %>%
  filter(id %in% colnames(asv_dat)[4:ncol(asv_dat)])

# Identifying positive controls in the samples and linking to metadata
pos <- data.frame(
  id = names(asv_dat)[grepl("pos", names(asv_dat), ignore.case = TRUE)]
)
# Identify which pos IDs are missing from metadat
missing_pos <- pos$id[!pos$id %in% metadat$id]
# Create a minimal metadata frame for them
new_rows <- data.frame(id = missing_pos)
# Bind them to metadat
metadat <- bind_rows(metadat, new_rows)

# Set mock samples and load Zymo data
mock_samples<-metadat %>%
  filter(grepl("pos",id)) %>%
  select(id) %>%
  pull()
zymo_fasta<-readDNAStringSet("../data/processed/zymobiomics_D6311/ssrRNAs/all_ssrRNAs.fasta")
zymo_conc_exp<-fread("../data/processed/zymobiomics_D6311/theoretical_concentrations.tsv",sep=",",header=TRUE,data.table=FALSE)

# Search exact sequence of each ASV against the Zymo mock community members
asv_zymo<-NULL
for (i in 1:nrow(asv_dat)){
  cat("Searching ASV",i,"in the Zymo mock community...")
  if (length(which(grepl(asv_dat$seq[i],zymo_fasta)))>0){
    tmp<-data.frame(asv_id=asv_dat$asv_id[i],zymo_member=names(zymo_fasta)[which(grepl(asv_dat$seq[i],zymo_fasta))]) 
    asv_zymo<-asv_zymo %>%
      bind_rows(tmp)
  }
  cat("DONE\n")
}

asv_zymo<-asv_zymo %>%
  separate(zymo_member,into=c("zymo_genus","zymo_spp",NA,NA),sep="_") %>%
  unite(zymo_member,zymo_genus,zymo_spp,sep=" ") %>%
  distinct()

asv_zymo %>%
  arrange(zymo_member)

# Match the ASVs to the Zymo members and join the expected abundances
zymo_conc_obs<-asv_dat %>%
  select(asv_id,all_of(mock_samples)) %>%
  column_to_rownames(var="asv_id") %>%
  decostand(method="total",MARGIN=2) %>% as.data.frame() %>%
  rownames_to_column(var="asv_id") %>%
  pivot_longer(-asv_id,names_to="sample",values_to="abund_obs") %>%
  full_join(asv_zymo,by="asv_id") %>%
  select(-asv_id) %>%
  group_by(zymo_member,sample) %>%
  dplyr::summarise(abund_obs=sum(abund_obs)) %>%
  mutate(zymo_member=ifelse(is.na(zymo_member),"All not mock members",zymo_member)) %>%
  mutate(percent_obs=abund_obs*100) %>%
  left_join(zymo_conc_exp,by="zymo_member") %>%
  mutate(conc_16S_exp=ifelse(is.na(conc_16S),0,conc_16S)) %>%
  mutate(conc_16S_and_18S_exp=ifelse(is.na(conc_16S_and_18S),0,conc_16S_and_18S))

# Plot observed vs. expected abundances for all samples
ggplot(data=zymo_conc_obs,aes(x=conc_16S_exp,y=percent_obs,colour=zymo_member)) +
  geom_abline() +
  geom_point(size=3,alpha=0.7) +
  facet_wrap(~sample) +
  scale_x_continuous(trans="log10",breaks=c(0.0001,0.01,0.1,1,100),labels=paste(c("0.0001","0.01","0.1","1","100"),"%",sep="")) +
  scale_y_continuous(trans="log10",breaks=c(0.0001,0.01,0.1,1,100),labels=paste(c("0.0001","0.01","0.1","1","100"),"%",sep="")) +
  xlab("Theoretical 16S concentration in Zymo mock community\n(% abundance)") +
  ylab("Observed % abundance") +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave(filename="../data/processed/04_mock-observed-vs-expected-abundances.png",width=10,height=5,dpi=500)