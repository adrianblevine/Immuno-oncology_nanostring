
library(here)
suppressPackageStartupMessages({
library(tidyverse)
library(Biobase)
library(ggpubr)
require(janitor)
})

# TIS weights file, with columns for GENE and Weights
TIS_FILE <- here("ref/example_tis_weights_file.csv")

# expression set of normalized gene counts with HGNC gene identifiers saved
# as minimal row data, which must match GENE column in TIS_FILE
ESET_RDS <- here("results/eset_norm.rds")

eset <- readRDS(ESET_RDS)
edata <- as_tibble(exprs(eset)) 

# load weights 
weights <- read_csv(TIS_FILE) 

# filter fdata for endogenous probes and remove unnecessary columns
fdata <- fData(eset) %>% 
    janitor::clean_names() %>% 
    filter(code_class=="Endogenous") %>% 
    # fill missing HGNC symbol for IGL gene to avoid error later
    mutate(hgnc_symbol = coalesce(hgnc_symbol, gene_name)) %>% 
    dplyr::select(hgnc=hgnc_symbol, probe=probe_id)

# select TIS genes and confirm proper ordering wrt weights
dat <- as.data.frame(edata)
rownames(dat) <- fdata$hgnc
dat <- dat[weights$GENE,]
stopifnot(identical(rownames(dat), weights$GENE))

# calculate TIS as matrix multiplication of log2 transformed gene counts 
# with weights
tis <- t(log2(dat)) %*% weights$Weights

# format as tibble
pdata <- pData(eset) %>% 
    janitor::clean_names() %>% 
    dplyr::select(file_name, study_id, sample_type=sample_type_2)

# join with annotation from pdata
scores <- tibble(file_name = names(dat), TIS = tis[,1]) %>% 
    left_join(pdata, by="file_name")

# save TIS scores
write_tsv(scores, here("results/TIS_scores.tsv"))

# save example boxplot
ggboxplot(scores, x="sample_type", y="TIS", color="sample_type",
          add="jitter", palette="nejm", legend="none",
          xlab="Sample Type", ylab="Tumor Inflammation Signature",
          title="Example TIS boxplot",
          subtitle="Note: for workflow testing only")
ggsave(here("results/TIS_boxplot.pdf"), h=4,w=4)
