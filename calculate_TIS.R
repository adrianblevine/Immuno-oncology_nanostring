
library(tidyverse)
library(Biobase)

# TIS weights file, with columns for GENE and Weights
TIS_FILE <- "/path/to/tis_weights.csv"

# expression set of normalized gene counts with HGNC gene identifiers saved
# as minimal row data, which must match GENE column in TIS_FILE
ESET_RDS <- "/path/to/expression_set.rds"

eset <- readRDS(ESET_RDS)
fdata <- fData(eset)
pdata <- pData(eset)
edata <- as_tibble(exprs(eset)) 

# load weights 
weights <- read_csv(TIS_FILE) 

# filter fdata for endogenous probes and remove unnecessary columns
fdata <- filter(fdata, CodeClass=="Endogenous") %>% 
    janitor::clean_names() %>% 
    dplyr::select(hgnc=hgnc, probe=probe)

# select TIS genes and confirm proper ordering wrt weights
dat <- as.data.frame(edata)
rownames(dat) <- fdata$hgnc
dat <- dat[weights$GENE,]
assert(identical(rownames(dat), weights$GENE))

# calculate TIS as matrix multiplication of log2 transformed gene counts 
# with weights
tis <- t(log2(dat)) %*% weights$Weights

# format as tibble
scores <- tibble(study_id = names(dat), TIS = tis[,1])



