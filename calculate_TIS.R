
library(tidyverse)
library(Biobase)

TIS_FILE <- "/path/to/tis_weights.csv"
ESET_RDS <- "/path/to/expression_set.rds"

eset <- readRDS(ESET_RDS)
fdata <- fData(eset)
pdata <- pData(eset)
edata <- exprs(eset) %>% as_tibble()


# load weights and make table
weights <- read_csv(TIS_FILE) %>%
    mutate(GENE = map_chr(GENE, function(x) str_replace(x, "\\.", "-")))

# filter fdata for endogenous probes and remove unnecessary columns
fdata <- filter(fdata, CodeClass=="Endogenous") %>% 
    clean_names() %>% 
    dplyr::select(hgnc=hgnc, probe=probe)

# select TIS genes and confirm proper ordering wrt weights
dat <- as.data.frame(norm_)
rownames(dat) <- fdata$hgnc
dat <- dat[weights$GENE,]
assert(identical(rownames(dat), weights$GENE))

# calculate TIS as matrix multiplication with weights
tis <- t(log2(dat)) %*% weights$Weights

# format as tibble
scores <- tibble(study_id = names(dat), TIS = tis[,1])



