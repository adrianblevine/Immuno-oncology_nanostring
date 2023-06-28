# Loads data from RCC files, normalizes samples with standard 
# NCounter method that I coded, and compares to the NanoStringNorm 
# package as a sanity check.
# Prints tables summarizing the most important categories
# Calculates TIS and adds to pdata
# Saves as csv files, as well as an expression set
# 
# 2022-April-13: also saves refernce set for clinical pipeline

# Preamble ----------------------------------------------------------------

library(janitor)
library(tidyverse)
library(gt)


#library(NanoStringNorm)
library(clusterProfiler) # need for entrez gene identifiers
library(org.Hs.eg.db) # need for entrez gene identifiers

source(here("R/fxns_misc.R"))    
source(here("R/fxns_norm.R"))

OUTP <- here("results/12-13")
dir.create(OUTP, showWarnings=F)

### hard coded variables #################################################

ESET_RDS <- here("data/eset_raw.rds")

LABEL_FILE <- here("data/nanostring_log.tsv")



# load data ---------------------------------------------------------------

eset <- readRDS(ESET_RDS)
fdata <- fData(eset)
pdata <- pData(eset)
edata <- exprs(eset) %>% as_tibble()

# make reference set for pipeline -----------------------------------------

fdata_ref <- dplyr::select(fdata, -c(Accession))
pdata_ref <- dplyr::select(pdata, FileName, study_id,
                           sample_type, subtype, molecular, sample_id)
edata_ref <- edata
rownames(edata_ref) <- fdata_ref$probe
colnames(edata_ref) <- pdata_ref$study_id

eset_ref <- make_expression_set(edata_ref, pdata_ref, fdata_ref)
saveRDS(eset_ref, here("data/pipeline_ref_eset.rds"))


# normalize  -----------------------------------------------------------

norm_ <- normalize_rna(edata, fdata) %>% as_tibble()

# NanoStringNorm (to compare) ---------------------------------------------
# verify that my normalization functions produce the same result
# as the (now unsupported) NanoStringNorm packages (which is the
# same as the default NSolver normalization)

# 
# anno <-  fData(eset) %>% 
#     dplyr::rename("Name" = "probe", "Code.Class" = "CodeClass")
# rownames(anno) <- anno$Name
# 
# x <- as.data.frame(edata)
# rownames(x) <- anno$Name
# 
# nsnorm <- NanoStringNorm(
#         x = x,
#         anno = anno,
#         CodeCount = 'geo.mean',
#         Background = "none",
#         SampleContent = 'housekeeping.geo.mean',
#         OtherNorm = "none",
#         round.values = TRUE,
#         is.log = FALSE,
#         take.log = F,
#         return.matrix.of.endogenous.probes = TRUE,
#         verbose=F
#     ) %>% as_tibble()
# 
# assert( identical(nsnorm, norm_))




