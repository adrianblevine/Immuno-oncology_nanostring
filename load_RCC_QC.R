
library(tidyverse)
library(NanoStringQCPro)
require(Biobase)

# all RCC files should be in a single folder
rcc_files <- dir("/path/to/rcc/dir", full.names = TRUE)

# RLF file is NanoString panel specific and is required to load RCC files
RLF_FILE <- "/path/to/file.rlf"

# optional additional annotation data
PDATA_FILE <- "/path/to/pdate_file.tsv"

# QC cutoff for geo mean of housekeeping genes
HK_CUTOFF <- 100 

# FUNCTIONS ---------------------------------------------------------------

geo_mean <- function(x) {
    exp(sum(log(x))/length(x))
}

# load RCC files ----------------------------------------------------------

# extraPdata: files must be tab- separated and contain a column labelled 
# "FileName" whose values correspond exactly to the rcc filenames 
# (including .RCC extension) 
# Note: pseudocount of 1 is added to all measurements to enable subsequent 
# log transformation of the data

rcc_set <- newRccSet(rcc_files,
                     rlf = RLF_FILE,
                     extraPdata = c(LABEL_FILE),
                     addEgAnnotations = TRUE, 
                     dropPdataCols = c("FileVersion", "SoftwareVersion", 
                                       "Owner", "SystemAPF", 
                                       "ScannerID", "CartridgeBarcode"), 
)

# QC plots ---------------------------------------------------------------

# filter for Housekeeping genes
hk <- rcc_set[fData(rcc_set)$CodeClass == "Housekeeping",]

# calculate geometric mean
hk <- apply(exprs(hk), 2, geo_mean)

# format as tibble and flag samples before specified cutoff
hk <- tibble("sample" = names(hk),
             "hk_gm" = hk,
             "flag" = map_dbl(hk, function(x) ifelse(x < HK_CUTOFF, 1, 0)))
hk

# QC report ---------------------------------------------------------------
### optional step to generate more comprehensive QC report

# rcc_set_pp <- preprocRccSet(rcc_set,
#                             doBackground = T, bgReference="negatives", 
#                             doPosCtrlNorm = TRUE, pcnSummaryFunction = geo_mean,
#                             doContentNorm = TRUE, normMethod = "housekeeping",
#                             normSummaryFunction = geo_mean, hkgenes = NULL,
#                             hkfeatures = NULL,
#                             doPresAbs = TRUE, paStringency = 2,
#                             quietly = FALSE)
# 
# # Results for intermediate preprocessing steps are included in additional 
# # matrices in assayData
# ls(assayData(rcc_set_pp))
# 
# # settings for each preprocessing step are recorded in correspondingly 
# # named elements in the preprocessing slot
# preproc(rcc_set_pp)
# 
# qc_example_rccSet <- makeQCReport(rcc_set_pp, 
#                                   outputBaseName = "nanostringqcpro_report",
#                                   outputDir = OUTP,
#                                   preprocOverride = F)

# save expression set ---------------------------------------------------

edata <- as.matrix(exprs(rcc_set))
pdata <- as.data.frame(pData(rcc_set)) %>% 
    Biobase::AnnotatedDataFrame()

fdata <- as.data.frame(fData(rcc_set)) %>% 
    Biobase::AnnotatedDataFrame()

assert(identical(colnames(edata), rownames(pdata)) )

eset <- Biobase::ExpressionSet(assayData = edata,
                               phenoData = pdata,
                               featureData = fdata)

idx <- hk$flag != 1
eset_pass <- eset[,idx]

saveRDS(eset_pass, "/path/to/eset_raw.rds")

