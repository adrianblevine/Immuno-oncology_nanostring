# Need to remove additional files (not in log) before running this

library(here)
library(formattable)
library(gt)
library(janitor)
library(tidyverse)

library(NanoStringQCPro)

source(here("R/fxns_misc.R"))    
source(here("R/fxns_norm.R"))    

OUTP <- here("results/12-13/qc")
dir.create(OUTP, showWarnings=F)

### hard coded variables #################################################

LABEL_FILE <- here("data/nanostring_log.tsv")

RCC_DIR <- here("data/RCC_files")
rcc_files <- dir(RCC_DIR, full.names = TRUE)

RLF_FILE <- here("ref/Hawkins_Imm_1_C2903+A3378.rlf")

# generate QC report with nanostring QC pro - very helpful, but takes a while
QC_REPORT <- F
HK_CUTOFF <- 100 # QC cutoff for geo mean of housekeeping genes

# Load RCC files ----------------------------------------------------------

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

pdata <- pData(rcc_set)

hk <- rcc_set[fData(rcc_set)$CodeClass == "Housekeeping",]
hk <- apply(exprs(hk), 2, geo_mean)
hk <- tibble("sample" = names(hk),
             "hk_gm" = hk,
             "flag" = map_dbl(hk, function(x) ifelse(x < HK_CUTOFF, 1, 0)),
             "year" = pdata$sample_year,
             "sample_type" = pdata$sample_type)
write_tsv(hk, file.path(OUTP, "housekeeping_QC.tsv"))

ggplot(hk) + geom_boxplot(aes(x=sample_type, y=hk_gm))
ggsave(file.path(OUTP, "hk_genes_by_sample_type.png"))

### plot endogenous genes
endog <- rcc_set[fData(rcc_set)$CodeClass == "Endogenous",]
endog <- apply(exprs(endog), 2, geo_mean)
endog <- tibble("sample" = names(endog),
                "endog_gm" = endog)
df <- bind_cols(endog, sample_type=pData(rcc_set)$sample_type)
ggplot(df) + geom_boxplot(aes(x=sample_type, y=endog_gm))
ggsave(file.path(OUTP, "endog_genes_by_sample_type.png"))

### plot endogenous:HK ratio
vec <- endog$endog_gm/hk$hk_gm
df<- tibble("sample" = names(vec),
            "div" = vec) %>% 
    bind_cols(sample_type=pData(rcc_set)$sample_type)
ggplot(df) + geom_boxplot(aes(x=sample_type, y=div))
ggsave(file.path(OUTP, "endog_hk_ratio_by_sample_type.png"))

### HK genes by sample age
df <- filter(hk, !is.na(year))
ggplot(df) + geom_point(aes(x=year, y=hk_gm))
ggsave(file.path(OUTP, "hk_gm_by_sample_year.png"))

### breakdown of samples failing qc
df <- filter(hk, flag==1)
df <- table(df$sample_type) %>% data.frame()
colnames(df) <- c("sample type", "N")
pdf(file.path(OUTP, "qc_failure_table.pdf"))
gridExtra::grid.table(df)
dev.off()

# QC report ---------------------------------------------------------------

rcc_set_pp <- preprocRccSet(rcc_set,
                            doBackground = T, bgReference="negatives", 
                            doPosCtrlNorm = TRUE, pcnSummaryFunction = geo_mean,
                            doContentNorm = TRUE, normMethod = "housekeeping",
                            normSummaryFunction = geo_mean, hkgenes = NULL,
                            hkfeatures = NULL,
                            doPresAbs = TRUE, paStringency = 2,
                            quietly = FALSE)

# Results for intermediate preprocessing steps are included in additional 
# matrices in assayData
ls(assayData(rcc_set_pp))

# settings for each preprocessing step are recorded in correspondingly 
# named elements in the preprocessing slot
preproc(rcc_set_pp)

#rcc_set = getBackground(rcc_set, bgReference="negatives")
#contentNorm(rccSet=bgcorr_example_rccSet, method="housekeeping") >


if (QC_REPORT) {
    qc_example_rccSet <- makeQCReport(rcc_set_pp, 
                                      outputBaseName = "nanostringqcpro_report",
                                      outputDir = OUTP,
                                      preprocOverride = F,
                                      
    )
}

# Tidy pdata and consolidate labels --------------------------------------

edata <- exprs(rcc_set) %>% as_tibble()

# need Accession for NanoStringNorm
fdata <- fData(rcc_set) %>% as_tibble() %>% 
    dplyr::select(probe = GeneName, CodeClass, Accession, hgnc=Hgnc_Symbol)

cols <- names(pData(rcc_set)[13:length(pData(rcc_set))])

pdata <- pData(rcc_set) %>% as_tibble() %>% 
    na_if("") %>% 
    mutate_if(is.character, str_squish) %>% 
    dplyr::select(FileName, cols) %>%
    dplyr::select(FileName, study_id, sample_type, subtype, 
                  molecular, run=run_1, sample_year, sample_id)# %>% 
#    mutate(study_id = make_clean_names(study_id, case="all_caps")) %>% 

names(edata) <- pdata$study_id

### make sure names line up
vec <- map_chr(names(edata), function(x) str_split(x, "\\.")[[1]][1])
setdiff(pdata$study_id, vec)


# summarize data ----------------------------------------------------------

tab <- dplyr::count(pdata, sample_type) 
write_tsv(tab, file.path(OUTP, "breakdown_sample_type.tsv"))
tab <- gt(tab) %>% 
    tab_style(style = list(cell_text(weight = "bold")),
          locations = cells_column_labels())
gtsave(tab, file.path(OUTP, "breakdown_sample_type.png"))

for (type in c("LGG", "HGG", "Ependymoma", "Medulloblastoma", "MMRD")) {
    df <- filter(pdata, sample_type==type)
    tab <- dplyr::count(df, subtype)
    write_tsv(tab, file.path(OUTP, paste0("breakdown_subtype_", type, ".tsv")))
    tab <- gt(tab) %>% 
        tab_header(title=type) %>% 
        tab_style(style = list(cell_text(weight = "bold")),
                  locations = cells_column_labels())
    gtsave(tab, file.path(OUTP, paste0("breakdown_subtype_", type, ".png")))
}

# save --------------------------------------------------------------------

# only probe without Hgnc label is IGL
idx <- is.na(fdata$hgnc)
fdata$hgnc[idx] <- fdata$probe[idx]

eset <- make_expression_set(edata,pdata,fdata)

idx <- hk$flag != 1
eset_pass <- eset[,idx]

write_rds(eset_pass, here("data/eset_raw.rds"))


