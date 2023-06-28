
library(tidyverse)
library(testit)
library(Biobase)

ESET_RDS <- "/path/to/eset_raw.rds"

eset <- readRDS(ESET_RDS)
fdata <- fData(eset)
pdata <- pData(eset)
edata <- exprs(eset) %>% as_tibble()

# FUNCTIONS -----------------------------------------------------------

separate_controls <- function(edata, fdata, column="Code.Class"){
    positive <- filter(edata, fdata[column] == "Positive") 
    hk <- filter(edata, fdata[column] == "Housekeeping") 
    endogenous <- filter(edata, fdata[column] == "Endogenous")
    negative <- filter(edata, fdata[column] == "Negative") 
    dat <- list(pos = positive, neg = negative, 
                hk = hk, endog = endogenous)
    return(dat)
}

geo_mean <- function(x) exp(sum(log(x))/length(x))

geo_mean_scale <- function(dat, ctrl) {
    geo_means <- apply(ctrl, 2, geo_mean)
    av_geo <- mean(geo_means)
    fac <- av_geo/geo_means
    norm_ <- mapply(`*`, dat, fac) %>% as_tibble()
    return(norm_)
}

normalize_rna <- function(edata, fdata, round=T) {
    n_hk <- sum(fdata$CodeClass == "Housekeeping") 
    n_endog <- sum(fdata$CodeClass == "Endogenous")
    dat <- separate_controls(edata, fdata, column="CodeClass")
    
    # KEY IS TO APPLY POSITIVE CONTROL SCALING TO THE HOUSEKEEPING GENES
    dat1 <- rbind(dat$endog, dat$hk)
    norm1 <- geo_mean_scale(dat1, dat$pos)
    dat2 <- norm1[1:n_endog, ]
    # use scaled HK genes for sample content normalization
    hk <- norm1[(n_endog + 1):(n_endog+n_hk), ]
    norm2 <- geo_mean_scale(dat2, hk)
    if (round) norm2 <- round(norm2)
    return(norm2)
}

# normalize  -----------------------------------------------------------

norm_ <- normalize_rna(edata, fdata) %>% 
    as_tibble()

# then can make and save expression set as shown in load_RCC_QC.R using the
# pdata and fdata from original expression set, but removing the housekeeping
# genes from fdata

# NanoStringNorm (to compare) ---------------------------------------------
# Can also verify that these normalization functions produce the same result
# as the NanoStringNorm packages (which is the same as the default 
# NSolver normalization)
# This is no longer supported and can be installed from the archive at 
# https://cran.r-project.org/src/contrib/Archive/NanoStringNorm/
# 
# library(NanoStringNorm) # 
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




