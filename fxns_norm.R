
library(tidyverse)
library(testit)

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
    ### see https://stackoverflow.com/questions/51110216/how-to-multiply-each-column-by-each-scalar-in-r
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