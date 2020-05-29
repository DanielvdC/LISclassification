# Packages ----------------------------------------------------------------

library(readr)
library(stringr)
library(eeguana)
library(osfr)
library(caret)
library(dplyr)
library(eegkit)
library(purrr)
library(devtools)
library(tidyr)
library(tibble)
library(factoextra)
library(Metrics)
library(data.table)

# Variable and path definition -------------------------------------------------

pre_files <- 'INPUT_FILES'
N400_c <- c('Cz', 'CP1', 'CP2', 'P3', 'Pz', 'P4', 'POz')
N400_t <- c(250, 600)
N400_c_txt <- paste0(N400_c, collapse = ", ")

# Preprocessing and PCA per .RDS file ------------------------------------------

list.files(pre_files, 'RDS')

eeg_data_ica <- map_dfr(
  list.files(pre_files, "RDS"),
  function(f) {
    RDSfile <- paste0(pre_files,"/",f)
    message(sprintf("Loading RDS file %s...", RDSfile))
    eeg <- readRDS(RDSfile)
    eeg <- eeg %>% eeg_events_to_NA(.type=='artifact')
    eeg <- eeg %>% select(one_of(N400_c))
    eeg <- eeg_downsample(eeg, q = 4, n = 8, ftype = 'iir', multiple_times = FALSE)
    eeg <- eeg %>% dplyr::filter(between(as_time(.sample, "milliseconds"), N400_t[1], N400_t[2]))
    eeg <- eeg %>% dplyr::filter(constraint == 'Constraining' & region %in% regions)
    temp <- as_tibble(eeg$.signal[,3:9])
    eeg <- c(eeg$.signal, eeg$.segments)
    princ <- prcomp(temp)
    pca <- predict(princ, newdata = temp)[,1]
    eeg['PC1'] <- as_tibble(pca)
    eeg <- eeg[-c(3:32, 34:46)]
    eeg <- eeg[c(1,2, 4, 3)]
    eeg$.sample <- as.factor(eeg$.sample)
    eeg <- as.data.frame(eeg)
    saveRDS(eeg, file = RDSfile)
  }
)

# Merge preprocessed files together in 1 dataframe -----------------------------
# Without subject numbering
c <- 0
for (file in list.files('PREPROCESSED_FILES/')){
  if (c == 0){
    RDSfile <- paste0('PREPROCESSED_FILES/', file)
    df <- readRDS(RDSfile)
    df$.sample <- as.factor(df$.sample)
    df <- as.data.frame(df)
    c = c + 1
  }
  else{
    RDSfile <- paste0('PREPROCESSED_FILES/', file)
    df1 <- readRDS(RDSfile)
    df1$.sample <- as.factor(df1$.sample)
    df1 <- as.data.frame(df1)
    df <- rbind(df, df1)
  }
}

#With subject numbering
c = 0
for (file in list.files('PREPROCESSED_FILES/')){
  if (c == 0){
    RDSfile <- paste0('PREPROCESSED_FILES/', file)
    dataset <- readRDS(RDSfile)
    dataset['subject'] <- c
    c = c + 1
  }
  else {
    RDSfile <- paste0('PREPROCESSED_FILES/', file)
    temp_dataset <- readRDS(RDSfile)
    temp_dataset['subject'] <- c
    dataset <- rbind(dataset, temp_dataset)
    rm(temp_dataset)
    c = c + 1
  }
}

# Write complete dataframe to .csv file ----------------------------------------
write.csv(dataset, file = 'csv_name.csv')
eeg_data <- read.csv('csv_name.csv')
