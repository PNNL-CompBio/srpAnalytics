library(tidyverse)
library(data.table)

# Set the working directory 
setwd("~/Git_Repos/srpAnalytics/zfBmd/files/")

# Read all csv files 
files <- list.files(".", full.names = T) %>%
  .[!grepl(".xlsx", .)]

# Process datasets--------------------------------------------------------------

# Process files one at a time 
standard_clean <- function(x) {
  fread(x) %>%
    rename(chemical = `Chemical ID`, plate = `Plate ID`, well = WELL, conc = CONC) %>%
    select(-c(DATE, CASRN, comments)) %>%
    pivot_longer(5:ncol(.)) %>% 
    rename(endpoint = name) %>%
    mutate(conc = as.numeric(gsub("X", "", conc), plate = as.character(plate)))
}

## "./Single_CSV_10JUL29-01-018_03-07-2017.csv"
file1 <- standard_clean(files[1])
head(file1)
unique(file1$endpoint)

## "./Single_CSV_A090024_03-07-2017.csv"
file2 <- standard_clean(files[2])
head(file2)
unique(file2$endpoint)

## "./Single_CSV_Anderson2016_Samples_2018AUG29.csv"
file3 <- fread(files[3]) %>%
  rename(chemical = chemical.id) %>%
  select(-c(bottle.id, date)) %>%
  pivot_longer(5:ncol(.)) %>% 
  rename(endpoint = name)
unique(file3$endpoint)

## "./Single_CSV_File_125 PAH_FOR RELEASE_12-10-201567369.csv"
file4 <- fread(files[4]) %>%
  rename(chemical = `Chemical ID`, plate = `Barcode ID`, well = WELL, conc = CONC) %>%
  select(-c(CASRN, DATE)) %>%
  pivot_longer(5:ncol(.)) %>% 
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM|X", "", conc)))
unique(file4$endpoint)

## "./Single_CSV_File_138_PAH_FOR_RELEASE_2016DEC07.csv"
file5 <- fread(files[5]) %>%
  rename(chemical = `Chemical.ID`, conc = CONC, plate = Plate, well = WELL) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM|X", "", conc)))
unique(file5$endpoint)

## "./Single_CSV_File11-08-201591139_Ph2013.csv"
file6 <- fread(files[6]) %>%
  rename(plate = `Barcode ID`, chemical = `Chemical ID`, well = WELL, conc = CONC) %>%
  select(-c(CASRN, DATE)) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("X_TMT|X", "", conc)))
unique(file6$endpoint)

## "./Single_CSV_FSES_1400_Overlap_Toxcast_ZF_Raw_Data_2015OCT30.csv"
file7 <- fread(files[7]) %>%
  rename(plate = `Barcode.ID`, chemical = `Chemical.ID`, well = WELL, conc = CONC) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM", "", conc)))
unique(file7$endpoint)

## "./Single_CSV_FSES_19unique_chem_2018MAY22.csv"
file8 <- fread(files[8]) %>%
  rename(plate = `Plate.ID`, chemical = `Chemical.ID`, well = WELL, conc = CONC) %>%
  select(-c(CASRN, DATE, comments)) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM", "", conc)))
unique(file8$endpoint)

# Combine datasets--------------------------------------------------------------

# bind them all together
all_data <- rbind(file1, file2, file3, file4, file5, file6, file7, file8)




