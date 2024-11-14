# This file contains some cleaning and checks before uploading raw data to figshare

library(tidyverse)
library(data.table)

#####################
## MORPHOLOGY DATA ##
#####################

# Set the working directory 
setwd("~/Git_Repos/srpAnalytics/zfBmd/files/")

# Read all csv files 
files <- list.files(".", full.names = T) %>%
  .[!grepl(".xlsx", .)]

# Process datasets--------------------------------------------------------------

#' Converts MM/DD/YYYY to YYYYMMDD
date_fixer <- function(x) {
  y <- strsplit(x, "/") %>% unlist()
  month <- y[1]
  day <- y[2]
  year <- y[3]
  year <- ifelse(nchar(year) == 2, paste0("20", year), year)
  month <- ifelse(nchar(month) == 1, paste0("0", month), month)
  day <- ifelse(nchar(day) == 1, paste0("0", day), day)
  return(paste0(year, month, day))
}

## "./Single_CSV_10JUL29-01-018_03-07-2017.csv"
file1 <- fread(files[1]) %>%
  rename(plate.id = `Plate ID`, chemical.id = `Chemical ID`, conc = CONC, 
         well = WELL, date = DATE) %>%
  select(-c(CASRN, comments)) %>%
  mutate(date = map_chr(date, date_fixer)) %>%
  pivot_longer(6:ncol(.)) %>% 
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("X", "", conc), 
         plate.id = as.character(plate.id)),
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
head(file1)
unique(file1$endpoint)

## "./Single_CSV_A090024_03-07-2017.csv"
file2 <- fread(files[2]) %>%
  rename(plate.id = `Plate ID`, chemical.id = `Chemical ID`, conc = CONC, 
         well = WELL, date = DATE) %>%
  select(-c(CASRN, comments)) %>%
  mutate(date = map_chr(date, date_fixer)) %>%
  pivot_longer(6:ncol(.)) %>% 
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("X", "", conc), 
                           plate.id = as.character(plate.id)),
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
head(file2)
unique(file2$endpoint)

## "./Single_CSV_Anderson2016_Samples_2018AUG29.csv"
file3 <- fread(files[3]) %>%
  rename(plate.id = plate) %>%
  pivot_longer(7:ncol(.)) %>% 
  rename(endpoint = name)
unique(file3$endpoint)

## "./Single_CSV_File_125 PAH_FOR RELEASE_12-10-201567369.csv"
file4 <- fread(files[4]) %>%
  rename(chemical.id = `Chemical ID`, plate.id = `Barcode ID`, well = WELL, conc = CONC,
         date = DATE) %>%
  select(-CASRN) %>%
  mutate(date = gsub("-", "", date)) %>%
  pivot_longer(6:ncol(.)) %>% 
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM|X", "", conc)),
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file4$endpoint)

## "./Single_CSV_File_138_PAH_FOR_RELEASE_2016DEC07.csv"
file5 <- fread(files[5]) %>%
  rename(chemical.id = `Chemical.ID`, conc = CONC, plate.id = Plate, well = WELL) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM|X", "", conc)),
         date = "20161207",
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file5$endpoint)

## "./Single_CSV_File11-08-201591139_Ph2013.csv"
file6 <- fread(files[6]) %>%
  rename(plate.id = `Barcode ID`, chemical.id = `Chemical ID`, well = WELL, conc = CONC,
         date = DATE) %>%
  select(-c(CASRN)) %>%
  pivot_longer(6:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("X_TMT|X", "", conc)),
         date = gsub("-", "", date),
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file6$endpoint)

## "./Single_CSV_FSES_1400_Overlap_Toxcast_ZF_Raw_Data_2015OCT30.csv"
file7 <- fread(files[7]) %>%
  rename(plate.id = `Barcode.ID`, chemical.id = `Chemical.ID`, well = WELL, conc = CONC) %>%
  pivot_longer(5:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM", "", conc)),
         date = "20151030",
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file7$endpoint)

## "./Single_CSV_FSES_19unique_chem_2018MAY22.csv"
file8 <- fread(files[8]) %>%
  rename(plate.id = `Plate.ID`, chemical.id = `Chemical.ID`, well = WELL, conc = CONC,
         date = DATE) %>%
  select(-c(CASRN, comments)) %>%
  pivot_longer(6:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM", "", conc)),
         bottle.id = NA,
         date = map_chr(date, date_fixer)) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file8$endpoint)

## "./Single_Tab_19missing_samples_ 04-01-2019 180923.csv"
file9 <- fread(files[9]) %>%
  pivot_longer(7:ncol(.)) %>%
  rename(endpoint = name, plate.id = plate) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file9$endpoint)

## "./Single_Tab_5_new_PAH  09-24-2018 144441.csv" 
file10 <- fread(files[10]) %>%
  pivot_longer(7:ncol(.)) %>%
  rename(endpoint = name, plate.id = plate) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file10$endpoint)

## "./Single_Tab_File_3PAH_Supermix3_03-07-2017 223014.csv"
file11 <- fread(files[11]) %>%
  rename(plate.id = `Plate ID`, chemical.id = `Chemical ID`, well = WELL, conc = CONC,
         date = DATE) %>%
  select(-c(CASRN, comments)) %>%
  pivot_longer(6:ncol(.)) %>%
  rename(endpoint = name) %>%
  mutate(conc = as.numeric(gsub("uM|X", "", conc)),
         date = map_chr(date, date_fixer),
         bottle.id = NA) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file11$endpoint)

## "./Single_Tab_Seal_Coat  10-02-2018 225135.csv"
file12 <- fread(files[12]) %>%
  pivot_longer(7:ncol(.)) %>%
  rename(endpoint = name, plate.id = plate) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file12$endpoint)

## "./Single_Tab_St Mary Bioremediation  10-17-2018 150650.csv"
file13 <- fread(files[13]) %>%
  pivot_longer(7:ncol(.)) %>%
  rename(endpoint = name, plate.id = plate) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file13$endpoint)

## "./Single_Tab_Thermal Remediatione  09-19-2018 174508.csv"
file14 <- fread(files[14]) %>%
  pivot_longer(7:ncol(.)) %>%
  rename(endpoint = name, plate.id = plate) %>%
  select(chemical.id, bottle.id, conc, plate.id, well, date, endpoint, value)
unique(file14$endpoint) 

##  "./Tanguay Phase 4 zf 104alkyl PAH morphology data PNNL 2023OCT05.csv"
#file15 <- fread(files[15]) 

# Combine datasets--------------------------------------------------------------

# bind them all together
all_data <- rbind(file1, file2, file3, file4, file5, file6, file7, file8, file9,
                  file10, file11, file12, file13, file14)

# Run checks
all_data <- all_data %>% 
  mutate(date = ifelse(date == "120720", "20120720", date),
         date = ifelse(date == "110307", "20110307", date))

# Check dates 
date_check <- all_data$date %>% unique() %>% as.Date("%Y%m%d")
min(date_check)
max(date_check)

# Check other variables
unique(all_data$chemical.id) # There can be no NA's in this category
unique(all_data$bottle.id)
unique(all_data$conc)
unique(all_data$plate.id)
unique(all_data$well)
unique(all_data$endpoint)
unique(all_data$value) # Only 0, 1, and NA allowed

# Remove any NA concentrations
all_data <- all_data %>%
  filter(!is.na(conc))

fwrite(all_data, "~/Downloads/Zfish_Morphology_Legacy_2011-2018.csv", quote = F, row.names = F)

################
## MAPPING ID ##
################

# Load zebrafish data 
zfish <- fread("~/Downloads/Zfish_Morphology_Legacy_2011-2018.csv") %>%
  select(chemical.id, bottle.id) %>%
  unique()

# Read and concatenate mapping files
mapfiles <- list.files("~/Downloads/mappings/", full.names = T)

all_maps <- rbind( 
  fread(mapfiles[1]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    select(-`Project Name`) %>% mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[2]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    select(-`Project Name`) %>% mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[3]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[4]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[5]) %>% rename(chemical.name = TestSubstance.ChemicalName, chemical.id = Chemical.ID, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[6]),
  fread(mapfiles[7]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[8]) %>% select(-project_name),
  fread(mapfiles[9]),
  fread(mapfiles[10]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = Chemical.ID, casrn = TestSubstance_CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[11]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = Chemical.ID, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[12]),
  fread(mapfiles[13]),
  fread(mapfiles[14]),
  fread(mapfiles[15]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn),
  fread(mapfiles[16]) %>% rename(chemical.name = TestSubstance_ChemicalName, chemical.id = `Chemical ID`, casrn = CASRN) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn), 
  fread(mapfiles[17]) %>% rename(chemical.name = ChemicalName, casrn = CAS, chemical.id = `Zebrafish ID`) %>%
    mutate(bottle.id = NA) %>% select(bottle.id, chemical.id, chemical.name, casrn)
) %>%
  filter(chemical.id %in% zfish$chemical.id) %>%
  mutate(chemical.name = gsub(",", "_", chemical.name),
         chemical.name = gsub("?", "", chemical.name, fixed = T),
         casrn = trimws(casrn),
         casrn = gsub("/", "-", casrn)) %>%
  unique() 

# Remove naming inconsistencies
all_maps <- all_maps %>%
  filter(chemical.name %in% c("9-aminophenanthrene", "PAH Mix", "Naphtho[12-b]fluoranthene",
                              "Supermix 10 [50x]", "_A130485___") == FALSE)

# Find and select first duplicate
dupes <- table(all_maps$chemical.id, dnn = "chemical.id") %>% data.frame() %>% filter(Freq > 1)

cleaned_maps <- rbind(
  all_maps %>%
    filter(chemical.id %in% dupes$chemical.id) %>%
    arrange(chemical.id, casrn) %>%
    group_by(chemical.id) %>%
    slice_head(n = 1),
  all_maps %>%
    filter(chemical.id %in% dupes$chemical.id == FALSE) 
)

# Make sure there is none missing! 
zfish$chemical.id[zfish$chemical.id %in% cleaned_maps$chemical.id == FALSE]

# Last checks
cleaned_maps <- cleaned_maps %>%
  mutate(casrn = ifelse(casrn == "", NA, casrn))

cleaned_maps$bottle.id %>% unique()
cleaned_maps$chemical.id %>% unique()
cleaned_maps$chemical.name %>% unique()
cleaned_maps$casrn %>% unique()

sum(is.na(cleaned_maps$casrn))

fwrite(cleaned_maps, "~/Downloads/Zfish_Chemical_Mappings_Legacy_2011-2018.csv", quote = F, row.names = F)





