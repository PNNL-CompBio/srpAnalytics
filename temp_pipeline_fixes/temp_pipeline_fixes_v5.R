library(tidyverse)
library(data.table)
library(xlsx)

######################
## SAMPLES PIPELINE ##
######################

# USE VERSION 4 SAMPLE DATA-----------------------------------------------------

## SAMPLES METADATA ##
samples_final = fread("~/Downloads/samples_mlb.csv") %>%
  select(Sample_ID, ClientName, SampleNumber, date_sampled, sample_matrix,
         technology, projectName, SampleName, LocationLat, LocationLon,
         LocationName, projectLink, AlternateName) %>%
  filter(projectName != "NULL") 

samples_final %>%
  fwrite("~/Downloads/all_srp_data/srpCompendiumV5/samples.txt", quote = F, row.names = F, sep = "\t")

## SAMPLE MEASUREMENTS: SAMPLES TO CHEMICALS ## 

# Read chemical metadata
chemical_metadata = read.xlsx("ChemicalMetadata_Ver_1_2_DD.xlsx", 1) %>%
  mutate(cas = gsub("`", "", cas))

# Load the accounted for 
samplesToChemicals_accounted = rbind(
  fread("~/Downloads/all_srp_data/cleaned_data_v2/EnvSampler_Sample_20230101_FigShare263116-30545993_StHelensStudyAirData.txt"),
  fread("~/Downloads/all_srp_data/cleaned_data_v2/EnvSampler_Sample_20210427_FigShare263116-31173814_FSESDataForPNNL.txt"),
  fread("~/Downloads/all_srp_data/cleaned_data_v2/EnvSampler_Sample_20230101_FigShare263116-30545993_StHelensStudyAirData.txt")
) %>%
  select(-Sample_ID) %>%
  left_join(samples_final %>% select(SampleNumber, Sample_ID)) %>%
  filter(!is.na(Sample_ID)) %>%
  select(-Chemical_ID) %>%
  left_join(
    chemical_metadata %>% select(chemical_id, cas) %>% rename(Chemical_ID = chemical_id, cas_number = cas)
  ) %>%
  select(Sample_ID, Chemical_ID, measurement_value, measurement_value_qualifier,
         measurement_value_unit, measurement_value_molar, environmental_concentration,
         environmental_concentration_qualifier, environmental_concentration_unit,
         environmental_concentration_molar, environmental_concentration_molar_unit,
         test_template, test_method) %>%
  filter(!is.na(Chemical_ID)) # Remove mixtures and pavement studies

# Find any missing data that is not a mixture or a pavement measurement 
samplesToChemicals_unaccounted = fread("~/Downloads/all_srp_data/31197685/samplesToChemicals.csv") %>%
  filter(Sample_ID %in% samplesToChemicals_accounted$Sample_ID == FALSE) %>%
  filter(Sample_ID %in% samples_final$Sample_ID) %>%
  select(-Chemical_ID) %>%
  left_join(
    chemical_metadata %>% select(chemical_id, cas) %>% rename(Chemical_ID = chemical_id, cas_number = cas)
  ) %>%
  mutate(
    Chemical_ID = ifelse(cas_number == "1506-02-1, 21145-77-7", 10340, Chemical_ID),
    Chemical_ID = ifelse(cas_number == "202-33-5 & 199-54-2", 10901, Chemical_ID),
    Chemical_ID = ifelse(cas_number == "205-83-4 & 238-04-0", 2327, Chemical_ID),
    Chemical_ID = ifelse(cas_number == "575-43-9 & 575-41-7", 66, Chemical_ID),
    Chemical_ID = ifelse(cas_number == "77392-71-3 , 198-55-0", 3757, Chemical_ID)
  ) %>%
  select(Sample_ID, Chemical_ID, measurement_value, measurement_value_qualifier,
         measurement_value_unit, measurement_value_molar, environmental_concentration,
         environmental_concentration_qualifier, environmental_concentration_unit,
         environmental_concentration_molar, environmental_concentration_molar_unit,
         test_template, test_method)

# Bind the two
samplesToChemicals_final = rbind(samplesToChemicals_accounted, samplesToChemicals_unaccounted) %>%
  filter(Sample_ID %in% samples_final$Sample_ID)

# Check chemicals
all(unique(samplesToChemicals_final$Chemical_ID) %in% chemical_metadata$chemical_id)

samplesToChemicals_final %>%
  fwrite("~/Downloads/all_srp_data/srpCompendiumV5/samplesToChemicals.txt", quote = F, row.names = F, sep = "\t")

########################
## CHEMICALS PIPELINE ##
########################

# Extract chemicals
chemical_metadata %>%
  select(-chemical_class) %>%
  rename(
    Chemical_ID = chemical_id, cas_number = cas, DTXCID = dtxcid, 
    PREFERRED_NAME = preferredName, INCHIKEY = inchikey, 
    SMILES = smiles, MOLECULAR_FORMULA = molFormula,
    AVERAGE_MASS = averageMass, chemical_class = chem_class_WEB
  ) %>%
  select(Chemical_ID, cas_number, DTXCID, PREFERRED_NAME, INCHIKEY,
         SMILES, MOLECULAR_FORMULA, AVERAGE_MASS, chemical_class, chemDescription) %>%
  fwrite("~/Downloads/all_srp_data/srpCompendiumV5/chemicals.txt", quote = F, row.names = F, sep = "\t")

#############################
## BENCHMARK DOSE PIPELINE ##
#############################

# Chemicals---------------------------------------------------------------------

# Samples-----------------------------------------------------------------------




