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

EM = read.xlsx("EndpointMetadata_Ver_1_1.xlsx", 1) %>% 
  select(2:4) %>%
  setNames(c("End_Point", "End_Point_Type", "End_Point_Name"))

chemical_fix_bmd = function(data) {
  data %>%
    left_join(EM) %>%
    select(Chemical_ID, Model, BMD10, BMD10_Flag, BMD50, BMD50_Flag, Min_Dose,
           Max_Dose, AUC_Norm, End_Point_Name, End_Point_Type, DataQC_Flag)
}
chemical_fix_dose = function(data)  {
  data %>% 
    left_join(EM) %>% 
    select(-End_Point) %>%
    mutate(Chemical_ID = as.numeric(Chemical_ID))
}

# Chemicals---------------------------------------------------------------------

# BMDs
zebrafishChemBMDs = bind_rows(
  fread("../zfBmd/outputs/continuous/continuous_BenchmarkDose.csv") %>%
    rename(cas = Chemical_ID) %>%
    left_join(chemical_metadata %>% select(cas, chemical_id) %>% rename(Chemical_ID = chemical_id)) %>%
    chemical_fix_bmd(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_BMDs_BC.csv") %>%
    chemical_fix_bmd(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_BMDs_LPR.csv") %>%
    chemical_fix_bmd()
)
fwrite(zebrafishChemBMDs, "~/Downloads/all_srp_data/srpCompendiumV5/zebrafishChemBMDs.txt",
       quote = F, row.names = F, sep = "\t")

# Dose
zebrafishChemDoseResponse = bind_rows(
  fread("../zfBmd/outputs/continuous/continuous_Dose.csv") %>% 
    rename(cas = Chemical_ID) %>%
    left_join(chemical_metadata %>% select(cas, chemical_id) %>% rename(Chemical_ID = chemical_id)) %>%
    chemical_fix_dose(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_Dose_BC.csv") %>% chemical_fix_dose(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_Dose_LPR.csv") %>% chemical_fix_dose()
) %>%
  select(Chemical_ID, Dose, Response, CI_Lo, CI_Hi, End_Point_Name, End_Point_Type)
fwrite(zebrafishChemDoseResponse, "~/Downloads/all_srp_data/srpCompendiumV5/zebrafishChemDoseResponse.txt",
       quote = F, row.names = F, sep = "\t")

# Coords
zebrafishChemXYCoords = bind_rows(
  fread("../zfBmd/outputs/continuous/continuous_Fits.csv") %>%
    rename(cas = Chemical_ID) %>%
    left_join(chemical_metadata %>% select(cas, chemical_id) %>% rename(Chemical_ID = chemical_id)) %>%
    chemical_fix_dose(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_Fits_BC.csv") %>% chemical_fix_dose(),
  fread("../zfBmd/outputs/zebrafish/zebrafish_chem_Fits_LPR.csv") %>% chemical_fix_dose()
) %>%
  select(Chemical_ID, X_vals, Y_vals, End_Point_Name, End_Point_Type) %>%
  filter(!is.na(X_vals) & !is.na(Y_vals))
fwrite(zebrafishChemXYCoords, "~/Downloads/all_srp_data/srpCompendiumV5/zebrafishChemXYCoords.txt",
       quote = F, row.names = F, sep = "\t")

# Samples-----------------------------------------------------------------------

# Convert sample identifiers
id_converter = fread("~/Downloads/all_srp_data/metadata/SampleMetadata_Ver_1_1.csv") %>%
  rename(Wrong_ID = Sample_ID) %>%
  full_join(samples_final %>% select(Sample_ID, SampleNumber), by = "SampleNumber")
  
# BMDs
zebrafishSampBMDs = rbind(
  fread("../zfBmd/outputs/sample/sample_chem_BMDs_BC.csv"),
  fread("../zfBmd/outputs/sample/sample_chem_BMDs_LPR.csv")
) %>%
  rename(Wrong_ID = Chemical_ID) %>%
  left_join(id_converter, by = "Wrong_ID") %>%
  filter(!is.na(Sample_ID)) %>%
  left_join(EM) %>%
  select(Sample_ID, Model, BMD10, BMD10_Flag, BMD50, BMD50_Flag, Min_Dose,
         Max_Dose, AUC_Norm, End_Point_Name, End_Point_Type, DataQC_Flag)
fwrite(zebrafishSampBMDs, "~/Downloads/all_srp_data/srpCompendiumV5/zebrafishSampBMDs.txt",
       quote = F, row.names = F, sep = "\t")



