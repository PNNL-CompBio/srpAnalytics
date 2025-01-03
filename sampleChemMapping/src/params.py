from pandas import isna
# =========================================================
# Define Columns / Table Schemas
# =========================================================

REQUIRED_SAMPLE_COLUMNS = [
    "ClientName",
    "SampleNumber",
    "date_sampled",
    "sample_matrix",
    "technology",
    "projectName",
    "SampleName",
    "LocationLat",
    "projectLink",
    "LocationLon",
    "LocationName",
    "LocationAlternateDescription",
    "AlternateName",
    "cas_number",
    "date_sample_start",
    "measurement_value",
    "measurement_value_qualifier",
    "measurement_value_unit",
    "measurement_value_molar",
    "measurement_value_molar_unit",
    "environmental_concentration",
    "environmental_concentration_qualifier",
    "environmental_concentration_unit",
    "environmental_concentration_molar",
    "environmental_concentration_molar_unit",
]

# @TODO: need to rename the water columns
# new_sample_columns = {
#     'environmental_concentration': 'water_concentration',
#     'environmental_concentration_qualifier': 'water_concentration_qualifier',
#     'environmental_concentration_unit': 'water_concentration_unit',
#     'environmental_concentration_molar': 'water_concentration_molar',
#     'environmental_concentration_molar_unit': 'water_concentration_molar_unit'
# }


REQUIRED_COMPTOX_COLUMNS = [
    "INPUT",
    "DTXSID",
    "PREFERRED_NAME",
    "INCHIKEY",
    "SMILES",
    "MOLECULAR_FORMULA",
    "AVERAGE_MASS",
    "PUBCHEM_DATA_SOURCES",
]

SAMPLE_CHEM_COLUMNS = [
    "Sample_ID",
    "Chemical_ID",
    "measurement_value",
    "measurement_value_qualifier",
    "measurement_value_unit",
    "measurement_value_molar",
    "measurement_value_molar_unit",
    "environmental_concentration",
    "environmental_concentration_qualifier",
    "environmental_concentration_unit",
    "environmental_concentration_molar",
    "environmental_concentration_molar_unit",
]

SAMPLE_COLUMNS = [
    "Sample_ID",
    "ClientName",
    "SampleNumber",
    "date_sampled",
    "sample_matrix",
    "technology",
    "projectName",
    "SampleName",
    "LocationLat",
    "projectLink",
    "LocationLon",
    "LocationName",
    "LocationAlternateDescription",
    "AlternateName",
    "date_sample_start",
]

REQUIRED_CHEM_COLUMNS = [
    "Chemical_ID",
    "cas_number",
    "DTXSID",
    "PREFERRED_NAME",
    "INCHIKEY",
    "SMILES",
    "MOLECULAR_FORMULA",
    "AVERAGE_MASS",
    "PUBCHEM_DATA_SOURCES",
    "chem_source",
    "chemical_class",
    "chemDescription",
]

REQUIRED_MAPPING_COLUMNS = [
    "SampleNumber",
    "Sample_ID",
    "zf_lims_id",
]  # added these to complete mapping

# required from bmd calculation
# add in Sample_ID or Chemical_ID depending on table
REQUIRED_BMD_COLUMNS = {
    "bmd": [
        "Chemical_ID",
        "End_Point",
        "Model",
        "BMD10",
        "BMD50",
        "Min_Dose",
        "Max_Dose",
        "AUC_Norm",
        "DataQC_Flag",
        "BMD_Analysis_Flag",
    ],
    "doseRep": ["Chemical_ID", "End_Point", "Dose", "Response", "CI_Lo", "CI_Hi"],
    "fitVals": ["Chemical_ID", "End_Point", "X_vals", "Y_vals"],
}

# Define MASV mapping params
MASV_SOURCE = [
    "pharmacological",
    "personalCare",
    "industrial",
    "pulpAndPaper",
    "pestProduct",
    "natural",
    "pestRodenticide",
    "consumerProduct",
    "pestHerbicide",
    "pestInsecticide",
    "pestFungicide",
    "pestGeneral",
    "flameRetardant",
]

MASV_CC = [
    "PAH",
    "OPAH",
    "PCB",
    "PBB",
    "PBDE",
    "deuterated",
    "dioxinsAndFurans",
    "haloEthers",
    "OPFR",
    "phenol",
    "aniline",
    "uncategorized",
]

QC_FLAGS = {0: "Poor", 1: "Poor", 4: "Moderate", 5: "Moderate"}
# BMD_FLAGS ?
# QC_SCALE = {0: "Poor", 1: "Poor", 2: "Good", 3: "Good", 4: "Moderate", 5: "Poor"}


# Map the newClass values based on conditions provided
def map_classification(x):
    if x in [
        "industrial",
        "industrial; aniline",
        "Industrial",
        "industrial; consumerProduct; phenol",
        "industrial; consumerProduct; aniline",
        "industrial; phenol",
    ]:
        return "Industrial"
    elif x in ["PAH; industrial", "PAH"]:
        return "PAH"
    elif x in [
        "personalCare; personalCare; natural; natural; consumerProduct; consumerProduct",
        "personalCare; natural; consumerProduct",
        "personalCare; natural",
        "pharmacological; personalCare; industrial; natural; consumerProduct",
    ]:
        return "Natural"
    elif x == "pestFungicide":
        return "Fungicide"
    elif isna(x) or x == "NA":
        return "Unclassified"
    else:
        return x
