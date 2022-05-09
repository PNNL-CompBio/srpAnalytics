#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os, random, sys, time

def format(mfile):
    """
    This module subsets morphology data to relevant columns, calculates sums (total, number of NA, and sum)
    for downstream dose response calculations, removes plates with >50% missingness in their baseline, rolls up
    to the chemical id and concentration level, and adds missing endpoints. 
    
    Attributes
    ----
    mfile: Path to the morphological file as a string. Required.
    
    """
    
    ##############################################
    ## READ FILE AND SUBSET TO RELEVANT COLUMNS ##
    ##############################################

    # Read morphology file 
    df_morph = pd.read_csv('./test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv', header = 0)

    # List relevant column names
    relevant_columns = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']

    # The input file must absolutely have these columns, no exceptions 
    if all(col in df_morph.columns for col in relevant_columns) == False:
        sys.exit(print("The input file", mfiles, "must have the columns:", ', '.join(relevant_columns)))

    # Keep only relevant columns
    df_morph = df_morph.loc[:,relevant_columns]

    ###########################################
    ## CALCULATE VARIABLES FOR DOSE RESPONSE ##
    ###########################################

    # Create groups of each chemical id, concentration, and plate id
    plate_groups = df_morph.drop(["well"], 1).groupby(by = ["chemical.id", "conc", "plate.id", "endpoint"], as_index = False)

    # Get the number of samples per group
    num_tot_samples = plate_groups.size().rename(columns = {"size": "num.tot"})

    # Get the number of non-na samples per groups
    num_nonna = plate_groups.count().rename(columns = {"value": "num.nonna"})

    # Get the number affected
    num_affected = plate_groups.sum().rename(columns = {"value": "num.affected"})

    # Merge to create missingness dataframe
    plate_groups = pd.merge(pd.merge(num_tot_samples, num_nonna), num_affected)

    # Create IDs of chemical.id, plate.id, and endpoint in plate_groups 
    ids = []
    for row in range(len(plate_groups)):
        ids.append(str(plate_groups["chemical.id"][row]) + " " + str(plate_groups["plate.id"][row]) + " " + str(plate_groups["endpoint"][row]))
    plate_groups["ids"] = ids

    #####################################################################
    ## REMOVE VARIABLES WITH HIGH MISSINGNESS IN BASELINE MEASUREMENTS ##
    #####################################################################

    # Identify 0 (baseline) concentrations with high missingness (greater than 50% missing or less than 50% non-missing)
    missingness = plate_groups.loc[plate_groups["conc"] == 0]
    missingness["keep"] = missingness["num.nonna"] / missingness["num.tot"] > 0.5 # TODO: Add a report of what was removed --> txt file "nothing removed"

    # Identify plates to keep 
    tokeep = missingness.loc[missingness["keep"]]["ids"].tolist()
    plate_groups = plate_groups[plate_groups["ids"].isin(tokeep)]

    # Stop if everything gets removed
    if len(plate_groups) == 0:
        sys.exit("Everything was removed with the 50% missingness filter")

    #######################################
    ## REGROUP WITHOUT PLATE IDS AND SUM ##
    #######################################

    # First, remove plate.id and ids column
    chemical_groups = plate_groups.drop(columns = ["plate.id", "ids"])

    # Group by chemical.id, concentration, and endpoint. Then, sum the results. 
    chemical_groups = chemical_groups.groupby(by = ["chemical.id", "conc", "endpoint"]).sum().reset_index()

    ##################################
    ## SUBSET TO RELEVANT ENDPOINTS ##
    ##################################

    # List the relevant endpoints, which is different for BRAIN samples  
    if "BRAI" in list(chemical_groups["endpoint"].unique()):
        relevant_endpoints = ['AXIS', 'BRAI', 'CFIN', 'CIRC', 'DNC_', 'DP24', 'EYE_', 'JAW_', 'MO24', 
                              'MORT', 'NC24', 'NC__', 'OTIC', 'PE__', 'PFIN', 'PIG_', 'SM24', 'SNOU', 
                              'SOMI', 'SWIM', 'TRUN', 'TR__', 'YSE_']
    else:
        relevant_endpoints = ['AXIS', 'BRN_', 'CRAN', 'DNC_', 'DP24', 'EDEM', 'LTRK', 'MO24', 'MORT', 
                              'MUSC', 'NC__', 'SKIN','SM24', 'TCHR']

    # Subset down to the relevant endpoints 
    chemical_groups = chemical_groups[chemical_groups["endpoint"].isin(relevant_endpoints)]

    ###########################
    ## ADD MISSING ENDPOINTS ##
    ###########################

    def new_endpoint(endpoints, new_name):
        """
        Generate a new endpoint which is a sum of existing endpoints.

        Attributes
        ----
        endpoints: A list of column names, as strings, to sum.
        new_name: The name of the new endpoint. 

        """
        sub_df = chemical_groups[chemical_groups["endpoint"].isin(endpoints)]
        sub_df["endpoint"] = new_name
        sub_df = sub_df.groupby(by = ["chemical.id", "conc", "endpoint"]).sum().reset_index()
        return(sub_df)

    # New endpoints to add is a smaller list if the sample is not from BRAIN
    if "BRAI" in list(chemical_groups["endpoint"].unique()):

        chemical_groups = pd.concat(
            [new_endpoint(['MO24','DP24','SM24','NC24'], 'ANY24'),
             new_endpoint(['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', 'PE__', 'BRAI', 
                          'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC', 'TRUN', 'SWIM', 'NC__', 'TR__', 
                          'ANY24'], 'ANY120'),
             new_endpoint(['MO24','MORT'], 'TOT_MORT'),
             new_endpoint(['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', 'PE__', 
                          'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__', 'TR__'], 'ALL_BUT_MORT'),
             new_endpoint(['BRAI','OTIC','PFIN'], 'BRN_'),
             new_endpoint(['EYE_', 'SNOU', 'JAW_'], 'CRAN'),
             new_endpoint(['YSE_','PE__'], 'EDEM'),
             new_endpoint(['TRUN','CFIN'], 'LTRK'),
             new_endpoint(['CIRC','SWIM','SOMI'], 'MUSC'),
             new_endpoint(['PIG_'], 'SKIN'),
             new_endpoint(['TR__'], 'TCHR'),
             chemical_groups]
        )

    else:

        chemical_groups = pd.concat(

            # 1. Add any effect at 24hrs (combination of MO24, DP24 and SM24) 
            [new_endpoint(['MO24','DP24','SM24'], 'ANY24'),

            # 2. Any effect within 5 days (combination of all measurements at both time points)
            new_endpoint(['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24'], 'ANY120'),

            # 3. Total mortality (MO24 + MORT) 
            new_endpoint(['MO24','MORT'], 'TOT_MORT'),

            # 4. Any effect except mortality (#2 minus MO24 and MORT)
            new_endpoint(['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', 'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR'], 'ALL_BUT_MORT'),

            # Add original dataframe
            chemical_groups]
        )

    ############################
    ## RETURN FORMATTED TABLE ##
    ############################
    return(chemical_groups)

def main():
    args = sys.argv[0:]
    mfile = args[1]
    format(mfile)
    
if __name__ == "__main__":
    main()
