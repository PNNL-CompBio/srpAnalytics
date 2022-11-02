#!/usr/bin/env python
# coding: utf-8
    
import pandas as pd
import numpy as np
import sys
from icecream import ic

def pre_launch_cleaning(thePaths, theDataType):
    """
    Compiles csvs and removes duplicates from datasets. 

    If theData is "morphology", the dataframe is simply uniqued.
    If theData is "lpr", NA measurements are converted to 0, NA concentrations are removed, 
        and an average is taken over duplicates. 
    """
    # Create a vector to hold data
    theData = []

    # Read all files 
    for thePath in thePaths:
        theData.append(pd.read_csv(thePath))

    # Combine dataframes 
    theData = pd.concat(theData)

    # Now, do filetype specific cleaning
    if theDataType == "morphology":

        # Set the plate.id to a string for correct uniquing
        theData["plate.id"] = theData["plate.id"].astype("string") 
        
        # Remove duplicates in the data
        theData = theData.drop_duplicates()

        # Return results
        return(theData)

    if theDataType == "lpr":

        # Convert NA measurements to 0
        theData["value"] = theData["value"].fillna(0)

        # Remove NA concentrations
        theData = theData.dropna(subset = ["conc"])

        # Average values in duplicates
        theData = theData.groupby(by = ["chemical.id", "conc", "plate.id", "well", "variable"]).agg({'value':'mean'}).reset_index()

        # Return results
        return(theData)

def format_morpho_input(morpho_data):
    """
    Formats morphology file for the zfBMD pipeline. This involves:
        1. checking input columns and endpoints
        2. setting endpoints for fish that have experienced mortality or in the 'do not count' column as NA
        3. add missing endpoints with a bit-wise 'or'
        4. calculating the number of samples, non-na samples, and number effected
    
    Returns
    ----
    1. a DataFrame ready for dose response calculations and containing the following columns: 
    chemical.id, conc, endpoint, num.tot, num.nonna, num.affected, id, frac.affected

    2. all unique endpoints found in this file

    3. a list of wells that experience mortality at 5 days

    4. a list of wells that experience mortality at 24 hours
    
    Parameters
    ----
    morph_path: a string indicating the path to the morphology file
    """

    ##############################################
    ## READ FILE AND SUBSET TO RELEVANT COLUMNS ##
    ##############################################

    # Read morphology file 
    df_morph = morpho_data.reset_index()

    # List relevant column names
    relevant_columns = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']

    # The input file must absolutely have these columns, no exceptions 
    if all(col in df_morph.columns for col in relevant_columns) == False:
        sys.exit(print("The input morphology file must have the columns:", ', '.join(relevant_columns)))

    # Keep only relevant columns
    df_morph = df_morph.loc[:,relevant_columns]

    ##################################
    ## SUBSET TO RELEVANT ENDPOINTS ##
    ##################################

    # List the relevant endpoints, which is different for BRAIN samples  
    if "BRAI" in list(df_morph["endpoint"].unique()):
        relevant_endpoints = ['AXIS', 'BRAI', 'BRN_', 'CRAN', 'DNC_', 'CFIN', 'CIRC', 'DNC_', 
                              'DP24', 'EYE_', 'EDEM', 'JAW_', 'LTRK', 'MO24', 'MORT', 'MUSC', 
                              'NC24', 'NC__', 'OTIC', 'PE__', 'PFIN', 'PIG_', 'SKIN', 
                              'SM24', 'SNOU', 'SOMI', 'SWIM', 'TCHR', 'TRUN', 'TR__', 'YSE_']
    else:
        relevant_endpoints = ['AXIS', 'BRN_', 'CRAN', 'DNC_', 'DP24', 'EDEM', 'LTRK', 'MO24', 'MORT', 
                            'MUSC', 'NC__', 'SKIN','SM24', 'TCHR']

    # Get all endpoints 
    theEndpoints = df_morph['endpoint'].unique().tolist()

    # Get endpoints that are not expected 
    unexpected = [end for end in theEndpoints if end not in relevant_endpoints]

    # Print a message for missing endpoints 
    if (len(unexpected) != 0):
        print("The following endpoints were discovered and are unexpected:", unexpected)

    # Subset down to the relevant endpoints 
    df_morph = df_morph[df_morph["endpoint"].isin(relevant_endpoints)]

    #################################
    ## CONVERT DO NOT COUNTS TO NA ##
    #################################

    # Generate well id for subsetting endpoints
    df_morph["well.id"] = df_morph["chemical.id"].astype(str) + " " + df_morph["conc"].astype(str) + " " + df_morph["plate.id"].astype(str) + " " + df_morph["well"].astype(str)

    # If there's a do not count category, remove data 
    if ("DNC_" in theEndpoints):

        # Let the user know that data have been removed
        print("A do not count category was detected and those wells were changed to NA.")

        # Wells to remove
        toRmWells = df_morph[(df_morph["endpoint"] == "DNC_") & (df_morph["value"] == 1)]["well.id"]

        # Make the wells NA
        df_morph[df_morph["well.id"].isin(toRmWells)]["value"] = np.nan

        # And now remove the DNC category entirely
        df_morph = df_morph[df_morph["endpoint"] != "DNC_"]

    #############################
    ## CONVERT MORTALITY TO NA ##
    #############################

    # Convert endpoints affected by nor
    if "MORT" in theEndpoints:

        # Let users know that the affected organisms were changed to NA
        print("Mortality data at 5 days was detected, and affected endpoints in wells were changed to NA.")

        # Pull out the wells with mortatility 
        MortWells = df_morph[(df_morph["endpoint"] == "MORT") & (df_morph["value"] == 1)]["well.id"]

        # Determine non-24 hour endpoints without MORT
        non24 = [end for end in relevant_endpoints if end not in ["DP24", "MO24", "SM24", "MORT"]]

        # Change affected values to NA
        df_morph[df_morph["well.id"].isin(MortWells) & df_morph["endpoint"].isin(non24)]["value"] = np.nan

    # Remove mortality at 24 hour
    if "MO24" in theEndpoints:

        # Let users know that the affected organisms were changed to NA
        print("Mortality data at 24 hours was detected and affected endpoints in wells were changed to NA.")

        # Pull out the wells with mortality at 24 hours 
        Mort24Wells = df_morph[(df_morph["endpoint"] == "MO24") & (df_morph["value"] == 1)]["well.id"]

        # All but the MO24 endpoints should be NA when MO24 == 1
        nonMO24 = [end for end in relevant_endpoints if end != "MO24"]

        # Change affected values to NA
        df_morph[df_morph["well.id"].isin(Mort24Wells) & df_morph["endpoint"].isin(nonMO24)]["value"] = np.nan

    ###########################
    ## ADD MISSING ENDPOINTS ##
    ###########################

    def new_endpoint(endpoints, new_name):
        """
        Generate a new endpoint which is a binary "or" of other endpoints,
        meaning that if there is a 1 in any of the other endpoints, the 
        resulting endpoint is a 1. Otherwise, it is 0 unless the other 
        endpoints are all NA. Then the final value is NA.
        
        Parameters
        ----
        endpoints: A list of column names, as strings, to binary "or". 
        new_name: The name of the new endpoint. 
        
        """
        sub_df = df_morph[df_morph["endpoint"].isin(endpoints)]
        sub_df["endpoint"] = new_name
        sub_df = sub_df.groupby(by = ["chemical.id", "conc", "plate.id", "well", "endpoint"], as_index = False).sum()
        sub_df['value'].values[sub_df['value'] > 1] = 1 
        return(sub_df)
        
    # New endpoints to add is a smaller list if the sample is not from BRAIN
    if "BRAI" in list(df_morph["endpoint"].unique()):
        
        df_morph = pd.concat(
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
            df_morph]
        )
        
    else:
        
        df_morph = pd.concat(

            # 1. Add any effect at 24hrs (combination of MO24, DP24 and SM24) 
            [new_endpoint(['MO24','DP24','SM24'], 'ANY24'),

            # 2. Any effect within 5 days (combination of all measurements at both time points)
            new_endpoint(['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24'], 'ANY120'),

            # 3. Total mortality (MO24 + MORT) 
            new_endpoint(['MO24','MORT'], 'TOT_MORT'),

            # 4. Any effect except mortality (#2 minus MO24 and MORT)
            new_endpoint(['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', 'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR'], 'ALL_BUT_MORT'),
            
            # Add original dataframe
            df_morph]
        )
        
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
    plate_groups["ids"] = plate_groups["chemical.id"].astype(str) + " " + plate_groups["plate.id"].astype(str) + " " + plate_groups["endpoint"].astype(str)

    #########################################################################
    ## REMOVE VARIABLES WITH HIGH MEASURED EFFECT IN BASELINE MEASUREMENTS ##
    #########################################################################

    # Identify 0 (baseline) concentrations with high measured effect
    missingness = plate_groups.loc[plate_groups["conc"] == 0]

    # Remove if more than 50% is 1's 
    missingness["keep"] = (missingness["num.nonna"] / missingness["num.tot"]) > 0.5

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
    chemical_groups = plate_groups.drop(columns = ["plate.id", "ids", "well.id"])

    # Group by chemical.id, concentration, and endpoint. Then, sum the results. 
    dose_response = chemical_groups.groupby(by = ["chemical.id", "conc", "endpoint"]).sum().reset_index()

    ############################
    ## RETURN FORMATTED TABLE ##
    ############################
    return([dose_response, theEndpoints, MortWells, Mort24Wells])

def format_lpr_input(lpr_data, theEndpoints, MortWells, Mort24Wells):
    """
    Formats an LPR file for the zfBMD pipeline. This involves:
        1. removing endpoints for fish that have experienced mortality
        2. summing counts to the minute
        3. removing plates with dominant outliers at concentration 0
        4. converting the continuous response variable to a dichotomous variable
    
    Returns
    ----
    a DataFrame ready for dose response calculations and containing the following columns: 
    chemical.id, conc, endpoint, num.tot, num.nonna, num.affected, id, frac.affected
    
    Parameters
    ----
    lpr_path: a string indicating the path to the lpr file 
    theEndpoints: all unique endpoints in the morphology file
    MortWells: a list of wells with mortality from format_morpho_input 
    Mort24Wells: a list of wells with mortality at 24 hours from format_morpho_input
    """

    ##########################################
    ## READ FILE, KEEP ALL RELEVANT COLUMNS ##
    ##########################################

    # Read LPR data
    df_LPR = lpr_data.reset_index()
    df_LPR_data = lpr_data.reset_index()

    # List relevant column names
    relevant_LPR_columns = ['chemical.id', 'conc', 'plate.id', 'well', 'variable', 'value']

    # The input file must absolutely have these columns, no exceptions 
    if all(col in df_LPR.columns for col in relevant_LPR_columns) == False:
        sys.exit(print("The input LPR file must have the columns:", ', '.join(relevant_LPR_columns)))

    # Determine if the LPR data has been taken in 6 second intervals, or 1 minute intervals 
    maxTime = max([int(x.replace("t", "")) for x in df_LPR["variable"].unique().tolist()])
    if (maxTime == 239):
        
        # Let the user know what timepoints are being calculated
        print("Max dose time is 240. The assumed time interval is 6 seconds.")

        # Add time 
        df_LPR["timepoint"] = np.floor((pd.to_numeric(df_LPR["variable"].replace("t", "", regex = True))) / 10)

    else:
        sys.exit(print("The max time interval of", maxTime, "is not known."))

    ###################################################
    ## CONFIRM FISH THAT EXPERIENCE MORTALITY ARE NA ##
    ###################################################

    ic()

    # Make an ID column
    df_LPR["well.id"] = df_LPR["chemical.id"].astype(str) + " " + df_LPR["conc"].astype(str) + " " + df_LPR["plate.id"].astype(str) + " " + df_LPR["well"].astype(str)

    ic()

    # Remove wells that experience mortality
    if "MORT" in theEndpoints:
        df_LPR[df_LPR["well.id"].isin(MortWells)]["value"] = np.nan
    if "MO24" in theEndpoints:
        df_LPR[df_LPR["well.id"].isin(Mort24Wells)]["value"] = np.nan

    ####################################
    ## SUM INTENSITIES AND GET COUNTS ##
    ####################################

    ic()

    # Group LPR data by chemical id, concentration, plate id, well and timepoint
    df_LPR_grouped = df_LPR.groupby(by = ["chemical.id", "conc", "plate.id", "well", "timepoint"], as_index = False)

    # Get the sum and NA values
    df_LPR_sum = df_LPR_grouped.sum().rename(columns = {"value":"Sum"})
    df_LPR_NA = df_LPR_grouped.apply(lambda df: any(np.isnan(df["value"]))).rename(columns = {None:"Remove"})["Remove"]

    # Remove the NA values 
    df_LPR_sum["Remove"] = df_LPR_NA
    df_LPR_sum = df_LPR_sum[df_LPR_sum["Remove"] == False].drop(columns = ["Remove"])

    #######################################
    ## ADD ENDPOINT NUMERIC CALCULATIONS ##
    #######################################

    ic()

    # Group data frames for endpoint calculations 
    df_LPR_sum_grouped = df_LPR_sum.groupby(by = ["chemical.id", "conc", "plate.id", "well"])
    LPR_Endpoints = df_LPR_sum_grouped.apply(lambda df: df[df["timepoint"] == 3]["Sum"]).reset_index().rename(columns = {"Sum":"Time3"}).drop(columns = "level_4")
    LPR_Endpoints["Time2"] = df_LPR_sum_grouped.apply(lambda df: df[df["timepoint"] == 2]["Sum"]).reset_index()["Sum"]
    LPR_Endpoints["MOV1"] = LPR_Endpoints["Time3"] - LPR_Endpoints["Time2"]
    LPR_Endpoints = LPR_Endpoints.drop(["Time2", "Time3"], axis = 1)

    # Define a function to sum values
    def sum_endpoints(num_list):
        return(df_LPR_sum_grouped.apply(lambda df: df[df["timepoint"].isin(num_list)]["Sum"]).reset_index()[["Sum"]])

    # Add MOV2 - MOV4
    LPR_Endpoints["MOV2"] = sum_endpoints([9]) - sum_endpoints([8])
    LPR_Endpoints["MOV3"] = sum_endpoints([15]) - sum_endpoints([14])
    LPR_Endpoints["MOV4"] = sum_endpoints([21]) - sum_endpoints([20])

    # Add AUC1 - AUC4
    LPR_Endpoints["AUC1"] = sum_endpoints([3,4,5]) - sum_endpoints([0,1,2])
    LPR_Endpoints["AUC2"] = sum_endpoints([9,10,11]) - sum_endpoints([6,7,8])
    LPR_Endpoints["AUC3"] = sum_endpoints([15,16,17]) - sum_endpoints([12,13,14])
    LPR_Endpoints["AUC4"] = sum_endpoints([21,22,23]) - sum_endpoints([18,19,20])

    # Create a function to transform the numeric to dichotomous
    def to_dichotomous(toCalculate, endpoint):
        '''

        A function that converts the continuous LPR response to dichotomous

        Parameters
        ----
        toCalculate: a pandas dataframe with chemical.id, concentration, plate.id, and the endpoint
        endpoint: the endpoint as a string
        '''

        # Pull plate and chemical groups concentration at 0
        plateGroups = toCalculate[toCalculate["conc"] == 0].groupby(["chemical.id", "plate.id"])

        # Pull quartile calculations
        rangeValues = plateGroups.apply(lambda df: df[endpoint].quantile(0.25)).reset_index().rename(columns = {0:"Q1"})
        rangeValues["Q3"] = plateGroups.apply(lambda df: df[endpoint].quantile(0.75)).reset_index()[0]

        # Add IQR and lower and upper bonds. 
        rangeValues["IQR"] = rangeValues["Q3"] - rangeValues["Q1"]
        rangeValues["Low"] = rangeValues["Q1"] - (1.5 * rangeValues["IQR"])
        rangeValues["High"] = rangeValues["Q3"] + (1.5 * rangeValues["IQR"])

        # Add range values
        toCalculate = toCalculate.merge(rangeValues[["chemical.id", "plate.id", "Low", "High"]])

        # The number affected is the number of hypoactive (negative) or hyperactive (outside of 1.5 * IQR)
        toCalculate["num.affected"] = (toCalculate[endpoint] < 0) | (toCalculate[endpoint] < toCalculate["Low"]) | (toCalculate[endpoint] > toCalculate["High"])

        # Group by chemical id and concentration
        preEndpoint = toCalculate.drop(columns = [endpoint, "Low", "High", "plate.id"]).groupby(["chemical.id", "conc"])

        # Caculate the endpoint and add other columns
        Endpoint = preEndpoint.sum().reset_index()
        Endpoint["num.nonna"] = preEndpoint.size().reset_index()[0]
        Endpoint["endpoint"] = endpoint
        Endpoint["num.tot"] = df_LPR_data[df_LPR_data["variable"] == "t0"].groupby(["chemical.id", "conc"]).size().reset_index()[0]

        # Order and return endpoint
        Endpoint = Endpoint[["chemical.id", "conc", "endpoint", "num.tot", "num.nonna", "num.affected"]]
        return(Endpoint)

    ic()

    dose_response = pd.concat([
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "MOV1"]], "MOV1"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "MOV2"]], "MOV2"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "MOV3"]], "MOV3"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "MOV4"]], "MOV4"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "AUC1"]], "AUC1"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "AUC2"]], "AUC2"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "AUC3"]], "AUC3"),
        to_dichotomous(LPR_Endpoints[["chemical.id", "conc", "plate.id", "AUC4"]], "AUC4")
    ])

    return(dose_response)