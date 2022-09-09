#!/usr/bin/env python
# coding: utf-8

import scipy.stats as stats
import numpy as np

def generate_BMD_flags(dose_response):
    '''This function performs feasibility analysis
    for dose response data. The value returned is a
    flag indicating data quality as defined below:

    0: Not enough dose groups for BMD analysis. BMD analysis not performed. BAD.
    1: No trend detected in dose-response data. BMD analysis not performed. BAD. 
    2: Good dose-response data. BMD analysis is performed. GOOD.
    3: Dose-response data quality poor. BMD analysis might be unreliable. GOOD. 
    4: Data resolution poor. BMD analysis might be unreliable. MODERATE.
    5: Negative correlation. 

    '''

    ########################
    ## CALCULATE RESPONSE ##
    ########################

    # Add an id column
    dose_response["ids"] = dose_response["chemical.id"].astype(str) + " " + dose_response["endpoint"].astype(str)

    # Response is the number affected over the number of embryos (non-na)
    dose_response["frac.affected"] = dose_response["num.affected"] / dose_response["num.nonna"]

    ####################################
    ## GENERATE QUALITY CONTROL FLAGS ##
    ####################################

    # Count the number of unique concentrations per chemical id and endpoint pairing 
    BMD_Flags = dose_response.groupby(["chemical.id", "endpoint", "ids"])["conc"].nunique().reset_index().rename(columns = {"conc": "num.conc"})

    # Add a flag category
    BMD_Flags["flag"] = None

    # If there are less than 3 concentrations, the BMD flag is 0 - not enough dose groups
    BMD_Flags["flag"].values[BMD_Flags["num.conc"] < 3] = 0

    # Change dose response 
    dose_response["conc"] = dose_response["conc"].astype('float') + 1e-15

    # Calculate the spearman correlation
    spear = dose_response[["chemical.id", "endpoint", "conc", "frac.affected"]].groupby(["chemical.id", "endpoint"]).corr(method = "spearman").unstack().iloc[:,1].reset_index()
    spear = spear.set_axis(["chemical.id", "endpoint", "spearman"], axis = 1)

    # Merge spearman to the BMD_Flags dataframe
    BMD_Flags = BMD_Flags.merge(spear)

    # If spearman correlation is below 0.2 or is NaN, the Flag is 1 - no strong trend 
    BMD_Flags["flag"].values[(BMD_Flags["spearman"] < 0.2) | (BMD_Flags["spearman"].isna())] = 1

    # If the correlation is above 0.25 or below 0.8, run a t-test
    BMD_Flags["run.ttest"] = (BMD_Flags["spearman"] >= 0.20) & (BMD_Flags["spearman"] <= 0.80)

    # If the correlation is above 0.8, assign the flag at 2 - good
    BMD_Flags["flag"].values[BMD_Flags["spearman"] > 0.8] = 2

    # Run the t-test only where indicated 
    ttest = dose_response[dose_response["ids"].isin(BMD_Flags[BMD_Flags["run.ttest"]]["ids"].to_list())][["ids", "frac.affected"]]
    ttest = ttest.groupby("ids").apply(lambda df: stats.ttest_1samp(np.diff(df["frac.affected"]), 0)[1]).reset_index().rename(columns = {0:"ttest.pval"})

    # Merge ttest results 
    BMD_Flags = BMD_Flags.merge(ttest, on = "ids", how = "outer")

    # A p-value of less than 0.05 gets a flag of 2 - good, from 0.05 to 0.32 gets a 3 - unreliable, 
    # and greater than 0.32 gets very unreliable. 
    BMD_Flags["flag"].values[BMD_Flags["ttest.pval"] <= 0.05] = 2
    BMD_Flags["flag"].values[(BMD_Flags["ttest.pval"] > 0.05) & (BMD_Flags["ttest.pval"] <= 0.32)] = 3
    BMD_Flags["flag"].values[BMD_Flags["ttest.pval"] > 0.32] = 4

    return(BMD_Flags)