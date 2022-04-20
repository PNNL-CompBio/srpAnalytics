#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os, random, sys, time

def format(mfile):
    """
    Pivot the morphological file to a wide format. Returns a dataframe.
    Additional endpoints are calculated, and total endpoints vary depending on 
    whether the inputted data is brain or not. 
    
    Attributes
    ----
    mfile: Path to the morphological file as a string. Required.
    
    """
    
    ###################################
    ## PIVOT WIDER AND CLEAN RESULTS ##
    ###################################

    # Read morphology file 
    df_morph = pd.read_csv('./test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv', header = 0)

    # List relevant column names
    relevant_columns = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']

    # The input file must absolutely have these columns, no exceptions 
    if all(col in df_morph.columns for col in relevant_columns) == False:
        sys.exit(print("The input file", mfiles, "must have the columns:", ', '.join(relevant_columns)))

    # Keep only relevant columns
    df_morph = df_morph.loc[:,columns_to_keep]

    # Pivot long format table wider, and sum up the ones 
    pivot_wider = df_morph.pivot_table(index = ['chemical.id', 'plate.id', 'well', 'conc'], columns = ['endpoint'], 
                                       values = ['value'], aggfunc = 'first')

    # Convert to a clean pandas dataframe 
    pivot_wider = pd.DataFrame(pivot_wider.to_records())

    # Clean up the column names 
    new_colnames = []

    # Fix column names with the extra 'value' annotation. Otherwise, ignore.
    for name in pivot_wider.columns:
        if "value" in name:
            new_colnames.append(re.sub("\\(value,|\\)", "", re.sub("'", "", name)).strip())
        else:
            new_colnames.append(name)

    # Rename columns
    pivot_wider = pivot_wider.set_axis(new_colnames, axis = 1, inplace = False)

    # Get the relevant endpoints, which is different for brain samples 
    if "BRAI" in pivot_wider.columns:
        relevant_endpoints = ['AXIS', 'BRAI', 'CFIN', 'CIRC', 'DNC_', 'DP24', 'EYE_', 'JAW_', 'MO24', 
                              'MORT', 'NC24', 'NC__', 'OTIC', 'PE__', 'PFIN', 'PIG_', 'SM24', 'SNOU', 
                              'SOMI', 'SWIM', 'TRUN', 'TR__', 'YSE_']
    else:
        relevant_endpoints = ['AXIS', 'BRN_', 'CRAN', 'DNC_', 'DP24', 'EDEM', 'LTRK', 'MO24', 'MORT', 
                              'MUSC', 'NC__', 'SKIN','SM24', 'TCHR']

    # Add relevant endpoints to the list of relevant columns, with value and endpints removed
    relevant_columns.remove("endpoint")
    relevant_columns.remove("value")
    for x in relevant_endpoints:
        relevant_columns.append(x)

    # Subset by those relevant columns 
    pivot_wider = pivot_wider.loc[:,relevant_columns]

    ###########################
    ## ADD MISSING ENDPOINTS ##
    ###########################

    if 'BRAI' not in pivot_wider.columns: 

        # 1. Add any effect at 24hrs (combination of MO24, DP24 and SM24) 
        pivot_wider['ANY24'] = pivot_wider[['MO24','DP24','SM24']].sum(axis=1,skipna=True,min_count=1)

        # 2. Any effect within 5 days (combination of all measurements at both time points)
        pivot_wider['ANY120'] = pivot_wider[['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 
                                             'NC__', 'SKIN', 'TCHR', 'ANY24']].sum(axis=1,skipna=True,min_count=1)

        # 3. Total mortality (MO24 + MORT) 
        pivot_wider['TOT_MORT'] = pivot_wider[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)

        # 4. Any effect except mortality (#2 minus MO24 and MORT)
        pivot_wider['ALL_BUT_MORT'] = pivot_wider[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', 'LTRK', 'MUSC', 'NC__', 
                                             'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)

    else:
        pivot_wider['ANY24'] = pivot_wider[['MO24','DP24','SM24','NC24']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['ANY120'] = pivot_wider[['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', 
                                       'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC', 
                                       'TRUN', 'SWIM', 'NC__', 'TR__', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['TOT_MORT'] = pivot_wider[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['ALL_BUT_MORT'] = pivot_wider[['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', 
                                             'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__', 
                                             'TR__']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['BRN_'] = pivot_wider[['BRAI','OTIC','PFIN']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['CRAN'] = pivot_wider[['EYE_', 'SNOU', 'JAW_']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['EDEM'] = pivot_wider[['YSE_','PE__']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['LTRK'] = pivot_wider[['TRUN','CFIN']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['MUSC'] = pivot_wider[['CIRC','SWIM','SOMI']].sum(axis=1,skipna=True,min_count=1)
        pivot_wider['SKIN'] = pivot_wider[['PIG_']]
        pivot_wider['TCHR'] = pivot_wider[['TR__']]


    return(pivot_wider)

def main():
    args = sys.argv[0:]
    mfile = args[1]
    format(mfile)
    
if __name__ == "__main__":
    main()
