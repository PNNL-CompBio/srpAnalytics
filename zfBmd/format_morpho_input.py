#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os, random, sys, time

def format(mfile):
    """
    Pivot the morphological file to a wide format. Returns the dataframe. 
    
    Attributes
    ----
    mfile: Path to the morphological file as a string. Required.
    
    """

    # Read morphology file 
    df_morph = pd.read_csv(complete_file_path, header = 0)
    
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
        
    # Add relevant endpoints to the list of relevant columns 
    for x in relevant_endpoints:
        relevant_columns.append(x)
            
    # Subset by those relevant columns 
    pivot_wider = pivot_wider.loc[:,relevant_columns]
    
    # Rename the earlier columns
    return(pivot_wider)