#!/usr/bin/env python
# coding: utf-8

######################
## IMPORT LIBRARIES ##
######################

# Import python libraries 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time
import argparse
import tarfile
import re

# Import zfBMD specific functions 
import bmd_analysis_morpho as bmd
import bmd_analysis_LPR_7_PAH_t0_t239 as bmd_LPR
import format_LPR_input as format_LPR
import format_morpho_input as format_morpho

###########################
## COLLECT CLI ARGUMENTS ##
###########################

parser = argparse.ArgumentParser('Run the QC and BMD analysis as well as join with \
extract data to store in SRP data analytics portal')

parser.add_argument('--morpho', dest='morpho',\
                    help='Comma-delimited list of morphological files to be processed. Required.',\
                    default=None)
parser.add_argument('--LPR', dest='lpr', \
                    help='Comma-delimited list of LPR-related files to be processed containing the same samples as morpho. Optional.',\
                    default=None)
parser.add_argument('--test', dest='test',\
                    help='Set this flag to test code with internal files.',\
                    action='store_true', default=False)


parser.add_argument('--test-extract', dest='test_extract',\
                    help='Set this flag to run morpho test code with extract data',\
                    action='store_true', default=False) # David: Do we need this? 

##############################
## DEFINE AND RUN FUNCTIONS ##
##############################

def run_lpr_on_file(lpr_file, morph_file, full_devel='full'):
    """
    Reformats inputted morphology (required) or lpr (optional) files and runs main zfBMD pipeline
    
    Attributes
    ----
    lpr_file: a string indicating the path to the lpr file (optional)
    morph_file: a string indicating the path to the morphology file (required) 
    full_devel: a variable to delete, most likely. Hard-coded to "full" 
    """
    
    # Generate a new name for the LPR file in a wider format 
    LPR_input_csv_file_name_wide = lpr_file[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"
    
    # Reformat LPR if necessary
    chem_ind = None
    if not os.path.exists(LPR_input_csv_file_name_wide):
       chem_ind = format_LPR.format(lpr_file, full_devel, LPR_input_csv_file_name_wide)
    
    # Generate a new name for the morphological file in wider format
    morpho_input_csv_file_name_wide = morph_file[:-4] + "_wide_DNC_0_"+full_devel+".csv"
    
    # Reformat Morphological file if necessary 
    if not os.path.exists(morpho_input_csv_file_name_wide):
        res0 = format_morpho.format(morph_file, full_devel, morpho_input_csv_file_name_wide, chem_ind)

    # Take reformatted files and run them through the main pipeline 
    res = bmd_LPR.runBmdPipeline(morpho_input_csv_file_name_wide, LPR_input_csv_file_name_wide, full_devel)
    return res

def run_morpho_on_file(morph_file, full_devel='full'):
    """
    Reformats inputted morpohology file - will be combined with above 
    """

    morpho_input_csv_file_name_wide = morph_file[:-4] + "_wide_DNC_0_"+full_devel+".csv"

    chem_ind = None
    if not os.path.exists(morpho_input_csv_file_name_wide):
        res0 = format_morpho.format(morph_file, full_devel, morpho_input_csv_file_name_wide, chem_ind)

    res = bmd.runBmdPipeline(morpho_input_csv_file_name_wide, full_devel)
    return res

def merge_data_and_write_files(file_dict, path = "/tmp"): # David: Move output path to command line arguments? 
    """
    Takes a dictionary of outputted data.frames, merges them, and returns three csv files.
    
    | File Name | Description, per each chemical ID and endpoint
    |-----------|-----------------------------------------------------------------------------------------------|
    | New BMDS  | The best fitting model is reported with AUC and BMD measurements                              |
    | New Fits  | The fitted values are reported with 7 predicted measurements between each actual measurement  |
    | New Dose  | Dose and reponses are written with confidence intervals                                       |

    Attributes
    ------
    file_dict: A python dictionary of data.frames to merge. The order is bmds, fits, and dose. Required.
    path : A string for the output folder of the three files. Defaults to "/tmp". 
    """

    # TODO: Change to just write from the dictionary. An extra list is not necessary. 
    bmds = []
    fits = []
    dose = []
    for dataset, filelist in file_dict.items():
        bmds.append(filelist[0])
        fits.append(filelist[1])
        dose.append(filelist[2])

    # Concatenate all files together 
    pd.concat([pd.read_csv(f) for f in bmds]).to_csv(path+'/new_bmds.csv') # David: What do we want the names to be? 
    pd.concat([pd.read_csv(f) for f in fits]).to_csv(path+'/new_fits.csv')
    pd.concat([pd.read_csv(f) for f in dose]).to_csv(path+'/new_dose.csv')
    return [path+'/new_bmds.csv', path+'/new_fits.csv', path+'/new_dose.csv']


def main():
    """
    Runs the zfBMD pipeline, where ...
    """

    # Start a timer for timestamping purposes 
    start_time = time.time()
    
    # Parse inputted arguments from the command line 
    args = parser.parse_args()

    # Initialize a dictionary to hold outputted files 
    files = dict()

    # If there is no morpholoigcal data, stop code and write warning. 
    if args.morpho is None and args.test == False:
        sys.exit("--morpho cannot be blank, since morphological data is required to run zfBMD.")
    else if args.test == False:
        mfiles = args.morpho.split(",")
    
    # LPR data is optional, so if it's included, split the arguments. 
    if args.lpr is not None and args.test == False:
        lfiles = args.lpr.split(',')

    fd = 'full' # TODO: Remove hard-coded variable 
    
    # If test mode has been activated, run with both morpho and lpr file in test_files. 
    if args.test:
        fd = 'devel' # TODO: Double check that this does anything 
        lfiles = ['/zfBmd/test_files/7_PAH_zf_LPR_data_2021JAN11_3756.csv']
        mfiles = ['/zfBmd/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv']
    

    if len(lfiles) > 0:
        if len(lfiles) != len(mfiles):
            sys.exit("Cannot calculate LPR without morphological files, please re-run with --morpho argument")
        else:
            print("Calculating LPR endpoints for ", str(len(lfiles)), " LPR files")
            for i in range(len(lfiles)):
                fname = lfiles[i]
                files[fname] = run_lpr_on_file(fname, mfiles[i], fd)
    elif len(mfiles) > 0:
        print("Calculating morphological endpoints for "+str(len(mfiles))+' files')
        for f in mfiles:
            files[f] = run_morpho_on_file(f, fd)

    
    # Merge each BMDS, Fits, and Dose data.frame together, and write to output files 
    merge_data_and_write_files(files, "/tmp")

    
    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print("zfBMD pipeline completed. It took:" + str(time_took))

if __name__ == "__main__":
    main()
