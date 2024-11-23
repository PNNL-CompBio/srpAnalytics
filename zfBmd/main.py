#!/usr/bin/env python
# coding: utf-8

######################
## IMPORT LIBRARIES ##
######################

# Import python libraries 
import pandas as pd
import sys
import argparse
from bmdrc.BinaryClass import BinaryClass
from bmdrc.LPRClass import LPRClass

# Import specific support_functions for this main pipeline function
from support_functions import combine_datasets, preprocess_morpho, preprocess_lpr, run_filters, write_outputs 

# Example commands

## morphology only: python3 main.py --morpho test_files/test_morphology.csv
## morphology & lpr: python3 main.py --morpho test_files/test_morphology.csv --LPR test_files/test_behavioral.csv --both True

###########################
## COLLECT CLI ARGUMENTS ##
###########################

parser = argparse.ArgumentParser('Run the QC and BMD analysis for the SRP analytics compendium')

parser.add_argument('--morpho', dest = 'morpho', nargs = "+", \
                    help = 'Pathway to the morphological file to be processed. Required. \
                            Assumed and needed column names are: chemical.id, conc, plate.id, well, variable, value. \
                            Data assumed to be in long format.',\
                    default = None)
parser.add_argument('--LPR', dest = 'lpr', nargs = "+", \
                    help = 'Pathway to the light photometer response (LPR) file to be processed containing the same \
                            samples as morpho. Optional. Unless both is True, only LPR data will be returned. \
                            Assumed format is long. Required columns are: chemical.id, conc, plate.id, well, variable, value.',\
                    default = None)
parser.add_argument('--both', dest = 'both', \
                    help = 'Return both morpho and LPR endpoints. Optional. Default is True.',\
                    default = True)
parser.add_argument('--output', dest = 'output', \
                    help = 'The output folder for files. Default is current directory.',\
                    default = '.')
parser.add_argument('--test', dest = 'test',\
                    help = 'Set this flag to test code with internal files. Default is False.',\
                    action = 'store_true', default=False)
parser.add_argument('--report', dest ='rep', \
                    help = 'Generate a report in the output folder. Default is True',
                    action = 'store_true', default = True)
                    
##############################
## DEFINE AND RUN FUNCTIONS ##
##############################

def main():
    """
    Gather the input arguments and run the zfBMD pipeline which: 

    1. formats input data
    2. calculate benchmark doses 
    3. select and run models

    Returns
    ----
    """

    # Parse inputted arguments from the command line 
    args = parser.parse_args()

    # Pull arguments
    morpho_paths = args.morpho
    lpr_paths = args.lpr

    # If there is no morphological data and this is not the test mode, stop. 
    if args.morpho is None and args.test == False:
        sys.exit("--morpho cannot be blank, since morphological data is required to run zfBMD.")
    
    # Users must supply both morpho and lpr data if both is selected
    if args.both and args.lpr is None and args.test == False:
        sys.exit("--lpr cannot be blank if args.both is True.")

    # Load test data if test is true 
    if args.test == True:
        morpho_paths = './test_files/test_morphology.csv'
        lpr_paths = './test_files/test_behavioral.csv'

    ### 0. Pre-launch combination-------------------------------------------------------------

    print("...Concatenating datasets")

    # Morpho data always needs to be provided. Concatenate datasets together. 
    morpho_data = combine_datasets(morpho_paths)

    # LPR data is optional. Concatenate datasets together. 
    if (args.lpr is not None or (args.test and args.both)):
        print("...LPR data detected")
        lpr_data = combine_datasets(lpr_paths)

    ### 1. Input Data Modules--------------------------------------------------------------------
    
    print("...Formatting morphology data")
    BC = BinaryClass(df = morpho_data, chemical = "chemical.id", concentration = "conc", \
                     plate = "plate.id", well = "well", endpoint = "endpoint", value = "value", \
                     format = "long")

    if (args.lpr is not None or (args.test and args.both)):
        print("...Formatting LPR data")
        LPR = LPRClass(df = lpr_data, chemical = "chemical.id", concentration = "conc", 
                       plate = "plate.id", well = "well", time = "variable", value = "value", 
                       cycle_length = 20.0, cycle_cooldown = 10.0, starting_cycle = "light")

    ### 2. Pre-Processing modules-----------------------------------------------------------------

    if (args.lpr is None or args.both):
        print("...Pre-Processing morphology data")
        preprocess_morpho(BC)

    if (args.lpr is not None or (args.test and args.both)):
        print("...Pre-Processing LPR data")
        preprocess_lpr(LPR, BC)

    ### 3. Filtering Modules----------------------------------------------------------------------

    if (args.lpr is None or args.both):
        print("...Filtering morphology data")
        run_filters(BC)

    if (args.lpr is not None or (args.test and args.both)):
        print("...Filtering LPR data")
        run_filters(LPR)

    ### 4. Model Fitting Modules------------------------------------------------------------------

    if (args.lpr is None or args.both):
        print("...Fitting models to morphology data")
        BC.fit_models(diagnostic_mode = True)

    if (args.lpr is not None or (args.test and args.both)):
        print("...Fitting models to LPR data")

        # Subset down to AUC2 and MOV2 and rename them "AUC" and "MOV"
        LPR.plate_groups = LPR.plate_groups[LPR.plate_groups[LPR.endpoint].isin(["AUC2", "MOV2"])]
        LPR.plate_groups[LPR.endpoint] = LPR.plate_groups[LPR.endpoint].str.replace('\d+', '', regex = True)
        LPR.plate_groups["bmdrc.Endpoint.ID"] = LPR.plate_groups[LPR.chemical].astype(str) + " " + LPR.plate_groups[LPR.endpoint].astype(str)
        LPR.fit_models(diagnostic_mode = True)

    ### 5. Format and export outputs------------------------------------------------------------

    print("...Exporting Results")

    if args.lpr is None and args.both == False:

        # Output files
        write_outputs(BC, "BC")

    elif args.lpr is not None and args.both == False:

        # Output files
        write_outputs(LPR, "LPR")

    elif args.both:

        # Output both files
        write_outputs(BC, "BC")
        write_outputs(LPR, "LPR")

if __name__ == "__main__":
    main()
