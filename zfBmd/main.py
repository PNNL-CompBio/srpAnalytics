#!/usr/bin/env python
# coding: utf-8

######################
## IMPORT LIBRARIES ##
######################

# Import python libraries 
import pandas as pd
import sys
import argparse

# Import zfBMD specific functions 
from format_binary_data import format_morpho_input
from format_binary_data import format_lpr_input 
from calculate_BMD_flags import generate_BMD_flags
from select_and_run_models import model_fitting
from select_and_run_models import export_BMDs
from select_and_run_models import export_fits
from select_and_run_models import export_doses

###########################
## COLLECT CLI ARGUMENTS ##
###########################

parser = argparse.ArgumentParser('Run the QC and BMD analysis for the SRP analytics compendium')

parser.add_argument('--morpho', dest='morpho',\
                    help='Pathway to the morphological file to be processed. Required.',\
                    default=None)
parser.add_argument('--LPR', dest='lpr', \
                    help='Pathway to the light photometer response (LPR) file to be processed containing the same \
                          samples as morpho. Optional. Unless both is True, only LPR data will be returned.',\
                    default=None)
parser.add_argument('--both', dest='both', \
                    help='Return both morpho and LPR endpoints. Optional. Default is False.',\
                    default=False)
parser.add_argument('--output', dest = 'output', \
                    help = 'The output folder for files. Default is current directory.',\
                    default = '.')
parser.add_argument('--test', dest='test',\
                    help='Set this flag to test code with internal files. Default is False.',\
                    action='store_true', default=False)
                    
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
    morpho_path = args.morpho
    lpr_path = args.lpr

    # If there is no morphological data and this is not the test mode, stop. 
    if args.morpho is None and args.test == False:
        sys.exit("--morpho cannot be blank, since morphological data is required to run zfBMD.")
    
    # Users must supply both morpho and lpr data if both is selected
    if args.both and args.lpr is None and args.test == False:
        sys.exit("--lpr cannot be blank if args.both is True.")

    # Load test data if test is true 
    if args.test == True:
        morpho_path = './test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv'
        lpr_path = './test_files/7_PAH_zf_LPR_data_2021JAN11_3756.csv'
    
    ### 1. Format data--------------------------------------------------------------------------
    print("...Formatting morphology data")
    #dose_response, theEndpoints, MortWells, Mort24Wells = format_morpho_input(morpho_path)

    #dose_response.to_csv("./Results/dose_response.csv", index = False)
    #pd.DataFrame(theEndpoints, columns = ["endpoints"]).to_csv("./Results/theEndpoints.csv", index = False)
    #pd.DataFrame(MortWells, columns = ["theWells"]).to_csv("./Results/MortWells.csv", index = False)
    #pd.DataFrame(Mort24Wells, columns = ["theWells"]).to_csv("./Results/Mort24Wells.csv", index = False)

    dose_response = pd.read_csv("./Results/dose_response.csv")
    theEndpoints = pd.read_csv("./Results/theEndpoints.csv")["endpoints"].to_list()
    MortWells = pd.read_csv("./Results/MortWells.csv")["theWells"].to_list()
    Mort24Wells = pd.read_csv("./Results/Mort24Wells.csv")["theWells"].to_list()

    if (args.lpr is not None or (args.test and args.both)):
        print("...Formatting LPR data")
        #lpr_dose_response = format_lpr_input(lpr_path, theEndpoints, MortWells, Mort24Wells)
        #lpr_dose_response.to_csv("./Results/lpr_dose_response.csv", index = False)
        lpr_dose_response = pd.read_csv("./Results/lpr_dose_response.csv")

    ### 2. Calculate dose response--------------------------------------------------------------

    if (args.lpr is None or args.both):
        print("...Generating morphology flags")
        BMD_Flags = generate_BMD_flags(dose_response)
        BMD_Flags = pd.read_csv("./Results/BMD_Flags.csv")

    if (args.lpr is not None or (args.test and args.both)):
        print("...Generating LPR flags")
        lpr_BMD_Flags = generate_BMD_flags(lpr_dose_response)

    ### 3. Select and run models----------------------------------------------------------------

    if (args.lpr is None or args.both):
        print("...Fitting models for morphology data")
        model_selection, lowqual_model, BMD_Flags, model_results = model_fitting(dose_response, BMD_Flags)
    
    if args.lpr is not None or (args.test and args.both):
        print("...Fitting models for LPR data")
        lpr_model_selection, lpr_lowqual_model, lpr_BMD_Flags, lpr_model_results = model_fitting(lpr_dose_response, lpr_BMD_Flags)

    ### 4. Format and export outputs------------------------------------------------------------

    print("...Exporting Results")

    if args.lpr is None and args.both == False:

        # Benchmark Dose #
        BMDS_Final = export_BMDs(dose_response, BMD_Flags, model_selection, lowqual_model)
        BMDS_Final_Clean = BMDS_Final.drop("ids", axis = 1)
        BMDS_Final_Clean.to_csv(args.output + "/new_bmds.csv", index = False)

        # Fits #
        export_fits(model_results, dose_response, BMDS_Final).to_csv(args.output + "/new_fits.csv", index = False)

        # Doses # 
        export_doses(dose_response).to_csv(args.output + "/new_dose.csv", index = False)

    elif args.lpr is not None and args.both == False:

        # Benchmark Dose #
        lpr_BMDS_Final = export_BMDs(lpr_dose_response, lpr_BMD_Flags, lpr_model_selection, lpr_lowqual_model)
        lpr_BMDS_Final.to_csv(args.output + "/new_bmds.csv", index = False)

        # Fits #
        export_fits(lpr_model_results, lpr_dose_response, lpr_BMDS_Final).to_csv(args.output + "/new_fits.csv", index = False)
        
        # Doses #
        export_doses(lpr_dose_response).to_csv(args.output + "/new_dose.csv", index = False)

    elif args.both:

        # Benchmark Dose #
        BMDS_Final = export_BMDs(dose_response, BMD_Flags, model_selection, lowqual_model)
        lpr_BMDS_Final = export_BMDs(lpr_dose_response, lpr_BMD_Flags, lpr_model_selection, lpr_lowqual_model)
        pd.concat([BMDS_Final, lpr_BMDS_Final]).to_csv(args.output + "/new_bmds.csv", index = False)

        # Fits #
        pd.concat(
            [export_fits(model_results, dose_response, BMDS_Final),
             export_fits(lpr_model_results, lpr_dose_response, lpr_BMDS_Final)]
        ).to_csv(args.output + "/new_fits.csv", index = False)
        
        # Doses # 
        pd.concat(
            [export_doses(dose_response),
             export_doses(lpr_dose_response)]
        ).to_csv(args.output + "/new_dose.csv", index = False)

if __name__ == "__main__":
    main()
