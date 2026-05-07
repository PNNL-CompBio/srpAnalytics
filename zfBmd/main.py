#!/usr/bin/env python
# coding: utf-8

######################
## IMPORT LIBRARIES ##
######################
import argparse
import os
import sys
from glob import glob

import pandas as pd
from bmdrc.BinaryClass import BinaryClass
from bmdrc.LPRClass import LPRClass

# Import specific support_functions for this main pipeline function
from support_functions import (
    combine_datasets,
    preprocess_morpho,
    run_filters,
    write_outputs,
)

# Example commands

## morphology only: python3 main.py --morpho test_files/test_morphology.csv
## lpr only: python3 main.py --lpr test_files/test_behavioral.csv
## morphology & lpr: python3 main.py --morpho test_files/test_morphology.csv --lpr test_files/test_behavioral.csv --both True

## python3 main.py --morpho files/Tanguay_Phase_4_zf_104alkyl_PAH_morphology_data_PNNL_2023OCT05.csv files/Zfish_Morphology_Legacy_2011-2018.csv --lpr files/Tanguay_Phase_4_zf_104alkyl_PAH_LPR_data_PNNL_2023OCT05.csv

## python3 main.py --lpr files/Tanguay_Phase_4_zf_104alkyl_PAH_LPR_data_PNNL_2023OCT05.csv

###########################
## COLLECT CLI ARGUMENTS ##
###########################
parser = argparse.ArgumentParser(
    "Run the QC and BMD analysis for the SRP analytics compendium"
)

parser.add_argument(
    "--morpho",
    dest="morpho",
    nargs="+",
    help=(
        "Pathway to the morphological file to be processed. "
        "Assumed format is long and required column names are: "
        "chemical.id, conc, plate.id, well, variable, value."
    ),
    default=None,
)
parser.add_argument(
    "--lpr",
    dest="lpr",
    nargs="+",
    help=(
        "Pathway to the light photometer response (LPR) file to be "
        "processed. Assumed format is long. Required columns are: "
        "chemical.id, conc, plate.id, well, variable, value."
    ),
    default=None,
)
parser.add_argument(
    "--output_dir",
    dest="output_dir",
    help="The output folder for files. Default is current directory.",
    default=".",
)
parser.add_argument(
    "--prefix",
    "--file_prefix",
    "--output_prefix",
    dest="file_prefix",
    help="Filename prefix for saved output files.",
    default="zebrafish",
)
parser.add_argument(
    "--report",
    dest="rep",
    help="Generate a report in the output folder. Default is True",
    action="store_true",
    default=True,
)


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

    ### 0. Pre-launch combination-------------------------------------------------------------
    if args.morpho is not None:
        print("...Concatenating morpho datasets")
        morpho_data = combine_datasets(morpho_paths)

    if args.lpr is not None:
        print("...Concatenating LPR datasets")
        lpr_data = combine_datasets(lpr_paths)

    ### 1. Input Data Modules--------------------------------------------------------------------
    if args.morpho is not None:
        print("...Formatting morphology data")
        BC = BinaryClass(
            df=morpho_data,
            chemical="chemical.id",
            concentration="conc",
            plate="plate.id",
            well="well",
            endpoint="endpoint",
            value="value",
            format="long",
        )

    if args.lpr is not None:
        print("...Formatting LPR data")
        LPR = LPRClass(
            df=lpr_data,
            chemical="chemical.id",
            concentration="conc",
            plate="plate.id",
            well="well",
            time="variable",
            value="value",
            cycle_length=20.0,
            cycle_cooldown=10.0,
            starting_cycle="light",
        )

    ### 2. Pre-Processing modules-----------------------------------------------------------------
    if args.morpho is not None:
        print("...Pre-Processing morphology data")
        preprocess_morpho(BC)

    # LPR data has MORT and MO24 fish set to NA

    ### 3. Filtering Modules----------------------------------------------------------------------
    if args.morpho is not None:
        print("...Filtering morphology data")
        run_filters(BC)

    if args.lpr is not None:
        print("...Filtering LPR data")
        run_filters(LPR)

    ### 4. Model Fitting Modules------------------------------------------------------------------
    if args.morpho is not None:
        print("...Fitting models to morphology data")
        BC.fit_models(diagnostic_mode=True)

    if args.lpr is not None:
        print("...Fitting models to LPR data")

        # Subset down to AUC2 and MOV2 and rename them "AUC" and "MOV"
        LPR.plate_groups = LPR.plate_groups[
            LPR.plate_groups[LPR.endpoint].isin(["AUC2", "MOV2"])
        ]
        LPR.plate_groups[LPR.endpoint] = LPR.plate_groups[LPR.endpoint].str.replace(
            "\d+", "", regex=True
        )
        LPR.plate_groups["bmdrc.Endpoint.ID"] = (
            LPR.plate_groups[LPR.chemical].astype(str)
            + " "
            + LPR.plate_groups[LPR.endpoint].astype(str)
        )
        LPR.fit_models(diagnostic_mode=True)

    ### 5. Format and export outputs------------------------------------------------------------
    print("...Exporting Results")

    # Export chemical data
    if args.morpho is not None:
        write_outputs(BC, "BC", args.output_dir, f"{args.file_prefix}_chem")
    if args.lpr is not None:
        write_outputs(LPR, "LPR", args.output_dir, f"{args.file_prefix}_chem")

    # Export sample data
    print("...Checking output")
    for fclass in ["BMDs", "Doses", "Fits"]:
        output_files = glob(f"{args.output_dir}/{args.file_prefix}_{fclass}_*.csv")
        for of in output_files:
            cmd = f"linkml-validate --schema srpAnalytics.yaml {of} --target-class zf{fclass}"
            print(cmd)
            os.system(cmd)


if __name__ == "__main__":
    main()
