#############
## IMPORTS ##
#############
import os
from typing import Optional, Union

import numpy as np
import pandas as pd
from bmdrc.BinaryClass import BinaryClass
from bmdrc.LPRClass import LPRClass

# Define type for data classes
BmdrcDataClass = Union[BinaryClass, LPRClass]

#########################
## COMBINE DATA.FRAMES ##
#########################

## TO FIX: log probit erroring out on some curves (how to fix infinity estimate)


# A support function to concatenate datasets together
def combine_datasets(thePaths):

    # Create a vector to hold data
    theData = []

    # Read all files
    for thePath in thePaths:
        theData.append(pd.read_csv(thePath))

    # Combine dataframes
    theData = pd.concat(theData)

    # Replace spaces
    theData["chemical.id"] = theData["chemical.id"].astype(str)
    theData["chemical.id"] = theData["chemical.id"].str.replace(" ", "_")

    return theData


######################################
## PRE-PROCESSING SUPPORT FUNCTIONS ##
######################################
# All pre-processing required for morphology data
def preprocess_morpho(BC):

    ## Endpoint Check--------------------------------------------------------------------------------------

    # Extract out the endpoints
    theEndpoints = BC.df[BC.endpoint].unique().tolist()

    # List the relevant endpoints, which is different for BRAIN samples
    if "BRAI" in theEndpoints:
        relevant_endpoints = [
            "AXIS",
            "BRAI",
            "BRN_",
            "CRAN",
            "DNC_",
            "CFIN",
            "CIRC",
            "DNC_",
            "DP24",
            "EYE_",
            "EDEM",
            "JAW_",
            "LTRK",
            "MO24",
            "MORT",
            "MUSC",
            "NC24",
            "NC__",
            "OTIC",
            "PE__",
            "PFIN",
            "PIG_",
            "SKIN",
            "SM24",
            "SNOU",
            "SOMI",
            "SWIM",
            "TCHR",
            "TRUN",
            "TR__",
            "YSE_",
        ]
    else:
        relevant_endpoints = [
            "AXIS",
            "BRN_",
            "CRAN",
            "DNC_",
            "DP24",
            "EDEM",
            "LTRK",
            "MO24",
            "MORT",
            "MUSC",
            "NC__",
            "SKIN",
            "SM24",
            "TCHR",
        ]

    # Get endpoints that are not expected
    unexpected = [end for end in theEndpoints if end not in relevant_endpoints]

    # Print a message for missing endpoints
    if len(unexpected) != 0:
        print(
            "......Pre-Processing: The following endpoints were discovered and are unexpected:",
            unexpected,
        )

    # Subset down to the relevant endpoints
    BC.df = BC.df[BC.df[BC.endpoint].isin(relevant_endpoints)]

    ## Convert do not counts to NA-------------------------------------------------------------------------

    # If there's a do not count category, remove data
    if "DNC_" in theEndpoints:

        # Let the user know that data have been removed
        print(
            "......Pre-Processing: A do not count category was detected and those wells were changed to NA."
        )

        # Set wells with a "do not count" to NA
        BC.set_well_to_na(endpoint_name="DNC_", endpoint_value=1)

        # Remove the DNC category
        BC.remove_endpoints("DNC_")

    ## Convert Mortality to NA------------------------------------------------------------------------------

    # Convert wells affected by mortality at 5 days to NA
    if "MORT" in theEndpoints:

        # Let users know that the affected organisms were changed to NA
        print(
            "......Pre-Processing: Mortality data at 5 days was detected, and affected endpoints in wells were changed to NA."
        )

        # Set wells with mortality at 5 days to NA, with the exception of the 24 hour timepoints
        BC.set_well_to_na(
            endpoint_name="MORT",
            endpoint_value=1,
            except_endpoint=["DP24", "MO24", "SM24", "MORT"],
        )

    # Convert wells affected by mortality at 24 hours to NA
    if "MO24" in theEndpoints:

        # Let users know that the affected organisms were changed to NA
        print(
            "......Pre-Processing: Mortality data at 24 hours was detected and affected endpoints in wells were changed to NA."
        )

        # Set wells at mortality at 24 hours to NA
        BC.set_well_to_na(
            endpoint_name="MO24", endpoint_value=1, except_endpoint="MO24"
        )

    ## Make new endpoints------------------------------------------------------------------------------------

    # Let users know that endpoints are being added
    print("......Pre-Processing: Adding additional endpoints")

    if "BRAI" in theEndpoints:

        # Remove any endpoints taht are calculated by us
        BC.df = BC.df[
            np.logical_not(
                BC.df[BC.endpoint].isin(
                    [
                        "ANY24",
                        "ANY120",
                        "TOT_MORT",
                        "ALL_BUT_MORT",
                        "BRN_",
                        "CRAN",
                        "EDEM",
                        "LTRK",
                        "MUSC",
                        "SKIN",
                        "TCHR",
                    ]
                )
            )
        ]

        # Define new endpoints in a dictionary
        EndpointDictionary = {
            "ANY24": ["MO24", "DP24", "SM24", "NC24"],
            "ANY120": [
                "MORT",
                "YSE_",
                "AXIS",
                "EYE_",
                "SNOU",
                "JAW_",
                "OTIC",
                "PE__",
                "BRAI",
                "SOMI",
                "PFIN",
                "CFIN",
                "PIG_",
                "CIRC",
                "TRUN",
                "SWIM",
                "NC__",
                "TR__",
                "ANY24",
            ],
            "TOT_MORT": ["MO24", "MORT"],
            "ALL_BUT_MORT": [
                "DP24",
                "SM24",
                "NC24",
                "YSE_",
                "AXIS",
                "EYE_",
                "SNOU",
                "JAW_",
                "OTIC",
                "PE__",
                "BRAI",
                "SOMI",
                "PFIN",
                "CFIN",
                "PIG_",
                "CIRC",
                "TRUN",
                "SWIM",
                "NC__",
                "TR__",
            ],
            "BRN_": ["BRAI", "OTIC", "PFIN"],
            "CRAN": ["EYE_", "SNOU", "JAW_"],
            "EDEM": ["YSE_", "PE__"],
            "LTRK": ["TRUN", "CFIN"],
            "MUSC": ["CIRC", "SWIM", "SOMI"],
            "SKIN": ["PIG_"],
            "TCHR": ["TR__"],
        }

        # Add new endpoints
        BC.combine_and_create_new_endpoints(EndpointDictionary)

        # Remove renamed endpoints
        # BC.remove_endpoints(["PIG_", "TR__"])

    else:

        # Define new endpoints in a dictionary
        EndpointDictionary = {
            "ANY24": ["MO24", "DP24", "SM24"],
            "ANY120": [
                "AXIS",
                "BRN_",
                "CRAN",
                "EDEM",
                "LTRK",
                "MORT",
                "MUSC",
                "NC__",
                "SKIN",
                "TCHR",
                "ANY24",
            ],
            "TOT_MORT": ["MO24", "MORT"],
            "ALL_BUT_MORT": [
                "AXIS",
                "BRN_",
                "CRAN",
                "DP24",
                "EDEM",
                "LTRK",
                "MUSC",
                "NC__",
                "SKIN",
                "SM24",
                "TCHR",
            ],
        }

        # Add new endpoints
        BC.combine_and_create_new_endpoints(EndpointDictionary)


#################################
## FILTERING SUPPORT FUNCTIONS ##
#################################


# Run filters on either the BC or LPR object
def run_filters(obj):

    # Negative control filter, default 50%
    obj.filter_negative_control(apply=True, diagnostic_plot=False)

    # Minimum concentration filter, default is 3
    obj.filter_min_concentration(apply=True, diagnostic_plot=False)

    # Correlation score filter
    obj.filter_correlation_score(apply=True, diagnostic_plot=False)


###################
## WRITE OUTPUTS ##
###################
def write_outputs(obj: BmdrcDataClass, tag: str, outname: Optional[str]):
    """Write all output files: BMDs, Dose, Fits, and Report markdown.

    Parameters
    ----------
    obj : BmdrcDataClass
        `bmdrc.BinaryClass.BinaryClass` or `bmdrc.LPRClass.LPRClass`
    tag : str
        Filename tag
    outname : Optional[str]
        Name of output directory
    """
    if outname is None:
        out = "."
    else:
        out = str(outname)

    # Output bmds: Chemical_ID, End_Point, Model, BMD10, BMDL, BMD50, AUC, Min_Dose, Max_Dose, AUC_Norm,
    # DataQC_Flag, BMD_Analysis_Flag, BMD10_Flag, BMD50_Flag, ids
    obj.output_benchmark_dose(f"{out}/new_BMDS_{tag}.csv")

    # Output dose: Chemical_ID, End_Point, Dose, num.affected, num.nonna, ids, CI_Lo, CI_Hi
    obj.output_dose_table(f"{out}/new_Dose_{tag}.csv")

    # Output fits: Chemical_ID, End_Point, X_vals, Y_vals
    obj.output_fits_table(f"{out}/new_Fits_{tag}.csv")

    # Output reports
    obj.report(f"{out}/new_report_{tag}/")
