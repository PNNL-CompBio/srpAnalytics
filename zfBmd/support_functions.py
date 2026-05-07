#############
## IMPORTS ##
#############
import os
from typing import Optional, Union

import numpy as np
import pandas as pd
from bmdrc.BinaryClass import BinaryClass
from bmdrc.LPRClass import LPRClass
from params import (
    BRAIN_ENDPOINT_DICT,
    BRAIN_ENDPOINTS,
    CALC_ENDPOINTS,
    ENDPOINT_DICT,
    ENDPOINTS,
)

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
def preprocess_morpho(BC: BinaryClass):
    """Run all pre-processing required for morphology data.

    Parameters
    ----------
    BC : BinaryClass
        Loaded morphology data object.
    """
    # Get relevant endpoints (different for brain)
    endpoints = BC.df[BC.endpoint].unique().tolist()
    relevant = BRAIN_ENDPOINTS if "BRAI" in endpoints else ENDPOINTS
    BC.df = BC.df[BC.df[BC.endpoint].isin(relevant)]

    # Find unexpected endpoints
    unexpected = [e for e in endpoints if e not in relevant]
    if len(unexpected) != 0:
        print(
            "......Pre-Processing: The following endpoints "
            f"were discovered and are unexpected: {unexpected}",
        )

    # Convert "do not count"s to N/A and remove
    if "DNC_" in endpoints:
        print(
            "......Pre-Processing: A do not count category was "
            "detected and those wells were changed to NA."
        )
        BC.set_well_to_na(endpoint_name="DNC_", endpoint_value=1)
        BC.remove_endpoints("DNC_")

    # Convert wells affected by mortality at 5 days to N/A
    if "MORT" in endpoints:
        print(
            "......Pre-Processing: Mortality data at 5 days was "
            "detected, and affected endpoints in wells were changed to NA."
        )

        # Set wells with mortality at 5 days to N/A, minus the 24 hour timepoints
        BC.set_well_to_na(
            endpoint_name="MORT",
            endpoint_value=1,
            except_endpoint=["DP24", "MO24", "SM24", "MORT"],
        )

    # Convert wells affected by mortality at 24 hours to N/A
    if "MO24" in endpoints:
        print(
            "......Pre-Processing: Mortality data at 24 hours was "
            "detected and affected endpoints in wells were changed to NA."
        )
        BC.set_well_to_na(
            endpoint_name="MO24", endpoint_value=1, except_endpoint="MO24"
        )

    # Create additional endpoints
    print("......Pre-Processing: Adding additional endpoints")
    if "BRAI" in endpoints:
        # Remove any endpoints that are calculated by us
        BC.df = BC.df[np.logical_not(BC.df[BC.endpoint].isin(CALC_ENDPOINTS))]

        # Define new endpoints to be added
        endpoint_dict = BRAIN_ENDPOINT_DICT
        BC.combine_and_create_new_endpoints(endpoint_dict)

        # Remove renamed endpoints
        # BC.remove_endpoints(["PIG_", "TR__"])

    else:
        # Define new endpoints to be added
        endpoint_dict = ENDPOINT_DICT
        BC.combine_and_create_new_endpoints(endpoint_dict)


#################################
## FILTERING SUPPORT FUNCTIONS ##
#################################
# Run filters on either the BC or LPR object
def run_filters(dataset: BmdrcDataClass):

    # Negative control filter, default 50%
    dataset.filter_negative_control(apply=True, diagnostic_plot=False)

    # Minimum concentration filter, default is 3
    dataset.filter_min_concentration(apply=True, diagnostic_plot=False)

    # Correlation score filter
    dataset.filter_correlation_score(apply=True, diagnostic_plot=False)


###################
## WRITE OUTPUTS ##
###################
def write_outputs(
    dataset: BmdrcDataClass,
    tag: str,
    save_to: Optional[str],
    data_type: str = "zebrafish",
):
    """Write all output files: BMDs, Dose, Fits, and Report markdown.

    Parameters
    ----------
    obj : BmdrcDataClass
        `bmdrc.BinaryClass.BinaryClass` or `bmdrc.LPRClass.LPRClass`
    tag : str
        Filename tag
    save_to : Optional[str]
        Name of output directory
    data_type : str
        Data type, e.g. zebrafish or human, by default "zebrafish"
    """
    save_to = "." if save_to is None else str(save_to)

    # Output bmds: Chemical_ID, End_Point, Model, BMD10, BMDL, BMD50, AUC, Min_Dose, Max_Dose, AUC_Norm,
    # DataQC_Flag, BMD_Analysis_Flag, BMD10_Flag, BMD50_Flag, ids
    dataset.output_benchmark_dose(f"{save_to}/{data_type}_BMDs_{tag}.csv")

    # Output dose: Chemical_ID, End_Point, Dose, num.affected, num.nonna, ids, CI_Lo, CI_Hi
    dataset.output_dose_table(f"{save_to}/{data_type}_Dose_{tag}.csv")

    # Output fits: Chemical_ID, End_Point, X_vals, Y_vals
    dataset.output_fits_table(f"{save_to}/{data_type}_Fits_{tag}.csv")

    # Output reports
    dataset.report(f"{save_to}/{data_type}_report_{tag}/")
