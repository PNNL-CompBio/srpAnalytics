"""map_samples_to_chemicals: Python version of mapSamplesToChems.R

Original rewrite using generative AI; modified by @christinehc
"""

# =========================================================
# Imports
# =========================================================
import os
import sys
from argparse import ArgumentParser
from typing import Optional

import numpy as np
import pandas as pd

from numpy.typing import ArrayLike

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.format import rename_duplicates, snakeify, snakeify_all_columns
from src.mapping import rename_chemical_class
from src.metadata import build_chem_metadata, get_endpoint_metadata
from src.params import (
    MASV_CC,
    MASV_SOURCE,
    QC_FLAGS,
    REQUIRED_BMD_COLUMNS,
    REQUIRED_SAMPLE_COLUMNS,
    SAMPLE_CHEM_COLUMNS,
    SAMPLE_COLUMNS,
)
from src.tables import sample_id_master_table

# These pathways refer to absolute pathways in the docker image
# setting these three parameters, can be appended
# data_dir = 'https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data'
# OUT_DIR = "/tmp"
OUT_DIR = "."

# Set CompTox API key
CTX_API_KEY = "5aded20c-9485-11ef-87c3-325096b39f47"

# =========================================================
# Functions
# =========================================================


###################################
# Metadata Collection
###################################


# FIXME: Seems unused?
#'getNewChemicalClass
#'This file reads in the files and processes the chemical class names
#'to be friendly for the website
#' @depracated As we move to the new MASV classes
#'@param data.dir
#'@return data.frame
def get_new_chemical_class(data_dir: str) -> pd.DataFrame:
    # Read in the Excel sheet for PAHs
    pahs = (
        pd.read_excel(
            os.path.join(data_dir, "PAH_and_1530_SRP_Summary.xlsx"), sheet_name=3
        )
        .loc[:, ["casrn"]]
        .rename(columns={"casrn": "cas_number"})
        .assign(classification="PAH")
    )

    # Add extra data for particular CAS numbers
    extras = pd.DataFrame(
        {
            "cas_number": ["3074-03-01", "7496-02-08", "6373-11-01"],
            "classification": "PAH",
        }
    )

    # Read in the Excel sheet for non-PAHs
    non_pahs = (
        pd.read_excel(
            os.path.join(data_dir, "PAH_and_1530_SRP_Summary.xlsx"), sheet_name=4
        )
        .loc[:, ["casrn", "classification"]]
        .rename(columns={"casrn": "cas_number"})
    )

    # Combine all the dataframes
    full_class = pd.concat([pahs, extras, non_pahs], ignore_index=True)
    full_class["newClass"] = full_class["classification"].apply(rename_chemical_class)

    return full_class


def process_fses(filename: str, snake_case: bool = True) -> pd.DataFrame:
    # Replace invalid values with nulls for filtering
    df = pd.read_csv(filename)[REQUIRED_SAMPLE_COLUMNS].replace(
        {"BLOD": "0", "NULL": "0", "nc:BDL": "0"}
    )
    if snake_case:
        df = snakeify_all_columns(df)

    # Remove null and invalid entries
    df = df[
        (df["sample_number"].notna())
        & (df["cas_number"].notna())
        & (~df["measurement_value"].isin(["0", np.nan]))
        & (~df["measurement_value_molar"].isin(["0", np.nan]))
    ]

    # Format FSES location data
    df["location_lon"] = pd.to_numeric(df["location_lon"], errors="coerce")

    # Only allow negative longitudes
    # Note: Our data is already all negative, so unnecessary?
    # df["location_lon"] = np.where(
    #     df["location_lon"].gt(0), -df["location_lon"], df["location_lon"]
    # )
    return df


def build_sample_data(
    fses_files: list[str],
    chem_metadata: pd.DataFrame,
    sample_id_file: str,
    sample_mapping: Optional[str] = None,
) -> pd.DataFrame:
    """Select relevant data from curated tables.

    Parameters
    ----------
    fses_files : list[str]
        List of FSES files from the Barton lab with sample info
    chem_metadata : pd.DataFrame
        Chemical metadata table containing identifier mapping
    sample_id_file : str
        /path/to/sample_id_file (CSV file)
    sample_mapping : Optional[str], optional
        Mapping file containing new data, by default None

    Returns
    -------
    pd.DataFrame
        _description_
    """
    data = []
    # Read and preprocess each file in fses_files
    for f in fses_files:
        # Replace invalid values with nulls for filtering
        sample = process_fses(f)
        data.append(sample)

    # Concatenate samples and merge with chemical metadata
    data = pd.concat(data).merge(
        chem_metadata[["chemical_id", "cas_number", "average_mass"]],
        on="cas_number",
        how="left",
    )

    # Merge with sample IDs
    ids = sample_id_master_table(data["sample_number"], sample_id_file)
    data = data.merge(ids, on="sample_number", how="left").drop_duplicates()

    # Rename duplicate sample names as sample:01, :02, etc.
    data["sample_name"] = rename_duplicates(data)

    # Merge with sample name remappings if provided
    if sample_mapping is not None:
        sample_remap_cols = [
            "sample_id",
            "new_project_name",
            "new_sample_name",
            "new_location_name",
        ]
        sample_name_remap = (
            snakeify_all_columns(pd.read_excel(sample_mapping, sheet_name=0))
            .rename(columns={"project_name": "new_project_name"})
            .loc[:, sample_remap_cols]
            .drop_duplicates()
        )

        data = data.merge(
            sample_name_remap, on="sample_id", how="left"
        ).drop_duplicates()

        # Fill in NAs with new values from remapping table
        nas = data["project_name"].isna()
        for col in ["project_name", "location_name", "sample_name"]:
            data.loc[nas, col] = data.loc[nas, f"new_{col}"]

    # Drop unnecessary columns and rows with missing cas_number
    data = data.drop(
        columns=[
            "new_project_name",
            "new_sample_name",
            "new_location_name",
            "average_mass",
        ]
    ).dropna(subset=["cas_number"])

    return data


def combine_v2_chemical_endpoint_data(
    bmd_files: list[str],
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_details: Optional[pd.DataFrame] = None,
):
    """Combine chemical endpoint data

    Parameters
    ----------
    bmd_files : list[str]
        List of BMD files
    is_extract : bool, optional
        True if data is for extracts, by default False
    chem_data : Optional[pd.DataFrame], optional
        Tabulated chemical data, by default None
    endpoint_details : Optional[pd.DataFrame], optional
        Tabulated endpoint data, by default None

    Returns
    -------
    _type_
        _description_
    """
    print(f"Combining bmd files: {', '.join(bmd_files)}")

    # Read and concatenate the specified columns from all provided BMD files into one df
    cols = REQUIRED_BMD_COLUMNS["bmd"]
    df = pd.concat([pd.read_csv(bmd_file).loc[:, cols] for bmd_file in bmd_files])
    df = snakeify_all_columns(df)  # .rename(columns={"chemical_id": "tmp_id"})

    # Remove duplicate entries (based on 'chemical_id' and 'end_point')
    print(df.columns)
    df = df.drop_duplicates(subset=["chemical_id", "end_point"])

    # Process extracts
    if is_extract:
        chem_data = snakeify_all_columns(chem_data)
        print(chem_data.columns)

        # Split `sample_id` column on '-', take first two parts, and rename split cols
        chem_data[["tmp_id", "sub"]] = chem_data["sample_id"].str.split(
            "-", expand=True, n=1
        )
        print(chem_data.columns)
        # df["chemical_id"] = df["chemical_id"].astype(str)
        df = df.merge(
            chem_data[["chemical_id", "sample_id"]], on="chemical_id", how="outer"
        )

        # Fill NAs and format sample names
        df["sample_id"] = df["sample_id"].fillna(df["chemical_id"])
        df["sample_name"] = df["sample_id"].apply(
            lambda x: f"Sample {x}"
            if pd.isnull(df.get("sample_name"))
            else df["sample_name"]
        )

        # Fill missing 'end_point' values with "NoData", merge with `endpoint_details` on 'end_point',
        # drop unused columns, remove rows with missing 'sample_id', drop duplicates,
        # and fill remaining missing 'location_name' values with "None"
        print(endpoint_details.columns)
        # print(endpoint_details.head(2))
        df = (
            df.fillna({"end_point": "NoData"})
            .merge(endpoint_details, on="end_point", how="right")
            .drop(columns=["end_point", "chemical_id"])
            .dropna(subset=["sample_id"])
            .drop_duplicates()
            .fillna({"location_name": "None"})
        )

    # Similar operations for non-extract case:
    # Fill missing 'end_point' values with "NoData", merge with `endpoint_details` on 'end_point'
    # drop unused column, remove rows with missing 'cas_number', drop duplicates,
    # and fill remaining missing 'chemical_class' values with "Unclassified"
    else:
        df = (
            df.fillna({"end_point": "NoData"})
            .merge(endpoint_details, on="end_point", how="right")
            .drop(columns="end_point")
            .dropna(subset=["cas_number"])
            .drop_duplicates()
            .fillna({"chemical_class": "Unclassified"})
        )

    # Grade data QC on scale and replace NULL values in model col
    df["data_qc_flag"] = df["data_qc_flag"].map(QC_FLAGS).fillna("Good")
    df["model"] = df["model"].str.replace("NULL", "None")

    return df


def combine_chemical_data(
    bmd_files: list[str],
    data_type: str = "fit",
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_details: Optional[pd.DataFrame] = None,
):
    print(f'Combining {data_type} files: {", ".join(bmd_files)}')

    col_type = (
        "fitVals"
        if data_type == "fit"
        else "doseRep"
        if data_type == "dose"
        else ValueError("Invalid data_type. Must be 'fit' or 'dose'.")
    )
    cols = REQUIRED_BMD_COLUMNS[col_type]

    df = pd.concat([pd.read_csv(file)[cols] for file in bmd_files])
    df = snakeify_all_columns(df)

    df["combined"] = df[["chemical_id", "end_point"]].astype(str).agg(" ".join, axis=1)
    df = df.drop_duplicates(subset=["combined"], keep="last")

    if data_type == "fit":
        df = df[df["x_vals"] != "NULL"]
        df["x_vals"] = pd.to_numeric(df["x_vals"])
        df["y_vals"] = pd.to_numeric(df["y_vals"])

    df = df.merge(
        snakeify_all_columns(endpoint_details),
        on="end_point",
        how="right",
    )
    df = df.drop(columns=["end_point", "description"]).drop_duplicates()

    if is_extract:
        chem_data[["tmp_id", "sub"]] = chem_data["sample_id"].str.split(
            "-", expand=True
        )
        chem_data = chem_data[["sample_id", "tmp_id"]].drop_duplicates()

        df["tmp_id"] = df["chemical_id"].astype(str)
        df = df.drop(columns="chemical_id").merge(chem_data, on="tmp_id", how="left")

        df["sample_id"] = df["sample_id"].fillna(df["tmp_id"])
        df = df.drop(columns="tmp_id")

    return df.drop_duplicates()


def _flatten_class_df(
    df: pd.DataFrame,
    keep_cols: ArrayLike,
    drop_cols: ArrayLike,
    var_name: str,
    id_cols: ArrayLike = ["CASNumber", "ParameterName"],
) -> pd.DataFrame:
    """Flatten and format long-form presence/absence class data.

    Parameters
    ----------
    df : pd.DataFrame
        Long-form class dataframe
    keep_cols : ArrayLike
        List of columns to keep (e.g. MASV_CC or MASV_SOURCE)
    drop_cols : ArrayLike
        List of columns to drop (e.g. MASV_CC or MASV_SOURCE)
    var_name : str
        Variable name (e.g. "chem_source", "chem_class", ...)
    id_cols : ArrayLike, optional
        List of chemical ID cols, by default ["CASNumber", "ParameterName"]

    Returns
    -------
    pd.DataFrame
        Formatted dataframe
    """
    df = df.drop(columns=drop_cols)  # drop non-source cols
    df = df.melt(  # reorg by cas, param name, source/class, and pos/neg
        id_vars=id_cols,
        value_vars=keep_cols,
        var_name=var_name,
        value_name="posNeg",
    )

    # Remove null pos/negs, then remove pos/neg col
    df = df[df["posNeg"] != "NULL"]  # remove nulls
    df = df.drop(columns="posNeg")

    # For sources, combine all pesticides into 1 category
    if "source" in var_name:
        df[var_name] = df[var_name].str.replace(r"^pest.*", "pesticide", regex=True)

    # Compile all values for each compound
    df = df.groupby(id_cols)[var_name].apply(lambda x: ";".join(x)).reset_index()
    return df


def masv_chem_class(
    class_file: str,
    save_to: str = "MASV_classAndSource.csv",
    id_cols: ArrayLike = ["CASNumber", "ParameterName"],
) -> pd.DataFrame:
    """Reads full MASV class annotations and assigns values to chemicals.

    Parameters
    ----------
    class_file : str
        Filename of class annotation data.
    save_to : str, optional
        Desired output filename, by default "MASV_classAndSource.csv"
    id_cols : ArrayLike, optional
        List of chemical ID cols, by default ["CASNumber", "ParameterName"]

    Returns
    -------
    pd.DataFrame
        Combined chemical class information
    """
    data = pd.read_excel(class_file, sheet_name=0)

    # Reorganize data by source
    sources = data.copy()
    sources = _flatten_class_df(
        sources,
        keep_cols=MASV_SOURCE,
        drop_cols=MASV_CC,
        var_name="chem_source",
        id_cols=id_cols,
    )

    # Reorganize data by class
    classes = data.copy()
    classes = _flatten_class_df(
        classes,
        keep_cols=MASV_CC,
        drop_cols=MASV_SOURCE,
        var_name="chemical_class",
        id_cols=id_cols,
    )

    # Combine source and class information
    combined = sources.merge(classes, on=id_cols, how="outer")
    combined = combined.fillna({"chemical_class": "Unclassified"})
    combined = combined.rename(columns={c: snakeify(c) for c in id_cols})

    combined.to_csv(save_to, index=False)
    return combined


# =========================================================
# Command Line Interface (CLI)
# =========================================================


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-s",
        "--sample",
        dest="is_sample",
        action="store_true",
        default=False,
        help="Flag to indicate file is for samples",
    )
    parser.add_argument(
        "-c",
        "--chemical",
        dest="is_chem",
        default=False,
        action="store_true",
        help="Flag to indicate file is for chemicals",
    )
    parser.add_argument(
        "-d",
        "--drc_files",
        dest="dose_response",
        default="",
        help="Dose response curve file",
    )
    parser.add_argument(
        "-p",
        "--sample_id",
        dest="sample_id_file",
        default="",
        help="Sample mapping file location",
    )
    parser.add_argument(
        "-i",
        "--chem_id",
        dest="chem_id_file",
        default="",
        help="Chemical ID file location",
    )
    parser.add_argument(
        "-e",
        "--ep_map",
        dest="endpoint_mapping_file",
        default="",
        help="Endpoint naming file location",
    )
    parser.add_argument(
        "-l",
        "--chem_class",
        dest="chem_class_file",
        default="",
        help="Chemical class file location",
    )
    parser.add_argument(
        "-f",
        "--sample_files",
        dest="sample_files",
        default="",
        help="Comma delimited list of FSES files to merge",
    )
    parser.add_argument(
        "-y",
        "--chem_desc",
        dest="chem_desc",
        default="",
        help="Descriptions of chemicals",
    )
    parser.add_argument(
        "-m",
        "--sample_map",
        dest="sample_map",
        default="",
        help="File that maps sample locations",
    )
    parser.add_argument(
        "--metadata",
        dest="metadata",
        default="https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data/srp_build_files.csv",
        help="Metadata file location (i.e. srp_build_files.csv)",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        default="/tmp/",
        help="File that maps sample locations",
    )

    args = parser.parse_args()
    print(args)

    chem_class = masv_chem_class(args.chem_class_file)
    chem_metadata = build_chem_metadata(args.metadata)

    sample_files_list = args.sample_files.split(",")
    chem_sample = build_sample_data(
        sample_files_list, chem_metadata, args.sample_id_file, args.sample_map
    )

    endpoint_details = get_endpoint_metadata(
        args.endpoint_mapping_file
    ).drop_duplicates()

    if args.is_sample or args.is_chem:
        all_files = args.dose_response.split(",")
        bmd_files = [file for file in all_files if "bmd" in file]
        dose_files = [file for file in all_files if "dose" in file]
        fit_files = [file for file in all_files if "fit" in file]

        chem_data = chem_sample if args.is_sample else chem_metadata

        bmds = (
            combine_v2_chemical_endpoint_data(
                bmd_files,
                is_extract=args.is_sample,
                chem_data=chem_data,
                endpoint_details=endpoint_details,
            )
            .dropna(subset=["bmd_analysis_flag"])
            .query("bmd_analysis_flag != 'NA'")
        )
        curves = combine_chemical_data(
            fit_files,
            data_type="fit",
            is_extract=args.is_sample,
            chem_data=chem_data,
            endpoint_details=endpoint_details,
        )
        dose_reps = combine_chemical_data(
            dose_files,
            data_type="dose",
            is_extract=args.is_sample,
            chem_data=chem_data,
            endpoint_details=endpoint_details,
        ).dropna(subset=["dose"])

        if args.is_sample:
            bmds.to_csv(
                os.path.join(args.output_dir, "zebrafishSampBMDs.csv"),
                index=False,
                quotechar='"',
            )
            curves.to_csv(
                os.path.join(args.output_dir, "zebrafishSampXYCoords.csv"),
                index=False,
                quotechar='"',
            )
            dose_reps.to_csv(
                os.path.join(args.output_dir, "zebrafishSampDoseResponse.csv"),
                index=False,
                quotechar='"',
            )

        if args.is_chem:
            nas = bmds[bmds["Chemical_ID"].isna()]["Chemical_ID"]
            to_remove = set(nas) - set(chem_sample["chemical_id"])
            bmds = bmds[~bmds["Chemical_ID"].isin(to_remove)]
            curves = curves[~curves["Chemical_ID"].isin(to_remove)]
            dose_reps = dose_reps[~dose_reps["Chemical_ID"].isin(to_remove)]

            bmds.to_csv(
                os.path.join(args.output_dir, "zebrafishChemBMDs.csv"),
                index=False,
                quotechar='"',
            )
            curves.to_csv(
                os.path.join(args.output_dir, "zebrafishChemXYCoords.csv"),
                index=False,
                quotechar='"',
            )
            dose_reps.to_csv(
                os.path.join(args.output_dir, "zebrafishChemDoseResponse.csv"),
                index=False,
                quotechar='"',
            )

    else:
        chem_metadata.to_csv(
            os.path.join(args.output_dir, "chemicals.csv"), index=False, quotechar='"'
        )
        chem_sample[SAMPLE_COLUMNS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "samples.csv"), index=False, quotechar='"'
        )
        chem_sample[SAMPLE_CHEM_COLUMNS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "sampleToChemicals.csv"),
            index=False,
            quotechar='"',
        )


main()
