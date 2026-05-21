"""map_samples_to_chemicals: Python version of mapSamplesToChems.R

Original rewrite using generative AI; modified by @christinehc
"""

# =========================================================
# Imports
# =========================================================
import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from dotenv import load_dotenv
from numpy.typing import ArrayLike
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.data import FigshareDataLoader, load_figshare_url
from src.format import rename_duplicates
from src.manifest import DataManifest
from src.mapping import rename_chemical_class
from src.metadata import build_chem_metadata, format_endpoint_metadata
from src.params import MANIFEST_FILEPATH, MASV_CC, MASV_SOURCE
from src.schema import combine_schema_cols, get_cols_from_schema
from src.tables import sample_id_master_table

# =========================================================
# Setup and Parameters
# =========================================================
SAMPLE_COLS = get_cols_from_schema("samples")
SAMP2CHEM_COLS = get_cols_from_schema("samplesToChemicals")
FSES_COLS = combine_schema_cols("samples", "samplesToChemicals")

# These pathways refer to absolute pathways in the docker image
# setting these three parameters, can be appended
# data_dir = 'https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data'
OUTPUT_DIR = os.getenv("OUTPUT_DIR")
OUTPUT_DIR = OUTPUT_DIR if OUTPUT_DIR is not None else "tmp"
# OUTPUT_DIR = "."

# Set CompTox API key
load_dotenv()
CTX_API_KEY = os.getenv("CTX_API_KEY")

# Figshare file loader
loader = FigshareDataLoader(
    Path(OUTPUT_DIR) / ".figshare_cache",  # api_token=FIGSHARE_API_TOKEN
)
manifest = DataManifest(MANIFEST_FILEPATH)


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


def load_clean_fses(filename: str, loader: FigshareDataLoader) -> pd.DataFrame:
    """Clean FSES input files.

    Steps:
        1. Remove invalid and null entries
        2. Format location data and enforces negative longitudes
        3. Remove entries with no data

    Parameters
    ----------
    filename : str
        Path to FSES file
    loader : FigshareDataLoader
        Figshare data loader object

    Returns
    -------
    pd.DataFrame
        Cleaned FSES data
    """
    # Parse CSV files and HTTPS (figshare) separately
    if os.path.splitext(filename)[1] == ".csv":
        fses = pd.read_csv(filename, dtype={"Sample_ID": str, "Chemical_ID": str})
    elif os.path.splitext(filename)[1] == "":
        fses = load_figshare_url(loader, filename)

    # Clean entries
    fses = fses[FSES_COLS].replace({"BLOD": "0", "NULL": "0", "nc:BDL": "0"})
    fses = fses.drop(columns=["Chemical_ID"])  # Drop chem ID (for later merge)

    # Remove null and invalid entries
    fses = fses[
        (fses["SampleNumber"].notna())
        & (fses["SampleNumber"] != "None")
        & (fses["cas_number"] != "NULL")
        & (~fses["measurement_value_molar"].isin(["0"]))
        & (~fses["measurement_value"].isin(["0", "NULL", ""]))
    ]

    # Format FSES location data and require negative longitudes
    # NOTE: Our data is already all negative, so unnecessary?
    fses["LocationLon"] = pd.to_numeric(fses["LocationLon"], errors="coerce")
    fses["LocationLon"] = np.where(
        fses["LocationLon"].gt(0), -fses["LocationLon"], fses["LocationLon"]
    )

    # Drop entries with no data
    fses = fses.dropna(
        subset=[
            "SampleNumber",
            "SampleName",
            "date_sampled",
            "LocationLat",
            "LocationLon",
            "measurement_value",
            "environmental_concentration",
        ],
        how="all",
    )
    return fses


def build_sample_data(
    fses_files: list[str],
    loader: FigshareDataLoader,
    chem_metadata: pd.DataFrame,
    sample_id_file: str,
    sample_mapping: str = None,
):
    """Select relevant data from curated tables.

    Parameters
    ----------
    fses_files : list[str]
        List of FSES files from the Barton lab with sample info
    loader : FigshareDataLoader
        Figshare data loader object (for loading FSES files)
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
    # Read and process all FSES files
    data = pd.concat(
        [load_clean_fses(f, loader) for f in fses_files], ignore_index=True
    )

    # Add chemical metadata and sample IDs
    chem_metadata = chem_metadata[["Chemical_ID", "cas_number", "averageMass"]]
    chem_metadata = chem_metadata.drop_duplicates()
    data = data.merge(chem_metadata, on="cas_number", how="left")

    # Get sample IDs
    sample_ids = sample_id_master_table(data["SampleNumber"], sample_id_file)
    data = (
        data.merge(sample_ids, on="SampleNumber", how="left", suffixes=("_x", ""))
        .drop_duplicates()
        .drop(columns=["Sample_ID_x"])
    )

    # NOTE: Moved to end
    # Rename duplicate sample names as sample:01, :02, etc.
    # data["SampleName"] = rename_duplicates(data, col="SampleName")

    # Fill in missing concentration data
    blanks_mask = data["measurement_value_molar"] == ""
    if blanks_mask.any():
        data.loc[blanks_mask, "measurement_value_molar"] = (
            pd.to_numeric(data.loc[blanks_mask, "measurement_value"])
            * 1000
            / pd.to_numeric(data.loc[blanks_mask, "averageMass"])
        )

    # Remove chem data without CAS and clean up duplicates
    data = data[data["cas_number"].notna()].drop_duplicates()

    # Merge with sample name remappings if provided
    if sample_mapping is not None:
        sample_remap_cols = [
            "Sample_ID",
            "ProjectName",
            "NewSampleName",
            "NewLocationName",
        ]
        remap = pd.read_excel(sample_mapping, sheet_name=0)
        remap = remap[sample_remap_cols].drop_duplicates()
        data = data.merge(remap, on="Sample_ID", how="left").drop_duplicates()

        # Fill in NAs with new values from remapping table
        missing = data["projectName"].isna()
        if missing.any():
            # tqdm.write("Missing values found!")
            remap_col_names = {
                "projectName": "ProjectName",
                "LocationName": "NewLocationName",
                "SampleName": "NewSampleName",
            }
            for old_col, new_col in remap_col_names.items():
                data.loc[missing, old_col] = data.loc[missing, new_col]

            # tqdm.write"columns after missing filled in:", data.columns)

        # Drop unnecessary columns and rows with missing CAS
        data = data.drop(
            columns=[
                "ProjectName",
                "NewSampleName",
                "NewLocationName",
                "averageMass",
            ]
        )
        data = data.dropna(subset=["cas_number"])
        data = data[data["cas_number"] != "N/A"].drop_duplicates()
        data["cas_number"] = data["cas_number"].astype(str)  # force str

    # Rename duplicate sample names as sample:01, :02, etc.
    data["SampleName"] = rename_duplicates(data, col="SampleName")
    return data


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
    save_to: str = os.path.join(OUTPUT_DIR, "MASV_classAndSource.csv"),
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

    combined.to_csv(save_to, index=False)
    return combined


# =========================================================
# Command Line Interface (CLI)
# =========================================================
def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--sample_id_file",
        "--sample_id",
        "-p",
        dest="sample_id_file",
        default="",
        help="Sample mapping file location",
    )
    parser.add_argument(
        "--sample_map",
        "--sample_map_file",
        "-m",
        dest="sample_map",
        default="",
        help="File that maps sample locations",
    )
    parser.add_argument(
        "--chemical_id",
        "--chem_id",
        "-i",
        dest="chem_id_file",
        default="",
        help="Chemical ID file location",
    )
    parser.add_argument(
        "--endpoint_map",
        "--ep_map",
        "-e",
        dest="endpoint_mapping_file",
        default="",
        help="Endpoint naming file location",
    )
    parser.add_argument(
        "--chemical_class",
        "--chem_class",
        "--chemical_class_file",
        "--chem_class_file",
        "-l",
        dest="chem_class_file",
        default="",
        help="Chemical class file location",
    )
    parser.add_argument(
        "--sample_files",
        "-f",
        dest="sample_files",
        default="",
        help="Comma delimited list of FSES files to merge",
    )
    parser.add_argument(
        "--chemical_description",
        "--chemical_desc",
        "--chem_desc",
        "-y",
        dest="chem_desc",
        default="",
        help="Descriptions of chemicals",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        dest="output_dir",
        default=OUTPUT_DIR,
        help="File that maps sample locations",
    )
    args = parser.parse_args()

    # -----------------
    # Chem Metadata
    # -----------------
    tqdm.write("Getting chemical metadata...")
    chem_class = masv_chem_class(args.chem_class_file)
    chem_ids = load_figshare_url(loader, args.chem_id_file)

    if not os.path.exists(os.path.join(args.output_dir, "chem_metadata.tsv")):
        tqdm.write("No metadata found. Building metadata...")

        try:
            chem_metadata = build_chem_metadata(
                chem_ids,
                save_to=os.path.join(args.output_dir, "chem_metadata.tsv"),
            )
        except requests.exceptions.ConnectionError as e:
            # Get pre-built chemical metadata
            chem_metadata_file = manifest.get(name="metadata_2026-01")
            chem_metadata = load_figshare_url(loader, chem_metadata_file)

    else:
        tqdm.write("Metadata found! Reading from previous file...")
        chem_metadata = pd.read_csv(
            os.path.join(args.output_dir, "chem_metadata.tsv"), sep="\t"
        )
    tqdm.write("Done!")

    # -----------------
    # Sample Data
    # -----------------
    tqdm.write("Getting sample data...")
    sample_files_list = args.sample_files.split(",")
    chem_sample = build_sample_data(
        sample_files_list, loader, chem_metadata, args.sample_id_file, args.sample_map
    )
    tqdm.write("Done!")

    # -----------------
    # Endpoint Metadata
    # -----------------
    tqdm.write("Getting endpoint details...")
    endpoint_metadata = load_figshare_url(
        loader, args.endpoint_mapping_file, sheet_name=3
    )
    endpoint_metadata = format_endpoint_metadata(endpoint_metadata)
    tqdm.write("Done!")

    # else:
    tqdm.write("Saving data...")
    chem_metadata.fillna("NULL").to_csv(
        os.path.join(args.output_dir, "chemicals.csv"), index=False, quotechar='"'
    )
    chem_sample[SAMPLE_COLS].drop_duplicates().fillna("NULL").to_csv(
        os.path.join(args.output_dir, "samples.csv"), index=False, quotechar='"'
    )
    chem_sample[SAMP2CHEM_COLS].drop_duplicates().fillna("NULL").to_csv(
        os.path.join(args.output_dir, "samplesToChemicals.csv"),
        index=False,
        quotechar='"',
    )
    tqdm.write("Done!")


if __name__ == "__main__":
    main()
