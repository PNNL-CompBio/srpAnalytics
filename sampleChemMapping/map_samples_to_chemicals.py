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
from typing import Optional

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.data import FigshareDataLoader
from src.format import rename_duplicates
from src.mapping import rename_chemical_class
from src.metadata import build_chem_metadata, get_endpoint_metadata
from src.params import MASV_CC, MASV_SOURCE
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
OUTPUT_DIR = "/tmp"
# OUTPUT_DIR = "."

# Set CompTox API key
CTX_API_KEY = "5aded20c-9485-11ef-87c3-325096b39f47"

# Figshare file loader
loader = FigshareDataLoader(
    Path(OUTPUT_DIR) / ".figshare_cache",  # api_token=FIGSHARE_API_TOKEN
)


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


def load_clean_fses(filename: str) -> pd.DataFrame:
    """Load and clean FSES input files.

    Steps:
        1. Remove invalid and null entries
        2. Format location data and enforces negative longitudes
        3. Remove entries with no data

    Parameters
    ----------
    filename : str
        Path to FSES file

    Returns
    -------
    pd.DataFrame
        Cleaned FSES data
    """
    # Parse CSV files and HTTPS (figshare) separately
    if os.path.splitext(filename)[1] == ".csv":
        # Replace invalid values with nulls for filtering
        fses = pd.read_csv(filename, dtype={"Sample_ID": str, "Chemical_ID": str})
    elif os.path.splitext(filename)[1] == "":
        file_id = filename.split("/")[-1]
        _ = loader.load_data(file_id)
        fname = loader.get_file_path(file_id).as_posix()
        fses = pd.read_csv(fname)

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
    chem_metadata: pd.DataFrame,
    sample_id_file: str,
    sample_mapping: str = None,
):
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
    # Read and process all FSES files
    data = pd.concat([load_clean_fses(f) for f in fses_files], ignore_index=True)

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


def combine_chemical_endpoint_data(
    bmd_files: list[str],
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_details: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Combine chemical endpoint data.

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
    pd.DataFrame
        _description_
    """
    tqdm.write(f"Combining bmd files: {', '.join(bmd_files)}")

    # Read and concatenate the specified columns from all BMD files
    cols = get_cols_from_schema("BMDs")
    files = [pd.read_csv(file)[cols] for file in bmd_files]
    df = pd.concat(files)

    # Remove duplicates
    df = df.drop_duplicates(subset=["Chemical_ID", "End_Point"])

    # For extracts, create unique sample IDs
    if is_extract:
        # Split `sample_id` column on '-', take first two parts, and rename split cols
        sd_samp = chem_data.copy()
        split_result = sd_samp["Sample_ID"].str.split("-", expand=True)
        sd_samp["tmp_id"] = split_result[0]
        sd_samp = sd_samp.dropna(subset=["tmp_id"])

        # Prepare for join
        full_bmd = df.copy()
        full_bmd["tmp_id"] = full_bmd["Chemical_ID"].astype(str)
        full_bmd = full_bmd.drop(columns=["Chemical_ID"]).dropna(subset=["tmp_id"])
        full_bmd = pd.merge(
            full_bmd,
            sd_samp[["Sample_ID", "tmp_id"]].drop_duplicates(),
            on="tmp_id",
            how="outer",
        )

        # Fill missing sample IDs with tmp_id
        nas = full_bmd["Sample_ID"].isna()
        full_bmd.loc[nas, "Sample_ID"] = full_bmd.loc[nas, "tmp_id"]

        # Format sample names
        new_nas = (
            full_bmd["SampleName"].isna()
            if "SampleName" in full_bmd
            else pd.Series(True, index=full_bmd.index)
        )
        if any(new_nas):
            full_bmd.loc[new_nas, "SampleName"] = "Sample " + full_bmd.loc[
                new_nas, "Sample_ID"
            ].astype(str)

        # For extracts:
        # 1. Missing endpoint values -> "NoData"
        # 2. Merge endpoint details into dataframe
        # 3. Drop unused column (not needed after merge)
        # 4. Remove rows with missing sample ID
        # 5. Remove duplicate rows
        # 6. Missing location name -> "None"
        full_bmd = full_bmd.fillna({"End_Point": "NoData"})
        full_bmd = full_bmd.merge(endpoint_details, on="End_Point", how="right")
        full_bmd = full_bmd.drop(columns=["End_Point", "tmp_id"])
        full_bmd = full_bmd[full_bmd["Sample_ID"].notna()]
        full_bmd = full_bmd.drop_duplicates()
        full_bmd = full_bmd.fillna({"LocationName": "None"})

    # For non-extracts:
    # 1. Missing endpoint values -> "NoData"
    # 2. Merge endpoint details into dataframe
    # 3. Drop unused column (not needed after merge)
    # 4. Remove rows with missing CAS RN
    # 5. Remove duplicate rows
    # 6. Missing chemical class -> "Unclassified"
    else:
        full_bmd = df.copy()
        full_bmd = full_bmd.fillna({"End_Point": "NoData"})
        full_bmd = full_bmd.merge(endpoint_details, on="End_Point", how="right")
        full_bmd = full_bmd.drop(columns=["End_Point"])
        full_bmd = full_bmd.drop_duplicates()
        full_bmd = full_bmd.fillna({"chemical_class": "Unclassified"})

    # Handle QC flags
    full_bmd = full_bmd.rename(columns={"DataQC_Flag": "qc_num"})
    full_bmd["DataQC_Flag"] = full_bmd["qc_num"].apply(
        lambda x: "Poor" if x in [0, 1] else "Moderate" if x in [4, 5] else "Good"
    )
    full_bmd["Model"] = full_bmd["Model"].str.replace("NULL", "None")
    full_bmd = full_bmd.drop(columns=["qc_num"])

    return full_bmd


def combine_chemical_data(
    bmd_files: list[str],
    data_type: str = "Fits",
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_details: Optional[pd.DataFrame] = None,
):
    """_summary_

    Parameters
    ----------
    bmd_files : list[str]
        List of BMD files
    data_type : str, optional
        Type of data, by default "Fits"
        Options are: "Dose"/"dose" or "Fits"/"fits"/"fit"
    is_extract : bool, optional
        If True, applies to extracts (sample) data, by default False
    chem_data : Optional[pd.DataFrame], optional
        _description_, by default None
    endpoint_details : Optional[pd.DataFrame], optional
        _description_, by default None

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    ValueError
        _description_
    """
    # Determine columns based on data type
    if "fit" in data_type.lower():
        data_type = "Fits"
    elif "dose" in data_type.lower():
        data_type = "Dose"
    else:
        raise ValueError(
            f"Invalid data_type {data_type} " "(must contain 'fit' or 'dose)."
        )

    # cols = REQUIRED_BMD_COLUMNS[col_type]
    cols = get_cols_from_schema(data_type)

    # Read all files into df (and throw error for invalid)
    files = list()
    for file in bmd_files:
        try:
            df = pd.read_csv(file)[cols]
            files.append(df)
        except Exception as e:
            tqdm.write(f"Error reading {file}: {e}")
    if not files:
        tqdm.write("No valid files found")
        return pd.DataFrame()
    df = pd.concat(files, ignore_index=True)

    # Create combined identifier
    df["combined"] = df["Chemical_ID"].astype(str) + " " + df["End_Point"]

    # Process fitness data
    if data_type == "fit":
        df = df[df["X_vals"] != "NULL"]
        df["X_vals"] = pd.to_numeric(df["X_vals"])
        df["Y_vals"] = pd.to_numeric(df["Y_vals"])

    df = df.drop_duplicates(subset=["combined"], keep="last")

    # Combine endpoint metadata with chemical data
    df = df.merge(endpoint_details, on="End_Point", how="right")
    df = df.drop(columns=["End_Point", "Description"]).drop_duplicates()

    # Process for extracts if needed
    if is_extract:
        # Use temp (split) IDs for joining
        tmp_chem_data = chem_data.copy()
        split_result = tmp_chem_data["Sample_ID"].str.split("-", expand=True)
        tmp_chem_data["tmp_id"] = split_result[0]
        tmp_chem_data = tmp_chem_data[["Sample_ID", "tmp_id"]].drop_duplicates()
        df["tmp_id"] = df["Chemical_ID"].astype(str)
        df = df.drop(columns=["Chemical_ID"])
        df = df.merge(tmp_chem_data, on="tmp_id", how="left")

        # Fill missing sample IDs with temp ID
        mask = df["Sample_ID"].isna()
        df.loc[mask, "Sample_ID"] = df.loc[mask, "tmp_id"]

        # Delete temp ID column
        df = df.drop(columns=["tmp_id"])

    # For non-extracts, keep only chemicals that are in sample data
    elif not is_extract and chem_data is not None:
        df = df[df["Chemical_ID"].isin(chem_data["Chemical_ID"])]

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
        default=OUTPUT_DIR,
        help="File that maps sample locations",
    )

    args = parser.parse_args()

    tqdm.write("Getting chemical metadata...")
    chem_class = masv_chem_class(args.chem_class_file)
    if not os.path.exists(os.path.join(args.output_dir, "chem_metadata.tsv")):
        tqdm.write("No metadata found. Building metadata...")
        chem_metadata = build_chem_metadata(
            args.metadata, save_to=os.path.join(args.output_dir, "chem_metadata.tsv")
        )
    else:
        tqdm.write("Metadata found! Reading from previous file...")
        chem_metadata = pd.read_csv(
            os.path.join(args.output_dir, "chem_metadata.tsv"), sep="\t"
        )
    tqdm.write("Done!")

    tqdm.write("Getting sample data...")
    sample_files_list = args.sample_files.split(",")
    chem_sample = build_sample_data(
        sample_files_list, chem_metadata, args.sample_id_file, args.sample_map
    )
    tqdm.write("Done!")

    tqdm.write("Getting endpoint details...")
    endpoint_details = get_endpoint_metadata(
        args.endpoint_mapping_file
    ).drop_duplicates()
    tqdm.write("Done!")

    if args.is_sample or args.is_chem:
        all_files = args.dose_response.split(",")
        bmd_files = [file for file in all_files if "BMDs" in file]
        dose_files = [file for file in all_files if "Dose" in file]
        fit_files = [file for file in all_files if "Fits" in file]

        chem_data = chem_sample if args.is_sample else chem_metadata

        tqdm.write("Combining chemical endpoint data...")
        bmds = (
            combine_chemical_endpoint_data(
                bmd_files,
                is_extract=args.is_sample,
                chem_data=chem_data,
                endpoint_details=endpoint_details,
            )
            .dropna(subset=["BMD_Analysis_Flag"])
            .query("BMD_Analysis_Flag != 'NA'")
        )
        tqdm.write("Done!")

        tqdm.write("Combining chemical fitness data...")
        curves = combine_chemical_data(
            fit_files,
            data_type="fit",
            is_extract=args.is_sample,
            chem_data=chem_data,
            endpoint_details=endpoint_details,
        )
        tqdm.write("Done!")

        tqdm.write("Combining chemical dose response data...")
        dose_reps = combine_chemical_data(
            dose_files,
            data_type="dose",
            is_extract=args.is_sample,
            chem_data=chem_data,
            endpoint_details=endpoint_details,
        ).dropna(subset=["Dose"])
        tqdm.write("Done!")

        if args.is_sample:
            tqdm.write("Saving sample data to CSV...")
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
            tqdm.write("Done!")

        if args.is_chem:
            tqdm.write("Removing invalid chemical IDs...")
            nas = bmds[bmds["Chemical_ID"].isna()]["Chemical_ID"]
            to_remove = set(nas) - set(chem_sample["Chemical_ID"])
            bmds = bmds[~bmds["Chemical_ID"].isin(to_remove)]
            curves = curves[~curves["Chemical_ID"].isin(to_remove)]
            dose_reps = dose_reps[~dose_reps["Chemical_ID"].isin(to_remove)]
            tqdm.write("Done!")

            tqdm.write("Saving chemical data to CSV...")
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
            tqdm.write("Done!")

    else:
        tqdm.write("Saving data...")
        chem_metadata.to_csv(
            os.path.join(args.output_dir, "chemicals.csv"), index=False, quotechar='"'
        )
        chem_sample[SAMPLE_COLS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "samples.csv"), index=False, quotechar='"'
        )
        chem_sample[SAMP2CHEM_COLS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "samplesToChemicals.csv"),
            index=False,
            quotechar='"',
        )
        tqdm.write("Done!")


if __name__ == "__main__":
    main()
