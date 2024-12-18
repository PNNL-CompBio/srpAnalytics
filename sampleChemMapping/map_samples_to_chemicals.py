"""map_samples_to_chemicals: Python version of mapSamplesToChems.R

Original rewrite using generative AI; modified by @christinehc
"""

# =========================================================
# Imports
# =========================================================
import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd

from numpy.typing import ArrayLike
from src.format import snakeify, snakeify_all_columns
from src.mapping import rename_chemical_class
from src.metadata import build_chem_metadata, get_endpoint_metadata
from src.params import (
    MASV_CC,
    MASV_SOURCE,
    REQUIRED_BMD_COLUMNS,
    REQUIRED_SAMPLE_COLUMNS,
    SAMPLE_CHEM_COLUMNS,
    SAMPLE_COLUMNS,
)
from src.tables import sample_id_master_table

# These pathways refer to absolute pathways in the docker image
# setting these three parameters, can be appended
# data_dir = 'https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data'
out_dir = "/tmp/"
# out_dir = "./"

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


#' buildSampleData - takes the curated information and selects the data we need
#' @param data.dir
#' @return data.frame
# fses_files, #files from barton that contain sample info
# chemMeta, #metadata for chemicals including identifier mapping
# sampIds, #new ids for samples
# sampMapping  ##mapping for sample names to clean up
def build_sample_data(
    fses_files: ArrayLike, chem_metadata: pd.DataFrame, sample_ids, samp_mapping=None
) -> pd.DataFrame:
    sample_data = []

    # Read and process each FSES file
    for f in fses_files:
        # Replace invalid values with nulls for filtering
        tmp = pd.read_csv(f)[REQUIRED_SAMPLE_COLUMNS].replace(
            {"BLOD": "0", "NULL": "0", "nc:BDL": "0"}
        )
        tmp = snakeify_all_columns(tmp)

        # Remove null and invalid entries
        tmp = tmp[
            (tmp["sample_number"].notna())
            & (tmp["cas_number"].notna())
            & (~tmp["measurement_value"].isin(["0", np.nan]))
            & (~tmp["measurement_value_molar"].isin(["0", np.nan]))
        ]

        # Format FSES location data
        tmp["location_lon"] = pd.to_numeric(tmp["location_lon"], errors="coerce")

        # Only allow negative longitudes
        # Note: Our data is already all negative, so unnecessary?
        # tmp["location_lon"] = np.where(
        #     tmp["location_lon"].gt(0), -tmp["location_lon"], tmp["location_lon"]
        # )
        sample_data.append(tmp)

    # Add metadata to sample chemical table
    sample_data = pd.concat(sample_data).merge(
        chem_metadata[["chemical_id", "cas_number", "average_mass"]],
        on="cas_number",
        how="left",
    )
    # Enforce snake_case for all col names
    sample_data = snakeify_all_columns(sample_data)

    # Add sample ID information to sample data table
    sample_data = sample_data.merge(
        sample_id_master_table(sample_data["sample_number"], sample_ids),
        on="sample_number",
        how="left",
    ).drop_duplicates()

    # Find all duplicate sample names
    all_samp_names = (
        sample_data[["sample_id", "sample_name"]]
        .drop_duplicates()
        .assign(is_dupe=lambda df: df.duplicated("sample_name"))
    )
    dupe_samp_names = all_samp_names.query("is_dupe").sort_values("sample_name")
    num_dupes = (
        dupe_samp_names.groupby("sample_name")
        .size()
        .reset_index(name="nid")
        .merge(dupe_samp_names)
    )

    new_names = (
        num_dupes[["sample_name", "nid"]]
        .drop_duplicates()
        .assign(
            new_name=lambda df: df.apply(
                lambda row: ":".join(
                    [f"{row['sample_name']}:{i+1}" for i in range(int(row["nid"]))]
                ),
                axis=1,
            )
        )
        .drop(columns="nid")
        .rename(columns={"SampleName": "old_sn"})
    )

    full_rep = pd.concat([dupe_samp_names, new_names], axis=1)[
        ["Sample_ID", "new_name"]
    ]
    new_samp_names = pd.concat(
        [
            all_samp_names.query("~is_dupe").drop("is_dupe", axis=1),
            full_rep.rename(columns={"new_name": "SampleName"}),
        ]
    )
    sample_data = (
        sample_data.drop("SampleName", axis=1)
        .merge(new_samp_names, on="Sample_ID", how="left")
        .copy()
    )

    sample_name_remap = (
        pd.read_excel(samp_mapping, sheet_name=0)
        .loc[:, ["Sample_ID", "ProjectName", "NewSampleName", "NewLocationName"]]
        .drop_duplicates()
    )

    sample_data = sample_data.merge(
        sample_name_remap, on="Sample_ID", how="left"
    ).drop_duplicates()

    nas = sample_data["projectName"].isna()
    sample_data.loc[nas, "projectName"] = sample_data.loc[nas, "ProjectName"]
    sample_data.loc[nas, "LocationName"] = sample_data.loc[nas, "NewLocationName"]
    sample_data.loc[nas, "SampleName"] = sample_data.loc[nas, "NewSampleName"]

    sample_data = sample_data.drop(
        columns=["ProjectName", "NewSampleName", "NewLocationName", "AVERAGE_MASS"]
    ).dropna(subset=["cas_number"])

    return sample_data


def combine_v2_chemical_endpoint_data(
    bmd_files, is_extract=False, samp_chem=None, endpoint_details=None
):
    print(f"Combining bmd files: {', '.join(bmd_files)}")
    cols = REQUIRED_BMD_COLUMNS["bmd"]
    mid_bmd = pd.concat([pd.read_csv(bmd_file).loc[:, cols] for bmd_file in bmd_files])

    dupes = mid_bmd[["Chemical_ID", "End_Point"]].duplicated().values
    if dupes.any():
        mid_bmd = mid_bmd[~dupes]

    if is_extract:
        sd_samp = (
            samp_chem["Sample_ID"]
            .str.split("-", expand=True)
            .iloc[:, 0:2]
            .rename(columns={0: "tmp_id", 1: "sub"})
            .merge(samp_chem[["Sample_ID"]].drop_duplicates(), on="sub", how="left")
        )

        full_bmd = (
            mid_bmd.assign(tmp_id=mid_bmd["Chemical_ID"].astype(str))
            .drop("Chemical_ID", axis=1)
            .merge(sd_samp, on="tmp_id", how="left")
        )

        nas = full_bmd["Sample_ID"].isna()
        full_bmd.loc[nas, "Sample_ID"] = full_bmd.loc[nas, "tmp_id"]

        new_nas = full_bmd["SampleName"].isna()
        if new_nas.any():
            full_bmd.loc[new_nas, "SampleName"] = (
                "Sample " + full_bmd.loc[new_nas, "Sample_ID"]
            )

        full_bmd = (
            full_bmd.fillna({"End_Point": "NoData"})
            .merge(endpoint_details, on="End_Point", how="right")
            .drop(columns=["End_Point", "tmp_id"])
            .dropna(subset=["Sample_ID"])
            .drop_duplicates()
            .fillna({"LocationName": "None"})
        )
    else:
        full_bmd = (
            mid_bmd.fillna({"End_Point": "NoData"})
            .merge(endpoint_details, on="End_Point", how="right")
            .drop(columns="End_Point")
            .dropna(subset=["cas_number"])
            .drop_duplicates()
            .fillna({"chemical_class": "Unclassified"})
        )

    full_bmd = (
        full_bmd.rename(columns={"DataQC_Flag": "qc_num"})
        .assign(
            DataQC_Flag=full_bmd["qc_num"].replace(
                {0: "Poor", 1: "Poor", 4: "Moderate", 5: "Moderate"}, regex=True
            )
        )
        .assign(Model=lambda df: df["Model"].str.replace("NULL", "None"))
        .drop(columns="qc_num")
    )

    return full_bmd


def combine_chemical_fit_data(
    bmd_files, is_extract=False, samp_chem=None, endpoint_details=None
):
    print(f"Combining fit files: {', '.join(bmd_files)}")
    cols = REQUIRED_BMD_COLUMNS["fitVals"]
    files = [pd.read_csv(bmd_file).loc[:, cols] for bmd_file in bmd_files]

    chem_eps = [
        set(
            file.assign(combined=lambda df: df["Chemical_ID"] + df["End_Point"])[
                "combined"
            ]
        )
        for file in files
    ]
    new_chem_eps = chem_eps.copy()

    if len(chem_eps) > 1:
        for i in range(1, len(new_chem_eps)):
            previous_eps = set.union(*chem_eps[:i])
            new_chem_eps[i] -= previous_eps

    fixed_files = []
    for i, (file, new_eps) in enumerate(zip(files, new_chem_eps)):
        fixed_files.append(
            file.assign(combined=lambda df: df["Chemical_ID"] + df["End_Point"])
            .query("combined in @new_eps")
            .replace({"X_vals": "NULL"}, np.nan)
            .assign(
                X_vals=lambda df: df["X_vals"].astype(float),
                Y_vals=lambda df: df["Y_vals"].astype(float),
            )
        )

    mid_bmd = pd.concat(fixed_files)
    full_bmd = (
        mid_bmd.merge(endpoint_details, on="End_Point", how="right")
        .drop_duplicates()
        .drop(columns=["End_Point", "Description"])
    )

    if is_extract:
        sd_samp = (
            samp_chem["Sample_ID"]
            .str.split("-", expand=True)
            .iloc[:, 0:2]
            .rename(columns={0: "tmp_id", 1: "sub"})
            .merge(samp_chem[["Sample_ID"]].drop_duplicates(), on="sub", how="left")
        )

        full_bmd = (
            full_bmd.assign(tmp_id=full_bmd["Chemical_ID"].astype(str))
            .drop("Chemical_ID", axis=1)
            .merge(sd_samp, on="tmp_id", how="left")
        )

        nas = full_bmd["Sample_ID"].isna()
        full_bmd.loc[nas, "Sample_ID"] = full_bmd.loc[nas, "tmp_id"]
        full_bmd = full_bmd.drop(columns=["tmp_id"])

    return full_bmd


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
    class_file,
    save_to: str = "MASV_classAndSource.csv",
    id_cols: ArrayLike = ["CASNumber", "ParameterName"],
):
    """Reads full MASV class annotations and assigns values to chemicals.

    Parameters
    ----------
    class_file : _type_
        _description_
    save_to : str, optional
        _description_, by default "MASV_classAndSource.csv"
    id_cols : ArrayLike, optional
        _description_, by default ["CASNumber", "ParameterName"]

    Returns
    -------
    _type_
        _description_
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


def combine_chemical_dose_data(
    bmd_files, is_extract=False, samp_chem=None, endpoint_details=None
):
    print(f'Combining dose response files: {", ".join(bmd_files)}')
    cols = REQUIRED_BMD_COLUMNS["doseRep"]
    files = [pd.read_csv(bmd_file).loc[:, cols] for bmd_file in bmd_files]

    chem_eps = [
        set(
            file.assign(combined=lambda df: df["Chemical_ID"] + df["End_Point"])[
                "combined"
            ]
        )
        for file in files
    ]
    new_chem_eps = chem_eps.copy()

    if len(chem_eps) > 1:
        for i in range(1, len(new_chem_eps)):
            previous_eps = set.union(*chem_eps[:i])
            new_chem_eps[i] -= previous_eps

    fixed_files = []
    for i, (file, new_eps) in enumerate(zip(files, new_chem_eps)):
        fixed_files.append(
            file.assign(combined=lambda df: df["Chemical_ID"] + df["End_Point"])
            .query("combined in @new_eps")
            .loc[:, cols]
        )

    mid_bmd = pd.concat(fixed_files)
    full_bmd = (
        mid_bmd.merge(endpoint_details, on="End_Point", how="right")
        .drop_duplicates()
        .drop(columns=["End_Point", "Description"])
    )

    if is_extract:
        sd_samp = (
            samp_chem["Sample_ID"]
            .str.split("-", expand=True)
            .iloc[:, 0:2]
            .rename(columns={0: "tmp_id", 1: "sub"})
            .merge(samp_chem[["Sample_ID"]].drop_duplicates(), on="sub", how="left")
        )

        full_bmd = (
            full_bmd.assign(tmp_id=full_bmd["Chemical_ID"].astype(str))
            .drop("Chemical_ID", axis=1)
            .merge(sd_samp, on="tmp_id", how="left")
            .assign(Sample_ID=lambda df: df["Sample_ID"].fillna(df["tmp_id"]))
        )

    return full_bmd


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
        dest="dose_res_stat",
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
        dest="samp_map",
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

    chem_class = masv_chem_class(args.chem_class_file)
    chem_meta = build_chem_metadata(args.metadata)

    sample_files_list = args.sample_files.split(",")
    samp_chem = build_sample_data(
        sample_files_list, chem_meta, args.sample_id_file, args.samp_map
    )

    endpoint_details = get_endpoint_metadata(
        args.endpoint_mapping_file
    ).drop_duplicates()

    if args.is_sample or args.is_chem:
        all_files = args.dose_res_stat.split(",")
        bmd_files = [file for file in all_files if "bmd" in file]
        dose_files = [file for file in all_files if "dose" in file]
        fit_files = [file for file in all_files if "fit" in file]

        meta_file = samp_chem if args.is_sample else chem_meta

        bmds = (
            combine_v2_chemical_endpoint_data(
                bmd_files,
                is_extract=args.is_sample,
                samp_chem=meta_file,
                endpoint_details=endpoint_details,
            )
            .dropna(subset=["BMD_Analysis_Flag"])
            .query("BMD_Analysis_Flag != 'NA'")
        )
        curves = combine_chemical_fit_data(
            fit_files,
            is_extract=args.is_sample,
            samp_chem=meta_file,
            endpoint_details=endpoint_details,
        )
        dose_reps = combine_chemical_dose_data(
            dose_files,
            is_extract=args.is_sample,
            samp_chem=meta_file,
            endpoint_details=endpoint_details,
        ).dropna(subset=["Dose"])

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
            nas = bmds["Chemical_ID"].isna()
            to_remove = set(nas) - set(samp_chem["Chemical_ID"])
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
        chem_meta.to_csv(
            os.path.join(args.output_dir, "chemicals.csv"), index=False, quotechar='"'
        )
        samp_chem[SAMPLE_COLUMNS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "samples.csv"), index=False, quotechar='"'
        )
        samp_chem[SAMPLE_CHEM_COLUMNS].drop_duplicates().to_csv(
            os.path.join(args.output_dir, "sampleToChemicals.csv"),
            index=False,
            quotechar='"',
        )


main()
