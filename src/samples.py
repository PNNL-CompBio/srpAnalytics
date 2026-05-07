"""samples.py: Handle zebrafish sample files.

author(s): @christinehc, @sgosline
"""

# =========================================================
# Imports
# =========================================================
from typing import Optional

import pandas as pd
from tqdm import tqdm

from .schema import get_cols_from_schema


# =========================================================
# Functions
# =========================================================
def combine_chemical_data(
    bmd_files: list[str],
    data_type: str = "Fits",
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_metadata: Optional[pd.DataFrame] = None,
):
    """Combine chemical data with endpoint metadata.

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
    endpoint_metadata : Optional[pd.DataFrame], optional
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

    # Remove End_Point_Name (currently added later)
    cols = get_cols_from_schema(data_type)

    # Remove End_Point_Name (currently added later)
    cols = [c for c in cols if c != "End_Point_Name"]

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
    df = df.merge(endpoint_metadata, on="End_Point", how="right")
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

        # Delete temp ID column and move sample ID to front
        df = df.drop(columns=["tmp_id"])
        cols = ["Sample_ID"] + [col for col in df.columns if col != "Sample_ID"]
        df = df[cols]

    # For non-extracts, keep only chemicals that are in sample data
    elif not is_extract and chem_data is not None:
        df = df[df["Chemical_ID"].isin(chem_data["Chemical_ID"])]

    return df.drop_duplicates()


def combine_chemical_endpoint_data(
    bmd_files: list[str],
    is_extract: bool = False,
    chem_data: Optional[pd.DataFrame] = None,
    endpoint_metadata: Optional[pd.DataFrame] = None,
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
    endpoint_metadata : Optional[pd.DataFrame], optional
        Tabulated endpoint data, by default None

    Returns
    -------
    pd.DataFrame
        Combined chemical endpoint data
    """
    tqdm.write(f"Combining bmd files: {', '.join(bmd_files)}")

    # Read and concatenate the specified columns from all BMD files
    cols = get_cols_from_schema("BMDs")
    cols = [
        c for c in cols if c != "End_Point_Name"
    ]  # Remove End_Point_Name (currently added later)
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
        full_bmd = full_bmd.merge(endpoint_metadata, on="End_Point", how="right")
        full_bmd = full_bmd.drop(columns=["End_Point", "tmp_id"])
        full_bmd = full_bmd[full_bmd["Sample_ID"].notna()]
        full_bmd = full_bmd.drop_duplicates()
        full_bmd = full_bmd.fillna({"LocationName": "None"})

        # Move sample ID to front
        cols = ["Sample_ID"] + [col for col in full_bmd.columns if col != "Sample_ID"]
        full_bmd = full_bmd[cols]

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
        full_bmd = full_bmd.merge(endpoint_metadata, on="End_Point", how="right")
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
