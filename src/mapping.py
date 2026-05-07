"""mapping.py: Functions to handle mapping operations.

NOTE: This may be superseded/replaced by manifest.py in future updates

author(s): @christinehc
"""

# =========================================================
# Imports
# =========================================================
import pandas as pd


# =========================================================
# Functions
# =========================================================
# @FIXME: Possibly delete? Seems unused
def load_mapping_reference(
    data_dir: str = "https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data",
    filename: str = "srp_build_files.csv",
) -> pd.DataFrame:
    """Get online data files
    Note: most data files are stored on github and can be edited there, then this will be rerun

    Parameters
    ----------
    data_dir : _type_, optional
        _description_, by default "https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data"
    filename : str, optional
        _description_, by default "srp_build_files.csv"

    Returns
    -------
    pd.DataFrame
        _description_
    """
    return pd.read_csv(f"{data_dir}/{filename}")


# @FIXME: Possibly delete? Seems unused
def get_mapping_file(df: pd.DataFrame, map_type: str, return_first: bool = True) -> str:
    """Get mapping filename from master list.

    Parameters
    ----------
    df : pd.DataFrame
        Mapping file reference table
    map_type : str
        Type of mapping file
    return_first : bool, optional
        If True, returns only first result.
        If False, returns comma-separated list of results.
            by default True

    Returns
    -------
    pd.DataFrame
        /path/to/map_type/file
    """
    if return_first:
        return list(df.loc[df.name == map_type].location)[0]
    return ",".join(list(df.loc[df.data_type == map_type].location))


# Map the newClass values based on conditions provided
def rename_chemical_class(x):
    if x in [
        "industrial",
        "industrial; aniline",
        "Industrial",
        "industrial; consumerProduct; phenol",
        "industrial; consumerProduct; aniline",
        "industrial; phenol",
    ]:
        return "Industrial"
    elif x in ["PAH; industrial", "PAH"]:
        return "PAH"
    elif x in [
        "personalCare; personalCare; natural; natural; consumerProduct; consumerProduct",
        "personalCare; natural; consumerProduct",
        "personalCare; natural",
        "pharmacological; personalCare; industrial; natural; consumerProduct",
    ]:
        return "Natural"
    elif x == "pestFungicide":
        return "Fungicide"
    elif pd.isna(x) or x == "NA":
        return "Unclassified"
    else:
        return x
