# =========================================================
# Imports
# =========================================================

import pandas as pd

from .format import snakeify
from numpy.typing import ArrayLike

# #################################
# Master ID tables
#
# Note: The database requires Sample_ID and Chemical_ID be unique.
# They are in some files but not others, so tables are
# automatically updated below.
# #################################


def chem_id_master_table(df: pd.DataFrame, cas_ids: ArrayLike) -> pd.DataFrame:
    """Generate master table for chemical ID

    Note: database requires Sample_ID and  Chemical_ID be unique. They are in some files but not others
    Thus, the tables below are automatically updated

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of chemical ID information (map_type="chemID")
    cas_ids : ArrayLike
        CAS IDs

    Returns
    -------
    pd.DataFrame
        Chemical ID master table
    """
    # Clean up input data and reduce # columns
    cols = ["cas_number", "zf.cid", "Chemical_ID", "chemical_class"]
    df = df[cols].drop_duplicates().dropna(subset=["cas_number"])

    # For compounds missing chemical IDs, assign new chem ID
    missing = set(cas_ids) - set(df["cas_number"])
    if missing:
        print(f"Missing {len(missing)} chemical IDs; adding them now...")
        max_id = df["Chemical_ID"].max() + 1

        # Create table entries for missing chemical IDs
        missing_df = pd.DataFrame(
            {
                "cas_number": list(missing),
                "zf.cid": [""] * len(missing),
                "Chemical_ID": range(max_id, max_id + len(missing)),
                "chemical_class": [""] * len(missing),
            }
        )

        # Append missing rows and format output col names
        df = pd.concat([df, missing_df])
        df.columns = [snakeify(c) for c in df.columns]

        # TODO: how do we update with new ids?
        # write.csv(newMap,paste0(data.dir,'chemicalIdMapping.csv'),row.names = F)
    return df


def sample_id_master_table(
    existing_sample_numbers: ArrayLike, smap: str
) -> pd.DataFrame:
    """Generate master table for sample ID

    Note: database requires Sample_ID and  Chemical_ID be unique. They are in some files but not others
    Thus, the tables below are automatically updated

    Parameters
    ----------
    existing_sample_numbers : ArrayLike
        _description_
    smap : str
        _description_

    Returns
    -------
    pd.DataFrame
        Sample ID master table
    """
    map_df = pd.read_csv(smap).loc[:, ["Sample_ID", "SampleNumber"]].drop_duplicates()

    missing = set(existing_sample_numbers) - set(map_df["SampleNumber"])
    if missing:
        print(f"Missing {len(missing)} sample IDs; adding them now...")
        max_id = map_df["Sample_ID"].astype(float).max(skipna=True) + 1
        missing_df = pd.DataFrame(
            {
                "Sample_ID": range(int(max_id), int(max_id) + len(missing)),
                "SampleNumber": list(missing),
            }
        )
        map_df = pd.concat([map_df, missing_df])

    return map_df


# ##################################
# Duplicate Removal
#
# There are two causes of duplicates in the data
# 1) Samples that are evaluated from multiple files:
#    Duplicated samples need to be removed.
# 2) Multiple chemical ids mapping to a single CAS ID:
#    A single chemical ID must be selected.
# ##################################


def remove_chem_id_duplicates(
    df: pd.DataFrame,
    chem_ids: pd.DataFrame,
    cols: ArrayLike = ["AUC_Norm", "End_Point_Name", "Model"],
) -> pd.DataFrame:
    """Remove duplicates from anything with Chemical_ID in the field

    Parameters
    ----------
    df : pd.DataFrame
        Table to be de-duplicated
    chem_ids : pd.DataFrame
        Mapping of chemical_IDs to CAS ids
    cols : ArrayLike, optional
        List of columns that make the data distinct,
            by default ["AUC_Norm", "End_Point_Name", "Model"]

    Returns
    -------
    pd.DataFrame
        De-duplicated table
    """
    group_by_cas = chem_ids.groupby("cas_number")["Chemical_ID"].nunique()
    dupes = group_by_cas[group_by_cas > 1].index

    chem_ids = (
        chem_ids[chem_ids["cas_number"].isin(dupes)]
        .loc[:, ["cas_number", "Chemical_ID"]]
        .assign(in_dataset=lambda df: df["Chemical_ID"].isin(df["Chemical_ID"]))
    )

    dist_tab = df.merge(chem_ids, on="Chemical_ID", how="left").drop_duplicates(
        subset=["cas_number"] + cols
    )
    dupe_counts = (
        dist_tab.groupby(["cas_number", "End_Point_Name", "Model"])
        .size()
        .reset_index(name="n_ids")
    )

    dupe_cas = dupe_counts[dupe_counts["n_ids"] > 1]["cas_number"].unique()
    return dupe_cas