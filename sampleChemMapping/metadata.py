# =========================================================
# Imports
# =========================================================


import ctxpy as ctx
import numpy as np
import pandas as pd

from .mapping import get_mapping_file
from .format import format_cas, snakeify
from .tables import chem_id_master_table

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


def query_comptox_by_cas(
    df: pd.DataFrame,
    keep_cols: list[str] = ["cas_number", "chemical_id", "chemical_class"],
    data_cols: list[str] = ["preferredName", "smiles", "dtxsid", "dtxcid"],
) -> pd.DataFrame:
    """Perform CompTox API search function by CAS.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing chemical data
    keep_cols : list[str], optional
        Columns from input dataframe to keep,
            by default ["cas_number", "chemical_id", "chemical_class"]
    data_cols : list[str], optional
        Fields from CompTox search to keep as data columns,
            by default ["preferredName", "smiles", "dtxsid", "dtxcid"]

    Returns
    -------
    pd.DataFrame
        DataFrame with selected chemical data and CompTox results
    """
    # Initialize CompTox API
    chem = ctx.Chemical(x_api_key=CTX_API_KEY)

    data = {c.lower(): [] for c in keep_cols}
    comptox_fields = {snakeify(c): [] for c in data_cols}
    data = {**data, **comptox_fields}

    for i, cas in enumerate(df["cas_number"].values):
        # Fix CAS erroneously converted to Excel dates
        cas = format_cas(cas)

        data["cas_number"] += [cas]
        for col in keep_cols:
            data[col.lower()] += [df[col].values[i]]

        try:
            details = chem.search(by="equals", word=cas)[0]

            for c in data_cols:
                data[snakeify(c)] += [details[c]]

        except (KeyError, TypeError):
            for c in data_cols:
                data[snakeify(c)] += [np.nan]

    data = pd.DataFrame(data)
    return data


def query_comptox_by_dtxsid(
    df: pd.DataFrame, data_cols: list[str] = ["averageMass", "inchikey", "molFormula"]
) -> pd.DataFrame:
    """Perform CompTox API detail query by DTXSID.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing DTXSIDs
    data_cols : list[str], optional
        Fields from CompTox details to keep as data column,
            by default ["averageMass", "inchikey", "molFormula"]

    Returns
    -------
    pd.DataFrame
        DataFrame with desired CompTox results

    Raises
    ------
    ValueError
        Raised if DTXSID fails the CTX API search
    """
    # Initialize CompTox API
    chem = ctx.Chemical(x_api_key=CTX_API_KEY)

    result = {snakeify(c): [] for c in data_cols}

    for sid in df["dtxsid"]:
        if isinstance(sid, str):
            try:
                details = chem.details(by="dtxsid", word=sid)
                for c in data_cols:
                    result[snakeify(c)] += [details[c]]
            except:
                raise ValueError(f"Encountered CTX error with DTXSID {sid}.")

        else:
            for c in data_cols:
                result[snakeify(c)] += [np.nan]

    result = pd.DataFrame(result)
    return pd.concat([df, result], axis=1)


def query_comptox(
    df: pd.DataFrame,
    keep_cols: list[str] = ["cas_number", "chemical_id", "chemical_class"],
    cas_data_cols: list[str] = ["preferredName", "smiles", "dtxsid", "dtxcid"],
    dtxsid_data_cols: list[str] = ["averageMass", "inchikey", "molFormula"],
) -> pd.DataFrame:
    """Get CompTox API results from CAS + DTXSID queries.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing chemical data
    keep_cols : list[str], optional
        Columns from input dataframe to keep,
            by default ["cas_number", "chemical_id", "chemical_class"]
    cas_data_cols : list[str], optional
        Fields from CompTox CAS search to keep as data columns,
            by default ["preferredName", "smiles", "dtxsid", "dtxcid"]
    dtxsid_data_cols : list[str], optional
        Fields from CompTox DTXSID details to keep as data column,
            by default ["averageMass", "inchikey", "molFormula"]

    Returns
    -------
    pd.DataFrame
        DataFrame with desired CompTox results
    """
    df = query_comptox_by_cas(df, keep_cols=keep_cols, data_cols=cas_data_cols)
    df = query_comptox_by_dtxsid(df, data_cols=dtxsid_data_cols)

    return df


def get_chem_metadata(
    filename: str,
    keep_cols: list[str] = ["cas_number", "chemical_id", "chemical_class"],
    cas_data_cols: list[str] = ["preferredName", "smiles", "dtxsid", "dtxcid"],
    dtxsid_data_cols: list[str] = ["averageMass", "inchikey", "molFormula"],
) -> pd.DataFrame:
    """Get chemical metadata, which is stored in `data.dir`

    Parameters
    ----------
    filename : str
        Path to master table of data files and types
            (i.e. srp_build_files.csv)
    keep_cols : list[str], optional
        Columns to keep from chem ID mapping file,
            by default ["cas_number", "chemical_id", "chemical_class"]
    cas_data_cols : list[str], optional
        Columns to keep from CAS search,
            by default ["preferredName", "smiles", "dtxsid", "dtxcid"]
    dtxsid_data_cols : list[str], optional
        Columns to keep from detailed CAS search,
            by default ["averageMass", "inchikey", "molFormula"]

    Returns
    -------
    pd.DataFrame
        _description_
    """
    # Get chem ID file
    mappings = pd.read_csv(filename)
    chem_ids = get_mapping_file(mappings, "chemID")

    # Get chemical metadata from CompTox dashboard
    metadata = (
        query_comptox(
            chem_ids,
            keep_cols=keep_cols,
            cas_data_cols=cas_data_cols,
            dtxsid_data_cols=dtxsid_data_cols,
        )
        .replace(["NA", "N/A"], np.nan)
        .dropna(subset=["cas_number"])
    )

    # Create table linking CAS to Chem IDs; drop duplicate entries
    chem_ids = chem_id_master_table(chem_ids, metadata["cas_number"])
    chem_ids = chem_ids[["cas_number", "chemical_id"]].drop_duplicates()

    # Add chem metadata information to Chem ID table
    metadata = metadata.merge(chem_ids, on="cas_number", how="left")

    # Clean up metadata table
    metadata["PREFERRED_NAME"] = metadata["PREFERRED_NAME"].str.replace(
        r"^-$", "Chemical name unknown", regex=True
    )  # replace "-" entries with more informative blurb
    metadata = metadata.fillna(
        {"newClass": "Unclassified", "PREFERRED_NAME": "Chemical name unknown"}
    ).drop(columns=["ParameterName"])  # fill unknowns with NaN

    # If no CAS, remove chemicals from final table
    nocas = metadata["cas_number"].str.contains("NOCAS", na=False)
    if nocas.any():
        print(f"removing {nocas.sum()} chems with no cas")
        metadata = metadata[~nocas]
    return metadata


def get_endpoint_metadata(filename: str) -> pd.DataFrame:
    """_summary_

    Parameters
    ----------
    filename : str
        Path to endpoint file

    Returns
    -------
    pd.DataFrame
        _description_
    """
    endpoint_details = (
        pd.read_excel(filename, sheet_name=3)
        .rename(
            columns={
                "Abbreviation": "End_Point",
                "Simple name (<20char)": "End_Point_Name",
                "Ontology Link": "endPointLink",
            }
        )
        .assign(
            IncludeInPortal="No",
            End_Point="NoData",
            End_Point_Name=None,
            Description="No data",
            endPointLink="",
        )
        .pipe(
            lambda df: df.append(
                {
                    "IncludeInPortal": "No",
                    "End_Point": "NoData",
                    "End_Point_Name": None,
                    "Description": "No data",
                    "endPointLink": "",
                },
                ignore_index=True,
            )
        )
    )
    return endpoint_details
