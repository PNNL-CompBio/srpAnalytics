import re
from numpy import nan
from pandas import DataFrame, read_csv, to_numeric
from .params import REQUIRED_SAMPLE_COLUMNS

RENAME = {"casrn": "cas_number"}


def format_cas(cas: str) -> str:
    if isinstance(cas, str):
        if "/" not in cas:
            return cas
        cas = cas.split("/")
        return f"{cas[2]}-{int(cas[0]):02d}-{cas[1]}"
    return cas


def snakeify(name: str) -> str:
    """
    Convert a camelCase string to snake_case.

    Written by AI Incubator

    Parameters:
        name (str): The camelCase string.

    Returns:
        str: The snake_case string.
    """
    # Convert camelCase to snake_case using regular expressions
    snake_case_name = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    snake_case_name = re.sub("([a-z0-9])([A-Z])", r"\1_\2", snake_case_name).lower()

    return snake_case_name


def snakeify_all_columns(df: DataFrame, rename: dict[str, str] = RENAME) -> DataFrame:
    """Convert all dataframe camelCase column names to snake_case.

    Parameters
    ----------
    df : DataFrame
        Input dataframe
    rename : dict[str, str], optional
        Column mappings for manual renaming, by default RENAME
            RENAME = {"casrn": "cas_number"}

    Returns
    -------
    DataFrame
        DataFrame with snake_case column names
    """
    df = df.rename(columns=RENAME)
    df.columns = [snakeify(c) for c in df.columns]
    return df


def chunker(seq: DataFrame, size: int) -> DataFrame:
    """Splits input data into chunks of specified size.

    Borrowed from https://stackoverflow.com/a/25701576

    Parameters
    ----------
    seq : DataFrame
        _description_
    size : int
        _description_

    Returns
    -------
    DataFrame
        _description_
    """
    return (seq.iloc[pos : pos + size] for pos in range(0, len(seq), size))


def process_fses(filename: str, snake_case: bool = True) -> DataFrame:
    # Replace invalid values with nulls for filtering
    df = read_csv(filename)[REQUIRED_SAMPLE_COLUMNS].replace(
        {"BLOD": "0", "NULL": "0", "nc:BDL": "0"}
    )
    if snake_case:
        df = snakeify_all_columns(df)

    # Remove null and invalid entries
    df = df[
        (df["sample_number"].notna())
        & (df["cas_number"].notna())
        & (~df["measurement_value"].isin(["0", nan]))
        & (~df["measurement_value_molar"].isin(["0", nan]))
    ]

    # Format FSES location data
    df["location_lon"] = to_numeric(df["location_lon"], errors="coerce")

    # Only allow negative longitudes
    # Note: Our data is already all negative, so unnecessary?
    # df["location_lon"] = np.where(
    #     df["location_lon"].gt(0), -df["location_lon"], df["location_lon"]
    # )
    return df


def rename_duplicates(data: DataFrame, col: str = "sample_name") -> list[str]:
    """Rename duplicate name entries with index, i.e. sample:01, etc.

    Parameters
    ----------
    data : DataFrame
        Table of samples
    col : str, optional
        Column containing duplicates, by default "sample_name"

    Returns
    -------
    list[str]
        Values of data[col] with sample:01, sample:02, etc.
            for duplicate entries.
            Note: non-duplicate values are unchanged.
    """
    is_duplicate = data[col].duplicated(keep=False)

    new_names, counts = list(), dict()

    # If duplicate, rename as name:01, etc., otherwise keep name
    for s, d in zip(data[col], is_duplicate):
        c = counts.get(s, 0) + 1
        counts[s] = c
        if d:
            new_names.append(f"{s}:{c:02d}")
        else:
            new_names.append(s)

    return new_names
