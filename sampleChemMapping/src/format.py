import re
from pandas import DataFrame

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
