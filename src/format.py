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
    snake_case_name = snake_case_name.replace("__", "_")

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
    samples = data[["Sample_ID", col]].drop_duplicates()
    samples["is_duplicate"] = samples[col].duplicated(keep=False)
    duplicates = samples[samples["is_duplicate"]]
    print("len duplicates", len(duplicates))

    if len(duplicates) == 0:
        return data

    # If duplicate, rename as name:01, etc., otherwise keep name
    new_names = {}
    for sname in duplicates[col].unique():
        sample_ids = duplicates[duplicates[col] == sname]["Sample_ID"].tolist()
        for i, sid in enumerate(sample_ids, 1):
            new_names[sid] = f"{sname}:{i:02d}"

    # Merge back into original data
    sid2sname = {
        sid: (new_names[sid] if sid in new_names.keys() else sname)
        for sid, sname in zip(samples["Sample_ID"], samples[col])
    }
    return [sid2sname[sid] for sid in data["Sample_ID"]]
