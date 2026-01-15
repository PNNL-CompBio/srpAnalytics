"""schema.py: Functionality related to parsing LinkML schemas.

author(s): @christinehc
"""

# =========================================================
# Imports
# =========================================================
from os.path import exists

from linkml_runtime.utils.schemaview import SchemaView

# =========================================================
# Setup and Parameters
# =========================================================
SCHEMA_FILE = "../srpAnalytics.yaml"
ZEBRAFISH_DTYPE_TO_SUFFIX = {"bmd": "BMDs", "dose": "Dose", "fit": "Fits"}


# =========================================================
# Functions
# =========================================================
def map_zebrafish_data_to_schema(sample_type: str, data_type: str):
    stype = {"chemical": "Chem", "extract": "Samp"}
    return f"zebrafish{stype[sample_type]}{ZEBRAFISH_DTYPE_TO_SUFFIX[data_type]}"


def get_slots_from_schema(classname: str, filename: str = SCHEMA_FILE) -> list[str]:
    """Get all slots for a given class in schema.

    Note: Does not include attributes.

    Parameters
    ----------
    classname : str
        Name of class
    filename : str, optional
        Path to schema file, by default 'srpAnalytics.yaml'

    Returns
    -------
    list[str]
        List of class slots

    Raises
    ------
    FileNotFoundError
        If `filename` is not found.
    """
    if not exists(filename):
        raise FileNotFoundError(f"Schema file '{filename}' was not found.")

    view = SchemaView(filename)
    return view.get_class(classname).slots


def get_cols_from_schema(classname: str, filename: str = SCHEMA_FILE) -> list[str]:
    """Get columns for a given class in schema.

    Note: Combines all slots and attributes for a given class
        to give a full list of fields for the class in question.

    Parameters
    ----------
    classname : str
        Name of class
    filename : str, optional
        Path to schema file, by default 'srpAnalytics.yaml'

    Returns
    -------
    list[str]
        List of column names

    Raises
    ------
    FileNotFoundError
        If `filename` is not found.
    """
    if not exists(filename):
        raise FileNotFoundError(f"Schema file '{filename}' was not found.")

    view = SchemaView(filename)
    return view.class_slots(classname)


def combine_schema_cols(
    classname_a: str, classname_b: str, filename: str = SCHEMA_FILE
) -> list[str]:
    """Get combined columns for 2 classes in a schema.

    Parameters
    ----------
    classname_a : str
        Name of first class
    classname_b : str
        Name of second class
    filename : str, optional
        PPath to schema file, by default 'srpAnalytics.yaml'

    Returns
    -------
    list[str]
        List of column names
    """
    cols_a = get_cols_from_schema(classname_a, filename=filename)
    cols_b = get_cols_from_schema(classname_b, filename=filename)
    return list(dict.fromkeys(cols_a + cols_b))
