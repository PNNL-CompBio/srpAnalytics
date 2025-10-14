"""exposome_summary_stats: Build exposome data and get summary statistics."""

# ===============================================
#  IMPORTS
# ===============================================
import json
import re
import sys
from os.path import join

import pandas as pd
import requests

# ===============================================
#  CONFIG
# ===============================================
OUTPUT_DIR = "/tmp/"

PROJ2NAME = {
    "ADIPO": "Human adipocyte cell lines",
    "HEPG2": "Human Hepg2 cell lines",
    "MCF10A": "Human MCF10A cell lines",
    "TG-GATEs": "Human TG-GATEs",
}


# ===============================================
#  FUNCTIONS
# ===============================================
def check_args():  # Check for the path to the chemical mapping file
    """Enforce that the minimum number of arguments is called."""
    if len(sys.argv) < 2:
        print("Need to call script with path to chemical mapping file.")
        sys.exit()


def _load_projects(
    url: str = "https://montilab.bu.edu/Xposome-API/projects?all=Yes",
    verbose: bool = True,
) -> list[str]:
    """Load list of projects from exposome API.

    Parameters
    ----------
    url : str, optional
        URL to projects, by default "https://montilab.bu.edu/Xposome-API/projects?all=Yes"
    verbose : bool, optional
        If True, prints stdout, by default True

    Returns
    -------
    list[str]
        List of projects

    Raises
    -------
    HTTPError
        If invalid URL, terminates.
    """
    res = requests.get(url)
    res.raise_for_status()  # Raises HTTPError for invalid
    projects = json.loads(res.json()[0])

    if verbose:
        print(f"Data from {len(projects)} projects found.")

    return projects


def _load_chemicals(project: str) -> pd.DataFrame:
    """Load chemicals associated with project.

    Parameters
    ----------
    project : str
        Project (e.g. ADIPO, HEPG2, etc.)

    Returns
    -------
    pd.DataFrame
        Tabulated chemical information with the following columns:
            Chemical_ID: str
            Project: str
            Chemical_Name: str
            CAS: str
    """
    url = f"https://montilab.bu.edu/Xposome-API/chemicals?projects={project}&chemical_ids=all"
    res = requests.get(url)
    if res.status_code == 200:
        data = pd.DataFrame(json.loads(res.json()[0]))
        data = data.rename(columns={"Chemical_Id": "Chemical_ID"})
        return data.dropna(subset=["Chemical_ID"])
    return


def format_concentration(
    df: pd.DataFrame, column_name: str = "Concentration"
) -> pd.DataFrame:
    """Split concentration column into separate numeric and unit columns,
    only recognizing units 'uM' or 'mg/kg' from a larger string.

    Written using AI incubator 01/2025

    Parameters
    ----------
    df : pandas.DataFrame
        Data containing the concentration column.
    column_name : str
        Name of the column to be split (default is 'Concentration').

    Returns
    -------
    pandas.DataFrame
        Data with separate 'Concentration' and 'Unit' columns.
    """

    def extract_unit_and_numeric(c):
        unit_match = re.search(r"(\d*\.?\d+)\s*(uM|mg/kg)", c)
        if unit_match:
            return float(unit_match.group(1)), unit_match.group(2)
        return None, None

    # Apply the extraction function to split concentration and unit
    df["Concentration"], df["Unit"] = zip(
        *df[column_name].apply(extract_unit_and_numeric)
    )

    # Remove rows where units are not recognized
    df = df.dropna(subset=["Unit", "Concentration"])

    # Rearrange columns to place 'Concentration' and 'Unit' adjacently
    cols = [c for c in df.columns if c not in ["Concentration", "Unit"]]
    df = df[cols + ["Concentration", "Unit"]]

    return df


def summarize(
    df: pd.DataFrame,
    id_var: str = "Gene",
    prefix: str = "ModZScore ",
    value_name: str = "ModZScore",
) -> pd.DataFrame:
    """Summarize/reshape dataframe by certain statistics.

    e.g. the default parameters would take a dataframe of the form:

    | Gene   | ModZScore 0.1 | ModZScore 1.0 | ModZScore 10.0 | OtherColumn |
    |--------|---------------|---------------|----------------|-------------|
    | GATA1  | 0.5           | 1.2           | 2.7            | value1      |
    | SOX2   | -0.3          | 0.8           | 1.9            | value2      |
    | PAX6   | 0.2           | 0.4           | 1.5            | value3      |

    and reshape to

    | Gene   | Concentration | ModZScore |
    |--------|--------------|-----------|
    | GATA1  | 0.1          | 0.5       |
    | GATA1  | 1.0          | 1.2       |
    | GATA1  | 10.0         | 2.7       |
    | SOX2   | 0.1          | -0.3      |
    | SOX2   | 1.0          | 0.8       |
    | SOX2   | 10.0         | 1.9       |
    | PAX6   | 0.1          | 0.2       |
    | PAX6   | 1.0          | 0.4       |
    | PAX6   | 10.0         | 1.5       |

    Note the example above was generated using AI on 10/2025.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    id_var : str, optional
        Column name for summary statistic of interest, by default "Gene"
    prefix : str, optional
        Column prefix for summary statistics to group, by default "ModZScore "
    value_name : str, optional
        Name of grouped statistic, by default "ModZScore"

    Returns
    -------
    pd.DataFrame
        Reshaped dataframe.
    """
    if " " not in prefix:
        prefix = f"{prefix} "

    df = df.melt(
        id_vars=[id_var],
        value_vars=[c for c in df.columns if prefix in c],
        var_name="Concentration",
        value_name=value_name,
    )
    df["Concentration"] = [c.lstrip(prefix) for c in df["Concentration"]]
    return df


def getGoTerms(chemical_id: str, project: str) -> pd.DataFrame:
    """Retrieve Gene Ontology (GO) terms via the BU Xposome portal.

    Parameters
    ----------
    chemical_id : str
        Chemical identifier (e.g. CAS number)
    project : str
        Name of project by which to query API

    Returns
    -------
    pd.DataFrame
        Dataframe containing GO summarized info
    """
    url = f"https://montilab.bu.edu/Xposome-API/gs_enrichment?project={project}&chemical_id={chemical_id}"
    res = requests.get(url)
    if res.status_code != 200:
        return pd.DataFrame()

    # Load GO terms enrichment summary statistics data
    summary = pd.json_normalize(res.json())

    if len(summary) > 0:
        # Reshape into 3 cols: GenesetName, Concentration, and GSScore
        summary = summarize(
            summary, id_var="GenesetName", prefix="GS Score ", value_name="GSScore"
        )
        summary[["Condition", "Concentration"]] = summary["Concentration"].str.split(
            "_", expand=True
        )
        summary = format_concentration(summary)  # format conc + unit cols

    display_link = f"https://montilab.bu.edu/Xposome/?page={project}&tab=chemical_explorer&chemical_id={chemical_id}&stat=gs_enrichment"

    summary["Project"] = project
    summary["cas_number"] = chemical_id
    summary["Link"] = display_link
    return summary


def getGenes(chemical_id, project):
    """Retrieve and gene expression data from the BU Xposome portal.

    Parameters
    ----------
    chemical_id : str
        Chemical identifier (e.g. CAS number)
    project : str
        Name of project (e.g. "TG-GATEs", "MCF10A")

    Returns
    -------
    pd.DataFrame
        Processed dataframe containing summarized gene expression
        data. Returns an empty DataFrame if no data is found or
        if an error occurs.

    Notes
    -----
    The function applies project-specific processing:
    - For TG-GATEs: Parses "Condition_High_Concentration" format
    - For MCF10A: Parses "Condition_Concentration" format
    - For other projects: Assigns "WT" as default condition
    """
    url = f"https://montilab.bu.edu/Xposome-API/gene_expression?project={project}&chemical_id={chemical_id}&landmark=FALSE&do.scorecutoff=FALSE"

    try:
        res = requests.get(url)

        # Load gene expression summary statistics data
        summary = pd.DataFrame(json.loads(res.json()[0]))

        if summary.empty:
            return pd.DataFrame()

        # Reshape into 3 cols: Gene, Concentration, and ModZScore
        summary = summarize(summary)

        if project == "TG-GATEs":
            summary[["Condition", "High", "Concentration"]] = summary[
                "Concentration"
            ].str.split("_", expand=True)
            summary.drop(columns=["High"], inplace=True)
        elif project == "MCF10A":
            summary[["Condition", "Concentration"]] = summary[
                "Concentration"
            ].str.split("_", expand=True)
        else:
            summary["Condition"] = "WT"

        summary = format_concentration(summary)  # format conc + unit cols

        display_link = f"https://montilab.bu.edu/Xposome/?page={project}&tab=chemical_explorer&chemical_id={chemical_id}&stat=gene_expression"

        summary["Project"] = project
        summary["cas_number"] = chemical_id
        summary["Link"] = display_link
        return summary

    except (requests.RequestException, json.JSONDecodeError, KeyError) as e:
        print(f"Error processing {chemical_id} in {project}: {e}")
        return pd.DataFrame()


# ===============================================
#  RUN SCRIPT
# ===============================================
if __name__ == "__main__":
    check_args()

    # Load all chemicals and projects
    chems = pd.read_csv(sys.argv[1], encoding="utf-8-sig").dropna(
        subset=["Chemical_ID"]
    )
    projects = _load_projects()

    genes, gos = list(), list()
    for proj in projects:
        c = _load_chemicals(proj)

        # Check for empty data
        if c is None or c.empty:
            print(f"No chemicals found for project {proj}")
            continue

        overlap = set(chems["cas_number"]).intersection(set(c["CAS"]))
        print(f"Found {len(overlap)} CAS ids in common in project {proj}")

        gg = pd.concat([getGenes(chem, proj) for chem in overlap], ignore_index=True)
        # gt = pd.concat([getGoTerms(chem, proj) for chem in overlap], ignore_index=True)
        genes.append(gg)  # gos.append(gt)

    # Combine all gene data and include friendly project names
    genes = pd.concat(genes, ignore_index=True)
    genes["project_id"] = genes["Project"]
    genes["Project"] = [PROJ2NAME[p] for p in genes["project_id"]]
    # gos = pd.concat(gos, ignore_index=True)

    # Get significant genes
    genes = genes.loc[genes["ModZScore"].abs() > 1.63].copy()
    chems = chems[["cas_number", "Chemical_ID"]].drop_duplicates()

    print(genes.columns)

    genes = (
        genes.groupby(["Project", "cas_number", "Concentration", "Link", "Condition"])
        .agg(nGenes=("Gene", "nunique"))
        .reset_index()
        .merge(chems, on="cas_number")
    )

    # Enforce schema
    genes = genes.rename(columns={"Concentration": "concentration"})

    genes.to_csv(join(OUTPUT_DIR, "exposomeGeneStats.csv"), index=False)
