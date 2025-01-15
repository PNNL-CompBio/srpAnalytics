"""exposome_summary_stats: Converted from R to Python with AI Incubator."""

# ===============================================
#  IMPORTS
# ===============================================
import json
import re
import requests
import sys
from os.path import abspath, dirname, join
import pandas as pd

sys.path.append(dirname(dirname(abspath(__file__))))
from src.format import snakeify_all_columns

# ===============================================
#  CONFIG
# ===============================================
OUTPUT_DIR = "."

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
    """Load chemicals associated with project

    Parameters
    ----------
    project : str
        Project (e.g. ADIPO, HEPG2, etc.)

    Returns
    -------
    pd.DataFrame
        Tabulated chemical information with the following columns:
            Chemical_Id: str
            Project: str
            Chemical_Name: str
            CAS: str
    """
    url = f"https://montilab.bu.edu/Xposome-API/chemicals?projects={project}&chemical_ids=all"
    res = requests.get(url)
    if res.status_code == 200:
        return pd.DataFrame(json.loads(res.json()[0]))
    return


def format_concentration(
    df: pd.DataFrame, column_name: str = "Concentration"
) -> pd.DataFrame:
    """
    Splits a concentration column into separate numeric and unit columns,
    only recognizing units 'uM' or 'mg/kg' from a larger string.

    Parameters:
    - df: pandas DataFrame containing the concentration column.
    - column_name: str, name of the column to be split (default is 'Concentration').

    Returns:
    - DataFrame with separate 'Concentration' and 'Unit' columns.

    Written by AI incubator 01/2025
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


def getGoTerms(chemical_id, project):
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

    summary["project"] = project
    summary["cas_number"] = chemical_id
    summary["link"] = display_link
    return snakeify_all_columns(summary)


def getGenes(chemical_id, project):
    url = f"https://montilab.bu.edu/Xposome-API/gene_expression?project={project}&chemical_id={chemical_id}&landmark=FALSE&do.scorecutoff=FALSE"
    res = requests.get(url)
    if res.status_code != 200:
        return pd.DataFrame()

    # Load gene expression summary statistics data
    summary = pd.DataFrame(json.loads(res.json()[0]))

    if len(summary) > 0:
        # Reshape into 3 cols: Gene, Concentration, and ModZScore
        summary = summarize(summary)

        if proj == "TG-GATEs":
            summary[["Condition", "High", "Concentration"]] = summary[
                "Concentration"
            ].str.split("_", expand=True)
            summary.drop(columns=["High"], inplace=True)
        elif proj == "MCF10A":
            summary[["Condition", "Concentration"]] = summary[
                "Concentration"
            ].str.split("_", expand=True)
        else:
            summary["Condition"] = "WT"

        summary = format_concentration(summary)  # format conc + unit cols

    display_link = f"https://montilab.bu.edu/Xposome/?page={project}&tab=chemical_explorer&chemical_id={chemical_id}&stat=gene_expression"

    summary["project"] = project
    summary["cas_number"] = chemical_id
    summary["link"] = display_link
    return snakeify_all_columns(summary)


# ===============================================
#  RUN SCRIPT
# ===============================================
if __name__ == "__main__":
    check_args()

    # Load all chemicals and projects
    chems = snakeify_all_columns(pd.read_csv(sys.argv[1], encoding="utf-8-sig"))
    projects = _load_projects()

    genes, gos = list(), list()
    for proj in projects:
        c = _load_chemicals(proj)

        overlap = set(chems["cas_number"]).intersection(set(c["CAS"]))
        print(f"Found {len(overlap)} CAS ids in common in project {proj}")

        gg = pd.concat([getGenes(chem, proj) for chem in overlap], ignore_index=True)
        # gt = pd.concat([getGoTerms(chem, proj) for chem in overlap], ignore_index=True)
        genes.append(gg)  # gos.append(gt)

    # Combine all gene data and include friendly project names
    genes = pd.concat(genes, ignore_index=True)
    genes["project_id"] = genes["project"]
    genes["project"] = [PROJ2NAME[p] for p in genes["project_id"]]
    # gos = pd.concat(gos, ignore_index=True)

    # Get significant genes
    genes = genes.loc[genes["mod_z_score"].abs() > 1.63].copy()
    chems = chems[["cas_number", "chemical_id"]].drop_duplicates()

    genes = (
        genes.groupby(["project", "cas_number", "concentration", "link", "condition"])
        .agg(n_genes=("gene", "nunique"))
        .reset_index()
        .merge(chems, on="cas_number")
    )

    genes.to_csv(join(OUTPUT_DIR, "exposomeGeneStats.csv"), index=False)
