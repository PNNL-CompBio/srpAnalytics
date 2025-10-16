"""
Build script moved to python for better extendability and interoperability.

"""

import argparse
import os
from typing import Optional

import pandas as pd
from src.mapping import get_mapping_file, load_mapping_reference
from tqdm import tqdm

# DEFINE OUTPUT DIRECTORY
OUTPUT_DIR = "/tmp"  # "/tmp"


# TODO: Write this function
def fitCurveFiles(morpho_behavior_tuples):
    """
    get new curve fits, list of tuples of morpho/behavior pairs
    """


def combineFiles(location_list: pd.DataFrame, ftype: str) -> pd.DataFrame:
    """Combine files and account for duplicates and desired schema.

    Parameters
    ----------
    location_list : pd.DataFrame
        _description_
    ftype : str
        File type, one of ["bmd", "fit", "dose"]

    Returns
    -------
    pd.DataFrame
        Table of concatenated files with duplicates removed
    """
    dflist = []
    required_cols = {
        "bmd": [
            "Chemical_ID",
            "End_Point",
            "Model",
            "BMD10",
            "BMD50",
            "Min_Dose",
            "Max_Dose",
            "AUC_Norm",
            "DataQC_Flag",
            "BMD_Analysis_Flag",
        ],  # ,"BMD10_Flag","BMD50_Flag{"),
        "dose": ["Chemical_ID", "End_Point", "Dose", "Response", "CI_Lo", "CI_Hi"],
        "fit": ["Chemical_ID", "End_Point", "X_vals", "Y_vals"],
    }

    tqdm.write(f"Concatenating {ftype}...")
    for loc in location_list.location:
        f = pd.read_csv(loc)[required_cols[ftype]]
        dflist.append(f)

    # Check for empty list
    if not dflist:
        tqdm.write("Warning: No valid files found for concatenation")
        return pd.DataFrame(columns=required_cols[ftype])

    df = pd.concat(dflist).drop_duplicates()
    return df


def runSampMap(
    is_sample: bool = False,
    drcfiles: list = [],
    sid: str = "",
    smap: str = "",
    cid: str = "",
    emap: str = "",
    cclass: str = "",
    fses: str = "",
    descfile: str = "",
    output_dir: str = OUTPUT_DIR,
) -> list[str]:
    """Run sample-to-chemical mapping.

    Parameters
    ----------
    is_sample : bool, optional
        If True, runs sample mapping mode; else, runs
        chemical mapping mode, by default False
    drcfiles : list, optional
        List of dose-response curve files to process, by default []
    sid : str, optional
        Sample ID, by default ""
    smap : str, optional
        /path/to/sample_mapping_file, by default ""
    cid : str, optional
        Chemical ID, by default ""
    emap : str, optional
        /path/to/endpoint_mapping_file, by default ""
    cclass : str, optional
        /path/to/chemical_class_file, by default ""
    fses : str, optional
        /path/to/sample_files, by default ""
    descfile : str, optional
        /path/to/chemical_description_file, by default ""
    output_dir : str, optional
        Directory to save output, by default OUTPUT_DIR (='/tmp')

    Returns
    -------
    list[str]
        List of paths to the generated output files, including:
        - Core data files
            - samples.csv
            - chemicals.csv
            - samplesToChemicals.csv
        - Zebrafish files for both chemical and sample measurements
            - zebrafish{Samp,Chem}XYCoords.csv
            - zebrafish{Samp,Chem}DoseResponse.csv
            - zebrafish{Samp,Chem}BMDs.csv)
    """
    import subprocess

    drc = ",".join(drcfiles)
    args = (
        f"--sample_id={sid} "
        f"--sample_map={smap} "
        f"--chem_id={cid} "
        f"--ep_map={emap} "
        f"--chem_class={cclass} "
        f"--sample_files={fses} "
        f"--chem_desc={descfile} "
        f"--output_dir={output_dir} "
    )
    if is_sample:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py --sample --drc_files={drc} {args}"
    elif len(drcfiles) > 0:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py --chemical --drc_files={drc} {args}"
    else:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py {args}"

    process = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    for line in process.stdout.splitlines():
        if line.strip():
            tqdm.write(line)
    # tqdm.write("\nRunning sample mapping with the following parameters:\n")
    # tqdm.write(f"{cmd}\n")
    # os.system(cmd)
    # tqdm.write(f"ls -la {output_dir} \n")
    # os.system(f"ls -la {output_dir}")
    ##now we validate the files that came out.
    dblist = [
        os.path.join(output_dir, "samples.csv"),
        os.path.join(output_dir, "chemicals.csv"),
        os.path.join(output_dir, "samplesToChemicals.csv"),
    ]
    for ftype in ["XYCoords.csv", "DoseResponse.csv", "BMDs.csv"]:
        dblist.append(os.path.join(output_dir, f"zebrafishChem{ftype}"))
        dblist.append(os.path.join(output_dir, f"zebrafishSamp{ftype}"))
    return dblist
    # runSchemaCheck(dblist)


def runExposome(chem_id_file: str) -> list[str]:
    """Pull exposome data.

    Parameters
    ----------
    chem_id_file : str
        Path to file containing chemical IDs for which to pull
        exposome data

    Returns
    -------
    list[str]
        List containing path to output exposomeGeneStats.csv file
    """
    cmd = f"python exposome/exposome_summary_stats.py {chem_id_file}"
    tqdm.write(cmd)
    os.system(cmd)
    return [os.path.join(OUTPUT_DIR, "exposomeGeneStats.csv")]


def runExpression(gex: str, chem: str, ginfo: str) -> list[str]:
    """Parse gene expression data using R.

    Parameters
    ----------
    gex : str
        Path to gene expression data
    chem : str
        Path to chemical data file
    ginfo : str
        Path to gene info file

    Returns
    -------
    list[str]
        List of these three output files:
            - "{OUTPUT_DIR}/srpDEGPathways.csv": Enriched pathways
                in differentially expressed genes
            - "{OUTPUT_DIR}/srpDEGStats.csv" : Summary statistics
                for differentially expressed genes
            - "{OUTPUT_DIR}/allGeneEx.csv" : All gene expression data
        Note that OUTPUT_DIR = "/tmp" by default.
    """
    cmd = f"Rscript zfExp/parseGexData.R {gex} {chem} {ginfo}"
    tqdm.write(cmd)
    os.system(cmd)
    return [
        os.path.join(OUTPUT_DIR, "srpDEGPathways.csv"),
        os.path.join(OUTPUT_DIR, "srpDEGStats.csv"),
        os.path.join(OUTPUT_DIR, "allGeneEx.csv"),
    ]


def runSchemaCheck(dbfiles: list[Optional[str]] = []):
    """Validate database files against schema using LinkML.

    Parameters
    ----------
    dbfiles : list[Optional[str]], optional
        List of database files, by default []
    """
    ##TODO: make this work with internal calls
    for filename in dbfiles:
        classname = os.path.basename(filename).split(".")[0]
        cmd = f"linkml-validate --schema srpAnalytics.yaml {filename} --target-class {classname}"
        tqdm.write(cmd)
        os.system(cmd)


def main():
    """Run data processing and analytics pipeline for Superfund data.

    This is the main entrypoint for the Superfund data processing
    pipeline designed to run inside a Docker container. The pipeline
    processes various chemical and sample data files, performs
    benchmark dose calculations, maps samples to chemicals, runs
    exposome analyses, and processes gene expression data. The
    workflow is modular, allowing specific components to be executed
    at a time depending on the supplied command line arguments.

    Workflow Components:
    -------------------
    1. Data Preparation:
       - Loads mapping reference data
       - Identifies morphology and behavior data pairs for chemicals
       - Retrieves various mapping files (sample IDs, chemical IDs, endpoints, etc.)

    2. Benchmark Dose (BMD) Analysis:
       - Calculates dose-response curves and benchmark doses for chemical exposures
       - Combines results across different sample types (chemical, extract) and data types (BMD, fit, dose)

    3. Sample-Chemical Mapping:
       - Links samples to chemicals using various reference files
       - Validates outputs against schema definitions

    4. Exposome Analysis:
       - Processes exposome data for chemicals to identify environmental exposures

    5. Gene Expression Analysis:
       - Processes differential gene expression data associated with chemical exposures
       - Performs pathway analysis on differentially expressed genes

    Command-line Arguments:
    ----------------------
    `--bmd` : Re-run benchmark dose calculation and dependent commands
    `--samps` : Re-run sample-chemical mapping
    `--expo` : Re-run exposome sample collection
    `--geneEx` : Re-run gene expression generation

    Outputs:
    --------
    Various CSV files stored in OUTPUT_DIR, including:
    - Core data
        - samples.csv
        - chemicals.csv
        - samplesToChemicals.csv
    - Zebrafish assay data
        - zebrafish{Chem,Samp}BMDs.csv
        - zebrafish{Chem,Samp}DoseResponse.csv
        - zebrafish{Chem,Samp}XYCoords.csv
    - exposomeGeneStats.csv (exposome analysis)
    - srpDEGPathways.csv, srpDEGStats.csv, allGeneEx.csv (gene expression results)

    Notes:
    ------
    - Intermediate files are created during processing and removed after use
    - All outputs are validated against the LinkML schema definitions
    - Progress is tracked using tqdm progress bars and informative messages
    """
    df = load_mapping_reference()

    ####
    # file parsing - collects all files we might need for the tool below
    ####
    ##first find the morphology and behavior pairs for chemical sources
    chemdf = df.loc[df.sample_type == "chemical"]
    morph = chemdf.loc[chemdf.data_type == "morphology"]
    beh = chemdf.loc[chemdf.data_type == "behavior"]
    tupes = []
    for n in morph.name:
        tupes.append(
            [morph.loc[morph.name == n].location, beh.loc[beh.name == n].location]
        )

    ##now map sample information
    sid = get_mapping_file(df, "sampId")
    cid = get_mapping_file(df, "chemId")
    cclass = get_mapping_file(df, "class1")
    emap = get_mapping_file(df, "endpointMap")
    fses = get_mapping_file(df, "sample", return_first=False)
    descfile = get_mapping_file(df, "chemdesc")
    smap = get_mapping_file(df, "sampMap")
    gex1 = get_mapping_file(df, "expression", return_first=False)
    ginfo = get_mapping_file(df, "geneInfo")

    ###now we can call individiual commands
    parser = argparse.ArgumentParser(
        "Pull files from github list of files and call appropriate command"
    )
    parser.add_argument(
        "--bmd",
        dest="bmd",
        action="store_true",
        default=False,
        help="Re-run benchmark dose calculation and dependent commands",
    )
    parser.add_argument(
        "--samps",
        dest="samps",
        action="store_true",
        default=False,
        help="Re run sample-chem mapping",
    )
    parser.add_argument(
        "--expo",
        dest="expo",
        action="store_true",
        default=False,
        help="Re run exposome sample collection",
    )
    parser.add_argument(
        "--geneEx",
        dest="geneEx",
        action="store_true",
        default=False,
        help="Re run gene expression generation",
    )

    args = parser.parse_args()

    ##call bmdrc on all morphology/behavior pairs for sample sources
    if args.bmd:
        tqdm.write("Re-running benchmark dose collection...")
        newbmds, newfits, newdoses = [], [], []
        fitCurveFiles()

    # ------------------------------------------------------------------------
    # Benchmark Dose (BMD) Calculation / Sample-Chem Mapping (SAMPS) Workflows
    # ------------------------------------------------------------------------
    if args.bmd or args.samps:  ### need to rerun samples if we have created new bmds
        # add chemical BMDS, fits, curves to existing data
        chem_files, sample_files = [], []

        # Define files and set progress bar incrementes for concatenating each
        sample_type = ["chemical", "extract"]
        data_type = ["bmd", "fit", "dose"]
        total_iterations = len(sample_type) * len(data_type)
        progress_bar = tqdm(total=total_iterations, desc="Combining files")

        for st in sample_type:
            tqdm.write(f"Processing {st} samples...")

            for dt in data_type:
                fdf = combineFiles(
                    df.loc[df.sample_type == st].loc[df.data_type == dt], dt
                )
                fname = os.path.join(OUTPUT_DIR, f"tmp_{st}_{dt}.csv")
                fdf.to_csv(fname, index=False)
                if st == "chemical":
                    chem_files.append(fname)
                else:
                    sample_files.append(fname)
                progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Combining files... Done!")
        progress_bar.close()

        # Define fixed params for sample mapping
        all_res = list()
        sampmap_args = {
            "sid": sid,
            "smap": smap,
            "cid": cid,
            "emap": emap,
            "cclass": cclass,
            "fses": fses,
            "descfile": descfile,
        }

        # Iterate through sampMap params
        sampmap_params = [
            {"is_sample": True, "drcfiles": sample_files},
            {"is_sample": False, "drcfiles": chem_files},
            {"is_sample": False, "drcfiles": []},
        ]
        progress_bar = tqdm(
            range(len(sampmap_params)),
            desc="Running sample mapping",
        )

        # Perform sample mapping
        for smp in sampmap_params:
            smpargs = {**smp, **sampmap_args}
            tqdm.write(
                f"Performing sample mapping with the following parameters: {smpargs}"
            )
            res = runSampMap(**smpargs)
            all_res.extend(res)
            progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Running sample mapping... Done!")
        progress_bar.close()

        # Collect all unique files and remove temp files
        all_res = list(set(all_res))
        for f in sample_files + chem_files:
            # tqdm.write(f"Filename: {f}")
            # os.system(f"head {f}")
            os.system(f"rm {f}")

        # Validate schema
        runSchemaCheck(res)

    # -----------------
    # Exposome Workflow
    # -----------------
    if args.expo:
        res = runExposome(cid)
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(res)

    # ------------------------
    # Gene Expression Workflow
    # ------------------------
    if args.geneEx:
        if not os.path.exists(os.path.join(OUTPUT_DIR, "chemicals.csv")):
            runSampMap(
                is_sample=False,
                drcfiles=[],
                sid=sid,
                smap=smap,
                cid=cid,
                emap=emap,
                cclass=cclass,
                fses=fses,
                descfile=descfile,
                output_dir=OUTPUT_DIR,
            )

        res = runExpression(gex1, os.path.join(OUTPUT_DIR, "chemicals.csv"), ginfo)
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(res)


if __name__ == "__main__":
    main()
