"""build_script.py: Build SRP database from raw files.

authors: @sgosline, @christinehc
"""

# =========================================================
# Imports
# =========================================================
import argparse
import itertools
import os
import subprocess
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from src.data import FigshareDataLoader
from src.manifest import DataManifest
from src.params import MANIFEST_FILEPATH
from src.schema import (
    ZEBRAFISH_DTYPE_TO_SUFFIX,
    get_cols_from_schema,
    map_zebrafish_data_to_schema,
)
from tqdm import tqdm

# =========================================================
# Setup/Parameters
# =========================================================
OUTPUT_DIR = "tmp"  # "./tmp"


manifest = DataManifest(MANIFEST_FILEPATH)
loader = FigshareDataLoader(
    Path(OUTPUT_DIR) / ".figshare_cache",  # api_token=FIGSHARE_API_TOKEN
)


# =========================================================
# Functions
# =========================================================
def fitCurveFiles(
    morpho_filename: Union[list, str, None] = None,
    lpr_filename: Union[list, str, None] = None,
    output_dir: str = OUTPUT_DIR,
    file_prefix: str = "zebrafish",
):
    """Create benchmark dose response curve fit files.

    Parameters
    ----------
    morpho_filename : Optional[str]
        Path to file containing morphology data
    lpr_filename : Optional[str]
        Path to file containing behavioral data
    output_dir : str, optional
        Path to which to save output files, by default OUTPUT_DIR
    file_prefix : str, optional
        Prefix for output filenames, by default "zebrafish"

    Raises
    ------
    ValueError
        If one of `morpho_filename` and `lpr_filename` is not provided.
    subprocess.CalledProcessError
        If zfBmd/main.py is not executed successfully
            (i.e. process exits with non-zero return code)
    Exception
        If an Exception is raised while executing zfBmd/main.py

    Examples
    --------
    >>> import pandas as pd
    >>> fitCurveFiles(
        "/path/to/morphology.csv", "/path/to/behavioral.csv", output_dir, file_prefix
        )
    # Creates 3 files: "{file_prefix}_chem_{[BMDs, Dose, Fits]}.csv"
    #   in the output directory

    """
    args = ""
    if morpho_filename is None and lpr_filename is None:
        raise ValueError(
            "At least one of `morpho_filename` and `lpr_filename`"
            " must be provided by the user."
        )

    # Construct flexible shell command from input
    for cli_arg, filename in zip(
        ["--morpho", "--lpr"], [morpho_filename, lpr_filename]
    ):
        if filename:
            if isinstance(filename, list):
                figshare_id_list = [f.split("/")[-1] for f in filename]
                for fid in figshare_id_list:
                    _ = loader.load_data(fid)
                figshare_files = " ".join(
                    [loader.get_file_path(fid).as_posix() for fid in figshare_id_list]
                )
                args = f"{args} {cli_arg} {figshare_files}"
            if isinstance(filename, str) or isinstance(filename, Path):
                figshare_id = filename.split("/")[-1]
                _ = loader.load_data(figshare_id)
                args = f"{args} {cli_arg} {loader.get_file_path(figshare_id)}"

    cmd = f"python -u zfBmd/main.py {args} --output_dir {output_dir} --prefix {file_prefix}"
    # tqdm.write(cmd)

    try:
        process = subprocess.run(cmd, capture_output=True, text=True, shell=True)

        # Verify successful command execution
        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=process.returncode,
                cmd=cmd,
                output=process.stdout,
                stderr=process.stderr,
            )

        # Show command line logging messages
        for line in process.stdout.splitlines():
            if line.strip():
                tqdm.write(line)
    except Exception as e:
        tqdm.write(f"An error occurred while trying to run the command: {str(e)}")
        raise e

    # When complete, clear data cache
    loader.clear_cache()


def combineZebrafishFiles(
    data_files: list[str], sample_type: str, data_type: str
) -> pd.DataFrame:
    """Combine preprocessed zebrafish sample files.

    Parameters
    ----------
    data_files : list[str]
        List of data files to concatenate
    sample_type : str
        Sample type, one of ["chemical", "extract"]
    data_type : str
        File type, one of ["bmd", "dose", "fit"]

    Returns
    -------
    pd.DataFrame
        Table of concatenated files with duplicates removed
    """
    required_cols = get_cols_from_schema(
        map_zebrafish_data_to_schema(sample_type, data_type),
    )

    tqdm.write(f"Concatenating {sample_type} {data_type}s...")
    if len(data_files) != 0:
        df = pd.concat([pd.read_csv(f) for f in data_files], ignore_index=True)
        df = df[required_cols].drop_duplicates()
        return df

    # If no files found, return empty df
    tqdm.write("Warning: No valid files found for concatenation")
    return pd.DataFrame(columns=required_cols[data_type])


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

    try:
        process = subprocess.run(cmd, capture_output=True, text=True, shell=True)

        # Verify successful command execution
        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=process.returncode,
                cmd=cmd,
                output=process.stdout,
                stderr=process.stderr,
            )

        # Show command line logging messages
        for line in process.stdout.splitlines():
            if line.strip():
                tqdm.write(line)
    except Exception as e:
        tqdm.write(f"An error occurred while trying to run the command: {str(e)}")
        raise e

    # Validate sample, chem, and mapping files
    dblist = [
        os.path.join(output_dir, "samples.csv"),
        os.path.join(output_dir, "chemicals.csv"),
        os.path.join(output_dir, "samplesToChemicals.csv"),
    ]
    return dblist


def runExposome(
    chem_id_file: str,
    output_dir: str = OUTPUT_DIR,
) -> list[str]:
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
    return [os.path.join(output_dir, "exposomeGeneStats.csv")]


def runExpression(
    gex: str,
    chem: str,
    ginfo: str,
    output_dir: str = OUTPUT_DIR,
) -> list[str]:
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
        os.path.join(output_dir, "srpDEGPathways.csv"),
        os.path.join(output_dir, "srpDEGStats.csv"),
        os.path.join(output_dir, "allGeneEx.csv"),
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


# =========================================================
# Command Line Parser
# =========================================================
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
        - zebrafish_BMDs_{BC,LPR}_{Chem,Samp}.csv
        - zebrafish_Dose_{BC,LPR}_{Chem,Samp}.csv
        - zebrafish_Fits_{BC,LPR}_{Chem,Samp}XYCoords.csv
    - exposomeGeneStats.csv (exposome analysis)
    - srpDEGPathways.csv, srpDEGStats.csv, allGeneEx.csv (gene expression results)

    Notes:
    ------
    - Intermediate files are created during processing and removed after use
    - All outputs are validated against the LinkML schema definitions
    - Progress is tracked using tqdm progress bars and informative messages
    """
    # ----------------------------
    # Command Line Argument Parser
    # ----------------------------
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
    parser.add_argument(
        "--output_dir",
        dest="output_dir",
        default=OUTPUT_DIR,
        help="Directory to store output files (default: '/tmp')",
    )
    args = parser.parse_args()

    # ---------------------------
    # File Parsing and Collection
    # ---------------------------
    # Map sample information
    sid = manifest.get(name="sampId")  # get_mapping_file(df, "sampId")
    cid = manifest.get(name="chemId")
    cclass = manifest.get(name="class1")
    emap = manifest.get(name="endpointMap")
    fses = manifest.get(data_type="sample", return_first=False)
    descfile = manifest.get(name="chemdesc")
    smap = manifest.get(name="sampMap")
    gex1 = manifest.get(data_type="expression", return_first=False)
    ginfo = manifest.get(name="geneInfo")

    # ------------------------------------------------------------------------
    # Benchmark Dose (BMD) Calculation / Sample-Chem Mapping (SAMPS) Workflows
    # ------------------------------------------------------------------------
    if args.bmd or args.samps:  ### need to rerun samples if we have created new bmds
        # Add chemical BMDS, fits, curves to existing data
        # sample_files, chem_files = [], []

        # Find morphology and behavior pairs for chemical extracts
        zebrafish_chem_files = manifest.get(
            data_type=["morphology", "behavior"],
            sample_type="chemical",
            version=4,
            return_first=False,
        )

        # Get zebrafish sample data
        zebrafish_samp_files = manifest.get(
            data_type=["bmd", "dose", "fit"],
            sample_type="extract",
            return_first=False,
            version=4,
        )

        # Define files and set progress bar increments for concatenating each
        total_iterations = len(zebrafish_chem_files) * len(zebrafish_samp_files)
        progress_bar = tqdm(total=total_iterations, desc="Combining files")

        # Process chemical files (using BMDRC) and collect output files
        tqdm.write(
            "Fitting benchmark dose response curves for zebrafish chemical extracts..."
        )
        # for chem_data in zebrafish_chem_files:
        #     morph_data, beh_data = chem_data
        fitCurveFiles(
            morpho_filename=[f[0] for f in zebrafish_chem_files],
            lpr_filename=[f[1] for f in zebrafish_chem_files],
            output_dir=args.output_dir,
            file_prefix="zebrafish",
        )
        fitted_chem_files = [
            os.path.join(args.output_dir, f"zebrafish_chem_{f}_{d}.csv")
            for f, d in itertools.product(["BMDs", "Dose", "Fits"], ["BC", "LPR"])
        ]

        # Process sample files (using preprocessed data)
        tqdm.write("Combining data for zebrafish sample extracts...")
        fitted_sample_files = list()
        for dtype, samp_data in zip(["bmd", "dose", "fit"], zip(*zebrafish_samp_files)):
            tqdm.write("Processing extracts data...")
            d = ZEBRAFISH_DTYPE_TO_SUFFIX[dtype]

            # TODO: fix this
            combined = combineZebrafishFiles(
                data_files=samp_data, sample_type="extract", data_type=dtype
            )
            combined_filename = os.path.join(
                args.output_dir, f"zebrafish_sample_{d}.csv"
            )
            combined.to_csv(combined_filename, index=False)
            fitted_sample_files.append(combined_filename)
            progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Combining files... Done!")
        progress_bar.close()

        # TODO: Add LinkML validation
        # for ftype in ["XYCoords.csv", "DoseResponse.csv", "BMDs.csv"]:
        #     dblist.append(os.path.join(output_dir, f"zebrafish_chem_{ftype}"))
        #     dblist.append(os.path.join(output_dir, f"zebrafish_samp_{ftype}"))

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
            "output_dir": args.output_dir,
        }

        # Iterate through sampMap params
        sampmap_params = [
            {"is_sample": True, "drcfiles": fitted_sample_files},
            {"is_sample": False, "drcfiles": fitted_chem_files},
            {"is_sample": False, "drcfiles": []},
        ]
        progress_bar = tqdm(
            range(len(sampmap_params)),
            desc="Running sample mapping",
        )

        # Perform sample mapping
        for smp in sampmap_params:
            smpargs = {**smp, **sampmap_args}
            res = runSampMap(**smpargs)
            all_res.extend(res)
            progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Running sample mapping... Done!")
        progress_bar.close()

        # Collect all unique files and remove temp files
        all_res = list(set(all_res))
        # for f in fitted_sample_files + fitted_chem_files:
        #     os.system(f"rm {f}")

        # Validate schema
        runSchemaCheck(res)

    # -----------------
    # Exposome Workflow
    # -----------------
    if args.expo:
        res = runExposome(cid, output_dir=args.output_dir)
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(res)

    # ------------------------
    # Gene Expression Workflow
    # ------------------------
    if args.geneEx:
        if not os.path.exists(os.path.join(args.output_dir, "chemicals.csv")):
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
                output_dir=args.output_dir,
            )

        res = runExpression(
            gex1,
            os.path.join(args.output_dir, "chemicals.csv"),
            ginfo,
            output_dir=args.output_dir,
        )
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(res)


if __name__ == "__main__":
    main()
