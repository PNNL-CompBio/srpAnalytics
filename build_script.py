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
import sys
import traceback
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from src.data import FigshareDataLoader, figshare_url_to_id, load_figshare_url
from src.manifest import DataManifest
from src.params import MANIFEST_FILEPATH
from src.samples import combine_chemical_data, combine_chemical_endpoint_data
from src.schema import get_cols_from_schema, map_zebrafish_data_to_schema
from tqdm import tqdm

# =========================================================
# Setup/Parameters
# =========================================================
OUTPUT_DIR = "tmp"  # "./tmp"

manifest = DataManifest(MANIFEST_FILEPATH)
loader = FigshareDataLoader(
    Path(OUTPUT_DIR) / ".figshare_cache", api_token=os.getenv("FIGSHARE_API_TOKEN")
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
    morpho_filename : Union[list, str, None]
        Path or list of paths to file(s) containing morphology data
    lpr_filename : Optional[str]
        Path or list of paths to file(s) containing behavioral data
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
                figshare_id_list = [figshare_url_to_id(f) for f in filename]
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
        process = subprocess.run(cmd, text=True, shell=True)  # capture_output=True,

        # Verify successful command execution
        if process.returncode != 0:
            sys.stderr.write("=== SUBPROCESS STDOUT ===\n")
            sys.stderr.write(process.stdout if process.stdout else "(empty)\n")
            sys.stderr.write("=== SUBPROCESS STDERR ===\n")
            sys.stderr.write(process.stderr if process.stderr else "(empty)\n")
            sys.stderr.flush()
            raise subprocess.CalledProcessError(
                returncode=process.returncode,
                cmd=cmd,
                output=process.stdout,
                stderr=process.stderr,
            )
        for line in process.stdout.splitlines():
            if line.strip():
                tqdm.write(line)

        # Show command line logging messages
        for line in process.stdout.splitlines():
            if line.strip():
                tqdm.write(line)
    except Exception as e:
        tqdm.write(f"An error occurred while trying to run the command: {str(e)}")
        raise e

    # When complete, clear data cache
    loader.clear_cache()


# Combine all zebrafish files. Includes both chem and sample data
def combineZebrafishFiles(
    data_files: list[str],
    sample_type: str,
    data_type: str,
    ids: pd.DataFrame,
    sample_id_map: Optional[pd.DataFrame] = None,
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
    ids : pd.DataFrame
        Full list of IDs for data type
            - Loaded from chemicals.csv for chemicals data
            - Loaded from samples.csv for sample/extracts data
    sample_id_map : Optional[pd.DataFrame]
        DataFrame containing Sample_ID and SampleNumber column
            mappings for all samples, by default None
        NOTE: This parameter will remain optional/unused until
            the pipeline actually processes sample data
            (currently it does not)

    Returns
    -------
    pd.DataFrame
        Table of concatenated files with duplicates removed
    """
    required_cols = get_cols_from_schema(
        map_zebrafish_data_to_schema(sample_type, data_type),
    )

    # Remove End_Point_Name (currently added later)
    required_cols = [
        c
        for c in required_cols
        if c not in ["Chemical_ID", "Sample_ID", "End_Point_Name"]
    ]
    id_col = ["Chemical_ID"] if sample_type == "extract" else None
    required_cols = id_col + required_cols

    tqdm.write(f"Concatenating {sample_type} {data_type} files...")
    if len(data_files) != 0:
        df = pd.concat([pd.read_csv(f) for f in data_files], ignore_index=True)
        df = df[required_cols].drop_duplicates()

        # Process extracts
        if sample_type == "extract":
            # Use temp (split) IDs for joining
            tmp_ids = ids.copy()
            split_ids = tmp_ids["Sample_ID"].str.split("-", expand=True)
            tmp_ids["tmp_id"] = split_ids[0]
            tmp_ids = tmp_ids[["Sample_ID", "tmp_id"]].drop_duplicates()
            df["tmp_id"] = df["Chemical_ID"].astype(str)
            df = df.drop(columns=["Chemical_ID"])
            df = df.merge(tmp_ids, on="tmp_id", how="left")

            # Fill missing sample IDs with temp ID
            mask = df["Sample_ID"].isna()
            df.loc[mask, "Sample_ID"] = df.loc[mask, "tmp_id"]

            # Delete temp ID column and move sample ID to front
            df = df.drop(columns=["tmp_id"])
            cols = ["Sample_ID"] + [col for col in df.columns if col != "Sample_ID"]
            df = df[cols]

        # For non-extracts, keep only chemicals that are in sample data
        elif sample_type != "extract" and ids is not None:
            df = df[df["Chemical_ID"].isin(ids["Chemical_ID"])]

        return df.drop_duplicates()

    # If no files found, return empty df
    tqdm.write("Warning: No valid files found for concatenation")
    return pd.DataFrame(columns=required_cols[data_type])


def combineZebrafishSampleFiles(
    bmd_files: list[str],
    dose_files: list[str],
    fit_files: list[str],
    chem_data: pd.DataFrame,
    endpoint_metadata: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
) -> list[str]:
    """Combine preprocessed zebrafish sample files into final output files.

    Parameters
    ----------
    bmd_files : list[str]
        List of BMD files
    dose_files : list[str]
        List of dose response files
    fit_files : list[str]
        List of fit files
    chem_data : pd.DataFrame
        Sample data containing Sample_ID mappings
    endpoint_metadata : pd.DataFrame
        Endpoint metadata for merging
    output_dir : str, optional
        Output directory, by default OUTPUT_DIR

    Returns
    -------
    list[str]
        List of paths to generated output files
    """
    output_files = []

    # BMDs
    tqdm.write("Combining BMD data for zebrafish sample extracts...")
    bmds = (
        combine_chemical_endpoint_data(
            bmd_files,
            is_extract=True,
            chem_data=chem_data,
            endpoint_metadata=endpoint_metadata,
        )
        .dropna(subset=["BMD_Analysis_Flag"])
        .query("BMD_Analysis_Flag != 'NA'")
    )
    bmd_output = os.path.join(output_dir, "zebrafishSampBMDs.csv")
    bmds.fillna("NULL").to_csv(bmd_output, index=False, quotechar='"')
    output_files.append(bmd_output)

    # XYCoords/Fits
    tqdm.write("Combining fits data for zebrafish sample extracts...")
    curves = combine_chemical_data(
        fit_files,
        data_type="fit",
        is_extract=True,
        chem_data=chem_data,
        endpoint_metadata=endpoint_metadata,
    )
    fits_output = os.path.join(output_dir, "zebrafishSampXYCoords.csv")
    curves.fillna("NULL").to_csv(fits_output, index=False, quotechar='"')
    output_files.append(fits_output)

    # Dose Response
    tqdm.write("Combining dose response data for zebrafish sample extracts...")
    dose_reps = combine_chemical_data(
        dose_files,
        data_type="dose",
        is_extract=True,
        chem_data=chem_data,
        endpoint_metadata=endpoint_metadata,
    ).dropna(subset=["Dose"])
    dose_output = os.path.join(output_dir, "zebrafishSampDoseResponse.csv")
    dose_reps.fillna("NULL").to_csv(dose_output, index=False, quotechar='"')
    output_files.append(dose_output)
    return output_files


def runSampMap(
    sample_id_file: str = "",
    sample_map_file: str = "",
    chemical_id: str = "",
    endpoint_map: str = "",
    chem_class_file: str = "",
    fses_files: str = "",
    chem_desc_file: str = "",
    output_dir: str = OUTPUT_DIR,
) -> list[str]:
    """Run sample-to-chemical mapping.

    Parameters
    ----------
    sample_id_file : str, optional
        File location for Sample ID mapping, by default ""
    sample_map_file : str, optional
        /path/to/sample_mapping_file, by default ""
    chemical_id : str, optional
        Chemical ID, by default ""
    endpoint_map : str, optional
        /path/to/endpoint_mapping_file, by default ""
    chem_class_file : str, optional
        /path/to/chemical_class_file, by default ""
    fses_files : str, optional
        /path/to/sample_files, by default ""
    chem_desc_file : str, optional
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
    args = (
        f"--sample_id_file={sample_id_file} "
        f"--sample_map={sample_map_file} "
        f"--chemical_id={chemical_id} "
        f"--endpoint_map={endpoint_map} "
        f"--chemical_class={chem_class_file} "
        f"--sample_files={fses_files} "
        f"--chemical_description={chem_desc_file} "
        f"--output_dir={output_dir} "
    )
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
        traceback.print_exception(e)
        raise e

    # TODO: Validate sample, chem, and mapping files

    # Return output files
    output_files = (
        os.path.join(output_dir, "samples.csv"),
        os.path.join(output_dir, "chemicals.csv"),
        os.path.join(output_dir, "samplesToChemicals.csv"),
    )
    return output_files


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


def runSchemaCheck(
    files: list[Optional[str]] = None, classes: Union[str, list[str], None] = None
):
    """Validate database files against schema using LinkML.

    Parameters
    ----------
    files : list[Optional[str]], optional
        List of database files, by default []
    classes : Union[str, list[str], None]
        Class name or list of class names, optional, by default None
        If one class name is supplied, it is assumed to apply to all
            files in the corresponding list.
        If a list is supplied, it is assumed to correspond with the
            files in the file list
            (i.e. [class_a, class_b] maps to [file_a, file_b])
        If no class name(s) is/are supplied, class is automatically
            determined via the filename.
    """
    if files is None:
        files = []

    # If classname = str, apply same classname to all
    if isinstance(classes, str) or classes is None:
        classes = [classes] * len(files)

    # Ensure class list matches file list (if list supplied)
    if isinstance(classes, list) and len(classes) != len(files):
        raise ValueError(
            "Classnames must correspond to filenames "
            f"({len(files)} files supplied with {len(classes)}"
            " classnames)."
        )

    ##TODO: make this work with internal calls
    for filename, classname in zip(files, classes):
        if classname is None:
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
    tqdm.write("Retrieving files from manifest...")
    sample_id_file = manifest.get(name="sampId")  # get_mapping_file(df, "sampId")
    chemical_id = manifest.get(name="chemId", version=4)
    chem_class_file = manifest.get(name="class1")
    endpoint_map = manifest.get(name="endpointMap", version=4)
    fses_files = manifest.get(
        data_type="sample", return_first=False, version=4
    )  # use new files
    chem_desc_file = manifest.get(name="chemdesc")
    sample_map_file = manifest.get(name="sampMap")
    gex1 = manifest.get(data_type="expression", return_first=False)
    ginfo = manifest.get(name="geneInfo")

    # Run sample-to-chemical mapping
    tqdm.write("Running sample/chemical mapping...")
    sampmap_args = {
        "sample_id_file": sample_id_file,
        "sample_map_file": sample_map_file,
        "chemical_id": chemical_id,
        "endpoint_map": endpoint_map,
        "chem_class_file": chem_class_file,
        "fses_files": fses_files,
        "chem_desc_file": chem_desc_file,
        "output_dir": args.output_dir,
    }
    samples_file, chemicals_file, samples_to_chemicals_file = runSampMap(**sampmap_args)

    # ------------------------------------------------------------------------
    # Benchmark Dose (BMD) Calculation / Sample-Chem Mapping (SAMPS) Workflows
    # ------------------------------------------------------------------------
    if args.bmd or args.samps:  ### need to rerun samples if we have created new bmds
        # Add chemical BMDS, fits, curves to existing data
        # sample_files, chem_files = [], []

        # Find morphology data for chemical extracts
        zebrafish_chem_morpho = manifest.get(
            data_type="morphology",  # ["morphology", "behavior"]
            sample_type="chemical",
            version=4,
            return_first=False,
        )

        # Get zebrafish chemical LPR data (pre-processed)
        zebrafish_chem_lpr = manifest.get(
            data_type=["bmd", "dose", "fit"],
            sample_type="chemical",
            return_first=False,
            version=4,
        )

        # Get zebrafish sample data
        zebrafish_samp_files = manifest.get(
            data_type=["bmd", "dose", "fit"],
            sample_type="extract",
            return_first=False,
            version=4,
        )

        # Define files and set progress bar increments for concatenating each
        total_iterations = 3
        progress_bar = tqdm(total=total_iterations, desc="Combining files")

        # Process chemical files (using BMDRC) and collect output files
        tqdm.write(
            "Fitting benchmark dose response curves for zebrafish chemical extracts..."
        )
        fitted_chem_files = [
            os.path.join(args.output_dir, f"zebrafish_chem_{f}_{d}.csv")
            for f, d in itertools.product(["BMDs", "Dose", "Fits"], ["BC", "LPR"])
        ]
        if not os.path.exists(fitted_chem_files[0]):  # Skip if files exist
            fitCurveFiles(
                morpho_filename=zebrafish_chem_morpho,  # [f[0] for f in zebrafish_chem_files],
                lpr_filename=None,  # [f[1] for f in zebrafish_chem_files],
                output_dir=args.output_dir,
                file_prefix="zebrafish",
            )

        # Process LPR and add to chem files
        for file_set in zip(*zebrafish_chem_lpr):
            tmp = list()
            for f in file_set:
                fid = f.split("/")[-1]
                _ = loader.load_data(fid)
                fname = loader.get_file_path(fid).as_posix()
                ftype = os.path.splitext(os.path.basename(fname))[0].split("_")[-1]
                tmp.append(pd.read_csv(fname))
            tmp = pd.concat(tmp, ignore_index=True)
            tmp.to_csv(
                os.path.join(args.output_dir, f"zebrafish_chem_{ftype}_LPR.csv"),
                index=False,
            )

        # Load endpoints and strip odd trailing spaces
        endpoint_names = load_figshare_url(
            loader, endpoint_map, sheet_name="Dictionary"
        )
        for c in endpoint_names.columns:
            endpoint_names[c] = [
                v.strip() if isinstance(v, str) else v for v in endpoint_names[c]
            ]

        # Combine LPR and BC data, add endpoint names, save
        temp_files = list()
        for ftype in ["BMDs", "Dose", "Fits"]:
            morpho_file = os.path.join(
                args.output_dir, f"zebrafish_chem_{ftype}_BC.csv"
            )
            behavior_file = os.path.join(
                args.output_dir, f"zebrafish_chem_{ftype}_LPR.csv"
            )
            tmp = pd.concat(
                [
                    pd.read_csv(morpho_file),
                    pd.read_csv(behavior_file),
                ],
                ignore_index=True,
            )
            tmp = (
                pd.merge(
                    tmp,
                    endpoint_names[["Abbreviation", "Simple name (<20char)"]],
                    how="left",
                    left_on="End_Point",
                    right_on="Abbreviation",
                )
                .rename(columns={"Simple name (<20char)": "End_Point_Name"})
                .drop(columns=["Abbreviation"])
            )
            tmp.to_csv(
                os.path.join(args.output_dir, f"zebrafishChem{ftype}.csv"), index=False
            )

            # Track files for later deletion
            temp_files.append(morpho_file)
            temp_files.append(behavior_file)

        # Process sample files (using preprocessed data)
        tqdm.write("Combining data for zebrafish sample extracts...")
        fitted_sample_files = list()
        for dtype, sample_data in zip(
            ["BMDs", "Dose", "Fits"], zip(*zebrafish_samp_files)
        ):
            tqdm.write("Processing extracts data...")
            samples = pd.read_csv(samples_file)

            # Combine zebrafish files
            combined = combineZebrafishFiles(
                data_files=sample_data,
                sample_type="extract",
                data_type=dtype,
                ids=samples,
            )
            combined = (
                pd.merge(
                    combined,
                    endpoint_names[["Abbreviation", "Simple name (<20char)"]],
                    how="left",
                    left_on="End_Point",
                    right_on="Abbreviation",
                )
                .rename(columns={"Simple name (<20char)": "End_Point_Name"})
                .drop(columns=["Abbreviation"])
            )
            combined_filename = os.path.join(
                args.output_dir, f"zebrafishSamp{dtype}.csv"
            )
            combined.to_csv(combined_filename, index=False)
            fitted_sample_files.append(combined_filename)
            progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Combining files... Done!")
        progress_bar.close()

        # Define fixed params for sample mapping
        all_results = list()

        # Collect all unique files and remove temp files
        all_results = list(set(all_results))
        for f in temp_files:
            os.remove(f)

        # Clean up separate LPR/BC files
        # os.remove(morpho_file)
        # os.remove(behavior_file)

        # Validate schema
        # TODO: fix schema check for combined files
        runSchemaCheck(all_results)
        for ftype in ["BMDs", "Dose", "Fits"]:
            runSchemaCheck(
                [os.path.join(args.output_dir, f"zebrafishChem{ftype}.csv")],
                classes=[
                    map_zebrafish_data_to_schema(
                        sample_type="chemical", data_type=ftype
                    )
                ],
            )
            runSchemaCheck(
                [os.path.join(args.output_dir, f"zebrafishSamp{ftype}.csv")],
                classes=[
                    map_zebrafish_data_to_schema(sample_type="extract", data_type=ftype)
                ],
            )

    # -----------------
    # Exposome Workflow
    # -----------------
    if args.expo:
        figshare_id = figshare_url_to_id(chemical_id)
        _ = loader.load_data(figshare_id)
        chem_id_map_file = loader.get_file_path(figshare_id).as_posix()
        result = runExposome(chem_id_map_file, output_dir=args.output_dir)
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(result)

    # ------------------------
    # Gene Expression Workflow
    # ------------------------
    if args.geneEx:
        # if not os.path.exists(os.path.join(args.output_dir, "chemicals.csv")):
        #     runSampMap(
        #         is_sample=False,
        #         dose_response_files=[],
        #         sample_id_file=sample_id_file,
        #         sample_map_file=sample_map_file,
        #         chemical_id=chemical_id,
        #         endpoint_map=endpoint_map,
        #         chem_class_file=chem_class_file,
        #         fses_files=fses_files,
        #         chem_desc_file=chem_desc_file,
        #         output_dir=args.output_dir,
        #     )

        result = runExpression(
            gex1,
            os.path.join(args.output_dir, "chemicals.csv"),
            ginfo,
            output_dir=args.output_dir,
        )
        # for f in res:
        #     tqdm.write(f"Filename: {f}")
        #     os.system(f"head {f}")
        runSchemaCheck(result)


if __name__ == "__main__":
    main()
