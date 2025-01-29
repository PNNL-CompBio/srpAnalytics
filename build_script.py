"""
Build script moved to python for better extendability and interoperability.

"""

import argparse
import os
from typing import Optional
import pandas as pd
from tqdm import tqdm
from src.mapping import load_mapping_reference, get_mapping_file

# DEFINE OUTPUT DIRECTORY
OUTPUT_DIR = "."  # "/tmp"


# TODO: Write this functions
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
    """
    run sample mapping
    """
    print("drcfiles:", drcfiles)
    drc = ",".join(drcfiles)
    args = (
        f" --sample_id={sid}"
        f" --sample_map={smap}"
        f" --chem_id={cid}"
        f" --ep_map={emap}"
        f" --chem_class={cclass}"
        f" --sample_files={fses}"
        f" --chem_desc={descfile}"
        f" --sample_map={smap}"
        f" --output_dir={output_dir}"
    )
    if is_sample:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py --sample --drc_files={drc} {args}"
    elif len(drcfiles) > 0:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py --chemical --drc_files={drc} {args}"
    else:
        cmd = f"python sampleChemMapping/map_samples_to_chemicals.py {args}"

    tqdm.write("\nRunning sample mapping with the following parameters:\n")
    tqdm.write(f"{cmd}\n")
    os.system(cmd)
    tqdm.write("ls -la . \n")
    os.system(f"ls -la {output_dir}")
    ##now we validate the files that came out.
    dblist = [
        os.path.join(output_dir, "samples.csv"),
        os.path.join(output_dir, "chemicals.csv"),
        os.path.join(output_dir, "sampleToChemicals.csv"),
    ]
    for ftype in ["XYCoords.csv", "DoseResponse.csv", "BMDs.csv"]:
        dblist.append(os.path.join(output_dir, f"zebrafishChem{ftype}"))
        dblist.append(os.path.join(output_dir, f"zebrafishSamp{ftype}"))
    return dblist
    # runSchemaCheck(dblist)


def runExposome(chem_id_file):
    """
    run exposome data pull
    """
    # cmd = f"Rscript exposome/exposome_summary_stats.R {chem_id_file}"
    cmd = f"python exposome/exposome_summary_stats.py {chem_id_file}"
    tqdm.write(cmd)
    os.system(cmd)
    return [os.path.join(OUTPUT_DIR, "exposomeGeneStats.csv")]


def runExpression(gex, chem, ginfo):
    """
    run expression parsing
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
    """
    run schema checking
    """
    ##TODO: make this work with internal calls
    for filename in dbfiles:
        classname = os.path.basename(filename).split(".")[0]
        cmd = f"linkml-validate --schema srpAnalytics.yaml {filename} --target-class {classname}"
        tqdm.write(cmd)
        os.system(cmd)


def main():
    """
    this wrapping script is placed into every docker image to pull the files
    from the repo and initiate the appropriate call to the underlying code.
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
            sampmap_args = {**smp, **sampmap_args}
            res = runSampMap(**sampmap_args)
            all_res.extend(res)
            progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Running sample mapping... Done!")
        progress_bar.close()

        # Collect all unique files
        all_res = list(set(all_res))
        for f in sample_files + chem_files:
            os.system(f"rm {f}")

        # Validate schema
        runSchemaCheck(res)

    # -----------------
    # Exposome Workflow
    # -----------------
    if args.expo:
        res = runExposome(cid)
        runSchemaCheck(res)

    # ------------------------
    # Gene Expression Workflow
    # ------------------------
    if args.geneEx:
        if not os.path.exists(os.path.join(OUTPUT_DIR, "chemicals.csv")):
            runSampMap(False, [], sid, cid, emap, cclass, fses, descfile)

        res = runExpression(gex1, os.path.join(OUTPUT_DIR, "chemicals.csv"), ginfo)
        runSchemaCheck(res)


main()
