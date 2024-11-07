"""
Build script moved to python for better extendability and interoperability.

"""

import argparse
import os
from typing import Optional
import pandas as pd
from tqdm import tqdm


def collectFiles(
    data_dir: str = "https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data",
    filename: str = "srp_build_files.csv",
) -> pd.DataFrame:
    """Collect files to be fed into each module upon build file update.

    Parameters
    ----------
    data_dir : str, optional
        _description_, by default "https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data"
    filename : str, optional
        _description_, by default "srp_build_files.csv"

    Returns
    -------
    pandas.DataFrame
        _description_
    """
    df = pd.read_csv(os.path.join(data_dir, filename))
    return df


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
    fulldf = pd.concat(dflist)
    fulldf = fulldf.drop_duplicates()
    return fulldf


def runSampMap(
    is_sample: bool = False,
    drcfiles: list = [],
    smap: str = "",
    cid: str = "",
    emap: str = "",
    cclass: str = "",
    ctfile: str = "",
    fses: str = "",
    descfile: str = "",
) -> list[str]:
    """
    run sample mapping
    """
    drc = ",".join(drcfiles)
    args = (
        f" --sampId={smap}"
        f" --chemid={cid}"
        f" --epMap={emap}"
        f" --chemClass={cclass}"
        f" --compToxFile={ctfile}"
        f" --sampleFiles={fses}"
        f" --chemDesc={descfile}"
        f" --sampMap={smap}"
    )
    if is_sample:
        cmd = f"Rscript sampleChemMapping/mapSamplesToChems.R --sample --drcFiles={drc} {args}"
    elif len(drcfiles) > 0:
        cmd = f"Rscript sampleChemMapping/mapSamplesToChems.R --chemical --drcFiles={drc} {args}"
    else:
        cmd = f"Rscript sampleChemMapping/mapSamplesToChems.R {args}"

    tqdm.write("\nRunning sample mapping with the following parameters:\n")
    tqdm.write(f"{cmd}\n")
    # os.system(cmd)
    tqdm.write("ls -la . \n")
    ##now we validate the files that came out.
    dblist = ["/tmp/samples.csv", "/tmp/chemicals.csv", "/tmp/sampleToChemicals.csv"]
    for ftype in ["XYCoords.csv", "DoseResponse.csv", "BMDs.csv"]:
        dblist.append(f"/tmp/zebrafishChem{ftype}")
        dblist.append(f"/tmp/zebrafishSamp{ftype}")
    return dblist
    # runSchemaCheck(dblist)


def runExposome(chem_id_file):
    """
    run exposome data pull
    """
    cmd = f"Rscript exposome/exposome_summary_stats.R {chem_id_file}"
    tqdm.write(cmd)
    os.system(cmd)
    return ["/tmp/exposomeGeneStats.csv"]


def runExpression(gex, chem, ginfo):
    """
    run expression parsing
    """
    cmd = f"Rscript zfExp/parseGexData.R {gex} {chem} {ginfo}"
    tqdm.write(cmd)
    os.system(cmd)
    return ["/tmp/srpDEGPathways.csv", "/tmp/srpDEGStats.csv", "/tmp/allGeneEx.csv"]


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
    df = collectFiles()

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
    sid = list(df.loc[df.name == "sampId"].location)[0]
    cid = list(df.loc[df.name == "chemId"].location)[0]
    cclass = list(df.loc[df.name == "class1"].location)[0]
    emap = list(df.loc[df.name == "endpointMap"].location)[0]
    fses = ",".join(list(df.loc[df.data_type == "sample"].location))
    ctfile = list(df.loc[df.name == "compTox"].location)[0]
    descfile = list(df.loc[df.name == "chemdesc"].location)[0]
    smap = list(df.loc[df.name == "sampMap"].location)[0]
    gex1 = ",".join(list(df.loc[df.data_type == "expression"].location))
    ginfo = list(df.loc[df.name == "geneInfo"].location)[0]

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
        chemfiles = []
        sampfiles = []

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
                fname = f"/tmp/tmp_{st}_{dt}.csv"
                fdf.to_csv(fname, index=False)
                if st == "chemical":
                    chemfiles.append(fname)
                else:
                    sampfiles.append(fname)
                progress_bar.update(1)

        # Update progress bar after completion
        progress_bar.set_description("Combining files... Done!")
        progress_bar.close()

        # Define fixed params for sample mapping
        all_res = list()
        sampmap_args = {
            "smap": smap,
            "cid": cid,
            "emap": emap,
            "cclass": cclass,
            "ctfile": ctfile,
            "fses": fses,
            "descfile": descfile,
        }

        # Iterate through sampMap params
        sampmap_params = [
            {"is_sample": True, "drcfiles": sampfiles},
            {"is_sample": False, "drcfiles": chemfiles},
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
        for f in sampfiles + chemfiles:
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
        if not os.path.exists("/tmp/chemicals.csv"):
            runSampMap(False, [], smap, cid, emap, cclass, ctfile, fses, descfile)

        res = runExpression(gex1, "/tmp/chemicals.csv", ginfo)
        runSchemaCheck(res)


main()
