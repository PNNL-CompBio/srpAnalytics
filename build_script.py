"""
Build script moved to python for better extendability and interoperability.

"""

import os
import pandas as pd
import argparse


def collectFiles(
    data_dir="https://raw.githubusercontent.com/PNNL-CompBio/srpAnalytics/main/data",
    filename="srp_build_files.csv",
):
    """
    every time the build file is updated, this script will collect the files and return
    a dictionary of files to be fed into each module
    """
    df = pd.read_csv(data_dir + "/" + filename)
    return df


def fitCurveFiles(morpho_behavior_tuples):
    """
    get new curve fits, list of tuples of morpho/behavior pairs
    """


def combineFiles(location_list, ftype):
    """
    helper function to combine duplicates
    """
    dflist = []
    required_columns = {
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

    print("concatenating " + ftype)
    for loc in location_list.location:
        f = pd.read_csv(loc)[required_columns[ftype]]
        dflist.append(f)
    fulldf = pd.concat(dflist)
    fulldf = fulldf.drop_duplicates()

    return fulldf.drop_duplicates()


def runSampMap(
    is_sample=False,
    drcfiles=[],
    smap="",
    cid="",
    emap="",
    cclass="",
    ctfile="",
    fses="",
    desfile="",
):
    """
    run sample mapping
    """
    if is_sample:
        cmd = (
            "Rscript sampleChemMapping/mapSamplesToChems.R --sample --drcFiles="
            + ",".join(drcfiles)
            + " --sampId="
            + smap
            + " --chemId="
            + cid
            + " --epMap="
            + emap
            + " --chemClass="
            + cclass
            + " --compToxFile="
            + ctfile
            + " --sampleFiles="
            + fses
            + " --chemDesc="
            + desfile
            + " --sampMap="
            + smap
        )
    elif len(drcfiles) > 0:
        cmd = (
            "Rscript sampleChemMapping/mapSamplesToChems.R --chemical --drcFiles="
            + ",".join(drcfiles)
            + " --sampId="
            + smap
            + " --chemId="
            + cid
            + " --epMap="
            + emap
            + " --chemClass="
            + cclass
            + " --compToxFile="
            + ctfile
            + " --sampleFiles="
            + fses
            + " --chemDesc="
            + desfile
            + " --sampMap="
            + smap
        )
    else:
        cmd = (
            "Rscript sampleChemMapping/mapSamplesToChems.R --sampId="
            + smap
            + " --chemId="
            + cid
            + " --epMap="
            + emap
            + " --chemClass="
            + cclass
            + " --compToxFile="
            + ctfile
            + " --sampleFiles="
            + fses
            + " --chemDesc="
            + desfile
            + " --sampMap="
            + smap
        )

    print(cmd)
    os.system(cmd)
    print("ls -la .")
    ##now we validate the files that came out.
    dblist = ["/tmp/samples.csv", "/tmp/chemicals.csv", "/tmp/sampleToChemicals.csv"]
    for ftype in ["XYCoords.csv", "DoseResponse.csv", "BMDs.csv"]:
        dblist.append("/tmp/zebrafishChem" + ftype)
        dblist.append("/tmp/zebrafishSamp" + ftype)
    return dblist
    # runSchemaCheck(dblist)


def runExposome(chem_id_file):
    """
    run exposome data pull
    """
    cmd = "Rscript exposome/exposome_summary_stats.R " + chem_id_file
    print(cmd)
    os.system(cmd)
    return ["/tmp/exposomeGeneStats.csv"]


def runExpression(gex, chem, ginfo):
    """
    run expression parsing
    """
    cmd = "Rscript zfExp/parseGexData.R " + gex + " " + chem + " " + ginfo
    print(cmd)
    os.system(cmd)
    return ["/tmp/srpDEGPathways.csv", "/tmp/srpDEGStats.csv", "/tmp/allGeneEx.csv"]


def runSchemaCheck(dbfiles=[]):
    """
    run schema checking
    """
    ##TODO: make this work with internal calls
    for filename in dbfiles:
        classname = os.path.basename(filename).split(".")[0]
        cmd = (
            "linkml-validate --schema srpAnalytics.yaml "
            + filename
            + " --target-class "
            + classname
        )
        print(cmd)
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
    desfile = list(df.loc[df.name == "chemdesc"].location)[0]
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
        print("Re-running benchmark dose collection")
        newbmds, newfits, newdoses = [], [], []
        fitCurveFiles()

    if args.bmd or args.samps:  ### need to rerun samples if we have created new bmds
        # add chemical BMDS, fits, curves to existing data
        chemfiles = []
        sampfiles = []
        # print(fses)
        for st in ["chemical", "extract"]:
            for dt in ["bmd", "fit", "dose"]:
                fdf = combineFiles(
                    df.loc[df.sample_type == st].loc[df.data_type == dt], dt
                )
                fname = "/tmp/tmp_" + st + "_" + dt + ".csv"
                fdf.to_csv(fname, index=False)
                if st == "chemical":
                    chemfiles.append(fname)
                else:
                    sampfiles.append(fname)
        res1 = runSampMap(
            True, sampfiles, smap, cid, emap, cclass, ctfile, fses, desfile
        )
        res2 = runSampMap(
            False, chemfiles, smap, cid, emap, cclass, ctfile, fses, desfile
        )
        res3 = runSampMap(False, [], smap, cid, emap, cclass, ctfile, fses, desfile)
        res = res1 + res2 + res3
        res = list(set(res))
        for f in sampfiles + chemfiles:
            os.system("rm " + f)
        ##now we run validation
        runSchemaCheck(res)
    if args.expo:
        res = runExposome(cid)
        runSchemaCheck(res)
    if args.geneEx:
        if not os.path.exists("/tmp/chemicals.csv"):
            runSampMap(False, [], smap, cid, emap, cclass, ctfile, fses, desfile)

        res = runExpression(gex1, "/tmp/chemicals.csv", ginfo)
        runSchemaCheck(res)


main()
