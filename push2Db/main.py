#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time
import argparse
import tarfile
import re
from ingest import pull_raw_data, test_connection
import validate as valid

OUT_FOLDER='/tmp'
IF_EXITS='replace' # options: "append", "replace", "fail"
DB='develop' # options: "develop", "production"


parser = argparse.ArgumentParser('Take the files and moved to a database')

parser.add_argument('--update-db', dest='update_db', action='store_true', \
                    help='Include --update-db if you want to update the database',\
                    default=False)


def merge_files(path, file_dict):
    """
    merge_files takes a dictionary of files and joints them to a single file to
    added to the next step of the algorithm

    Attributes
    ------
    path : str
    file_dict: dict
    """

    ## three lists of files to collect
    bmds = []
    fits = []
    dose = []
    for dataset, filelist in file_dict.items():
        bmds.append(filelist[0])
        fits.append(filelist[1])
        dose.append(filelist[2])

    ##concatenate all the files together
    pd.concat([pd.read_csv(f) for f in bmds]).to_csv(path+'/new_bmds.csv')
    pd.concat([pd.read_csv(f) for f in fits]).to_csv(path+'/new_fits.csv')
    pd.concat([pd.read_csv(f) for f in dose]).to_csv(path+'/new_dose.csv')
    return [path+'/new_bmds.csv',path+'/new_fits.csv',path+'/new_dose.csv']


def run_lpr_on_file(lpr_file, morph_file, full_devel='full'):
    """
    runs LPR code on a file
    Attributes
    ----
    unformatted_file: str
    """
    LPR_input_csv_file_name_wide = lpr_file[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"

    command = "python3 /srpAnalytics/format_LPR_input.py " + str(lpr_file) + " " + str(full_devel)

    if not os.path.exists(LPR_input_csv_file_name_wide):
        print(command)
        res0 = os.system(command)

    #print ("LPR_input_csv_file_name_wide:" + str(LPR_input_csv_file_name_wide))

    #print ("morpho_input_csv_file_name:" + str(morph_file))
    #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall.csv
    morpho_input_csv_file_name_wide = morph_file[:-4] + "_wide_DNC_0.csv"

    if not os.path.exists(morpho_input_csv_file_name_wide):
        command = "python3 /srpAnalytics/format_morpho_input.py " + str(morph_file) + " " + str(full_devel)
        print(command)
        res0 = os.system(command)
            #time.sleep(20)


    res = bmd_LPR.runBmdPipeline(morpho_input_csv_file_name_wide, \
                                             LPR_input_csv_file_name_wide, full_devel)
    return res

def run_morpho_on_file(morph_file, full_devel='full'):
    """
    formats and runs morphological BMD on file
    """
    print("morpho_input_csv_file_name:" + str(morph_file))
    #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall.csv
    command = "python3 /srpAnalytics/format_morpho_input.py " + str(morph_file) + " " + str(full_devel)
    print(command)
    os.system(command)

    morpho_input_csv_file_name_wide = morph_file[:-4] + \
        "_wide_DNC_0.csv"
    print("morpho_input_csv_file_name_wide:" + str(morpho_input_csv_file_name_wide))
    res = bmd.runBmdPipeline(morpho_input_csv_file_name_wide, \
                                             full_devel)
    return res

def build_db_with_files(fdict):
    """
    Builds database from dictoary of of filelists
    """
    merged_files = merge_files(os.getcwd(), fdict)
    command = "Rscript /srpAnalytics/buildv1database.R "

    #if args.isSample:
    #        command = command+'--samples  '
    #    else:
    command = command+'--chemicals '
    if len(merged_files) == 3:
        command = command + ','.join(merged_files)
        print(command)
        res0 = os.system(command)
        for m in merged_files:
            res0 = os.system('rm '+m)


def main():
    """
    main method for command line
    """
    start_time = time.time()
    args = parser.parse_args()
    #print(args)
    #flist = args.files.split(',')
    #print(flist)

    ##collecting a list of files to add to DB
    files = dict()

    if args.lpr is None:
        lfiles = ''
    else:
        lfiles = args.lpr.split(',')
    if args.morpho is None:
        mfiles = ''
    else:
        mfiles = args.morpho.split(',')

    print(lfiles)
    print(mfiles)
    fd = 'full'
    if args.test_lpr:
        print("Testing LPR code\n")
        lfiles = ['/srpAnalytics/test_files/7_PAH_zf_LPR_data_2021JAN11_3756.csv']
        mfiles = ['/srpAnalytics/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv']
        fd = 'devel'
 #       files['test'] = run_lpr_on_file(test_lpr, test_morph, 'devel')
    elif args.test_morpho:
        mfiles = ['/srpAnalytics/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv']
        print("Testing morphological code\n")
        fd = 'devel'
#        files['test'] = run_morpho_on_file(test_morph, 'devel')

    if len(lfiles) > 0:
        if len(lfiles) != len(mfiles):
            print("Cannot calculate LPR without morphological files, please re-run with --morpho argument")
            sys.exit()
        else:
            print('Calculating LPR endpoints for '+str(len(lfiles))+' LPR files')
            for i in range(len(lfiles)):
                fname = lfiles[i]
                files[fname] = run_lpr_on_file(fname, mfiles[i], fd)
    elif len(mfiles) > 0:
        print("Calculating morphological endpoints for "+str(len(mfiles))+' files')
        for f in mfiles:
            files[f] = run_morpho_on_file(f, fd)

#    if args.morpho == "":

    ##here we run the gene data collection
    if args.get_genes:
        print("Collecting data from BU")
        command = 'Rscript /srpAnalytics/exposome_summary_stats.R'
        os.system(command)

    ##here we run R code to merge all the files togeter
    if len(files) == 0:
        print("Testing database rebuild")
        command = "Rscript /srpAnalytics/buildv1database.R"
        os.system(command)
    else:
        print("building database with new files:")
        print(files)
        build_db_with_files(files)

    allfiles = ['/tmp/'+a for a in os.listdir('/tmp') if 'csv' in a]
    print(allfiles)
    if args.validate:
        print("Validating existing files for database ingest")
        ##get files
        for fval in allfiles:
            valid.verify(pd.read_csv(fval, quotechar='"', quoting=1), re.sub('.csv', '', os.path.basename(fval)))
        ##validate
    allfiles = allfiles+['/srpAnalytics/README.md']

    print('Now zipping up '+str(len(allfiles))+' files')
    tar = tarfile.open("/tmp/srpAnalyticsCompendium.tar.gz", "w:gz")
    for fname in allfiles:
        tar.add(fname)
    tar.close()

    if args.update_db:
        print('Saving to {}...'.format(DB))
        pull_raw_data(folder=OUT_FOLDER, if_exists=IF_EXITS, database=DB)
        print('Finished saving to database.')
  #  else: # if not saving to database, check connection to DB is okay
  #      print("Testing connection to database...", end='')
  #      okay, error = test_connection(database=DB)
  #      if okay:
  #          print('Connection OK')
  #      else:
  #          print('Connection failed, {}'.format(error))
    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print ("Done, it took:" + str(time_took))


if __name__ == "__main__":
    main()
