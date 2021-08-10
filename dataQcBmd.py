#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time
import argparse
import tarfile
from ingest import pull_raw_data, test_connection

OUT_FOLDER='/tmp'
IF_EXITS='replace' # options: "append", "replace", "fail"
DB='develop' # options: "develop", "production"


sys.path.insert(0, './qc_BMD')

from qc_BMD import bmd_analysis_morpho as bmd
from qc_BMD import bmd_analysis_LPR_7_PAH_t0_t239 as bmd_LPR

parser = argparse.ArgumentParser('Run the QC and BMD analysis as well as join with \
extract data to store in SRP data analytics portal')

#parser.add_argument('--label', dest='label', help='Label to store data', \
#                    default='newdata')
parser.add_argument('--isSample', dest='isSample', action='store_true',\
                    default=False, help='Set this flag if we are processing a sample not a chemical')
parser.add_argument('files', nargs='?', default='',\
                    help='Morphological files for regular BMD input or LPR (with --LPR option)')
parser.add_argument('--devel', dest='devel',\
                    help='Set this flag to run test code instead of full analysis',\
                    action='store_true', default=False)
parser.add_argument('--LPR', dest='LPR', \
                    help='If this tag is added, then LPR is calculated.',\
                    action='store_true', default=False)
parser.add_argument('--update-db', dest='update_db', action='store_true', help='Include --update-db if you want to update the database', default=False)

############ (developer comment)
# for morphological data, only morphological data is needed as input
# for LPR processing, both morphological data and LPR data are needed as inputs


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

if __name__ == "__main__":
    """
    main method for command line
    """
    start_time = time.time()
    args = parser.parse_args()
    flist = args.files.split(',')
    #print(flist)

    files = dict()
    if flist[0] == '':
        print("No new files, just re-building archive")
        command = "Rscript /srpAnalytics/buildv1database.R"
        os.system(command)
    else:
        for morpho_input_csv_file_name in flist:
            if args.devel:
                full_devel = "devel"
            else:
                full_devel = "full"
            print ("full_devel:" + str(full_devel)) # devel



            ########## <begin> tall format (Oregon state original) -> wide format (so that BMD can be calculated)
            print ("morpho_input_csv_file_name:" + str(morpho_input_csv_file_name))
            #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall.csv

            command = "python3 /srpAnalytics/format_morpho_input.py " + \
                    str(morpho_input_csv_file_name) + " " + str(full_devel)
            print(command)
            os.system(command)
            #time.sleep(20)

            morpho_input_csv_file_name_wide = morpho_input_csv_file_name[:-4] + \
                "_wide_DNC_0.csv"
            print ("morpho_input_csv_file_name_wide:" + str(morpho_input_csv_file_name_wide))
            # actual file is not saved here, but it is ok to be used at following procedures


            print ("args.LPR:" + str(args.LPR)) # True
            if (args.LPR == True):
                # for LPR reformatting (tall->wide), both morphological and LPR is needed

                LPR_input_csv_file_name = morpho_input_csv_file_name.replace("morphology", "LPR")

                command = "python3 /srpAnalytics/format_LPR_input.py " + \
                    str(LPR_input_csv_file_name) + " " + str(full_devel)
                print(command)
                os.system(command)

                LPR_input_csv_file_name_wide = LPR_input_csv_file_name[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"

                print ("LPR_input_csv_file_name_wide:" + str(LPR_input_csv_file_name_wide))

                #devel
                print ("press enter to continue")
                bypass_can = input()
                #devel

                #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall_wide_t0_t239_devel.csv
            ########### <end> tall format (Oregon state original) -> wide format (so that BMD can be calculated)



            ########### <begin> BMD calculation
            if (args.LPR == True):
                files[morpho_input_csv_file_name] = bmd_LPR.runBmdPipeline(morpho_input_csv_file_name_wide, \
                                             LPR_input_csv_file_name_wide, \
                                             full_devel)
            else: # if (args.LPR == False):
                files[morpho_input_csv_file_name] = bmd.runBmdPipeline(morpho_input_csv_file_name_wide, \
                                             full_devel)
            ########### <end> BMD calculation



        merged_files = merge_files(os.getcwd(), files)
        command = "Rscript /srpAnalytics/buildv1database.R "

        if args.isSample:
            command = command+'--samples  '
        else:
            command = command+'--chemicals '
        if len(merged_files) == 3:
            command = command + ','.join(merged_files)
            print(command)
            os.system(command)
            for m in merged_files:
                os.system('rm '+m)

         #wd <- paste0(getwd(),'/')
         ##UPDATE TO PYTHON     allfiles<-paste0(wd, c('README.md',list.files(path='.')[grep('csv',list.files(path='.'))]))
        allfiles = ['README.md'] + [a for a in os.listdir('/tmp') if 'csv' in a]
        print(allfiles)
        print('Now zipping up'+str(len(allfiles))+'files')
        tar = tarfile.open("/tmp/srpAnalyticsCompendium.tar.gz", "w:gz")
        for fname in allfiles:
            #os.system('mv '+fname+' /tmp')
            tar.add('/tmp/'+fname)
        tar.close()

    if args.update_db:
        print('Saving to {}...'.format(DB))
        pull_raw_data(folder=OUT_FOLDER, if_exists=IF_EXITS, database=DB)
        print('Finished saving to database.')
    else: # if not saving to database, check connection to DB is okay
        print("Testing connection to database...", end='')
        okay, error = test_connection(database=DB)
        if okay:
            print('Connection OK')
        else:
            print('Connection failed, {}'.format(error))
    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print ("Done, it took:" + str(time_took))
