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


##impor  BMD files from directory
from qc_BMD import bmd_analysis_morpho as bmd
from qc_BMD import bmd_analysis_LPR_7_PAH_t0_t239 as bmd_LPR

parser = argparse.ArgumentParser('Run the QC and BMD analysis as well as join with \
extract data to store in SRP data analytics portal')

#parser.add_argument('--label', dest='label', help='Label to store data', \
#                    default='newdata')
#parser.add_argument('files', nargs='?', default='',\
#                    help='Morphological files for regular BMD input or LPR (with --LPR option)')
parser.add_argument('--morpho',dest='morpho',\
                    help='Comma-delimited list of morphological files to be processed',\
                    default='')
parser.add_argument('--LPR', dest='lpr', \
                    help='Comma-delimited list of LPR-related files to be processed. MUST correspond to similar files in the morpho argument',\
                    default='')
parser.add_argument('--test-lpr', dest='test_lpr',\
                    help='Set this flag to run LPR test code instead of full analysis',\
                    action='store_true', default=False)
parser.add_argument('--test-morpho', dest='test_morpho',\
                    help='Set this flag to run LPR test code instead of full analysis',\
                    action='store_true', default=False)

parser.add_argument('--validate', dest='validate', \
                    help='If this tag is added, then we validate existing files',\
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


def run_lpr_on_file(lpr_file,morph_file, full_devel='full'):
    """
    runs LPR code on a file
    Attributes
    ----
    unformatted_file: str
    """

    command = "python3 /srpAnalytics/format_LPR_input.py " + \
        str(lpr_file) + " " + str(full_devel)
    print(command)
    os.system(command)

    LPR_input_csv_file_name_wide = lpr_file[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"

    #print ("LPR_input_csv_file_name_wide:" + str(LPR_input_csv_file_name_wide))

    #print ("morpho_input_csv_file_name:" + str(morph_file))
    #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall.csv
    command = "python3 /srpAnalytics/format_morpho_input.py " + \
                    str(morph_file) + " " + str(full_devel)
    print(command)
    os.system(command)
            #time.sleep(20)

    morpho_input_csv_file_name_wide = morph_file[:-4] + "_wide_DNC_0.csv"

    res = bmd_LPR.runBmdPipeline(morpho_input_csv_file_name_wide, \
                                             LPR_input_csv_file_name_wide, full_devel)
    return res

def run_morpho_on_file(morph_file,full_devel='full'):
    """
    formats and runs morphological BMD on file
    """
    print ("morpho_input_csv_file_name:" + str(morph_file))
    #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall.csv
    command = "python3 /srpAnalytics/format_morpho_input.py " + \
                    str(morph_file) + " " + str(full_devel)
    print(command)
    os.system(command)
            #time.sleep(20)

    morpho_input_csv_file_name_wide = morph_file[:-4] + \
        "_wide_DNC_0.csv"
    print ("morpho_input_csv_file_name_wide:" + str(morpho_input_csv_file_name_wide))
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
        os.system(command)
        for m in merged_files:
            os.system('rm '+m)


def main():
    """
    main method for command line
    """
    start_time = time.time()
    args = parser.parse_args()
    #flist = args.files.split(',')
    #print(flist)

    ##collecting a list of files to add to DB
    files = dict()

    mfiles = []
    lfiles = []
    if args.morpho!='':
        mfiles = args.morpho.split(',')
        print("Calculating morphological endpoints for "+str(len(mfiles))+' files')
        for f in mfiles:
            files[f] = run_morpho_on_file(f)
    if args.lpr!="":
        lfiles = args.lpr.split(',')
        if len(lfiles)!=len(mfiles):
            print("Cannot calculate LPR without morphological files, please re-run with --morpho argument")
            sys.exit()
        else:
            print( 'Calculating LPR endpoints for '+str(len(lfiles))+'LPR files')
            for i in range(len(files)):
                fname=lfiles[i]
                files[fname] = run_lpr_on_file(fname,mfiles[i])
    if args.morpho=="":
        if args.test_lpr:
            print("Testing LPR code")
            test_lpr = '/srpAnalytics/test_files/7_PAH_zf_LPR_data_2020NOV11_tall.csv'
            test_morph = '/srpAnalytics/test_files/7_PAH_zf_morphology_data_2020NOV11_tall.csv'
            res = run_lpr_on_file(test_lpr,test_morph,'devel')
        elif args.test_morpho:
            test_morph = '/srpAnalytics/test_files/7_PAH_zf_morphology_data_2020NOV11_tall.csv'
            print("Testing morphological code")
            res = run_morpho_on_file(test_morph,'devel')
    if len(files)==0:
        print("Testing database rebuild")
        command = "Rscript /srpAnalytics/buildv1database.R"
        os.system(command)

    else:
        print("building database with new files")
        build_db_with_files(files)
    if args.validate:
        print("Validating existing files")
        ##get files
        ##validate
#    else:
#        for morpho_input_csv_file_name in flist:
           #ull_devel = "full"
           # print ("full_devel:" + str(full_devel)) # devel
            ########## <begin> tall format (Oregon state original) -> wide format (so that BMD can be calculated)


            # actual file is not saved here, but it is ok to be used at following procedures


 #           print ("args.LPR:" + str(args.LPR)) # True
       #     if (args.LPR == True):
                # for LPR reformatting (tall->wide), both morphological and LPR is needed

        #          bypass_can = input()
                #devel

                #to_be_processed/7_PAH_zf_LPR_data_2021JAN11_tall_wide_t0_t239_devel.csv
            ########### <end> tall format (Oregon state original) -> wide format (so that BMD can be calculated)



            ########### <begin> BMD calculation
         #   if (args.LPR == True):
          #                                            full_devel)
           # else: # if (args.LPR == False):
            #    files[morpho_input_csv_file_name] =
            ########### <end> BMD calculation




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
