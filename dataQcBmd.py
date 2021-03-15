#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time
import argparse
sys.path.insert(0, './qc_BMD')
#from qc_BMD import bmd_analysis_full as bmd
from qc_BMD import bmd_analysis_02 as bmd
from qc_BMD import wide2dicho_LPR_7_PAH_t0_t239_display_to_print_03 as bmd_LPR

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
parser.add_argument('--LPR', dest='LPR', type=os.path.abspath,\
                    help='LPR input csv file, needed for LPR data processing')

############ developer comment:
# for morphological data, only morphological data is needed as input
# for LPR processing, both morphological data and LPR data re needed as inputs


def merge_files(path,file_dict):
    '''
    merge_files takes a dictionary of files and joints them to a single file to
    added to the next step of the algorithm
    '''

    bmds = []
    fits = []
    dose = []
    for dataset, filelist in file_dict.items():
        bmds.append(filelist[0])
        fits.append(filelist[1])
        dose.append(filelist[2])
    pd.concat([pd.read_csv(f) for f in bmds]).to_csv(path+'/new_bmds.csv')
    pd.concat([pd.read_csv(f) for f in fits]).to_csv(path+'/new_fits.csv')
    pd.concat([pd.read_csv(f) for f in dose]).to_csv(path+'/new_dose.csv')
    return [path+'/new_bmds.csv',path+'/new_fits.csv',path+'/new_dose.csv']

if __name__ == "__main__":
    start_time = time.time()
    args = parser.parse_args()
    flist = args.files.split(',')
    print(flist)
    files = dict()
    if flist[0] == '':
        print("No new files, just re-building archive")
        command = "Rscript /srpAnalytics/buildv1database.R"
        os.system(command)
    else:
        for input_csv_file_name in flist:
            print ("input_csv_file_name:" + str(input_csv_file_name))

            #  if args.isSample:
            ##this doesn't exist yet - we need to read in sample information and merge with
            ##existing dose response data
            #      command = "Rscript 04_mergeWithChems.R "+str(input_csv_file_name)
            #      print(command)
            #      os.system(command)
            #  else:
            if args.devel:
                full_devel = "devel"
            else:
                full_devel = "full"

            if args.LPR == None:
                command = "python3 /srpAnalytics/01_reformat_df_data.py " + \
                    str(input_csv_file_name) + " " + str(full_devel)
                print(command)
                os.system(command)
            else:
                # for LPR reformatting both morphological and LPR is needed
                command = "python3 /srpAnalytics/01_reformat_df_data.py " + \
                    str(input_csv_file_name) + " " + str(full_devel)
                print(command)
                os.system(command)

                command = "python3 /srpAnalytics/01_reformat_df_LPR_data.py " + \
                    str(args.LPR) + " " + str(full_devel)
                print(command)
                os.system(command)

            output_complete_file_path = input_csv_file_name[:-4] + \
                "_wide_DNC_0.csv"
            # actual file is not saved here, but it is ok to be used at following procedures

            if args.LPR != None:
                output_complete_file_path_LPR = args.LPR[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"

            if args.LPR == None:
                files[input_csv_file_name] = bmd.runBmdPipeline(output_complete_file_path, full_devel)
            else:
                files[input_csv_file_name] = bmd_LPR.runBmdPipeline(output_complete_file_path, \
                                               output_complete_file_path_LPR, \
                                               full_devel)

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

    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print ("Done, it took:" + str(time_took))
