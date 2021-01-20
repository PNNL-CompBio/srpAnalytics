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

parser = argparse.ArgumentParser('Run the QC and BMD analysis as well as join with \
    extract data to store in SRP data analytics portal')

#parser.add_argument('--label', dest='label', help='Label to store data', \
#                    default='newdata')
parser.add_argument('--isSample', dest='isSample', action='store_true',\
                    default=False, help='Set this flag if we \
                    are processing a sample not a chemical')
parser.add_argument('files', nargs='?', default='')
parser.add_argument('--devel', dest='devel',\
                    help='Set this flag to run test code instead of full analysis',\
                    action='store_true', default=False)

if __name__ == "__main__":
    args = parser.parse_args()
    flist = args.files.split(',')
    print(flist)

    if len(flist)==0:
        print("No new files, just re-building archive")
        command = "Rscript /srpAnalytics/03_mergeWithExtracts.R"
        os.system(command)
    else:
        for input_csv_file_name in flist:
            #  if args.isSample:
            ##this doesn't exist yet - we need to read in sample information and merge with
            ##existing dose response data
            #      command = "Rscript 04_mergeWithChems.R "+str(input_csv_file_name)
            #      print(command)
            #      os.system(command)
            #  else:
            command = "python3 /srpAnalytics/01_reformat_df_data.py " + str(input_csv_file_name)
            print(command)
            os.system(command)
            output_complete_file_path = input_csv_file_name[:-4] + "_wide_DNC_0.csv"

            #command = "python3 /srpAnalytics/qc_BMD/02_bmd_analysis.py " + \
                #    str(output_complete_file_path)
            if args.devel:
                fulldev = "devel"
            else:
                fulldev = "full"
            files = bmd.runBmdPipeline(output_complete_file_path, fulldev)

            command = "Rscript /srpAnalytics/03_mergeWithExtracts.R "
            if args.isSample:
                command = command+'--samples  '
            else:
                command = command+'--chemicals '
            if len(files) == 3:
                command = command + ','.join(files)
                print(command)
                os.system(command)
