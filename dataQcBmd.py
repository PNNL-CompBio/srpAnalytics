#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time
import argparse

parser = argparse.ArgumentParser('Run the QC and BMD analysis as well as join with \
    extract data to store in SRP data analytics portal')

parser.add_argument('--label', dest='label', help='Label to store data', \
                    default='newdata')
parser.add_argument('--isSample', dest='isSample', action='store_true',\
                    help='Set this flag if we are processing a sample not a chemical')
parser.add_argument('files', nargs='?', default='')

if __name__ == "__main__":
    args = parser.parse_args()
    flist = args.files.split(',')
    print(flist)
#    args = sys.argv[0:]
#    if (len(args) != 2):
#        print ("Specify <tall_raw.csv> \n")
#        print ("For tutorial, python run_all.py input/7_PAH_zf_morphology_data_2020NOV11_tall.csv\n")
#        exit(1)

    for input_csv_file_name in flist:
        if args.isSample:
            ##this doesn't exist yet - we need to read in sample information and merge with 
            ##existing dose response data
            command = "Rscript 04_mergeWithChems.R "+str(input_csv_file_name)
            print(command)
            os.system(command)
        else:
            command = "python3 01_reformat_zf_data.py " + str(input_csv_file_name)
            print(command)
            os.system(command)

            output_complete_file_path = input_csv_file_name[:-4] + "_wide_DNC_0.csv"

            command = "python3 qc_BMD/02_bmd_analysis.py " + str(output_complete_file_path)
            print(command)
            os.system(command)
            #added in Rscript
            command = "Rscript 03_mergeWithExtracts.R"
            print(command)
            os.system(command)
