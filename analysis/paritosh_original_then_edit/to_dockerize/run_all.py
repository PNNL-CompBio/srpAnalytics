#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time

if (__name__ == "__main__"):
    args = sys.argv[0:]
    if (len(args) != 2):
        print ("Specify tall_raw.csv \n")
        print ("For example, python run_all.py input/7_PAH_zf_morphology_data_2020NOV11_tall.csv\n")
        exit(1)

    input_csv_file_name = args[1]

    command = "python preprocessing/DN_reformat_most_7_PAHs.py " + str(input_csv_file_name)
    print (command)
    os.system(command)
    
    output_complete_file_path = input_csv_file_name[:-4] + "_wide_DNC_0.csv"
    
    command = "python qc_BMD/Main_Wrapper_BMD_Analysis_7_PAHs.py " + str(output_complete_file_path)
    print (command)
    os.system(command)
    