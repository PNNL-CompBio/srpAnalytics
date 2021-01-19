#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time

from datetime import datetime as dt

today = dt.now()  
time_now_date = today.strftime('%Y_%m_%d')


if (__name__ == "__main__"):
    args = sys.argv[0:]
    if (len(args) != 2):
        print ("Specify <tall_raw.csv> \n")
        print ("For tutorial, python run_all.py input/7_PAH_zf_morphology_data_2021JAN11_3756_only_for_devel.csv\n")
        exit(1)

    input_csv_file_name = args[1]
    
    command = "python 1_preprocessing/reformat_most_7_PAHs_arg_added.py " + str(input_csv_file_name)
    print (command)
    os.system(command)
    
       
    input_csv_file_name_base = os.path.basename(input_csv_file_name)
    output_file_name = input_csv_file_name_base[:-4] + "_ran_at_" + str(time_now_date) + '.csv'
    output_file_name_w_path = os.path.join("input", output_file_name)
    
    
    command = "python 2_qc_BMD/Main_Wrapper_BMD_Analysis_7_PAHs.py " + str(output_file_name_w_path)
    print (command)
    os.system(command)
    