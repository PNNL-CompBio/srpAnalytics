#!/usr/bin/env python
# coding: utf-8

#7_PAH_zf_morphology
import sys

import numpy as np
import pandas as pd
import os, sys, time
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

starting_dir = os.getcwd()
print(starting_dir)

import BMD_BMDL_estimation as bmdest
import generate_dose_response as gdr
import Plot_Save as ps

def main():
    args = sys.argv[0:]
    complete_file_path = args[1]
    print("Running command line bmd analysis on "+complete_file_path)
    full_devel = args[2]
    files = runBmdPipeline(complete_file_path, full_devel)
    print(files)

def runBmdPipeline(complete_file_path, full_devel):
    ##read in the data
    df_morpho = pd.read_csv(complete_file_path, header = 0)


    ##let's collect the file names to return
    filenames = []
    test_data_sim = 0
    if(test_data_sim == 0):
        # Add aggregate endpoints
        # 1. Any effect at 24hrs (combination of MO24, DP24 and SM24) >> 'ANY24'
        df_morpho['ANY24'] = df_morpho[['MO24','DP24','SM24']].sum(axis=1,skipna=True,min_count=1)

        # 2. Any effect within 5 days (combination of all measurements at both time points)
        df_morpho['ANY120'] = df_morpho[['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24']].sum(axis=1,skipna=True,min_count=1)

        # 3. Total mortality (MO24 + MORT) >> 'TOT_MORT'
        df_morpho['TOT_MORT'] = df_morpho[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)

        # 4. Any effect except mortality (#2 minus MO24 and MORT) >> 'ALL_BUT_MORT'
        df_morpho['ALL_BUT_MORT'] = df_morpho[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', \
                                               'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)


    if os.path.isdir("output") is False:
        os.mkdir("output")

    output_folder = os.path.join(starting_dir, "output")
    os.chdir(output_folder)

    start_time = time.time()

    # full -> 17 (without DNC) unlike phase_I_II (18 endpoints), 7_PAH lacks NC24
    if full_devel == "full":
        end_points = ['ANY24', 'ANY120', 'AXIS', 'ALL_BUT_MORT', 'BRN_',\
                      'CRAN', 'DP24', 'EDEM', 'LTRK', 'MO24', 'MORT', 'MUSC',\
                      'NC__', 'SKIN', 'SM24', 'TCHR', 'TOT_MORT']
        chemical_id_from_here = np.unique(df_morpho['chemical.id'])
    else:
        end_points = ['ANY24']
        chemical_id_from_here = [3756]

    for chemical_id in chemical_id_from_here:
        #print("chemical_id:" + str(chemical_id))
        for end_point in end_points:
            os.chdir(output_folder)
            # subset original dataframe for a user-specified chemical and end_point pair
            df_morpho_end_point_chemical_id = df_morpho.loc[df_morpho['chemical.id'] == \
                chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]

            # Binarize end-point hits (Values > 1 are forced to 1)
            end_point_hits = df_morpho_end_point_chemical_id[end_point]
            end_point_hits.loc[end_point_hits > 0] = 1

            dose_response = gdr.gen_dose_response(df_morpho_end_point_chemical_id, end_point)

            qc_flag = gdr.BMD_feasibility_analysis(dose_response)

            test_dose_response = gdr.reformat_dose_response(dose_response)

            ##TODO: currently we format the files in individual methods;
            ## we should do that here, then pass the file in from here
            if(qc_flag in [0, 1]):
                # No BMD analysis required. Generate report and exit
                filenames = ps.save_results_poor_data_or_no_convergence(test_dose_response,\
                                                            qc_flag, str(chemical_id),\
                                                            end_point, None)
            else:
                # Fit dose response models
                model_predictions = bmdest.analyze_dose_response_data(test_dose_response)
                # Select best model
                selected_model_params = bmdest.select_model(model_predictions)
                # Check if unique model is found
                unique_model_flag = selected_model_params['no_unique_model_found_flag']
                if(unique_model_flag == 0):
                    # Generate report
                    # print(f"test_dose_response.dose:{test_dose_response.dose}")
                    # print(f"test_dose_response.dose[-1:]:{test_dose_response.dose[-1:]}")
                    # print(f"qc_flag:{qc_flag}")
                    # print(f"model_predictions:{model_predictions}")
                    # print(f"selected_model_params:{selected_model_params}")
                    # print(f"chemical_id:{chemical_id}")
                    # print(f"end_point:{end_point}")

                    filenames = ps.save_results_good_data_unique_model(test_dose_response,\
                                                           qc_flag, model_predictions,\
                                                           selected_model_params,\
                                                           str(chemical_id), end_point)
                else:
                    bmd_analysis_flag = selected_model_params['model_select_flag']
                    if(bmd_analysis_flag == 1):
                        filenames = ps.save_results_poor_data_or_no_convergence(test_dose_response, \
                                                                    qc_flag, str(chemical_id),\
                                                                    end_point, selected_model_params)
                    else:
                        filenames = ps.save_results_good_data_nounique_model(test_dose_response, qc_flag,\
                                                                 model_predictions, selected_model_params, \
                                                                 str(chemical_id), end_point)

    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print("Done, it took:"+str(time_took))
        # for all combinations of 342 chemicals and 18 endpoints, 4 minutes took for qc only
        # for all combinations of 342 chemicals and 18 endpoints, 104~165 minutes took for qc and bmd report

    os.chdir(starting_dir)
    #print (df_morpho_end_point_chemical_id)
    np.asarray(df_morpho_end_point_chemical_id['plate.id'])
    full_paths = [output_folder+'/'+f for f in filenames]
    return full_paths

if __name__ == '__main__':
    main()
