#!/usr/bin/env python
# coding: utf-8

#7_PAH_zf_morphology
import sys

import numpy as np
import pandas as pd
import os, random, sys, time
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
    df_morph = pd.read_csv(complete_file_path, header = 0)


    ##let's collect the file names to return
    filenames = []
    test_data_sim = 0
    if(test_data_sim == 0):
        # Add aggregate endpoints
        # 1. Any effect at 24hrs (combination of MO24, DP24 and SM24) >> 'ANY24'
       # df_morph['ANY24'] = df_morph[['MO24','DP24','SM24']].sum(axis=1,skipna=True,min_count=1)

        # 2. Any effect within 5 days (combination of all measurements at both time points)
        #df_morph['ANY120'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24']].sum(axis=1,skipna=True,min_count=1)

        # 3. Total mortality (MO24 + MORT) >> 'TOT_MORT'
        #df_morph['TOT_MORT'] = df_morph[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)

        # 4. Any effect except mortality (#2 minus MO24 and MORT) >> 'ALL_BUT_MORT'
        #df_morph['ALL_BUT_MORT'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', \
                                               #'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)
        # Add aggregate endpoints for 7 PAH

        ## JUSTIFICATION: 7 PAH dataset doesn't have "BRAI" endpoint.
        ## On the other hand, extract/phase I,II have "BRAI" endpoint.
##this was added ina merge conflict, i hope it's ok (Sara)
        if 'BRAI' not in df_morph.columns: # as 7 PAH
            # 1. Any effect at 24hrs (combination of MO24, DP24 and SM24) >> 'ANY24'
            df_morph['ANY24'] = df_morph[['MO24','DP24','SM24']].sum(axis=1,skipna=True,min_count=1)
            # 2. Any effect within 5 days (combination of all measurements at both time points)
            df_morph['ANY120'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24']].sum(axis=1,skipna=True,min_count=1)

            # 3. Total mortality (MO24 + MORT) >> 'TOT_MORT'
            df_morph['TOT_MORT'] = df_morph[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)

            # 4. Any effect except mortality (#2 minus MO24 and MORT) >> 'ALL_BUT_MORT'
            df_morph['ALL_BUT_MORT'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM', \
                                                   'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)



        else: # as extract
            df_morph['ANY24'] = df_morph[['MO24','DP24','SM24','NC24']].sum(axis=1,skipna=True,min_count=1)
            df_morph['ANY120'] = df_morph[['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', \
                                                               'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC', \
                                                               'TRUN', 'SWIM', 'NC__', 'TR__', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
            df_morph['TOT_MORT'] = df_morph[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)
            df_morph['ALL_BUT_MORT'] = df_morph[['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC', \
                                                               'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__', \
                                                               'TR__']].sum(axis=1,skipna=True,min_count=1)
            df_morph['BRN_'] = df_morph[['BRAI','OTIC','PFIN']].sum(axis=1,skipna=True,min_count=1)
            df_morph['CRAN'] = df_morph[['EYE_', 'SNOU', 'JAW_']].sum(axis=1,skipna=True,min_count=1)
            df_morph['EDEM'] = df_morph[['YSE_','PE__']].sum(axis=1,skipna=True,min_count=1)
            df_morph['LTRK'] = df_morph[['TRUN','CFIN']].sum(axis=1,skipna=True,min_count=1)
            df_morph['MUSC'] = df_morph[['CIRC','SWIM','SOMI']].sum(axis=1,skipna=True,min_count=1)
            df_morph['SKIN'] = df_morph[['PIG_']]
            df_morph['TCHR'] = df_morph[['TR__']]

### END Of merge conflict

    if os.path.isdir("output") is False:
        os.mkdir("output")

    output_folder = os.path.join(starting_dir, "output")
    os.chdir(output_folder)

    start_time = time.time()

    # full -> 17 (without DNC) unlike phase_I_II (18 endpoints), 7_PAH lacks NC24
    if full_devel == "full":
        if 'BRAI' not in df_morph.columns: # as 7 PATH
            end_points = ['ANY24', 'ANY120', 'AXIS', 'ALL_BUT_MORT', 'BRN_',\
                      'CRAN', 'DP24', 'EDEM', 'LTRK', 'MO24', 'MORT', 'MUSC',\
                      'NC__', 'SKIN', 'SM24', 'TCHR', 'TOT_MORT']
        else:
            end_points = ['ANY24','ANY120','AXIS','ALL_BUT_MORT','BRN_','CRAN','DP24','EDEM','LTRK','MO24','MORT','MUSC','NC24','NC__','SKIN','SM24','TCHR','TOT_MORT']
        chemical_id_from_here = np.unique(df_morph['chemical.id'])
    else:
        end_points = ['ANY24']
        chemical_id_from_here = random.sample(set(np.unique(df_morph['chemical.id'])), 1)

    for chemical_id in chemical_id_from_here:
        #print("chemical_id:" + str(chemical_id))
        for end_point in end_points:
            os.chdir(output_folder)
            # subset original dataframe for a user-specified chemical and end_point pair
            df_morph_end_point_chemical_id = df_morph.loc[df_morph['chemical.id'] == \
                chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]

            # Binarize end-point hits (Values > 1 are forced to 1)
            end_point_hits = df_morph_end_point_chemical_id[end_point]
            end_point_hits.loc[end_point_hits > 0] = 1

            dose_response = gdr.gen_dose_response(df_morph_end_point_chemical_id, end_point)

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
    #print (df_morph_end_point_chemical_id)
    np.asarray(df_morph_end_point_chemical_id['plate.id'])
    full_paths = [output_folder+'/'+f for f in filenames]
    return full_paths

if __name__ == '__main__':
    main()
