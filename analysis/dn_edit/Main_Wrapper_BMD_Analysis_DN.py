#!/usr/bin/env python
# coding: utf-8

# In[1]:

import datetime
import numpy as np
import pandas as pd
import os
from scipy import stats
from matplotlib import pyplot as plt

import generate_dose_response as gdr
import BMD_BMDL_estimation as bmdest
import Plot_Save as ps

import warnings
warnings.filterwarnings('ignore')


# In[2]:


os.getcwd()

#complete_file_path = '/Users/kimd999/Dropbox/script/python/toxicology/DN_try/Phase_I_II.csv'
complete_file_path = '/Users/kimd999/Dropbox/script/python/srpAnalytics/analysis/dn_edit/Phase_I_II.csv'

morphological_data = pd.read_csv(complete_file_path, header = 0)
#display(morphological_data.head())
#display(morphological_data.columns)
#display(np.unique(morphological_data.well))

### clear existing result csv files.
if (os.path.isfile("bmd_vals.csv") == True):
    os.remove("bmd_vals.csv")
if (os.path.isfile("fit_vals.csv") == True):
    os.remove("fit_vals.csv")
if (os.path.isfile("dose_response_vals.csv") == True):
    os.remove("dose_response_vals.csv")


# In[3]:


test_data_sim = 0
if(test_data_sim == 0):
    # Add aggregate endpoints
    # 1. Any effect at 24hrs (combination of MO24, DP24 and SM24) >> 'ANY24'
    # 2. Any effect within 5 days (combination of all measurements at both time points)
    # 3. Total mortality (MO24 + MORT) >> 'TOT_MORT'
    # 4. Any effect except mortality (#2 minus MO24 and MORT) >> 'ANY_MORT'
    # Add new endpoints
    # BRAIN	OTIC	PFIN >> 'BRN_'
    # EYE	SNOUT	JAW >> 'CRAN'
    # YSE	PE >> 'EDEM'
    # TRUNK	CFIN >> 'LTRK'
    # CIRC	SWIM	SOMITE >> 'MUSC'
    # PIG_ >> 'SKIN'
    # TR_ >> 'TCHR'
    morphological_data['ANY24'] = morphological_data[['MO24','DP24','SM24','NC24']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['ANY120'] = morphological_data[['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC',                                                        'TRUN', 'SWIM', 'NC__', 'TR__', 'DNC_', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['TOT_MORT'] = morphological_data[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['ANY_MORT'] = morphological_data[['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__',                                                        'TR__', 'DNC_']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['BRN_'] = morphological_data[['BRAI','OTIC','PFIN']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['CRAN'] = morphological_data[['EYE_', 'SNOU', 'JAW_']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['EDEM'] = morphological_data[['YSE_','PE__']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['LTKR'] = morphological_data[['TRUN','CFIN']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['MUSC'] = morphological_data[['CIRC','SWIM','SOMI']].sum(axis=1,skipna=True,min_count=1)
    morphological_data['SKIN'] = morphological_data[['PIG_']]
    morphological_data['TCHR'] = morphological_data[['TR__']]


# In[10]:

today = datetime.datetime.today()
date_based_working_folder = today.strftime('%Y%m%d') + "_" + str(today.hour) \
                            + "_" + str(today.minute) + "_" + str(today.second)
output_folder = os.path.join("output", date_based_working_folder)
output_folder_abs_path = os.path.abspath(output_folder)
print ("output_folder_abs_path:" + str(output_folder_abs_path))
if (os.path.isdir(output_folder_abs_path) == False):
    os.mkdir(output_folder_abs_path)
    
# Specify end_point and chemical of interest
# *********************************************
# Perform a check of the existence of "essential" column labels
# *********************************************
#end_points = ['ANY24','ANY120','TOT_MORT','ANY_MORT','BRN_','CRAN','EDEM','LTKR','MUSC','SKIN','TCHR']
end_points = ['AXIS','NC__','MO24','DP24','SM24','MORT']
#end_points = ['ANY24']
#for chemical_id in np.unique(morphological_data['chemical.id']):
for chemical_id in [53, 54]:
    print(chemical_id)
    for end_point in end_points:
        print(end_point)
        # subset original dataframe for a user-specified chemical and end_point pair
        morphological_data_end_point_chemical_id = morphological_data.loc[morphological_data['chemical.id'] == chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]

        # Binarize end-point hits (Values > 1 are forced to 1)
        end_point_hits = morphological_data_end_point_chemical_id[end_point]
        end_point_hits.loc[end_point_hits > 0] = 1

        dose_response = gdr.gen_dose_response(morphological_data_end_point_chemical_id, end_point)
        qc_flag = gdr.BMD_feasibility_analysis(dose_response)
        test_dose_response = gdr.reformat_dose_response(dose_response)
        print("test_dose_response"+str(test_dose_response))

        if(qc_flag in [0, 1]):
            # No BMD analysis required. Generate report and exit
            print ("qc_flag in [0, 1]")
            ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point, output_folder_abs_path, None)
        else:
            # Fit dose response models
            model_predictions = bmdest.analyze_dose_response_data(test_dose_response)
            # Select best model
            selected_model_params = bmdest.select_model(model_predictions)
            # Check if unique model is found
            unique_model_flag = selected_model_params['no_unique_model_found_flag']
            if(unique_model_flag == 0):
                # Generate report
                print(test_dose_response.dose[-1:])
                ps.save_results_good_data_unique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), output_folder_abs_path, end_point)
            else:
                bmd_analysis_flag = selected_model_params['model_select_flag']
                if(bmd_analysis_flag == 1):
                    print ("bmd_analysis_flag = 1")
                    exit(1)
                    ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point, output_folder_abs_path, selected_model_params)
                else:
                    print ("save_results_good_data_nounique_model")
                    exit(1)
                    ps.save_results_good_data_nounique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), output_folder_abs_path, end_point)


# In[ ]:


#display(morphological_data.head())
morphological_data.columns


# In[ ]:


test_dose_response.dose


# In[ ]:


test_dose_response.dose.iloc[0]+test_dose_response.dose.iloc[1]


# In[ ]:


dose_response['num_affect']/dose_response['num_embryos']


# In[ ]:




