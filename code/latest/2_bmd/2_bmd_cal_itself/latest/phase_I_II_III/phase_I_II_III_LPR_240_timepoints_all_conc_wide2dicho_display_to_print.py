#!/usr/bin/env python
# coding: utf-8

# In[45]:


# [overall explanation about this code]
# just calculate BMD (already input file is prepared)
########## e.g. filter phase I, II LPR against phase I, II morpho
########## filter phase III LPR against phase III morpho
########## merge them


# In[46]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, random, shutil, time
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import sys

#mac
#util_path = "/Users/kimd999/research/script_not_in_dropbox/srpAnalytics/code/latest/util"

#constance
#'''
args = sys.argv[0:]
py_file = args[0]
py_file_wo_path = os.path.basename(py_file)

code_location = os.path.dirname(os.path.abspath(py_file))
index_of_latest = code_location.index('latest')
util_path = os.path.join(code_location[:index_of_latest], "latest", "util")
print ("util_path:"+ str(util_path))
#'''

sys.path.insert(0, util_path)


# In[47]:


starting_dir = os.getcwd()
print (starting_dir)


# In[56]:


# mac       - phase I, II, III - LPR - after_merging - 240 timepoints in min
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II_III/LPR/input/phase_I_II_III_LPR_240_timepoints_morpho_filtered_min.csv'
# 383 unique chemicals

# constance
complete_file_path = '/people/kimd999/tox/phase_I_II_III/LPR/input/phase_I_II_III_LPR_240_timepoints_morpho_filtered_min.csv'
# 383 unique chemicals

df_lpr_min = pd.read_csv(complete_file_path, header = 0)
df_lpr_min.head()
print(df_lpr_min.shape)

print(len(np.unique(df_lpr_min['chemical.id'])))


# In[58]:


## Calculate MOV, AUC for all chemical concentrations

delta_mov_auc = df_lpr_min[['chemical.id', 'conc', 'plate.id', 'well']].copy()

trans_points = [2,8,14,20] # "official"

end_points = ['MOV', 'AUC']

num_light = 3 # seems reasonable since interval between middle points of each peak ~= 6
num_dark  = 3

for trans_index, trans_point in enumerate(trans_points):
#    print ("\ntrans_index:" + str(trans_index))
    print ("trans_point:" + str(trans_point))
    
    for just_index, end_point in enumerate(end_points):
        if (end_point == 'MOV'):
            delta_mov_auc['MOV' + str(trans_index + 1)] = df_lpr_min['t' + str(trans_point + 1)] - df_lpr_min['t' + str(trans_point)]
        else:
            delta_mov_auc['AUC' + str(trans_index + 1)]             = sum(df_lpr_min['t' + str(trans_point + 1 + index_count)]                   for index_count in range(num_dark))             - sum(df_lpr_min['t' + str(trans_point - index_count)]                   for index_count in range(num_light))

print(delta_mov_auc.head())

print(delta_mov_auc.shape)

#cwd = os.getcwd()
#print (cwd)
#delta_mov_auc.to_csv("delta_mov_auc.csv", index=False)


# In[51]:


import generate_dose_response as gdr
import BMD_BMDL_estimation as bmdest
import Plot_Save as ps


# In[52]:


# This box is essential for BMD calculation
# Rename column headers to make it compatible with earlier data received from Lisa
delta_mov_auc.rename(columns={"chemical.id": "Chemical.ID", "conc": "CONC", "plate.id": "Plate", "well": "WELL"}, inplace = True)
print(delta_mov_auc.head())
#display(delta_mov_auc.tail())


# In[54]:


# Calculate BMD finally
start_time = time.time()
os.chdir(starting_dir)

if (os.path.isdir("output") == True):
    shutil.rmtree("output")
os.mkdir("output")

output_folder = os.path.join(starting_dir, "output")
os.chdir(output_folder)

full_devel = "full"
#full_devel = "devel"

if (full_devel == "full"):
    chemical_id_from_here = np.unique(delta_mov_auc['Chemical.ID'])
    end_points_from_here = ['MOV1','AUC1','MOV2','AUC2','MOV3','AUC3','MOV4','AUC4']
else:
    chemical_id_from_here = [53]
    end_points_from_here = ['MOV1']

#report = True
report = False

total_number_of_chemicals_to_processed = len(chemical_id_from_here)
number_of_chemicals_processed = 0

for chemical_id in chemical_id_from_here:
    print("\nchemical_id:" + str(chemical_id))
    for end_point in end_points_from_here:
        if (report): print("end_point:" + str(end_point))
        # subset original dataframe for a user-specified chemical and end_point pair
        delta_mov_auc_end_point_chemical_id = delta_mov_auc.loc[delta_mov_auc['Chemical.ID'] == chemical_id,['Chemical.ID', 'CONC', 'Plate', 'WELL', end_point]]
        #print("delta_mov_auc_end_point_chemical_id:\n"+str(delta_mov_auc_end_point_chemical_id))
        #print("type(delta_mov_auc_end_point_chemical_id):\n"+str(type(delta_mov_auc_end_point_chemical_id)))
        #print("type(end_point):\n"+str(type(end_point)))

        dose_response = gdr.gen_dose_response_behavior(delta_mov_auc_end_point_chemical_id, end_point)
        if (report): print("dose_response:\n"+str(dose_response))
        qc_flag = gdr.BMD_feasibility_analysis(dose_response)
        print ("qc_flag:"+str(qc_flag))
        test_dose_response = gdr.reformat_dose_response(dose_response)
        #test_dose_response = dose_response
        if(qc_flag in [0, 1]):
            # No BMD analysis required. Generate report and exit
            ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point, None)
        else:
            # Fit dose response models
            model_predictions = bmdest.analyze_dose_response_data(test_dose_response)
            # Select best model
            selected_model_params = bmdest.select_model(model_predictions)
            # Check if unique model is found
            unique_model_flag = selected_model_params['no_unique_model_found_flag']
            if(unique_model_flag == 0):
                # Generate report
                ps.save_results_good_data_unique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), end_point)
            else:
                bmd_analysis_flag = selected_model_params['model_select_flag']
                if(bmd_analysis_flag == 1):
                    ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point, selected_model_params)
                else:
                    ps.save_results_good_data_nounique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), end_point)
                    
    number_of_chemicals_processed += 1
    print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_processed)
    print(print_this)

cwd = os.getcwd()
print (cwd)

time_filename = 'running_time.txt'
f_time = open(time_filename, 'w')
f_time.write(str(time_took))
f_time.close()

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took))
# 1 chemical (3756) and 2 endpoints (['MOV1','AUC1']), 140 seconds took
# 7 chemicals and 2 endpoints (['MOV1','AUC1']), 6 minutes took
# [mac] 186 chemicals and 2 endpoints (['MOV1','AUC1']), 6 hrs took
# [constance] 186 chemicals and 2 endpoints (['MOV1','AUC1']), 4.5 hrs took

