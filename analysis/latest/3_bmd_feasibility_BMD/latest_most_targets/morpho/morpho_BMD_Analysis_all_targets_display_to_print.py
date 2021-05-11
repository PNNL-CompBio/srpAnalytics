#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os, random, shutil, sys, time
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import sys

#mac
#util_path = "/Users/kimd999/research/script_not_in_dropbox/srpAnalytics/analysis/latest/util"

#constance
#'''
args = sys.argv[0:]
py_file = args[0]

code_location = os.path.dirname(os.path.abspath(py_file))
index_of_latest = code_location.index('latest')
util_path = os.path.join(code_location[:index_of_latest], "latest", "util")
#'''

sys.path.insert(0, util_path)


# In[2]:


starting_dir = os.getcwd()
print (starting_dir)


# In[3]:


# mac       - 7 PAH - full
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/7_PAH/morpho/input/wide/7_PAH_zf_morphology_data_2021JAN11_wide_made_in_2021_01_19_DNC_0.csv'

# mac       - extracts - full
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/extracts/input/wide/extracts_morpho_PNNL_10-28-2020_005958_wide_DNC_0.csv'
# 42 unique chemical IDs

# mac       - phase I && II - morpho - devel
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/input/morpho/devel/wide/Tanguay_Phase_3_zf_morphology_data_PNNL_2021MAR23_devel_w_23_endpoints_wide_DNC_0_devel.csv'

# mac       - phase III - morpho - full - 23 endpoints
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/input/morpho/full/wide/Tanguay_Phase_3_zf_morphology_data_PNNL_2021MAR23_full_w_23_endpoints_wide_DNC_0_full.csv'
# 118 chemicals

# constance - 7 PAH - full
#complete_file_path = '/people/kimd999/tox/7_PAH/input/7_PAH_zf_morphology_data_2020NOV11_wide_DNC_0.csv'

# constance - extracts - full
#complete_file_path = '/qfs/people/kimd999/tox/extracts/input/tall_extracts_morphology_PNNL_10-28-2020_005958_wide_DNC_0.csv'

# constance - phase III - full
complete_file_path = '/people/kimd999/tox/phase_III/input/Tanguay_Phase_3_zf_morphology_data_PNNL_2021MAR23_full_w_23_endpoints_wide_DNC_0_full.csv'

df_morph = pd.read_csv(complete_file_path, header = 0)
pd.set_option('display.max_columns', None)
print(df_morph.head())
#display(np.unique(morpho_data.well))


# In[4]:


unique_chemical_id_s = np.unique(df_morph['chemical.id'])
print (len(unique_chemical_id_s))
print(df_morph.columns)
print(len(df_morph.columns))


# In[5]:


# Add super-endpoints
if ("PAH" in complete_file_path):
    # 1. Any effect at 24hrs (combination of MO24, DP24 and SM24) >> 'ANY24'
    df_morph['ANY24'] = df_morph[['MO24','DP24','SM24']].sum(axis=1,skipna=True,min_count=1)
    
    # 2. Any effect within 5 days (combination of all measurements at both time points)
    df_morph['ANY120'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'EDEM', 'LTRK', 'MORT', 'MUSC', 'NC__', 'SKIN', 'TCHR', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
    
    # 3. Total mortality (MO24 + MORT) >> 'TOT_MORT'
    df_morph['TOT_MORT'] = df_morph[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)

    # 4. Any effect except mortality (#2 minus MO24 and MORT) >> 'ALL_BUT_MORT'
    df_morph['ALL_BUT_MORT'] = df_morph[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM',                                                              'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)
else:
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
    df_morph['ANY24'] = df_morph[['MO24','DP24','SM24','NC24']].sum(axis=1,skipna=True,min_count=1)
    df_morph['ANY120'] = df_morph[['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC',                                                        'TRUN', 'SWIM', 'NC__', 'TR__', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
    df_morph['TOT_MORT'] = df_morph[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)
    df_morph['ALL_BUT_MORT'] = df_morph[['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__',                                                        'TR__']].sum(axis=1,skipna=True,min_count=1)
    df_morph['BRN_'] = df_morph[['BRAI','OTIC','PFIN']].sum(axis=1,skipna=True,min_count=1)
    df_morph['CRAN'] = df_morph[['EYE_', 'SNOU', 'JAW_']].sum(axis=1,skipna=True,min_count=1)
    df_morph['EDEM'] = df_morph[['YSE_','PE__']].sum(axis=1,skipna=True,min_count=1)
    df_morph['LTRK'] = df_morph[['TRUN','CFIN']].sum(axis=1,skipna=True,min_count=1)
    df_morph['MUSC'] = df_morph[['CIRC','SWIM','SOMI']].sum(axis=1,skipna=True,min_count=1)
    df_morph['SKIN'] = df_morph[['PIG_']]
    df_morph['TCHR'] = df_morph[['TR__']]
print(df_morph.columns)
print(len(df_morph.columns))


# In[6]:


print(df_morph.head())


# In[7]:


if (os.path.isdir("output") == True):
    shutil.rmtree("output")
os.mkdir("output")

output_folder = os.path.join(starting_dir, "output")
os.chdir(output_folder)

if (os.path.isdir("report") == False):
    os.mkdir("report")
    
df_morph_filename = os.path.join("report", 'df_morpho_after_merging_endpoints.csv')
df_morph.to_csv(df_morph_filename, index=False)


# In[8]:


import generate_dose_response as gdr

import BMD_BMDL_estimation as bmdest
import Plot_Save as ps


# In[9]:


# BMD calculation
start_time = time.time()

os.chdir(starting_dir)

os.chdir(output_folder)

bmd_feasibility_flag_filename = os.path.join("report", 'bmd_feasibility_flag.csv')
print ("bmd_feasibility_flag_filename:" + str(bmd_feasibility_flag_filename))

bmd_feasibility_flag_file_out = open(bmd_feasibility_flag_filename, "w")

write_this = "bmd_feasibility_flag\n"
bmd_feasibility_flag_file_out.write(write_this)

full_devel = "full"
#full_devel = "devel"

chemical_id_from_here = np.unique(df_morph['chemical.id'])

if (full_devel == "full"):
    if ("PAH" in complete_file_path):
        end_points = ['ANY24','ANY120','AXIS','ALL_BUT_MORT','BRN_','CRAN','DP24','EDEM','LTRK','MO24','MORT','MUSC','NC__', 'SKIN','SM24','TCHR','TOT_MORT']
    else: # full_oregon_state_request -> 18 (without DNC)
        end_points = ['ANY24','ANY120','AXIS','ALL_BUT_MORT','BRN_','CRAN','DP24','EDEM',                  'LTRK','MO24','MORT','MUSC','NC24','NC__','SKIN','SM24','TCHR','TOT_MORT']
else:
    end_points = ['ANY24']
    choose_this_number = min(len(chemical_id_from_here), 1)
    randomly_chosen = random.sample(set(chemical_id_from_here), choose_this_number)
    chemical_id_from_here = []
    for i in range(len(randomly_chosen)):
        chemical_id_from_here.append(randomly_chosen[i])
    
total_number_of_chemicals_to_processed = len(chemical_id_from_here)
number_of_chemicals_processed = 0

for chemical_id in chemical_id_from_here:
    print("\nchemical_id:" + str(chemical_id))
    for end_point in end_points:
        os.chdir(output_folder)
        # subset original dataframe for a user-specified chemical and end_point pair
        df_per_chemical = df_morph.loc[df_morph['chemical.id'] == chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]
        
        # Binarize end-point hits (Values > 1 are forced to 1)
        end_point_hits = df_per_chemical[end_point]
        end_point_hits.loc[end_point_hits > 0] = 1
                  
        dose_response = gdr.gen_dose_response(df_per_chemical, end_point)
        
        bmd_feasibility_flag = gdr.BMD_feasibility_analysis(dose_response)
        bmd_feasibility_flag_file_out.write(str(bmd_feasibility_flag)+"\n")

        #print ("dose_response:" + str(dose_response))        
        '''dose  num_affect  frac_affect  num_embryos  tot_wells
        0   0.0         0.0     0.000000         26.0       32.0
        1   0.1         1.0     0.032258         31.0       32.0
        2   0.5         1.0     0.062500         16.0       32.0
        '''
        
        test_dose_response = gdr.reformat_dose_response(dose_response)
        #print ("test_dose_response:" + str(test_dose_response))
        '''dose  num_affected  total_num
        0   0.0           0.0       26.0
        1   0.1           1.0       31.0
        2   0.5           1.0       16.0
        '''
        
        #bmd_feasibility_flag_folder = "bmd_feasibility_" + str(bmd_feasibility_flag)
        #if (os.path.isdir(str(bmd_feasibility_flag_folder)) == False):
        #    os.mkdir(str(bmd_feasibility_flag_folder))
        #os.chdir(str(bmd_feasibility_flag_folder))

        if(bmd_feasibility_flag in [0, 1]):
            # No BMD analysis required. Generate report and exit
            ps.save_results_poor_data_or_no_convergence(test_dose_response, bmd_feasibility_flag, str(chemical_id), end_point, None)
        else:
            # Fit dose response models
            model_predictions = bmdest.analyze_dose_response_data(test_dose_response)
            #print ("model_predictions"+str(model_predictions))
            
            # Select best model
            selected_model_params = bmdest.select_model(model_predictions)
            # Check if unique model is found
            unique_model_flag = selected_model_params['no_unique_model_found_flag']
            if(unique_model_flag == 0):
                # Generate report
                #print(test_dose_response.dose[-1:])
                
                ps.save_results_good_data_unique_model(test_dose_response, bmd_feasibility_flag, model_predictions, selected_model_params, str(chemical_id), end_point)
            else:
                bmd_analysis_flag = selected_model_params['model_select_flag']
                if(bmd_analysis_flag == 1):
                    ps.save_results_poor_data_or_no_convergence(test_dose_response, bmd_feasibility_flag, str(chemical_id), end_point, selected_model_params)
                else:
                    ps.save_results_good_data_nounique_model(test_dose_response, bmd_feasibility_flag, model_predictions, selected_model_params, str(chemical_id), end_point)
                    
    number_of_chemicals_processed += 1
    print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_processed)
    print(print_this)

bmd_feasibility_flag_file_out.close()
end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("BMD calculation is done, it took:"+str(time_took)) 
# (old) for all combinations of 342 chemicals and 18 endpoints, 4 minutes took for bmd_feasibility only
# (old) for all combinations of 342 chemicals and 18 endpoints, 104~165 minutes took for bmd_feasibility and bmd report
# (04/27/2021) for all combinations of 335 chemicals and 18 endpoints, 15.5 hrs took for bmd_feasibility and bmd report in constance

os.chdir(output_folder)
time_filename = os.path.join("report", 'running_time.txt')
f_time = open(time_filename, 'w')
f_time.write(str(time_took))
f_time.close() 

