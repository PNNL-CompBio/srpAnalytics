#!/usr/bin/env python
# coding: utf-8

# In[15]:


from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, random, shutil, sys, time
from scipy import stats
import seaborn as sns

#mac
#util_path = "/Users/kimd999/Dropbox/script/python/srpAnalytics/code/latest/util"

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

import generate_dose_response as gdr
import BMD_BMDL_estimation as bmdest
import Plot_Save as ps

import warnings
warnings.filterwarnings('ignore')


# In[2]:


starting_dir = os.getcwd()
print (starting_dir)


# In[3]:


#mac
#wide_file_path = '/Users/kimd999/research/projects/Katrina/per_each_data/extracts/input/latest/wide/ZF_120_SRP_Sample_Extracts_with_Dilution_Factors_2021MAY21_wide_DNC_0.csv'

# constance
wide_file_path = '/people/kimd999/tox/extracts/input/ZF_120_SRP_Sample_Extracts_with_Dilution_Factors_2021MAY21_wide_DNC_0.csv'

df = pd.read_csv(wide_file_path, header = 0)
pd.set_option('display.max_columns', None)
print(df.head())
print(df.columns)
print(np.unique(df.well))


# In[4]:


#np.sum(df['MO24'] == 1)


# In[5]:


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
    df['ANY24'] = df[['MO24','DP24','SM24','NC24']].sum(axis=1,skipna=True,min_count=1)
    df['ANY120'] = df[['MORT', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC',                                                        'TRUN', 'SWIM', 'NC__', 'TR__', 'ANY24']].sum(axis=1,skipna=True,min_count=1)
    df['TOT_MORT'] = df[['MO24','MORT']].sum(axis=1,skipna=True,min_count=1)
    df['ALL_BUT_MORT'] = df[['DP24','SM24','NC24', 'YSE_', 'AXIS', 'EYE_', 'SNOU', 'JAW_', 'OTIC',                                                        'PE__', 'BRAI', 'SOMI', 'PFIN', 'CFIN', 'PIG_', 'CIRC','TRUN', 'SWIM', 'NC__',                                                        'TR__']].sum(axis=1,skipna=True,min_count=1)
    df['BRN_'] = df[['BRAI','OTIC','PFIN']].sum(axis=1,skipna=True,min_count=1)
    df['CRAN'] = df[['EYE_', 'SNOU', 'JAW_']].sum(axis=1,skipna=True,min_count=1)
    df['EDEM'] = df[['YSE_','PE__']].sum(axis=1,skipna=True,min_count=1)
    df['LTRK'] = df[['TRUN','CFIN']].sum(axis=1,skipna=True,min_count=1)
    df['MUSC'] = df[['CIRC','SWIM','SOMI']].sum(axis=1,skipna=True,min_count=1)
    df['SKIN'] = df[['PIG_']]
    df['TCHR'] = df[['TR__']]


# In[6]:


print(df.head())


# In[7]:


print(df.columns)
print(len(df.columns))


# In[8]:


os.chdir(starting_dir)

if (os.path.isdir("output") == False):
    os.mkdir("output")

output_folder = os.path.join(starting_dir, "output")
os.chdir(output_folder)

if (os.path.isdir("report") == False):
    os.mkdir("report")
    

# <begin> just checking purpose
df_filename = os.path.join("report", 'df_after_aggregating_endpoints.csv')
df_file_out = open(df_filename, "w")
df.to_csv(df_filename, index=False)
df_file_out.close()
# <end> just checking purpose


# In[9]:


print(df.shape)
df_final = df.dropna(how='any')


# In[10]:


print(df_final.shape)
print(df_final.head())


# In[ ]:


# BMD analysis

start_time = time.time()

#full_devel = "full"
full_devel = "devel"

if (full_devel == "full"):
    chemical_id_from_here = np.unique(df['chemical.id'])
    
    # deal 18 end_points
    end_points = ['ANY24','ANY120','AXIS','ALL_BUT_MORT','BRN_','CRAN','DP24','EDEM','LTRK','MO24','MORT','MUSC','NC__','NC24', 'SKIN','SM24','TCHR','TOT_MORT']
else:
    chemical_id_from_here = [101]
    end_points = ['ANY24']

    
total_number_of_chemicals_to_process = len(chemical_id_from_here)
number_of_chemicals_processed = 0

for chemical_id in chemical_id_from_here:
    print(f"\nchemical_id:{chemical_id}")
        
    for end_point in end_points:
        #print(f"end_point:{end_point}")
        os.chdir(output_folder)
    
        # subset original dataframe for a user-specified chemical and end_point pair
        df_end_point_chemical_id = df.loc[df['chemical.id'] == chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]
        #print (f"df_end_point_chemical_id:{df_end_point_chemical_id}")
        
        # Binarize end-point hits (Values > 1 are forced to 1)
        end_point_hits = df_end_point_chemical_id[end_point]
        end_point_hits.loc[end_point_hits > 0] = 1
        
        dose_response = gdr.gen_dose_response(df_end_point_chemical_id, end_point)
        #print (f"dose_response:{dose_response}")
        
        qc_flag = gdr.BMD_feasibility_analysis(dose_response)
        
        test_dose_response = gdr.reformat_dose_response(dose_response)
        #print (f"test_dose_response:{test_dose_response}")
        
    
#        qc_flag_folder = "qc_" + str(qc_flag)
 #       if (os.path.isdir(str(qc_flag_folder)) == False):
  #          os.mkdir(str(qc_flag_folder))
   #     os.chdir(str(qc_flag_folder))

        
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
                #print(test_dose_response.dose[-1:])
                ps.save_results_good_data_unique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), end_point)
            else:
                bmd_analysis_flag = selected_model_params['model_select_flag']
                if(bmd_analysis_flag == 1):
                    ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point, selected_model_params)
                else:
                    ps.save_results_good_data_nounique_model(test_dose_response, qc_flag, model_predictions, selected_model_params, str(chemical_id), end_point)

    number_of_chemicals_processed += 1
    print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_process)
    print(print_this)
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    
end_time = time.time()

time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took)) 
# for all combinations of 42 chemicals and 18 endpoints, <212 minutes took for qc and bmd

f_time = open('running_time.txt', 'w')
f_time.write(str(time_took))
f_time.close()


# In[ ]:


a=b


# In[ ]:


#print (df_end_point_chemical_id)
np.asarray(df_end_point_chemical_id['plate.id'])


# In[ ]:


#[np.unique(df_end_point_chemical_id['chemical.id'].values()),1]


# In[ ]:


os.chdir(starting_dir)

#qc_flag_filename="/Users/kimd999/research/projects/toxicity/result/old_Phase_I_II/newest_criteria_no_avg/report/qc_flag.csv"
qc_flag_filename = os.path.join("report", 'qc_flag.csv')
print ("qc_flag_filename:"+str(qc_flag_filename))
qc_flag_data = pd.read_csv(qc_flag_filename, index_col=None)
#display(qc_flag_data.head())
ds = pd.Series({"Column": qc_flag_data["qc_flag"]})
plt.figure(figsize=(8,4))
sns.countplot(x="Column", data=ds)
plt.show()


# In[ ]:


'''os.chdir(starting_dir)

sns.set_theme(style="whitegrid")
print ("array_filename:"+str(array_filename))
array_report_data = pd.read_csv(array_filename, index_col=None)
print(array_report_data.head())
#ax = sns.barplot(x="end_point", y="len_test_dose_response", data=array_report_data)

ds = pd.Series({"Column": array_report_data["len_test_dose_response"]})
plt.figure(figsize=(8,4))
plt.xlabel("leng")
sns.countplot(x="Column", data=ds)
plt.show()
print ("done")'''


# In[ ]:


test_dose_response.dose


# In[ ]:


#test_dose_response.dose.iloc[0]+test_dose_response.dose.iloc[1]


# In[ ]:


dose_response['num_affect']/dose_response['num_embryos']

