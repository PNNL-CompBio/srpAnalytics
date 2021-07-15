#!/usr/bin/env python
# coding: utf-8

# In[1]:


#7_PAH_zf_morphology
import numpy as np
import pandas as pd
import os, random, shutil, sys, time
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

#mac
#util_path = "/Users/kimd999/Dropbox/script/python/srpAnalytics/code/latest/util"
#sys.path.insert(0, util_path)

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

# In[2]:



# In[3]:


import BMD_BMDL_estimation as bmdest
import generate_dose_response as gdr
import Plot_Save as ps


# In[4]:


# old - used to work
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/7_PAH/01_11_2021/input/wide/7_PAH_zf_morphology_data_2021JAN11_wide_made_in_2021_01_19_DNC_0.csv'

# mac
#complete_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/7_PAH/morpho/input/wide/7_PAH_zf_morphology_data_2020NOV11_wide_made_in_2021_07_14.csv'

# constance
complete_file_path = '/people/kimd999/tox/7_PAH/input/7_PAH_zf_morphology_data_2020NOV11_wide_made_in_2021_07_14.csv'

dir_of_inputfile = os.path.dirname(complete_file_path)
os.chdir(dir_of_inputfile)

df_morpho = pd.read_csv(complete_file_path, header = 0)
pd.set_option('display.max_columns', None)
print(df_morpho.head())
print(df_morpho.columns)
#display(np.unique(df_morpho.well))


# In[5]:


#np.sum(morphological_data['MO24'] == 1)


# In[6]:


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
    df_morpho['ALL_BUT_MORT'] = df_morpho[['AXIS', 'BRN_', 'CRAN', 'DP24', 'EDEM',                                            'LTRK', 'MUSC', 'NC__', 'SKIN', 'SM24', 'TCHR']].sum(axis=1,skipna=True,min_count=1)


# In[7]:


print(df_morpho.head())
print(df_morpho.tail())


# In[8]:


if (os.path.isdir("output") == True):
    shutil.rmtree("output")
os.mkdir("output")

output_folder = os.path.join(dir_of_inputfile, "output")
os.chdir(output_folder)

if (os.path.isdir("report") == True):
    shutil.rmtree("output")
os.mkdir("report")
    
# df_morpho_filename = os.path.join("report", 'df_morpho.csv')
# df_morpho_file_out = open(df_morpho_filename, "w")
# df_morpho.to_csv(df_morpho_filename, index=False)
# df_morpho_file_out.close()


# In[ ]:


# Goal: calculate BMD
# Specify end_point and chemical of interest
# Perform a check of the existence of "essential" column labels
os.chdir(output_folder)
start_time = time.time()


'''
qc_flag_filename = os.path.join("report", 'qc_flag.csv')
qc_flag_file_out = open(qc_flag_filename, "w")

write_this = "qc_flag\n"
qc_flag_file_out.write(write_this)

erased_morphological_data_end_point_chemical_id_filename = os.path.join("report", 'erased_morphological_data_end_point_chemical_id.csv')
erased_morphological_data_end_point_chemical_id_file = open(erased_morphological_data_end_point_chemical_id_filename, "w")
write_this="chemical_id,plate_id,end_point\n"
erased_morphological_data_end_point_chemical_id_file.write(write_this)
erased_morphological_data_end_point_chemical_id_file.close()


erased_morphological_data_end_point_chemical_id_filename_0p25_erased = erased_morphological_data_end_point_chemical_id_filename[:-4] + '_0p25_erased.csv'
erased_morphological_data_end_point_chemical_id_file_0p25_erased = open(erased_morphological_data_end_point_chemical_id_filename_0p25_erased, "w")
write_this="chemical_id,end_point,dose\n"
erased_morphological_data_end_point_chemical_id_file_0p25_erased.write(write_this)
erased_morphological_data_end_point_chemical_id_file_0p25_erased.close()


erased_morphological_data_end_point_chemical_id_filename_0p25_kept = erased_morphological_data_end_point_chemical_id_filename[:-4] + '_0p25_kept.csv'
erased_morphological_data_end_point_chemical_id_file_0p25_kept = open(erased_morphological_data_end_point_chemical_id_filename_0p25_kept, "w")
write_this="chemical_id,end_point,dose\n"
erased_morphological_data_end_point_chemical_id_file_0p25_kept.write(write_this)
erased_morphological_data_end_point_chemical_id_file_0p25_kept.close()
'''

full_devel = "full"
#full_devel = "devel"

# full -> 17 (without DNC) unlike phase_I_II (18 endpoints), 7_PAH lacks NC24
if (full_devel == "full"):
    end_points = ['ANY24','ANY120','AXIS','ALL_BUT_MORT','BRN_','CRAN','DP24','EDEM','LTRK','MO24','MORT','MUSC','NC__', 'SKIN','SM24','TCHR','TOT_MORT']
else: # full_devel = "devel"
    end_points = ['ANY24','CRAN']

    
if (full_devel == "full"):
    # all chemicals
    chemical_id_from_here = np.unique(df_morpho['chemical.id'])
else: # full_devel = "devel"
    chemical_id_from_here = [3756]

for chemical_id in chemical_id_from_here:
    print("chemical_id:" + str(chemical_id))

    for end_point in end_points:
        os.chdir(output_folder)
        # subset original dataframe for a user-specified chemical and end_point pair
        df_morpho_end_point_chemical_id = df_morpho.loc[df_morpho['chemical.id'] == chemical_id,['chemical.id', 'conc', 'plate.id', 'well', end_point]]
        
        # Binarize end-point hits (Values > 1 are forced to 1)
        end_point_hits = df_morpho_end_point_chemical_id[end_point]
        end_point_hits.loc[end_point_hits > 0] = 1
        
        #print (str(morphological_data_end_point_chemical_id))
   #     morphological_data_end_point_chemical_id.to_csv('morpho.csv', index=False)

#        f_end_point = open('end_point.txt', 'w')
 #       f_end_point.write(str(end_point))
  #      f_end_point.close()
                  
        dose_response = gdr.gen_dose_response(df_morpho_end_point_chemical_id, end_point)
        
        qc_flag = gdr.BMD_feasibility_analysis(dose_response)
        # qc_flag = gdr.BMD_feasibility_analysis_qc_1(dose_response)
        # qc_flag_file_out.write(str(qc_flag)+"\n")
        
        test_dose_response = gdr.reformat_dose_response(dose_response)
        
#        write_this = str(chemical_id) + "," + str(end_point) + "," + str(len(test_dose_response)) + "\n"
 #       print ("write_this:"+str(write_this))
  #      f_out.write(write_this)
    
        #qc_flag_folder = "qc_" + str(qc_flag)
        #if (os.path.isdir(str(qc_flag_folder)) == False):
        #    os.mkdir(str(qc_flag_folder))
        #os.chdir(str(qc_flag_folder))

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
                print(test_dose_response.dose[-1:])
                ps.save_results_good_data_unique_model(test_dose_response, qc_flag, model_predictions, selected_model_params,                                                        str(chemical_id), end_point)
            else:
                bmd_analysis_flag = selected_model_params['model_select_flag']
                if(bmd_analysis_flag == 1):
                    ps.save_results_poor_data_or_no_convergence(test_dose_response, qc_flag, str(chemical_id), end_point,                                                                 selected_model_params)
                else:
                    ps.save_results_good_data_nounique_model(test_dose_response, qc_flag, model_predictions,                                                              selected_model_params, str(chemical_id), end_point)
#test_dose_f_out.close()
#f_out.close()
#qc_flag_file_out.close()
end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("BMD calculation is done, it took:"+str(time_took)) 
# took 21 minutes for all 7 PAHs and all endpoints

os.chdir(output_folder)
      
time_filename = os.path.join("report", 'running_time.txt')
f_time = open(time_filename, 'w')
f_time.write(str(time_took))
f_time.close()


# In[ ]:


#print (morphological_data_end_point_chemical_id)
np.asarray(morphological_data_end_point_chemical_id['plate.id'])


# In[ ]:


#[np.unique(morphological_data_end_point_chemical_id['chemical.id'].values()),1]


# In[ ]:


os.chdir(output_folder)

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

