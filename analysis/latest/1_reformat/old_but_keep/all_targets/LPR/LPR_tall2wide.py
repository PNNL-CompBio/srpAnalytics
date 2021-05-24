#!/usr/bin/env python
# coding: utf-8

# ### Reformat LPR behavioral data to have t0-t239
# ### While reformatting, divide data into 240 and 15 timepoints sets respectively

# In[31]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, random, time
from datetime import datetime

import warnings
warnings.filterwarnings('ignore')


# In[32]:


starting_dir = os.getcwd()
print (starting_dir)


# In[33]:


# mac - phase I & II - LPR - before merging - chemical 53 only
complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/latest/after_merging/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_merged.csv'

# mac - phase III - full
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/input/original/behavior/LPR/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'
# -> 240 unique variables and 215 unique chemical IDs


# In[34]:


df_behav = pd.read_csv(complete_input_file_path, header = 0)
df_behav = df_behav.rename({'endpoint': 'timepoint'}, axis=1)
df_behav = df_behav.rename({'variable': 'timepoint'}, axis=1)

df_behav['chemical.id'] = df_behav['chemical.id'].astype(str)
# this recasting is needed for "df_select_1846 = df_select.loc[df_select['chemical.id'] == '1846',:]" later

display(df_behav.head())
display(df_behav.tail())


# In[41]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'timepoint', 'value']
df_select = df_behav.loc[:,columns_to_keep]

display(df_select.head())
display("number of unique chemical.id:" + str(len(np.unique(df_select['chemical.id']))))
display("number of unique timepoints: " + str(len(np.unique(df_select['timepoint']))))


# In[36]:


nan = df_select[df_select['value'].isna()]
display(nan.head())

# [phase III] there is no nan in 'chemical.id', 'conc', 'plate.id', 'well', 'variable'

# don't drop na now for easier proceesing for now


# ### Transpose time points 

# In[44]:


# old using groupby
#'''

start_time = time.time()
       
def reformat(chemical_index, df_select, df_reformatted_240_timepoints, df_reformatted_15_timepoints):
    df_per_chemical = df_select.loc[df_select['chemical.id'] == chemical_index,:]
    #display (df_per_chemical.head())

    # Append chemical_plate_well as a unique identifier
    # takes long time (~1 min)
    df_per_chemical.insert(0, 'chemical_plate_well', df_per_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    for cpw in np.unique(df_per_chemical.chemical_plate_well):
        #print (str(cpw))
        df_per_cpw = df_per_chemical.loc[df_per_chemical.chemical_plate_well == cpw,:]
        df_per_cpw_grouped = df_per_cpw.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in df_per_cpw_grouped:
            concat_this = pd.DataFrame(
                    {
                    'chemical.id': np.unique(df_per_cpw['chemical.id']),
                    'plate.id': np.unique(df_per_cpw['plate.id']),
                    'well': np.unique(df_per_cpw['well']),
                    'chemical_plate_well': np.unique(df_per_cpw['chemical_plate_well']),
                    'conc': np.unique(df_per_cpw['conc'])
                    })
            
            timepoints_15 = False # init
            # rename timepoint columns if this is for 15 endpoints
            for time_point in np.arange(len(np.unique(group.timepoint))):
                if (len(np.unique(group.timepoint)) == 15):
                    timepoints_15 = True
                    time_point = time_point + 3
                timepoint = 't'+ str(time_point)
                concat_this = pd.concat([concat_this, pd.DataFrame({timepoint: df_per_cpw.value[df_per_cpw.timepoint == timepoint].values})],axis = 1)

            if (timepoints_15 == False):
                df_reformatted_240_timepoints = pd.concat([df_reformatted_240_timepoints, concat_this])
            else:
                df_reformatted_15_timepoints = pd.concat([df_reformatted_15_timepoints, concat_this])

    return df_reformatted_240_timepoints, df_reformatted_15_timepoints
########### end of def reformat(chemical_index, behav_select, df_reformatted):


df_reformatted_240_timepoints = pd.DataFrame()
df_reformatted_15_timepoints = pd.DataFrame()

full_devel = "full"
#full_devel = "devel"

chemical_id_from_here = np.unique(df_behav['chemical.id'])

if (full_devel == "devel"):
    randomly_chosen = random.sample(set(chemical_id_from_here), 2)
    chemical_id_from_here = []
    for i in range(len(randomly_chosen)):
        chemical_id_from_here.append(randomly_chosen[i])

#chemical_id_from_here = ['1030', '1119']
chemical_id_from_here = ['53']
# 1119 chemical.id ->  15 timepoints
# 1030 chemical.id -> 240 timepoints

total_number_of_chemicals_to_processed = len(chemical_id_from_here)
number_of_chemicals_processed = 0

f_report = open("unique_chemical_ids_per_run.txt", "w")

for chemical_index in chemical_id_from_here:
    print_this = "\nchemical_index:" + str(chemical_index) + "\n"
    print(print_this)
    f_report.write(print_this)
    
    df_reformatted_240_timepoints, df_reformatted_15_timepoints     = reformat(chemical_index, df_select, df_reformatted_240_timepoints, df_reformatted_15_timepoints)
    
    number_of_chemicals_processed += 1
    print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_processed)
    print(print_this)
    
    write_this = 'number of unique chemical.id:' + str(len(np.unique(df_reformatted_240_timepoints['chemical.id']))) + "\n"
    print(write_this)
    f_report.write(write_this)

    write_this = 'unique chemical.id:' + str(np.unique(df_reformatted_240_timepoints['chemical.id'])) + "\n\n"
    print(write_this)
    f_report.write(write_this)
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Transposing timepoints is done. It took " + str(time_took)) 
# took 75 seconds in pnnl laptop for 1 chemical
# took 5~7 hrs in pnnl laptop for 196 chemicals
#'''

f_report.close()


# In[ ]:


display(df_reformatted_240_timepoints.head())

chemical_id_is_NaN = df_reformatted_240_timepoints['chemical.id'].isnull().values.any()

display(chemical_id_is_NaN)
row_has_NaN = is_NaN.any(axis=1)

display (len(np.unique(df_reformatted_240_timepoints['chemical.id'])))
display     (np.unique(df_reformatted_240_timepoints['chemical.id']))


'''
is_NaN = df_reformatted_240_timepoints.isnull()
row_has_NaN = is_NaN.any(axis=1)
rows_with_NaN = df_reformatted_240_timepoints[row_has_NaN]

display(rows_with_NaN.head())
display(rows_with_NaN.tail())
'''


# In[ ]:


display(df_reformatted_15_timepoints.head())


# In[ ]:


display (len(np.unique(df_reformatted_240_timepoints['chemical.id'])))
display (np.unique(df_reformatted_240_timepoints['chemical.id']))


# In[ ]:


reformatted_data_filename = str(complete_input_file_path[:-4]) + "_wide_t0_t239_" + str(full_devel) + ".csv"
display ("reformatted_data_filename:", reformatted_data_filename)
df_reformatted_240_timepoints.to_csv(reformatted_data_filename, index=False)

reformatted_data_filename = str(complete_input_file_path[:-4]) + "_wide_t3_t17_" + str(full_devel) + ".csv"
display ("reformatted_data_filename:", reformatted_data_filename)
df_reformatted_15_timepoints.to_csv(reformatted_data_filename, index=False)


# In[ ]:




