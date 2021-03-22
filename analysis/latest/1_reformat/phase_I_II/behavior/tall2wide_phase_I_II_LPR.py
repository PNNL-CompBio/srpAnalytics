#!/usr/bin/env python
# coding: utf-8

# In[1]:


##### overall explanation
# 1. reformat LPR data to have t0-t239
# 2. reformat EPR data to have t1-t49

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, time

import warnings
warnings.filterwarnings('ignore')

# ### Reformat LPR behavioral data

# In[17]:


complete_input_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/input/behavioral/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX.csv'
behav_all_data = pd.read_csv(complete_input_file_path, header = 0)
#behav_all_data.head()

# In[33]:


#len(behav_all_data)

# In[4]:


TX_bottle_id_behav_all_data = behav_all_data[behav_all_data['bottle.id'].str.contains("^TX")]
TX_bottle_id_behav_all_data.head()

# In[29]:


len(TX_bottle_id_behav_all_data)

# In[5]:


##### (start) replace plate.id in behavioral data with plate.id in morphological data (common data will be bottle.id)

complete_morph_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/input/morphology/after_Lisa_cleanup/tall/zf_morphology_data_335_chemicals_2020DEC16.csv'
morph_all_data = pd.read_csv(complete_morph_file_path, header = 0)
#morph_all_data.head()

# In[32]:


#len(morph_all_data)

# In[6]:


TX_bottle_id_morph_all_data = morph_all_data[morph_all_data['bottle.id'].str.contains("^TX")]
TX_bottle_id_morph_all_data.head()

# In[24]:


len(TX_bottle_id_morph_all_data)

# In[15]:


#new_behav_all_data = TX_bottle_id_behav_all_data.loc[TX_bottle_id_morph_all_data['bottle.id'],['chemical.id', 'bottle.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']]

# In[ ]:


##### (end) replace plate.id in behavioral data with plate.id in morphological data (common data will be bottle.id)

# In[11]:


#new_behav_all_data = behav_all_data.loc[behav_all_data['bottle.id'] == "TX002271",['chemical.id', 'bottle.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']]

# In[12]:


#new_behav_all_data.head()

# In[9]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
behav_all_data_select = behav_all_data.loc[:,columns_to_keep]
behav_all_data_select.head()

# In[10]:


display(np.unique(behav_all_data['endpoint']))

# In[11]:


display(np.unique(behav_all_data_select['endpoint']))

# In[16]:


'''
reformat_data = pd.DataFrame()

for chemical_index in np.unique(behav_all_data['chemical.id']):
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print(chemical_index)
    # Append chemical_plate_well as a unique identifier
    behav_data_chemical.insert(0, 'chemical_plate_well', behav_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    for cpw in np.unique(behav_data_chemical.chemical_plate_well):
        temp_df = behav_data_chemical.loc[behav_data_chemical.chemical_plate_well == cpw,:]
        temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in temp_df_grouped:
            if(len(group.endpoint) == 15):
                temp = pd.DataFrame(
                        {
                        'chemical.id': np.unique(temp_df['chemical.id']),
                        'plate.id': np.unique(temp_df['plate.id']),
                        'well': np.unique(temp_df['well']),
                        'chemical_plate_well': np.unique(temp_df['chemical_plate_well']),
                        'conc': np.unique(temp_df['conc']),
                        't3': temp_df.value[temp_df.endpoint == 't3'].values,
                        't4': temp_df.value[temp_df.endpoint == 't4'].values,
                        't5': temp_df.value[temp_df.endpoint == 't5'].values,
                        't6': temp_df.value[temp_df.endpoint == 't6'].values,
                        't7': temp_df.value[temp_df.endpoint == 't7'].values,
                        't8': temp_df.value[temp_df.endpoint == 't8'].values,
                        't9': temp_df.value[temp_df.endpoint == 't9'].values,
                        't10': temp_df.value[temp_df.endpoint == 't10'].values,
                        't11': temp_df.value[temp_df.endpoint == 't11'].values,
                        't12': temp_df.value[temp_df.endpoint == 't12'].values,
                        't13': temp_df.value[temp_df.endpoint == 't13'].values,
                        't14': temp_df.value[temp_df.endpoint == 't14'].values,
                        't15': temp_df.value[temp_df.endpoint == 't15'].values,
                        't16': temp_df.value[temp_df.endpoint == 't16'].values,
                        't17': temp_df.value[temp_df.endpoint == 't17'].values,
                        }
                            )
                reformat_data = pd.concat([reformat_data, temp])'''

# In[17]:


#reformat_data.head()

# In[18]:


#reformat_data.shape

# In[19]:


#reformat_data.to_csv('Phase_I_II_t3_t17_LPR.csv',index=False)

# In[ ]:


###### replace plate id with what is listed in the morphology data to the LPR.
###### for your example, TP967-E9-P1, will go in place of plate.id “1” in the LPR file.

# read morpho data and match/join


# ### t0_t239 time points

# In[ ]:


start_time = time.time()

reformat_data = pd.DataFrame()
max_time = 240
    
full_devel = "full"
#full_devel = "devel"

if (full_devel == "full"):
    # all chemicals
    chemical_id_from_here = np.unique(behav_all_data['chemical.id'])
else: # full_devel = "devel"
    chemical_id_from_here = np.unique([53,54])

    
for chemical_index in chemical_id_from_here:
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print("chemical_index:" + str(chemical_index))
    # Append chemical_plate_well as a unique identifier
    behav_data_chemical.insert(0, 'chemical_plate_well', behav_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    for cpw in np.unique(behav_data_chemical.chemical_plate_well):
        temp_df = behav_data_chemical.loc[behav_data_chemical.chemical_plate_well == cpw,:]
        #display(temp_df.head())
        temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
        #display(temp_df_grouped.head())
        for name, group in temp_df_grouped:
            if(len(group.endpoint) == 240):
                temp = pd.DataFrame(
                        {
                        'chemical.id': np.unique(temp_df['chemical.id']),
                        'plate.id': np.unique(temp_df['plate.id']),
                        'well': np.unique(temp_df['well']),
                        'chemical_plate_well': np.unique(temp_df['chemical_plate_well']),
                        'conc': np.unique(temp_df['conc'])
                        })
                #print(temp.head())
                # Append additonal columns corresponding to time points
                for time_point in np.arange(max_time):
                    end_point = 't'+ str(time_point)
                    temp = pd.concat([temp, pd.DataFrame({end_point: temp_df.value[temp_df.endpoint == end_point].values})],axis = 1)
                #print(temp.head())
                #print(reformat_data)
                reformat_data = pd.concat([reformat_data, temp])
                #print(reformat_data)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took)) 
# took 2 minutes for 2 chemicals

# In[13]:


display(reformat_data)

# In[18]:


reformatted_data_filename = str(complete_input_file_path[:-4]) + "_wide_t0_t239_" + str(full_devel) + ".csv"
reformatted_data.to_csv(reformatted_data_filename,index=False)

# In[19]:


len(np.unique(reformat_data['chemical.id']))
