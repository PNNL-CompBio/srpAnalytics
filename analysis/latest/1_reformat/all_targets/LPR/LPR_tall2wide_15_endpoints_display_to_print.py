#!/usr/bin/env python
# coding: utf-8

# ### Reformat LPR behavioral data to have t0-t239 or t0-t14

# In[9]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, random, time

import warnings
warnings.filterwarnings('ignore')


# In[10]:


starting_dir = os.getcwd()
print (starting_dir)


# In[13]:


# mac - phase I & II - LPR - 240 endpoints
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/after_Lisa_fix/full/240_endpoints/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_full_w_240_endpoints.csv'
# -> 196 unique chemical IDs

# mac - phase I & II - LPR - 15 endpoints
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/after_Lisa_fix/full/15_endpoints/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_full_w_15_endpoints.csv'
# -> 148 unique chemical IDs

# mac - phase III - full
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/input/original/behavior/LPR/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'
# -> 240 unique variables and 215 unique chemical IDs

# constance - phase I & II - LPR - 15 endpoints
complete_input_file_path= '/people/kimd999/tox/phase_I_II/LPR/input/tall/after_Lisa_fix/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_full_w_15_endpoints.csv'

# constance - phase III - full
#complete_input_file_path= '/people/kimd999/tox/phase_III/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'


# In[14]:


df_behav = pd.read_csv(complete_input_file_path, header = 0)
df_behav = df_behav.rename({'endpoint': 'variable'}, axis=1)

print(df_behav.head())


# In[17]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'variable', 'value']
df_select = df_behav.loc[:,columns_to_keep]

df_select['chemical.id'] = df_select['chemical.id'].astype(str)
# this recasting is needed for "df_select_1846 = df_select.loc[df_select['chemical.id'] == '1846',:]" later

print("number of unique variables:" + str(len(np.unique(df_select['variable']))))
print("number of chemical IDs:" + str(len(np.unique(df_select['chemical.id']))))


# In[21]:


# recast some types

start_time = time.time()
       
def fix(chemical_index, df_select, df_fixed):

    df_per_chemical = df_select.loc[df_select['chemical.id'] == str(chemical_index),:]
    #print (df_per_chemical)

    try: # for chemical.id = 1846 
        df_per_chemical['plate.id'] = df_per_chemical['plate.id'].astype(float)
        df_per_chemical['plate.id'] = df_per_chemical['plate.id'].astype(int)    
    except: # for plate.id like TPP166-811...
        df_per_chemical['plate.id'] = df_per_chemical['plate.id'].astype(str)

    df_fixed = pd.concat([df_fixed, df_per_chemical])
    return df_fixed
########### end of def fix(chemical_index, df_select, df_fixed):

df_fixed = pd.DataFrame()

full_devel = "full"
#full_devel = "devel"

chemical_id_from_here = np.unique(df_behav['chemical.id'])

if (full_devel == "devel"):
    randomly_chosen = random.sample(set(chemical_id_from_here), 1)
    chemical_id_from_here = []
    for i in range(len(randomly_chosen)):
        chemical_id_from_here.append(randomly_chosen[i])

#chemical_id_from_here = ['471']

total_number_of_chemicals_to_processed = len(chemical_id_from_here)
number_of_chemicals_processed = 0

for chemical_index in chemical_id_from_here:
    print("chemical_index:" + str(chemical_index))
    df_fixed = fix(chemical_index, df_select, df_fixed)
    
    #number_of_chemicals_processed += 1
    #print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_processed)
    #print(print_this)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("recast some types is done. It took :"+str(time_took)) 
# took 12 seconds for minutes for 148 chemicals


# In[22]:


print(df_fixed.head())


# In[ ]:


### Transpose time points 

start_time = time.time()
       
def reformat(chemical_index, df_fixed, reformatted):
    df_per_chemical = df_fixed.loc[df_fixed['chemical.id'] == str(chemical_index),:]
    #print (df_per_chemical)

    # Append chemical_plate_well as a unique identifier
    # takes long time (1 min?)
    df_per_chemical.insert(0, 'chemical_plate_well', df_per_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    for cpw in np.unique(df_per_chemical.chemical_plate_well):
        #print (str(cpw))
        df_per_cpw = df_per_chemical.loc[df_per_chemical.chemical_plate_well == cpw,:]
        df_per_cpw_grouped = df_per_cpw.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in df_per_cpw_grouped:
            #print ("name:"+str(name))
            #print ("group:"+str(group))
            concat_this = pd.DataFrame(
                    {
                    'chemical.id': np.unique(df_per_cpw['chemical.id']),
                    'plate.id': np.unique(df_per_cpw['plate.id']),
                    'well': np.unique(df_per_cpw['well']),
                    'chemical_plate_well': np.unique(df_per_cpw['chemical_plate_well']),
                    'conc': np.unique(df_per_cpw['conc'])
                    })
            
            # rename endpoint columns if this is for 15 endpoints
            for time_point in np.arange(len(np.unique(group.variable))):
                #print ("np.unique(group.variable):"+str(np.unique(group.variable)))
                if (len(np.unique(group.variable)) == 15):
                    time_point = time_point + 3
                variable = 't'+ str(time_point)
                #print ("\nvariable:"+str(variable))
                concat_this = pd.concat([concat_this, pd.DataFrame({variable: df_per_cpw.value[df_per_cpw.variable == variable].values})],axis = 1)
            reformatted = pd.concat([reformatted, concat_this])
    return reformatted
########### end of def reformat(chemical_index, df_fixed, reformatted):


reformatted = pd.DataFrame()

full_devel = "full"
#full_devel = "devel"

chemical_id_from_here = np.unique(df_behav['chemical.id'])

if (full_devel == "devel"):
    randomly_chosen = random.sample(set(chemical_id_from_here), 1)
    chemical_id_from_here = []
    for i in range(len(randomly_chosen)):
        chemical_id_from_here.append(randomly_chosen[i])

#chemical_id_from_here = ['1846']
#chemical_id_from_here = ['471']

total_number_of_chemicals_to_processed = len(chemical_id_from_here)
number_of_chemicals_processed = 0

for chemical_index in chemical_id_from_here:
    print("chemical_index:" + str(chemical_index))
    reformatted = reformat(chemical_index, df_fixed, reformatted)
    
    number_of_chemicals_processed += 1
    print_this = str(number_of_chemicals_processed) + " chemicals processed out of " + str(total_number_of_chemicals_to_processed)
    print(print_this)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("transposing time points is done. It took :"+str(time_took)) 
# took 5 minutes for 148 chemicals


# In[28]:


print(reformatted.head())


# In[12]:


'''
nan = df_select[df_select['value'].isna()]
print(nan.head())

# [phase III] there is no nan in 'chemical.id', 'conc', 'plate.id', 'well', 'variable'

print("before dropna, len(behav_select):"+str(len(df_select)))
df_select = df_select.dropna(how='any')
print("after dropna,  len(behav_select):"+str(len(df_select)))

print(df_select.head())

# phase I & II -> dropped some
# phase III    -> dropped many
'''


# In[ ]:


''' # basic check of variable #
full_devel = "full"
#full_devel = "devel"

if (full_devel == "full"):
    chemical_id_from_here = np.unique(behav_select['chemical.id'])
else: # full_devel = "devel"
    chemical_id_from_here = np.unique([234])
    
for chemical_index in chemical_id_from_here:
    behav_per_chemical = behav_select.loc[behav_select['chemical.id'] == chemical_index,:]
    print("chemical_index:" + str(chemical_index))

    variables = np.unique(behav_per_chemical['variable'])
    print("variables:" + str(variables))
#    variable_splited = variable.split("t")
    
 #   print (min(variable_splited[1]))
    var_len = len(np.unique(behav_per_chemical['variable']))
    if (var_len != 15):
        display("number of variable:" + str(len(np.unique(behav_per_chemical['variable']))))
'''


# In[69]:


'''df_fixed_1846 = df_fixed.loc[df_fixed['chemical.id'] == '1846',:]
print(df_fixed_1846.head())
#df_select_1846.to_csv("df_select_1846.csv",index=False)
#display(np.unique(df_select_1846['plate.id']))
print(df_fixed_1846['plate.id'])'''


# In[71]:


'''df_select_471 = df_fixed.loc[df_fixed['chemical.id'] == '471',:]
print(df_select_471.head())
#df_select_1846.to_csv("df_select_1846.csv",index=False)
print(np.unique(df_select_471['plate.id']))
'''


# In[79]:


reformatted_data_filename = str(complete_input_file_path[:-9]) + "_wide_t3_t17_" + str(full_devel) + ".csv"
print ("reformatted_data_filename:", reformatted_data_filename)
reformatted.to_csv(reformatted_data_filename,index=False)


# In[80]:


print(len(np.unique(reformatted['chemical.id'])))
#display(len(np.unique(reformatted_w_non_240_endpoints['chemical.id'])))


# In[ ]:




