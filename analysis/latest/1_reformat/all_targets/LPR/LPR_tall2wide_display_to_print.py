#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, time

import warnings
warnings.filterwarnings('ignore')

# In[2]:


starting_dir = os.getcwd()
print (starting_dir)

# ### Reformat LPR behavioral data

# In[3]:


# mac - phase I & II - full
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall.csv'

# mac - phase I & II - 240 endpoints
complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/tall/bifurcated/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_w_240_endpoints.csv'

# mac - phase I & II - 15 endpoints
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_I_II/input/LPR/tall/bifurcated/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_w_15_endpoints.csv'

# mac - phase III - full
#complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/input/original/behavior/LPR/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'


# constance - phase I & II - full
#complete_input_file_path= '/qfs/people/kimd999/tox/phase_I_II/LPR/input/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_full_w_15_endpoints.csv'
#complete_input_file_path= '/qfs/people/kimd999/tox/phase_I_II/LPR/input/tall/344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_full_w_240_endpoints.csv'

# constance - phase III - full
#complete_input_file_path= '/people/kimd999/tox/phase_III/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'

# In[4]:


behav_all_data = pd.read_csv(complete_input_file_path, header = 0)
behav_all_data = behav_all_data.rename({'endpoint': 'variable'}, axis=1)
print(behav_all_data.head())

# In[5]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'variable', 'value']
behav_select = behav_all_data.loc[:,columns_to_keep]

print("number of unique variable:" + str(len(np.unique(behav_select['variable']))))
print("number of chemical_id:" + str(len(np.unique(behav_select['chemical.id']))))
# phase I && II
# 240 variables -> 196 chemicals
# 15  variables -> 148 chemicals

# phase III
# 240 variables -> 215 chemicals

# In[6]:


nan = behav_select[behav_select['value'].isna()]
print(nan.head())

# [phase III] there is no nan in 'chemical.id', 'conc', 'plate.id', 'well', 'variable'

# In[7]:


print("before dropna, len(behav_select):"+str(len(behav_select)))
behav_select = behav_select.dropna(how='any')
print("after dropna,  len(behav_select):"+str(len(behav_select)))

print(behav_select.head())

# phase I & II -> dropped some
# phase III    -> dropped many

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

# ### Transpose time points 

# In[8]:


start_time = time.time()
       
def reformat(chemical_index, behav_select, reformatted):
    behav_per_chemical = behav_select.loc[behav_select['chemical.id'] == chemical_index,:]
    #print (behav_per_chemical)

    # Append chemical_plate_well as a unique identifier
    # takes long time (1 min?)
    behav_per_chemical.insert(0, 'chemical_plate_well', behav_per_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    for cpw in np.unique(behav_per_chemical.chemical_plate_well):
        #print (str(cpw))
        per_cpw = behav_per_chemical.loc[behav_per_chemical.chemical_plate_well == cpw,:]
        per_cpw_grouped = per_cpw.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in per_cpw_grouped:
            final = pd.DataFrame(
                    {
                    'chemical.id': np.unique(per_cpw['chemical.id']),
                    'plate.id': np.unique(per_cpw['plate.id']),
                    'well': np.unique(per_cpw['well']),
                    'chemical_plate_well': np.unique(per_cpw['chemical_plate_well']),
                    'conc': np.unique(per_cpw['conc'])
                    })
            # Append additonal columns corresponding to time points
            for time_point in np.arange(len(np.unique(group.variable))):
                if (len(np.unique(group.variable)) == 15):
                    time_point = time_point + 3
                variable = 't'+ str(time_point)
                #print ("\nvariable:"+str(variable))
                final = pd.concat([final, pd.DataFrame({variable: per_cpw.value[per_cpw.variable == variable].values})],axis = 1)
            reformatted = pd.concat([reformatted, final])
    return reformatted
########### end of def reformat(chemical_index, behav_select, reformatted):

reformatted = pd.DataFrame()

#full_devel = "full"
full_devel = "devel"

if (full_devel == "full"):
    chemical_id_from_here = np.unique(behav_all_data['chemical.id'])
else: # full_devel = "devel"
    #chemical_id_from_here = [53,54]
    chemical_id_from_here = [414]

for chemical_index in chemical_id_from_here:
    print("chemical_index:" + str(chemical_index))
    reformatted = reformat(chemical_index, behav_select, reformatted)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took)) 
# took 2 minutes for 2 chemicals

# In[ ]:


# ds = pd.Series({"Column": len_group_variable})
# plt.figure(figsize=(8,4))
# sns.countplot(x="Column", data=ds)
# plt.show()

# In[51]:


print(reformatted.head())

# In[52]:


reformatted_data_filename = str(complete_input_file_path[:-9]) + "_wide_t0_t239_" + str(full_devel) + "_414.csv"
print ("reformatted_data_filename:", reformatted_data_filename)
reformatted.to_csv(reformatted_data_filename,index=False)

# 15 variables
#reformatted_data_filename = str(complete_input_file_path[:-9]) + "_wide_t0_t239_" + str(full_devel) + "_w_non_240_endpoints.csv"
#reformatted_w_non_240_endpoints.to_csv(reformatted_data_filename,index=False)

# In[40]:


print(len(np.unique(reformatted['chemical.id'])))
#display(len(np.unique(reformatted_w_non_240_endpoints['chemical.id'])))

# In[ ]:



