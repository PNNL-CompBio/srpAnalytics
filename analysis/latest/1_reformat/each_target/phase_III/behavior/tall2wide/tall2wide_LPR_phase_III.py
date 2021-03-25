#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, time

import warnings
warnings.filterwarnings('ignore')

# In[10]:


starting_dir = os.getcwd()
print (starting_dir)

# ### Reformat LPR behavioral data

# In[11]:


#mac
complete_input_file_path= '/Users/kimd999/research/projects/toxicity/per_each_data/phase_III/original/behavior/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'

#constance
#complete_input_file_path= '/people/kimd999/tox/phase_III/Tanguay_Phase_3_zf_LPR_data_PNNL_2021MAR23.csv'

behav_all_data = pd.read_csv(complete_input_file_path, header = 0)
display(len(np.unique(behav_all_data['chemical.id'])))
display(behav_all_data.head())
#display(behav_all_data.tail())

# In[12]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'variable', 'value']
behav_all_data_select = behav_all_data.loc[:,columns_to_keep]
behav_all_data_select.head()

# In[13]:


display("chemical_id_count:"+str(len(np.unique(behav_all_data_select['chemical.id']))))

# In[14]:


#behav_all_data_select = behav_all_data_select.dropna(subset=['chemical.id'])
### no row is dropped for this phase III data

#display("chemical_id_count:"+str(len(np.unique(behav_all_data_select['chemical.id']))))

# ### Transpose time points

# In[ ]:


start_time = time.time()

reformatted_w_240_variables = pd.DataFrame()
reformatted_w_non_240_variables = pd.DataFrame()
    
len_group_variable = []
    
full_devel = "full"
#full_devel = "devel"

if (full_devel == "full"):
    # all chemicals
    chemical_id_from_here = np.unique(behav_all_data['chemical.id'])
else: # full_devel = "devel"
    chemical_id_from_here = np.unique([234])
    
for chemical_index in chemical_id_from_here:
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print("chemical_index:" + str(chemical_index))
    
    print("len(behav_data_chemical):" + str(len(behav_data_chemical)))
    
    
    # Append chemical_plate_well as a unique identifier
    behav_data_chemical.insert(0, 'chemical_plate_well', behav_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))
    
    #output_filename = str(chemical_index) + ".csv"
    #behav_data_chemical.to_csv(output_filename,index=False)
    
    for cpw in np.unique(behav_data_chemical.chemical_plate_well):
        per_cpw = behav_data_chemical.loc[behav_data_chemical.chemical_plate_well == cpw,:]

        per_cpw_grouped = per_cpw.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in per_cpw_grouped:
            
            #display("group.variable", group.variable.head())
            len_group_variable.append(str(len(group.variable)))
            temp = pd.DataFrame(
                    {
                    'chemical.id': np.unique(per_cpw['chemical.id']),
                    'plate.id': np.unique(per_cpw['plate.id']),
                    'well': np.unique(per_cpw['well']),
                    'chemical_plate_well': np.unique(per_cpw['chemical_plate_well']),
                    'conc': np.unique(per_cpw['conc'])
                    })
            #display("before concat:", temp.head())
            
            # Append additonal columns corresponding to time points
            for time_point in np.arange(len(group.variable)):
                end_point = 't'+ str(time_point)
                #print ("\nend_point:"+str(end_point))
                #print ("pd.DataFrame({end_point: per_cpw.value[per_cpw.variable == end_point]}):"\
                #   +str(pd.DataFrame({end_point: per_cpw.value[per_cpw.variable == end_point]})))
                #print ("pd.DataFrame({end_point: per_cpw.value[per_cpw.variable == end_point].values}):"\
                #   +str(pd.DataFrame({end_point: per_cpw.value[per_cpw.variable == end_point].values})))
                temp = pd.concat([temp, pd.DataFrame({end_point: per_cpw.value[per_cpw.variable == end_point].values})],axis = 1)
            #display("after concat:", temp.head())


            if(len(group.variable) == 240): # because we are dealing t0-t239
                reformatted_w_240_variables = pd.concat([reformatted_w_240_variables, temp])
            else: # len(group.variable) != 240
                reformatted_w_non_240_variables = pd.concat([reformatted_w_non_240_variables, temp])
             

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took)) 
# took 2 minutes for 2 chemicals

# In[ ]:


# ds = pd.Series({"Column": len_group_variable})
# plt.figure(figsize=(8,4))
# sns.countplot(x="Column", data=ds)
# plt.show()

# In[ ]:


display(reformatted_w_240_endpoints)
display(reformatted_w_non_240_endpoints)

# In[ ]:


reformatted_data_filename = str(complete_input_file_path[:-9]) + "_wide_t0_t239_" + str(full_devel) + "_w_240_endpoints.csv"
reformatted_w_240_endpoints.to_csv(reformatted_data_filename,index=False)

reformatted_data_filename = str(complete_input_file_path[:-9]) + "_wide_t0_t239_" + str(full_devel) + "_w_non_240_endpoints.csv"
reformatted_w_non_240_endpoints.to_csv(reformatted_data_filename,index=False)

# In[ ]:


#display(len(np.unique(reformatted_w_240_endpoints['chemical.id'])))
display(len(np.unique(reformatted_w_non_240_endpoints['chemical.id'])))

# In[ ]:


reformatted_w_non_240_endpoints_155 = reformatted_w_non_240_endpoints[reformatted_w_non_240_endpoints['chemical.id']==155]

display(reformatted_w_non_240_endpoints_155.head())
display(reformatted_w_non_240_endpoints_155.tail())

output_filename = str(complete_input_file_path[:-4]) + "_155_chemical_id.csv"
print ("output_filename:"+str(output_filename))
reformatted_w_non_240_endpoints_155.to_csv(output_filename, index=False)

# In[ ]:


reformatted_w_non_240_endpoints_163 = reformatted_w_non_240_endpoints[reformatted_w_non_240_endpoints['chemical.id']==163]

display(reformatted_w_non_240_endpoints_163.head())
display(reformatted_w_non_240_endpoints_163.tail())

output_filename = str(complete_input_file_path[:-4]) + "_163_chemical_id.csv"
print ("output_filename:"+str(output_filename))
reformatted_w_non_240_endpoints_163.to_csv(output_filename, index=False)

# In[ ]:



