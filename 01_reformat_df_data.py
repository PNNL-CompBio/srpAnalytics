#!/usr/bin/env python
# coding: utf-8

# In[16]:


# To deal with 7_PAH_zf_morphology_data_2020NOV11.csv

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time

args = sys.argv[0:]


#complete_file_path = '/Users/kimd999/Dropbox/script/python/srpAnalytics/analysis/paritosh_original_then_edit/to_dockerize/data/7_PAH_zf_morphology_data_2020NOV11_tall.csv'
complete_file_path = args[1]
print ("complete_file_path:" + str(complete_file_path))

morph_all_data = pd.read_csv(complete_file_path, header = 0)


# In[18]:


morph_all_data.head(2400).tail


# In[19]:


np.sum(morph_all_data['value'] == 1)


# In[20]:


np.sum(morph_all_data['value'] == 0) # -> 82%


# In[21]:


np.sum(np.isnan(morph_all_data['value']))


# In[22]:


morph_all_data.shape


# In[23]:


# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
morph_all_data_select = morph_all_data.loc[:,columns_to_keep]
morph_all_data_select.head()


# In[14]:


start_time = time.time()

reformat_data = pd.DataFrame()
total_number_of_unique_chemicals = 0
total_number_of_chemical_plate_well = 0
for chemical_index in np.unique(morph_all_data['chemical.id']):
    print("chemical_index:" + str(chemical_index))
    total_number_of_unique_chemicals += 1
    morph_data_chemical = morph_all_data_select.loc[morph_all_data['chemical.id'] == chemical_index,:]

#    if (chemical_index != 3756):
 #       continue

    # Append chemical_plate_well as a unique identifier
    morph_data_chemical.insert(0, 'chemical_plate_well', morph_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))

    for cpw in np.unique(morph_data_chemical.chemical_plate_well):
        total_number_of_chemical_plate_well += 1
        temp_df = morph_data_chemical.loc[morph_data_chemical.chemical_plate_well == cpw,:]
        temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
        for name, group in temp_df_grouped:
            try:
#            if(len(group.endpoint) == 14):
                temp = pd.DataFrame( {
                        'chemical.id': np.unique(temp_df['chemical.id']),
                        'plate.id': np.unique(temp_df['plate.id']),
                        'well': np.unique(temp_df['well']),
                        'chemical_plate_well': np.unique(temp_df['chemical_plate_well']),
                        'conc': np.unique(temp_df['conc']),
                        'AXIS': temp_df.value[temp_df.endpoint == 'AXIS'].values,
                        'BRN_': temp_df.value[temp_df.endpoint == 'BRN_'].values,
                        'CRAN': temp_df.value[temp_df.endpoint == 'CRAN'].values,
                        'DNC_': temp_df.value[temp_df.endpoint == 'DNC_'].values,
                        'DP24': temp_df.value[temp_df.endpoint == 'DP24'].values,
                        'EDEM': temp_df.value[temp_df.endpoint == 'EDEM'].values,
                        'LTRK': temp_df.value[temp_df.endpoint == 'LTRK'].values,
                        'MO24': temp_df.value[temp_df.endpoint == 'MO24'].values,
                        'MORT': temp_df.value[temp_df.endpoint == 'MORT'].values,
                        'MUSC': temp_df.value[temp_df.endpoint == 'MUSC'].values,
                        'NC__': temp_df.value[temp_df.endpoint == 'NC__'].values,
                        'SKIN': temp_df.value[temp_df.endpoint == 'SKIN'].values,
                        'SM24': temp_df.value[temp_df.endpoint == 'SM24'].values,
                        'TCHR': temp_df.value[temp_df.endpoint == 'TCHR'].values,  }  )
                reformat_data = pd.concat([reformat_data, temp])
            except:
                print ("len(group.endpoint) != 14")
                print ("chemical_plate_well:" + str(cpw))
print ("total_number_of_unique_chemicals:" + str(total_number_of_unique_chemicals))
print ("total_number_of_chemical_plate_well:" + str(total_number_of_chemical_plate_well))
end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took))


# In[15]:


pd.set_option('display.max_columns', None)
print (reformat_data.head())
print ("reformat_data.shape:" + str(reformat_data.shape))


# In[33]:


reformat_data_DNC_0 = pd.DataFrame()


# In[34]:


reformat_data_DNC_0 = reformat_data.loc[reformat_data['DNC_'] == 0.0]
print ("reformat_data_DNC_0.shape:" + str(reformat_data_DNC_0.shape))


# In[ ]:


'''
nan_count = 0
zero_count = 0
one_count = 0
for (columnName, columnData) in reformat_data.iteritems():
    if (columnName == "chemical.id") or (columnName == "plate.id") or (columnName == "well") or (columnName == "chemical_plate_well") or (columnName == "conc"):
        continue
#    print('Colunm Name : ', columnName)
#    print('Column Contents : ', columnData.values[i])
#    print (len(columnData.values))
    for i in range(len(columnData.values)):
        if (str(columnData.values[i]) == "nan"):
            nan_count += 1
        elif (str(columnData.values[i]) == "0.0"):
            zero_count += 1
        elif (str(columnData.values[i]) == "1.0"):
            one_count += 1
#        print('Column Contents : ', columnData.values[i])
print ("nan_count:" + str(nan_count))
print ("zero_count:" + str(zero_count))
print ("one_count:" + str(one_count))
'''


# In[35]:


pd.set_option('display.max_columns', None)


# In[39]:

output_complete_file_path = complete_file_path[:-4] + "_wide_DNC_0.csv"

#reformat_data.to_csv('7_PAH_zf_morphology_data_2020NOV11_wide_DNC_0_1.csv',index=False)
#reformat_data_DNC_0.to_csv('7_PAH_zf_morphology_data_2020NOV11_wide_DNC_0.csv',index=False)
reformat_data_DNC_0.to_csv(output_complete_file_path,index=False)


# In[37]:


len(np.unique(reformat_data['chemical.id']))


# ### Reformat behavioral data

# In[ ]:


#complete_file_path = '/Users/pand381/Downloads/Phase 1 and 2 for PNNL database 2020JUNE21/344 zf LPR data phase 1 and 2 - 2020JUNE25.csv'
#behav_all_data = pd.read_csv(complete_file_path, header = 0)


# In[ ]:


#behav_all_data.head()


# In[ ]:


# Keep only relevant columns
#columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
#behav_all_data_select = behav_all_data.loc[:,columns_to_keep]
#behav_all_data_select.head()


# In[ ]:


'''reformat_data = pd.DataFrame()

for chemical_index in np.unique(behav_all_data['chemical.id']):
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print("chemical_index:" + str(chemical_index))
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


# In[ ]:


#reformat_data.head()


# In[ ]:


#reformat_data.shape


# In[ ]:


#reformat_data.to_csv('Phase_I_II_t3_t17_LPR.csv',index=False)


# ### t0_t239 time points

# In[ ]:


'''reformat_data = pd.DataFrame()

#for chemical_index in np.unique([54]):
for chemical_index in np.unique(behav_all_data['chemical.id']):
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print(chemical_index)
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
                max_time = 240
                for time_point in np.arange(max_time):
                    end_point = 't'+ str(time_point)
                    temp = pd.concat([temp, pd.DataFrame({end_point: temp_df.value[temp_df.endpoint == end_point].values})],axis = 1)
                #print(temp.head())
                #print(reformat_data)
                reformat_data = pd.concat([reformat_data, temp])
                #print(reformat_data)'''


# In[ ]:


#display(reformat_data)


# In[ ]:


#reformat_data.to_csv('Phase_I_II_t0_t239_LPR.csv',index=False)


# In[ ]:


#len(np.unique(reformat_data['chemical.id']))


# ## Reformat EPR data

# In[ ]:


#complete_file_path = '/Users/pand381/Downloads/Phase 1 and 2 for PNNL database 2020JUNE21/344 zf EPR data phase 1 and 2 - 2020JUNE25.csv'
#behav_all_data = pd.read_csv(complete_file_path, header = 0)


# In[ ]:


#behav_all_data.head()


# In[ ]:


# Keep only relevant columns
#columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
#behav_all_data_select = behav_all_data.loc[:,columns_to_keep]
#behav_all_data_select.head()


# In[ ]:


'''reformat_data = pd.DataFrame()

#for chemical_index in np.unique([54]):
for chemical_index in np.unique(behav_all_data['chemical.id']):
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]
    print(chemical_index)
    # Append chemical_plate_well as a unique identifier
    behav_data_chemical.insert(0, 'chemical_plate_well', behav_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))

    for cpw in np.unique(behav_data_chemical.chemical_plate_well):
        temp_df = behav_data_chemical.loc[behav_data_chemical.chemical_plate_well == cpw,:]
        #display(temp_df.head())
        temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
        #display(temp_df_grouped.head())
        for name, group in temp_df_grouped:
            if(len(group.endpoint) == 49):
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
                max_time = 49
                for time_point in np.arange(max_time):
                    end_point = 't'+ str(time_point+1)
                    temp = pd.concat([temp, pd.DataFrame({end_point: temp_df.value[temp_df.endpoint == end_point].values})],axis = 1)
                #print(temp.head())
                #print(reformat_data)
                reformat_data = pd.concat([reformat_data, temp])
                #print(reformat_data)'''


# In[ ]:


#display(reformat_data)


# In[ ]:


#reformat_data.to_csv('Phase_I_II_t1_t49_EPR.csv',index=False)
