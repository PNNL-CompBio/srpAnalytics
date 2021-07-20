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

print ("reformat for morphological data starts----------------------")

#complete_file_path = '/Users/kimd999/Dropbox/script/python/srpAnalytics/analysis/paritosh_original_then_edit/to_dockerize/data/7_PAH_zf_morphology_data_2020NOV11_tall.csv'
complete_file_path = args[1]
print ("complete_file_path:" + str(complete_file_path))

full_devel = args[2]




morph_all_data = pd.read_csv(complete_file_path, header = 0)


morph_all_data.head(2400).tail


np.sum(morph_all_data['value'] == 1)


np.sum(morph_all_data['value'] == 0) # -> 82%


np.sum(np.isnan(morph_all_data['value']))


morph_all_data.shape



# Keep only relevant columns
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
morph_all_data_select = morph_all_data.loc[:,columns_to_keep]
morph_all_data_select.head()


start_time = time.time()

reformat_data = pd.DataFrame()
total_number_of_unique_chemicals = 0
total_number_of_chemical_plate_well = 0

if (full_devel == "full"):
    # all chemicals
    chemical_id_from_here = np.unique(morph_all_data['chemical.id'])
else: # full_devel = "deve01l"
    chemical_id_from_here = np.unique([3756])

#for chemical_index in np.unique(morph_all_data['chemical.id']):
for chemical_index in chemical_id_from_here:

    print("chemical_index:" + str(chemical_index))
    total_number_of_unique_chemicals += 1
    morph_data_chemical = morph_all_data_select.loc[morph_all_data['chemical.id'] == chemical_index,:]

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


pd.set_option('display.max_columns', None)
print (reformat_data.head())
print ("reformat_data.shape:" + str(reformat_data.shape)) #(288, 19)

reformat_data_DNC_0 = pd.DataFrame()

reformat_data_DNC_0 = reformat_data.loc[reformat_data['DNC_'] == 0.0]
print ("reformat_data_DNC_0.shape:" + str(reformat_data_DNC_0.shape)) #(287, 19)


output_complete_file_path = complete_file_path[:-4] + "_wide_DNC_0.csv"
print ("output_complete_file_path:" + str(output_complete_file_path))
reformat_data_DNC_0.to_csv(output_complete_file_path,index=False)
#a=b
print ("Reformat of morpho data is done")
