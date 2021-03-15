#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time

# ### Reformat LPR behavioral data

args = sys.argv[0:]

print ("reformat_LPR starts----------------------")

starting_dir = os.getcwd()

#complete_input_file_path = '/Users/kimd999/research/projects/toxicity/per_each_data/7_PAH/01_11_2021/input/tall/7_PAH_zf_LPR_data_2021JAN11.csv'
complete_input_file_path = args[1]
behav_all_data = pd.read_csv(complete_input_file_path, header = 0)

full_devel = args[2]

behav_all_data.head()

behav_all_data = behav_all_data.dropna()



#behav_all_data_certain_TX = behav_all_data.loc[behav_all_data['bottle.id'] == "TX002271",['chemical.id', 'bottle.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']]
#behav_all_data_certain_TX.head()


# Keep only relevant columns
behav_all_data = behav_all_data.rename(columns = {"variable":"endpoint"})
columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
behav_all_data_select = behav_all_data.loc[:,columns_to_keep]
behav_all_data_select.head()


print(np.unique(behav_all_data_select['endpoint']))
print(len(np.unique(behav_all_data_select['endpoint'])))


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


start_time = time.time()

reformatted_data = pd.DataFrame()
max_time = 240

#full_devel = "full"
#full_devel = "devel"

if (full_devel == "full"):
    # all chemicals
    chemical_id_from_here = np.unique(behav_all_data['chemical.id'])
else: # full_devel = "devel"
    chemical_id_from_here = np.unique([3756])

for chemical_index in chemical_id_from_here:
    print("chemical_index:" + str(chemical_index))
    behav_data_chemical = behav_all_data_select.loc[behav_all_data['chemical.id'] == chemical_index,:]

    # Append chemical_plate_well as a unique identifier
    behav_data_chemical.insert(0, 'chemical_plate_well', behav_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))

    for cpw in np.unique(behav_data_chemical.chemical_plate_well):
        temp_df = behav_data_chemical.loc[behav_data_chemical.chemical_plate_well == cpw,:]
        #print("temp_df.head():\n" + str(temp_df.head()))
        temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
        #print("temp_df_grouped.head():\n" + str(temp_df_grouped.head()))
        for name, group in temp_df_grouped:
            #print ("len(group.endpoint):" + str(len(group.endpoint)))
            if(len(group.endpoint) == max_time):
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
                #print(reformatted_data)
                reformatted_data = pd.concat([reformatted_data, temp])
                #print(reformatted_data)

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds"
print ("Done, it took:"+str(time_took))
# for 2 chemicals -> took 2 minutes
# for all 7 chemicals -> took 6 minutes

print("reformatted_data:\n" + str(reformatted_data))



print ("starting_dir: " + str(starting_dir)) # /srpAnalytics

cwd = os.getcwd()
print ("cwd: " + str(cwd)) # /srpAnalytics

output_complete_file_path = complete_input_file_path[:-4] + "_wide_t0_t239_" + str(full_devel) + ".csv"
print ("output_complete_file_path after reformat:" + str(output_complete_file_path))
reformatted_data.to_csv(output_complete_file_path,index=False)

print ("output_complete_file_path existence in same py code:" + str(os.path.isfile(output_complete_file_path)))
