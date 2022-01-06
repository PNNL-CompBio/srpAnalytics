#!/usr/bin/env python
# coding: utf-8

# In[16]:


# To deal with 7_PAH_zf_morphology_data_2020NOV11.csv

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, random, sys, time



def format(complete_file_path, full_devel, output_complete_file_path, chem_ind=None):

    df_morph = pd.read_csv(complete_file_path, header = 0)

    # Keep only relevant columns
    columns_to_keep = ['chemical.id', 'conc', 'plate.id', 'well', 'endpoint', 'value']
    df_morph_select = df_morph.loc[:,columns_to_keep]

    ##not all plate ids are integers
    df_morph_select['plate.id'] = [str(a) for a in df_morph_select['plate.id']]
    df_morph_select.head()

    start_time = time.time()

    df_reformatted = pd.DataFrame()
    total_number_of_unique_chemicals = 0
    total_number_of_chemical_plate_well = 0

    if (full_devel == "full"):
        # all chemicals
        chemical_id_from_here = pd.unique(df_morph['chemical.id'])
    elif (chem_ind is not None): # full_devel = "deve01l"
        chemical_id_from_here = chem_ind
    else:
        chemical_id_from_here = random.sample(set(pd.unique(df_morph['chemical.id'])), 1)

    ##here we iterate through every chemical
    for chemical_index in chemical_id_from_here:

        print("chemical_index:" + str(chemical_index))
        total_number_of_unique_chemicals += 1
        morph_data_chemical = df_morph_select.loc[df_morph['chemical.id'] == chemical_index,:]

        # Append chemical_plate_well as a unique identifier
        morph_data_chemical.insert(0, 'chemical_plate_well', morph_data_chemical.loc[:,['chemical.id','plate.id', 'well']].apply(lambda x: '_'.join(x.map(str)), axis = 1))

        for cpw in np.unique(morph_data_chemical.chemical_plate_well):
            total_number_of_chemical_plate_well += 1
            temp_df = morph_data_chemical.loc[morph_data_chemical.chemical_plate_well == cpw,:]

            temp_df_grouped = temp_df.groupby(['chemical.id', 'plate.id', 'well'])
            for name, group in temp_df_grouped:

                ## JUSTIFICATION: 7 PAH dataset doesn't have "BRAI" endpoint.
                ## On the other hand, extracts/phase I,II have "BRAI" endpoint.
                if 'BRAI' not in temp_df.endpoint.values: # as 7 PAH
                    try:
                        #            if(len(group.endpoint) == 14):
                        temp = pd.DataFrame( {
                            'chemical.id': pd.unique(temp_df['chemical.id']),
                            'plate.id': pd.unique(temp_df['plate.id']),
                            'well': pd.unique(temp_df['well']),
                            'chemical_plate_well': pd.unique(temp_df['chemical_plate_well']),
                            'conc': pd.unique(temp_df['conc']),
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
                            'TCHR': temp_df.value[temp_df.endpoint == 'TCHR'].values,
                        }  )
                        df_reformatted = pd.concat([df_reformatted, temp])
                    except:
                        print ("len(group.endpoint) != 14")
                        print ("chemical_plate_well:" + str(cpw))
                else: #as extracts
                    temp = pd.DataFrame(
                        {
                        'chemical.id': pd.unique(temp_df['chemical.id']),
                        'plate.id': pd.unique(temp_df['plate.id']),
                        'well': pd.unique(temp_df['well']),
                        'chemical_plate_well': pd.unique(temp_df['chemical_plate_well']),
                        'conc': pd.unique(temp_df['conc']),
                        'AXIS': temp_df.value[temp_df.endpoint == 'AXIS'].values,
                        'BRAI': temp_df.value[temp_df.endpoint == 'BRAI'].values,
                        'CFIN': temp_df.value[temp_df.endpoint == 'CFIN'].values,
                        'CIRC': temp_df.value[temp_df.endpoint == 'CIRC'].values,
                        'DNC_': temp_df.value[temp_df.endpoint == 'DNC_'].values,
                        'DP24': temp_df.value[temp_df.endpoint == 'DP24'].values,
                        'EYE_': temp_df.value[temp_df.endpoint == 'EYE_'].values,
                        'JAW_': temp_df.value[temp_df.endpoint == 'JAW_'].values,
                        'MO24': temp_df.value[temp_df.endpoint == 'MO24'].values,
                        'MORT': temp_df.value[temp_df.endpoint == 'MORT'].values,
                        'NC24': temp_df.value[temp_df.endpoint == 'NC24'].values,
                        'NC__': temp_df.value[temp_df.endpoint == 'NC__'].values,
                        'OTIC': temp_df.value[temp_df.endpoint == 'OTIC'].values,
                        'PE__': temp_df.value[temp_df.endpoint == 'PE__'].values,
                        'PFIN': temp_df.value[temp_df.endpoint == 'PFIN'].values,
                        'PIG_': temp_df.value[temp_df.endpoint == 'PIG_'].values,
                        'SM24': temp_df.value[temp_df.endpoint == 'SM24'].values,
                        'SNOU': temp_df.value[temp_df.endpoint == 'SNOU'].values,
                        'SOMI': temp_df.value[temp_df.endpoint == 'SOMI'].values,
                        'SWIM': temp_df.value[temp_df.endpoint == 'SWIM'].values,
                        'TRUN': temp_df.value[temp_df.endpoint == 'TRUN'].values,
                        'TR__': temp_df.value[temp_df.endpoint == 'TR__'].values,
                        'YSE_': temp_df.value[temp_df.endpoint == 'YSE_'].values,
                        })
                    df_reformatted = pd.concat([df_reformatted, temp])
    print ("total_number_of_unique_chemicals:" + str(total_number_of_unique_chemicals))
    print ("total_number_of_chemical_plate_well:" + str(total_number_of_chemical_plate_well))
    end_time = time.time()
    time_took = str(round((end_time-start_time), 1)) + " seconds"
    print ("Done, it took:"+str(time_took))


    pd.set_option('display.max_columns', None)
    print (df_reformatted.head())
    #print ("df_reformatted.shape:" + str(df_reformatted.shape)) #(288, 19)

    df_reformatted_DNC_0 = pd.DataFrame()

    df_reformatted_DNC_0 = df_reformatted.loc[df_reformatted['DNC_'] == 0.0]
    print ("df_reformatted_DNC_0.shape:" + str(df_reformatted_DNC_0.shape)) #(287, 19)

    df_reformatted_DNC_0.to_csv(output_complete_file_path,index=False)

    print ("Reformat of morpho data is done")
    return chem_ind


def main():
    args = sys.argv[0:]

    print ("reformat for morphological data starts----------------------")

    #complete_file_path = '/Users/kimd999/Dropbox/script/python/srpAnalytics/analysis/paritosh_original_then_edit/to_dockerize/data/7_PAH_zf_morphology_data_2020NOV11_tall.csv'

    complete_file_path = args[1]
    #(from dataQcBmd.py) /srpAnalytics/test_files/7_PAH_zf_morphology_data_2020NOV11_tall.csv
    print ("complete_file_path:" + str(complete_file_path))

    full_devel = args[2]
    output_complete_file_path = complete_file_path[:-4] + "_wide_DNC_0.csv"
    print ("output_complete_file_path:" + str(output_complete_file_path))
    ind = format(complete_file_path, full_devel,output_complete_file_path)



if __name__ == "__main__":
    main()
