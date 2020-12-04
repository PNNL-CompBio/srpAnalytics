import numpy as np
import pandas as pd
import os, time
from matplotlib import pyplot as plt
import seaborn as sns


##### <begin> sort by 'Chemical_ID', 'End_Point'
complete_file_path = 'working_example/input/bmd_vals_all_qc_analyzed.csv'
#complete_file_path = 'devel.csv'
morphological_data = pd.read_csv(complete_file_path, header = 0)

sorted_morphological_data = morphological_data.sort_values(['Chemical_ID', 'End_Point'], ascending=[1, 1])

sorted_morphological_data = sorted_morphological_data[['Chemical_ID', 'End_Point', 'BMD10_call']]
#print (sorted_morphological_data)

filename_sorted_df = 'sorted.csv'
sorted_morphological_data.to_csv(filename_sorted_df, index=False)
##### <end> sort by 'Chemical_ID', 'End_Point'


##### <begin> prepare heatmap array
heatmap_array = []
for i in range(0, len(np.unique(sorted_morphological_data['Chemical_ID']))+1):
  new = []
  for j in range(0, len(np.unique(sorted_morphological_data['End_Point']))+1):
      new.append('')
  heatmap_array.append(new)
print (len(heatmap_array))
print (len(heatmap_array[0]))

Chemical_ID_index = 1
old_Chemical_ID = -999
End_Point_index = 1 # initial
f_in = open (filename_sorted_df, "r")
for line in f_in:
    #print ("\nline:" + str(line))
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "Chemical_ID"):
      continue
    if (old_Chemical_ID == -999):
      heatmap_array[Chemical_ID_index][0] = str(Chemical_ID)
    elif (old_Chemical_ID != Chemical_ID):
      Chemical_ID_index += 1
      End_Point_index = 1 # re-initial      
      heatmap_array[Chemical_ID_index][0] = str(Chemical_ID)
    
    End_Point = splited_line[1]
    
    heatmap_array[0][End_Point_index] = str(End_Point)
    
    BMD10_call = splited_line[2].rstrip()
    #print ("Chemical_ID_index:" + str(Chemical_ID_index))
    #print ("End_Point_index:" + str(End_Point_index))
    heatmap_array[Chemical_ID_index][End_Point_index] = str(BMD10_call)
    #print ("heatmap_array:" + str(heatmap_array))
    
    if (old_Chemical_ID != Chemical_ID):
      old_Chemical_ID = Chemical_ID
      
    End_Point_index += 1
f_in.close()

DF = pd.DataFrame(heatmap_array)
DF.to_csv("ready_for_heatmap.csv", index=False, header=False)
##### <end> prepare heatmap




'''
df_unique = pd.DataFrame()
for Chemical_ID in np.unique(sorted_morphological_data['Chemical_ID']):
  sorted_morphological_data_per_chemical = sorted_morphological_data.loc[sorted_morphological_data['Chemical_ID'] == Chemical_ID,['Chemical_ID', 'End_Point', 'BMD10_call']]
  #print (str(sorted_morphological_data))
  for End_Point in np.unique(sorted_morphological_data_per_chemical['End_Point']):
    sorted_morphological_data_per_chemical_endpoint = sorted_morphological_data_per_chemical.loc[sorted_morphological_data_per_chemical['End_Point'] == End_Point,['Chemical_ID', 'End_Point', 'BMD10_call']]
    #print (str(sorted_morphological_data_per_chemical_endpoint))
    df_unique = df_unique.append(sorted_morphological_data_per_chemical_endpoint, ignore_index=True)
    print (str(df_unique))

filename_unique = "sorted_morphological_data_per_chemical_endpoint.csv"
if (os.path.isfile(filename_unique) == True):
  os.remove(filename_unique)
df_unique.to_csv(filename_unique, index=False, header=False)



##### <begin> prepare heatmap array
heatmap_array = []
for i in range(0, len(sorted_morphological_data)+1):
  new = []
  for j in range(0, len(sorted_morphological_data)+1):
      new.append('')
  heatmap_array.append(new)

Chemical_ID_index = 1
End_Point_index = 1
f_in = open (filename_sorted_df, "r")
for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "Chemical_ID"):
      continue
    heatmap_array[Chemical_ID_index][0] = str(Chemical_ID)
    
    End_Point = splited_line[1]
    heatmap_array[0][End_Point_index] = str(End_Point)
    
    BMD10_call = splited_line[2]
    heatmap_array[Chemical_ID_index][End_Point_index] = str(BMD10_call)
    
    Chemical_ID_index += 1
    End_Point_index += 1
f_in.close()

DF = pd.DataFrame(heatmap_array)
DF.to_csv("arr.csv", index=False, header=False)
##### <end> prepare heatmap
'''


'''
##### <begin> prepare heatmap array
f_in = open (filename_sorted_df, "r")
filename_ready_for_heatmap = filename_sorted_df[:-4] + "_read_for.csv"
f_out = open (filename_sorted_df, "w")
f_out.write(",\n") # just for very first cell
for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    f_out.write(str(Chemical_ID) + ",")
    End_Point = splited_line[1]
    BMD10_call = splited_line[2]
f_in.close()
f_out.close()

##### <end> prepare heatmap
'''