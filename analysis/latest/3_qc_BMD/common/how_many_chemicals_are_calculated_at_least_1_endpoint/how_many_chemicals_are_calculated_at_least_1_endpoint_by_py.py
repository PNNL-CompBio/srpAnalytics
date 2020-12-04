import csv, os, sys
import numpy as np
import pandas as pd

#### <begin> make bmd_call
filename_w_bmd_call = "/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/output/old_less_detailed_report/all_qc/bmd_vals_all_qc.csv"
#filename_w_bmd_call = "/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/output/old_less_detailed_report/head.csv"
f_in = open (filename_w_bmd_call, "r")

filename_w_bmd_call_analyzed = filename_w_bmd_call[:-4] + "_analyzed.csv"
f_out = open (filename_w_bmd_call_analyzed, "w")

for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    BMD10 = splited_line[3]
    if (Chemical_ID == "Chemical_ID"):
        write_this = line.rstrip() + ",BMD10_call\n"
    else:
        try:
            print (type(float(BMD10))) # essential
            write_this = line.rstrip() + ",1\n"
        except:
            write_this = line.rstrip() + ",0\n"            
    f_out.write(write_this)
f_in.close()
f_out.close()
#### <end> make bmd_call


# <begin> count how many chemicals are calculated at least 1 endpoint
chemical_w_at_least_1_bmd10 = 0
chemical_w_no_bmd10 = 0
df_per_Chemical_ID_w_no_bmd10 = pd.DataFrame()

filename_chemical_w_no_bmd10 = "chemical_w_no_bmd10.csv"
if (os.path.isfile(filename_chemical_w_no_bmd10) == True):
    os.remove(filename_chemical_w_no_bmd10)
f_out = open (filename_chemical_w_no_bmd10, "a")
header_written = False # initial

df_w_bmd_call = pd.read_csv(filename_w_bmd_call_analyzed, header = 0)
for Chemical_ID in np.unique(df_w_bmd_call['Chemical_ID']):
    #print ("Chemical_ID:" + str(Chemical_ID))
    df_per_Chemical_ID = df_w_bmd_call.loc[df_w_bmd_call['Chemical_ID'] == Chemical_ID,['Chemical_ID', 'End_Point', 'BMD10_call']]
    
    calculated_at_least_1_bmd10 = False # initial
    for End_Point in np.unique(df_per_Chemical_ID['End_Point']):
        if (calculated_at_least_1_bmd10 == True):
            break
        df_per_Chemical_ID_End_Point = df_per_Chemical_ID.loc[df_per_Chemical_ID['End_Point'] == End_Point,['Chemical_ID', 'End_Point', 'BMD10_call']]
        #print ("df_per_Chemical_ID_End_Point:\n" + str(df_per_Chemical_ID_End_Point))
        for (columnName, columnData) in df_per_Chemical_ID_End_Point.iteritems():
            #print('Colunm Name : ', columnName) 
            #print('Column Contents : ', columnData.values)
            if (columnName == "BMD10_call"):
                if (columnData.values == 1):
                    calculated_at_least_1_bmd10 = True
                    break
    
    ## assessment is done
    if (calculated_at_least_1_bmd10 == True):
        chemical_w_at_least_1_bmd10 += 1
    else:
        df_write_this = df_w_bmd_call.loc[df_w_bmd_call['Chemical_ID'] == Chemical_ID]
        if (header_written == False):
            header_written = True
            df_write_this.to_csv(filename_chemical_w_no_bmd10, mode='a', header=True, index=False)
        else:
            df_write_this.to_csv(filename_chemical_w_no_bmd10, mode='a', header=False, index=False)
        
        chemical_w_no_bmd10 += 1
    calculated_at_least_1_bmd10 = False # re-initial
f_out.close()


print ("chemical_w_at_least_1_bmd10:" + str(chemical_w_at_least_1_bmd10))
print ("chemical_w_no_bmd10:" + str(chemical_w_no_bmd10))
# <end> count how many chemicals are calculated at least 1 endpoint            
    