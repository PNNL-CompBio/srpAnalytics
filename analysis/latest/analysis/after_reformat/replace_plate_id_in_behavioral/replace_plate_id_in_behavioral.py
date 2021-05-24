# purpose: replace plate.id in behavior data with bottle.id in morphology data

import subprocess, time

def show_time(process, time_start, time_end):
    time_took = str(process) + " finished in "
    if (round((time_end-time_start)/60, 1) < 1):
      time_took = time_took + str(round((time_end-time_start), 1)) + " seconds "
    elif (round((time_end-time_start)/60/60, 1) < 1):
      time_took = time_took + str(round((time_end-time_start)/60, 1)) + " minutes "
    else:
      time_took = time_took + str(round((time_end-time_start)/60/60, 1)) + " hours "
    time_took = time_took + "(wall clock)"
    return time_took
############### end of show_time function

start_time = time.time()

# used to work
#beha_full_csv = "344_zf_LPR_data_phase_1_and_2_-_2020JUNE25.csv"

beha_full_csv = "344_zf_LPR_data_phase_1_2_2020JUNE25_updated_plate_id_for_TX_tall_fixed_full_15_timepoints_wide_t3_t17_full.csv"


beha_full_in = open (beha_full_csv, "r")

beha_full_csv_after_replacement = beha_full_csv[:-4] + "_replaced.csv"


for line in beha_full_in:
    beha_full_out = open (beha_full_csv_after_replacement, "a")
    splited_line = line.split(",")
    bottle_id = splited_line[1]
    #print ("bottle_id:" + str(bottle_id))
    if (bottle_id == "\"bottle.id\""):
        beha_full_out.write(line)
        continue
    #print ("bottle_id[:3]:" + str(bottle_id[:3]))
    if (bottle_id[:3] != "\"TX"):
        beha_full_out.write(line)
        continue
    old_plate_id = splited_line[3]
    
    #command = "grep -m 1 " + str(bottle_id) + " zf_morphology_data_335_chemicals_2020DEC16.csv"
    command = "grep -m 1 " + str(bottle_id) + " zf_morphology_data_335_chemicals_2020DEC16_fixed.csv"
    
    #print ("\ncommand:" + str(command))
    subprocessed = subprocess.check_output(command, shell=True)
    #print ("subprocessed:" + str(subprocessed))
    subprocessed = subprocessed.decode('UTF-8')
    splited_subprocessed = subprocessed.split(",")
    new_plate_id = splited_subprocessed[3]
    #print ("old_plate_id:" + str(old_plate_id))
    #print ("new_plate_id:" + str(new_plate_id))
    
    new_line = ''
    for i in range (len(splited_line)):
        if (i == 0):
            new_line = str(splited_line[i])
        elif (i != 3):
            new_line = new_line + "," + str(splited_line[i])
        else:
            new_line = new_line + "," + str(new_plate_id)
    
    #print ("line:" + str(line))
    #print ("new_line:" + str(new_line))
    beha_full_out.write(new_line)
    beha_full_out.close()
    
beha_full_in.close()


end_time = time.time()
write_this = show_time("replacement of plate id", start_time, end_time)
'''
In constance, it took 8 hrs to transform
344_zf_LPR_data_phase_1_and_2_-_2020JUNE25.csv
into
344_zf_LPR_data_phase_1_and_2_-_2020JUNE25_replaced.csv

If I had used wide format, it would have taken much shorter time.
'''