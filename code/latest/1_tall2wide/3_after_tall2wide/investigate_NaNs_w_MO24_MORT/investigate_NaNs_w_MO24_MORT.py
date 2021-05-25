#f_in = open ("head.csv", "r")
f_in = open ("/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/input/wide/344_zf_morphology_data_phase_1_and_2_-_2020JUNE25_wide_DNC_0.csv", "r")

nan_except_MO24 = 0
nan_except_MORT_SM24_DP24_NC24_MO24 = 0

filename_chemical_plate_well_MO24 = "has_at_least_1_non_nan_response_other_than_MO24_even_when_MO24_eq_1.csv"
f_MO24 = open (filename_chemical_plate_well_MO24, "w")

filename_chemical_plate_well_MORT = "has_at_least_1_non_nan_response_other_than_24series_and_MORT_even_when_MORT_eq_1.csv"
f_MORT = open (filename_chemical_plate_well_MORT, "w")
for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "chemical.id"):
        f_MO24.write(str(line))
        f_MORT.write(str(line))
    chemical_plate_well = splited_line[3]
    MO24 = splited_line[13]
    if (MO24 == "1.0") or (MO24 == "1"):
        empty_response = 0
        non_nan_response = 0
        for i in range(len(splited_line)):
            if i < 5: # not endpoints
                continue
            elif i ==9: # DNC_
                continue
            elif i ==13: # MO24
                continue
            response = splited_line[i].rstrip()
            if (response == ''):
                empty_response += 1
            else:
                #print ("\n[dealing MO24] Not empty!")
                #print ("i=" + str(i) + ". Its response=" + str(response))
                non_nan_response += 1
        if (non_nan_response != 0):
            #print ("There is at least 1 non-nan response (other than MO24) even when MO24=1. Its chemical_plate_well: " + str(chemical_plate_well))
            f_MO24.write(str(line))
        nan_except_MO24 = nan_except_MO24 + empty_response
    
    MORT = splited_line[14]
    if (MORT == "1.0") or (MORT == "1"):
        empty_response = 0
        non_nan_response = 0
        for i in range(len(splited_line)):
            if i < 5: # not endpoints
                continue
            elif i ==9: # DNC_
                continue
            elif i ==10: # DP24
                continue
            elif i ==13: # MO24
                continue
            elif i ==14: # MORT
                continue
            elif i ==15: # NC24
                continue
            elif i ==21: # SM24
                continue
            response = splited_line[i].rstrip()
            if (response == ''):
                empty_response += 1
            else:
                print ("\n[dealing MORT] Not empty!")
                print ("i=" + str(i) + ". Its response=" + str(response))
                non_nan_response += 1
        if (non_nan_response != 0):
            print ("There is at least 1 non-nan response (other than 24 series and MORT) even when MORT=1. Its chemical_plate_well: " + str(chemical_plate_well))
            f_MORT.write(str(line))
        nan_except_MORT_SM24_DP24_NC24_MO24 = nan_except_MORT_SM24_DP24_NC24_MO24 + empty_response

f_MO24.close()
f_MORT.close()
print ("\nNumber of nan entries that are not MO24 (nan_except_MO24):" + str(nan_except_MO24) )
print ("Number of nan entries that are not MORT_SM24_DP24_NC24_MO24 (nan_except_MORT_SM24_DP24_NC24_MO24):" + str(nan_except_MORT_SM24_DP24_NC24_MO24) )
