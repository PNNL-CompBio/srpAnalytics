f_in = open ("working_example/input/344_zf_morphology_data_phase_1_and_2_-_2020JUNE25_wide_DNC_0_head.csv", "r")
#f_in = open ("/Users/kimd999/research/projects/toxicity/per_each_data/Phase_I_II/input/wide/344_zf_morphology_data_phase_1_and_2_-_2020JUNE25_wide_DNC_0.csv", "r")
empty_response = 0
filled_response = 0
for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "chemical.id"):
        continue
    for i in range(len(splited_line)):
        if i < 5: # not endpoints
            continue
        elif i ==9: # DNC_
            continue
        response = splited_line[i]
        if (response == ''):
            empty_response += 1
        else:
            filled_response += 1
print ("empty_response:" + str(empty_response) + ", filled_response:" + str(filled_response))
