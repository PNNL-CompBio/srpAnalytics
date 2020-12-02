#f_in = open ("morphological_data_no_old_endpoints_DNC_head.csv", "r")

full_csv = "morphological_data_no_old_endpoints_DNC.csv"
f_in = open (full_csv, "r")

empty_response = 0
filled_response = 0
total_response = 0
zero_response = 0
one_response = 0

for line in f_in:
    num_endpoint = 0
    #print ("line:" + str(line))
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "chemical.id"):
        continue
    for i in range(len(splited_line)):
        if i < 5: # not endpoints
            continue
        num_endpoint += 1
        response = splited_line[i]
        #print ("response:" + str(response))
        total_response += 1
        if (response == "0") or (response == "0.0"):
            zero_response += 1
            filled_response += 1
        elif (response == "1") or (response == "1.0"):
            one_response += 1
            filled_response += 1
        else:
            empty_response += 1

print ("num_endpoint:" + str(num_endpoint))
print ("total_response:" + str(total_response))
print ("empty_response:" + str(empty_response) + ", filled_response:" + str(filled_response))
print ("zero_response:" + str(zero_response) + ", one_response:" + str(one_response))
