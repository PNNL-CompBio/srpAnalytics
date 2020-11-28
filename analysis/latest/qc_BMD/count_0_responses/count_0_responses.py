import csv, os, sys

dic_1st = {}
#f_in = open ("dose_response_vals_qc_1_new_head.csv", "r")
f_in = open ("dose_response_vals_qc_1.csv", "r")
for line in f_in:
    splited_line = line.split(",")
    Chemical_ID = splited_line[0]
    if (Chemical_ID == "Chemical_ID"):
        continue
    End_Point = splited_line[1]
    combi = (str(Chemical_ID), str(End_Point), "Dose_count")
    if combi not in dic_1st.keys():
        #print ("combi not in dic_1st.keys()")
        Dose_count = 0
        zero_Response = 0
    Dose_count += 1
    dic_1st[str(Chemical_ID), str(End_Point), "Dose_count"] = Dose_count
    
    Response = splited_line[3]
    if (float(Response) < 0.00001):
        zero_Response += 1
    dic_1st[str(Chemical_ID), str(End_Point), "zero_Response"] = zero_Response
f_in.close()
print ("dic_1st:\n"+str(dic_1st))

report_filename = 'Dose_count_zero_Response.csv'
report_file_out = open(report_filename, "w")
report_file_out.write("Dose_count,zero_Response,one_response,zero_Response/Dose_count\n")

for key, value in dic_1st.items():
    if ((key[2]) == "Dose_count"):
        Dose_count = value
        write_this = str(value) + ","
    elif ((key[2]) == "zero_Response"):
        zero_Response = value
        write_this = write_this + str(zero_Response) + "," + str(Dose_count-zero_Response) \
                    + "," + str(round(100*(zero_Response/Dose_count),0)) + "%\n"
        report_file_out.write(write_this)
        write_this = ''
report_file_out.close()

dic_2nd_to_merge = {}
report_file_in = open(report_filename, "r")
for line in report_file_in:
    splited_line = line.split(",")
    Dose_count = splited_line[0]
    if (Dose_count == "Dose_count"):
        continue
    
    zero_Response = int(splited_line[1])
    combi = (str(Dose_count), "zero_Response")
    if combi not in dic_2nd_to_merge.keys():
        dic_2nd_to_merge[str(Dose_count), "zero_Response"] = zero_Response
    else:
        dic_2nd_to_merge[str(Dose_count), "zero_Response"] = dic_2nd_to_merge[str(Dose_count), "zero_Response"]+zero_Response
        combi = (str(Dose_count), "zero_Response")
    
    one_response = int(splited_line[2])
    combi = (str(Dose_count), "one_response")
    if combi not in dic_2nd_to_merge.keys():
        dic_2nd_to_merge[str(Dose_count), "one_response"] = one_response
    else:
        dic_2nd_to_merge[str(Dose_count), "one_response"] = dic_2nd_to_merge[str(Dose_count), "one_response"]+one_response
print ("dic_2nd_to_merge:\n"+str(dic_2nd_to_merge))