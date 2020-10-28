'''
f = open("guru.txt","a+")
for i in range(2):
	f.write("Append")
f.close()
'''

#'''
final_count = 0
f_out = open("qc_1_2_case.txt")
for line in f_out:
	final_count = str(line)
f_out.close()

f_out = open("qc_1_2_case.txt", 'a+')
f_out.write(str(final_count))
                        #write_this = int(line)+1
                        #f_out.write(str(write_this))
f_out.close()
#'''
