import glob, os, random, shutil, sys, time

if (__name__ == "__main__"):
    start_time = time.time()
    args = sys.argv[1:]
    print (len(args))
    if (len(args) < 1):
        print ("Specify py filename \n")
        print ("For example, python display_to_print.py wide2dicho_LPR_7_PAH_t0_t239.py\n")
        exit(1)
    
    args = sys.argv[1:]
    input_py_file_name = args[0]
    
    new_py_file_name = input_py_file_name[:-3] + "_display_to_print.py"
    if (os.path.isfile(new_py_file_name) == True):
        os.remove(new_py_file_name)
    new_py_file = open(new_py_file_name, 'w')
    
    input_file = open(input_py_file_name, 'r')
    for line in input_file:
        if (line[:7] == "display"):
            new_line = "print" + line[7:]
            new_py_file.write(new_line)
        else:
            new_py_file.write(line)
    input_file.close()
    new_py_file.close()

end_time = time.time()
time_took = str(round((end_time-start_time), 1)) + " seconds "

print ("display_to_print" + str(time_took))
