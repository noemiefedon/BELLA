import os
import fnmatch

def find_in_files(directory, find, filePattern):
    print(find)
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            with open(filepath) as f:
                s = f.readlines()
                for ind_line, line in enumerate(s):
                    if find in line:
                        print(" Found in " + filepath \
                              + " at line " + str(ind_line + 1))

directory = r'C:\BELLA'
find = "calc_penalty_oopo_ss("
find_in_files(directory, find, '*.py')