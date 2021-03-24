import os
import fnmatch

def findReplace(directory, find, replace, filePattern):
    for path, dirs, files in os.walk(os.path.abspath(directory)):
        for filename in fnmatch.filter(files, filePattern):
            filepath = os.path.join(path, filename)
            with open(filepath) as f:
                s = f.read()
            s = s.replace(find, replace)
            with open(filepath, "w") as f:
                f.write(s)

find = "from src.LAYLA_V02"
replace = "from src.LAYLA_V02"

#findReplace(r'C:\BELLA', find, replace, '*.py')
#findReplace(r'C:\LAYLA', find, replace, '*.py')
#findReplace(r'C:\RELAY', find, replace, '*.py')
findReplace(r'C:\BELLA', find, replace, '*.py')