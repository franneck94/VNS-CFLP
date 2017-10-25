import os
from math import isclose

def roundoff_error(exact, approximate):
    if (exact == 0.0 or approximate == 0.0):
        return abs(exact + approximate)
    else:
        return abs(approximate/exact - 1.0)

def float_equal(float1, float2, epsilon=2.0e-9):
    return (roundoff_error(float1, float2) < epsilon)

dir_ = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm"
#dir_ = "C:/Users/Administrator/Dropbox/Bachelorarbeit/Programm"
path = dir_+"/Testdaten"
dest = dir_+"/Ergebnisse"

sol_dirs = [d for d in os.listdir(dest) if "cap" in d]

others = []
y = []
x = []

for sol_dir in sol_dirs:
    sols = [s for s in os.listdir(dest+"/"+sol_dir+"/") if "_solutioncap" in s and "Filter" not in s]
    for set_ in sols: 
        print(set_)
        try:
            with open(dest+"/"+sol_dir+"/"+set_) as f:
                for line in f:
                    if line.startswith("y"):
                        val = line.split(" ")[-1]
                        if float(val) != 0.0:
                            y.append(line)
                    elif line.startswith("x"):
                        val = line.split(" ")[-1]
                        if float(val) != 0.0 and not float_equal(0.0, float(val)):
                            x.append(line)
                    else:
                        others.append(line)

                with open(dest+"/"+sol_dir+"/"+set_, "w") as f:
                    for line in others:
                        f.write(line)
                    f.write("\n")
                    for line in y:
                        f.write(line)
                    f.write("\n")
                    for line in x:
                        f.write(line)
                others = []
                y = []
                x = []
        except:
            pass