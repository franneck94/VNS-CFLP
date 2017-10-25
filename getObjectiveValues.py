import os
from itertools import product

dir_ = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm"
#dir_ = "C:/Users/Administrator/Dropbox/Bachelorarbeit/Programm"
path = dir_+"/Testdaten"
dest = dir_+"/Ergebnisse"

sol_dirs = [file for file in os.listdir(dest) if "cap" in file 
            and file not in ["cap10", "cap11", "cap12", "cap13", "cap14", "cap15"]]
combs = [[2, 3], [30, 40], [0, 1, 2, 3]]

fx = []
fx_all = []
sol = []
min_el, max_el = 10000000000, 0

for sol_dir in sol_dirs:
    gurobi_sol = "_solution"+sol_dir+"TimeLimit.sol"
    print("OPEN: "+ dest+"/"+sol_dir)
    for act_comb in ["3_40_0", "3_40_2"]:
        try:
            with open(dest+"/"+sol_dir+"/"+act_comb+".txt") as f:
                for line in f:
                    if "f(x)" in line:
                        val = float(line.split(" ")[-1])
                        fx.append(val)
                        sol.append(act_comb)
        except:
            pass
    with open(dest+"/"+sol_dir+"/"+gurobi_sol) as f:
        line = f.readline()
        val = float(line.split(" ")[-1])
        fx.append(val)
        sol.append("Gurobi (30s)")

    # Compute vals to Min
    min_el = min(fx)
    fx[:] = [x / min_el for x in fx]
    fx[:] = [round(x, 3) for x in fx]

    # Print out Vals
    for val, comb in zip(fx, sol):
        print(val, "     -     " ,comb)
    fx_all.append(fx)
    fx = []
    sol = []

# sums = [sum(i) for i in zip(*fx_all)]
# for s in sums:
#     print(round(s, 3))
# print(min(sums))