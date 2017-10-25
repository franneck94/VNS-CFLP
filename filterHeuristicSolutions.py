import os

path = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Ergebnisse"
test_sets = [file for file in os.listdir(path) if "cap" in file]

elements = []
solution_val = 0.0
optimal_solution_val = 0.0

for testset in test_sets:
    sol_file = path+"/"+testset+"/"+"_solution"+testset+".sol"
    with open(sol_file) as f:
        for line in f:
            if "# Objective value = " in line:
                optimal_solution_val = float(line.split(" ")[-1])
    files = [file for file in os.listdir(path+"/"+testset) if "solution" not in file]
    for file in files:
        with open(path+"/"+testset+"/"+file, 'r') as f:
            for line in f:
                if "Soltuion:" in line:
                    solution_val = float(line.split(" ")[-1])
        if solution_val != 0.0:
            gap = (solution_val - optimal_solution_val) / optimal_solution_val
            if gap > 0.02:
                try:
                    if os.path.exists(path+"/"+testset+"/"+file):
                        os.remove(path+"/"+testset+"/"+file)
                except:
                    print("Cant remove")
        with open(path+"/"+testset+"/"+file, "a") as f:
            f.write("\n")
            f.write("GAP : "+ str(gap))
        solution_val = 0.0