from gurobipy import *
import os
import numpy as np

dir_ = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm"
#dir_ = "C:/Users/Administrator/Dropbox/Bachelorarbeit/Programm"
path = dir_+"/Testdaten"
dest = dir_+"/Ergebnisse"

test_set = [test for test in os.listdir(path) if "capb" in test or "capc" in test]
#test_set = ["cap16"]

for testset in test_set:
    print("Compute ",testset)
    demand, capacity, fixedCosts = [], [], []

    # Warehouse demand in thousands of units
    with open(path+"/"+testset+"/dj.txt", "r") as f:
        for line in f:
            demand.append(float(line.split(" ")[-1]))

    # Plant capacity in thousands of units
    with open(path+"/"+testset+"/bi.txt", "r") as f:
        for line in f:
            capacity.append(float(line.split(" ")[-1]))

    # Fixed costs for each plant
    with open(path+"/"+testset+"/fi.txt", "r") as f:
        for line in f:
            fixedCosts.append(float(line.split(" ")[-1]))

    m_locNum = len(capacity)
    m_cusNum = len(demand)
    cij = np.zeros(m_locNum * m_cusNum)

    # Transportation costs per thousand units
    with open(path+"/"+testset+"/cij.txt", "r") as f:
        for i,line in zip(range(m_locNum * m_cusNum) ,f):
            cij[i] = float(line.split(" ")[-1])

    shape = (m_locNum, m_cusNum)
    cij = np.reshape(cij, shape)

    # Range of plants and warehouses
    facility = range(len(capacity))
    customer = range(len(demand))

    m = Model()

    m.setParam('TimeLimit', 30)

    y = m.addVars(facility, vtype=GRB.BINARY, obj=fixedCosts, name="y")
    x = m.addVars(facility, customer, obj=cij, name="x")
    # Production constraints
    # Note that the right-hand limit sets the production to zero if the plant
    # is closed
    m.addConstrs( (x[i, j] - y[i]) <= 0 for i in facility for j in customer )
 
    m.addConstrs( sum(x[i, j] for i in facility) == 1  for j in customer )
    m.addConstrs( sum(demand[j] * x[i,j] for j in customer) <= capacity[i] * y[i] for i in facility )

    obj = sum(x[i,j] * cij[i, j] for i in facility for j in  customer) + sum(y[i] * fixedCosts[i] for i in facility)

    m.setObjective(obj ,GRB.MINIMIZE)

    m.optimize()

    # Print solution
    print('\nTOTAL COSTS: %g' % m.objVal)
    runtime = m.runtime
    m.write(dest+"/"+testset+"/_solution"+testset+"TimeLimit.sol")
    # Append Runtime to File
    with open(dest+"/"+testset+"/_solution"+testset+"TimeLimit.sol", "a") as f:
        end_line = "Runtime: "+str(runtime)
        f.write(end_line)