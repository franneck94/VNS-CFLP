import os
import errno
import itertools
import numpy as np
import random

directory = 'C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten/'
from_ =  'C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten/'

m = 0
n = 0

demand, capacity, fixedCosts = [], [], []

######## READ IN CAPA DATA ##########

# Warehouse demand in thousands of units
with open(from_+"/"+"capa"+"/dj.txt", "r") as f:
    for line in f:
        demand.append(float(line.split(" ")[-1]))

# Plant capacity in thousands of units
with open(from_+"/"+"capa"+"/bi.txt", "r") as f:
    for line in f:
        capacity.append(float(line.split(" ")[-1]))

# Fixed costs for each plant
with open(from_+"/"+"capa"+"/fi.txt", "r") as f:
    for line in f:
        fixedCosts.append(float(line.split(" ")[-1]))

m_locNum = len(capacity)
m_cusNum = len(demand)
cij = np.zeros(m_locNum * m_cusNum)

# Transportation costs per thousand units
with open(from_+"/"+"capa"+"/cij.txt", "r") as f:
    for i,line in zip(range(m_locNum * m_cusNum) ,f):
        cij[i] = float(line.split(" ")[-1])

shape = (m_locNum, m_cusNum)
cij = np.reshape(cij, shape)

###### START CREATING TESTSETS ########
# 19, 20, 21, 50, 51, 52, 53, 54, 55, 56

l = ["g"]

to_create = ["cap"+str(number) for number in l]

for testset, number in zip(to_create, l):

    m = 100
    n = 1000

    capacity = 8000.0

    with open(directory + testset + '/fi.txt', 'w') as f:
        for i in range(m):
            if i < m-1:
                f.write(str(random.choice(fixedCosts))+ "\n")
            else:
                f.write(str(random.choice(fixedCosts)))
    with open(directory + testset + '/bi.txt', 'w') as f:
        for i in range(m):
            if i < m-1:
                f.write(str(capacity)+ "\n")
            else:
                f.write(str(capacity))
    with open(directory + testset + '/dj.txt', 'w') as f:
        for i in range(n):
            if i < n-1:
                f.write(str(random.choice(demand))+ "\n")
            else:
                f.write(str(random.choice(demand)))
    with open(directory + testset + '/cij.txt', 'w') as f:
        for i in range(m * n):
            if i < m*n -1:
                f.write(str(random.choice(random.choice(cij)))+ "\n")
            else:
                f.write(str(random.choice(random.choice(cij))))