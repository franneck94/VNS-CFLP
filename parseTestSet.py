import os
import errno
import itertools

directory = 'C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten/Raw DataSet/'
# listdir = [file for file in os.listdir(directory) if file not in ['capa.txt', 'capb.txt', 'capc.txt']]

# for d in listdir:
#     print('Opening dir: ', directory+'/'+d)
#     with open(directory+'/'+d) as f:
#         firstLine = True
#         numLocs = 0
#         numCusts = 0
#         text  = ''
#         f_i = []
#         d_j = []
#         b_i = []
#         first_line = f.readline()
#         numLocs = int(first_line.split()[0])
#         numCusts = int(first_line.split()[1])
#         c_ij = [[(0, 0, 0) for cus in range(numCusts)]  for loc in range(numLocs) ]
#         for i, line in enumerate(f):
#             if i < numLocs:
#                 b_i.append((i, float(line.split()[0])))
#                 f_i.append((i, float(line.split()[1])))
#             else:
#                 for number in line.split():
#                     text += ' ' + number
#         text = text[1:]
#         for index, item in enumerate(text.split(' ')):
#             if index % (numLocs+1) == 0:
#                 d_j.append((index, float(item)))
#         text = [val for index, val in enumerate(text.split(' ')) if  index % (numLocs+1) != 0]
#         for customer in range(numCusts):
#             firstLine = True
#             actual_allocating_costs = []
#             for index, val in enumerate(text[:numLocs]):
#                 c_ij[index][customer] = (index, customer, float(val))
#             if len(text) > numLocs:
#                 text = text[numLocs:]

# directory = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten"
# for d in listdir:
#     files_to_save = ['cij.txt', 'dj.txt', 'bi.txt', 'fi.txt']
#     for file, data in zip(files_to_save, [c_ij, d_j, b_i, f_i]):
#         print(directory+'/'+d.split('.')[0]+'/'+file)
#         os.makedirs(os.path.dirname(directory+'/'+d.split('.')[0]+'/'+file), exist_ok=True)
#         with open(directory+'/'+d.split('.')[0]+'/'+file, "w") as f:
#             if len(data[0]) == 2:
#                 for val in data:
#                     f.write(' '.join(map(str, val))+'\n')
#             else:
#                 for val_i in data:
#                     for val_j in val_i:
#                         f.write(' '.join(map(str, val_j)) + '\n')    
# 

correct = [56, 94, 89]
listdir = [file for file in os.listdir(directory) if file in ['capa.txt', 'capb.txt', 'capc.txt']]   
capacity_amount = 8000.0

for c, d in zip(range(listdir), listdir):
    print('Opening dir: ', directory+d)
    with open(directory+'/'+d) as f:
        #Init Vars
        firstLine = True
        numLocs = 0
        numCusts = 0
        text  = ''
        f_i = []
        d_j = []
        b_i = []
        first_line = f.readline()
        numLocs = int(first_line.split()[0])
        numCusts = int(first_line.split()[1])
        c_ij = []
        
        # Start Parsing
        for i, line in enumerate(f):
            if i <= numLocs and i > 0:
                b_i.append(capacity_amount)
                f_i.append(line.split(" ")[2])
            else:
                for number in line.split(" "):
                    text += " " + number
        
        text_list = [item for index, item in enumerate(text.split()) if item != " " and item != "\n" and item != "capacity"]
        text_list[c] = correct[c]

        for item,counter in zip(text_list, range(len(text_list))):
            if counter % (numLocs+1) == 0:
                d_j.append(item)
            else:
                if " " in item:
                    c_ij.append(item.split(" ")[-1])
                else:
                    c_ij.append(item)
            

directory = "C:/Users/Jan/Dropbox/Bachelorarbeit/Programm/Testdaten"
for d in listdir:
    files_to_save = ['cij.txt', 'dj.txt', 'bi.txt', 'fi.txt']
    for file, data in zip(files_to_save, [c_ij, d_j, b_i, f_i]):
        print(directory+'/'+d.split('.')[0]+'/'+file)
        os.makedirs(os.path.dirname(directory+'/'+d.split('.')[0]+'/'+file), exist_ok=True)
        with open(directory+'/'+d.split('.')[0]+'/'+file, "w") as f:
            for val in data:
                f.write(str(val)+'\n')