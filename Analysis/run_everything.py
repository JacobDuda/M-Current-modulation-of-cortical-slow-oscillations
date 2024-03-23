import subprocess
import numpy as np


#code used for coarse search
'''
#these are all the values that we will test the parameters at
vals = [0,0.33,0.66,0.99]

#this array will contain all the combinations we will run through
cs = np.zeros((64,4))
#loop through all possible combinations
counter = 0
for i in range(4):
    for j in range(4):
        for k in range(4):
            cs[counter] = np.array((vals[i],vals[j],vals[k],counter))
            counter += 1
'''

#code used for fine search
'''
#these are all the values that we will test the parameters at
Mvals = [0,0.1,0.2,0.3]
Hvals = [0.7,0.8,0.9,1]
Nvals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7]

#this array will contain all the combinations we will run through
cs = np.zeros((128,4))
#loop through all possible combinations
counter = 0


for i in range(4):
    for j in range(4):
        for k in range(8):
            cs[counter] = np.array((Mvals[i],Hvals[j],Nvals[k],counter))
            counter += 1


for c in cs[62:]:
    command  = r"./exec "+str(c[0])+" "+str(c[1])+" "+str(c[2])+ " 123456 "+str(int(c[3]))
    print("running "+command)
    results = subprocess.call(command)
    print("Finished "+str(c[3]+1)+" simulations...")
'''

#code used for individual parameter exploration
vals = np.arange(0,1.1,0.1)

cs = np.zeros((11,4))

for i in range(11):
    cs[i] = np.array((1,1,vals[i],i))

for c in cs:
    command  = str("./exec "+str(np.round(c[0],1))+" "+str(np.round(c[1],1))+" "+str(np.round(c[2],1))+ " 123456 "+str(int(c[3])))
    print("running "+command)
    results = subprocess.call(command)
    print("Finished "+str(c[3]+1)+" simulations...")



