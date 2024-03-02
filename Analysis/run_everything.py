import subprocess
import numpy as np

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

for c in cs[42:]:
    command  = r"./exec "+str(c[0])+" "+str(c[1])+" "+str(c[2])+ " 123456 "+str(int(c[3]))
    print("running "+command)
    results = subprocess.call(command)
    print("Finished "+str(c[3]+1)+" simulations...")

    


