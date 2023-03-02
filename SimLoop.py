import sim_modules as sm
import numpy as np
import matplotlib.pyplot  as plt
import sys
import os

# Large scale iteration code

participant_vec = [100] #75, 100]
timesteps = [200] #[50, 100, 200, 400]
ampContempVec = [0.1]# [0.01, 0.05, 0.1 ]
ampLaggedVec=[0.1] #[0.01, 0.05, 0.1]
ampMeasureVec=[1.0] #[0.1, 0.5, 1, 1.5]
maskVec=[True] #, False]


f = open('input_sample.py', 'r')
lines = f.readlines()

for part in participant_vec:
    for t in timesteps:
        for ampC in ampContempVec:
            for ampL in ampLaggedVec:
                for ampM in ampMeasureVec:
                    for mask in maskVec:
                        print('writing for ', part,t,ampC,ampL,ampM,mask) 
                        g = open('input.py','w')
                        for line in lines:
                            if (line.find('num_participants=100')>-1):
                                newline=line.replace('num_participants=100', 'num_participants=%s'%part)
                                g.write(newline)
                                #print(newline)
                            elif (line.find('steps=100')>-1):
                                newline=line.replace('steps=100','steps=%s'%t)
                                g.write(newline)
                                #print(newline)
                            elif (line.find('ampContemp = 0.1')> -1):
                                newline=line.replace('ampContemp = 0.1','ampContemp = %s'%str(ampC))
                                g.write(newline)
                                #print(newline)
                            elif (line.find('ampLagged = 0.1')>-1):
                                newline=line.replace('ampLagged = 0.1','ampLagged = %s'%str(ampL))
                                g.write(newline)
                                #print(newline)
                            elif (line.find('ampMeasure = 0.1')>-1):
                                newline=line.replace('ampMeasure = 0.1','ampMeasure = %s'%str(ampM))
                                g.write(newline)
                                #print(newline)
                            elif (line.find('maskZero=True')>-1):
                                newline=line.replace('maskZero=True','maskZero = %s'%str(mask))
                                g.write(newline)
                            else: 
                                g.write(line)
                                #print('did not work')
                        g.close()
                        os.system('python SimGenerator.py')




                            
