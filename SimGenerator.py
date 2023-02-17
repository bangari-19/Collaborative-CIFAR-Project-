import sim_modules as sm
import numpy as np
import matplotlib.pyplot  as plt
import sys 
import input as input 
import os

size= input.size
pathin= input.pathin
pathout=input.pathout
cols = 2*size
print('reading from ', input.num_participants, input.steps, input.ampContemp, input.ampLagged, input.ampMeasure, input.maskZero)


for j in range(input.num_iterations):
    savepath=pathout+'sims/numParticipants_%s/steps_%s/Contemp_%s_amp_%s/Lagged_%s_amp_%s/Measure_%s_amp_%s/mask_%s/rep_%i/'%(input.num_participants,str(input.steps),input.covContempName,str(input.ampContemp),input.covLaggedName, str(input.ampLagged),input.measurecovName, str(input.ampMeasure),str(input.maskZero),j)
    try:
        os.makedirs(savepath)
        print('making %s'%savepath)
    except:
        print('directory %s exists'%savepath)

    for number in range(input.num_participants):
        csvnum = number+1
        file = pathin+'%i.csv'%csvnum
        data = np.loadtxt(file, skiprows=1, usecols=range(0,cols), delimiter=',')
        matContemp = data[:,size:] #same day (contemporaneous) beta values
        matLagged = data[:,:size] #lagged beta values
        
        if input.covContempName=='randn':
            covContemp = np.random.randn(size,size)
        else:
            covContemp = np.ones(size)

        if input.covLaggedName=='randn':
            covLagged = np.random.randn(size,size)
        else:
            covLagged =np.ones(size)

        if input.measurecovName=='diag':
            measureCov= input.ampMeasure*np.eye(size)
        
        if input.maskZero:
            maskContemp =  sm.make_mask(matContemp, contemp=True)
            maskLagged = sm.make_mask(matLagged, contemp=False)
        else:
            maskContemp =  np.ones((input.size, input.size))
            np.fill_diagonal(maskContemp,0) # still make diagonal zeros
            maskLagged = np.ones((input.size, input.size))

    

        samples = sm.generate_timeseries(input.start, input.steps, input.ampContemp, matContemp, covContemp, input.ampLagged, matLagged, covLagged,measureCov)

        if input.clip_samples:
        # Checking the clipping
            samples_clip = sm.clip_timeseries(samples, input.clip_indices, input.clip_mins, input.clip_maxs)
        
        savefile = savepath+'ind_%i.txt'%(csvnum)
        np.savetxt(savefile,samples, delimiter=',')