import sim_modules as sm
import numpy as np
import matplotlib.pyplot  as plt
import sys 
import input as input 
import os

size= input.size
path= input.path
cols = 2*size

savepath=path+'sims/numParticipants_%s/mask_%s/'%(input.num_participants,str(input.maskZero))
try:
    os.makedirs(savepath)
except:
    print('directory %s exists'%savepath)


for number in range(input.num_participants)[0:2]:
    
    csvnum = number+1
    file = path+'%i.csv'%csvnum
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
    
    for j in range(input.num_iterations):
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
        
        savefile = savepath+'Person_%i_Series_t%s_covContemp%s_contempAmp%2.2f_covLagged%s_laggedAmp_%2.2f_measureCov%s_measureAmp%2.2f_iter%i.txt'%(csvnum,input.steps,input.covContempName,input.ampContemp,input.covLaggedName,input.ampLagged,input.measurecovName,input.ampMeasure,j)
        np.savetxt(savefile,samples, delimiter=',')