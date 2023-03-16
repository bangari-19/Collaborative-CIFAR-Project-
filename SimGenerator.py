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
print('Sim Directory is %s'%input.pathin)
print('reading from ', input.num_participants, input.steps, input.ampContemp, input.ampLagged, input.ampMeasure, input.maskZero)
printpath=pathout+'sims/numParticipants_%s'%input.num_participants
print('Main dir is %s'%printpath)

for j in range(input.num_iterations):
    savepath=pathout+'sims/numParticipants_%s/steps_%s/Contemp_%s_amp_%s/Lagged_%s_amp_%s/Measure_%s_amp_%s/mask_%s/rep_%i/'%(input.num_participants,str(input.steps),input.covContempName,str(input.ampContemp),input.covLaggedName, str(input.ampLagged),input.measurecovName, str(input.ampMeasure),str(input.maskZero),j)
    try:
        os.makedirs(savepath)
        if input.debug:
            print('making %s'%savepath)
    except:
        if input.debug:
            print('directory %s exists'%savepath)
        pass

    for number in range(input.num_participants):
        
        csvnum = number+1
        if input.debug:
            print('Writing sims for participant %i and iteration %i'%(csvnum,j))
        file = pathin+'%i.csv'%csvnum
        data = np.loadtxt(file, skiprows=1, usecols=range(0,cols), delimiter=',')
        matContemp = data[:,size:] #same day (contemporaneous) beta values
        matLagged = data[:,:size] #lagged beta values
        
        # sm.save_matrix(matContemp, 'matContempOrig.csv')
        # sm.save_matrix(matLagged, 'matLaggedOrig.csv')
        # sm.plot_matrix(matLagged, True, 'matLaggedOrig.png')
        # sm.plot_matrix(matContemp, True, 'matContempOrig.png')
        
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

        # sm.plot_matrix(maskLagged, True, 'maskLagged.png')
        # sm.plot_matrix(maskContemp, True, 'maskContemp.png')

        samples = sm.generate_timeseries(input.start, input.steps, input.ampContemp, matContemp, covContemp, input.ampLagged, matLagged, covLagged,measureCov,input.save)

        if input.clip_samples:
        # Checking the clipping
            samples_clip = sm.clip_timeseries(samples, input.clip_indices, input.clip_mins, input.clip_maxs)
            samples=samples_clip
            if input.debug:
                print('saved clip')
        
        if input.clip_outliers:
            if input.debug:
                print('Clipping samples')
            samples_clip = sm.clip_outliers(samples, input.clip_sigma,input.ampMeasure, input.debug)
            samples=samples_clip
            if input.debug:
                print('saved clip')

        savefile = savepath+'ind_%i.txt'%(csvnum)
        np.savetxt(savefile,samples, delimiter=',')