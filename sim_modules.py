#module file containing all bits needed in python
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

import random
random.seed(10)
import os
from scipy.stats import norm

plt.close()
def plot_matrix(mat, save=False, name=None):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc('image', cmap='cool')
    plt.figure()
    plt.imshow(mat)
    plt.colorbar()
    if save:
        plt.savefig(name)
    plt.close()

def save_matrix(mat, name):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np

    np.savetxt(name,mat, delimiter=',')
# separate out the contemporaneous diagonal masking

def make_mask(mat, contemp):
    size = np.shape(mat)[0]
    mask = np.ones((size,size))

    # Also only want to iterate over contemporaneous variables that are not set by themselves
    if contemp:
        np.fill_diagonal(mask,0)
     # Only want to iterate on entries that are non-zero
    
    inds = np.where(mat==0)
    mask[inds]=0
    
    return mask

def coeff_draw_from_cov(amplitude, mat,cov,mask):
    "draw a coefficient matrix from the mean matrix and the covariance"
    import numpy as np
    import matplotlib.pyplot as plt
    size=len(cov)

    plotmat=False
    # noRand=True # are we adding random noise to matrices
    
    # if noRand:
    #     ## If making the coefficient matrix the same as the one given by the sims team:
    #     sample_step = 0*np.matmul(cov,np.dot(amplitude,np.random.randn(size,size))) 
    # else:
    sample_step = np.matmul(cov,np.dot(amplitude,np.random.randn(size,size))) 

    # generate a sample of the coefficient matrix based on the structure. 
    
    # make a draw of the covariance matrix
    stepmat = mat + sample_step  # additive rather than multiplicative to match the original simulations
    #plot_matrix(sample_step,True, 'sample_stepmat%2.2f.png'%amplitude)
    
    if plotmat:
        plot_matrix(mat,True, 'origmat%2.2f.png'%amplitude)
        plot_matrix(stepmat,True, 'stepmat%2.2f.png'%amplitude)
    # plot_matrix(mask,True, 'mask%2.2f.png'%amplitude)
    # noisy_mat=np.zeros((size,size))
    # for i in range(size):
    #     for j in range(size):
    #         noisy_mat[i,j] = stepmat[i,j]*mask[i,j]
    noisy_mat = np.multiply(stepmat,mask)# mask out the relevant terms
    if plotmat:
        plot_matrix(noisy_mat,True, 'maskedmat%2.2f.png'%amplitude)
    return noisy_mat

def generate_timeseries(start,len,contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=False):
    ''' Using Y = Ylag*lagmat + Y*contempmat
    so Y= (Ylag*lagmat)*(I-contempmat)^-1
    '''
    import numpy as np
    size=np.shape(start)[0]
    #print(size,len)
    samples = np.zeros((len+100,size))
    samples[0,:] = start
    lagmask = make_mask(lagmat, contemp=False)
    contempmask = make_mask(contempmat, contemp=True)

    noisy_lagmat = coeff_draw_from_cov(lagamp, lagmat,lagcov, lagmask)
    noisy_contempmat = coeff_draw_from_cov(contempamp, contempmat,contempcov, contempmask)
    
    if save:
            save_matrix(noisy_lagmat, 'noisy_lagmat.csv')
            plot_matrix(noisy_lagmat, True, 'noisy_lagmat.png')
            save_matrix(noisy_contempmat, 'noisy_contempmat.csv')
            plot_matrix(noisy_contempmat, True, 'noisy_contempmat.png')

    #debug=True
    for i in range(1,len+100):
        samples[i,:] = np.matmul(np.matmul(samples[i-1,:], noisy_lagmat),np.linalg.inv(np.eye(size)-noisy_contempmat))+start
        samples[i,:]+= np.random.multivariate_normal(np.zeros(size),measurecov) # adding additional measurement noise
    
    samples-=start
    samples=samples[101:,:]
    return samples

def clip_timeseries(timeseries, indices_to_clip, min_vec, max_vec):

    ts = timeseries.copy()


    for i in indices_to_clip:
        min = min_vec[i]
        max = max_vec[i]
        min_inds = np.where(timeseries[:,i]<min)[0]
        ts[min_inds,i] = min
        max_inds = np.where(timeseries[:,i]>max)[0]
        ts[max_inds,i]=max

    return ts
    #clipping all to be within 0-5

    # indices = np.where(samples[:,:]<-10)[0]
    # samples[indices] = -10
    # indices = np.where(samples[:,:]>10)[0]
    # samples[indices] = 10
def clip_outliers(timeseries, sigma, measure_amp,debug,start, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov):

    #ts = timeseries.copy()
    ts = timeseries
    for i in range(np.shape(timeseries)[1]):
        
        min = -sigma*np.sqrt(measure_amp) #np.mean(timeseries[:,i])-sigma*np.sqrt(measure_amp)
        max = sigma*np.sqrt(measure_amp) #np.mean(timeseries[:,i])+sigma*np.sqrt(measure_amp)
        #print(min,max, i, sigma*np.sqrt(measure_amp),measure_amp)
        min_inds = np.where(timeseries[:,i]<min)[0]

        if debug:
            if len(min_inds)>0: print('bad vals min', timeseries[min_inds,i], 'inds', min_inds)
        countmax=2

        for minind in min_inds:
            count=0
            
            tmp=generate_timeseries(ts[minind-1,:],3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=False)[1]
            if debug:
                print('param', i, tmp[i],'tmp[i]','min= ', min,count, minind)

            while((ts[minind,i]< min) and count<countmax):
               
                tmp=generate_timeseries(ts[minind-1,:],3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=False)[1]

                if debug:
                    print(f" iteration {count} for param {i} and ind {minind} in min")
                    print(ts[minind-1:minind+2,i], 'param ',i, ' before')
                    #print(tmp[i],'tmp[i]',min, 'min after',count, countmax) # still a vector of 6 parameters
                

                ts[minind,i]=tmp[i]
                count+=1

                if debug:   
                    print(ts[minind-1:minind+2,i], 'after', count,countmax)
                    print('-----')

                if count==(countmax-1): 
                    ts[minind,i] = min*(1 - 0.01*np.random.rand()) #min here is less then zero
                    if debug:
                        print(ts[minind,i],'failing due to count')
                
            if debug:
                print('while loop complete')
                print(ts[minind,i],'new value for param', i)
                print('=======')

        max_inds = np.where(timeseries[:,i]>max)[0]
        
        if debug:
            if len(max_inds)>0: print('bad vals max', timeseries[max_inds,i], 'inds', max_inds)
        
        for maxind in max_inds:
            count=0
            
            tmp= generate_timeseries(ts[maxind-1,:],3, contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov, save=False)[1]
            
            if debug:
                print('param', i,tmp[i],'tmp[i]','max=',max, 'count=',count,maxind)

            while((ts[maxind,i] > max) and (count<countmax)):
                
                if debug:
                    print(f" iteration {count} for param {i} and ind {maxind} in max")
                    print(ts[maxind-1:maxind+2,i], 'ts before')
                    print(tmp[i], 'tmp[i]','max= ',max,  'count=', count)

                ts[maxind,i]= tmp[i]
                count+=1

                if debug:   
                    print(ts[maxind-1:maxind+2,i], 'after', count,countmax)
                    print('-----')


                if count==(countmax-1):
                    ts[maxind,i] = max*(1 - 0.01*np.random.rand())
                    if debug:
                        print(ts[maxind,i],'failing due to count')
                
            
            if debug:
                print('while loop complete')
                print(ts[maxind,i],'new value for param', i)
                print('=======')

        if debug:
            if (len(min_inds)==0 and len(max_inds)==0):
                print('no clipping needed for param %i'%i) 
            else:
                print('mean for param %i before clip: %2.4f, mean after clip %2.4f'%(i,np.mean(timeseries[:,i]),np.mean(ts[:,i])))
 
    return ts