#module file containing all bits needed in python
import numpy as np
from numpy.linalg import inv

import random
import os
from scipy.stats import norm

def plot_matrix(mat):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rc('image', cmap='cool')
    plt.figure()
    plt.imshow(mat)
    plt.colorbar()

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
    sample_step = np.matmul(cov,np.dot(amplitude,np.random.randn(size,size))) 
    # generate a sample of the coefficient matrix based on the structure. 
    # Should reproduce the R notebook perfectly if cov is uniform 1
    
    # make a draw of the covariance matrix
    stepmat = mat + sample_step  # additive rather than multiplicative to match the original simulations
    noisy_mat = np.multiply(mask,stepmat)# mask out the relevant terms
        
    return noisy_mat

def generate_timeseries(start,len,contempamp, contempmat, contempcov, lagamp, lagmat, lagcov, measurecov):
    ''' Using Y = Ylag*lagmat + Y*contempmat
    so Y= (Ylag*lagmat)*(I-contempmat)^-1
    '''
    import numpy as np
    size=np.shape(start)[0]
    #print(size,len)
    samples = np.zeros((len,size))
    samples[0,:] = start
    lagmask = make_mask(contempmat, contemp=False)
    contempmask = make_mask(lagmat, contemp=True)

    for i in range(1,len):
        
        noisy_lagmat = coeff_draw_from_cov(lagamp, lagmat,lagcov, lagmask)
        noisy_contempmat = coeff_draw_from_cov(contempamp, contempmat,contempcov, contempmask)
        samples[i,:] = np.matmul(np.matmul(samples[i-1,:], noisy_lagmat),np.linalg.inv(np.eye(size)-noisy_contempmat))+start
        samples[i,:]+=np.random.multivariate_normal(np.zeros(size),measurecov) # adding additional measurement noise
    
    samples-=start
    samples=samples[1:,:]
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