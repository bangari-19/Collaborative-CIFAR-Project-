import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import random

path = '/Users/mac/Desktop/Astro/CIFAR Project'
mean1, std1 = 2.4, 1.19 #mean and standard dev for variable 1 across all 100 individuals 
mean2, std2 = 3.02, 1.31
mean3, std3 = 5.01, 1.38
mean4, std4 = 4.05, 1.19
mean5, std5 = 3.21, 1.20
mean6, std6 = 4.23, 1.17

var1, var2, var3, var4, var5, var6 = np.array([std1, std2, std3, std4, std5, std6])**2

#functions used to generate simulated data for six variables and 100 individuals

#beta coefficients 
def get_betas(ind: str):
    path_betas = path+'/220429.FinalMatrices/'
    betas = np.loadtxt(path_betas+ind+'.csv', skiprows=1, usecols=range(12), delimiter=',')
    A = betas[:,6:] #contemporaneous coeffs
    B = betas[:,:6] #lagged coeffs
    return A, B

#covariance structures 
def get_cov(ind, alpha, cov_type):
    "Return the 6x6 covariance structure specified by cov_type. cov_type can be 'diagonal', 'AR1' or 'compound'."
    
    #Diagonal
    diagonal = alpha**2*np.diag([var1, var2, var3, var4, var5, var6])
    
    #AR1
    AR1_cov = np.zeros((6,6))
    
    #diagonal 0
    AR1_cov[0,0], AR1_cov[1,1], AR1_cov[2,2], AR1_cov[3,3], AR1_cov[4,4], AR1_cov[5,5] =  var1, var2, var3, var4, var5, var6 

    #diagonal 1
    AR1_cov[0,1], AR1_cov[1,2], AR1_cov[2,3], AR1_cov[3,4], AR1_cov[4,5] = alpha*np.array([std1*std2, std2*std3, std3*std4, std4*std5, std5*std6])
    AR1_cov[1,0], AR1_cov[2,1], AR1_cov[3,2], AR1_cov[4,3], AR1_cov[5,4] = AR1_cov[0,1], AR1_cov[1,2], AR1_cov[2,3], AR1_cov[3,4], AR1_cov[4,5]
    
    #diagonal 2
    AR1_cov[0,2], AR1_cov[1,3], AR1_cov[2,4], AR1_cov[3,5] = alpha**2*np.array([std1*std3, std2*std4, std3*std5, std4*std6])
    AR1_cov[2,0], AR1_cov[3,1], AR1_cov[4,2], AR1_cov[5,3] = AR1_cov[0,2], AR1_cov[1,3], AR1_cov[2,4], AR1_cov[3,5]
    
    #diagonal 3
    AR1_cov[0,3], AR1_cov[1,4], AR1_cov[2,5] = alpha**3*np.array([std1*std4, std2*std5, std3*std6])
    AR1_cov[3,0], AR1_cov[4,1], AR1_cov[5,2] = AR1_cov[0,3], AR1_cov[1,4], AR1_cov[2,5] 
    
    #diagonal 4
    AR1_cov[0,4], AR1_cov[1,5] = alpha**4*np.array([std1*std5, std2*std6])
    AR1_cov[4,0], AR1_cov[5,1] = AR1_cov[0,4], AR1_cov[1,5] 
    
    #diagonal 5
    AR1_cov[0,5] = AR1_cov[5,0] = alpha**5*(std1*std6)
    
    
    #Compound 
    cov = np.zeros((6,6))
    #diagonal 0
    cov[0, 0], cov[1, 1], cov[2, 2], cov[3, 3], cov[4, 4], cov[5, 5] =  var1, var2, var3, var4, var5, var6 

    #diagonal 1
    cov[0, 1], cov[1, 2], cov[2, 3], cov[3, 4], cov[4, 5] = alpha*np.array([std1*std2, std2*std3, std3*std4, std4*std5, std5*std6])
    cov[1,0], cov[2,1], cov[3,2], cov[4,3], cov[5,4] = cov[0,1], cov[1,2], cov[2,3], cov[3,4], cov[4,5] 
    
    #diagonal 2
    cov[0, 2], cov[1, 3], cov[2, 4], cov[3, 5] = alpha*np.array([std1*std3, std2*std4, std3*std5, std4*std6])
    cov[2,0], cov[3,1], cov[4,2], cov[5,3] = cov[0, 2], cov[1, 3], cov[2, 4], cov[3, 5]
    
    #diagonal 3
    cov[0, 3], cov[1, 4], cov[2, 5] = alpha*np.array([std1*std4, std2*std5, std3*std6])
    cov[3,0], cov[4,1], cov[5,2] = cov[0, 3], cov[1, 4], cov[2, 5] 
    
    #diagonal 4
    cov[0, 4], cov[1, 5] = alpha*np.array([std1*std5, std2*std6])
    cov[4,0], cov[5,1] = cov[0, 4], cov[1, 5]
    
    #diagonal 5
    cov[0, 5] = cov[5, 0] = alpha*(std1*std6)
    
    if cov_type == 'diagonal': 
        return diagonal
    if cov_type == 'AR1':
        return AR1_cov 
    if cov_type == 'compound': 
        return cov
    else: 
        print('cov_type should be one of the following: diagonal, AR1, compound')
        
#Simulated data using multivariate_normal method to generate noise
def get_matrix(ind, shape, covariance):
    """Return simulated data for variables in matrix of shape = (# of days, # of variables) for 
    individual ind with covariance structure given by covariance. Values less than zero are 
    clipped. """
    
    contemp_coeff, lagged_coeff = get_betas(ind)
    
    M = np.zeros(shape)
    noise = np.random.multivariate_normal(np.zeros(6), covariance)
    M[0,:] = noise
 
    for r in range(1, shape[0]): 
        noise = np.random.multivariate_normal(np.zeros(6), covariance)
        M[r] = np.matmul(np.matmul(np.linalg.inv(np.eye(6)-contemp_coeff),lagged_coeff),M[r-1]) 
        M[r] += np.matmul(np.linalg.inv(np.eye(6)-contemp_coeff),noise)
    M[:,:] += np.array([mean1, mean2, mean3, mean4, mean5, mean6])
    
    indices = np.where(M[:,:]<0)[0]
    M[indices,:] = 0
    
    return M

# def plot(num, alpha, cov_type):
#     """Plot data for a given alpha and covariance type."""
#     VR1, VR2, VR3, VR4, VR5, VR6 = np.loadtxt()
#     time = np.arange(0, days) 
#     plt.plot(time, VR1)
#     plt.plot(time, VR2)