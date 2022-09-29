import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import random

params = {'figure.titlesize': 15,
          'axes.labelsize': 10,
          'axes.titlesize': 12,
          'axes.linewidth': 1.5,
          'xtick.minor.visible': True,
          'xtick.major.width': 1.5,
          'xtick.direction': 'in',
          'ytick.direction': 'in',
          'xtick.major.size': 7,
          'xtick.minor.size': 5,
          'xtick.labelsize': 10,
          'ytick.minor.visible': True,
          'ytick.major.width': 1.5,
          'ytick.major.size': 7,
          'ytick.minor.size': 5,
          'ytick.labelsize': 10,
          'legend.fontsize': 10,
          'legend.labelspacing': 0.3,
          'font.family': 'sans-serif'}

plt.rcParams.update(params)

path = '/Users/mac/Desktop/Astro/CIFAR Project'

def get_data(ind: str):
    '''
    Load data for variables VR, SS, ANX, INT for individual ind
    
    Parameters
    ----------
    ind : str
        Individual number (01 - 13)
        
    Returns
    -------
    VR: array 
        Verbal recall (0 - 5)
    SS: array 
        Somatic sympton (1 - 5)
    ANX: array 
        Anxiety (1 - 5)
    INT: array 
        Intellectual interest (1 - 5)
    '''
    path_data = path+'/GIMME analyses/Cleaned GIMME Data/Control_Group/10' 
    VR, SS, ANX, INT = np.genfromtxt(path_data+ind+'.txt', skip_header=1, unpack=True) 
    return VR, SS, ANX, INT

def get_betas(ind: str):
    '''
    Load beta weights for individual ind
    
    Parameters
    ----------
    ind : str
        Individual number (01 - 13)
        
    Returns
    -------
    A: array 
        Contemporaneous coefficient matrix
    B: array 
        Lagged coefficient matrix
    '''
    path_betas = path+'/GIMME analyses/Cleaned GIMME Data/Control_Group_output/individual/10'
    betas = np.loadtxt(path_betas+ind+'Betas.csv', skiprows=1, usecols=range(1,9), delimiter=',')
    A = betas[:,4:]  # Contemporaneous coeffs
    B = betas[:,:4]  # Lagged coeffs
    return A, B

def covariance_matrix(ind: str, cov_type: str, alpha: float):
    '''
    Construct covariance matrix specified by cov_type and scaling alpha for individual ind.
    
    Parameters
    ----------
    ind: str
        Individual number (01 - 13)
    cov_type: str
        cov_type can be 'diagonal', 'AR1' or 'compound'
            - AR1 covariance structure is a first-order autoregressive structure with heterogenous variances
            - Compound covariance structure has heterogenous variances and constant correlation between elements
            - Diagonal covariance structure has heterogenous variances and zero correlation between elements 
    alpha: float
        Scaling parameters between elements 
    
    Returns
    -------
    diagonal : array 
        Diagonal covariance structure 
    AR1_cov : array 
        AR(1) heterogenous covariance structure
    cov : array 
        Compound symmetry: heterogenous covariance structure 
    '''
    # Load data for individual, calculate variance and standard deviation 
    VR, SS, ANX, INT = get_data(ind)
    var_vr, std_vr = np.nanvar(VR), np.nanstd(VR)  
    var_ss, std_ss = np.nanvar(SS), np.nanstd(SS)
    var_anx, std_anx = np.nanvar(ANX), np.nanstd(ANX)
    var_int, std_int = np.nanvar(INT), np.nanstd(INT)
    
    # Diagonal 
    diagonal = alpha**2*np.diag([var_vr, var_ss, var_anx, var_int])
    
    # AR1 
    AR1_cov = np.zeros((4,4))
    AR1_cov[0, 0], AR1_cov[1, 1], AR1_cov[2, 2], AR1_cov[3, 3] = var_vr, var_ss, var_anx, var_int  # Diagonal 0
    AR1_cov[0, 1], AR1_cov[1, 2], AR1_cov[2, 3] = alpha*np.array([std_vr*std_ss, std_ss*std_anx, std_anx*std_int])  # Diagonal 1
    AR1_cov[1, 0], AR1_cov[2, 1], AR1_cov[3, 2] = AR1_cov[0, 1], AR1_cov[1, 2], AR1_cov[2, 3]
    AR1_cov[0, 2], AR1_cov[1, 3] = alpha**2*np.array([std_vr*std_anx, std_ss*std_int])  # Diagonal 2
    AR1_cov[2, 0], AR1_cov[3, 1] = AR1_cov[0, 2], AR1_cov[1, 3]
    AR1_cov[0, 3] = AR1_cov[3, 0] = alpha**3*(std_vr*std_int)  # Diagonal 3
    
    # Compound 
    cov = np.zeros((4,4))
    cov[0, 0], cov[1, 1], cov[2, 2], cov[3, 3] = var_vr, var_ss, var_anx, var_int  # Diagonal 0
    cov[0, 1], cov[1, 2], cov[2, 3] = alpha*np.array([std_vr*std_ss, std_ss*std_anx, std_anx*std_int])  # Diagonal 1
    cov[1, 0], cov[2, 1], cov[3, 2] = cov[0, 1], cov[1, 2], cov[2, 3]
    cov[0, 2], cov[1, 3] = alpha*np.array([std_vr*std_anx, std_ss*std_int])  # Diagonal 2
    cov[2, 0], cov[3, 1] = cov[0, 2], cov[1, 3]
    cov[0, 3] = cov[3, 0] = alpha*(std_vr*std_int)  # Diagonal 3
    
    if cov_type == 'diagonal': 
        return diagonal
    if cov_type == 'AR1':
        return AR1_cov 
    if cov_type == 'compound': 
        return cov
    else: 
        print('Error: cov_type should be one of diagonal, AR1, or compound')
    
def get_matrix_multivar(ind: str, shape: tuple, covariance):
    '''
    Generate time series data for individual ind and store into matrix M  
    Data is simulated according to $M[r] = (I - A)^{-1}B + (I - A)^{-1}C$
        I is the identity matrix
        A is the contemperaneous coefficients matrix 
        B is the lagged coefficients matrix
        C is the noise matrix generated using the np.random.multivariate_normal method
    VR data is clipped between 0 - 5 and all other variable data 1 - 5
    
    Parameters
    ----------
    ind: str
        Individual number (01 - 13)
    shape: tuple
        Shape = (number of days, number of variables)
    covariance: array
        Covariance structure used to generate error 
    
    Returns
    -------
    M: array 
        Simulated data for variables VR, SS, ANX, INT
        Rows represent days, columns represent variables 
    '''
    # Load data for individual, calculate mean 
    # Load contemperaneous and lagged beta weights 
    VR, SS, ANX, INT = get_data(ind)
    contemp_coeff, lagged_coeff = get_betas(ind)
    mean_vr = np.nanmean(VR)  
    mean_ss = np.nanmean(SS)
    mean_anx = np.nanmean(ANX)
    mean_int = np.nanmean(INT)
    
    M = np.zeros(shape)
    noise = np.random.multivariate_normal(np.zeros(4), covariance)  # Using zero mean and covariance structure 
    M[0,:] = noise
 
    for r in range(1, shape[0]): 
        noise = np.random.multivariate_normal(np.zeros(4), covariance)
        M[r] = np.matmul(np.matmul(np.linalg.inv(np.eye(4)-contemp_coeff),lagged_coeff),M[r-1]) 
        M[r] += np.matmul(np.linalg.inv(np.eye(4)-contemp_coeff),noise)
    M[:,:] += np.array([mean_vr, mean_ss, mean_anx, mean_int])
    
    # Clipping VR (0-5)
    indices = np.where(M[:,0]<0)[0]
    M[indices,0] = 0
    indices = np.where(M[:,0]>5)[0]
    M[indices,0] = 5

    # Clipping other variables (1-5)
    indices = np.where(M[:,1:4]<1)[0]
    M[indices,1:4] = 1
    indices = np.where(M[:,1:4]>5)[0]
    M[indices,1:4] = 5

    return M

def get_matrix_vector(ind: str, shape: tuple, contemp_coeff, lagged_coeff, alpha: float):
    '''
    Generate time series data for individual ind and store into matrix M  
    Data is simulated according to $M[r] = (I - A)^{-1}B + (I - A)^{-1}C$
        I is the identity matrix
        A is the contemperaneous coefficients matrix 
        B is the lagged coefficients matrix
        C is the noise matrix generated using the np.random.randn method and scaled by alpha
    VR data is clipped between 0 - 5 and all other variable data 1 - 5
    
    Parameters
    ----------
    ind: str
        Individual number (01 - 13)
    shape: tuple
        Shape = (number of days, number of variables)
    alpha: float
        Scaling parameter 
    
    Returns
    -------
    M: array 
        Simulated data for variables VR, SS, ANX, INT
        Rows represent days, columns represent variables 
    '''
    VR, SS, ANX, INT = get_data(ind)
    std_vr, mean_vr = np.nanstd(VR), np.nanmean(VR)  
    std_ss, mean_ss = np.nanstd(SS), np.nanmean(SS)
    std_anx, mean_anx = np.nanstd(ANX), np.nanmean(ANX)
    std_int, mean_int = np.nanstd(INT), np.nanmean(INT)
    
    M = np.zeros(shape)
    noise = alpha*np.random.randn(4)*np.array([std_vr, std_ss, std_anx, std_int])
    M[0,:] = noise

    for r in range(1, shape[0]): 
        noise = alpha*np.random.randn(4)*np.array([std_vr, std_ss, std_anx, std_int])
        M[r] = np.matmul(np.matmul(np.linalg.inv(np.eye(4)-contemp_coeff),lagged_coeff),M[r-1])     
        M[r] += np.matmul(np.linalg.inv(np.eye(4)-contemp_coeff),noise)
    M[:,:] += np.array([mean_vr, mean_ss, mean_anx, mean_int])
    
    indices = np.where(M[:,0]<0)[0]
    M[indices,0] = 0
    indices = np.where(M[:,0]>5)[0]
    M[indices,0] = 5

    indices = np.where(M[:,1:4]<1)[0]
    M[indices,1:4] = 1
    indices = np.where(M[:,1:4]>5)[0]
    M[indices,1:4] = 5

    return M