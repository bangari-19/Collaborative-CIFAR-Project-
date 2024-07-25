#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 02:03:26 2024

@author: mac
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from SyncRNG import SyncRNG

# seed -> for single use
# rng -> for 'global' use, consistent rand. num. gen across multiple function calls
# set a universal seed or rng? will have to add 'rng' as parameter to all functions

universal_seed = 123
universal_rng = SyncRNG(seed=universal_seed)

def syncrng_box_muller(mu, sigma, n, seed=0, rng=None):
    """Generate n numbers from N(mu, sigma^2)"""
    rng = universal_rng if rng is None else rng

    two_pi = 2 * math.pi
    ngen = math.ceil(n / 2)
    out = [0.0] * 2 * ngen

    for i in range(ngen):
        u1 = 0.0
        u2 = 0.0

        while u1 == 0:
            u1 = rng.rand()
        while u2 == 0:
            u2 = rng.rand()

        mag = sigma * math.sqrt(-2.0 * math.log(u1))
        z0 = mag * math.cos(two_pi * u2) + mu
        z1 = mag * math.sin(two_pi * u2) + mu

        out[2*i] = z0
        out[2*i + 1] = z1

    return out[:n]

#plt.hist(syncrng_box_muller(mu=0, sigma=1, n=100), bins=10)

def generate_random_matrix(size, mu, sigma, rng=None):
    rng = universal_rng if rng is None else rng
    random_numbers = syncrng_box_muller(mu, sigma, size * size, rng=rng)
    return np.array(random_numbers).reshape(size, size)

def generate_mvn_samples(mean_vector, cov_matrix, rng=None):
    rng = universal_rng if rng is None else rng
    size = len(mean_vector)
    L = np.linalg.cholesky(cov_matrix).T
    Z = syncrng_box_muller(0, 1, size, rng=rng)
    mvn_samples = np.dot(L, Z) + mean_vector
    return mvn_samples

############################ EXAMPLES ############################
mean_vector = np.array([1, 3])
cov_matrix = np.array([[2, 0], [0, 2]])
n_samples = 1000
n_iterations = 3

mvn_sample = generate_mvn_samples(mean_vector, cov_matrix, rng=None)
print(mvn_sample)

random_sample = generate_random_matrix(2, 0, 1)
print(random_sample)

'''
# PLOT 1
for i in range(n_iterations):
    # Generate samples using generate_mvn_samples
    samples_custom = np.array([generate_mvn_samples(mean_vector, cov_matrix) for _ in range(n_samples)])

    # Generate samples using np.random.multivariate_normal
    samples_np = np.random.multivariate_normal(mean_vector, cov_matrix, n_samples)
    
    # Plot
    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot syncrng samples
    ax[0].scatter(samples_custom[:, 0], samples_custom[:, 1], alpha=0.5, label='Custom MVN', color='blue')
    ax[0].set_title(f'Custom MVN Samples - Iteration {i+1}')
    ax[0].set_xlabel('X')
    ax[0].set_ylabel('Y')
    ax[0].legend()
    
    # Plot numpy samples
    ax[1].scatter(samples_np[:, 0], samples_np[:, 1], alpha=0.5, label='NumPy MVN', color='orange')
    ax[1].set_title(f'NumPy MVN Samples - Iteration {i+1}')
    ax[1].set_xlabel('X')
    ax[1].set_ylabel('Y')
    ax[1].legend()
    
    # Plot both overlapped
    ax[2].scatter(samples_custom[:, 0], samples_custom[:, 1], alpha=0.5, label='Custom MVN', color='blue')
    ax[2].scatter(samples_np[:, 0], samples_np[:, 1], alpha=0.5, label='NumPy MVN', color='orange')
    ax[2].set_title(f'Overlapped MVN Samples - Iteration {i+1}')
    ax[2].set_xlabel('X')
    ax[2].set_ylabel('Y')
    ax[2].legend()
    
    plt.show()

'''
'''
# PLOT 2
fig, axes = plt.subplots(n_iterations, 3, figsize=(18, 6 * n_iterations))

for i in range(n_iterations):
    # Generate samples using generate_mvn_samples
    samples_custom = np.array([generate_mvn_samples(mean_vector, cov_matrix) for _ in range(n_samples)])

    # Generate samples using np.random.multivariate_normal
    samples_np = np.random.multivariate_normal(mean_vector, cov_matrix, n_samples)
    
    # Plot syncrng samples
    ax_custom = axes[i, 0]
    ax_custom.scatter(samples_custom[:, 0], samples_custom[:, 1], alpha=0.5, label='Custom MVN', color='blue')
    ax_custom.set_title(f'Custom MVN Samples - Iteration {i+1}')
    ax_custom.set_xlabel('X')
    ax_custom.set_ylabel('Y')
    ax_custom.legend()
    
    # Plot numpy samples
    ax_np = axes[i, 1]
    ax_np.scatter(samples_np[:, 0], samples_np[:, 1], alpha=0.5, label='NumPy MVN', color='orange')
    ax_np.set_title(f'NumPy MVN Samples - Iteration {i+1}')
    ax_np.set_xlabel('X')
    ax_np.set_ylabel('Y')
    ax_np.legend()
    
    # Plot both overlapped
    ax_overlap = axes[i, 2]
    ax_overlap.scatter(samples_custom[:, 0], samples_custom[:, 1], alpha=0.5, label='Custom MVN', color='blue')
    ax_overlap.scatter(samples_np[:, 0], samples_np[:, 1], alpha=0.5, label='NumPy MVN', color='orange')
    ax_overlap.set_title(f'Overlapped MVN Samples - Iteration {i+1}')
    ax_overlap.set_xlabel('X')
    ax_overlap.set_ylabel('Y')
    ax_overlap.legend()

plt.tight_layout()
plt.show()
'''
