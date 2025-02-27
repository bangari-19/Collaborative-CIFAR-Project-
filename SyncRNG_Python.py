import numpy as np
from SyncRNG import SyncRNG

# Set a universal seed at the top
# This seed is set in the beginning inside each controlled random function. It sets a deterministic algorithm
# This ensures that each function will produce the same output in both Python and R if it has the same inputs
universal_seed = 123
universal_rng = SyncRNG(seed=universal_seed)

def syncrng_box_muller(mu, sigma, n, rng=None):
    """
    Generate random samples from a normal distribution using the Box-Muller transform.

    Args:
        mu (float): Mean of the normal distribution.
        sigma (float): Standard deviation of the normal distribution.
        n (int): Number of random samples to generate.
        rng (SyncRNG): Instance of SyncRNG. If None, a universal_rng is used.

    Returns:
        np.ndarray: Array of n random samples from the specified normal distribution.
    """
    rng = rng or universal_rng
    two_pi = 2 * np.pi
    ngen = int(np.ceil(n / 2))
    out = np.zeros(2 * ngen)

    for i in range(ngen):
        u1 = 0.0
        u2 = 0.0
        while u1 == 0:
            u1 = rng.rand()
        while u2 == 0:
            u2 = rng.rand()

        mag = sigma * np.sqrt(-2.0 * np.log(u1))
        z0 = mag * np.cos(two_pi * u2) + mu
        z1 = mag * np.sin(two_pi * u2) + mu

        out[2 * i] = z0
        out[2 * i + 1] = z1

    return out[:n]

def generate_random_matrix(size, mu, sigma, rng=None):
    """
    Generate a square matrix with elements sampled from a normal distribution.

    Args:
        size (int): Dimension of the square matrix (number of rows and columns).
        mu (float): Mean of the normal distribution.
        sigma (float): Standard deviation of the normal distribution.
        rng (SyncRNG): Instance of SyncRNG. If None, a universal_rng is used.

    Returns:
        np.ndarray: size x size matrix with elements sampled from the specified normal distribution.
    """
    rng = rng or universal_rng
    random_numbers = syncrng_box_muller(mu, sigma, size * size, rng)
    return np.reshape(random_numbers, (size, size))

def generate_mvn_samples(mean_vector, cov_matrix, rng=None):
    """
    Generate a random vector from a multivariate normal distribution.

    Args:
        mean_vector (np.ndarray): Mean vector of the multivariate normal distribution.
        cov_matrix (np.ndarray): Covariance matrix of the multivariate normal distribution.
        rng (SyncRNG): Instance of SyncRNG. If None, a universal_rng is used.

    Returns:
        np.ndarray: Random vector sampled from the specified multivariate normal distribution.
    """
    rng = rng or universal_rng
    size = len(mean_vector)
    # Perform Cholesky decomposition
    L = np.linalg.cholesky(cov_matrix)  # Always returns lower triangular matrix
    Z = syncrng_box_muller(0, 1, size, rng)
    mvn_samples = np.dot(L, Z) + mean_vector
    return mvn_samples

# Check if output matches R 
if __name__ == "__main__":
    mu = 0
    sigma = 1
    size = 3
    mean_vector = np.array([1, 2, 3])
    cov_matrix = np.array([[1, 0.5, 0.3],
                           [0.5, 1, 0.4],
                           [0.3, 0.4, 1]])

    # Generate random numbers
    random_samples = syncrng_box_muller(mu, sigma, 10)
    print("Random Samples:", random_samples)

    # Generate a random matrix
    random_matrix = generate_random_matrix(size, mu, sigma)
    print("Random Matrix:\n", random_matrix)

    # Generate multivariate normal samples
    mvn_sample = generate_mvn_samples(mean_vector, cov_matrix)
    print("Multivariate Normal Sample:", mvn_sample)
