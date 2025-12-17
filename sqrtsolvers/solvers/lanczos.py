import numpy as np
from scipy.linalg import eigh_tridiagonal

def solve_lanczos(A, b, k=20):
    """
    Solves the system A^(1/2)x = b using the Lanczos algorithm to approximate A^(-1/2)b.
    This implementation uses full reorthogonalization to improve numerical stability.

    Args:
        A (np.ndarray): The matrix A.
        b (np.ndarray): The vector b.
        k (int, optional): The number of Lanczos iterations. Defaults to 20.

    Returns:
        np.ndarray: The solution vector x.
    """
    n = A.shape[0]
    q = b / np.linalg.norm(b)
    Q = np.zeros((n, k))
    alphas = np.zeros(k)
    betas = np.zeros(k - 1)

    for j in range(k):
        Q[:, j] = q
        v = A @ q
        alphas[j] = q @ v
        v = v - alphas[j] * q
        if j > 0:
            v = v - betas[j-1] * Q[:, j-1]

        # Reorthogonalization
        v = v - Q[:, :j+1] @ (Q[:, :j+1].T @ v)

        if j < k - 1:
            beta = np.linalg.norm(v)
            betas[j] = beta
            q = v / beta

    # Compute the eigendecomposition of T
    eigvals, eigvecs = eigh_tridiagonal(alphas, betas)

    # Approximate A^(-1/2)b
    # f(A)b approx Q f(T) Q^T b
    # A^(-1/2)b approx Q T^(-1/2) Q^T b
    # T^(-1/2) = V D^(-1/2) V^T
    # Q^T b = ||b|| * e_1

    T_inv_sqrt_b = eigvecs @ np.diag(1.0 / np.sqrt(np.maximum(eigvals, 1e-15))) @ eigvecs.T
    x = Q @ (T_inv_sqrt_b @ (np.linalg.norm(b) * np.eye(k, 1))).flatten()

    return x