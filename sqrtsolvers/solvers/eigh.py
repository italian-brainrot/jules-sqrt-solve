import numpy as np

def solve_eigh(A, b):
    """
    Solves the system A^(1/2)x = b using direct eigendecomposition.

    Args:
        A (np.ndarray): The matrix A.
        b (np.ndarray): The vector b.

    Returns:
        np.ndarray: The solution vector x.
    """
    L, Q = np.linalg.eigh(A)
    A_sqrt_inv = Q @ np.diag(1.0 / np.sqrt(L)) @ Q.T
    return A_sqrt_inv @ b
