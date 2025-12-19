import numpy as np

def inv_sqrt_ns(A, tol=1e-6, max_iter=100):
    """
    Computes the inverse square root of a matrix A using the Newton-Schulz method.

    Args:
        A (np.ndarray): The matrix A.
        tol (float, optional): The tolerance for convergence. Defaults to 1e-6.
        max_iter (int, optional): The maximum number of iterations. Defaults to 100.

    Returns:
        np.ndarray: The solution vector x.
    """
    n = A.shape[0]

    trace = np.trace(A)
    if trace > 0:
        # A better scaling factor based on the trace of A
        alpha = np.sqrt(trace / n)
        X = np.eye(n) / alpha

    else:
        # Fallback to the robust Frobenius norm
        frob_norm = np.sqrt(np.trace(A @ A.T))
        X = np.eye(n) / frob_norm

    I = np.eye(n)

    for _ in range(max_iter):
        X_old = X
        AX2 = A @ X @ X
        X = 0.5 * X @ (3 * I - AX2)
        if np.linalg.norm(X - X_old) < tol:
            break

    return X

def solve_ns(A, b, tol=1e-6, max_iter=100):
    """
    Solves the system A^(1/2)x = b using an iterative method
    based on the Newton-Schulz iteration for the inverse square root of A.

    Args:
        A (np.ndarray): The matrix A.
        b (np.ndarray): The vector b.
        tol (float, optional): The tolerance for convergence. Defaults to 1e-6.
        max_iter (int, optional): The maximum number of iterations. Defaults to 100.

    Returns:
        np.ndarray: The solution vector x.
    """
    A_inv_sqrt = inv_sqrt_ns(A, tol=tol, max_iter=max_iter)
    return A_inv_sqrt @ b
