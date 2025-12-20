import numpy as np
from scipy.linalg import eigh_tridiagonal

def lanczos_tridiag(A, v, k, full_reorthogonalization=True):
    """
    The Lanczos algorithm for symmetric matrices.
    """
    n = A.shape[0]
    k = min(k, n)
    Q = np.zeros((n, k))
    alphas = np.zeros(k)
    betas = np.zeros(k - 1)

    if np.linalg.norm(v) == 0:
        T = np.diag(alphas) + np.diag(betas, k=1) + np.diag(betas, k=-1)
        return np.zeros((n, k)), T

    q = v / np.linalg.norm(v)

    for j in range(k):
        Q[:, j] = q
        w = A @ q
        alphas[j] = np.dot(w, q)

        # Reorthogonalization
        if full_reorthogonalization:
            w = w - alphas[j] * q
            if j > 0:
                w = w - betas[j-1] * Q[:, j-1]
            w = w - Q[:, :j+1] @ (Q[:, :j+1].T @ w)
        else:
             w = w - alphas[j] * Q[:, j] - (betas[j-1] * Q[:, j-1] if j > 0 else 0)


        if j < k - 1:
            betas[j] = np.linalg.norm(w)
            if betas[j] < 1e-10:
                # Lanczos breakdown, terminate early
                k = j + 1
                alphas = alphas[:k]
                betas = betas[:k-1]
                Q = Q[:, :k]
                break
            q = w / betas[j]

    T = np.diag(alphas) + np.diag(betas, k=1) + np.diag(betas, k=-1)
    return Q, T


def solve_lanczos(A, b, k=10):
    """
    Solves A^(1/2)x = b using the Lanczos algorithm to approximate A^(-1/2)b.
    """
    Q, T = lanczos_tridiag(A, b, k)

    eigvals, eigvecs = eigh_tridiagonal(np.diag(T), np.diag(T, k=1))

    # Filter out small or negative eigenvalues
    eigvals[np.abs(eigvals) < 1e-12] = 1e-12
    eigvals[eigvals < 0] = 1e-12

    T_inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T

    # The vector b is aligned with the first Lanczos vector q_1 with a scaling factor
    k_actual = T.shape[0]
    e1 = np.zeros(k_actual)
    e1[0] = 1.0
    x = Q @ T_inv_sqrt @ e1 * np.linalg.norm(b)
    return x
