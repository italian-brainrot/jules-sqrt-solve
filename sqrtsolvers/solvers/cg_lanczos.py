import numpy as np

def lanczos_algorithm(A, v, k):
    """
    The Lanczos algorithm for symmetric matrices with full reorthogonalization.
    """
    n = A.shape[0]
    k = min(k, n)
    Q = np.zeros((n, k + 1))
    alphas = np.zeros(k)
    betas = np.zeros(k)

    if np.linalg.norm(v) == 0:
        return np.zeros((n, k)), np.zeros((k, k))

    q = v / np.linalg.norm(v)
    Q[:, 0] = q

    for j in range(k):
        w = A @ Q[:, j]
        alphas[j] = np.dot(w, Q[:, j])
        w = w - alphas[j] * Q[:, j]
        if j > 0:
            w = w - betas[j-1] * Q[:, j-1]

        for i in range(j + 1):
            w = w - np.dot(w, Q[:, i]) * Q[:, i]

        betas[j] = np.linalg.norm(w)

        if betas[j] < 1e-10:
            k = j + 1
            break

        Q[:, j+1] = w / betas[j]

    T = np.diag(alphas[:k]) + np.diag(betas[:k-1], k=1) + np.diag(betas[:k-1], k=-1)
    return Q[:, :k], T


def apply_A_sqrt(A, p, k=10):
    """
    Computes A^(1/2)p using the Lanczos algorithm.
    """
    if np.linalg.norm(p) < 1e-10:
        return np.zeros_like(p)

    Q, T = lanczos_algorithm(A, p, k)

    eigvals, eigvecs = np.linalg.eigh(T)
    eigvals = np.maximum(eigvals, 1e-12)
    T_sqrt = eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T

    k_actual = T.shape[0]
    p_in_Q_basis = np.zeros(k_actual)
    p_in_Q_basis[0] = np.linalg.norm(p)

    return Q @ T_sqrt @ p_in_Q_basis

def solve_cg_lanczos(A, b, max_iter=100, tol=1e-6, lanczos_k=10):
    """
    Solves A^(1/2)x = b using the conjugate gradient method, where the
    matrix-vector products with A^(1/2) are approximated using the
    Lanczos algorithm.
    """
    def matvec(p):
        return apply_A_sqrt(A, p, k=lanczos_k)

    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rs_old = np.dot(r, r)

    if np.sqrt(rs_old) < tol:
        return x

    for i in range(max_iter):
        Ap = matvec(p)
        alpha = rs_old / np.dot(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rs_new = np.dot(r, r)

        if np.sqrt(rs_new) < tol:
            break

        p = r + (rs_new / rs_old) * p
        rs_old = rs_new

    return x
