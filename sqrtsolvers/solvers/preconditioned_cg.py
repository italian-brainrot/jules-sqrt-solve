import numpy as np
import numpy as np
from .cg_lanczos import apply_A_sqrt
from .chebyshev import chebyshev_poly
from .lanczos import lanczos_tridiag

def solve_pcg_chebyshev(A, b, x0=None, max_iter=100, tol=1e-6, lanczos_k=30, chebyshev_k=10):
    """
    Solves A^(1/2)x = b using the preconditioned conjugate gradient method
    with a Chebyshev polynomial preconditioner.
    """
    if x0 is None:
        x = np.zeros_like(b)
    else:
        x = x0.copy()

    r = b - apply_A_sqrt(A, x, k=lanczos_k)

    # Preconditioner based on Chebyshev approximation of A^(-1/2)
    def preconditioner(r_vec):
        if np.linalg.norm(r_vec) < 1e-12:
            return np.zeros_like(r_vec)
        _, T = lanczos_tridiag(A, r_vec, k=lanczos_k, full_reorthogonalization=False)
        eigvals = np.linalg.eigvalsh(T)
        lambda_min = max(np.min(eigvals), 1e-12)
        lambda_max = np.max(eigvals)

        return chebyshev_poly(A, r_vec, lambda z: 1/np.sqrt(z), chebyshev_k, lambda_min, lambda_max)


    z = preconditioner(r)
    p = z.copy()
    rs_old = np.dot(r, z)

    if np.sqrt(np.dot(r, r)) < tol:
        return x

    for i in range(max_iter):
        Ap = apply_A_sqrt(A, p, k=lanczos_k)

        alpha = rs_old / np.dot(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap

        if np.sqrt(np.dot(r, r)) < tol:
            break

        z = preconditioner(r)
        rs_new = np.dot(r, z)
        p = z + (rs_new / rs_old) * p
        rs_old = rs_new

    return x
