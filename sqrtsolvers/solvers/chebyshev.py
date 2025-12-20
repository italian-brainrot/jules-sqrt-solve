import numpy as np
from numpy.polynomial.chebyshev import Chebyshev

def chebyshev_poly(A, b, f, k, lambda_min, lambda_max):
    """
    Computes p(A)b where p is the best Chebyshev polynomial approximation
    of f(z) in the interval [lambda_min, lambda_max].
    """
    # Map for the operator
    alpha = (lambda_max + lambda_min) / 2.0
    beta = (lambda_max - lambda_min) / 2.0

    def mapped_A(v):
        return (A @ v - alpha * v) / beta

    # Get the Chebyshev polynomial coefficients
    pts = np.linspace(lambda_min, lambda_max, 100)
    cheb_poly = Chebyshev.fit(pts, f(pts), deg=k, domain=[lambda_min, lambda_max])
    c = cheb_poly.coef

    # Clenshaw's algorithm
    w0 = b
    w1 = mapped_A(w0)
    x = c[0] * w0 + c[1] * w1

    for i in range(2, k + 1):
        w2 = 2 * mapped_A(w1) - w0
        x += c[i] * w2
        w0, w1 = w1, w2

    return x

def solve_chebyshev(A, b, k=20, lanczos_k=10):
    """
    Solves the system A^(1/2)x = b using Chebyshev polynomial approximation.
    """
    from sqrtsolvers.solvers.lanczos import lanczos_tridiag
    _, T = lanczos_tridiag(A, b, k=lanczos_k)
    eigvals = np.linalg.eigvalsh(T)
    lambda_min = max(np.min(eigvals), 1e-12)
    lambda_max = np.max(eigvals)

    return chebyshev_poly(A, b, lambda z: 1/np.sqrt(z), k, lambda_min, lambda_max)
