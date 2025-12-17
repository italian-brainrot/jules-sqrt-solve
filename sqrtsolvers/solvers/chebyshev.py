import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from scipy.linalg import eigh_tridiagonal

def solve_chebyshev(A, b, k=20, lanczos_k=10):
    """
    Solves the system A^(1/2)x = b using Chebyshev polynomial approximation.

    Args:
        A (np.ndarray): The matrix A.
        b (np.ndarray): The vector b.
        k (int, optional): The degree of the Chebyshev polynomial. Defaults to 20.
        lanczos_k (int, optional): The number of Lanczos iterations to estimate the spectral range.

    Returns:
        np.ndarray: The solution vector x.
    """
    n = A.shape[0]

    # 1. Estimate the spectral range of A using Lanczos
    q = b / np.linalg.norm(b)
    Q = np.zeros((n, lanczos_k))
    alphas = np.zeros(lanczos_k)
    betas = np.zeros(lanczos_k - 1)

    for j in range(lanczos_k):
        Q[:, j] = q
        v = A @ q
        alphas[j] = q @ v
        v = v - alphas[j] * q
        if j > 0:
            v = v - betas[j-1] * Q[:, j-1]

        if j < lanczos_k - 1:
            beta = np.linalg.norm(v)
            betas[j] = beta
            if beta < 1e-12:
                break
            q = v / beta

    eigvals, _ = eigh_tridiagonal(alphas, betas)
    lambda_min = np.min(eigvals)
    lambda_max = np.max(eigvals)

    # 2. Find the best polynomial approximation of z^(-1/2)
    # The domain for the Chebyshev fit needs to be [-1, 1]
    # We map [lambda_min, lambda_max] to [-1, 1]
    # Let z = (lambda_max - lambda_min)/2 * t + (lambda_max + lambda_min)/2
    # Then t = (2z - lambda_max - lambda_min) / (lambda_max - lambda_min)

    # We want to approximate f(z) = 1/sqrt(z)

    # We create a set of points in the spectral range to fit the polynomial
    pts = np.linspace(lambda_min, lambda_max, 100)

    # The function to approximate
    f = lambda z: 1.0 / np.sqrt(z)

    # Get the Chebyshev polynomial
    cheb_poly = Chebyshev.fit(pts, f(pts), deg=k, domain=[lambda_min, lambda_max])

    # 3. Use Clenshaw algorithm to compute p(A)b
    # The Clenshaw algorithm for p(A)b is given by:
    # y_0 = c_k * b
    # y_1 = c_{k-1} * b + 2 * (A - alpha*I) / beta * y_0
    # y_j = c_{k-j} * b + 2 * (A - alpha*I) / beta * y_{j-1} - y_{j-2} for j=2,...,k
    # x = y_k

    # The Chebyshev recurrence is T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)
    # Our polynomial is p(z) = sum_{i=0 to k} c_i T_i(t), where t is the mapped variable

    # We can use a simpler three-term recurrence for applying the Chebyshev polynomial
    # to the vector b

    # Map for the operator
    alpha = (lambda_max + lambda_min) / 2.0
    beta = (lambda_max - lambda_min) / 2.0

    def mapped_A(v):
        return (A @ v - alpha * v) / beta

    c = cheb_poly.coef

    w0 = b
    w1 = mapped_A(w0)
    x = c[0] * w0 + c[1] * w1

    for i in range(2, k + 1):
        w2 = 2 * mapped_A(w1) - w0
        x += c[i] * w2
        w0, w1 = w1, w2

    return x
