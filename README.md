# Novel Solver for A^(1/2)x = b

This repository contains the implementation and benchmarking of solvers for the matrix equation A^(1/2)x = b, where A is a symmetric positive definite matrix.

## Solvers

Currently four solvers are implemented in this repository:

1.  **Eigendecomposition**: This method computes the eigendecomposition of the matrix A to find its square root and then solves the system directly. It is accurate and robust but can be computationally expensive for large matrices.

2.  **Lanczos**: This is an iterative method for finding eigenvalues and eigenvectors of a symmetric matrix. In this repository, the Lanczos algorithm is used to approximate A^(-1/2)b.

3.  **Newton-Schulz**: This is an iterative solver based on the Newton-Schulz iteration for finding the inverse square root of the matrix A. The iteration is given by:
    X_{k+1} = 0.5 * X_k * (3I - A * X_k^2)

4.  **Chebyshev**: This method approximates the function f(z) = z^(-1/2) with a Chebyshev polynomial, and then computes the solution as a matrix-polynomial-vector product. It can be more stable than Newton-Schulz for ill-conditioned matrices.

## Benchmarking

The solvers are benchmarked on a range of symmetric positive definite matrices with varying sizes and condition numbers. The following metrics are used for comparison:

*   **Execution Time**: The wall-clock time taken to solve the system.
*   **Relative Error**: The 2-norm of the difference between the computed solution and the true solution, divided by the 2-norm of the true solution.

A warm-start scenario is also benchmarked if it is supported by the solver. We obtain a new matrix equation by adding low rank correction to A and perturbing b. We then use the solver to solve this equation, and use the solution as a warm start.

## Results

The benchmark results are summarized in the table below:

```
--------------------------------------------------------------------------------
Solver                         | Size       | Cond       | Time (s)        | Error          
--------------------------------------------------------------------------------
Eigendecomposition             | 10         | 1e+01      | 0.000151        | 1.611221e-15
Lanczos (k=5)                  | 10         | 1e+01      | 0.000506        | 9.869019e-03
Lanczos (k=50)                 | 10         | 1e+01      | 0.001248        | 1.622088e-15
Newton-Schulz (k=5)            | 10         | 1e+01      | 0.000145        | 2.106957e-05
Newton-Schulz (k=50)           | 10         | 1e+01      | 0.000137        | 3.177023e-16
Chebyshev (k=5)                | 10         | 1e+01      | 0.000762        | 1.184544e-02
Chebyshev (k=50)               | 10         | 1e+01      | 0.002025        | 2.045273e-15
Eigendecomposition             | 10         | 1e+04      | 0.000056        | 3.975556e-13
Lanczos (k=5)                  | 10         | 1e+04      | 0.000202        | 3.193981e-01
Lanczos (k=50)                 | 10         | 1e+04      | 0.001082        | 2.407978e-13
Newton-Schulz (k=5)            | 10         | 1e+04      | 0.000115        | 2.976397e-01
Newton-Schulz (k=50)           | 10         | 1e+04      | 0.001994        | nan
Chebyshev (k=5)                | 10         | 1e+04      | 0.000592        | 1.693419e+00
Chebyshev (k=50)               | 10         | 1e+04      | 0.001848        | 4.761975e-06
Eigendecomposition             | 10         | 1e+08      | 0.000058        | 1.238661e-10
Lanczos (k=5)                  | 10         | 1e+08      | 0.000239        | 1.043984e-01
Lanczos (k=50)                 | 10         | 1e+08      | 0.001086        | 1.991683e-09
Newton-Schulz (k=5)            | 10         | 1e+08      | 0.000116        | 1.039520e-01
Newton-Schulz (k=50)           | 10         | 1e+08      | 0.000921        | nan
Chebyshev (k=5)                | 10         | 1e+08      | 0.000556        | 1.796566e+02
Chebyshev (k=50)               | 10         | 1e+08      | 0.001862        | 4.239729e-04
Eigendecomposition             | 100        | 1e+01      | 0.001677        | 2.068418e-15
Lanczos (k=5)                  | 100        | 1e+01      | 0.000214        | 9.613844e-03
Lanczos (k=50)                 | 100        | 1e+01      | 0.001248        | 1.210597e-15
Newton-Schulz (k=5)            | 100        | 1e+01      | 0.001344        | 2.511794e-06
Newton-Schulz (k=50)           | 100        | 1e+01      | 0.001672        | 8.170532e-16
Chebyshev (k=5)                | 100        | 1e+01      | 0.000618        | 4.494523e-03
Chebyshev (k=50)               | 100        | 1e+01      | 0.001946        | 7.851619e-10
Eigendecomposition             | 100        | 1e+04      | 0.001550        | 3.404796e-14
Lanczos (k=5)                  | 100        | 1e+04      | 0.000217        | 1.694096e-01
Lanczos (k=50)                 | 100        | 1e+04      | 0.001252        | 6.454948e-06
Newton-Schulz (k=5)            | 100        | 1e+04      | 0.001149        | 1.439958e-01
Newton-Schulz (k=50)           | 100        | 1e+04      | 0.009406        | nan
Chebyshev (k=5)                | 100        | 1e+04      | 0.000594        | 1.608013e-01
Chebyshev (k=50)               | 100        | 1e+04      | 0.001979        | 1.322931e-01
Eigendecomposition             | 100        | 1e+08      | 0.001639        | 1.366208e-10
Lanczos (k=5)                  | 100        | 1e+08      | 0.000226        | 8.369837e-02
Lanczos (k=50)                 | 100        | 1e+08      | 0.001283        | 2.498959e-02
Newton-Schulz (k=5)            | 100        | 1e+08      | 0.001148        | 2.866966e-02
Newton-Schulz (k=50)           | 100        | 1e+08      | 0.009538        | nan
Chebyshev (k=5)                | 100        | 1e+08      | 0.000631        | 6.278616e-02
Chebyshev (k=50)               | 100        | 1e+08      | 0.001970        | 2.543364e-02
Eigendecomposition             | 1000       | 1e+01      | 0.185755        | 3.372289e-15
Lanczos (k=5)                  | 1000       | 1e+01      | 0.001154        | 1.029816e-02
Lanczos (k=50)                 | 1000       | 1e+01      | 0.006844        | 2.073080e-15
Newton-Schulz (k=5)            | 1000       | 1e+01      | 0.398956        | 2.613246e-06
Newton-Schulz (k=50)           | 1000       | 1e+01      | 0.427680        | 1.552004e-15
Chebyshev (k=5)                | 1000       | 1e+01      | 0.002795        | 4.798983e-03
Chebyshev (k=50)               | 1000       | 1e+01      | 0.008187        | 4.544999e-10
Eigendecomposition             | 1000       | 1e+04      | 0.168302        | 1.204612e-14
Lanczos (k=5)                  | 1000       | 1e+04      | 0.001172        | 1.078717e-01
Lanczos (k=50)                 | 1000       | 1e+04      | 0.006977        | 3.079505e-02
Newton-Schulz (k=5)            | 1000       | 1e+04      | 0.278953        | 6.642356e-02
Newton-Schulz (k=50)           | 1000       | 1e+04      | 2.662046        | nan
Chebyshev (k=5)                | 1000       | 1e+04      | 0.002847        | 9.285141e-02
Chebyshev (k=50)               | 1000       | 1e+04      | 0.008859        | 5.156652e-02
Eigendecomposition             | 1000       | 1e+08      | 0.167850        | 1.330017e-10
Lanczos (k=5)                  | 1000       | 1e+08      | 0.001108        | 1.149061e-01
Lanczos (k=50)                 | 1000       | 1e+08      | 0.007183        | 7.165214e-02
Newton-Schulz (k=5)            | 1000       | 1e+08      | 0.265549        | 7.750715e-02
Newton-Schulz (k=50)           | 1000       | 1e+08      | 2.603875        | nan
Chebyshev (k=5)                | 1000       | 1e+08      | 0.002762        | 9.905950e-02
Chebyshev (k=50)               | 1000       | 1e+08      | 0.008436        | 7.242280e-02
--------------------------------------------------------------------------------
```

### Analysis

*   **Eigendecomposition**: This method remains the most accurate and is very fast for the matrix sizes tested.
*   **Lanczos**: The Lanczos method provides a good approximation for well-conditioned matrices, but its accuracy degrades significantly as the condition number increases.
*   **Newton-Schulz**: This solver outperforms Lanczos slightly on very small matrices. Its performance degrades with increasing condition number, and it can be unstable for ill-conditioned matrices.
*   **Chebyshev**: The Chebyshev solver is a competitive alternative to the other iterative methods. It is significantly faster than Eigendecomposition and Newton-Schulz for larger matrices. Its accuracy is comparable to Lanczos, and it is more stable than Newton-Schulz for ill-conditioned matrices.
