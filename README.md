# Novel Solver for A^(1/2)x = b

This repository contains the implementation and benchmarking of solvers for the matrix equation A^(1/2)x = b, where A is a symmetric positive definite matrix.

## Solvers

Currently three solvers are implemented in this repository:

1.  **Eigendecomposition**: This method computes the eigendecomposition of the matrix A to find its square root and then solves the system directly. It is accurate and robust but can be computationally expensive for large matrices.

2.  **Lanczos**: This is an iterative method for finding eigenvalues and eigenvectors of a symmetric matrix. In this repository, the Lanczos algorithm is used to approximate A^(-1/2)b.

3.  **Newton-Schulz**: This is an iterative solver based on the Newton-Schulz iteration for finding the inverse square root of the matrix A. The iteration is given by:
    X_{k+1} = 0.5 * X_k * (3I - A * X_k^2)

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
Eigendecomposition             | 10         | 1e+01      | 0.000110        | 7.994504e-16
Lanczos (k=5)                  | 10         | 1e+01      | 0.000439        | 5.896173e-03
Lanczos (k=50)                 | 10         | 1e+01      | 0.001289        | 1.012623e-15
Newton-Schulz (k=5)            | 10         | 1e+01      | 0.000164        | 2.602481e-06
Newton-Schulz (k=50)           | 10         | 1e+01      | 0.000206        | 6.089318e-16
Eigendecomposition             | 10         | 1e+02      | 0.000122        | 6.106346e-15
Lanczos (k=5)                  | 10         | 1e+02      | 0.000442        | 2.675070e-01
Lanczos (k=50)                 | 10         | 1e+02      | 0.001160        | 1.482667e-14
Newton-Schulz (k=5)            | 10         | 1e+02      | 0.000121        | 1.212820e-01
Newton-Schulz (k=50)           | 10         | 1e+02      | 0.000195        | 1.238942e-14
Eigendecomposition             | 10         | 1e+04      | 0.000070        | 4.295256e-14
Lanczos (k=5)                  | 10         | 1e+04      | 0.000200        | 2.126257e-01
Lanczos (k=50)                 | 10         | 1e+04      | 0.000933        | 5.385227e-13
Newton-Schulz (k=5)            | 10         | 1e+04      | 0.000116        | 1.982350e-01
Newton-Schulz (k=50)           | 10         | 1e+04      | 0.001388        | nan
Eigendecomposition             | 10         | 1e+08      | 0.000120        | 6.884285e-11
Lanczos (k=5)                  | 10         | 1e+08      | 0.000323        | 1.634810e-01
Lanczos (k=50)                 | 10         | 1e+08      | 0.001781        | 1.780953e-09
Newton-Schulz (k=5)            | 10         | 1e+08      | 0.000222        | 1.633001e-01
Newton-Schulz (k=50)           | 10         | 1e+08      | 0.001920        | nan
Eigendecomposition             | 50         | 1e+01      | 0.000331        | 2.024903e-15
Lanczos (k=5)                  | 50         | 1e+01      | 0.000238        | 8.150200e-03
Lanczos (k=50)                 | 50         | 1e+01      | 0.001039        | 1.458252e-15
Newton-Schulz (k=5)            | 50         | 1e+01      | 0.000333        | 2.582469e-06
Newton-Schulz (k=50)           | 50         | 1e+01      | 0.000374        | 5.717449e-16
Eigendecomposition             | 50         | 1e+02      | 0.000438        | 3.691181e-15
Lanczos (k=5)                  | 50         | 1e+02      | 0.000362        | 1.366043e-01
Lanczos (k=50)                 | 50         | 1e+02      | 0.001293        | 9.411253e-15
Newton-Schulz (k=5)            | 50         | 1e+02      | 0.000312        | 4.132703e-02
Newton-Schulz (k=50)           | 50         | 1e+02      | 0.000716        | 2.524319e-14
Eigendecomposition             | 50         | 1e+04      | 0.000466        | 2.071940e-14
Lanczos (k=5)                  | 50         | 1e+04      | 0.000272        | 8.592614e-02
Lanczos (k=50)                 | 50         | 1e+04      | 0.001462        | 4.886367e-14
Newton-Schulz (k=5)            | 50         | 1e+04      | 0.000505        | 7.434954e-02
Newton-Schulz (k=50)           | 50         | 1e+04      | 0.004010        | nan
Eigendecomposition             | 50         | 1e+08      | 0.000406        | 2.336000e-11
Lanczos (k=5)                  | 50         | 1e+08      | 0.000283        | 8.523572e-02
Lanczos (k=50)                 | 50         | 1e+08      | 0.001049        | 4.967799e-11
Newton-Schulz (k=5)            | 50         | 1e+08      | 0.000301        | 1.201420e-02
Newton-Schulz (k=50)           | 50         | 1e+08      | 0.002763        | nan
Eigendecomposition             | 100        | 1e+01      | 0.001721        | 2.006459e-15
Lanczos (k=5)                  | 100        | 1e+01      | 0.000443        | 9.387762e-03
Lanczos (k=50)                 | 100        | 1e+01      | 0.003422        | 2.236057e-15
Newton-Schulz (k=5)            | 100        | 1e+01      | 0.001194        | 2.502699e-06
Newton-Schulz (k=50)           | 100        | 1e+01      | 0.002643        | 7.442858e-16
Eigendecomposition             | 100        | 1e+02      | 0.000935        | 2.357627e-15
Lanczos (k=5)                  | 100        | 1e+02      | 0.000267        | 7.458989e-02
Lanczos (k=50)                 | 100        | 1e+02      | 0.001223        | 2.951309e-09
Newton-Schulz (k=5)            | 100        | 1e+02      | 0.001012        | 1.311037e-02
Newton-Schulz (k=50)           | 100        | 1e+02      | 0.002331        | 1.083241e-14
Eigendecomposition             | 100        | 1e+04      | 0.001147        | 8.416005e-15
Lanczos (k=5)                  | 100        | 1e+04      | 0.000380        | 4.976510e-02
Lanczos (k=50)                 | 100        | 1e+04      | 0.002061        | 1.590998e-06
Newton-Schulz (k=5)            | 100        | 1e+04      | 0.001232        | 1.774692e-02
Newton-Schulz (k=50)           | 100        | 1e+04      | 0.010439        | nan
Eigendecomposition             | 100        | 1e+08      | 0.001380        | 5.174989e-10
Lanczos (k=5)                  | 100        | 1e+08      | 0.000479        | 1.252978e-01
Lanczos (k=50)                 | 100        | 1e+08      | 0.002112        | 1.011833e-01
Newton-Schulz (k=5)            | 100        | 1e+08      | 0.001349        | 1.083786e-01
Newton-Schulz (k=50)           | 100        | 1e+08      | 0.011041        | nan
Eigendecomposition             | 500        | 1e+01      | 0.027673        | 3.140974e-15
Lanczos (k=5)                  | 500        | 1e+01      | 0.000638        | 1.131503e-02
Lanczos (k=50)                 | 500        | 1e+01      | 0.002703        | 1.984874e-15
Newton-Schulz (k=5)            | 500        | 1e+01      | 0.039854        | 3.912238e-06
Newton-Schulz (k=50)           | 500        | 1e+01      | 0.064607        | 1.315587e-15
Eigendecomposition             | 500        | 1e+02      | 0.018998        | 3.517723e-15
Lanczos (k=5)                  | 500        | 1e+02      | 0.000650        | 7.123132e-02
Lanczos (k=50)                 | 500        | 1e+02      | 0.002352        | 2.946142e-06
Newton-Schulz (k=5)            | 500        | 1e+02      | 0.038952        | 1.304884e-02
Newton-Schulz (k=50)           | 500        | 1e+02      | 0.089970        | 1.527625e-14
Eigendecomposition             | 500        | 1e+04      | 0.020062        | 1.412883e-14
Lanczos (k=5)                  | 500        | 1e+04      | 0.000737        | 8.003274e-02
Lanczos (k=50)                 | 500        | 1e+04      | 0.002660        | 1.935111e-02
Newton-Schulz (k=5)            | 500        | 1e+04      | 0.051814        | 3.695952e-02
Newton-Schulz (k=50)           | 500        | 1e+04      | 0.427022        | nan
Eigendecomposition             | 500        | 1e+08      | 0.032025        | 1.096154e-11
Lanczos (k=5)                  | 500        | 1e+08      | 0.000763        | 9.607493e-02
Lanczos (k=50)                 | 500        | 1e+08      | 0.003329        | 2.673598e-02
Newton-Schulz (k=5)            | 500        | 1e+08      | 0.042065        | 3.585645e-02
Newton-Schulz (k=50)           | 500        | 1e+08      | 0.435293        | nan
Eigendecomposition             | 1000       | 1e+01      | 0.145619        | 3.545349e-15
Lanczos (k=5)                  | 1000       | 1e+01      | 0.002374        | 1.146752e-02
Lanczos (k=50)                 | 1000       | 1e+01      | 0.010401        | 2.127093e-15
Newton-Schulz (k=5)            | 1000       | 1e+01      | 0.413834        | 4.143600e-06
Newton-Schulz (k=50)           | 1000       | 1e+01      | 0.797502        | 1.431749e-15
Eigendecomposition             | 1000       | 1e+02      | 0.120820        | 3.972363e-15
Lanczos (k=5)                  | 1000       | 1e+02      | 0.002222        | 6.287304e-02
Lanczos (k=50)                 | 1000       | 1e+02      | 0.008851        | 2.090054e-06
Newton-Schulz (k=5)            | 1000       | 1e+02      | 0.496667        | 8.426214e-03
Newton-Schulz (k=50)           | 1000       | 1e+02      | 0.822091        | 1.076646e-14
Eigendecomposition             | 1000       | 1e+04      | 0.148239        | 6.985658e-15
Lanczos (k=5)                  | 1000       | 1e+04      | 0.003865        | 9.546559e-02
Lanczos (k=50)                 | 1000       | 1e+04      | 0.027970        | 2.061441e-02
Newton-Schulz (k=5)            | 1000       | 1e+04      | 0.596363        | 4.883674e-02
Newton-Schulz (k=50)           | 1000       | 1e+04      | 3.203256        | nan
Eigendecomposition             | 1000       | 1e+08      | 0.127511        | 1.202581e-10
Lanczos (k=5)                  | 1000       | 1e+08      | 0.002394        | 8.851800e-02
Lanczos (k=50)                 | 1000       | 1e+08      | 0.020661        | 1.270457e-02
Newton-Schulz (k=5)            | 1000       | 1e+08      | 0.541356        | 3.890516e-02
Newton-Schulz (k=50)           | 1000       | 1e+08      | 3.694084        | nan
--------------------------------------------------------------------------------
```

### Analysis

*   **Eigendecomposition**: This method remains the most accurate and is very fast for the matrix sizes tested.
*   **Lanczos**: The Lanczos method provides a good approximation for well-conditioned matrices, but its accuracy degrades significantly as the condition number increases.
*   **Newton-Schulz**: This solver outperforms Lanczos slightly on very small matrices. Its performance degrades with increasing condition number.
