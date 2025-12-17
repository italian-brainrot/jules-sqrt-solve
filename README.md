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
Eigendecomposition             | 10         | 1e+01      | 0.000150        | 1.161391e-15   
Lanczos (k=5)                  | 10         | 1e+01      | 0.000466        | 9.583825e-03   
Lanczos (k=50)                 | 10         | 1e+01      | 0.001660        | 1.559696e-15   
Newton-Schulz (k=5)            | 10         | 1e+01      | 0.000268        | 4.909819e-06   
Newton-Schulz (k=50)           | 10         | 1e+01      | 0.000272        | 5.118157e-16   
Eigendecomposition             | 10         | 1e+04      | 0.000113        | 2.324098e-14   
Lanczos (k=5)                  | 10         | 1e+04      | 0.000353        | 4.142411e-01   
Lanczos (k=50)                 | 10         | 1e+04      | 0.003323        | 4.112463e-13   
Newton-Schulz (k=5)            | 10         | 1e+04      | 0.000221        | 3.839981e-01   
Newton-Schulz (k=50)           | 10         | 1e+04      | 0.002784        | nan            
Eigendecomposition             | 10         | 1e+08      | 0.000118        | 4.159207e-09   
Lanczos (k=5)                  | 10         | 1e+08      | 0.000399        | 2.591260e-01   
Lanczos (k=50)                 | 10         | 1e+08      | 0.001951        | 2.164300e-09   
Newton-Schulz (k=5)            | 10         | 1e+08      | 0.000219        | 2.588782e-01   
Newton-Schulz (k=50)           | 10         | 1e+08      | 0.001642        | nan            
Eigendecomposition             | 100        | 1e+01      | 0.001516        | 2.005557e-15   
Lanczos (k=5)                  | 100        | 1e+01      | 0.000374        | 1.126066e-02   
Lanczos (k=50)                 | 100        | 1e+01      | 0.002500        | 1.826147e-15   
Newton-Schulz (k=5)            | 100        | 1e+01      | 0.002412        | 5.032733e-06   
Newton-Schulz (k=50)           | 100        | 1e+01      | 0.002613        | 9.446298e-16   
Eigendecomposition             | 100        | 1e+04      | 0.001340        | 5.993718e-15   
Lanczos (k=5)                  | 100        | 1e+04      | 0.000372        | 1.255828e-01   
Lanczos (k=50)                 | 100        | 1e+04      | 0.003615        | 5.569157e-06   
Newton-Schulz (k=5)            | 100        | 1e+04      | 0.002594        | 6.845138e-02   
Newton-Schulz (k=50)           | 100        | 1e+04      | 0.010947        | nan            
Eigendecomposition             | 100        | 1e+08      | 0.001297        | 3.542119e-10   
Lanczos (k=5)                  | 100        | 1e+08      | 0.000384        | 1.131044e-01   
Lanczos (k=50)                 | 100        | 1e+08      | 0.002252        | 9.372736e-02   
Newton-Schulz (k=5)            | 100        | 1e+08      | 0.002411        | 1.004419e-01   
Newton-Schulz (k=50)           | 100        | 1e+08      | 0.013175        | nan            
Eigendecomposition             | 1000       | 1e+01      | 0.243655        | 3.676470e-15   
Lanczos (k=5)                  | 1000       | 1e+01      | 0.003633        | 9.920668e-03   
Lanczos (k=50)                 | 1000       | 1e+01      | 0.019655        | 1.923925e-15   
Newton-Schulz (k=5)            | 1000       | 1e+01      | 0.618184        | 2.734460e-06   
Newton-Schulz (k=50)           | 1000       | 1e+01      | 0.779746        | 1.374514e-15   
Eigendecomposition             | 1000       | 1e+04      | 0.186703        | 7.635312e-15   
Lanczos (k=5)                  | 1000       | 1e+04      | 0.002818        | 1.125388e-01   
Lanczos (k=50)                 | 1000       | 1e+04      | 0.012387        | 2.378403e-02   
Newton-Schulz (k=5)            | 1000       | 1e+04      | 0.358648        | 5.879788e-02   
Newton-Schulz (k=50)           | 1000       | 1e+04      | 3.967155        | nan            
Eigendecomposition             | 1000       | 1e+08      | 0.124042        | 1.626992e-10   
Lanczos (k=5)                  | 1000       | 1e+08      | 0.002192        | 1.035749e-01   
Lanczos (k=50)                 | 1000       | 1e+08      | 0.007194        | 3.863894e-02   
Newton-Schulz (k=5)            | 1000       | 1e+08      | 0.518498        | 5.362316e-02   
Newton-Schulz (k=50)           | 1000       | 1e+08      | 3.261142        | nan            
--------------------------------------------------------------------------------
```

### Analysis

*   **Eigendecomposition**: This method remains the most accurate and is very fast for the matrix sizes tested.
*   **Lanczos**: The Lanczos method provides a good approximation for well-conditioned matrices, but its accuracy degrades significantly as the condition number increases.
*   **Newton-Schulz**: This solver outperforms Lanczos slightly on very small matrices. Its performance degrades with increasing condition number.
