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
Eigendecomposition             | 10         | 1e+01      | 0.000107        | 4.595680e-16   
Lanczos (k=5)                  | 10         | 1e+01      | 0.000444        | 9.669811e-03   
Lanczos (k=50)                 | 10         | 1e+01      | 0.001376        | 1.223367e-15   
Newton-Schulz (k=5)            | 10         | 1e+01      | 0.000240        | 3.204304e-05   
Newton-Schulz (k=50)           | 10         | 1e+01      | 0.000209        | 4.876220e-16   
Chebyshev (k=5)                | 10         | 1e+01      | 0.001062        | 1.497139e-02   
Chebyshev (k=50)               | 10         | 1e+01      | 0.002558        | 5.085634e-15   
Eigendecomposition             | 10         | 1e+04      | 0.000083        | 3.713821e-14   
Lanczos (k=5)                  | 10         | 1e+04      | 0.000259        | 2.685840e-01   
Lanczos (k=50)                 | 10         | 1e+04      | 0.001276        | 4.585072e-13   
Newton-Schulz (k=5)            | 10         | 1e+04      | 0.000166        | 2.505447e-01   
Newton-Schulz (k=50)           | 10         | 1e+04      | 0.001566        | nan            
Chebyshev (k=5)                | 10         | 1e+04      | 0.000769        | 1.343778e+00   
Chebyshev (k=50)               | 10         | 1e+04      | 0.002231        | 4.680878e-06   
Eigendecomposition             | 10         | 1e+08      | 0.000081        | 2.308165e-10   
Lanczos (k=5)                  | 10         | 1e+08      | 0.000250        | 1.133393e-01   
Lanczos (k=50)                 | 10         | 1e+08      | 0.001288        | 3.624546e-09   
Newton-Schulz (k=5)            | 10         | 1e+08      | 0.000176        | 1.132503e-01   
Newton-Schulz (k=50)           | 10         | 1e+08      | 0.001195        | nan            
Chebyshev (k=5)                | 10         | 1e+08      | 0.000772        | 2.081412e+02   
Chebyshev (k=50)               | 10         | 1e+08      | 0.002227        | 3.838228e-04   
Eigendecomposition             | 100        | 1e+01      | 0.001264        | 2.112677e-15   
Lanczos (k=5)                  | 100        | 1e+01      | 0.000324        | 9.073770e-03   
Lanczos (k=50)                 | 100        | 1e+01      | 0.001605        | 1.239504e-15   
Newton-Schulz (k=5)            | 100        | 1e+01      | 0.001222        | 1.904760e-06   
Newton-Schulz (k=50)           | 100        | 1e+01      | 0.001483        | 8.567713e-16   
Chebyshev (k=5)                | 100        | 1e+01      | 0.000910        | 3.771421e-03   
Chebyshev (k=50)               | 100        | 1e+01      | 0.002570        | 3.072078e-10   
Eigendecomposition             | 100        | 1e+04      | 0.001061        | 1.619330e-14   
Lanczos (k=5)                  | 100        | 1e+04      | 0.000307        | 1.113959e-01   
Lanczos (k=50)                 | 100        | 1e+04      | 0.001686        | 3.173401e-06   
Newton-Schulz (k=5)            | 100        | 1e+04      | 0.001080        | 7.577578e-02   
Newton-Schulz (k=50)           | 100        | 1e+04      | 0.011200        | nan            
Chebyshev (k=5)                | 100        | 1e+04      | 0.001063        | 9.422034e-02   
Chebyshev (k=50)               | 100        | 1e+04      | 0.002487        | 6.959557e-02   
Eigendecomposition             | 100        | 1e+08      | 0.001052        | 1.289877e-11   
Lanczos (k=5)                  | 100        | 1e+08      | 0.000318        | 1.047424e-01   
Lanczos (k=50)                 | 100        | 1e+08      | 0.001610        | 2.347537e-02   
Newton-Schulz (k=5)            | 100        | 1e+08      | 0.001142        | 3.863038e-02   
Newton-Schulz (k=50)           | 100        | 1e+08      | 0.010924        | nan            
Chebyshev (k=5)                | 100        | 1e+08      | 0.001037        | 8.661828e-02   
Chebyshev (k=50)               | 100        | 1e+08      | 0.002537        | 2.379472e-02   
Eigendecomposition             | 1000       | 1e+01      | 0.158845        | 3.543412e-15   
Lanczos (k=5)                  | 1000       | 1e+01      | 0.004145        | 1.143928e-02   
Lanczos (k=50)                 | 1000       | 1e+01      | 0.028217        | 1.895816e-15   
Newton-Schulz (k=5)            | 1000       | 1e+01      | 0.429936        | 3.458694e-06   
Newton-Schulz (k=50)           | 1000       | 1e+01      | 0.597620        | 1.346602e-15   
Chebyshev (k=5)                | 1000       | 1e+01      | 0.011259        | 5.130858e-03   
Chebyshev (k=50)               | 1000       | 1e+01      | 0.036660        | 2.266691e-11   
Eigendecomposition             | 1000       | 1e+04      | 0.209702        | 1.182377e-14   
Lanczos (k=5)                  | 1000       | 1e+04      | 0.002502        | 9.994802e-02   
Lanczos (k=50)                 | 1000       | 1e+04      | 0.012518        | 2.422650e-02   
Newton-Schulz (k=5)            | 1000       | 1e+04      | 0.443195        | 5.054681e-02   
Newton-Schulz (k=50)           | 1000       | 1e+04      | 4.097119        | nan            
Chebyshev (k=5)                | 1000       | 1e+04      | 0.012188        | 8.265538e-02   
Chebyshev (k=50)               | 1000       | 1e+04      | 0.046072        | 3.660731e-02   
Eigendecomposition             | 1000       | 1e+08      | 0.151663        | 1.709811e-10   
Lanczos (k=5)                  | 1000       | 1e+08      | 0.002596        | 8.168562e-02   
Lanczos (k=50)                 | 1000       | 1e+08      | 0.015394        | 3.435528e-02   
Newton-Schulz (k=5)            | 1000       | 1e+08      | 0.318573        | 3.859637e-02   
Newton-Schulz (k=50)           | 1000       | 1e+08      | 3.367172        | nan            
Chebyshev (k=5)                | 1000       | 1e+08      | 0.004304        | 6.345043e-02   
Chebyshev (k=50)               | 1000       | 1e+08      | 0.007866        | 3.472298e-02   
--------------------------------------------------------------------------------
```

### Analysis

*   **Eigendecomposition**: This method remains the most accurate and is very fast for the matrix sizes tested.
*   **Lanczos**: The Lanczos method provides a good approximation for well-conditioned matrices, but its accuracy degrades significantly as the condition number increases.
*   **Newton-Schulz**: This solver outperforms Lanczos slightly on very small matrices. Its performance degrades with increasing condition number, and it can be unstable for ill-conditioned matrices.
*   **Chebyshev**: The Chebyshev solver is a competitive alternative to the other iterative methods. It is significantly faster than Eigendecomposition and Newton-Schulz for larger matrices. Its accuracy is comparable to Lanczos, and it is more stable than Newton-Schulz for ill-conditioned matrices.
