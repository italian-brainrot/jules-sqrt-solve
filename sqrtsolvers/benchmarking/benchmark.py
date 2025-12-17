import time
from collections.abc import Callable
from typing import NamedTuple

import numpy as np


def generate_spd_matrix(n, cond):
    """
    Generates a symmetric positive definite matrix of size n x n with a given condition number.

    Args:
        n (int): The size of the matrix.
        cond (float): The condition number of the matrix.

    Returns:
        np.ndarray: The generated matrix A.
        np.ndarray: The matrix A^(1/2).
    """
    # Generate a random orthogonal matrix
    Q, _ = np.linalg.qr(np.random.randn(n, n))

    # Generate eigenvalues with a specified condition number
    L = np.linspace(1, cond, n)

    A = Q @ np.diag(L) @ Q.T
    A_sqrt = Q @ np.diag(np.sqrt(L)) @ Q.T
    return A, A_sqrt

class Solver(NamedTuple):
    name: str
    func: Callable
    warm_start: bool

def run_benchmark(solvers: list[Solver], matrix_sizes: list[int], condition_numbers: list[float]):
    """
    Runs the benchmark for the given solvers, matrix sizes, and condition numbers.

    Args:
        solvers (dict): A dictionary of solver names to solver functions.
        matrix_sizes (list): A list of matrix sizes to test.
        condition_numbers (list): A list of condition numbers to test.

    Returns:
        list: A list of dictionaries containing the benchmark results.
    """
    results = []
    for n in matrix_sizes:
        for cond in condition_numbers:
            print(f"Running benchmark for size={n}, cond={cond:.0e}")
            A, A_sqrt = generate_spd_matrix(n, cond)
            x_true = np.random.randn(n)
            b = A_sqrt @ x_true

            # generate data for warm starting
            P = np.random.randn(n, max(n // 10, 1))
            A2 = A + (P @ P.T) * np.linalg.norm(A, ord='fro') * 0.2
            L_A2, Q_A2 = np.linalg.eigh(A2)
            L_A2 = np.abs(L_A2) + 1e-12
            x2_true = x_true + np.random.randn(n) * 0.2
            A2_sqrt = Q_A2 @ np.diag(np.sqrt(L_A2)) @ Q_A2.T
            b2 = A2_sqrt @ x2_true

            for solver in solvers:

                x0 = None
                if solver.warm_start:
                    x0 = solver.func(A2, b2, x0=None)

                start_time = time.perf_counter()

                if solver.warm_start:
                    x = solver.func(A, b, x0=x0)
                else:
                    x = solver.func(A, b)

                end_time = time.perf_counter()

                error = np.linalg.norm(x - x_true) / np.linalg.norm(x_true)

                results.append({
                    'solver': solver.name,
                    'size': n,
                    'condition_number': cond,
                    'time': end_time - start_time,
                    'error': error,
                })
    return results

def print_results(results):
    """
    Prints the benchmark results in a formatted table.

    Args:
        results (list): The list of benchmark results.
    """
    print("-" * 80)
    print(f"{'Solver':<30} | {'Size':<10} | {'Cond':<10} | {'Time (s)':<15} | {'Error':<15}")
    print("-" * 80)
    for res in results:
        print(f"{res['solver']:<30} | {res['size']:<10} | {res['condition_number']:<10.0e} | {res['time']:<15.6f} | {res['error']:<15.6e}")
    print("-" * 80)
