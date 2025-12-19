from sqrtsolvers.benchmarking.benchmark import run_benchmark, print_results, plot_results, Solver
from sqrtsolvers.solvers import eigh, lanczos, chebyshev, newton_schulz

def main():
    """
    Main function to run the benchmarks.
    """
    solvers = [
        Solver("Eigendecomposition", lambda A, b: eigh.solve_eigh(A, b), warm_start=False),
        Solver("Lanczos (k=10)", lambda A, b: lanczos.solve_lanczos(A, b, k=10), warm_start=False),
        Solver("Lanczos (k=100)", lambda A, b: lanczos.solve_lanczos(A, b, k=100), warm_start=False),
        Solver("Chebyshev", lambda A, b: chebyshev.solve_chebyshev(A, b), warm_start=False),
        Solver("Newton-Schulz", lambda A, b: newton_schulz.solve_ns(A, b), warm_start=False),
    ]


    matrix_sizes = [10, 100, 1000]
    condition_numbers = [1e1, 1e4, 1e8]

    results = run_benchmark(solvers, matrix_sizes, condition_numbers)
    print_results(results)
    plot_results(results)

if __name__ == '__main__':
    main()
