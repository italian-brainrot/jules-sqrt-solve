from sqrtsolvers.benchmarking.benchmark import run_benchmark, print_results, Solver
from sqrtsolvers.solvers import eigh, lanczos, newton_schulz

def main():
    """
    Main function to run the benchmarks.
    """
    solvers = [
        Solver("Eigendecomposition", lambda A, b: eigh.solve_eigh(A, b), warm_start=False),
        Solver("Lanczos (k=5)", lambda A, b: lanczos.solve_lanczos(A, b, k=5), warm_start=False),
        Solver("Lanczos (k=50)", lambda A, b: lanczos.solve_lanczos(A, b, k=50), warm_start=False),
        Solver("Newton-Schulz (k=5)", lambda A, b: newton_schulz.solve_ns(A, b, max_iter=5), warm_start=False),
        Solver("Newton-Schulz (k=50)", lambda A, b: newton_schulz.solve_ns(A, b, max_iter=50), warm_start=False),
    ]


    matrix_sizes = [10, 50, 100, 500, 1000]
    condition_numbers = [1e1, 1e2, 1e4, 1e8]

    results = run_benchmark(solvers, matrix_sizes, condition_numbers)
    print_results(results)

if __name__ == '__main__':
    main()
