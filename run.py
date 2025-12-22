from sqrtsolvers.benchmarking.benchmark import run_benchmark, print_results, plot_results, Solver
from sqrtsolvers.solvers import eigh, lanczos, chebyshev, newton_schulz, cg_lanczos, preconditioned_cg, pcg_chebyshev_tuned

def main():
    """
    Main function to run the benchmarks.
    """
    solvers = [
        Solver("Eigendecomposition", lambda A, b: eigh.solve_eigh(A, b), warm_start=False),
        Solver("Lanczos (k=10)", lambda A, b: lanczos.solve_lanczos(A, b, k=10), warm_start=False),
        Solver("Lanczos (k=100)", lambda A, b: lanczos.solve_lanczos(A, b, k=100), warm_start=False),
        Solver("Chebyshev", lambda A, b: chebyshev.solve_chebyshev(A, b), warm_start=False),
        Solver("Newton-Schulz (k=1)", lambda A, b: newton_schulz.solve_ns(A, b, max_iter=1), warm_start=False),
        Solver("Newton-Schulz (k=3)", lambda A, b: newton_schulz.solve_ns(A, b, max_iter=3), warm_start=False),
        Solver("CG Lanczos (k=10)", lambda A, b: cg_lanczos.solve_cg_lanczos(A, b, max_iter=10), warm_start=False),
        Solver("CG Lanczos (k=30)", lambda A, b: cg_lanczos.solve_cg_lanczos(A, b, max_iter=30), warm_start=False),
        Solver("PCG-Chebyshev", lambda A, b, x0=None: preconditioned_cg.solve_pcg_chebyshev(A, b, x0=x0, max_iter=10, lanczos_k=10, chebyshev_k=10), warm_start=True),
        Solver("PCG-Chebyshev-Tuned", lambda A, b, x0=None: pcg_chebyshev_tuned.solve_pcg_chebyshev_tuned(A, b, x0=x0, max_iter=10), warm_start=True),
    ]


    matrix_sizes = [10, 100, 1000]
    condition_numbers = [1e1, 1e4, 1e8]

    results = run_benchmark(solvers, matrix_sizes, condition_numbers)
    print_results(results)
    plot_results(results)

if __name__ == '__main__':
    main()
