# Assignment 5 – LAPACK (HPC)

This repository contains solutions for Assignment 5 of the HPC Tools course.

## Part 1: Iterative Improvement
- Solves complex linear systems using LU decomposition (LAPACK)
- Applies iterative refinement
- Records number of iterations required for convergence for n = 2 to 500

## Part 2: Determinant Minimization
- Computes determinant using LU factorization
- Uses bisection method to find the value of t that minimizes det(A(t))
- Tested for odd matrix sizes n = 3, 5, 7, 9

## Implementation Notes
- All numerical computations were performed on the Seagull HPC system
- LAPACK (Fortran interface) was used for all linear algebra operations
- Plotting was performed locally due to lack of visualization tools on the HPC system

## Files
- `iterative_improvement.c`
- `determinant_bisection.c`
- `iterations_vs_n.dat`
- `determinant_results.dat`
