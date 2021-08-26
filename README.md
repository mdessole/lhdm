# Lawson-Hanson algorithm with Deviation Maximization pivoting

Lawson-Hanson active set algorithm with Deviation Maximization pivoting implemented in C. A matlab wrapper is provided. 
NNLS problem: min ||A*X-B||, s.t. X>=0. 
Implementation used in

> Monica Dessole, Marco Dell'Orto, Fabio Marcuzzi "**[The Lawson-Hanson Algorithm with Deviation Maximization: Finite Convergence and Sparse Recovery](https://arxiv.org/abs/2108.05345)**", Preprint, 2021.

Compile mex function from matlab command line

```console
mex -v -R2017b lhdm.c -lmwblas -lmwlapack
```

## Dependencies

- Matlab 2017 or higher
- LAPACK
- BLAS

