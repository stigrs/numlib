# Benchmarks

The benchmarks show that Numlib can compete with Armadillo when a fast BLAS library is used. BLAS was not used for the
DAXPY operation. In this case Armadillo performs significantly better than Numlib. However, Numlib is able to compete 
with valarray on Windows. 

_Note:_ All benchmarks were done with OPENBLAS_NUM_THREADS = 1.

#### DAXPY operation 

* Linux (g++ 7.3.0):

| Size                 | 10   | 100 | 1000 | 10000 | 100000 |
|----------------------|------|-----|------|-------|--------|
| T_numlib/T_armadillo | 5.2  | 3.0 | 12.2 | 9.2   | 8.0    |
| T_numlib/T_valarray  | 14.2 | 4.5 | 14.0 | 10.2  | 14.4   |

* Windows (MSVC 2017):

| Size                 | 10   | 100 | 1000 | 10000 | 100000 |
|----------------------|------|-----|------|-------|--------|
| T_numlib/T_armadillo | 4.0  | 3.5 | 13.0 | 9.3   | 4.1    |
| T_numlib/T_valarray  | 2.7  | 1.0 | 1.0  | 1.7   | 1.7    |

#### DGEMM operation 

* Linux (g++ 7.3.0):

| Size                 | 10 x 5 | 100 x 50 | 1000 x 500 | 10000 x 5000 |
|----------------------|--------|----------|------------|--------------|
| T_numlib/T_armadillo | 0.2    | 0.85     | 0.53       | 1.05         |

* Windows (MSVC 2017):

| Size                 | 10 x 5 | 100 x 50 | 1000 x 500 | 
|----------------------|--------|----------|------------|
| T_numlib/T_armadillo | 0.03   | 0.88     | 0.98       | 

#### DGEMV operation

* Linux (g++ 7.3.0):

| Size                 | 10 x 5 | 100 x 50 | 1000 x 500 | 10000 x 5000 |
|----------------------|--------|----------|------------|--------------|
| T_numlib/T_armadillo | 0.21   | 0.79     | 0.92       | 0.86         |

* Windows (MSVC 2017):

| Size                 | 10 x 5 | 100 x 50 | 1000 x 500 | 10000 x 5000 |
|----------------------|--------|----------|------------|--------------|
| T_numlib/T_armadillo | 0.014  | 0.27     | 0.70       | 0.93         |

