GPML
====

Extensions for GPML in Matlab
- GPML extensions for sparse covariance function and derivative observations
- To see the details, run test/test_all_cov.m

Features
1. Sparse covariance function (covSparseiso)
  (ref: A. Melkumyan and F. Ramos, "A Sparse Covariance Function for Exact Gaussian Process Inference in Large Datasets," IJCAI 2009.)

2. Derivative observations
2.1 cov*DiffUnstable
: They consist of multiplications of partial derivatives.
: They produce NaN when a partial derivative includes divide-by-zero. (Thus it is called unstable.)
: They are used for validation of cov*Diff with ignoring divide-by-zero.
: Do not use them for your application.

2.2 cov*Diff
: They reduce complex multiplications and additions to simpler equations.
: They care divide-by-zero with limits. (Thus it is stable.)
: Use them for your application.
  

