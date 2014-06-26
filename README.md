GPML
====

Extensions for GPML in Matlab
- GPML extensions for sparse covariance function and derivative observations
- To see the details, run test/test_all_cov.m

Features
----
1. Sparse covariance function
  - **covSparseiso.m**
  - reference : [A. Melkumyan and F. Ramos, "A Sparse Covariance Function for Exact Gaussian Process Inference in Large Datasets," Proceedings of International Joint Conferences on Artificial Intelligence, pp. 1936-1942, 2009.](http://ijcai.org/papers09/Papers/IJCAI09-320.pdf)

2. Derivative observations (*: SEiso, Materniso3, and Sparseiso)
  1. **cov*DiffUnstable.m**
    - They consist of multiplications of partial derivatives.
    - They produce NaN when a partial derivative includes divide-by-zero. (Thus it is called unstable.)
    - They are used for validation of cov*Diff with ignoring divide-by-zero.
    - Do not use them for your application.

  2. **cov*Diff.m**
    - They reduce complex multiplications and additions to simpler equations.
    - They care divide-by-zero with limits. (Thus it is stable.)
    - Use them for your application.
  

