# SSM 1.0.4
* created sensitivity_permutations.R which will hold functions for permutation
  tests of FANOVA variances.
* constructH computes the H matrix used in the theta = Hy formulation.
* constructDelta computes the squared basis norms assuming Legendre polynomials
  are used.
* constructQ uses construct H and constructDelta to compute Q in the y^TQy
  formulation of FANOVA variance.

# SSM 1.0.3
* added sanity check test that ssm do actually interpolate at design points using   1 and 2 d models.

# SSM 1.0.2

* now exports new.distance (because I need to use it during simulations)

# SSM 1.0.1

* Defines S4 class SSM object
* Provides plot and predict methods for SSM object
* fit.ssm returns an SSM object representing a fitted model
* sensitivity.plot plots the Sobol or Total indices for an SSM object
