# PGDoubling-quad

This is a fork of [pgdoubling](https://bitbucket.org/fph/pgdoubling/) in which we introduce the ability to work with *quad-bases*, that is, using a triple (A,B,C) \in ℂ^{(n-k) × k} × ℂ^{(n-k) × t} × ℂ^{r × k} to represent the quasidefinite matrix
```
X=[ -C^T*C  A^T   ]
  [   A     B*B^T ]
```
appearing in the context of principal pivot transforms and permuted Riccati bases of Lagrangian subspaces.

The theory behind these algorithms is described in [Poloni, Strabić, *Principal pivot transforms of quasidefinite matrices and semidefinite Lagrangian subspaces* (provisional title, still to be submitted)].

# Usage

A *quad basis* is a structure containing the boolean vector `quad.v`, which is to be interpreted in the same way as `sym.v` in a symbasis from `pgdoubling`, and the matrix `quad.X=[C *; A B]`. The first block of this partition corresponds to indices for which `quad.v==true`, while the second block corresponds to indices for which `quad.v==false`. The block labelled `*` is arbitrary (and the algorithms are free to modify it for performance reasons).

The function `formQuadBasis` can be used to obtain a quadBasis from a triple $(A,B,C)$ coming from a control theory problem. Then, the function `optimizeQuadBasis` can be used to obtain an alternative one with a smaller $|X|$. Conversion functions `symBasisFromQuadBasis`, `hamiltonianPencilFromSymBasis` and `symplecticSubspaceFromSymBasis` (slightly misnamed, since it produces a *Lagrangian* subspace) are provided to switch between the different formats used in the algorithm.

The functions `reduceChuLM07CAREX` and `getFigures` are used to generate the numerical experiments in the papers.

# Authors

Federico Poloni <fpoloni@di.unipi.it> and Nataša Strabić <natasa.strabic@manchester.ac.uk>.
