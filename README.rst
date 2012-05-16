Purpose
=======
This package uses the algorithms described in [MehP12,MehP12_ppt] to compute Lagrangian invariant subspaces of Hamiltonian and "extended Hamiltonian" pencils. You can use it to solve algebraic Riccati equations, Lur'e equations/singular control problems, and finding the optimal gamma in H_infinity control.

Terminology
===========
A *CanBasis* is a basis matrix for an n-dimensional subspace in the form ``P*[I(n);X]``, where P is a permutation matrix and X is arbitrary. Usually the routines return CanBases in which each element of X is smaller (in modulus) than a certain configurable threshold.

A *SymBasis* is the structured version of a CanBasis, but for a Lagrangian subspace. P is a "symplectic swap matrix" (a symplectic analogous of a permutation), and X is Hermitian.

You can represent subspaces, matrices and matrix pencils through canBases and symBases. See the basic theory in [Meh12]. (In the following, all will be clearer if you have read this paper.)

Documentation
=============
All functions contain help; type ``help functionName`` in Matlab to access it. You can look in the "test" functions for some more usage examples. In the following there is a quick introduction.

Basic usage: matrix equations
=============================
::

[X,Y,U,V]=solveCARE(A,G,Q,...options...)

Returns the symmetric (semi)stabilizing solution of a CARE ``A'X+XA+Q=XGX``. Also returns the (semi)stabilizing solution of the "dual CARE" ``YA'+AY+YQY=G``. More important, U and V are basis matrices for the (semi)stabilizing and (semi)antistabilizing Lagrangian subspaces of the Hamiltonian, i.e., subspaces ``U=[U1;U2]`` and ``V=[V1;V2]`` such that ``X=U2/U1``, ``Y=V1/V2``. **You should always try to restate your equations in order to use U and V rather than X and Y**, as X and Y may be ill-conditioned while U and V are tame.

Recommended options for robustness: ``solveCARE(A,G,Q,'type','sign','safer',true)``. If your problem has large size, you may consider increasing the threshold with ``(...'threshold',sqrt(n),'diagonalThreshold',sqrt(n))``. The impact of the thresholds hasn't been investigated thoroughly yet.

::

[X,Y,U,V]=solveECARE(A,B,Q,R,S,...options...)

solve a (possibly singular) generic control problem in the form::

  [ 0  I 0]    [0  A  B]
  [ -I 0 0]s - [A' Q  S],
  [ 0  0 0]    [B' S' R]

returning the (semi)stabilizing X and Y. Following the more standard notation of this problem, U and V now are bases for the subspaces [X;I] (stabilizing) and [I;Y] (antistabilizing). Again, you should work with them rather than with X and Y.

::

 gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol,...options...);

run the gamma iteration [BenBMX07,MehP12_ppt] on a H_inf control problem, with tolerance tol.

Some more matrix equations could be added --- I expect these techniques to work also for X=Q+A*inv(X)*A', X=Q-A*inv(X)*A', and 0=P+QY+RY^2 (with some modifications).

Basic usage: building blocks
============================

``subspace2CanBasis`` takes a subspace and returns X and a permutation P defining its CanBasis.
``symplecticSubspace2SymBasis`` takes a Lagrangian subspace and returns X and a vector v defining its SymBasis. The "real" P is obtained by ``rowSwap(eye(2*n),v,'N')``, if needed.

You may also be interested in ``symplecticPencil2SymBasis``, ``hamiltonianPencil2SymBasis``, ``evenPencil2SymBasis`` (actually only works for even pencils in the form above, a generic even pencil does not have a symBasis).

Some random stuff that an end user might find interesting: ``checkCAREInvariantSubspaceResidual`` (computes and returns several error measures), ``Hamiltonian2RiccatiCoefficients``, ``hamiltonian`` (should be named ``RiccatiCoefficients2Hamiltonian``), ``pirq`` ((symplectic swap P)*R*Q factorization), ``ChuLM07Carex`` (returns the experiments used in [ChuLM07]), ``randomLagrangianSubspace``.

Feedback
========
The code is quite new at the moment; I am happy to hear about end users' experiences (I tried it and it works great; I tried it and it doesn't work) and to receive feedback and bug reports.
Also, if you use this code in a paper, please cite us when it is relevant. In the academia we basically live off citation counts.

License
=======
Check ``COPYING.txt`` to see if you qualify to use the library for free.

FAQ
===
Q: why are the function names so long?
A: I find it more expressive than the default 'shorten everything' (and slightly inconsistent) Matlab syntax. If you use Matlab's editor, you can configure it for tab completion and they will be fast to type as if they were short. If you use some other editor, well, I hope it supports tab completion. After all, you're using it because it's better than Matlab's builtin one, right?

Q: why is it so slow?
A: mainly, because of the \Pi RQ factorization (which should be rewritten in a MEX file, but MEX files are a mess) and because of option parsing (which could be removed, at the expenses of configurability and parameter-tweaking-ability). In general, some more performance could be squeezed out, the functions are not aggressively optimized.

Q: why do you roll your own option parsing functions? OOP is slow as hell in Matlab!
A: I haven't found yet an option parser which satisfy some basic requirements, namely: (1) can pass options unchanged to subroutines (2) complains at the end if and only if there are unused parameters (3) has an acceptable syntax which does not require me to specify option names in 18 different places. If you know of one, let me know.

Q: does it work when there is a mass matrix E?
A: no, and it would require significant changes to the algorithm to support it. We'll be working on that hopefully.

References
==========

[MehP12]
  Volker Mehrmann, Federico Poloni "Doubling Algorithms With Permuted Lagrangian Graph Bases". To appear in SIMAX.

[Meh12_ppt]
  Volker Mehrmann, Federico Poloni "Robust control with doubling and permuted Lagrangian bases" (provisional title). In preparation.

[ChuLM07]
  Chu, Delin; Liu, Xinmin; Mehrmann, Volker A numerical method for computing the Hamiltonian Schur form. Numer. Math. 105 (2007), no. 3, 375–412.

[BenBMX07]
  Benner, Peter; Byers, Ralph; Mehrmann, Volker; Xu, Hongguo A robust numerical method for the γ-iteration in H∞ control. Linear Algebra Appl. 425 (2007), no. 2-3, 548–570.

