Purpose
=======

This package uses the algorithms described in [MehP12, MehP13] to
compute Lagrangian invariant subspaces of Hamiltonian and "extended
Hamiltonian" pencils. The resulting algorithm is more accurate than
Matlab's `care` on several tricky test problems [carex].

It is the right choice for you if: 

* you want to compute a
**linear-quadratic regulator**, solve a **Riccati equation**, or a
generalized Riccati equation (those that correspond to a matrix pencil
in the form:
:::text
    [ 0  I 0]    [0  A  B]
    [ -I 0 0]s - [A' Q  S]
    [ 0  0 0]    [B' S' R]

); 

* your problem might also be a **singular control** problem (Lur'e
equations, singular `R`); there is also a function to compute the
optimal parameter gamma in **H\_infinity control**; 

* you are not
satisfied with Matlab's care, or wish to test a new algorithm; 

* your problem is **small-scale** (up to a few hundred degrees of freedom); 

* your system is **not a descriptor system** (i.e., no `E` coefficient,
just `\dot{x}=Ax`).

The current implementation is in **pure Matlab** (no mex-files), so it
is **slower** that `care`. Rewriting the more intensive parts in Fortran
or C will very likely result in a huge speedup.

Terminology
===========

A *CanBasis* is a basis matrix for an n-dimensional subspace in the form
`P*[I(n);X]`, where P is a permutation matrix and X is arbitrary.
Usually the routines return CanBases in which each element of X is
smaller (in modulus) than a certain configurable threshold.

A *SymBasis* is the structured version of a CanBasis for a Lagrangian
subspace. P is a "symplectic swap matrix" (a symplectic analogous of a
permutation), and X is Hermitian.

You can represent subspaces, matrices and matrix pencils through
canBases and symBases. See the basic theory in [MehP12]. (all will be
clearer if you have read this paper.)

Documentation
=============

All functions contain help; type `help functionName` in Matlab to access
it. You can look in the "test" functions for some more usage examples.
In the following there is a quick introduction.

Basic usage: matrix equations
=============================

```matlab
    [X,Y,U,V]=solveCARE(A,G,Q,...options...)
```
Returns the symmetric (semi)stabilizing solution of a CARE
`A'X+XA+Q=XGX`. Also returns the (semi)stabilizing solution of the "dual
CARE" `YA'+AY+YQY=G`. More important, U and V are basis matrices for the
(semi)stabilizing and (semi)antistabilizing Lagrangian subspaces of the
Hamiltonian, i.e., subspaces `U=[U1;U2]` and `V=[V1;V2]` such that
`X=U2/U1`, `Y=V1/V2`. **You should always try to restate your successive
computations in order to use U and V rather than X and Y**, as X and Y
may be ill-conditioned while U and V are tame.

Recommended options for robustness:
`solveCARE(A,G,Q,'type','sign','safer',true)`. If your problem has large
size, you may consider increasing the threshold with
`(...'threshold',sqrt(n),'diagonalThreshold',sqrt(n))`. The impact of
the thresholds hasn't been investigated thoroughly yet.

```matlab
[X,Y,U,V]=solveECARE(A,B,Q,R,S,...options...)
```
solve a (possibly singular) generic control problem in the form:

:::text
    [ 0  I 0]    [0  A  B]
    [ -I 0 0]s - [A' Q  S],
    [ 0  0 0]    [B' S' R]

returning the (semi)stabilizing X and Y. Following the more standard
notation of this problem, U and V now are bases for the subspaces [X;I]
(stabilizing) and [I;Y] (antistabilizing). Again, you should work with
them rather than with X and Y.
```matlab
    gamma=gammaIteration(A1,B1,B2,C1,C2,D11,D12,D21,D22,tol,...options...);
```
run the gamma iteration [BenBMX07,MehP13] on a H\_inf control problem,
with tolerance tol.

Some more matrix equations could be added --- These techniques work also
for DARE, X=Q+A\*inv(X)*A', X=Q-A*inv(X)\*A', and 0=P+QY+RY\^2 (with
some modifications).

Basic usage: building blocks
============================

`subspace2CanBasis` takes a subspace and returns X and a permutation P
defining its CanBasis. `symplecticSubspace2SymBasis` takes a Lagrangian
subspace and returns X and a vector v defining its SymBasis. The "real"
P is obtained by `rowSwap(eye(2*n),v,'N')`, if needed.

You may also be interested in `symplecticPencil2SymBasis`,
`hamiltonianPencil2SymBasis`, `evenPencil2SymBasis` (actually only works
for even pencils in the form above, a generic even pencil does not have
a symBasis).

Some random stuff that an end user might find interesting:
`checkCAREInvariantSubspaceResidual` (computes and returns several error
measures), `Hamiltonian2RiccatiCoefficients`, `hamiltonian` (should be
named `RiccatiCoefficients2Hamiltonian`), `pirq` ((symplectic swap
P)*R*Q factorization), `ChuLM07Carex` (returns the results of the
experiments shown in [ChuLM07]), `randomLagrangianSubspace`.

Feedback
========

The code is quite new at the moment; I am happy to hear about end users'
experiences ("I tried it and it works great; I tried it and it doesn't
work") and to receive feedback and bug reports. Also, if you use this
code in a paper, please cite us when it is relevant. In the academia we
basically live off citation counts.

License
=======

Check `COPYING.txt` to see if you qualify to use the library in other
projects for free.

FAQ
===

Q: why are the function names so long?

A: I find it more expressive than the default 'shorten everything' (and
slightly inconsistent) Matlab syntax. If you use Matlab's editor, you
can configure it for tab completion and they will be fast to type as if
they were short. If you use some other editor, well, I hope it supports
tab completion. After all, you're using it because it's better than
Matlab's builtin one, right?

Q: why is it so slow?

A: mainly, because of the Pi RQ factorization (which should be rewritten
in a MEX file, but MEX files are a mess) and because of option parsing
(which could be removed, at the expenses of configurability and
parameter-tweaking-ability). In general, some more performance could be
squeezed out, the functions are not aggressively optimized.

Q: why do you roll your own option parsing functions? OOP is slow as
hell in Matlab!

A: I haven't found yet an option parser that satisfies some basic
requirements, namely: (1) can pass options unchanged to subroutines (2)
complains at the end if there are unused parameters (3) has an
acceptable syntax which does not require me to specify option names in
18 different places. If you know of one, let me know.

Q: does it work when there is a mass matrix E?

A: no, and it would require significant changes to the algorithm to
support it. We'll be working on that hopefully.

References
==========

[MehP12]
:   Mehrmann, Volker; Poloni, Federico *Doubling algorithms with
    permuted Lagrangian graph bases.* SIAM J. Matrix Anal. Appl. 33
    (2012), no. 3, 780–805.

[MehP13]
:   Mehrmann, Volker; Poloni, Federico *Using permuted graph bases in H∞
    control.* Automatica J. IFAC 49 (2013), no. 6, 1790–1797.

[ChuLM07]
:   Chu, Delin; Liu, Xinmin; Mehrmann, Volker *A numerical method for
    computing the Hamiltonian Schur form.* Numer. Math. 105 (2007), no.
    3, 375–412.

[BenBMX07]
:   Benner, Peter; Byers, Ralph; Mehrmann, Volker; Xu, Hongguo *A robust
    numerical method for the γ-iteration in H∞ control.* Linear Algebra
    Appl. 425 (2007), no. 2-3, 548–570.

[carex]
:   Jörn Abels , Peter Benner *CAREX - A Collection of Benchmark Examples for Continuous-Time Algebraic Riccati Equations (Version 2.0) (1999)*,
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.40.4899


