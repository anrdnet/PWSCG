PWSCG
=====
This project contains an implementation of the PWSCG finite element method.

Plane Wave Semi-Continuous Galerkin is a finite element method using
complex exponential plane wave functions as basis functions instead of the
conventional polynomial basis functions.

Method is described in [this thesis](http://urn.nb.no/URN:NBN:no-49625).

### Dependencies:
 - GetFEM++ (home.gna.org/getfem/)
 - scons (for building) (www.scons.org)

### Build:
```
$ scons
```
Main binary will be available as `build/pwscg`.

### Usage:
```
$ pwcg {-p} {-z} {-n <subdiv>} {-s {<n>}+}+ -{k{,1,2,3},f,u}{r,i} <expr>

-p      Use linear polynomials instead of plane wave basis functions.
-z      Skip computation of error norms.
-n      Number of times to subdivide grid before output.
-s      Followed by a sequence of numbers, specifies the grid size of the
        simulation. Can be specified multiple times to run simulation several
        times on different grids.

-{k{,1,2,3},f,u}{r,i}
        Set the expression for `k`, `k1`, `k2`, `k3`, `f` or `u`. If suffix is
        `r` sets real part, if suffix is `i` sets imaginary part.
```

Any questions regarding this project can be directed to a@anrd.net
