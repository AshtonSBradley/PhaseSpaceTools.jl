# PhaseSpaceTools

[![Build Status](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t6i7kdnpffgnq6pg?svg=true)](https://ci.appveyor.com/project/AshtonSBradley/phasespacetools-jl)
[![Coverage Status](https://coveralls.io/repos/AshtonSBradley/PhaseSpaceTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)

Package for sampling some of the quantum initial states commonly encountered in quantum-optical and matter-wave bosonic systems. W and +P representations are implemented, being the most useful for dynamical simulations. Currently supports only single mode sampling.

## Install

```julia
Pkg.clone("https://github.com/AshtonSBradley/PhaseSpaceTools.jl.git")
Pkg.test("PhaseSpaceTools")
using PhaseSpaceTools
```
#### Coherent state
```julia
julia>
α = 1.0+im*2.0 #coherent state |α⟩
N = 1000 #number of samples
a,ā = coherent(α,N,dist="+P")
```
a special (trivial) case where the two phase space variables a and ā are complex conjugate in the +P representation.

#### Fock state
An approximate fock state sampler for W:
```julia
julia>
n = 100 #fock state |n⟩
N = 1000 #number of samples
a,ā = fock(n,N,dist="W")
```
provides a positive W approximation for large `n`.

See  `/examples/PhaseSpaceTools.ipynb` for more usage.
#### External links
[Numerical representation of quantum states in the positive-P and Wigner representations, Olsen, Bradley, Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)
