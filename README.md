# PhaseSpaceTools

[![Build Status](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl.svg?branch=master)](https://travis-ci.org/AshtonSBradley/PhaseSpaceTools.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t6i7kdnpffgnq6pg?svg=true)](https://ci.appveyor.com/project/AshtonSBradley/phasespacetools-jl)
[![Coverage Status](https://coveralls.io/repos/AshtonSBradley/PhaseSpaceTools.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl/coverage.svg?branch=master)](http://codecov.io/github/AshtonSBradley/PhaseSpaceTools.jl?branch=master)

Small package for sampling some of the quantum initial states commonly encountered in quantum-optical and matter-wave bosonic systems. Wigner (`W`) and positive-P (`+P`) representations are implemented, being the most useful for dynamical simulations. Currently supports only single mode sampling.

## Install

```julia
] add https://github.com/AshtonSBradley/PhaseSpaceTools.jl.git
using PhaseSpaceTools
```

## Usage

### States
`coherent`, `thermal`, `squeezed`, `fock`, `crescent`

#### Coherent state
A coherent state |α⟩ is sampled as
```julia
α = 1.0+im*2.0 #coherent state |α⟩
N = 1000 #number of samples
a,ā = coherent(α,N,dist="+P")
```
This is a special (trivial) case where the two phase space variables `a` and `ā` are complex conjugate, and non-stochastic in the `+P` representation.

#### Fock state
An approximate fock state sampler in the Wigner representation:
```julia
n = 100 #fock state |n⟩
N = 1000 #number of samples
a,ā = fock(n,N,dist="W")
```
Provides a positive `W` approximation that reproduces moments for large `n`.

## Examples

See  `/examples/PhaseSpaceTools.ipynb` for more usage.

## External links
[Numerical representation of quantum states in the positive-P and Wigner representations, M K Olsen, A S Bradley, Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

[Scipost Commentary with full erratum](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
