# PhaseSpaceTools.jl
A julia package for working with quantum phase space distributions.

# Overview
This package provides sampling methods for commonly used quantum states in various quantum phase-space representations, including Glauber-P, Positive-P, HusimiQ, and Wigner distributions.

There are also convenience methods for calculating operator averages from phase-space averages, and for sampling noises for solving SDEs in [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

# Installation
In the julia REPL, enter

`julia> ]add https://github.com/AshtonSBradley/PhaseSpaceTools.jl`

`Pkg> test PhaseSpaceTools`

# Usage Examples

## Coherent state
A trivial example is provided by the coherent state $|\alpha\rangle$. These are eigenstate of the Bose annihilation operator $a$, with commutation relation $[a,a^\dagger]=1$. They satisfy the eigenvalue equation

```math
 {\hat a}|\alpha\rangle = \alpha|\alpha\rangle
```

<<<<<<< HEAD
As eigenstate of the Bose annihilation operator ``a`` with commutator ``[a,a^\dagger]=1``, satisfying

```math
 {\hat a}|\alpha\rangle = \alpha|\alpha\rangle,
=======
where $\alpha\in ℂ$. The Fock basis $|n\rangle$ is often convenient to work in. Fock states are eigenstates of the number operator $n=a^† a$:

```math
a^† a|n\rangle = n|n\rangle.
```

In the Fock basis, the coherent states take the form

```math
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty\frac{\alpha^n}{n!}|n\rangle.
>>>>>>> c8343b6cbae37cd25cc96f4e9d67dae17649af6c
```

Coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamical problems to equivalent stochastic differential equations.

### Glauber-P
In the `glauberP` representation the state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```
### Positive-P
In the `positiveP` representation the simples way to sample the state is again as a point on the complex plane

```julia
α,ᾱ = positiveP(state,N)
```

## Fock state
A more difficult state to sample is the eigenstate of the number operator...



# Quantum states
Quantum states available in `PhaseSpaceTools.jl` are

## Coherent
```@docs
Coherent
```
## Fock
```@docs
Fock
```
## Crescent
```@docs
Crescent
```
## Squeezed
```@docs
Squeezed
```
## Thermal
```@docs
Thermal
```
## Bogoliubov
```@docs
Bogoliubov
```

# Distributions
The phase-space distributions that may be sampled for these states are

## Husimi-Q
```@docs
husimiQ
```
## Glauber-P
```@docs
glauberP
```
## Positive-P
```@docs
positiveP
```
## Wigner
```@docs
wigner
```
## Positive-W
```@docs
positiveW
```

# Noises

## randnc
```@docs
randnc
```

## realnoise
```@docs
realnoise
```
## realbridge
```@docs
realbridge
```

# Solving SDEs

## Complex noise example

## Real noise example

# Recovering normal order
A common approach to quantum phase space simulations involves performing truncated Wigner simulations, necessitating symmetric order of operator moments. We can easily recover normal order, the experimentally relevant form, using

```math
\langle (a^\dagger)^p a^q\rangle = \sum_{n=0}^{\textrm{min}(p,q)}
```

# Citing

# Index

```@index
```
