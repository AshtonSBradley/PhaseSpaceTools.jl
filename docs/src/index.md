# PhaseSpaceTools.jl
A julia package for working with quantum phase space distributions.

# Overview
This package provides sampling methods for commonly used quantum states in various quantum phase-space representations, including Glauber-P, Positive-P, Husimi-Q, and Wigner distributions.

There are also convenience methods for calculating operator averages from phase-space averages, and for sampling noises for solving SDEs using [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

# Installation
In the julia REPL

`julia> ]add https://github.com/AshtonSBradley/PhaseSpaceTools.jl`

`Pkg> test PhaseSpaceTools`

# Usage Examples

## Coherent state
A trivial example is provided by the coherent state $|\alpha\rangle$. These are eigenstate of the Bose annihilation operator ``a``, with commutation relation ``[a,a^\dagger]=1``. They satisfy the eigenvalue equation

```math
 a|\alpha\rangle = \alpha|\alpha\rangle
```

The Fock basis ``|n\rangle`` is often convenient to work in. Fock states are eigenstates of the number operator ``n=a^\dagger a``:

```math
a^\dagger a|n\rangle = n|n\rangle.
```

In the Fock basis, the coherent states take the form

```math
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty\frac{\alpha^n}{n!}|n\rangle.
```

Coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamics to equivalent stochastic differential equations.

### Glauber-P
In the Glauber-P representation the state may be sampled as a single point on the complex plane

```julia
N = 10000
state = Coherent(12.0)
α,ᾱ = glauberP(state,N)
```

In any other distribution the coherent state is more interesting to sample.
### Positive-P
In the Positive-P representation the simples way to sample the state is again as a point on the complex plane

```julia
α,ᾱ = positiveP(state,N)
```

Here ``α,ᾱ`` are once again complex conjugates.
However, we emphasize that this is no longer the case in any subsequent dynamical evolution of the ensemble of phase space points. Furthermore, unlike the Glauber-P distribution, the positive-P distribution for a particular quantum state is not unique.

In this package, we treat single mode problems, taking the general approach of sampling the Husimi-Q funcion, and then exploiting the simple convolutional relationship between the Husimi-Q and the positive-P representations:

```math
P(\alpha,\bar\alpha) = \frac{1}{4\pi^2}\langle (\alpha+\bar\alpha^*)/2|\rho|(\alpha+\bar\alpha^*)/2\rangle e^{-|\alpha-\bar\alpha^*|^2/4}
```
We use the transformation
```math
\mu = (\alpha+\bar\alpha^*)/2,\quad\quad\gamma = (\alpha-\bar\alpha^*)/2
```
with inverse
```math
\alpha = \mu + \gamma,\quad\quad \bar\alpha = \mu^* - \gamma^*.
```
The problem reduces to that of sampling the  distribution
```math
P(\mu,\gamma)=Q(\mu,\mu^*)e^{-|gamma|^2}{\pi}
```
where
```math
Q(\alpha,\alpha^*)\equiv \frac{\langle \alpha|\rho|\alpha\rangle}{\pi}
```
is the Husimi-Q function, and the remaining Gaussian is readily sampled.
## Fock state
A more difficult state to sample an eigenstate of the number operator.



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
