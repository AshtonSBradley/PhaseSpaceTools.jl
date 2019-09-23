# Solving SDEs

## Noises

### randnc
```@docs
randnc
```

### realnoise
```@docs
realnoise
```
### realbridge
```@docs
realbridge
```

## Complex noise example

## Real noise example

## Recovering normal order from Wigner ensemble averages
A common approach to quantum phase space simulations involves performing truncated Wigner simulations, allowing _symmetrically ordered_ operator averages.
However, in experiments, the physically measurable moments are usually those in _normal order_. We can recover normal order for an operator average of a given degree if we have access to all symmetrically ordered operator averages up to the same degree.

One can derive the following useful relationship between normal- and symmetrically-ordered operator averages

```math
\langle (a^\dagger)^p a^q\rangle = \sum_{n=0}^{\textrm{min}(p,q)}\frac{(-1)^nn!}{2^n}\binom{p}{n}\binom{q}{n}\langle \{(a^\dagger)^{p-n}a^{q-n}\}_\textrm{sym}\rangle
```

The relationship between normal ordered operator averages and ensemble averages of Wigner trajectories, denoted by $\overline{\left(\dots\right)}_W$ is

```math
\langle (a^\dagger)^p a^q\rangle = \sum_{n=0}^{\textrm{min}(p,q)}\frac{(-1)^nn!}{2^n}\binom{p}{n}\binom{q}{n}\overline{\left((\alpha^*)^{p-n}\alpha^{q-n}\right)}_W
```
