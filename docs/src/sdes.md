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

## Recovering normal order
A common approach to quantum phase space simulations involves performing truncated Wigner simulations, necessitating symmetric order of operator moments. We can easily recover normal order, the experimentally relevant form, using

```math
\langle (a^\dagger)^p a^q\rangle = \sum_{n=0}^{\textrm{min}(p,q)}
