"""
```julia
a,ā = coherent(β,N;dist)
```
Sample the phase-space distribution for a coherent state.

`β` is coherent (complex) amplitude.

`N` is number of samples.

`dist` is distrubtion. Can be either `W` or `+P` (default).

For a coherent state in +P, the distribution is just a point on the complex plane
at the position of the coherent amplitude.
"""
function coherent(β,N;dist="+P")
if dist=="+P"
    α = β*ones(N)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="W"
    α = β + (randn(N)+im*randn(N))/2
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
