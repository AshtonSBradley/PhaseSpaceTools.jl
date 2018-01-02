"""
```julia
a,ā = coherent(β,N;dist)
```
samples phase space distribution for coherent state

`β` is coherent amplitude (complex)

`N` is number of samples

`dist` is distrubtion. Can be either `+P` or `W`

For a coherent state in +P the distribution is just a point on the complex plane
at the position of the coherent amplitude.
Default (no value for dist) will give +P
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
