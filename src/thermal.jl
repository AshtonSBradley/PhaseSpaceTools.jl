"""
```julia
a,ā = thermal(β,n̄,N;dist)
```
Sample the phase-space distribution for a thermal state.

`β` is the coherent (complex) displacement.

`n̄` is the thermal average population.

`N` is the number of samples.

`dist` is the distribution. Can be `P` (default), `Q` or `W`
"""
function thermal(β,n̄,N;dist="P")
if dist=="P"
    α = β .+ (randn(N)+im*randn(N))*sqrt(n̄/2)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="W"
    α = β .+ (randn(N)+im*randn(N))*sqrt((n̄+.5)/2)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="Q"
    α = β .+ (randn(N)+im*randn(N))*sqrt((n̄+1.0)/2)
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
