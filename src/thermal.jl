"""
```julia
a,ā = thermal(β,n̄,N;dist)
```
samples phase space distribution for a thermal state:

`β` is complex displacement

`n̄` is thermal average population

`N` is number of samples

`dist` is the distribution. Can be `P`, `Q` or `W`
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
