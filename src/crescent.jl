"""
```julia
a,ā = crescent(β,ϵ,q,N;dist)
```
Sample phase space distribution for a crescent state

`β` is coherent amplitude

`ϵ` is squeeze paramter

`q` is shearing parameter

`N` is number of samples

`dist` is distribution. Can be either `Q` or `+P`
"""
function crescent(β,ϵ,q,N;dist="+P")
if dist=="+P"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    μ = β + (randn(N)*exp(-r)+im*randn(N)*exp(r))/sqrt(2)
    μ .*= exp.(im*q*randn(N))
    γ = (randn(N)+im*randn(N))/sqrt(2)
    α = μ + γ
    ᾱ = conj(μ - γ)
    return α,ᾱ
elseif dist=="Q"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β + (randn(N)*exp(-r)+im*randn(N)*exp(r))/sqrt(2)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α,ᾱ
elseif dist=="W"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β + 0.5*(randn(N)*exp(-r)+im*randn(N)*exp(r))*exp(-im*ϕ)
    α = α.*exp.(im*q*randn(N))
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
