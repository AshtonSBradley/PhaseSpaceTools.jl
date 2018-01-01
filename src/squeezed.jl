"""
```julia
a,ā = squeezed(β,ϵ,N;dist)
```
samples phase space distribution for a squeezed state

`β` is coherent amplitude

`ϵ` is (complex) squeezing parameter

`N` is number of samples

`dist` is the distribution; can be `W` or `+P`
"""
function squeezed(β,ϵ,N;dist="+P")
if dist=="+P"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    γ = (randn(N)+im*randn(N))/sqrt(2)
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(N)+im*sqrt(exp(r)*cosh(r)/2)*randn(N)
    α = β + exp(im*ϕ)*ν + γ
    ᾱ = conj(β) + exp(-im*ϕ)*conj(ν) - conj(γ)
    return α, ᾱ
elseif dist=="W"
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β + 0.5*(randn(N)*exp(-r)+im*randn(N)*exp(r))*exp(-im*ϕ)
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
