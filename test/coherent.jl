"""
    a,ā = coherent(β,N;dist=:posP)

Sample the phase-space distribution for a coherent state.

`β`: coherent (complex) amplitude.

`N`: number of samples.

`dist`: phase-space distribution, either `:W` or `:posP`.

For a coherent state in +P, the distribution is just a point on the complex plane
at the position of the coherent amplitude.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.
"""
function coherent(β,N;dist=:posP)
if dist==:posP
    α = β*ones(N)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist==:W
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist==:posW
    α = β .+ randnc(N)/sqrt(2)
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
