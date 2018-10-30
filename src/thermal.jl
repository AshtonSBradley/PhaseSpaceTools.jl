"""
    a,ā = thermal(β,n̄,N;dist="P")

Sample the phase-space distribution for a thermal state.

`β`: coherent (complex) displacement.

`n̄`: thermal average population.

`N`: number of samples.

`dist`: phase-space distribution. Can be `P`, `Q` or `W`.

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.
"""
function thermal(β,n̄,N;dist="P")
if dist=="P"
    α = β .+ sqrt(n̄)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="W"
    α = β .+ sqrt(n̄+.5)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
elseif dist=="Q"
    α = β .+ sqrt(n̄+1.0)*crandn(N)
    ᾱ = conj(α)
    return α, ᾱ
else error("distribution not implemented")
end
end
