#helper functions
"""
`a = crandn(N)`

Retunrs a vector of length `N` samples of complex random variates with mean zero and variance one:
```math
\\langle |a|^2\\rangle = 1; \\quad\\quad\\langle a^2\\rangle = \\langle (a^*) ^2\\rangle = 0,
```
Useful for creating a range of noises that show up in `CField` simulations.
"""
function crandn(N)
    return (randn(N).+im*randn(N))/sqrt(2)
end
