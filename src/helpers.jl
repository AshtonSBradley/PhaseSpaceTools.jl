#helper functions
"""
`a = crandn(N...)`

Returns an array of dimension specified by the length of the tuple `N`, containing samples of complex random variates with mean zero and variance one:
```math
\\langle |a|^2\\rangle = 1; \\quad\\quad\\langle a^2\\rangle = \\langle (a^*) ^2\\rangle = 0,
```
Useful for creating a range of noises that show up in `CField` simulations.

Usage:

`a = crandn(10)` returns a `10-element Array{Complex{Float64},1}`.

`a = crandn(50,100)` returns a `50x100-element Array{Complex{Float64},2}`.
"""
function crandn(args...)
    return randn(ComplexF64,args...)
end
