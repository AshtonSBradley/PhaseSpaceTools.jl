"""
    `a,ā = bogoliubov(u,v,n̄,N;dist="W")`

Sample the phase-space distribution for a Bogoliubov mode in the Wigner representation.

`u,v`: Bogoliubov amplitudes, with thermal population `n̄`.

`N`: number of samples.

`dist`: phase-space distribution, `W`.

Draws samples from the state

```math
\\alpha = \\beta u - \\beta^* v
```

where \$\beta\$, \$\beta^*\$ are sampled as thermal states with

```math
\\langle |\\beta|^2\\rangle = \bar n +\frac{1}{2}.
```

For standard `P,Q,W` distributions, `a` and `ā` are complex conjugate, while for `+P` etc,
`a` and `ā` are independent variables.
"""
function bogoliubov(u,v,n̄,N;dist="W")
    if dist=="W"
        b,b̄ = thermal(0.,n̄,N;dist="W")
        a = u*b + conj(v)*b̄ 
        ā = conj.(a)
    else error("distribution not implemented")
    end
    return a,ā
end
