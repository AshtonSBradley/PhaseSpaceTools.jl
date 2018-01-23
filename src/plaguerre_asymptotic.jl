"""
    plaguerre_asymptotic(x,n)

Define the radial phase-space probability distribution for sampling a fock state in the
`+W` representation. In the large `n` limit the distribution is sampled using the ansymptotic expansion
```math
P(x,n)=\\frac{1}{\\sqrt{\\pi}}\\exp{(-[x-\\sqrt{n+1}]^2)}
```
where
`x`: point to evaluate probability.

`n`: number of the fock state |n⟩.

"""
function plaguerre_asymptotic(x,n)
return exp(-(x-sqrt(n+1))^2)/sqrt(pi)
end
