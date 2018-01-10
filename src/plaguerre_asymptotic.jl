"""
```julia
P(x,n) = plaguerre_asymptotic(x,n)
```
Define the radial phase-space probability distribution for sampling a fock state in the
+P representation. In the large `n` limit the distribution is sampled using an ansymptotic expansion.

`x ` is the point to evaluate probability.

`n` is the number of the fock state |n⟩.

"""
function plaguerre_asymptotic(x,n)
return exp(-(x-sqrt(n+1))^2)/sqrt(pi)
end
