"""
```julia
P(x,n) = plaguerre_asymptotic(x,n)
```
Radial phase space probability distribution for fock state sampling in
+P representation, in the large `n` limit - asymptotic expansion.

`x ` is location to evaluate probability.

`n` is the number of the fock state.

"""
function plaguerre_asymptotic(x,n)
return exp(-(x-sqrt(n+1))^2)/sqrt(pi)
end
