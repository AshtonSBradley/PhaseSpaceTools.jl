"""
```julia
P(x,n) = plaguerre(x,n)
```
Radial phase space probability distribution for fock state sampling in
+P representation.

`x ` is location to evaluate probability.

`n` is the number of the fock state.

"""
function plaguerre(x,n)
    #helper function for +W Fock state
#P(x,n)=(4/3)*exp(-2*x^2/3)*x*laguerre(-4*x^2/3,n)/(3^n)
    return (4/3)*exp(-2*x^2/3 + log(x) + log(laguerren(-4*x^2/3,n)) - n*log(3))
end
