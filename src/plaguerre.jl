function plaguerre(x,n)
    #helper function for +W Fock state
#P(x,n)=(4/3)*exp(-2*x^2/3)*x*laguerre(-4*x^2/3,n)/(3^n)
    return (4/3)*exp(-2*x^2/3 + log(x) + log(laguerre(-4*x^2/3,n)) - n*log(3))
end
