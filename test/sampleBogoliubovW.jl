N = 100000
u = crandn(1)[1]
v = crandn(1)[1]

n̄ = 307
a,ā = bogoliubov(u,v,n̄,N)

# thermal mode population (zero coherent amplitude)
N̄ = real(mean(a.*ā))

# analytic form for Bogoliubov state:
nbog = (abs2.(u) + abs2.(v))*(n̄+.5)
