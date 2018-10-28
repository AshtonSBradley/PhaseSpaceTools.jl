N = 100000
u = crandn(1)[1]
v = crandn(1)[1]

n̄ = 307
a,ā = bogoliubov(u,v,n̄,N)

#test thermal mode population (zero coherent amplitude)
N̄ = real(mean(a.*ā))

#should agree with analytic form for quantum state:
nbog = (abs2.(u) + abs2.(v))*(n̄+.5)

isapprox(N̄,nbog,rtol=1e-2)
