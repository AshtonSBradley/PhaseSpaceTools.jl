# Bogoliubov
N = 100000
n̄ = 10

function randuv()
    u,v = randnc(2)
    nrm = abs2(u)-abs2(v)
    if nrm < 0
         new = [0 1;1 0]*[u; v]
         u,v = new
     end
    u /=sqrt(abs(nrm))
    v /=sqrt(abs(nrm))
    return u,v
end

u,v = randuv()
abs2(u)-abs2(v) ≈ 1.0

#thermal state
state = Bogoliubov(u,v,n̄)
a,ā = wigner(state,N)

# particle mode population
na = real(mean(a.*ā))-0.5
nth = (abs2(u)+abs2(v))*(n̄+0.5) - 0.5
@test isapprox(na,nth,rtol=1e-2)


#test vacuum limit
state = Bogoliubov(u,v,0)
a,ā = wigner(state,N)

na = real(mean(a.*ā))-0.5
nvac = abs2(v)
@test isapprox(na,nvac,atol=5e-2)
