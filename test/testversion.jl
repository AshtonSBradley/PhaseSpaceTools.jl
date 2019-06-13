#test
using Pkg
Pkg.activate(".")

using Test, LinearAlgebra, Revise, PhaseSpaceTools

# Coherent
N = 100000
α = 100
state = Coherent(α)

a,ā = wigner(state,N)
meana = mean(a)
n̄ = real(mean(@. ā*a)-.5)
Vn = mean(@. ā^2*a^2)-mean(@. a*ā)-n̄^2 |> abs

#test
@test norm(meana - α)/abs(α) < 0.01
@test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
@test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05

# scatter(a)

# Fock
n = 100
N = 100000
state = Fock(n)
a,ā = wigner(state,N)

meana = mean(a)
absa = abs(meana)
n̄ = mean(abs2.(a))-.5
Vn= abs(mean(abs2.(a).^2)-mean(abs2.(a))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

#test
@test absa < 0.1
@test n̄ - n < 0.01
@test Vn < 0.01
@test rel_num_var < 0.001

# scatter(a)


n = 99
N = 100000
state = Fock(n)
a,ā = positiveW(state,N)

meana = mean(a)
absa = abs(meana)
n̄ = real(mean(a.*ā)-.5)
Vn = abs(mean(a.*a.*ā.*ā)-abs(mean(a.*ā))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄)

#test
@test absa < 0.1
@test n̄ - n < 1
@test Vn < 50
@test rel_num_var < 0.1

# scatter(a)

#sample
n = 101
state = Fock(n)
N = 100000

a,ā = positiveW(state,N)

av_a = mean(a)
absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄)

#test
@test absa < 0.1
@test (n̄ - n)/n < 0.01
@test Vn < 30
@test rel_num_var < 0.1

# scatter(a)

#sample
β = 10
ϕ = π/16
r = 2
ϵ = r*exp(2*im*ϕ)
state = Squeezed(β,ϵ)
N = 10000

a,ā = positiveP(state,N)

av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā))
nbar = sinh(abs(ϵ)).^2 + abs2(β)

#test
@test norm(av_a - β)/abs(β) < 0.01
@test abs(n̄ - nbar)/abs(β)^2 < 0.01

# scatter(a)

# Squeezed
β = 10
ϕ = π/16
r = 1.5
ϵ = r*exp(2*im*ϕ)
state = Squeezed(β,ϵ)
N = 10000
a,ā = wigner(state,N)

av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
nbar = sinh(abs(ϵ)).^2+abs2(β)

#test
@test norm(av_a - β)/abs(β) < 0.01
@test abs(n̄ - nbar)/abs(β)^2 < 0.01

# scatter(a)

# Bogoliubov

#sample
N = 100000
n̄ = 10

function randuv()
    u,v = crandn(2)
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
abs2(u)+abs2(v)

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
