#Squeezed
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
@test norm(av_a - β)/abs(β) < 0.05
@test abs(n̄ - nbar)/abs(β)^2 < 0.05
