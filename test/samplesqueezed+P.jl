β = 10
ϕ = π/16
r = 2
ϵ = r*exp(2*im*ϕ)
N = 10000
a,ā = squeezed(β,ϵ,N;dist=:posP)

av_a = mean(a);absa = abs(av_a)
n̄ = real(mean(a.*ā))
nbar = sinh(abs(ϵ)).^2+abs2(β)
