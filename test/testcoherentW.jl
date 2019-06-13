#test W

α = 100*(randn()+im*randn())
N = 100000
a,ā = coherent(α,N;dist=:W)

av_a = mean(a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2);

@test norm(av_a - α)/abs(α) < 0.01
@test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
@test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05
