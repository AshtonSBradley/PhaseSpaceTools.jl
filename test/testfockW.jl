#test W

n = 100
N = 100000
a,ā = fock(n,N;dist=:W)

av_a = mean(a)
absa = abs(av_a)
n̄ = mean(abs2.(a))-.5
Vn= abs(mean(abs2.(a).^2)-mean(abs2.(a))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

@test absa < 0.1
@test n̄ - n < 0.01
@test Vn < 0.01
@test rel_num_var < 0.001
