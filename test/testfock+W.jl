#test +W

n = 99
N = 100000
a,ā = fock(n,N;dist=:posW)

av_a = mean(a)
absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-abs(mean(a.*ā))-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);

@test absa < 0.1
@test n̄ - n < 1
@test Vn < 50
@test rel_num_var < 0.1
