#test +W for large n

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
