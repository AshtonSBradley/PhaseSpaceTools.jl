#test W

n = 101
N = 100000
a,ā = fock(n,N;dist="+W")

av_a = mean(a)
absa = abs(av_a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2)
rel_num_var = sqrt(abs(Vn))/abs(n̄);
