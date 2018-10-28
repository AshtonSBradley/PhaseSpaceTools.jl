#sample Bogoliubov state

N = 100000
(u,v) = (0.1,0.4)
n̄ = 102
a,ā = bogoliubov(u,v,n̄,N)

av_a = mean(a)
n̄ = real(mean(a.*ā)-.5)
Vn= abs(mean(a.*a.*ā.*ā)-mean(a.*ā)-n̄.^2);

#test bog mode statistics 
abs2.(u)*n̄ + abs2.(v)*(n̄+1)
