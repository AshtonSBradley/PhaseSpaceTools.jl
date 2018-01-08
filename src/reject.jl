"""
```julia
samples = reject(P,w,N,Pmax)
```
Samples the distribution `P` over the window

`w=[w1,w2]`, where `P` is a function with `max(P)<=Pmax`.

`N` is number of samples

"""
function reject(P,w,N,Pmax)
    #Rejection sampling the probability distribution P over the window w=[w1,w2]
    #Pmax must be numerical with value max(P) <= Pmax

    samples = Array{Float64}(1)
    while length(samples) < N + 1
        y = w[1] + rand()*(w[2]-w[1])
        z = rand()*Pmax
        if z < P(y)
            push!(samples,y)
        end
    end
    return samples[2:end]
end
