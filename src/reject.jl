"""
```julia
samples = reject(P,w,N,Pmax)
```
Sample a probability distribution by rejection sampling over a finite window.

`P` is the probability distribution, normalized to 1.

`w=[w1,w2]` is the window for sampling.

`N` is number of samples.

`Pmax` is a numerical upper bound for `P`: `max(P) ≦ Pmax`.



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
