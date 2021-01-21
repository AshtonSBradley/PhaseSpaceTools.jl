### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ ccccd3d4-5b5c-11eb-2d48-8593366fb942
using PhaseSpaceTools, Plots, LaTeXStrings

# ╔═╡ b45b5d84-5b5c-11eb-2189-371f4b3bf28a
html"<button onclick=present()>Present</button>"

# ╔═╡ 9535e4b6-9308-4077-a6c6-f5860a8da4c2
md"""
# Sampling quantum states
Sampling single mode quantum states for use in quantum phase space simulations using the truncated Wigner and +P distributions. Similar methods can be applied to multimode simulations of Bose-Einstein condensates.

We use the methods presented in the article
- [Olsen, Bradley, Opt. Comm. 282 (2009) 3924-3929](https://doi.org/10.1016/j.optcom.2009.06.033)

with erratum [Olsen, Lewis-Swan, Bradley, Opt. Comm. 370 (2016) 327-328](https://doi.org/10.1016/j.optcom.2016.02.068)

A complete erratum can be found at [SciPost](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/).

# New States
We have add sampling for the additional distributions:
- Fock state for +W
- Thermal state for +P
- Bogoliubov states for P,Q,W

# Papers
For more information on Bose-Einstein condensate applications, see the review 

_Dynamics and statistical mechanics of ultra-cold Bose gases using c-field techniques_, P. B. Blakie, A. S. Bradley, M. J. Davis, R. J. Ballagh, and C. W. Gardiner, [Advances in Phyiscs 57, 363 (2008)](http://dx.doi.org/10.1080/00018730802564254)

For a recent paper using this package, see

_Steady states, squeezing, and entanglement in intracavity triplet down conversion_,
Mathew D. E. Denys, Murray K. Olsen, Luke S. Trainor, Harald G. L. Schwefel, and Ashton S. Bradley,  [Optics Communications 484 (2021) 126699](https://doi.org/10.1016/j.optcom.2020.126699), [arXiv](https://arxiv.org/abs/1907.09572)
"""

# ╔═╡ 2a3a9c5e-5b67-11eb-2250-41d162813a4a
md"# Examples"

# ╔═╡ bf17f9ca-8dcf-483d-a399-42bdb69f79f9
gr(size=(350,350),xlabel=L"\alpha_r",ylabel=L"\alpha_i",
	    xlims=(-20,20),ylims=(-20,20),ms=.2,
	    legend=false,grid=false,aspect_ratio=1,colorbar=true)

# ╔═╡ adce375f-31dd-4bbd-896a-6bdf823a4456
function show_meansW(N,State,p,Dist)
	s = State(p...)
	a,ā = Dist(s,N)
    n̄ = mean(a.*ā)-.5
    Vn= mean(a.^2 .*ā.^2)-mean(a.*ā)-n̄.^2
    Text("
	N = $N samples.\n  
    Averages:\n
	<â> = $(mean(a))\n
	<â⁺â> = $(n̄)\n
	V(n̂) = $(Vn)\n
	Relative width = $(sqrt(abs(Vn))/abs(n̄))")
end

# ╔═╡ eba146d8-5b5e-11eb-363c-1d84d1f72610
function sample_plot(N,State,p,Dist)
	s = State(p...)
	a,ā = Dist(s,N)
	scatter(real(a),imag(a),c=:red)
end

# ╔═╡ 86a58536-5b66-11eb-2e94-0f42924c2bee
md"# Coherent State
Let's look at the coherent state in different representations."

# ╔═╡ 2421fc20-5b5f-11eb-35c0-253c6b5420f6
sample_plot(500,Coherent,10,wigner)

# ╔═╡ e5fedb01-83be-4fe0-9072-be338998707a
show_meansW(100000,Coherent,10,wigner) # Poisson number statistics

# ╔═╡ 7b1a6556-5b66-11eb-21c0-53b40b983fab
md"# Thermal State"

# ╔═╡ 9d51cea0-5b5e-11eb-17f8-b9ec85cfd39a
sample_plot(500,Thermal,(0,10),wigner)

# ╔═╡ 130e58fa-72b1-4005-8e3b-c624da6db44f
show_meansW(100000,Thermal,(0,10),wigner)

# ╔═╡ 3a3f2a44-5b66-11eb-02df-ad97df0d48d0
md"# Squeezed state"

# ╔═╡ 10c936c9-452a-4c05-a5d7-6e2d1735d3c2
begin
	β = 10 
	ϕ = 0.
	r = 1
	ϵ = r*exp(2*im*ϕ)
	sample_plot(500,Squeezed,(β,ϵ),wigner)
end

# ╔═╡ b163aaf8-47c2-4d87-a92f-8703eb9f2d6d
show_meansW(100000,Squeezed,(β,ϵ),wigner)

# ╔═╡ 2c9b2906-5b66-11eb-2eb8-932c97b5e2b0
md"# Fock state"

# ╔═╡ 472762c2-5b60-11eb-361d-37601eaed6ff
sample_plot(500,Fock,100,wigner)

# ╔═╡ 37f1a96e-5b62-11eb-217f-9dc0fae96dd6
show_meansW(100000,Fock,100,wigner)

# ╔═╡ 4c12caae-5b62-11eb-1ad2-5fc5d8120fba
sample_plot(1000,Fock,100,positiveP)

# ╔═╡ c4a573b8-5b62-11eb-0726-3d42d1187b64
sample_plot(1000,Fock,100,positiveW)

# ╔═╡ e0583d34-5b62-11eb-10a8-f380c764daa7
sample_plot(1000,Fock,321,positiveP) #uses asymptotic expansion for n>320

# ╔═╡ 2df01300-5b63-11eb-0faf-9d8fb6d796c9
md"# Crescent state
A Kerr interaction will create shear for a coherent state. This can be modelled using a crescent state."

# ╔═╡ 3c8a60b4-5b63-11eb-13ef-037243089ec1
sample_plot(500,Crescent,(10,0,.5),wigner)

# ╔═╡ 6d01bc24-5b63-11eb-2fea-a3e993675754
sample_plot(500,Crescent,(10,0,.5),husimiQ)

# ╔═╡ 7befd068-5b63-11eb-15d8-a19df4d9da1b
sample_plot(500,Crescent,(10,0,.5),positiveP)

# ╔═╡ Cell order:
# ╟─b45b5d84-5b5c-11eb-2189-371f4b3bf28a
# ╟─9535e4b6-9308-4077-a6c6-f5860a8da4c2
# ╟─2a3a9c5e-5b67-11eb-2250-41d162813a4a
# ╠═ccccd3d4-5b5c-11eb-2d48-8593366fb942
# ╠═bf17f9ca-8dcf-483d-a399-42bdb69f79f9
# ╠═adce375f-31dd-4bbd-896a-6bdf823a4456
# ╠═eba146d8-5b5e-11eb-363c-1d84d1f72610
# ╟─86a58536-5b66-11eb-2e94-0f42924c2bee
# ╠═2421fc20-5b5f-11eb-35c0-253c6b5420f6
# ╠═e5fedb01-83be-4fe0-9072-be338998707a
# ╟─7b1a6556-5b66-11eb-21c0-53b40b983fab
# ╠═9d51cea0-5b5e-11eb-17f8-b9ec85cfd39a
# ╠═130e58fa-72b1-4005-8e3b-c624da6db44f
# ╟─3a3f2a44-5b66-11eb-02df-ad97df0d48d0
# ╠═10c936c9-452a-4c05-a5d7-6e2d1735d3c2
# ╠═b163aaf8-47c2-4d87-a92f-8703eb9f2d6d
# ╟─2c9b2906-5b66-11eb-2eb8-932c97b5e2b0
# ╠═472762c2-5b60-11eb-361d-37601eaed6ff
# ╠═37f1a96e-5b62-11eb-217f-9dc0fae96dd6
# ╠═4c12caae-5b62-11eb-1ad2-5fc5d8120fba
# ╠═c4a573b8-5b62-11eb-0726-3d42d1187b64
# ╠═e0583d34-5b62-11eb-10a8-f380c764daa7
# ╟─2df01300-5b63-11eb-0faf-9d8fb6d796c9
# ╠═3c8a60b4-5b63-11eb-13ef-037243089ec1
# ╠═6d01bc24-5b63-11eb-2fea-a3e993675754
# ╠═7befd068-5b63-11eb-15d8-a19df4d9da1b
