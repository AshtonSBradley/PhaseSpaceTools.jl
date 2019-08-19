var documenterSearchIndex = {"docs":
[{"location":"#PhaseSpaceTools.jl-1","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"A julia package for working with quantum phase space distributions.","category":"page"},{"location":"#Contents-1","page":"PhaseSpaceTools.jl","title":"Contents","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"","category":"page"},{"location":"#Overview-1","page":"PhaseSpaceTools.jl","title":"Overview","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"The main aim of this package is to provide sampling methods for commonly used quantum states in various quantum phase-space representations, including Glauber-P, Positive-P, HusimiQ, and Wigner distributions. There are also convenience methods for calculating operator averages from phase-space averages, and for sampling noises for solving SDEs in DifferentialEquations.jl.","category":"page"},{"location":"#Installation-1","page":"PhaseSpaceTools.jl","title":"Installation","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"In the julia REPL, enter","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"]add https://github.com/AshtonSBradley/PhaseSpaceTools.jl","category":"page"},{"location":"#Usage-Examples-1","page":"PhaseSpaceTools.jl","title":"Usage Examples","text":"","category":"section"},{"location":"#Coherent-state-1","page":"PhaseSpaceTools.jl","title":"Coherent state","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"A trivial example is provided by the coherent state","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"alpharangle = e^-alpha^22sum_n=0^inftyfracalpha^nnnrangle","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"As eigenstate of the Bose annihilation operator","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":" hat aalpharangle = alphaalpharangle","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"coherent states play a central role in quantum phase space methods, providing a mapping of many-body boson dynamical problems to equivalent stochastic differential equations.","category":"page"},{"location":"#Glauber-P-1","page":"PhaseSpaceTools.jl","title":"Glauber-P","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"In the glauberP representation the state may be sampled as a single point on the complex plane","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"N = 10000\nstate = Coherent(12.0)\na,ā = glauberP(state,N)","category":"page"},{"location":"#Positive-P-1","page":"PhaseSpaceTools.jl","title":"Positive-P","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"In the positiveP representation the simples way to sample the state is again as a point on the complex plane","category":"page"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"a,ā = positiveP(state,N)","category":"page"},{"location":"#Fock-state-1","page":"PhaseSpaceTools.jl","title":"Fock state","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"A more difficult state to sample is the eigenstate of the number operator","category":"page"},{"location":"#GlauberP-1","page":"PhaseSpaceTools.jl","title":"GlauberP","text":"","category":"section"},{"location":"#positiveP-1","page":"PhaseSpaceTools.jl","title":"positiveP","text":"","category":"section"},{"location":"#Quantum-states-1","page":"PhaseSpaceTools.jl","title":"Quantum states","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Quantum states available in PhaseSpaceTools.jl are","category":"page"},{"location":"#Coherent-1","page":"PhaseSpaceTools.jl","title":"Coherent","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Coherent","category":"page"},{"location":"#PhaseSpaceTools.Coherent","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Coherent","text":"Coherent(β)\n\nCreate a coherent state with complex amplitude β.\n\n\n\n\n\n","category":"type"},{"location":"#Fock-1","page":"PhaseSpaceTools.jl","title":"Fock","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Fock","category":"page"},{"location":"#PhaseSpaceTools.Fock","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Fock","text":"Fock(n)\n\nCreate a Fock state for particle number n.\n\n\n\n\n\n","category":"type"},{"location":"#Crescent-1","page":"PhaseSpaceTools.jl","title":"Crescent","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Crescent","category":"page"},{"location":"#PhaseSpaceTools.Crescent","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Crescent","text":"Crescent(β,ϵ,q)\n\nCreate a Cresecent state with parameters β, ϵ, q.\n\n\n\n\n\n","category":"type"},{"location":"#Squeezed-1","page":"PhaseSpaceTools.jl","title":"Squeezed","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Squeezed","category":"page"},{"location":"#PhaseSpaceTools.Squeezed","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Squeezed","text":"Squeezed(β,ϵ)\n\nCreate a Squeezed state with parameters β, ϵ.\n\n\n\n\n\n","category":"type"},{"location":"#Thermal-1","page":"PhaseSpaceTools.jl","title":"Thermal","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Thermal","category":"page"},{"location":"#PhaseSpaceTools.Thermal","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Thermal","text":"Thermal(β,n̄)\n\nCreate a Thermal state with parameters β, n̄.\n\n\n\n\n\n","category":"type"},{"location":"#Bogoliubov-1","page":"PhaseSpaceTools.jl","title":"Bogoliubov","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"Bogoliubov","category":"page"},{"location":"#PhaseSpaceTools.Bogoliubov","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.Bogoliubov","text":"Bogoliubov(u,v,n̄)\n\nCreate a Bogoliubov state with parameters u, v, n̄.\n\n\n\n\n\n","category":"type"},{"location":"#Phase-space-distributions-1","page":"PhaseSpaceTools.jl","title":"Phase-space distributions","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"The phase-space distributions supported for these states are","category":"page"},{"location":"#Husimi-Q-1","page":"PhaseSpaceTools.jl","title":"Husimi-Q","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"husimiQ","category":"page"},{"location":"#PhaseSpaceTools.husimiQ","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.husimiQ","text":"a,ā = husimiQ(state <: State,N)\n\nGenerate N samples from the Husimi-Q phase-space distribution for state.\n\nMoments of the Husimi-Q distribution generate quantum operator averages that are anti-normally ordered.\n\n\n\n\n\n","category":"function"},{"location":"#Glauber-P-2","page":"PhaseSpaceTools.jl","title":"Glauber-P","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"glauberP","category":"page"},{"location":"#PhaseSpaceTools.glauberP","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.glauberP","text":"a,ā = glauberP(state <: State,N)\n\nGenerate N samples from the Glauber-Sudarshan-P phase-space distribution for state.\n\nMoments of the Glauber-Sudarshan-P distribution generate quantum operator averages that are normally ordered.\n\n\n\n\n\n","category":"function"},{"location":"#Positive-P-2","page":"PhaseSpaceTools.jl","title":"Positive-P","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"positiveP","category":"page"},{"location":"#PhaseSpaceTools.positiveP","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.positiveP","text":"a,ā = positiveP(state <: State,N)\n\nGenerate N samples from the positive-P phase-space distribution for state.\n\nMoments of the positive-P distribution generate quantum operator averages that are normally ordered.\n\n\n\n\n\n","category":"function"},{"location":"#Wigner-1","page":"PhaseSpaceTools.jl","title":"Wigner","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"wigner","category":"page"},{"location":"#PhaseSpaceTools.wigner","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.wigner","text":"a,ā = wigner(state <: State,N)\n\nGenerate N samples from wigner phase-space distribution for state.\n\nMoments of the Wigner distribution generate symmetrically ordered quantum operator averages.\n\n\n\n\n\n","category":"function"},{"location":"#Positive-W-1","page":"PhaseSpaceTools.jl","title":"Positive-W","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"positiveW","category":"page"},{"location":"#PhaseSpaceTools.positiveW","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.positiveW","text":"a,ā = positiveW(state <: State,N)\n\nGenerate N samples from the positive-W phase-space distribution for state.\n\nImplemented states are\n\nFock(N)\n\nThe moments of the positive-W distribution generate quantum operator averages that are symmmetrically ordered.\n\n\n\n\n\n","category":"function"},{"location":"#Noises-1","page":"PhaseSpaceTools.jl","title":"Noises","text":"","category":"section"},{"location":"#randnc-1","page":"PhaseSpaceTools.jl","title":"randnc","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"randnc","category":"page"},{"location":"#PhaseSpaceTools.randnc","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.randnc","text":"a = randnc(N...)\n\nReturns an array of dimension length(N), containing samples of complex random variables with mean zero and variance one:\n\nlangle a^2rangle = 1 quadquadlangle a^2rangle = langle (a^*) ^2rangle = 0\n\nUseful for creating a range of noises that show up in phase-space simulations.\n\nUsage:\n\na = randnc(10) returns a 10-element Array{Complex{Float64},1}.\n\na = randnc(50,100) returns a 50x100-element Array{Complex{Float64},2}.\n\n\n\n\n\n","category":"function"},{"location":"#realnoise-1","page":"PhaseSpaceTools.jl","title":"realnoise","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"realnoise","category":"page"},{"location":"#PhaseSpaceTools.realnoise","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.realnoise","text":"Complex variable SDE's in DifferentialEquations.jl have complex noises by default dispatch on complex fields. This method defines a real noise for use with the library.\n\n\n\n\n\n","category":"function"},{"location":"#realbridge-1","page":"PhaseSpaceTools.jl","title":"realbridge","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"realbridge","category":"page"},{"location":"#PhaseSpaceTools.realbridge","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.realbridge","text":"Complex variable SDE's in DifferentialEquations.jl have complex noises by default dispatch on complex fields. This method defines a real brownian bridge for use with the library.\n\n\n\n\n\n","category":"function"},{"location":"#Solving-SDEs-1","page":"PhaseSpaceTools.jl","title":"Solving SDEs","text":"","category":"section"},{"location":"#Complex-noise-example-1","page":"PhaseSpaceTools.jl","title":"Complex noise example","text":"","category":"section"},{"location":"#Real-noise-example-1","page":"PhaseSpaceTools.jl","title":"Real noise example","text":"","category":"section"},{"location":"#Recovering-normal-order-1","page":"PhaseSpaceTools.jl","title":"Recovering normal order","text":"","category":"section"},{"location":"#Citing-1","page":"PhaseSpaceTools.jl","title":"Citing","text":"","category":"section"},{"location":"#Index-1","page":"PhaseSpaceTools.jl","title":"Index","text":"","category":"section"},{"location":"#","page":"PhaseSpaceTools.jl","title":"PhaseSpaceTools.jl","text":"","category":"page"}]
}
