[![DOI](https://zenodo.org/badge/863960464.svg)](https://doi.org/10.5281/zenodo.13912549)

# Simulation and analysis of Kuramoto-Sakaguchi populations interacting through complex networks


This repository contains the code accompaining the paper [From chimeras to extensive chaos in networks of heterogeneous
Kuramoto oscillator populations](https://arxiv.org/abs/2407.20408) by Pol Floriach, Jordi Garcia-Ojalvo, and Pau Clusella.

Two different set of codes are provided:

1. Some `auto-07p` files and instructions to obtain the main bifurcations displayied in the paper for a two-population model (Fig. 2). These is all provided in the `auto` folder. See the `TUTORIAL.md` in that folder.
2. The Julia code to reproduce most of the simulations in the paper. These functions are provided in the `KuramotoPopulationNetwork.jl` module in the `src` directory. In the following we show how to use the functions in that module to simulate the system using some examples.

	 
## Using the `KuramotoPopulationNetwork` module 

Open a Julia REPL in the directory, activate the project and load/install all required dependencies:

```julia
] activate .
] instantiate
using KuramotoPopulationNetwork
```

Apart from these packages, to run these examples we also need three other packages.
Load (or install) them with:

```
using Graphs, Plots, ColorSchemes
```

Then we can use `Graphs.jl` to generate a ring network with nearest neighbours:

```julia
N = 128;
g_ring = watts_strogatz(N, 2, 0.0);
```

Now we can reproduce of Fig. 5(a) in the paper:

```julia
t,R,ϕ,u0 = simulation(;K=7.0,p=0.9,α=1.2,gr=g_ring,trans=5e2,tmax=4e2,ic="homogeneous");
n = size(R,1)
range = n-10000:1:n;
heatmap(1:128,t[range],ϕ[range,:],c=:cyclic_mrybm_35_75_c68_n256)
heatmap(1:128,t[range],R[range,:],c=:linear_bmy_10_95_c78_n256)
```

We can compute and plot the full Lyapunov spectra for this case (this can take some minutes):

```julia
λ, = simulationLE(;K=7,p=0.9,α=1.2,gr=g_ring,trans=2e2,tmax=1e3,nLE=2*size(g_ring,1),ic=u0);
plot(λ)
```
With this data we can compute the attractor dimension and the dynamical entropy:

```julia
dimKY(λ)
entropy(λ)
```

Reproduction of Fig. 5(b) in the paper:
```julia
@time t,R,ϕ,u0 = simulation(;K=15.0,p=0.5,α=1.2,gr=g_ring,trans=5e2,tmax=4e2,ic="homogeneous");
n = size(R,1)
range = n-10000:1:n;
heatmap(1:128,t[range],ϕ[range,:],c=:cyclic_mrybm_35_75_c68_n256)
heatmap(1:128,t[range],R[range,:],c=:linear_bmy_10_95_c78_n256)
``` 

Simulation on the slow regime that appears for p=0.5. Notice we change `tmax` and `dt`:
```julia
t,R,ϕ,u0 = simulation(;K=20.0,p=0.5,α=1.2,gr=graph,trans=0e3,tmax=40e3,dt=1e-1,ic="homogeneous");
n = size(R,1)
range = n-1000:10:n;
heatmap(1:128,t[range],ϕ[range,:],c=:cyclic_mrybm_35_75_c68_n256)
heatmap(1:128,t[range],R[range,:],c=:linear_bmy_10_95_c78_n256)
```

Let's compute now the LE for a simulation of the ER network:
```julia
n = 128;
avk = 10;
g_er = erdos_renyi(n, avk*1.0/n);

t,R,ϕ,u0 = simulation(;K=14.0,p=0.9,α=1.2,gr=g_er,trans=5e2,tmax=4e2,ic="homogeneous");
λ, = simulationLE(;K=14,p=0.9,α=1.2,gr=g_er,trans=2e2,tmax=1e3,nLE=2*size(g_er,1),ic=u0);
```
and plot the results:

```julia
n = size(R,1)
trange = n-10000:1:n;
heatmap(1:128,t[trange],ϕ[trange,:],c=:cyclic_mrybm_35_75_c68_n256)
heatmap(1:128,t[trange],R[trange,:],c=:linear_bmy_10_95_c78_n256)
plot(λ)
```
