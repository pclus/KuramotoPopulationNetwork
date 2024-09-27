module KuramotoPopulationNetwork

using DifferentialEquations, Kronecker
using Graphs, LinearAlgebra, DynamicalSystems
using DelimitedFiles,SparseArrays,Statistics

export OAnet!, Jacobian_OAnet!, simulation, simulationLE, statsR, dimKY, entropy

# ------------------------------------------------
"""
        OAnet!(dv,v,par,t)

Coupled Ott-Antonsen equations for integration. 

        * `dv` and `v` are the velocity and state variables respectively.
        * `par` is a tuple that contains several parameters and auxiliary variables in order to improve performance:
                `K` is the global coupling strength of the system.
                `p` is the self-strength (μ in the paper).
                `Γ` is the heterogeneity of the system (Δ in the paper).
                `ω` is the mean natural frequency of the oscillators.
                `cosα` and `sinα`are the cosine and sine of the Sakaguchi phase shift `α`
                `list` is the adjacency list of the graph (see the documentation of Graphs.jl).
                `ideg` is a vector `ideg[i]=1/k[i]` where `k[i]` is the in-degree of node `i`.
                `sx,sy,A,B,C,D,P,Q` are N-dimensional auxiliary vectors used for performance.
                `N` is the number of nodes in the graph.
        * `t` is the time variable.

"""
function OAnet!(dv, v, par, t)
    K,p,Γ,ω,cosα,sinα, list,ideg,sx,sy,A,B,C,D,P,Q,N = par
    @inbounds x = @view v[1:2:end] 
    @inbounds y = @view v[2:2:end]
    @inbounds dx = @view dv[1:2:end]
    @inbounds dy = @view dv[2:2:end]
    @fastmath @. dx = -Γ*x - ω*y
    @fastmath @. dy =  ω*x - Γ*y

    @fastmath @. P =  (x*x) - (y*y)
    @fastmath @. Q =  2*x*y
    @fastmath @. A =  cosα*( 1-P ) + Q*sinα 
    @fastmath @. B =  sinα*( 1-P ) - Q*cosα 
    @fastmath @. C =  sinα*(-1-P ) - Q*cosα 
    @fastmath @. D =  cosα*( 1+P ) - Q*sinα 

    @fastmath @inbounds sx[:] = 0.5*K*[(1-p)*ideg[i]*sum(x[list[i]])+p*x[i] for i in 1:N]
    @fastmath @inbounds sy[:] = 0.5*K*[(1-p)*ideg[i]*sum(y[list[i]])+p*y[i] for i in 1:N]
    @fastmath @. dx += A*sx + B*sy
    @fastmath @. dy += C*sx + D*sy
end
# ------------------------------------------------

# ------------------------------------------------
"""
        Jacobian_OAnet!(M,v,par,t)

Analytical Jacobian of the system for the computation of Lyapunov exponents.
        * `M` is a matrix where the Jacobian is stored.
        * `v` is the vector of state variables.
        * `par` are the parameters of the system, as in `OAnet!`. 
        * `t` is the time variable.
"""
function Jacobian_OAnet!(M,v,par,t)
    K,p,Γ,ω,cosα,sinα, list,ideg,sx,sy,A,B,C,D,P,Q,N = par

    @inbounds x = @view v[1:2:end] 
    @inbounds y = @view v[2:2:end]

    d0 = diagind(M);
    d1 = diagind(M,1);
    ds1 = diagind(M,-1);
    @inbounds a = view(M,d0[1:2:end])
    @inbounds b = view(M,d1[1:2:end])
    @inbounds c = view(M,ds1[1:2:end])
    @inbounds d = view(M,d0[2:2:end])

    # block diagonal, for optimization xcos, xsin, ycos and ysin could be precomputed
    pre = 0.5*K*p
    @fastmath @. a = -Γ + 2.0*((-x*cosα+y*sinα)*sx + (-x*sinα-y*cosα)*sy) + A*pre
    @fastmath @. b = -ω + 2.0*(( y*cosα+x*sinα)*sx + ( y*sinα-x*cosα)*sy) + B*pre
    @fastmath @. c =  ω + 2.0*((-x*sinα-y*cosα)*sx + ( x*cosα-y*sinα)*sy) + C*pre
    @fastmath @. d = -Γ + 2.0*(( y*sinα-x*cosα)*sx + (-y*cosα-x*sinα)*sy) + D*pre

    for j in 1:N 
         for i in list[j]
           @inbounds pre = 0.5*K*ideg[i]*(1-p)
           @inbounds M[2*i-1, 2*j-1] = pre*A[i]
           @inbounds M[2*i-1, 2*j  ] = pre*B[i]
           @inbounds M[2*i  , 2*j-1] = pre*C[i]
           @inbounds M[2*i  , 2*j  ] = pre*D[i]
        end
    end
end
# ------------------------------------------------

# ------------------------------------------------
"""
        simulation(;kwargs...)

Performs a simulations of the system with the given parameters and network topology.
The (optional) function arguments are:
        * `K` the coupling of the system.
        * `p` the self-coupling strength (μ in the paper).
        * 'α' the Sakaguchi phase shift.
        * 'gr' a Graph object containing the network topology, see the documentation of Graphs.jl. The row-normalization of the inputs in the connectivity matrix is computed authomatically, thus the connectivity of `gr` does not need to be row-normalized.
        * 'tmax' simulation time after transient.
        * 'trans' transient simulation time.
        * 'ic' initial condition of the system. Can be:
                - A vector containing the state variables to initialize the system.
                - The String "random" for totally random initial values.
                - The String "antiphase" for simulations starting close to the antiphase state.
                - Any other string, which will cause the system to initialize close to the homogeneous state.
        * `dt` integration time step.


Returns:
        * A time vector.
        * A matrix `R` with the values of the Kuramoto order parameters at all times and nodes.
        * A matrix `Φ` with the values of the collective phase at all times and nodes.
        * The final state of the system `u0`.
"""
function simulation(;K=20.0,p=0.9,α=1.45,gr=erdos_renyi(20, 0.2),tmax=200.0,trans=0.0,ic="random",dt=0.01)
    N = length(gr.fadjlist)
    Γ = 1.0
    ω = tan(α)*(K*cos(α)-Γ)
    deg = degree(gr);
    ideg =  1.0./deg;
    par = K,p,Γ,ω,cos(α),sin(α),gr.fadjlist,ideg,zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),N;
    
    if typeof(ic)==String || length(ic)!=2*N
        if ic=="random"
            r0 = rand(N); # random r's
            ϕ0 = rand(N)*2*π;
        else
            d = 1-2*Γ/(K*cos(α))
            r0 = ( d<0 ? 0.0 : sqrt(d))
            r0 = @. r0*ones(N)+1e-3*randn(N);
            r0 = @. abs(r0);
            ϕ0 =1e-3*randn(N) # random perturbation
            if ic=="antiphase"
                @. ϕ0[1:2:N]=ϕ0[1:2:N]+π # random perturbation
            end
        end
        x0 = @. r0*cos(ϕ0);
        y0 = @. r0*sin(ϕ0);

        u0 = [ (i%2==1 ? x0[Int(ceil(i/2))] : y0[Int(i/2)] ) for i in 1:2*N];
    else
        u0 = ic 
    end


    ds=ContinuousDynamicalSystem(OAnet!, u0, par; 
        diffeq=(alg=RK4(),dt=dt,adaptive=false)) 
  
    ddt=(dt < 0.01 ? 0.01 : dt); # sample step (equivalent to "saveat")
    tr=trajectory(ds, tmax ;Δt=ddt,Ttr=trans) 
    R = [ @. sqrt(tr[1][:,2*i-1]^2+tr[1][:,2*i]^2) for i in 1:N ];
    ϕ = [ @. atan(tr[1][:,2*i],tr[1][:,2*i-1]) for i in 1:N ];

    R = reduce(hcat,R);
    ϕ = reduce(hcat,ϕ);
    return collect(tr[2]),R,ϕ,u0
end
# ------------------------------------------------

# ------------------------------------------------
"""
        simulationLE(;kwargs...)

Performs a simulations of the system in state and tangent space and provides the Lyapunov exponents.
All arguments are optional, and contains the same arguments as in `simulation` and two additional ones:
        * `nLE` number of exponents to be computed. 
        * `ns` time window between successive calls of the QR-orthogonalization for computation of the exponents.

Returns: 
        * `λ` the computed exponents.
        * `ds` the DynamicalSystem object at the end of the simulation.
"""
function simulationLE(;K=20.0,p=0.9,α=1.45,gr=erdos_renyi(20, 0.2),tmax=2e2,trans=2e2,nLE=1,ic="random",dt=0.01,ns=1)
    N = length(gr.fadjlist)
    Γ = 1.0
    ω = tan(α)*(K*cos(α)-Γ)
    deg = degree(gr);
    ideg =  1.0./deg;
    par = K,p,Γ,ω,cos(α),sin(α),gr.fadjlist,ideg,zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),N;
    
    if typeof(ic)==String || length(ic)!=2*N
        if ic=="random"
            r0 = rand(N); # random r's
            ϕ0 = rand(N)*2*π;
        
        else
            d = 1-2*Γ/(K*cos(α))
            r0 = ( d<0 ? 0.0 : sqrt(d))
            r0 = @. r0*ones(N)+1e-3*randn(N);
            r0 = @. abs(r0)
            ϕ0 = 1e-3*randn(N)
            if ic=="antiphase"
                @. ϕ0[1:2:N]=ϕ0[1:2:N]+π
            end
        end
        x0 = @. r0*cos(ϕ0);
        y0 = @. r0*sin(ϕ0);
        u0 = [ (i%2==1 ? x0[Int(ceil(i/2))] : y0[Int(i/2)] ) for i in 1:2*N];
    else
        u0 = ic 
    end

    # --------------------------------
    # Sparse Jacobian: at least for a ring of N=128 nodes the performance is not noticiably changed,
    # but for N = 256 it starts to improve. 
    if N>=256 
        Adj = adjacency_matrix(gr)+I(N)
        One = ones(2,2)*1.0
        J0 = Adj⊗One
        J0 = sparse(J0)
    else
        J0=zeros(2*N,2*N)
    end
    # --------------------------------
    
    ds=ContinuousDynamicalSystem(OAnet!, u0, par; 
        diffeq=(alg=RK4(),dt=dt,adaptive=false))
    tands = TangentDynamicalSystem(ds; k=nLE,J = Jacobian_OAnet!,J0=J0)
    
    Tdisp=tmax # computation time AFTER transient
    trans=trans # transient

    λ = lyapunovspectrum(tands, Int(1*Tdisp/ns); Δt = ns, Ttr=trans)
    return λ,ds
end
# ------------------------------------------------

# ------------------------------------------------
"""
        statsR(t,R)

Computes some statistical quantities from the matrix `R` returned by `simulation`.
"""
function statsR(t,R)
    n = size(R,2)
    nt = size(R,1)
    Rs = [ mean(R[:,i]) for i in 1:n ]      # mean over time for each population
    std_o_time = [ std(R[:,i]) for i in 1:n ]    # deviations over time (oscillations or not)
    std_o_space = [ std(R[j,:]) for j in 1:nt ]      # deviations over space (homogeneous or not)
    return mean(Rs),mean(std_o_time),mean(std_o_space)
end
# ------------------------------------------------

# ------------------------------------------------
# Kaplan-York dimension
"""
        dimKY(λ;tol=1e-4)

Computes the Kaplan-York dimension on the exponents in `λ`. 
`tol` specifies the threshold at which a exponent with |λ|<tol 
is considered to be zero.
"""
function dimKY(λ;tol=1e-4)
    n = length(λ)
    μ = [(abs(l) < tol ? 0.0 : l) for l in λ ]
    Sλ = [ sum(μ[1:i]) for i in 1:n]
    j = sum(Sλ .>= 0.0)
    if j>0 && j<n
        dim = j + Sλ[j]/abs(μ[j+1])
    else
        dim = 0
    end
    return dim
end
# ------------------------------------------------

# ------------------------------------------------
# Kolmogorov-Sinai entropy
"""
        entropy(λ)

Computes the Kolmogorov-Sinai entropy on the exponents given in `λ`
"""
function entropy(λ)
    j = sum(λ .> 0.0)
    return (j==0 ? 0 : sum(λ[1:j]))
end
# ------------------------------------------------

end #module
