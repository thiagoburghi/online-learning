using Plots: print
using DifferentialEquations, Random, Plots, LinearAlgebra, DelimitedFiles

# Models (include after defining the perturbation δ)
include("./HCO_kinetics.jl")
include("./HCO_odes.jl")

# General parameters
save_data = false
outinj = false
σnoise = 2.0
pcerror = 0.01   # percent error in true internal dynamics
pcêrror = 0.0   # percent error in observer internal dynamics

# Time parameters
Tfinal = 10000
tspan = (0.0,Tfinal)
dt_sim = 0.001
dt_plot = 0.1

# Observer parameters
α = 0.0025#0.0025
γ = 0.1#0.5
β = 0.0#0.0

# True system parameters
C1 = C2 = 1.0
(gNa1,gNa2) = (60.,60.)                 # Sodium current maximal conductance
(gK1,gK2) = (40.,40.)                   # Potassium current maximal conductance
(gL1,gL2) = (0.035,0.035)               # Leak maximal conductance
(gSyn1,gSyn2) = (4,4)  #(4,4)           # Synaptic maximal conductance 4
gCa0 = 0.11#0.1                          # Time-varying calcium
gCaf = 0.17
sig = x -> 1 ./ (1+exp(-x));
gCa = t -> gCa0 .+ (gCaf-gCa0)*sig((t-Tfinal/2)/1250)
(gCa1,gCa2) = (t -> gCa.(t), t -> gCa.(t))   
a1 = a2 = 2       
b1 = b2 = 0.1

# True system initial conditions 
# (chosen from steady state oscillations with gCa = 0.1)
V0 = v_ss
w₀ = w_ss 
x0 = [V0;w₀]

# Input currents:
Iapp1 = t -> -0.65 
Iapp2 = t -> -0.65

# Observer initial conditions
x̂01 = [-50. zeros(1,6)]
x̂02 = [-50. zeros(1,6)]
x̂0 = reshape([x̂01;x̂02],1,14)'
gNa1₀ = gNa2₀ = 80
gK1₀ = gK2₀ = 80
gCa1₀ = gCa2₀ = 1 
gL1₀ = gL2₀ = 1
gSyn1₀ = gSyn2₀ = 10
θ₀ = reshape([gNa1₀ gK1₀ gCa1₀ gL1₀ gSyn1₀; gNa2₀ gK2₀ gCa2₀ gL2₀ gSyn2₀],10,1)
P0 = (Matrix(I, 5, 5),Matrix(I, 5, 5))
ψ₀ = zeros(10,1)
z0 = [x0;x̂0;θ₀;reshape(P0[1],25,1);reshape(P0[2],25,1);ψ₀]

seed=101
# for seed = 11:19
# Flag for saving data to .txt files 
rng = MersenneTwister(seed)

# Internal parameters
η = (internal_parameters(pcerror,rng),internal_parameters(pcerror,rng))
η̂ = (internal_parameters(pcêrror,rng),internal_parameters(pcêrror,rng))

# ESTIMATE BOTH TRUE SYSTEM AND OBSERVER
p = ((Iapp1,Iapp2),(C1,C2),(gNa1,gNa2),(gK1,gK2),(gCa1,gCa2),(gL1,gL2),(gSyn1,gSyn2),
    (a1,a2),(b1,b2),(α,β,γ),(η,η̂),outinj,(σnoise,dt_sim))

##### ODE #####
# prob = ODEProblem(HCO_observer_ode,z0,tspan,p)
# sol = solve(prob,Rodas4(),saveat=0.1,reltol=1e-6,abstol=1e-6) #,reltol=1e-6,abstol=1e-6,AutoTsit5(Rosenbrock23()),Tsit5()

##### SDE #####
prob = SDEProblem(f_SDE!,g_SDE!,[z0;0;0],tspan,p,seed=123)
sol = solve(prob,EM(),dt=dt_sim,saveat=dt_plot,maxiters=1e7)
noise = sol[99:100,1,:]';

v = sol[1:2,1,:]';
w = sol[3:14,1,:]';
v₀ = v[end,:] 
w₀ = w[end,:] 
v̂ = sol[15:16,1,:]';
θ = sol[29:38,1,:]';

# NO OBSERVER
# prob = ODEProblem(HCO_ode,z0[1:14],tspan,p)
# sol = solve(prob,Rodas4(),saveat=0.1,reltol=1e-6,abstol=1e-6) #,reltol=1e-6,abstol=1e-6,AutoTsit5(Rosenbrock23()),Tsit5()
# v = sol[1:2,:]';
# w = sol[3:14,:]';

# prob = SDEProblem(f_SDE_no_observer!,g_SDE_no_observer!,z0[1:14],tspan,p)
# sol = solve(prob,EM(),dt=dt_sim,saveat=dt_plot)
# v = sol[1:2,:]';

if save_data
    writedlm(string("./../data/HCO_truevoltages_seed=",seed,"pcerror=",pcerror,".txt"),  hcat(sol.t[1:1:end],v[1:1:end,:]), " ")
    writedlm(string("./../data/HCO_parameters_seed=",seed,"pcerror=",pcerror,".txt"),  hcat(sol.t[1:1:end],θ[1:1:end,:],gCa.(sol.t)[1:1:end,:]), " ")
end

# VOLTAGES
p1=plot(sol.t, v[:,1],linewidth=1.5,legend=false)
p2=plot(sol.t, v[:,2],linewidth=1.5,legend=false)
voltages = plot(p1,p2,layout=(2,1),legend=false)
# png(string("./figs/voltages_seed=",seed))
# plot(sol.t,noise[:,1])

## PLOTS

# Na CONDUCTANCES 
plt1 = plot(sol.t,θ[:,1:2],labels=["ĝNa₁" "ĝNa₂"],linecolor=["red" "blue"])
plt1 = plot!([0,Tfinal],[gNa1,gNa1],linecolor="red",linestyle= :dash,labels="gNa₁")
plt1 = plot!([0,Tfinal],[gNa2,gNa2],linecolor="blue",linestyle= :dash,labels="gNa₂")

# K CONDUCTANCES 
plt2 = plot(sol.t,θ[:,3:4],labels=["ĝK₁" "ĝK₂"],linecolor=["red" "blue"]) 
plt2 = plot!([0,Tfinal],[gK1,gK1],linecolor="red",linestyle= :dash,labels="gK₁")
plt2 = plot!([0,Tfinal],[gK2,gK2],linecolor="blue",linestyle= :dash,labels="gK₂")

# Ca CONDUCTANCES 
plt3 = plot(sol.t,θ[:,5:6],labels=["ĝCa₁" "ĝCa₂"],linecolor=["red" "blue"])
plt3 = plot!(sol.t,[gCa1(sol.t)],linecolor="red",linestyle= :dash,labels="gCa₁")
plt3 = plot!(sol.t,[gCa2(sol.t)],linecolor="blue",linestyle= :dash,labels="gCa₂")
ylims!((-0.1,0.2))

# Leak CONDUCTANCES
plt4 = plot(sol.t,θ[:,7:8],labels=["ĝL₁" "ĝL₂"],linecolor=["red" "blue"]) 
plt4 = plot!([0,Tfinal],[gL1,gL1],linecolor="red",linestyle= :dash,labels="gL₁")
plt4 = plot!([0,Tfinal],[gL2,gL2],linecolor="blue",linestyle= :dash,labels="gL₂")

# # SYNAPTIC CONDUCTANCES
plt5 = plot(sol.t,θ[:,9:10],labels=["ĝSyn₁" "ĝSyn₂"],linecolor=["red" "blue"])
plt5 = plot!([0,Tfinal],[gSyn1,gSyn1],linecolor="red",linestyle= :dash,labels="gSyn₁")
plt5 = plot!([0,Tfinal],[gSyn2,gSyn2],linecolor="blue",linestyle= :dash,labels="gSyn₂")

firstpars = plot(p1,plt1,plt2,layout=(3,1),legend=true)
# png(string("./figs/firstpars_seed=",seed))
lastpars = plot(plt3,plt4,plt5,layout=(3,1),legend=true)
# png(string("./figs/lastpars_seed=",seed))
# end

## Validation
(ĝNa1,ĝNa2) = θ[end,1:2]                 # Sodium current maximal conductance
(ĝK1,ĝK2) = θ[end,3:4]                   # Potassium current maximal conductance
(ĝL1,ĝL2) = θ[end,7:8]               # Leak maximal conductance
(ĝSyn1,ĝSyn2) = θ[end,9:10] #(4,4)           # Synaptic maximal conductance 4

# ĝCa0 = 0.12
# ĝCa = t -> ĝCa0
(ĝCa1,ĝCa2) = (t -> θ[end,5], t -> θ[end,6])  

p̂ = ((Iapp1,Iapp2),(C1,C2),(ĝNa1,ĝNa2),(ĝK1,ĝK2),(ĝCa1,ĝCa2),(ĝL1,ĝL2),(ĝSyn1,ĝSyn2),
    (a1,a2),(b1,b2),(α,β,γ),(η̂,η̂),outinj,(σnoise,dt_sim))

prob = ODEProblem(HCO_ode,z0[1:14],tspan,p̂)
s̄ol = solve(prob,Rodas4(),saveat=0.1,reltol=1e-6,abstol=1e-6) #,reltol=1e-6,abstol=1e-6,AutoTsit5(Rosenbrock23()),Tsit5()
v̄ = s̄ol[1:2,:]';

p1=plot(s̄ol.t, v̄[:,1],linewidth=1.5,legend=false)
p2=plot(s̄ol.t, v̄[:,2],linewidth=1.5,legend=false)
voltages = plot(p1,p2,layout=(2,1),legend=false)