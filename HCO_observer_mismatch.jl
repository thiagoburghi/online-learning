using Plots: print
using DifferentialEquations, Random, Plots, LinearAlgebra, DelimitedFiles

# Flag for saving data to .txt files 
save_data = 0

# Select if true model (true) or adaptive observer (false) dynamics will be perturbed
perturbed = true

# Internal dynamics perturbation
# include("./HCO_perturbations.jl")
Random.seed!(2)
δₐ = 0.02                   # Perturbation amplitude
δᵦ = 0.02                   # Perturbation baseline
δ = δₐ*2*(0.5.-rand(Float64, (46,2)))
δ = δ .+ δᵦ*sign.(δ)        # perturbation = amplitude + sign(amplitude)*baseline
δ = (1 .+ δ)
if perturbed
    # true model is perturbed
    δ = hcat(δ, ones(46,2))
else
    # adaptive observer is perturbed
    δ = hcat(ones(46,2), δ)
end

# Models (include after defining the perturbation δ)
include("./multisines.jl")
include("./HCO_kinetics.jl")
include("./HCO_odes.jl")

# Definition of the timespan for numerical integration
Tfinal = 20000.
tspan = (0.0,Tfinal)

# Input currents:
Iapp1 = t -> -0.65 #.+ multisines(t;Tmax=10000,A=0.1,N=100)
Iapp2 = t -> -0.65 #.+ multisines(t;Tmax=10000,A=0.1,N=100)
#plot(0:0.01:10000,Iapp1(0:0.01:10000))

# Nominal parameters
C1 = C2 = 1.0
if perturbed
    (gNa1,gNa2) = (47.3, 57.4)      #(53.5,58.6) 1-2%   # Sodium current maximal conductance
    (gK1,gK2) = (33.1, 45.3)        #(36.5,42.6)     # Potassium current maximal conductance
    (gCa1,gCa2) = (t->0.08,t->0.18)      #(0.1,0.16)    # Calcium current maximal conductance; 
    (gL1,gL2) = (0.04, 0.03)        #(0.04,0.03)      # Leak maximal conductance
    (gSyn1,gSyn2) = (6.1, 2.4)      #(5.5,4.2)   # Synaptic maximal conductance
else
    (gNa1,gNa2) = (60.,60.)     # Sodium current maximal conductance
    (gK1,gK2) = (40.,40.)       # Potassium current maximal conductance
    (gCa1,gCa2) = (t->0.14,t->0.14)   # Calcium current maximal conductance; 
    (gL1,gL2) = (0.035,0.035)   # Leak maximal conductance
    (gSyn1,gSyn2) = (4,4)       # Synaptic maximal conductance
end

# Synaptic parameters
kfA1 = kfA2 = 2;        
krA1 = krA2 = 0.1;

# Observer parameters
α₁ = 0.0001 # 0.0001
α₂ = 1.
γ = 5
β = γ

# True system initial conditions
V01= -64
V02= -66.
x01 = [V01 0.5*ones(1,6)]   #gamNa_m(V01,1) gamNa_h(V01,1) gamK_m(V01,1) gamCa_m(V01,1) gamCa_h(V01,1) gamSyn(V01,1)]
x02 = [V02 0.5*ones(1,6)]   #gamNa_m(V02,2) gamNa_h(V02,2) gamK_m(V02,2) gamCa_m(V02,2) gamCa_h(V02,2) gamSyn(V02,2)]
x0 = reshape([x01;x02],1,14)'
# Observer initial conditions
x̂01 = [-50. zeros(1,6)]     #0.5*ones(1,6)]
x̂02 = [-50. zeros(1,6)]     #0.5*ones(1,6)]
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

## ESTIMATE BOTH TRUE SYSTEM AND OBSERVER
p = ((Iapp1,Iapp2),(C1,C2),(gNa1,gNa2),(gK1,gK2),(gCa1,gCa2),(gL1,gL2),(gSyn1,gSyn2),
    (kfA1,kfA2),(krA1,krA2),(α₁,α₂,β,γ))

prob = ODEProblem(HCO_observer_ode,z0,tspan,p)
sol = solve(prob,Tsit5(),saveat=0.1,maxiters=1e6) #,reltol=1e-6,abstol=1e-6,AutoTsit5(Rosenbrock23()),Tsit5()
v = sol[1:2,1,:]';
v̂ = sol[15:16,1,:]';
θ = sol[29:38,1,:]';

# VOLTAGES
p1=plot(sol.t, v[:,1],linewidth=1.5,legend=false)
p2=plot(sol.t, v[:,2],linewidth=1.5,legend=false)
voltages = plot(p1,p2,layout=(2,1),legend=false)

# Na CONDUCTANCES 
plt1 = plot(sol.t,θ[:,1:2],labels=["ĝNa₁" "ĝNa₂"],linecolor=["red" "blue"])
plt1 = plot!([0,Tfinal],[gNa1,gNa1],linecolor=["red"],linestyle= :dash,labels="gNa₁")
plt1 = plot!([0,Tfinal],[gNa2,gNa2],linecolor=["blue"],linestyle= :dash,labels="gNa₂")

# K CONDUCTANCES 
plt2 = plot(sol.t,θ[:,3:4],labels=["ĝK₁" "ĝK₂"],linecolor=["red" "blue"]) 
plt2 = plot!([0,Tfinal],[gK1,gK1],linecolor=["red"],linestyle= :dash,labels="gK₁")
plt2 = plot!([0,Tfinal],[gK2,gK2],linecolor=["blue"],linestyle= :dash,labels="gK₂")

# Ca CONDUCTANCES 
plt3 = plot(sol.t,θ[:,5:6],labels=["ĝCa₁" "ĝCa₂"],linecolor=["red" "blue"])
plt3 = plot!([0,Tfinal],[gCa1,gCa1],linecolor=["red"],linestyle= :dash,labels="gCa₁")
plt3 = plot!([0,Tfinal],[gCa2,gCa2],linecolor=["blue"],linestyle= :dash,labels="gCa₂")

# Leak CONDUCTANCES
plt4 = plot(sol.t,θ[:,7:8],labels=["ĝL₁" "ĝL₂"],linecolor=["red" "blue"]) 
plt4 = plot!([0,Tfinal],[gL1,gL1],linecolor=["red"],linestyle= :dash,labels="gL₁")
plt4 = plot!([0,Tfinal],[gL2,gL2],linecolor=["blue"],linestyle= :dash,labels="gL₂")

# # SYNAPTIC CONDUCTANCES
plt5 = plot(sol.t,θ[:,9:10],labels=["ĝSyn₁" "ĝSyn₂"],linecolor=["red" "blue"])
plt5 = plot!([0,Tfinal],[gSyn1,gSyn1],linecolor=["red"],linestyle= :dash,labels="gSyn₁")
plt5 = plot!([0,Tfinal],[gSyn2,gSyn2],linecolor=["blue"],linestyle= :dash,labels="gSyn₂")

if save_data == 1
    writedlm("./../data/HCO_parameters_mismatch.txt",  hcat(sol.t,θ), " ")
end
plot(plt1,plt2,plt5,layout=(3,1),legend=true)
# plot(plt3,plt4,layout=(2,1),legend=true)

## SIMULATE THE IDENTIFIED MODEL
Tfinal_id=10000.
tspan_id = (0.0,Tfinal_id)

Iapp_id = t -> -0.65
p_id = ((Iapp_id,Iapp_id),(C1,C2),θ[end,1:2],θ[end,3:4],(t->θ[end,5],t->θ[end,6]),θ[end,7:8],θ[end,9:10],(kfA1,kfA2),(krA1,krA2),"observer")
p_nom = ((Iapp_id,Iapp_id),(C1,C2),p[3],p[4],p[5],p[6],p[7],p[8],p[9],"nominal")

prob_id = ODEProblem(HCO_ode,x0,tspan_id,p_id)
prob_nom = ODEProblem(HCO_ode,x0,tspan_id,p_nom)

sol_nom = solve(prob_nom,Tsit5(),saveat=0.1)
sol_id = solve(prob_id,Tsit5(),saveat=0.1)
v_nom = sol_nom[1:2,1,:]';
v_id = sol_id[1:2,1,:]';

p1id=plot(sol_nom.t, v_nom[:,1],linewidth=1.5,legend=false)
p3id=plot(sol_id.t, v_id[:,1],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
p2id=plot(sol_nom.t, v_nom[:,2],linewidth=1.5,legend=false)
p4id=plot(sol_id.t, v_id[:,2],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
plot(p1id,p2id,p3id,p4id,layout=(4,1),legend=false)

if save_data == 1
    writedlm("./../data/HCO_voltages_mismatch.txt",  hcat(sol_nom.t,v_nom,v_id), " ")
end

# ## ESTIMATE TRUE SYSTEM ONLY
# # ODE parameters
# p = ((Iapp1,Iapp2),(C1,C2),(gNa1,gNa2),(gK1,gK2),(gCa1,gCa2),(gL1,gL2),(gSyn1,gSyn2),
#     (kfA1,kfA2),(krA1,krA2),"nominal")

# prob = ODEProblem(HCO_ode,x0,tspan,p)
# sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6,saveat=0.1,maxiters=1e6)

# # VOLTAGES
# p1=plot(sol.t, sol[1,:],linewidth=1.5,legend=false)
# p2=plot(sol.t, sol[2,:],linewidth=1.5,legend=false)
# voltages = plot(p1,p2,layout=(2,1),legend=false)