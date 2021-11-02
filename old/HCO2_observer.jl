using Plots: print
using DifferentialEquations, Random, Plots, LinearAlgebra

# Internal dynamics perturbations (open file to change)
Random.seed!(1)
include("./HCO2_perturbations.jl")

# Model equations
# This HCO model was taken from the upcoming book by Drion, Franci and Sepulchre
include("./HCO2_kinetics.jl")
include("./HCO2_odes.jl")
include("./HCO2_sdes.jl")

# Definition of the timespan for numerical integration
Tfinal = 15000.
tspan = (0.0,Tfinal)

# Input currents: 
Iapp1 = t -> -1.0 
Iapp2 = t -> -1.0 

# Noise paramters
Noiseσ1 = Noiseσ2 =     0.25        # originally 0.2
TauNoise1 = TauNoise2 = 5.0         # originally 2.0

# Passive properties
gLeak1=gLeak2= 0.1   # originally 0.1
C1=C2=1.0

# Maximal conductances
gNa1=gNa2=300.;     # Sodium current maximal conductance
gKD1=gKD2=120.     # Delayed-rectifier potassium current maximal conductance
gA1=gA2=2.;         # A-type potassium current maximal conductance
gT1=gT2=5.;         # T-type calcium current maximal conductance
gH1=gH2=0.1;        # H-current maximal conductance
gL1=gL2=1.;         # L-type calcium current maximal conductance
gKCa1=gKCa2=2.;     # Calcium-activated potassium current maximal conductance; 

# Synaptic parameters
tausyn1=20.0; tausyn2=20.0
gsyn1=15; gsyn2=15
Esyn1=VK; Esyn2=VK

# Observer parameters
α₁ = 0.0005 # 0.0005
α₂ = 1.
γ = 10.
β = γ

## Numerical integration
# ODE parameters
p= ((Iapp1,Iapp2),(gT1,gT2),(gKD1,gKD2),(gH1,gH2),(gNa1,gNa2),(gA1,gA2),(gLeak1,gLeak2),(gL1,gL2),
    (gKCa1,gKCa2),(C1,C2),(tausyn1,tausyn2),(gsyn1,gsyn2),(Esyn1,Esyn2),(α₁,α₂,β,γ),
    "output",(TauNoise1,TauNoise2),(Noiseσ1,Noiseσ2))

# True system initial conditions
V01=-60.
V02=-0.
x01 = [V01 m_inf(V01,1) h_inf(V01,1) mH_inf(V01,1) mt_inf(V01,1) ht_inf(V01,1) mA_inf(V01,1) hA_inf(V01,1) mKD_inf(V01,1) mL_inf(V01,1) -10*gL1*mL_inf(V01,1)*(V01-VCa) msyn_inf(V01,1)]
x02 = [V02 m_inf(V02,2) h_inf(V02,2) mH_inf(V02,2) mt_inf(V02,2) ht_inf(V02,2) mA_inf(V02,2) hA_inf(V02,2) mKD_inf(V02,2) mL_inf(V02,2) -10*gL2*mL_inf(V02,2)*(V02-VCa) msyn_inf(V02,2)]
x0 = reshape([x01;x02],1,24)'
# Observer initial conditions
x̂01 = [-50. 0.5*ones(1,9) 0 0.5]
x̂02 = [-50. 0.5*ones(1,9) 0 0.5]
x̂0 = reshape([x̂01;x̂02],1,24)'
θ₀ = η.*[gNa1 gH1 gT1 gA1 gKD1 gL1 gKCa1 gLeak1 gsyn1;
         gNa2 gH2 gT2 gA2 gKD2 gL2 gKCa2 gLeak2 gsyn2]
P0 = (Matrix(I, 9, 9),Matrix(I, 9, 9))
ψ₀ = zeros(18,1)
z0 = [x0;x̂0;reshape(θ₀,18,1);reshape(P0[1],81,1);reshape(P0[2],81,1);ψ₀]
# Noise initial condition
e0 = [0;0]

# SDES
# prob = SDEProblem(HCO_sde,HCO_noise,[x0;e0],tspan,p)
prob = SDEProblem(HCO_observer_sde,HCO_observer_noise,[z0;e0],tspan,p)
sol = solve(prob,saveat=0.1)
# ODES
# prob = ODEProblem(HCO_ode,x0,tspan,p)
# prob = ODEProblem(HCO_observer_ode,z0,tspan,p)
# prob = ODEProblem(HCO_observer_besancon,[x0;x̂0;reshape(θ₀,16,1)],tspan,p)
# sol = solve(prob,Tsit5(),saveat=0.1) #,reltol=1e-6,abstol=1e-6
 
θ = sol[49:66,1,:]';
ψ = sol[193:208,1,:]';
# P1 = sol[65:128,1,:]';

# VOLTAGES
p1=plot(sol.t, sol[1,:],linewidth=1.5,legend=false)
p1=plot!(sol.t, sol[25,:],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
p2=plot(sol.t, sol[2,:],linewidth=1.5,legend=false)
p2=plot!(sol.t, sol[26,:],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
voltages = plot(p1,p2,layout=(2,1),legend=false)

# noise = plot(sol.t,sol[25,:],linewidth=1.5,legend=false)
# Pmatrix = plot(sol.t,P1,legend=false)

# Na and K CONDUCTANCES 
plt1 = plot(sol.t,θ[:,1:2],labels=["ĝNa₁" "ĝNa₂"],linecolor=["red" "blue"])
plt1 = plot!(sol.t,θ[:,9:10],labels=["ĝK₁" "ĝK₂"],linecolor=["green" "orange"]) 
plt1 = plot!([0,Tfinal],[gNa1,gNa1],linecolor=["black"],linestyle= :dash,labels="gNa")
plt1 = plot!([0,Tfinal],[gKD1,gKD1],linecolor=["black"],linestyle= :dash,labels="gK")
# # # xlims!(1000,Tfinal)

# # SYNAPTIC CONDUCTANCES
# plt1 = plot(sol.t,θ[:,17:18],labels=["ĝSyn₁" "ĝSyn₂"],linecolor=["red" "blue"]) 
# plt1 = plot!([0,Tfinal],[gsyn1,gsyn2],linecolor=["black"],linestyle= :dash,labels="gNa")

# plt2 = plot(sol.t,θ[:,5:6],labels=["ĝT₁" "ĝT₂"],linecolor=["red" "blue"]) 
# plt2 = plot!(sol.t,θ[:,7:8],labels=["ĝA₁" "ĝA₂"],linecolor=["green" "orange"]) 
# plt2 = plot!([0,Tfinal],[gT1,gT1],linecolor=["black"],linestyle= :dash,labels="gT")
# plt2 = plot!([0,Tfinal],[gA1,gA1],linecolor=["black"],linestyle= :dash,labels="gA")
# # xlims!(1000,Tfinal)

# plt3 = plot(sol.t,θ[:,13:14],labels=["ĝKCa₁" "ĝKCa₂"],linecolor=["red" "blue"])
# plt3 = plot!(sol.t,θ[:,11:12],labels=["ĝL₁" "ĝL₂"],linecolor=["green" "orange"])
# plt3 = plot!([0,Tfinal],[gKCa1,gKCa1],linecolor=["black"],linestyle= :dash,labels="gKCa")
# plt3 = plot!([0,Tfinal],[gL1,gL1],linecolor=["black"],linestyle= :dash,labels="gL")
# # xlims!(1000,Tfinal)

# plt4 = plot(sol.t,θ[:,3:4],labels=["ĝH₁" "ĝH₂"],linecolor=["red" "blue"])
# plt4 = plot!(sol.t,θ[:,15:16],labels=["ĝLeak₁" "ĝLeak₂"],linecolor=["green" "orange"])
# plt4 = plot!([0,Tfinal],[gH1,gH1],linecolor=["black"],linestyle= :dash,labels="gH")
# plt4 = plot!([0,Tfinal],[gLeak1,gLeak1],linecolor=["black"],linestyle= :dash,labels="gLeak")
# xlims!(1000,Tfinal)

# p3 = plot(sol.t,P1,legend=false)
#savefig("myplot.png")

## SIMULATE THE IDENTIFIED MODEL
Tfinal_id=5000.
tspan_id = (0.0,Tfinal_id)

θ₀ = reshape(θ₀,1,18)
p_id= ((Iapp1,Iapp2),θ[end,5:6],θ[end,9:10],θ[end,3:4],θ[end,1:2],θ[end,7:8],θ[end,15:16],θ[end,11:12],θ[end,13:14],(C1,C2),(tausyn1,tausyn2),(gsyn1,gsyn2),(Esyn1,Esyn2),(α₁,α₂,β,γ))
p₀= ((Iapp1,Iapp2),θ₀[5:6],θ₀[9:10],θ₀[3:4],θ₀[1:2],θ₀[7:8],θ₀[15:16],θ₀[11:12],θ₀[13:14],(C1,C2),(tausyn1,tausyn2),θ₀[17:18],(Esyn1,Esyn2),(α₁,α₂,β,γ)) 

prob_nom = ODEProblem(HCO_ode,x0,tspan_id,p)
prob_id = ODEProblem(HCO_ode,x0,tspan_id,p_id)

sol_nom = solve(prob_nom,Tsit5(),saveat=0.1) 
sol_id = solve(prob_id,Tsit5(),saveat=0.1)

p1id=plot(sol_nom.t, sol_nom[1,:],linewidth=1.5,legend=false)
p1id=plot!(sol_id.t, sol_id[1,:],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
p2id=plot(sol_nom.t, sol_nom[2,:],linewidth=1.5,legend=false)
p2id=plot!(sol_id.t, sol_id[2,:],linewidth=1.5,legend=false,linestyle= :dash,linecolor="black")
plot(p1id,p2id,layout=(2,1),legend=false)