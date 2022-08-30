using Plots: print
using DifferentialEquations, Random, Plots, LinearAlgebra, DelimitedFiles
import Statistics: mean

# Flag for saving data to .txt files 
save_data = false
outinj = true

include("HH_odes.jl")

# True Parameters 
c = 1.
μ = [120.,36.,0.3]
θ = [μ[1]/c,μ[2]/c,μ[3]/c,1/c]
η = [-40,-62,-53]
E = (55.,-77.,-54.4)
Iapp = t -> 10 + sin(2*pi/10*t)
k = (9.0,-7.0,15.0)

# Observer parameters
α,γ,β = 0.1, 1.0, 1.0
#α,γ,β = 0.0001, 0.1, 1.0

# Integration parameters
dt_plot = 0.1 
dt_sim = 0.001
#Tfinal = 300.
Tfinal = 300.
tspan = (0.,Tfinal)

# Noise parameters
σnoise = 0 #0.05

# Initial conditions
x₀ = [-30 0.5 0.5 0.5];
x̂₀ = [-30 0.0 0.0 0.0];
θ̂₀ = [78 78 10 2];
η̂₀ = [-20 -20 -20];
P₀ = Matrix(I, 7, 7);
ψ₁₀ = [0 0 0 0 0 0 0];
ψ₂₀ = [0 0 0];
z₀ = [x₀ x̂₀ θ̂₀ η̂₀ reshape(P₀,1,49) ψ₁₀ ψ₂₀]

# System parameters
p = (Iapp,θ,η,E,(α,β,γ),k,(σnoise,dt_sim),outinj)

# Solve SDE problem (includes measurement noise)
# Important: Wₜ in the SDE is a Wiener process such that 
# Wₜ - Wₛ ∼ zero mean gaussian white noise with σ² = t-s
# When solving with Euler-Maruyama, dWₜ is approximated by ΔWₜ, which is 
# discrete-time zero-mean gaussian white noise with σ² = dt_sim
prob = SDEProblem(f_SDE!,g_SDE!,[z₀ 0],tspan,p)
sol = solve(prob,EM(),dt=dt_sim,saveat=dt_plot)
noise = sol[1,end,:]

# Solve ODE problem (no measurement noise)
# prob = ODEProblem(HH_ODE!,z₀,tspan,p)
# sol = solve(prob,Tsit5(),reltol=1e-6,abstol=1e-6,saveat=0.1,maxiters=1e5)

# Recover individual states
t = sol.t
v = sol[1,1,:];
w = sol[1,2:4,:];
v̂ = sol[1,5,:];
ŵ = sol[1,6:8,:];
θ̂ = sol[1,9:12,:];
η̂ = sol[1,13:15,:];
p = sol[1,15+1:15+49,:];
ψ₁ = sol[1,15+49+1:15+49+7,:];
ψ₂ = sol[1,15+49+7+1:15+49+7+3,:]

Vpower = 1/Tfinal*sum((v .- mean(v)).^2)*dt_plot
noisepower = 1/Tfinal*sum(noise.^2)*dt_plot
SNR = 10*log10(Vpower/noisepower)

if save_data
    writedlm(string("./../data/HH_voltages_a=",α,"_g=",γ,"_s=",σnoise,"_oi=",outinj,".txt"),  hcat(t[1:100:end],v[1:100:end],v̂[1:100:end]), " ")
    writedlm(string("./../data/HH_parameters_a=",α,"_g=",γ,"_s=",σnoise,"_oi=",outinj,".txt"),  hcat(t[1:100:end],θ̂[:,1:100:end]',η̂[:,1:100:end]'), " ")
end

## Plots
# pltv = plot(t,v)
# pltv = plot!(t,v̂,linecolor="red",linestyle= :dash)
# pltd = plot(t,permutedims(ŵ))
# pltd = plot(t,noise)
# pltP = plot(t,permutedims(p))

# gNa/c
plt1 = plot([0,Tfinal],[θ[1],θ[1]],linecolor="black",linestyle=:dash,labels="gNa/c")
plt1 = plot!(t,θ̂[1,:],linecolor="red")

# gK/c
plt2 = plot([0,Tfinal],[θ[2],θ[2]],linecolor="black",linestyle=:dash,labels="gK/c")
plt2 = plot!(t,θ̂[2,:],linecolor="red")

# gL/c
plt3 = plot([0,Tfinal],[θ[3],θ[3]],linecolor="black",linestyle=:dash,labels="gL/c")
plt3 = plot!(t,θ̂[3,:],linecolor="red")

plot(plt1,plt2,plt3,layout=(3,1))

## 1/c
plt8 = plot([0,Tfinal],[θ[4],θ[4]],linecolor="black",linestyle=:dash,labels="1/c")
plt8 = plot!(t,θ̂[4,:],linecolor="red")

### INTERNAL 
# 1 (m)
plt8 = plot([0,Tfinal],[η[1],η[1]],linecolor="black",linestyle=:dash,labels="vh_m")
plt8 = plot!(t,η̂[1,:],linecolor="red")

# 1 (h)
plt9 = plot([0,Tfinal],[η[2],η[2]],linecolor="black",linestyle=:dash,labels="vh_h")
plt9 = plot!(t,η̂[2,:],linecolor="red")

# 1 (n)
plt10 = plot([0,Tfinal],[η[3],η[3]],linecolor="black",linestyle=:dash,labels="vh_n")
plt10 = plot!(t,η̂[3,:],linecolor="red")

plot(plt8,plt9,plt10,layout=(3,1))