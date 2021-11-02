using Plots: print
using DifferentialEquations, Random, Plots, LinearAlgebra, DelimitedFiles

# Flag for saving data to .txt files 
save_data = false

# Models
include("./HCO_kinetics.jl")
include("./RB_odes.jl")

# Definition of the timespan for numerical integration
Tfinal = 4000.0
tspan = (0.0,Tfinal)
tplot = 0.0:0.1:Tfinal

# Input currents:
# T = period of the pulse train
# L = length of each pulse
# t₀ = time of the pulse within the first period
# a₀ = baseline input current
# a₁ = input current during pulse
function pulsetrain(t,t₀,L,a₀,a₁,T)
    for k=0:floor(Int,Tfinal/T)
        if (t > t₀+k*T) & (t <= t₀+L+k*T)
            return a₁
        end
    end
    return a₀
end
Iapp(t) = pulsetrain(t,500,500,-1,-2,1000)

# Desired Calcium maximal conductance
ḡCa(t) = 0.125

# True system parameters
C = 1.
gNa = 60.       # Sodium current maximal conductance
gK = 40.        # Potassium current maximal conductance
gCa = 0.18      # Calcium current maximal conductance;  0.14
gL = 0.035      # Leak maximal conductance

# True system initial conditions
V₀= -64
x₀ = [V₀ σNa_m(V₀) σNa_h(V₀)  σK_m(V₀) σCa_m(V₀) σCa_h(V₀)]

## OPEN-LOOP TRUE SYSTEM

p₀ = (Iapp,C,gNa,gK,gCa,gL)
prob₀ = ODEProblem(RB_ode!,x₀,tspan,p₀)
sol₀ = solve(prob₀,Tsit5(),saveat=0.1,reltol=1e-6,abstol=1e-6)
t = sol₀.t
v₀ = sol₀[1,1,:]

if save_data
    writedlm("./../data/RB_openloop.txt",  hcat(t,v₀,Iapp.(t)), " ")
end


plt1 = plot(t,v₀)
plt2 = plot(t,Iapp.(t),linecolor="gray")
plot(plt1,plt2,layout=(2,1))

## ADAPTIVE CONTROL

# Observer initial conditions
x̂₀ = [-50. zeros(1,5)]
gNa₀ = 80
gK₀ = 80
gCa₀ = 1
gL₀ = 1
θ̂₀ = [gNa₀ gK₀ gCa₀ gL₀]
P₀ = Matrix(I, 4, 4)
Ψ₀ = zeros(1,4)

# Observer parameters
α₁ = 0.001 # 0.0001
α₂ = 1.
γ = 5

# True system and observer initial conditions and parameters
z₀ = [x₀ x̂₀ θ̂₀ reshape(P₀,1,16) Ψ₀]
p = (Iapp,C,gNa,gK,gCa,gL,(α₁,α₂,γ),ḡCa)

prob = ODEProblem(RB_controlled_ode!,z₀,tspan,p)
sol = solve(prob,Tsit5(),saveat=0.1,reltol=1e-6,abstol=1e-6) #,reltol=1e-6,abstol=1e-6,AutoTsit5(Rosenbrock23()),Tsit5()

t = sol.t 
v = sol[1,1,:];
v̂ = sol[1,7,:];
θ̂ = sol[1,13:16,:]';
ĝNa = θ̂[:,1];       
ĝK = θ̂[:,2];       
ĝCa = θ̂[:,3];       
ĝL = θ̂[:,4];       

if save_data
    writedlm("./../data/RB_closedloop.txt",  hcat(t,v,Iapp.(t),ĝCa,ḡCa.(t)), " ")
end

# VOLTAGES
plt0 = plot(t, v,linewidth=1.5,legend=false)
plt0 = plot!(t, v̂,linewidth=1.5,linecolor="red",linestyle= :dash,legend=false)
# savefig(voltages,"./HCO.png")

# Na CONDUCTANCES 
plt1 = plot(t,ĝNa,labels="ĝNa",linecolor="red")
plt1 = plot!([0,Tfinal],[gNa,gNa],linecolor="black",linestyle= :dash,labels="gNa")

# K CONDUCTANCES 
plt2 = plot(t,ĝK,labels="ĝK",linecolor="red") 
plt2 = plot!([0,Tfinal],[gK,gK],linecolor="black",linestyle= :dash,labels="gK")

# Ca CONDUCTANCES 
plt3 = plot(t,ĝCa,labels="ĝCa",linecolor="red") 
plt3 = plot!([0,Tfinal],[gCa,gCa],linecolor="black",linestyle= :dash,labels="gCa")

plot(plt0,plt1,plt2,plt3,layout=(4,1))