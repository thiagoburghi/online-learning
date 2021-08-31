# using Plots
# using Random

# Perturbation amplitude
δₐ = 0.0
# Minimum perturbation 
δₘ = 0.0

# Kinetics perturbations
δ = δₐ*2*(0.5.-rand(Float64, (46,4)))
δ = δ .+ δₘ*sign.(δ)
δ = (1 .+ δ)
δ[:,1:2] = ones(46,2) # (no perturbation)

# plot(1:46,δ[:,3])

# Initial condition perturbations

ηₐ = 0.1
ηₘ = 0.2

# Initial parameter perturbations
η = ηₐ*2*(0.5.-rand(Float64, (2,9)))
η = η .+ ηₘ*sign.(η)
η = (1 .+ η)