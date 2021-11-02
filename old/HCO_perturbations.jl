# using Plots
# using Random

δₐ = 0.1
δₘ = 0.1
ηₐ = 0.1
ηₘ = 0.2

# Kinetics perturbations
δ = δₐ*2*(0.5.-rand(Float64, (46,2)))
δ = δ .+ δₘ*sign.(δ)
δ = (1 .+ δ)

if perturbed == 1
    # adaptive observer is perturbed
    δ = hcat(ones(46,2), δ)
else
    # true model is perturbed
    δ = hcat(δ, ones(46,2))
end

# plot(1:46,δ[:,3])

# Initial parameter perturbations
η = ηₐ*2*(0.5.-rand(Float64, (2,5)))
η = η .+ ηₘ*sign.(η)
η = (1 .+ η)