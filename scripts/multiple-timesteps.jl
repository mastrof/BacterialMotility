using BacterialMotility
using LinearAlgebra
using Distributions
using Plots

include("/home/riccardo/Templates/juliaplots_templates/plots_template.jl")


function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r.*sin.(θ), k .+ r.*cos.(θ)
end # function


Lbox = 1.2e3 # μm
d = 2 # dimensionality

# convenience functions
randposition() = (rand(d) .- 1/2) .* (Lbox/2)
function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)

custom_tumble!(bacterium) = tumble!(bacterium, Degenerate(π/4))
custom_revflick!(bacterium) = reverse_flick!(bacterium, Normal(π,π/8), Uniform(1π/8,3π/8))


function save_pos!(traj, bacteria, N)
    t = findfirst(i -> traj[i,:,:] == zeros(size(traj[1,:,:])), 1:N)
    for i in 1:length(bacteria)
        traj[t,i,:] .= bacteria[i].r
    end # for
end # function

species = ["0.25", "0.1", "0.2"]
Δt = [0.25, 0.1, 0.2]
s = properties.("IntegrationTimestep" .=> Δt)
T = 2.0
nsteps = round.(Int, T ./ Δt)
N = maximum(nsteps)
trajectories = zeros(N, 3, d)
population = [Bacterium{d}(id = species[i], r = randposition(), v = randvelocity(30.0), run! = run!, turn! = tumble!, state = s[i]) for i in 1:3]
num_bacteria = length(population)

callback() = save_pos!(trajectories, population, N)

multistep!(population, T, EmptyField; callback = callback)


speciescolor = Dict(species .=> 1:length(species))
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

plot(;plot_style(:Dark2)...)

plot!(circleShape(0,0,R), seriestype=:shape, lw=0.0,
      c=length(species)+1, lab=false, fillalpha=0.5)

plot!(trajectories[:,:,1], trajectories[:,:,2],# trajectories[:,:,3],
      lab=false, lw=1,
      linealpha=range(0.2, 1; length=nsteps), lc=linecolor)

for i in 1:length(species)
    plot!([0.0], [0.0], lab=species[i], lc=speciescolor[species[i]], legend=:topright, legendfontsize=5)
end # for

plot!(aspect_ratio=1)

