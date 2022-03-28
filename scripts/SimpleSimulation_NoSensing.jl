#=
  Simulation of two-dimensional bacterial motility without chemosensing

  Three different bacterial species (motile patterns):
  1) Run-tumble with default reorientation distribution (Uniform(0,2π))
  2) Run-tumble with custom reorientation distribution (Degenerate(π/4))
  2) Run-reverse-flick with default reorientation distributions (reverse:Degenerate(π), flick:Degenerate(π/2))
  3) Run-reverse-flick with custom reorientation distributions (reverse:Normal(π,π/8), flick:Uniform(3π/8,5π/8))
=#

using BacterialMotility
using LinearAlgebra
using Distributions
using Plots

nclude("/home/riccardo/Templates/juliaplots_templates/plots_template.jl")

Lbox = 1e3 # μm
d = 2 # dimensionality
U = 30.0 # μm/s, speed, same for all species

# convenience functions
randposition() = (rand(d) .- 1/2) .* (Lbox/2)
function randvelocity(d)
    x = rand(d) .- 1/2
    x ./ norm(x)
end # function
randvelocity() = randvelocity(d) .* U

custom_tumble!(bacterium) = tumble!(bacterium, Degenerate(π/4))
custom_revflick!(bacterium) = reverse_flick!(bacterium, Normal(π,π/8), Uniform(3π/8,5π/8))

# create population composed by different bacteria in different numbers
pop1 = [Bacterium{d}(id = "run-tumble", r = randposition(), v = randvelocity(), run! = run!, turn! = tumble!) for _ in 1:5]
pop2 = [Bacterium{d}(id = "run-tumble_2", r = randposition(), v = randvelocity(), run! = run!, turn! = custom_tumble!) for _ in 1:2]
pop3 = [Bacterium{d}(id = "run-reverse-flick", r = randposition(), v = randvelocity(), run! = run!, turn! = reverse_flick!) for _ in 1:3]
pop4 = [Bacterium{d}(id = "run-reverse-flick_2", r = randposition(), v = randvelocity(), run! = run!, turn! = custom_revflick!) for _ in 1:4]
population = vcat(pop1, pop2, pop3, pop4)
num_bacteria = length(population)

f = EmptyField() # no sensing behavior

dt = 0.05 # s
nsteps = 750
trajectories = zeros(nsteps, num_bacteria, d)
for t in 1:nsteps
    step!(population, dt, f)
    for n in 1:num_bacteria
        trajectories[t,n,:] .= population[n].r
    end # for
end # for

species = ["run-tumble", "run-tumble_2", "run-reverse-flick", "run-reverse-flick_2"]
speciescolor = Dict(species .=> 1:4)
specieslabels = ["Run-Tumble", "Custom Run-Tumble", "Run-Reverse-Flick", "Custom Run-Reverse-Flick"]
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

plot(;plot_style(:Dark2)...)
plot!(trajectories[:,:,1], trajectories[:,:,2], lab=false,
      lw=1, linealpha=range(0.2, 1; length=nsteps), lc=linecolor)
for i in 1:length(species)
    plot!([0.0], [0.0], lab=specieslabels[i], lc=speciescolor[species[i]], legend=:topright, legendfontsize=5)
end # for
plot!(ticks=false)
