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
d = 3 # dimensionality
U = 30.0 # μm/s, speed, same for all species

# convenience functions
randposition() = (rand(d) .- 1/2) .* (Lbox/2)
function randvelocity(d,U)
    x = rand(d) .- 1/2
    x ./ norm(x) .* U
end # function
randvelocity(U) = randvelocity(d,U)

custom_tumble!(bacterium) = tumble!(bacterium, Degenerate(π/4))
custom_revflick!(bacterium) = reverse_flick!(bacterium, Normal(π,π/8), Uniform(1π/8,3π/8))


function runsim!(traj, steps, population, f; save_every=1)
    i = 1
    for t in 1:nsteps
        step!(population, f)
        if (t-1) % save_every == 0
            for n in 1:num_bacteria
                traj[i,n,:] .= population[n].r
            end # for
            i += 1
        end # if
    end # for
end # function


R = 100.0
C = 6.0
Cbg = 0.0
f = Concentration_SteadyDiffusionSphericalSource_3D(R=R, C=C, Cbg=Cbg)
s = propertiesBrumley("ChemotacticPrecision" => 0.0)


# create population composed by different bacteria in different numbers
species = ["RT (Chemo,π/4)", "RRF (Brumley)", "RRF (No chemo)", "RRF (Boh)"]
pop1 = [BacteriumBrumley(id = species[1], r = randposition(), v = randvelocity(30.0), turn! = custom_tumble!, state=s) for _ in 1:3]
pop2 = [BacteriumBrumley(id = species[2], r = randposition(), v = randvelocity(46.5), turn! = reverse_flick!, state=s) for _ in 1:3]
pop3 = [Bacterium{3}(id = species[3], r = randposition(), v = randvelocity(46.5), run! = run!, turn! = reverse_flick!) for _ in 1:3]
pop4 = [Bacterium{3}(id = species[4], r = randposition(), v = randvelocity(26.5), run! = run!, turn! = custom_revflick!) for _ in 1:3]
population = vcat(pop1, pop2, pop3, pop4)
num_bacteria = length(population)


nsteps = 500
save_every = 1
trajectories = zeros(nsteps÷save_every, num_bacteria, d)

runsim!(trajectories, nsteps, population, f; save_every=save_every)


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

