using BacterialMotility
using LinearAlgebra
using Distributions
using Plots

include("/home/riccardo/Templates/juliaplots_templates/plots_template.jl")


function circleShape(h, k, r)
    θ = LinRange(0, 2π, 500)
    h .+ r.*sin.(θ), k .+ r.*cos.(θ)
end # function


Lbox = 900 # μm
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
custom_revflick!(bacterium) = reverse_flick!(bacterium, Normal(π,π/8), Normal(π/2,π/8))


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


R = 30.0
C = 2.0
Cbg = 0.01
f = Concentration_SteadyDiffusionSphericalSource_3D(R=R, C=C, Cbg=Cbg)
s_chemo = propertiesBrumley("ChemotacticPrecision" => 0.0)
s_nochemo = propertiesBrumley("ChemotacticPrecision" => 0.0, "ReceptorGain" => 0.0)

species = ["Random Walk", "Chemotaxis"]
pop1 = [BacteriumBrumley(id = species[1], r = randposition(), v = randvelocity(20.5), turn! = custom_revflick!, state = s_nochemo) for _ in 1:15]
pop2 = [BacteriumBrumley(id = species[2], r = randposition(), v = randvelocity(46.5), turn! = custom_revflick!, state = s_chemo) for _ in 1:15]
population = vcat(pop1, pop2)
num_bacteria = length(population)


nsteps = 500
save_every = 1
trajectories = zeros(nsteps÷save_every, num_bacteria, d)

runsim!(trajectories, nsteps, population, f; save_every=save_every)


speciescolor = Dict(species .=> 1:length(species))
linecolor = [speciescolor[bact.id] for _ in 1:1, bact in population]

function cfield(x,y)
    r = sqrt(x*x+y*y)
    r < R ? Cbg+C : Cbg+C*R/r
end # function
clims = (cfield(Lbox/3,Lbox/3), cfield(0,0))
xx = yy = -Lbox/3:4:Lbox/3
cc = cfield.(xx',yy)
for t in 2:2:nsteps
    p = plot(;plot_style(:Dark2)..., leg=false, axis=false, lims=(-Lbox/3,Lbox/3), aspect_ratio=1);
    heatmap!(p, xx, yy, cc, c=:bone, clims=clims)
    x = @view trajectories[1:t,:,1]
    xend = @view trajectories[t,:,1]
    y = @view trajectories[1:t,:,2]
    yend = @view trajectories[t,:,2]
    plot!(p, x, y, lab=false, lw=1, linealpha=range(0.1, 1; length=t), lc=linecolor)
    scatter!(p, xend, yend, lab=false, m=:c, msw=0, mc=linecolor|>permutedims)
    idx = lpad(t, 4, '0')
    savefig(p, "frame_$(idx).png")
end # for

#=
plot!(circleShape(0,0,R), seriestype=:shape, lw=0.0,
      c=length(species)+1, lab=false, fillalpha=0.5)

plot!(trajectories[:,:,1], trajectories[:,:,2],# trajectories[:,:,3],
      lab=false, lw=1,
      linealpha=range(0.2, 1; length=nsteps), lc=linecolor)

for i in 1:length(species)
    plot!([0.0], [0.0], lab=species[i], lc=speciescolor[species[i]], legend=:topright, legendfontsize=5)
end # for

plot!(axis=false, aspect_ratio=1, leg=false, lims=(-Lbox/3,Lbox/3))
=#
