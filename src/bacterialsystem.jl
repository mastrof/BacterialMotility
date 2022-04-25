export
    BacterialSystem,
    integrate!

@with_kw struct BacterialSystem
    clock = [0] # Vector{Int}
    field = EmptyField() # AbstractField
    boundary_conditions! = dummy
    callback_inner = dummy # Function
    callback_outer = dummy # Function
    population # AbstractVector{AbstractBacterium}
end # struct


function integrate!(bs::BacterialSystem, nsteps::Int; kwargs...)
    for n in 1:nsteps
        step!(bs.population, bs.field;
              bcs! = bs.boundary_conditions!, callback=bs.callback_inner, kwargs...)
        bs.clock[1] += 1
        bs.callback_outer(bs; kwargs...)
    end # for
end # function

function integrate!(bs::BacterialSystem, T::Float64; kwargs...)
    num_bacteria = length(bs.population)
    Δt = [b.state["IntegrationTimestep"] for b in bs.population]
    nsteps = T ./ Δt # steps for each bacterium to reach time T
    N = round(Int, maximum(nsteps))
    δ = nsteps ./ N
    accumulator = zeros(num_bacteria)
    for n in 1:N
        for i in 1:num_bacteria
            accumulator[i] += δ[i]
            if accumulator[i] ≥ 1
                step!(bs.population[i], bs.field;
                      bcs! = bs.boundary_conditions!, callback=bs.callback_inner, kwargs...)
                accumulator[i] -= 1
            end # if
        end # for
        bs.clock[1] += 1
        bs.callback_outer(bs; kwargs...)
    end # for
end # function
