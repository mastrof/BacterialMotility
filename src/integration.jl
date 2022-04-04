export
    step!, multistep!

@doc raw"""
    step!(bacterium::B, f::F; callback) where {B<:AbstractBacterium,F<:AbstractField}
    step!(population::T, f::F; callback_inner, callback_outer) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}

Perform a single integration step, composed as:
    1) run
    2) sensing
    3) reorientation attempt
    4) callback
"""
function step!(bacterium::B, f::F=EmptyField(); callback=dummy, kwargs...) where {B<:AbstractBacterium,F<:AbstractField}
    dt = bacterium.state["IntegrationTimestep"]
    bacterium.run!(bacterium, dt)
    bacterium.sense!(bacterium, f; kwargs...)
    ν = bacterium.state["ReorientationRate"]
    if rand() < ν * dt
        bacterium.turn!(bacterium)
    end # if
    callback(bacterium, f; kwargs...)
end # function

function step!(population::T, f::F=EmptyField(); callback=dummy, kwargs...) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}
    for bacterium in population
        step!(bacterium, f; callback=callback, kwargs...)
    end # for
end # function

