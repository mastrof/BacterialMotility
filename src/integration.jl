export
    step!

@doc raw"""
    step!(bacterium::B, f::F) where {B<:AbstractBacterium,F<:AbstractField}
    step!(bacteria::T, f::F) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}

Perform a single integration step, composed as:
    1) run
    2) sensing
    3) reorientation attempt
"""
function step!(bacterium::B, f::F=EmptyField) where {B<:AbstractBacterium,F<:AbstractField}
    dt = bacterium.state["IntegrationTimestep"]
    bacterium.run!(bacterium, dt)
    bacterium.sense!(bacterium, f)
    if rand() < bacterium.state["ReorientationRate"] * dt
        bacterium.turn!(bacterium)
    end # if
end # function


function step!(bacteria::T, f::F=EmptyField) where {B<:AbstractBacterium,T<:AbstractVector{B},F<:AbstractField}
    for bacterium in bacteria
        step!(bacterium, f)
    end # for
end # function

