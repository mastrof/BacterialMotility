export
    AbstractBacterium, Bacterium

@doc raw"""
    AbstractBacterium{D,T}

General interface for bacterium object.
    - `D`: dimensionality of the space in which the bacterium lives
    - `T`: type of the space in which the bacterium lives
"""
abstract type AbstractBacterium{D,T} end

@doc raw"""
    struct Bacterium{D,T} <: AbstractBacterium{D,T}

Concrete type for bacterium with sensing capabilities and one (or multiple) internal state variable.
    - `id::String`: identifier (e.g. species, strain...)
    - `r::MVector{D,T}`: position vector
    - `v::MVector{D,T}`: velocity vector
    - `run!`: function defining displacement, acts in place on r
    - `turn!`: function defining reorientation, acts in place on v
    - `sense!`: function defining sensing behavior
    - `state`: integration parameters, bacterium properties

All arguments and parameters are optional.
If parameters D and T are not given, they are assumed to be D=3 and T=Float64.
If r and v are not defined, they are initialized to zeros(T,D).
If run!, turn!, sense! are not defined, they are assigned to a dummy function that does nothing.
Some standard properties are provided by default for state.

run! and turn! should only take an instance of AbstractBacterium as input; sense! should take AbstractBacterium and AbstractField.
"""
@with_kw struct Bacterium{D,T} <: AbstractBacterium{D,T}
    id::String = "" # identifier
    r::MVector{D,T} = zeros(MVector{D,T}) # position vector
    v::MVector{D,T} = zeros(MVector{D,T})# velocity vector
    run! = dummy # displacement
    turn! = dummy # reorientation
    sense! = dummy # kinetic/tactic behavior
    state = properties()
end # state

Bacterium{D}(; kwargs...) where D = Bacterium{D,Float64}(; kwargs...)
Bacterium(; kwargs...) = Bacterium{3,Float64}(; kwargs...)
