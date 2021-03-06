export
    Bacterium

@doc raw"""
    struct Bacterium{D,T} <: AbstractBacterium{D,T}

Concrete type for bacterium with sensing capabilities and internal state variables.
    - `id::String`: identifier (e.g. species, strain...)
    - `r::MVector{D,T}`: position vector
    - `v::MVector{D,T}`: velocity vector
    - `run!`: function defining displacement, acts in place on r
    - `turn!`: function defining reorientation, acts in place on v
    - `sense!`: function defining sensing behavior
    - `state`: integration parameters, bacterium properties

All arguments and parameters are optional, except the dimensionality parameter D which must be specified.
If parameter T is not given, it is assumed to be T=Float64.
If r and v are not defined, they are initialized to zeros(T,D).
If run!, turn!, sense! are not defined, they are assigned to a dummy function that does nothing.
Some standard properties are provided by default for state (see `properties`).

run! and turn! should only take an instance of AbstractBacterium as input; sense! should take AbstractBacterium and AbstractField.
"""
@with_kw struct Bacterium{D,T} <: AbstractBacterium{D,T}
    id::String = "" # identifier
    r::MVector{D,T} = zeros(MVector{D,T}) # position vector
    v::MVector{D,T} = zeros(MVector{D,T})# velocity vector
    run! = run! # displacement
    turn! = dummy_turn # reorientation
    sense! = dummy_sense # kinetic/tactic behavior
    state = properties()
end # state

Bacterium{D}(; kwargs...) where D = Bacterium{D,Float64}(; kwargs...)
