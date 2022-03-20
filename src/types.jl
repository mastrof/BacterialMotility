export
    AbstractProperties, Properties, properties,
    AbstractBacterium, Bacterium


AbstractProperties = AbstractDict{String}
PropertyValue = Union{Real, Array}
Properties = Dict{String, PropertyValue}
DefaultProperties = Properties(
    "Radius" => 0.5, # μm, typical value for E.coli
    "ReorientationRate" => 1.0, # 1/s, typical value for E.coli
    "RunState" => 1, # currently only used for run-reverse-flick motility; 1=reverse 0=flick
    "RotationalDiffusivity" => 0.0 # rad²/s
)
properties(x...) = merge(DefaultProperties, Properties(x...))

function dummy(x...)
    ;
end

@doc raw"""
    AbstractBacterium{D,T,F<:Function,S<:AbstractProperties}

General interface for bacterium object.
    - `D`: dimensionality of the space in which the bacterium lives
    - `T`: type of the space in which the bacterium lives
"""
abstract type AbstractBacterium{D,T,F<:Function,S<:AbstractProperties} end

@doc raw"""
    struct Bacterium{D,T,F,S} <: AbstractBacterium{D,T,F,S}

Concrete type for bacterium with sensing capabilities and one (or multiple) internal state variable.
    - `id::String`: identifier (e.g. species, strain...)
    - `r::MVector{D,T}`: position vector
    - `v::MVector{D,T}`: velocity vector
    - `run!::F`: function defining displacement, acts in place on r
    - `turn!::F`: function defining reorientation, acts in place on v
    - `sense!::F`: function defining sensing behavior
    - `state::S`: properties, geometry, internal state variables

All arguments and parameters are optional.
If parameters D and T are not given, they are assumed to be D=3 and T=Float64. F<:Function and S<:AbstractProperties.
If r and v are not defined, they are initialized to zeros(T,D).
If run!, turn!, sense! are not defined, they are assigned to a dummy function that does nothing.
Some standard properties are provided by default for state.

Ideally, `run!` should take two parameters, bacterium (self) and an externally provided integration timestep;
turn! should only take bacterium (self) as input; sense! should take bacterium and an instance of AbstractField.
"""
@with_kw struct Bacterium{D,T,F,S} <: AbstractBacterium{D,T,F,S}
    id::String = "" # identifier
    r::MVector{D,T} = zeros(MVector{D,T}) # position vector
    v::MVector{D,T} = zeros(MVector{D,T})# velocity vector
    run!::F = dummy # displacement
    turn!::F = dummy # reorientation
    sense!::F = dummy # kinetic/tactic behavior
    state::S = properties()
end # state

Bacterium{D,T}(; kwargs...) where {D,T} = Bacterium{D,T,Function,Properties}(; kwargs...)
Bacterium{D}(; kwargs...) where D = Bacterium{D,Float64,Function,Properties}(; kwargs...)
Bacterium(; kwargs...) = Bacterium{3,Float64,Function,Properties}(; kwargs...)
