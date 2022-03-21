export
    AbstractProperties,
    Properties,
    properties

AbstractProperties = AbstractDict{String}
PropertyValue = Union{Real, Array}
Properties = Dict{String, PropertyValue}
DefaultProperties = Properties(
    "IntegrationTimestep" => 0.1, # s
    "Radius" => 0.5, # μm, typical value for E.coli
    "ReorientationRate" => 1.0, # 1/s, typical value for E.coli
    "RunState" => 1, # currently only used for run-reverse-flick motility; 1=reverse 0=flick
    "RotationalDiffusivity" => 0.0 # rad²/s
)

properties(x...) = merge(DefaultProperties, Properties(x...))
