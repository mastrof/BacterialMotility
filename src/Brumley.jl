export
    propertiesBrumley,
    affect_Brumley!,
    BacteriumBrumley

function propertiesBrumley(x...)
    DefaultPropertiesBrumley = Properties(
        "IntegrationTimestep" => 0.1, # s
        "SensoryTime" => 0.1, # s
        "ReceptorGain" => 50.0, # 1/μM
        "MotorGain" => 50.0, # 
        "ChemotacticPrecision" => 6.6, # 
        "InternalState" => 0.0, # 
        "AdaptationTime" => 1.3, # s
        "RunTimeUnbiased" => 0.45, # s
        "NutrientDiffusivity" => 608.0, # μm²s
        "ReorientationRate" => 1/0.45, # 1/s
        "RotationalDiffusivity" => 0.0349, # rad²/s
    )
    i = findfirst("SensoryTime" .== first.(x))
    j = findfirst("IntegrationTimestep" .== first.(x))
    if isnothing(i) && isnothing(j)
        properties(DefaultPropertiesBrumley..., x...)
    elseif isnothing(i)
        dt = x[j][2]
        properties(DefaultPropertiesBrumley..., x...,
                   "SensoryTime" => dt)
    elseif isnothing(j)
        T = x[i][2]
        properties(DefaultPropertiesBrumley..., x...,
                   "IntegrationTimestep" => T)
    else
        dt = x[j][2]
        T = x[i][2]
        if dt == T
            properties(DefaultPropertiesBrumley..., x...)
        else
            @warn "IntegrationTimestep must be equal to SensoryTime." *
                "\n" * "IntegrationTimestep will be set to the value $T."
            properties(DefaultPropertiesBrumley..., x...,
                       "IntegrationTimestep" => T)
        end # if
    end # if
end # function

function affect_Brumley!(bacterium, ϕ, ∇ϕ, ∂ₜϕ)
    r = bacterium.r
    v = bacterium.v
    Uᵣ = dot(v,r) / norm(r)
    dt = bacterium.state["IntegrationTimestep"]
    a = bacterium.state["Radius"]
    τ₀ = bacterium.state["RunTimeUnbiased"]
    tM = bacterium.state["AdaptationTime"]
    Dc = bacterium.state["NutrientDiffusivity"]
    κ = bacterium.state["ReceptorGain"]
    Γ = bacterium.state["MotorGain"]
    Π = bacterium.state["ChemotacticPrecision"]
    T = bacterium.state["SensoryTime"]
    S = bacterium.state["InternalState"]
    α = exp(-T/tM) # memory persistence factor
    T³ = T*T*T
    μ = Uᵣ*∇ϕ + ∂ₜϕ # mean concentration gradient at current position
    σ = Π * sqrt(3.0*ϕ / (π*a*Dc*T³)) # sensing noise at current position
    M = rand(Normal(μ,σ)) # gradient measurement
    S = (1-α)*κ*tM*M + α*S
    bacterium.state["InternalState"] = S
    ν = (1.0 + exp(-Γ*S)) / (2.0*τ₀)
    bacterium.state["ReorientationRate"] = min(ν, 1.0/dt)
end # function


@with_kw struct BacteriumBrumley{D} <: AbstractBacterium{D,Float64}
    id::String = ""
    r::MVector{D,Float64} = zeros(MVector{D,Float64})
    v::MVector{D,Float64} = zeros(MVector{D,Float64})
    run! = run!
    turn! = reverse_flick!
    sense! = (b, f; kwargs...) -> sense!(b, f, affect_Brumley!; kwargs...)
    state = propertiesBrumley()
end # struct

