export
    propertiesBrumley,
    affect_Brumley!,
    BacteriumBrumley

DefaultPropertiesBrumley = Properties(
    "ReceptorGain" => 50.0, # 1/μM
    "MotorGain" => 50.0, # 
    "ChemotacticPrecision" => 6.0, # 
    "InternalState" => 0.0, # 
    "SensoryTime" => 0.1, # s
    "AdaptationTime" => 1.3, # s
    "RunTimeUnbiased" => 0.45, # s
    "NutrientDiffusivity" => 608.0, # μm²s
    "ReorientationRate" => 1/0.45, # 1/s
)

propertiesBrumley(x...) = properties(DefaultPropertiesBrumley..., x...)

function affect_Brumley!(bacterium, ϕ, ∇ϕ)
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
    μ = Uᵣ * ∇ϕ # mean concentration gradient at current position
    σ = Π * sqrt(3.0*ϕ / (π*a*Dc*T³)) # sensing noise at current position
    M = rand(Normal(μ,σ)) # gradient measurement
    S = (1-α)*κ*tM*M + α*S
    bacterium.state["InternalState"] = S
    ν = (1.0 + exp(-Γ*S)) / (2.0*τ₀)
    bacterium.state["ReorientationRate"] = min(ν, 1.0/dt)
end # function


@with_kw struct BacteriumBrumley <: AbstractBacterium{3,Float64}
    id::String = ""
    r::MVector{3,Float64} = zeros(MVector{3,Float64})
    v::MVector{3,Float64} = zeros(MVector{3,Float64})
    run! = run!
    turn! = reverse_flick!
    sense! = (b,f) -> sense_f_∇f!(b, f, affect_Brumley!)
    state = propertiesBrumley()
end # struct

