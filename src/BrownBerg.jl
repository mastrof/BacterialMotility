export
    propertiesBrownBerg,
    affect_BrownBerg!,
    BacteriumBrownBerg

function propertiesBrownBerg(x...)
    DefaultPropertiesBrownBerg = Properties(
        "AdaptationTime" => 1.0, # s
        "ReceptorBindingConstant" => 100.0, # μM
        "MotorGain" => 660.0, # s
        "RunTimeUnbiased" => 0.67, # s
        "ReorientationRate" => 1.0/0.67, # 1/s
        "InternalState" => 0.0 #
    )
    properties(DefaultPropertiesBrownBerg..., x...)
end # function

function affect_BrownBerg!(bacterialsystem, i, ϕ, ∇ϕ, ∂ₜϕ)
    bacterium = bacterialsystem.population[i]
    r = bacterium.r
    v = bacterium.v
    dt = bacterium.state["IntegrationTimestep"]
    tM = bacterium.state["AdaptationTime"]
    τ₀ = bacterium.state["RunTimeUnbiased"]
    KD = bacterium.state["ReceptorBindingConstant"]
    α = bacterium.state["MotorGain"]
    S = bacterium.state["InternalState"] # weighted dPb/dt @ previous step
    
    dC_dt = dot(v,∇ϕ) + ∂ₜϕ # new gradient measurement
    M = KD / (KD + ϕ)^2 * dC_dt # dPb/dt from new measurement
    S = M*dt/tM + S*exp(-dt/tM) # weighted dPb/dt
    bacterium.state["InternalState"] = S
    logτ = log(τ₀) + α*S
    bacterium.state["ReorientationRate"] = min(exp(-logτ), 1.0/dt)
end # function

@with_kw struct BacteriumBrownBerg{D} <: AbstractBacterium{D,Float64}
    id::String = ""
    r::MVector{D,Float64} = zeros(MVector{D,Float64})
    v::MVector{D,Float64} = zeros(MVector{D,Float64})
    run! = run!
    turn! = tumble!
    sense! = (bs, i) -> sense!(bs, i, affect_BrownBerg!)
    state = propertiesBrownBerg()
end # struct
