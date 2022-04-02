export
    propertiesBrownBerg,
    memoryBrownBerg,
    affect_BrownBerg!,
    BacteriumBrownBerg

DefaultPropertiesBrownBerg = Properties(
    "AdaptationTime" => 1.0, # s
    "RunTimeUnbiased" => 0.67, # s
    "ReceptorBindingConstant" => 100.0, # μM
    "MotorGain" => 660.0 # s
)

propertiesBrownBerg(x...) = properties(DefaultPropertiesBrownBerg..., x...)

memoryBrownBerg(nsteps) = Dict(
    "NumSteps" => 0,
    "dPb_dt" => zeros(nsteps)
)

function affect_BrownBerg!(bacterium, ϕ, ∇ϕ)
    r = bacterium.r
    v = bacterium.v
    Uᵣ = dot(v,r) / norm(r)
    dt = bacterium.state["IntegrationTimestep"]
    tM = bacterium.state["AdaptationTime"]
    τ₀ = bacterium.state["RunTimeUnbiased"]
    KD = bacterium.state["ReceptorBindingConstant"]
    α = bacterium.state["MotorGain"]
    N = bacterium.memory["NumSteps"] + 1
    t = N*dt
    history_dPb_dt = bacterium.memory["dPb_dt"]
    dC_dt = Uᵣ * ∇ϕ
    dPb_dt = KD / (KD + ϕ)^2 * dC_dt
    history_dPb_dt[N] = dPb_dt
    weighted_dPb_dt = sum([history_dPb_dt[n] * exp(-(t-n*dt)/tM) for n in 1:N]) * dt/tM
    logτ = log(τ₀) + α*weighted_dPb_dt
    bacterium.state["ReorientationRate"] = min(exp(-logτ), 1.0/dt)
    bacterium.memory["NumSteps"] = N
end # function

@with_kw struct BacteriumBrownBerg{D} <: AbstractBacterium{D,Float64}
    id::String = ""
    r::MVector{D,Float64} = zeros(MVector{D,Float64})
    v::MVector{D,Float64} = zeros(MVector{D,Float64})
    run! = run!
    turn! = tumble!
    sense! = (b,f) -> sense_f_∇f!(b, f, affect_BrownBerg!)
    state = propertiesBrownBerg()
    memory = memoryBrownBerg(0)
end # struct
