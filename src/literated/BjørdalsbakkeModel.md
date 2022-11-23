```@meta
EditURL = "<unknown>/BjørdalsbakkeModel.jl"
```

````@example BjørdalsbakkeModel
include("../../src/CirculationModels.jl")

using ModelingToolkit, DifferentialEquations, Plots
using .CirculationModels
````

# Simploified two-chamber model

This follows Bjørdalsbakke et al.

````@example BjørdalsbakkeModel
include("BjørdalsbakkeParam.jl");

nstep = 1000
t = LinRange(0, τ, nstep)
#
# Start Modelling
kLV = 1 ./ maximum((t ./ Tau1fLV).^n1LV ./ (1 .+ (t ./ Tau1fLV).^n1LV) .* 1 ./ (1 .+ (t ./ Tau2fLV).^n2LV))

#
# Named Elements
@parameters t
@named LV = DHChamber(V₀ = 0.0, Eₘₐₓ=Eₘₐₓ, Eₘᵢₙ=Eₘᵢₙ, n₁=n1LV, n₂=n2LV, τ = τ, τ₁=Tau1fLV, τ₂=Tau2fLV, k = kLV, Eshift=0.0, Ev=Inf)

@named AV = ResistorDiode(R=Zao)
@named MV = ResistorDiode(R=Rmv)

@named Rs = Resistor(R=Rs)


@named Csa = Compliance(C=Csa)
@named Csv = Compliance(C=Csv)

@named ground = Ground()

# Connect the system
circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Csa.in)
    connect(Csa.out, Rs.in)
    connect(Rs.out, Csv.in)
    connect(Csv.out, MV.in)
    connect(MV.out, LV.in)
]


# Compose the whole ODE system
@named _circ_model = ODESystem(circ_eqs, t)

#
@named circ_model = compose(_circ_model,
                          [LV, AV, MV, Rs, Csa, Csv, ground])



# And simplify it
circ_sys = structural_simplify(circ_model)


prob = ODEProblem(circ_sys, [MCFP, MCFP, MCFP], (0.0, 20.0))

@time sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)


p1 = plot(sol, idxs=[LV.p,  Csa.in.p], tspan=(16 * τ, 17 * τ), xlabel = "Time [s]", ylabel = "Pressure [mmHg]",  hidexaxis = nothing) # Make a line plot
p2 = plot(sol, idxs=[LV.V], tspan=(16 * τ, 17 * τ),xlabel = "Time [s]", ylabel = "Volume [ml]",  linkaxes = :all)
p3 = plot(sol, idxs=[Csa.in.q,Csv.in.q], tspan=(16 * τ, 17 * τ),xlabel = "Time [s]", ylabel = "Flow rate [ml/s]", linkaxes = :all)
plot(p1, p2, p3, layout=@layout([a; b c]), legend = true)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

