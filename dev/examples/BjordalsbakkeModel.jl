# # Importing the required packages

using CirculatorySystemModels
using CirculatorySystemModels.ModelingToolkit
using CirculatorySystemModels.DifferentialEquations
using Plots
using DisplayAs

# # A simple single-chamber model
#
# ![Single chamber, closed-loop, lumped parameter model of the systemic circulation and the left ventricle. The circuit equivalent formulation of the model is depicted, with the pressures of each compartment, as well as most of the mechanical parameters. The model describes three compartments: the left ventricular, arterial and venous compartments. 𝑃𝑡ℎ is the intrathoracic pressure, 𝑃𝑙𝑣 is the left ventricular pressure and 𝐸𝑙𝑣(𝑡) indicates the left ventricular elastance function.](./BjordalsbakkeModelSketch.png)
#
#
# This follows Bjørdalsbakke et al.
#
# Bjørdalsbakke, N.L., Sturdy, J.T., Hose, D.R., Hellevik, L.R., 2022. Parameter estimation for closed-loop lumped parameter models of the systemic circulation using synthetic data. Mathematical Biosciences 343, 108731. https://doi.org/10.1016/j.mbs.2021.108731
#
#
# Changes from the published version above:
#
# - Capacitors are replaced by compliances. These are identical to capacitors, but have an additional parameter, the unstrained volume $V_0$, which allows for realistic blood volume modelling.
#   Compliances have an inlet and an oulet in line with the flow, rather than the ground connector of the dangling capacitor.
# - The aortic resistor is combined with the valve (diode) in the `ResistorDiode` element.
#
# [Jupyter Notebook](./BjordalsbakkeModel.ipynb)

# ## Define the parameters
#
# All the parameters are taken from table 1 of [Bjørdalsbakke2022].
#

# Cycle time in seconds
#
τ = 0.85

# Double Hill parameters for the ventricle
#
Eₘᵢₙ = 0.03
Eₘₐₓ = 1.5
n1LV    = 1.32;
n2LV    = 21.9;
Tau1fLV = 0.303 * τ;
Tau2fLV = 0.508 * τ

# Resistances and Compliances
#
Rs = 1.11
Csa = 1.13
Csv = 11.0

# Valve parameters
#
# Aortic valve basic
Zao = 0.033
# Mitral valve basic
Rmv = 0.006

# Inital Pressure (mean cardiac filling pressure)
MCFP = 7.0

# ## Calculating the additional `k` parameter
#
# The ventricle elastance is modelled as:
#
# $$E_{l v}(t)=\left(E_{\max }-E_{\min }\right) e(t)+E_{\min }$$
#
# where $e$ is a double-Hill function, i.e., two Hill-functions, which are multiplied by each other:
#
# $$e(\tau)= k \times \frac{\left(\tau / \tau_1\right)^{n_1}}{1+\left(\tau / \tau_1\right)^{n_1}} \times \frac{1}{1+\left(\tau / \tau_2\right)^{n_2}}$$
#
# and $k$ is a scaling factor to assure that $e(t)$ has a maximum of $e(t)_{max} = 1$:
#
# $$k = \max \left(\frac{\left(\tau / \tau_1\right)^{n_1}}{1+\left(\tau / \tau_1\right)^{n_1}} \times \frac{1}{1+\left(\tau / \tau_2\right)^{n_2}} \right)^{-1}$$
#

nstep = 1000
t = LinRange(0, τ, nstep)

kLV = 1 / maximum((t ./ Tau1fLV).^n1LV ./ (1 .+ (t ./ Tau1fLV).^n1LV) .* 1 ./ (1 .+ (t ./ Tau2fLV).^n2LV))


# ## Set up the model elements
#
# Set up time as a variable `t`
#
@variables t

# Heart is modelled as a single chamber (we call it `LV` for "Left Ventricle" so the model can be extended later, if required):
#
@named LV = DHChamber(V₀ = 0.0, Eₘₐₓ=Eₘₐₓ, Eₘᵢₙ=Eₘᵢₙ, n₁=n1LV, n₂=n2LV, τ = τ, τ₁=Tau1fLV, τ₂=Tau2fLV, k = kLV, Eshift=0.0)

# The two valves are simple diodes with a small resistance
# (resistance is needed, since perfect diodes would connect two elastances/compliances, which will lead to unstable oscillations):
#
@named AV = ResistorDiode(R=Zao) 
@named MV = ResistorDiode(R=Rmv) 

# The main components of the circuit are 1 resistor `Rs` and two compliances for systemic arteries `Csa`,
# and systemic veins `Csv` (names are arbitrary).
#
@named Rs = Resistor(R=Rs)

@named Csa = Compliance(C=Csa)
@named Csv = Compliance(C=Csv)

# We also need to define a base pressure level, which we use the `Ground` element for:
#
@named ground = Ground(P=0)

# ## Build the system
#
# ### Connections
#
# The system is built using the `connect` function. `connect` sets up the Kirchhoff laws:
#
# - pressures are the same in all connected branches on a connector
# - sum of all flow rates at a connector is zero
#
# The resulting set of Kirchhoff equations is stored in `circ_eqs`:
#
circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Csa.in)
    connect(Csa.out, Rs.in)
    connect(Rs.out, Csv.in)
    connect(Csv.out, MV.in)
    connect(MV.out, LV.in)
]

# ### Add the component equations
#
# In a second step, the system of Kirchhoff equations is completed by the component equations (both ODEs and AEs), resulting in the full, overdefined ODE set `circ_model`.
#
# _Note: we do this in two steps._
#
@named _circ_model = ODESystem(circ_eqs, t)

@named circ_model = compose(_circ_model,
                          [LV, AV, MV, Rs, Csa, Csv, ground])

# ### Simplify the ODE system
#
# The crucial step in any acausal modelling is the sympification and reduction of the OD(A)E system to the minimal set of equations. ModelingToolkit.jl does this in the `structural_simplify` function.
#
circ_sys = structural_simplify(circ_model)

# `circ_sys` is now the minimal system of equations. In this case it consists of 3 ODEs for the ventricular volume and the systemic and venous pressures.
#
# _Note: this reduces and optimises the ODE system. It is, therefore, not always obvious, which states it will use and which it will drop. We can use the `states` and `observed` function to check this. It is recommended to do this, since small changes can reorder states, observables, and parameters._
#
# States in the system are now:
states(circ_sys)

# Observed variables - the system will drop these from the ODE system that is solved, but it keeps all the algebraic equations needed to calculate them in the system object, as well as the `ODEProblem` and solution object - are:
observed(circ_sys)

# And the parameters (these could be reordered, so check these, too):
parameters(circ_sys)

# ### Define the ODE problem
#
# First defined initial conditions `u0` and the time span for simulation:
#
# _Note: the initial conditions are defined as a parameter map, rather than a vector, since the parameter map allows for changes in order. This map can include non-existant states (like `LV.p` in this case), which allows for exchanging the ventricle for one that's defined in terms of $dp/dt$)._

u0 = [
        LV.p => MCFP
        LV.V => MCFP/Eₘᵢₙ
        Csa.p => MCFP
        Csv.p => MCFP
        ]

# 
        
tspan = (0, 20)

# in this case we use the mean cardiac filling pressure as initial condition, and simulate 20 seconds.
#
# Then we can define the problem:
#
prob = ODEProblem(circ_sys, u0, tspan)

# ## Simulate
#
# The ODE problem is now in the MTK/DifferentialEquations.jl format and we can use any DifferentialEquations.jl solver to solve it:
#
sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12);

# ## Results
p1 = plot(sol, idxs=[LV.p,  Csa.in.p], tspan=(16 * τ, 17 * τ), xlabel = "Time [s]", ylabel = "Pressure [mmHg]",  hidexaxis = nothing) # Make a line plot
p2 = plot(sol, idxs=[LV.V], tspan=(16 * τ, 17 * τ),xlabel = "Time [s]", ylabel = "Volume [ml]",  linkaxes = :all)
p3 = plot(sol, idxs=[Csa.in.q,Csv.in.q], tspan=(16 * τ, 17 * τ),xlabel = "Time [s]", ylabel = "Flow rate [ml/s]", linkaxes = :all)
p4 = plot(sol, idxs=(LV.V, LV.p), tspan=(16 * τ, 17 * τ),xlabel = "Volume [ml]", ylabel = "Pressure [mmHg]", linkaxes = :all)

img = plot(p1, p2, p3, p4; layout=@layout([a b; c d]), legend = true)

img = DisplayAs.Text(DisplayAs.PNG(img))

img
