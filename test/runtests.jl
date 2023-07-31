##
using CirculatorySystemModels
using DifferentialEquations
using ModelingToolkit
using Test
using CSV
using DataFrames

##
# @testset "WK5" begin

# end

@testset "Shi Model" begin
    ##
    include("ShiParam.jl")

    ## Start Modelling
    @variables t

    ## Ventricles
    @named LV = ShiChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    @named LA = ShiChamber(V₀=v0_la, p₀=p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la / 2, τₑₚ=τpww_la, Eshift=τpwb_la)
    @named RV = ShiChamber(V₀=v0_rv, p₀=p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    # @named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)
    @named RA = ShiAtrium(V₀=v0_ra, p₀=1, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τpwb=τpwb_ra, τpww=τpww_ra) #, Ev=Inf)

    ## Valves as simple valves
    @named AV = OrificeValve(CQ=CQ_AV)
    @named MV = OrificeValve(CQ=CQ_MV)
    @named TV = OrificeValve(CQ=CQ_TV)
    @named PV = OrificeValve(CQ=CQ_PV)

    ####### Systemic Loop #######
    # Systemic Aortic Sinus ##
    @named SAS = CRL(C=Csas, R=Rsas, L=Lsas)
    # Systemic Artery ##
    @named SAT = CRL(C=Csat, R=Rsat, L=Lsat)
    # Systemic Arteriole ##
    @named SAR = Resistor(R=Rsar)
    # Systemic Capillary ##
    @named SCP = Resistor(R=Rscp)
    # Systemic Vein ##
    @named SVN = CR(R=Rsvn, C=Csvn)

    ####### Pulmonary Loop #######
    # Pulmonary Aortic Sinus ##
    @named PAS = CRL(C=Cpas, R=Rpas, L=Lpas)
    # Pulmonary Artery ##
    @named PAT = CRL(C=Cpat, R=Rpat, L=Lpat)
    # Pulmonary Arteriole ##
    @named PAR = Resistor(R=Rpar)
    # Pulmonary Capillary ##
    @named PCP = Resistor(R=Rpcp)
    # Pulmonary Vein ##
    @named PVN = CR(R=Rpvn, C=Cpvn)

    ##
    circ_eqs = [
        connect(LV.out, AV.in)
        connect(AV.out, SAS.in)
        connect(SAS.out, SAT.in)
        connect(SAT.out, SAR.in)
        connect(SAR.out, SCP.in)
        connect(SCP.out, SVN.in)
        connect(SVN.out, RA.in)
        connect(RA.out, TV.in)
        connect(TV.out, RV.in)
        connect(RV.out, PV.in)
        connect(PV.out, PAS.in)
        connect(PAS.out, PAT.in)
        connect(PAT.out, PAR.in)
        connect(PAR.out, PCP.in)
        connect(PCP.out, PVN.in)
        connect(PVN.out, LA.in)
        connect(LA.out, MV.in)
        connect(MV.out, LV.in)
    ]

    ## Compose the whole ODAE system
    @named _circ_model = ODESystem(circ_eqs, t)
    @named circ_model = compose(_circ_model,
        [LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

    ## And simplify it
    circ_sys = structural_simplify(circ_model)

    ## Setup ODE
    # Initial Conditions for Shi Valve
    #u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0,0, 0, 0, 0]
    # and for OrificeValve --- Commment this next line to use ShiValves
    u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

    prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-9, saveat=19:0.01:20)
    ShiSimpleSol = sol(19:0.01:20)

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiSimple.csv", DataFrame)

    @test SciMLBase.successful_retcode(sol)
    @test sum((ShiSimpleSol[LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiSimpleSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSol[RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiSimpleSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSol[LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiSimpleSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSol[RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiSimpleSol.u) ≈ 0 atol = 1e-3
    ##
end

##

@testset "Shi Model Complex" begin

##
    ##
    include("ShiParam.jl")
    ##
    ## Start Modelling
    @variables t

    ### Shi Heart (with AV stenosis: max AV opening angle = 40 degrees!)
    @named Heart = ShiHeart(τ=τ,
        LV_V₀=v0_lv, LV_p0=p0_lv, LV_Emin=Emin_lv, LV_Emax=Emax_lv, LV_τes=τes_lv, LV_τed=τed_lv, LV_Eshift=0.0,
        RV_V₀=v0_rv, RV_p0=p0_rv, RV_Emin=Emin_rv, RV_Emax=Emax_rv, RV_τes=τes_rv, RV_τed=τed_rv, RV_Eshift=0.0,
        LA_V₀=v0_la, LA_p0=p0_la, LA_Emin=Emin_la, LA_Emax=Emax_la, LA_τes=τpww_la / 2, LA_τed=τpww_la, LA_Eshift=τpwb_la,
        RA_V₀=v0_ra, RA_p0=p0_ra, RA_Emin=Emin_ra, RA_Emax=Emax_ra, RA_τes=τpww_ra / 2, RA_τed=τpww_ra, RA_Eshift=τpwb_ra,
        AV_CQ=CQ_AV, AV_Kp=Kp_av, AV_Kf=Kf_av, AV_Kb=0.0, AV_Kv=3.5, AV_θmax=40.0 * pi / 180, AV_θmin=5.0 * pi / 180,
        MV_CQ=CQ_MV, MV_Kp=Kp_mv, MV_Kf=Kf_mv, MV_Kb=0.0, MV_Kv=3.5, MV_θmax=75.0 * pi / 180, MV_θmin=5.0 * pi / 180,
        TV_CQ=CQ_TV, TV_Kp=Kp_tv, TV_Kf=Kf_tv, TV_Kb=0.0, TV_Kv=3.5, TV_θmax=75.0 * pi / 180, TV_θmin=5.0 * pi / 180,
        PV_CQ=CQ_PV, PV_Kp=Kp_pv, PV_Kf=Kf_pv, PV_Kb=0.0, PV_Kv=3.5, PV_θmax=75.0 * pi / 180, PV_θmin=5.0 * pi / 180,
    )
    ### Ventricles
    # @named LV = ShiChamber(V₀=v0_lv, p₀ = p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
    # # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    # @named LA = ShiChamber(V₀=v0_la, p₀ = p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la/2, τₑₚ=τpww_la, Eshift=τpwb_la)
    # @named RV = ShiChamber(V₀=v0_rv, p₀ = p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
    # # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    # # @named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)
    # @named RA = ShiAtrium(V₀=v0_ra, p₀ = 1, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τpwb=τpwb_ra, τpww=τpww_ra) #, Ev=Inf)

    # ## 4 Valves
    # @named AV = ShiValve(CQ=CQ_AV, Kp=Kp_av, Kf=Kf_av, Kb= 0.0, Kv=3.5, θmax=40.0 * pi / 180, θmin=5.0 * pi / 180)
    # @named MV = ShiValve(CQ=CQ_MV, Kp=Kp_mv, Kf=Kf_mv, Kb= 0.0, Kv=3.5, θmax=75.0 * pi / 180, θmin=5.0 * pi / 180)
    # @named TV = ShiValve(CQ=CQ_TV, Kp=Kp_tv, Kf=Kf_tv, Kb= 0.0, Kv=3.5, θmax=75.0 * pi / 180, θmin=5.0 * pi / 180)
    # @named PV = ShiValve(CQ=CQ_PV, Kp=Kp_pv, Kf=Kf_pv, Kb= 0.0, Kv=3.5, θmax=75.0 * pi / 180, θmin=5.0 * pi / 180)

    @named SystLoop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas,
        SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat,
        SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn)
    # ####### Systemic Loop #######
    # # Systemic Aortic Sinus ##
    # @named SAS = CRL(C=Csas, R=Rsas, L=Lsas)
    # # Systemic Artery ##
    # @named SAT = CRL(C=Csat, R=Rsat, L=Lsat)
    # # Systemic Arteriole ##
    # @named SAR = Resistor(R=Rsar)
    # # Systemic Capillary ##
    # @named SCP = Resistor(R=Rscp)
    # # Systemic Vein ##
    # @named SVN = CR(R=Rsvn, C=Csvn)

    @named PulmLoop = ShiPulmonaryLoop(PAS_C=Cpas, PAS_R=Rpas, PAS_L=Lpas,
        PAT_C=Cpat, PAT_R=Rpat, PAT_L=Lpat,
        PAR_R=Rpar, PCP_R=Rpcp, PVN_C=Cpvn, PVN_R=Rpvn)
    # ####### Pulmonary Loop #######
    # # Pulmonary Aortic Sinus ##
    # @named PAS = CRL(C=Cpas, R=Rpas, L=Lpas)
    # # Pulmonary Artery ##
    # @named PAT = CRL(C=Cpat, R=Rpat, L=Lpat)
    # # Pulmonary Arteriole ##
    # @named PAR = Resistor(R=Rpar)
    # # Pulmonary Capillary ##
    # @named PCP = Resistor(R=Rpcp)
    # # Pulmonary Vein ##
    # @named PVN = CR(R=Rpvn, C=Cpvn)

    ##
    circ_eqs = [
        connect(Heart.LHout, SystLoop.in)
        connect(SystLoop.out, Heart.RHin)
        connect(Heart.RHout, PulmLoop.in)
        connect(PulmLoop.out, Heart.LHin)
    ]
    # circ_eqs = [
    #     connect(Heart.LHout, SAS.in)
    #     connect(SAS.out, SAT.in)
    #     connect(SAT.out, SAR.in)
    #     connect(SAR.out, SCP.in)
    #     connect(SCP.out, SVN.in)
    #     connect(SVN.out, Heart.RHin)
    #     connect(Heart.RHout, PAS.in)
    #     connect(PAS.out, PAT.in)
    #     connect(PAT.out, PAR.in)
    #     connect(PAR.out, PCP.in)
    #     connect(PCP.out, PVN.in)
    #     connect(PVN.out, Heart.LHin)
    # ]

    ## Compose the whole ODAE system
    @named _circ_model = ODESystem(circ_eqs, t)
    @named circ_model = compose(_circ_model,
        [Heart, SystLoop, PulmLoop])
    # [Heart, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

    ## And simplify it
    circ_sys = structural_simplify(circ_model)

    ## Setup ODE
    # Initial Conditions for Shi Valve
    u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, 0, 0, 0, 0, 0, 0, 0, 0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

    prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-9, saveat=19:0.01:20)
    ShiComplexSol = sol(19:0.01:20)
    ##

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiComplex.csv", DataFrame)

    @test SciMLBase.successful_retcode(sol)
    @test sum((ShiComplexSol[Heart.LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiComplexSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSol[Heart.RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiComplexSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSol[Heart.LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiComplexSol.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSol[Heart.RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiComplexSol.u) ≈ 0 atol = 1e-3
    ##

##
end


##
@testset "Bjørdalsbakke" begin

##
    using ModelingToolkit, DifferentialEquations
    using CirculatorySystemModels

    # # A simple single-chamber model
    #
    # ![Single chamber, closed-loop, lumped parameter model of the systemic circulation and the left ventricle. The circuit equivalent formulation of the model is depicted, with the pressures of each compartment, as well as most of the mechanical parameters. The model describes three compartments: the left ventricular, arterial and venous compartments. 𝑃𝑡ℎ is the intrathoracic pressure, 𝑃𝑙𝑣 is the left ventricular pressure and 𝐸𝑙𝑣(𝑡) indicates the left ventricular elastance function.](./BjordalsbakkeModelSk
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
    n1LV = 1.32
    n2LV = 21.9
    Tau1fLV = 0.303 * τ
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
    # $k$ is a scaling factor to assure that $e(t)$ has a maximum of $e(t)_{max} = 1$:
    #
    # $$k = \max \left(\frac{\left(\tau / \tau_1\right)^{n_1}}{1+\left(\tau / \tau_1\right)^{n_1}} \times \frac{1}{1+\left(\tau / \tau_2\right)^{n_2}} \right)^{-1}$$ .
    #

    nstep = 1000
    t = LinRange(0, τ, nstep)

    kLV = 1 / maximum((t ./ Tau1fLV) .^ n1LV ./ (1 .+ (t ./ Tau1fLV) .^ n1LV) .* 1 ./ (1 .+ (t ./ Tau2fLV) .^ n2LV))


    # ## Set up the model elements
    #
    # Set up time as a parameter `t`
    #
    @parameters t

    # Heart is modelled as a single chamber (we call it `LV` for "Left Ventricle" so the model can be extended later, if required):
    #
    @named LV = DHChamber(V₀=0.0, Eₘₐₓ=Eₘₐₓ, Eₘᵢₙ=Eₘᵢₙ, n₁=n1LV, n₂=n2LV, τ=τ, τ₁=Tau1fLV, τ₂=Tau2fLV, k=kLV, Eshift=0.0, Ev=Inf)

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

    # `circ_sys` is now the minimal system of equations. In this case it consists of 3 ODEs for the three pressures.
    #
    # _Note: `structural_simplify` reduces and optimises the ODE system. It is, therefore, not always obvious, which states it will use and which it will drop. We can use the `states` and `observed` function to check this. It is recommended to do this, since small changes can reorder states, observables, and parameters._
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
    u0 = [MCFP, MCFP, MCFP]

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
    sol = solve(prob, Vern7(), reltol=1e-6, abstol=1e-9, saveat=(19:0.01:20))
    BBsol = sol(19:0.01:20)

    BBbench = CSV.read("BB.csv", DataFrame)

    @test SciMLBase.successful_retcode(sol)
    @test sum((BBsol[LV.V] .- BBbench[!, :LV_V]) ./ BBbench[!, :LV_V]) / length(BBsol.u) ≈ 0 atol = 1e-3
    @test sum((BBsol[Csa.p] .- BBbench[!, :Csa_p]) ./ BBbench[!, :Csa_p]) / length(BBsol.u) ≈ 0 atol = 1e-3
    @test sum((BBsol[Csv.p] .- BBbench[!, :Csv_p]) ./ BBbench[!, :Csv_p]) / length(BBsol.u) ≈ 0 atol = 1e-3

##
end




# Helpers to set up the tests
# df = DataFrame(
#     t=ShiSol.t,
#     LV_V=ShiSol[Heart.LV.V],
#     RV_V=ShiSol[Heart.RV.V],
#     LA_V=ShiSol[Heart.LA.V],
#     RA_V=ShiSol[Heart.RA.V],
#     SAS_S_p=ShiSol[SystLoop.SAS.C.p],
#     SAS_L_q=ShiSol[SystLoop.SAS.L.q],
#     SAT_S_p=ShiSol[SystLoop.SAT.C.p],
#     SAT_L_q=ShiSol[SystLoop.SAT.L.q],
#     SVN_C_p=ShiSol[SystLoop.SVN.C.p],
#     PAS_S_p=ShiSol[PulmLoop.PAS.C.p],
#     PAS_L_q=ShiSol[PulmLoop.PAS.L.q],
#     PAT_S_p=ShiSol[PulmLoop.PAT.C.p],
#     PAT_L_q=ShiSol[PulmLoop.PAT.L.q],
#     PVN_C_p=ShiSol[PulmLoop.PVN.C.p]
# )

# CSV.write("ShiComplex.csv", df)

##
