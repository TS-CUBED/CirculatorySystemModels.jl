##
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner
using Test
using CSV
using DataFrames
##
# @testset "WK5" begin

# end

@testset "Shi Model in V" begin
    ##
    include("ShiParam.jl")

    @independent_variables t

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

    ## Compose the whole ODE system
    @named _circ_model = System(circ_eqs, t)
    @named circ_model = compose(_circ_model,
        [LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

    ## And simplify it
    circ_sys = mtkcompile(circ_model)

    ## Setup ODE
    # Initial Conditions for Shi Valve
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0,0, 0, 0, 0]
    # and for OrificeValve --- Commment this next line to use ShiValves
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

    u0 = [
        LV.V => LV_Vt0
        RV.V => RV_Vt0
        LA.V => LA_Vt0
        RA.V => RA_Vt0
        SAS.C.p => pt0sas
        SAS.L.q => qt0sas
        SAT.C.p => pt0sat
        SAT.L.q => qt0sat
        SVN.C.p => pt0svn
        PAS.C.p => pt0pas
        PAS.L.q => qt0pas
        PAT.C.p => pt0pat
        PAT.L.q => qt0pat
        PVN.C.p => pt0pvn
    ]

    prob = ODEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time ShiSimpleSolV = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)
    # ShiSimpleSolV = ShiSimpleSolV(19:0.01:20)

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiSimple.csv", DataFrame)

    @test SciMLBase.successful_retcode(ShiSimpleSolV)
    @test sum((ShiSimpleSolV[LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
    @test sum((ShiSimpleSolV[RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiSimpleSolV.u) ≈ 0 atol = 1e-3
end

##
@testset "Shi Model in P" begin
    include("ShiParam.jl")

    ## Start Modelling
    @independent_variables t

    ## Ventricles
    @named LV = ShiChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0, inP=true)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    @named LA = ShiChamber(V₀=v0_la, p₀=p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la / 2, τₑₚ=τpww_la, Eshift=τpwb_la, inP=true)
    @named RV = ShiChamber(V₀=v0_rv, p₀=p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0, inP=true)
    # The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
    # @named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)
    @named RA = ShiAtrium(V₀=v0_ra, p₀=1, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τpwb=τpwb_ra, τpww=τpww_ra, inP=true) #, Ev=Inf)

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

    ## Compose the whole ODE system
    @named _circ_model = System(circ_eqs, t)
    @named circ_model = compose(_circ_model,
        [LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

    ## And simplify it
    circ_sys = mtkcompile(circ_model)

    ## Setup ODE
    # Initial Conditions for Shi Valve
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0,0, 0, 0, 0]
    # and for OrificeValve --- Commment this next line to use ShiValves
    # u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

    u0 = [
        LV.V => LV_Vt0
        RV.V => RV_Vt0
        LA.V => LA_Vt0
        RA.V => RA_Vt0
        SAS.C.p => pt0sas
        SAS.L.q => qt0sas
        SAT.C.p => pt0sat
        SAT.L.q => qt0sat
        SVN.C.p => pt0svn
        PAS.C.p => pt0pas
        PAS.L.q => qt0pas
        PAT.C.p => pt0pat
        PAT.L.q => qt0pat
        PVN.C.p => pt0pvn
    ]

    prob = ODEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time ShiSimpleSolP = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)
    # ShiSimpleSolP = ShiSimpleSolP(19:0.01:20)

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiSimple.csv", DataFrame)

    @test SciMLBase.successful_retcode(ShiSimpleSolP)
    @test sum((ShiSimpleSolP[LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiSimpleSolP.u) ≈ 0 atol = 2e-3
    @test sum((ShiSimpleSolP[RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiSimpleSolP.u) ≈ 0 atol = 2e-3
    @test sum((ShiSimpleSolP[LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiSimpleSolP.u) ≈ 0 atol = 2e-3
    @test sum((ShiSimpleSolP[RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiSimpleSolP.u) ≈ 0 atol = 2e-3
end

##

@testset "Shi Model Complex" begin
    include("ShiParam.jl")

    ## Start Modelling
    @independent_variables t

    ## Shi Heart (with AV stenosis: max AV opening angle = 40 degrees!)
    @named heart = ShiHeart(τ=τ,
        LV_V₀=v0_lv, LV_p₀=p0_lv, LV_Eₘᵢₙ=Emin_lv, LV_Eₘₐₓ=Emax_lv, LV_τ=τ, LV_τₑₛ=τes_lv, LV_τₑₚ=τed_lv, LV_Eshift=0.0,
        RV_V₀=v0_rv, RV_p₀=p0_rv, RV_Eₘᵢₙ=Emin_rv, RV_Eₘₐₓ=Emax_rv, RV_τ=τ, RV_τₑₛ=τes_rv, RV_τₑₚ=τed_rv, RV_Eshift=0.0,
        LA_V₀=v0_la, LA_p₀=p0_la, LA_Eₘᵢₙ=Emin_la, LA_Eₘₐₓ=Emax_la, LA_τ=τ, LA_τₑₛ=τpww_la / 2, LA_τₑₚ=τpww_la, LA_Eshift=τpwb_la,
        RA_V₀=v0_ra, RA_p₀=p0_ra, RA_Eₘᵢₙ=Emin_ra, RA_Eₘₐₓ=Emax_ra, RA_τ=τ, RA_τₑₛ=τpww_ra / 2, RA_τₑₚ=τpww_ra, RA_Eshift=τpwb_ra,
        AV_CQ=CQ_AV, AV_Kp=Kp_av, AV_Kf=Kf_av, AV_Kb=0.0, AV_Kv=3.5, AV_θmax=40.0 * pi / 180, AV_θmin=5.0 * pi / 180,
        MV_CQ=CQ_MV, MV_Kp=Kp_mv, MV_Kf=Kf_mv, MV_Kb=0.0, MV_Kv=3.5, MV_θmax=75.0 * pi / 180, MV_θmin=5.0 * pi / 180,
        TV_CQ=CQ_TV, TV_Kp=Kp_tv, TV_Kf=Kf_tv, TV_Kb=0.0, TV_Kv=3.5, TV_θmax=75.0 * pi / 180, TV_θmin=5.0 * pi / 180,
        PV_CQ=CQ_PV, PV_Kp=Kp_pv, PV_Kf=Kf_pv, PV_Kb=0.0, PV_Kv=3.5, PV_θmax=75.0 * pi / 180, PV_θmin=5.0 * pi / 180
    )
    @named syst_loop = ShiSystemicLoop(SAS_C=Csas, SAS_R=Rsas, SAS_L=Lsas,
        SAT_C=Csat, SAT_R=Rsat, SAT_L=Lsat,
        SAR_R=Rsar, SCP_R=Rscp, SVN_C=Csvn, SVN_R=Rsvn
    )
    @named pulm_loop = ShiPulmonaryLoop(PAS_C=Cpas, PAS_R=Rpas, PAS_L=Lpas,
        PAT_C=Cpat, PAT_R=Rpat, PAT_L=Lpat,
        PAR_R=Rpar, PCP_R=Rpcp, PVN_C=Cpvn, PVN_R=Rpvn
    )

    circ_eqs = [
        connect(heart.LHout, syst_loop.in)
        connect(syst_loop.out, heart.RHin)
        connect(heart.RHout, pulm_loop.in)
        connect(pulm_loop.out, heart.LHin)
    ]

    @named _circ_model = System(circ_eqs, t)
    @named circ_model = compose(_circ_model, [heart, syst_loop, pulm_loop])
    circ_sys = mtkcompile(circ_model)

    u0 = [
        heart.LV.V => LV_Vt0
        heart.RV.V => RV_Vt0
        heart.LA.V => LA_Vt0
        heart.RA.V => RA_Vt0
        heart.AV.θ => 5.0 * pi / 180 + 0.01
        heart.AV.ω => 0
        heart.MV.θ => 5.0 * pi / 180 + 0.01
        heart.MV.ω => 0
        heart.TV.θ => 5.0 * pi / 180 + 0.01
        heart.TV.ω => 0
        heart.PV.θ => 5.0 * pi / 180 + 0.01
        heart.PV.ω => 0
        syst_loop.SAS.C.p => pt0sas
        syst_loop.SAS.L.q => qt0sas
        syst_loop.SAT.C.p => pt0sat
        syst_loop.SAT.L.q => qt0sat
        syst_loop.SVN.C.p => pt0svn
        pulm_loop.PAS.C.p => pt0pas
        pulm_loop.PAS.L.q => qt0pas
        pulm_loop.PAT.C.p => pt0pat
        pulm_loop.PAT.L.q => qt0pat
        pulm_loop.PVN.C.p => pt0pvn
    ]

    prob = ODEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time ShiComplexSol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-9, saveat=19:0.01:20)
    # The callbacks prevent saveat from working as intended! So I need to interpolate the results:
    ShiComplexSolInt = ShiComplexSol(19:0.01:20)
    ##

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiComplex.csv", DataFrame)

    @test SciMLBase.successful_retcode(ShiComplexSol)
    @test sum((ShiComplexSolInt[heart.LV.V] .- ShiBench[!, :LV_V]) ./ ShiBench[!, :LV_V]) / length(ShiComplexSolInt.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSolInt[heart.RV.V] .- ShiBench[!, :RV_V]) ./ ShiBench[!, :RV_V]) / length(ShiComplexSolInt.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSolInt[heart.LA.V] .- ShiBench[!, :LA_V]) ./ ShiBench[!, :LA_V]) / length(ShiComplexSolInt.u) ≈ 0 atol = 1e-3
    @test sum((ShiComplexSolInt[heart.RA.V] .- ShiBench[!, :RA_V]) ./ ShiBench[!, :RA_V]) / length(ShiComplexSolInt.u) ≈ 0 atol = 1e-3
    ##
end


##
@testset "Bjørdalsbakke" begin

    ##
    using ModelingToolkit
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
    R_s = 1.11
    C_sa = 1.13
    C_sv = 11.0

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
    @independent_variables t

    # Heart is modelled as a single chamber (we call it `LV` for "Left Ventricle" so the model can be extended later, if required):
    #
    @named LV = DHChamber(V₀=0.0, Eₘₐₓ=Eₘₐₓ, Eₘᵢₙ=Eₘᵢₙ, n₁=n1LV, n₂=n2LV, τ=τ, τ₁=Tau1fLV, τ₂=Tau2fLV, k=kLV, Eshift=0.0, inP=true)

    # The two valves are simple diodes with a small resistance
    # (resistance is needed, since perfect diodes would connect two elastances/compliances, which will lead to unstable oscillations):
    #
    @named AV = ResistorDiode(R=Zao)
    @named MV = ResistorDiode(R=Rmv)

    # The main components of the circuit are 1 resistor `Rs` and two compliances for systemic arteries `Csa`,
    # and systemic veins `Csv` (names are arbitrary).
    #
    @named Rs = Resistor(R=R_s)

    @named Csa = Compliance(C=C_sa, inP=true, has_ep=true, has_variable_ep=true)
    # @named Csv = Compliance(C=C_sv, inV=true, has_ep=true)
    @named Csv = Elastance(E=1/C_sv, inP=false)

    # _Note: testing different parameters and formulations here. May need to do all of them separately to really do a unit test.
    #
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
        Csa.ep.p ~ 0
    ]

    # ### Add the component equations
    #
    # In a second step, the system of Kirchhoff equations is completed by the component equations (both ODEs and AEs), resulting in the full, overdefined ODE set `circ_model`.
    #
    # _Note: we do this in two steps._
    #
    @named _circ_model = System(circ_eqs, t)

    @named circ_model = compose(_circ_model,
        [LV, AV, MV, Rs, Csa, Csv])

    # ### Simplify the ODE system
    #
    # The crucial step in any acausal modelling is the sympification and reduction of the OD(A)E system to the minimal set of equations. ModelingToolkit.jl does this in the `mtkcompile` function.
    #
    circ_sys = mtkcompile(circ_model)

    # `circ_sys` is now the minimal system of equations. In this case it consists of 3 ODEs for the three pressures.

    #
    # _Note: `mtkcompile` reduces and optimises the ODE system. It is, therefore, not always obvious, which states it will use and which it will drop. We can use the `states` and `observed` function to check this. It is recommended to do this, since small changes can reorder states, observables, and parameters._
    #
    # States in the system are now:
    #unknowns(circ_sys)

    # Observed parameters - the system will drop these from the ODE system that is solved, but it keeps all the algebraic equations needed to calculate them in the system object, as well as the `ODEProblem` and solution object - are:
    #observed(circ_sys)

    # And the parameters (these could be reordered, so check these, too):
    #parameters(circ_sys)

    # ### Define the ODE problem
    #
    # First defined initial conditions `u0` and the time span for simulation:
    #
    u0 = [
        LV.p => MCFP
        Csa.p => MCFP
        Csv.p => MCFP
    ]

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
    @time BBsol = solve(prob, Vern7(), reltol=1e-6, abstol=1e-9, saveat=(19:0.01:20))
    # BBsol = sol(19:0.01:20)

    BBbench = CSV.read("BB.csv", DataFrame)

    @test SciMLBase.successful_retcode(BBsol)
    @test sum((BBsol[LV.V] .- BBbench[!, :LV_V]) ./ BBbench[!, :LV_V]) / length(BBsol.u) ≈ 0 atol = 1.5e-3
    @test sum((BBsol[Csa.p] .- BBbench[!, :Csa_p]) ./ BBbench[!, :Csa_p]) / length(BBsol.u) ≈ 0 atol = 1.5e-3
    @test sum((BBsol[Csv.p] .- BBbench[!, :Csv_p]) ./ BBbench[!, :Csv_p]) / length(BBsol.u) ≈ 0 atol = 1.5e-3

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
#     SAT_L_q=ShiSol[SystLoop.SAT.L.q] oj,
#     SVN_C_p=ShiSol[SystLoop.SVN.C.p],
#     PAS_S_p=ShiSol[PulmLoop.PAS.C.p],
#     PAS_L_q=ShiSol[PulmLoop.PAS.L.q],
#     PAT_S_p=ShiSol[PulmLoop.PAT.C.p],
#     PAT_L_q=ShiSol[PulmLoop.PAT.L.q],
#     PVN_C_p=ShiSol[PulmLoop.PVN.C.p]
# )

# CSV.write("ShiComplex.csv", df)

##
