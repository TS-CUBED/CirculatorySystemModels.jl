using CirculatorySystemModels
using DifferentialEquations
using ModelingToolkit
using Test
using CSV
using DataFrames

# @testset "WK5" begin

# end

@testset "Shi Model" begin
    ##
    include("ShiParam.jl")

    ## Start Modelling
    @variables t

    ### Ventricles
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
    @time sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9, saveat=19:0.01:20)
    ShiSol = sol

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiSimple.csv", DataFrame)

    @test ShiSol.retcode == :Success
    @test sum(ShiSol[LV.V] .- ShiBench[!, :LV_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[RV.V] .- ShiBench[!, :RV_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[LA.V] .- ShiBench[!, :LA_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[RA.V] .- ShiBench[!, :RA_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
end

@testset "Shi Model Complex" begin
    ##
    include("ShiParam.jl")

    ## Start Modelling
    @variables t

    ### Shi Heart
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
    u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0, 0, 0, 0, 0]

    prob = ODAEProblem(circ_sys, u0, (0.0, 20.0))
    ##
    @time sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9, saveat=19:0.01:20)
    ShiSol = sol
    ##

    ## Read benchmark data and compare
    ShiBench = CSV.read("ShiComplex.csv", DataFrame)

    @test ShiSol.retcode == :Success
    @test sum(ShiSol[Heart.LV.V] .- ShiBench[!, :LV_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[Heart.RV.V] .- ShiBench[!, :RV_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[Heart.LA.V] .- ShiBench[!, :LA_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
    @test sum(ShiSol[Heart.RA.V] .- ShiBench[!, :RA_V]) / length(ShiSol.u) ≈ 0 atol = 1e-3
end




# Helpers to set up the tests
# df = DataFrame(
#     LV_V=ShiSol[LV.V],
#     RV_V=ShiSol[RV.V],
#     LA_V=ShiSol[LA.V],
#     RA_V=ShiSol[RA.V],
#     SAS_S_p=ShiSol[SAS.C.p],
#     SAS_L_q=ShiSol[SAS.L.q],
#     SAT_S_p=ShiSol[SAT.C.p],
#     SAT_L_q=ShiSol[SAT.L.q],
#     SVN_C_p=ShiSol[SVN.C.p],
#     PAS_S_p=ShiSol[PAS.C.p],
#     PAS_L_q=ShiSol[PAS.L.q],
#     PAT_S_p=ShiSol[PAT.C.p],
#     PAT_L_q=ShiSol[PAT.L.q],
#     PVN_C_p=ShiSol[PVN.C.p]
# )

# CSV.write("ShiComplex.df", df)

# df2 = CSV.read("ShiSimple.df", DataFrame)
