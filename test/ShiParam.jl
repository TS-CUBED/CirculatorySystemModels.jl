τ = 1.0
Eshift=0.0
Ev=Inf
#### LV chamber parameters #### Checked
v0_lv = 5.0
p0_lv = 1.0
Emin_lv = 0.1
Emax_lv = 2.5
τes_lv = 0.3
τed_lv = 0.45
Eshift_lv = 0.0
#### RV Chamber parameters #### Checked
v0_rv = 10.0
p0_rv = 1.0
Emin_rv = 0.1
Emax_rv = 1.15
τes_rv = 0.3
τed_rv = 0.45
Eshift_rv = 0.0
### LA Atrium Parameters #### Checked
v0_la = 4.0
p0_la = 1.0
Emin_la = 0.15
Emax_la = 0.25
τpwb_la = 0.92
τpww_la = 0.09
τes_la = τpww_la/2
τed_la = τpww_la
Eshift_la = τpwb_la
### RA Atrium parameters #### Checked
v0_ra = 4.0
p0_ra = 1.0
Emin_ra = 0.15
Emax_ra = 0.25
τpwb_ra = 0.92
τpww_ra = 0.09
τes_ra = τpww_ra/2
τed_ra = τpww_ra
Eshift_ra = τpwb_ra
#### Valve parameters #### Checked
CQ_AV = 350.0
CQ_MV = 400.0
CQ_TV = 400.0
CQ_PV = 350.0
## Systemic Aortic Sinus #### Checked
Csas = 0.08
Rsas = 0.003
Lsas = 6.2e-5
pt0sas = 100.0
qt0sas = 0.0
## Systemic Artery #### Checked
Csat = 1.6
Rsat = 0.05
Lsat = 0.0017
pt0sat = 100.0
qt0sat = 0.0
## Systemic Arteriole #### Checked
Rsar = 0.5
## Systemic Capillary #### Checked 
Rscp = 0.52
## Systemic Vein #### Checked
Csvn = 20.5
Rsvn = 0.075
pt0svn = 0.0
qt0svn = 0.0
## Pulmonary Aortic Sinus #### Checked
Cpas = 0.18
Rpas = 0.002
Lpas = 5.2e-5
pt0pas = 30.0
qt0pas = 0.0
## Pulmonary Artery #### Checked
Cpat = 3.8
Rpat = 0.01
Lpat = 0.0017
pt0pat = 30.0
qt0pat = 0.0
## Pulmonary Arteriole #### Checked
Rpar = 0.05
## Pulmonary Capillary #### Checked
Rpcp = 0.25
## Pulmonary Vein #### Checked
Cpvn = 20.5
Rpvn = 0.006           # this was 0.006 originally and in the paper, seems to be wrong in the paper!
# CHANGED THIS IN THE CELLML MODEL AS WELL TO MATCH THE PAPER!!!!!
pt0pvn = 0.0
qt0pvn = 0.0
## KG diaphragm ## Not in cellML model
# left heart #
Kst_la = 2.5
Kst_lv = 20.0
Kf_sav = 0.0004
Ke_sav = 9000.0
M_sav = 0.0004
A_sav = 0.00047
# right heart # 
Kst_ra = 2.5
Kst_rv = 20.0
Kf_pav = 0.0004
Ke_pav = 9000.0
M_pav = 0.0004
A_pav = 0.00047
#
#### Diff valve params #### not in cellML model
Kp_av = 5500.0 # *  57.29578 # Shi Paper has values in radians!
Kf_av = 50.0
Kf_mv = 50.0
Kp_mv = 5500.0 # *  57.29578 
Kf_tv = 50.0
Kp_tv = 5500.0 # *  57.29578
Kf_pv = 50.0
Kp_pv = 5500.0 #*  57.29578
Kb_av = 2.0
Kv_av = 7.0
Kb_mv = 2.0
Kv_mv = 3.5
Kb_tv = 2.0
Kv_tv = 3.5
Kb_pv = 2.0
Kv_pv = 3.5
θmax_av = 75.0 * pi / 180
θmax_mv = 75.0 * pi / 180
θmin_av = 5.0 * pi / 180
θmin_mv = 5.0 * pi / 180
θmax_pv = 75.0 * pi / 180
θmax_tv = 75.0 * pi / 180
θmin_pv = 5.0 * pi / 180
θmin_tv = 5.0 * pi / 180

## pressure force and frictional force is the same for all 4 valves 

# Initial conditions #### Checked against cellML model

LV_Vt0 = 500
RV_Vt0 = 500
LA_Vt0 = 20
RA_Vt0 = 20

