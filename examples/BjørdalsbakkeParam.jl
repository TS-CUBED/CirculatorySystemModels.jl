# Heart rate and cycle time
HR = 70.58823529411765
τ = 60.0/HR
### Double Hill parameters as in Nikolai paper 
# Left Ventricle
Eₘᵢₙ = 0.03
Eₘₐₓ = 1.5
n1LV    = 1.32;
n2LV    = 21.9;
Tau1fLV = 0.303 * τ;
Tau2fLV = 0.508 * τ
## Model Parameters
Rs = 1.11
Csa = 1.13
Csv = 11.0

### Valve parameters ###
# Aortic valve basic 
Zao = 0.033
#Mitral valve basic 
Rmv = 0.006

#Inital Pressure
MCFP = 7.0