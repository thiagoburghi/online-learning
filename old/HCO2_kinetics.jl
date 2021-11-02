# The HCO model was taken from the upcoming book by Drion, Franci and Sepulchre

const VNa = 45.     # 45.
const VCa = 120.
const VK = -90.     #-90.
const VH= -43.
const Vleak = -55.

## Model gating functions

# Synaptic current
V12syn=-20.0*δ[45,:];
ksyn=4.0*δ[46,:];
msyn_inf(V::Float64,i)=1.0/( 1.0 + exp(-(V-V12syn[i])/ksyn[i]) )

# Na-current (m=activation variable, h=inactivation variable)
vh_α_m = 40*δ[1,:]
vh_β_m = 65*δ[2,:]
k_α_m = 10*δ[3,:]
k_β_m = 18*δ[4,:] 
alpha_m(V::Float64,i) = -0.025*(V+vh_α_m[i])/(exp(-(V+vh_α_m[i])/k_α_m[i]) - 1.0 )
beta_m(V::Float64,i) = exp(-(V+vh_β_m[i])/k_β_m[i])
m_inf(V::Float64,i) = alpha_m(V,i) / (alpha_m(V,i) + beta_m(V,i))
tau_m(V::Float64,i) = 1 / (alpha_m(V,i) + beta_m(V,i))

vh_α_h = 65*δ[5,:]
vh_β_h = 35*δ[6,:]
k_α_h = 20*δ[7,:]
k_β_h = 10*δ[8,:] 
alpha_h(V::Float64,i) = 0.0175*exp(-(V+vh_α_h[i])/k_α_h[i])
beta_h(V::Float64,i) = 0.25/(1.0 + exp(-(V+vh_β_h[i])/k_β_h[i]) )
h_inf(V::Float64,i) = alpha_h(V,i) / (alpha_h(V,i) + beta_h(V,i))
tau_h(V::Float64,i) = 1 / (alpha_h(V,i) + beta_h(V,i))

# KD-current (mKD=activation variable)
vh_α_mKD = 55*δ[9,:]
vh_β_mKD = 65*δ[10,:]
k_α_mKD = 10*δ[11,:]
k_β_mKD = 80*δ[12,:] 
KDshift=10.0
alpha_mKD(V::Float64,i) = 0.0025*(V+vh_α_mKD[i])/(1. - exp(-(V+vh_α_mKD[i])/k_α_mKD[i]) )
beta_mKD(V::Float64,i) = 0.03125*exp(-(V+vh_β_mKD[i])/k_β_mKD[i])
mKD_inf(V::Float64,i) = alpha_mKD(V-KDshift,i) / (alpha_mKD(V-KDshift,i) + beta_mKD(V-KDshift,i))
tau_mKD(V::Float64,i) = 1 / (alpha_mKD(V-KDshift,i) + beta_mKD(V-KDshift,i))

# H-current (mH=activation variable)
w_α_mH = 14.59*δ[13,:]
w_β_mH = 1.87*δ[14,:]
b_α_mH = 0.086*δ[15,:]
b_β_mH = 0.0701*δ[16,:]
alpha_mH(V::Float64,i)= exp(-w_α_mH[i]-(b_α_mH[i]*V))
beta_mH(V::Float64,i)= exp(-w_β_mH[i]+(b_β_mH[i]*V))
mH_inf(V::Float64,i)= alpha_mH(V,i) /(alpha_mH(V,i) + beta_mH(V,i))
tau_mH(V_taumH,i) = 1/(alpha_mH(V_taumH,i) + beta_mH(V_taumH,i))
# dmH_inf(V::Float64)=((((0 - (0.086 * 1)) * exp(-14.59 - 0.086*V)) * (exp(-14.59 - 0.086*V) + exp(-1.87 + 0.0701*V)) - exp(-14.59 - 0.086*V) * ((0 - (0.086 * 1)) * exp(-14.59 - 0.086*V) + (0.0701 * 1) * exp(-1.87 + 0.0701*V))) / (exp(-14.59 - 0.086*V) + exp(-1.87 + 0.0701*V)) ^ 2)

# A-current (mA=activation variable, hA=inactivation variable)
vh_∞_mA = 60*δ[17,:]
vh_τ_mA1 = 35.82*δ[18,:]
vh_τ_mA2 = 79.69*δ[19,:]
k_∞_mA = 8.5*δ[20,:]
k_τ_mA1 = 19.697*δ[21,:]
k_τ_mA2 = -12.7*δ[22,:]
mA_inf(V::Float64,i) = 1/(1+exp(-(V+vh_∞_mA[i])/k_∞_mA[i]))
tau_mA_temp(V::Float64,i) = 0.37 + 1/(exp((V+vh_τ_mA1[i])/k_τ_mA1[i])+exp((V+vh_τ_mA2[i])/k_τ_mA2[i]))
tau_mA(V::Float64,i) = tau_mA_temp(V,i)

vh_∞_hA = 78*δ[23,:]
vh_τ_hA1 = 46.05*δ[24,:]
vh_τ_hA2 = 238.4*δ[25,:]
k_∞_hA = 6*δ[26,:]
k_τ_hA1 = 5*δ[27,:]
k_τ_hA2 = -37.45*δ[28,:]
hA_inf_temp(V::Float64,i) = 1/(1+exp((V+vh_∞_hA[i])/k_∞_hA[i]))
hA_inf(V,i) = hA_inf_temp(V,i)
function tau_hA(V::Float64,i)
    if V < -63
        tau_hA = 1/(exp((V+vh_τ_hA1[i])/k_τ_hA1[i])+exp((V+vh_τ_hA2[i])/k_τ_hA2[i]))
    else
        tau_hA = 19
    end
    return tau_hA
end
#tau_hA(V::Float64)=50.

# T-type Ca-current (mt=activation variable, ht=inactivation variable)
vh_∞_mt = 57*δ[29,:]; 
vh_τ_mt1 = 131.6*δ[30,:];
vh_τ_mt2 = 16.8*δ[31,:];
k_∞_mt = 6.2*δ[32,:];
k_τ_mt1 = 16.7*δ[33,:]; 
k_τ_mt2 = 18.2*δ[34,:]; 
mt_inf(V::Float64,i) = 1/(1+exp(-(V+vh_∞_mt[i])/k_∞_mt[i]))
tau_mt(V::Float64,i) = 0.612 + 1/(exp(-(V+vh_τ_mt1[i])/k_τ_mt1[i])+exp((V+vh_τ_mt2[i])/k_τ_mt2[i]))*2

vh_∞_ht = 81*δ[35,:]; 
vh_τ_ht1 = 467*δ[36,:];
vh_τ_ht2 = 21.88*δ[37,:];
k_∞_ht = 4.03*δ[38,:];
k_τ_ht1 = 66.6*δ[39,:]; 
k_τ_ht2 = 10.2*δ[40,:];
ht_inf(V::Float64,i) = 1/(1+exp((V+vh_∞_ht[i])/k_∞_ht[i]))
function tau_ht(V::Float64,i)
    if V < -80
        tau_ht = exp((V+vh_τ_ht1[i])/k_τ_ht1[i])*2
    else
        tau_ht = (exp(-(V+vh_τ_ht2[i])/k_τ_ht2[i])+28)*2
    end
    return tau_ht
end

# L-type Ca-current (mL=activation variable) (from Drion2011)
vh_∞_mL = 55*δ[41,:]; 
vh_τ_mL = 45*δ[42,:];
k_∞_mL = 3*δ[43,:];
k_τ_mL = 400*δ[44,:]; 
mL_inf(V::Float64,i) = 1/(1+exp(-(V+vh_∞_mL[i])/k_∞_mL[i]))
tau_mL(V::Float64,i) = (72*exp(-(V+vh_τ_mL[i])^2/k_τ_mL[i])+6.)*2

# Intracellular calcium
ICa_pump(Ca::Float64)=0.1*Ca/(Ca+0.0001)

## Definition of stimulation function
heaviside(t)=(1+sign(t))/2
pulse(t,ti,tf)=heaviside(t-ti)-heaviside(t-tf)

## Model ordinary differential equations for current-clamp simulations
function HM_ODE(du,u,p,t)
    # Parameters
    Iapp=p[1]
    I1=p[2]
    I2=p[3]
    ti1=p[4]
    tf1=p[5]
    ti2=p[6]
    tf2=p[7]
    gT=p[8]
    gKD=p[9]
    gH=p[10]
    gNa=p[11]
    gA=p[12]
    gLeak=p[13]
    gL=p[14]
    gKCa=p[15]
    C=p[16]
    taunoise=p[18]
    Ain=p[19]
    Win=p[20]

    # Variables
    V=u[1]
    m=u[2]
    h=u[3]
    mH=u[4]
    mt=u[5]
    ht=u[6]
    mA=u[7]
    hA=u[8]
    mKD=u[9]
    mL=u[10]
    Ca=u[11]
    noise=u[12]

    # ODEs
    du[1]=1/C*(- gNa*m^3*h*(V-VNa) - gH*mH*(V-VH) - gT*mt^2*ht*(V-VCa) - gA*mA^4*hA*(V-VK) - gKD*mKD^4*(V-VK) - gLeak*(V-Vleak)- gL*mL*(V-VCa) - gKCa*(Ca/(15.0+Ca))^4*(V-VK) + Iapp +I1*pulse(t,ti1,tf1) +I2*pulse(t,ti2,tf2) + noise + Ain*sin(2*pi*Win*t))
    du[2]=1/tau_m(V)*(-m+m_inf(V))
    du[3]=1/tau_h(V)*(-h+h_inf(V))
    du[4]=(-mH+mH_inf(V))*1/tau_mH(V)
    du[5]=(-mt+mt_inf(V))*1/tau_mt(V)
    du[6]=(-ht+ht_inf(V))*1/tau_ht(V)
    du[7]=(-mA+mA_inf(V))*1/tau_mA(V)
    du[8]=(-hA+hA_inf(V))*1/tau_hA(V)
    du[9]=1/tau_mKD(V)*(-mKD+mKD_inf(V))
    du[10]=1/tau_mL(V)*(-mL+mL_inf(V))
    du[11]=(-0.1(gL*mL*(V-VCa)+0.0*ICa_pump(Ca))-0.01*Ca)/4
    du[12]=-noise/taunoise
end

function σ_HM(du,u,p,t)
    du[1]=0.0
    du[2]=0.0   #p[16]/tau_m(V)
    du[3]=0.0   #p[17]/tau_h(V)
    du[4]=0.0   #p[18]/tau_mH(V)
    du[5]=0.0   #p[19]/tau_mt(V)
    du[6]=0.0   #p[20]/tau_ht(V)
    du[7]=0.0   #p[21]/tau_mA(V)
    du[8]=0.0   #p[22]/tau_hA(V)
    du[9]=0.0   #p[23]/tau_mKD(V+15)
    du[10]=0.0  #p[24]/tau_mL(V)
    du[11]=0.0
    du[12]=p[17]
end
