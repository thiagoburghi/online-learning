# This model was adapted directly from the model in 
# "A positive feedback at the cellular level promotes robustness and modulation at the circuit level"
# By Dethier et al, 2015.

const ENa = 50.     
const EK = -80.     
const ECa = 120.
const EL = -49.
const ESyn = -80.

gam(v,c,vhalf,sig) = c./(1+exp((v+vhalf)/sig))
tau(v,c1,c2,vhalf,sig) = c1 + c2/(1+exp((v+vhalf)/sig))
tau2(v,c,vhalf,sig,c1,c2,vhalf2,sig2) = c/(1+exp((v+vhalf)/sig))*(c1 + c2/(1+exp((v+vhalf2)/sig2)))

### Na-current (m=activation variable, h=inactivation variable)
# m activation
cNa_m = 1;
vhalfNa_m = 35.5*δ[1,:];
sigNa_m = -5.29*δ[2,:];
gamNa_m(V,i) = gam(V,cNa_m,vhalfNa_m[i],sigNa_m[i])
# m time constant
c1Na_m = 1.32;
c2Na_m = -1.26;
vtauNa_m = 120*δ[3,:];
sigtauNa_m = -25*δ[4,:];
tauNa_m(V,i) = tau(V,c1Na_m,c2Na_m,vtauNa_m[i],sigtauNa_m[i])
# h activation
cNa_h = 1;
vhalfNa_h = 48.9*δ[5,:];
sigNa_h = 5.18*δ[6,:];
gamNa_h(V,i) = gam(V,cNa_h,vhalfNa_h[i],sigNa_h[i])
# h time constant
ctauNa_h = 0.67;
vtauNa_h = 62.9*δ[7,:];
sigtauNa_h = -10*δ[8,:];
c1Na_h = 1.5;
c2Na_h = 1;
vtau2Na_h = 34.9*δ[9,:]; 
sigtau2Na_h = 3.6*δ[10,:];
tauNa_h(V,i) = tau2(V,ctauNa_h,vtauNa_h[i],sigtauNa_h[i],c1Na_h,c2Na_h,vtau2Na_h[i],sigtau2Na_h[i])

### K-current (mK=activation variable)
# mK activation
cK_m = 1;
vhalfK_m = 12.3*δ[11,:];
sigK_m = -11.8*δ[12,:];
gamK_m(V,i) = gam(V,cK_m,vhalfK_m[i],sigK_m[i])
# mK time constant
c1K_m = 7.2;
c2K_m = -6.4;
vtauK_m = 28.3*δ[13,:];
sigtauK_m = -19.2*δ[14,:];
tauK_m(V,i) = tau(V,c1K_m,c2K_m,vtauK_m[i],sigtauK_m[i])

### Ca-current (mCa=activation variable, hCa = inactivation variable)
# mCa activation
cCa_m = 1;
vhalfCa_m = 67.1*δ[15,:];          #dethier: 57.1; follow-up paper: 67.1
sigCa_m = -7.2*δ[16,:];
gamCa_m(V,i)= gam(V,cCa_m,vhalfCa_m[i],sigCa_m[i])
# mCa time constant
c1Ca_m = 43.4;
c2Ca_m = -42.6;
vtauCa_m = 68.1*δ[17,:];
sigtauCa_m = -20.5*δ[18,:];
tauCa_m(V,i) = tau(V,c1Ca_m,c2Ca_m,vtauCa_m[i],sigtauCa_m[i])
# hCa activation
cCa_h = 1;
vhalfCa_h = 82.1*δ[19,:];
sigCa_h = 5.5*δ[20,:];
gamCa_h(V,i) = gam(V,cCa_h,vhalfCa_h[i],sigCa_h[i])
# hCa time constant
c1Ca_h = 140;                                               # dethier: 840 / me 140
c2Ca_h = -100;                                              # dethier: -718.4 / me -100
vtauCa_h = 55*δ[21,:];
sigtauCa_h = -16.9*δ[22,:];
tauCa_h(V,i) = tau(V,c1Ca_h,c2Ca_h,vtauCa_h[i],sigtauCa_h[i])

### Synaptic current
cSyn = 1;
vhalfSyn = 45*δ[23,:];
sigSyn = -2*δ[24,:];
gamSyn(V,i) = gam(V,cSyn,vhalfSyn[i],sigSyn[i])