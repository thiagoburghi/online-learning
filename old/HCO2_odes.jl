# The HCO model was taken from the upcoming book by Drion, Franci and Sepulchre

function HCO_ode(du,u,p,t)
    Iapp=p[1] 
    gT=p[2]
    gKD=p[3]
    gH=p[4]
    gNa=p[5]
    gA=p[6]
    gLeak=p[7]
    gL=p[8]
    gKCa=p[9]
    C=p[10]
    tausyn=p[11]
    gsyn=p[12]
    Esyn=p[13]

    V = u[1:2]
    m = u[3:4]
    h = u[5:6]
    mH = u[7:8]
    mt = u[9:10]
    ht = u[11:12]
    mA = u[13:14]
    hA = u[15:16]
    mKD = u[17:18]
    mL = u[19:20]
    Ca = u[21:22]
    msyn = u[23:24]

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*m[n]^3*h[n]*(V[n]-VNa) - gH[n]*mH[n]*(V[n]-VH) - gT[n]*mt[n]^2*ht[n]*(V[n]-VCa) - gA[n]*mA[n]^4*hA[n]*(V[n]-VK) + 
                        - gKD[n]*mKD[n]^4*(V[n]-VK) - gLeak[n]*(V[n]-Vleak)- gL[n]*mL[n]*(V[n]-VCa) - gKCa[n]*(Ca[n]/(15.0+Ca[n]))^4*(V[n]-VK) + 
                        - gsyn[n]*msyn[n]*(V[n]-Esyn[n]) + Iapp[n](t)  )
        du[2+n]=(-m[n]+m_inf(V[n],n))*1/tau_m(V[n],n)
        du[4+n]=(-h[n]+h_inf(V[n],n))*1/tau_h(V[n],n)
        du[6+n]=(-mH[n]+mH_inf(V[n],n))*1/tau_mH(V[n],n)
        du[8+n]=(-mt[n]+mt_inf(V[n],n))*1/tau_mt(V[n],n)
        du[10+n]=(-ht[n]+ht_inf(V[n],n))*1/tau_ht(V[n],n)
        du[12+n]=(-mA[n]+mA_inf(V[n],n))*1/tau_mA(V[n],n)
        du[14+n]=(-hA[n]+hA_inf(V[n],n))*1/tau_hA(V[n],n)
        du[16+n]=(-mKD[n]+mKD_inf(V[n],n))*1/tau_mKD(V[n],1)
        du[18+n]=(-mL[n]+mL_inf(V[n],n))*1/tau_mL(V[n],n)
        du[20+n]=(-0.01*Ca[n] - 0.1(gL[n]*mL[n]*(V[n]-VCa) + 0.0*ICa_pump(Ca[n])))/4
        du[22+n]=(-msyn[n]+msyn_inf(V[mod(2*n,3)],n))/tausyn[n]
    end
end

function HCO_observer_ode(du,u,p,t)
    Iapp=p[1] 
    gT=p[2]
    gKD=p[3]
    gH=p[4]
    gNa=p[5]
    gA=p[6]
    gLeak=p[7]
    gL=p[8]
    gKCa=p[9]
    C=p[10]
    tausyn=p[11]
    gsyn=p[12]
    Esyn=p[13]
    α₁,α₂,β,γ = p[14]
    mode = p[15]

    V = u[1:2];         V̂ = u[25:26];
    m = u[3:4];         m̂ = u[27:28];
    h = u[5:6];         ĥ = u[29:30];
    mH = u[7:8];        mĤ = u[31:32];
    mt = u[9:10];       mt̂ = u[33:34];
    ht = u[11:12];      ht̂ = u[35:36];
    mA = u[13:14];      mÂ = u[37:38];
    hA = u[15:16];      hÂ = u[39:40];
    mKD = u[17:18];     mKD̂ = u[41:42];
    mL = u[19:20];      mL̂ = u[43:44];
    Ca = u[21:22];      Câ = u[45:46];
    msyn = u[23:24];    msyn̂ = u[47:48];

    gNâ = u[49:50]
    gĤ = u[51:52]
    gT̂ = u[53:54]
    gÂ = u[55:56]
    gKD̂ = u[57:58]
    gL̂ = u[59:60]
    gKCâ = u[61:62]
    gLeak̂ = u[63:64]
    gsyn̂ = u[65:66]
    P1 = reshape(u[67:147],9,9);    P1 = (P1+P1')/2
    P2 = reshape(u[148:228],9,9);   P2 = (P2+P2')/2
    P = (P1,P2)
    ψ = reshape(u[229:246],2,9)

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*m[n]^3*h[n]*(V[n]-VNa) +
                        - gH[n]*mH[n]*(V[n]-VH) + 
                        - gT[n]*mt[n]^2*ht[n]*(V[n]-VCa) +
                        - gA[n]*mA[n]^4*hA[n]*(V[n]-VK) + 
                        - gKD[n]*mKD[n]^4*(V[n]-VK) +
                        - gL[n]*mL[n]*(V[n]-VCa) +
                        - gLeak[n]*(V[n]-Vleak) +
                        - gKCa[n]*(Ca[n]/(15.0+Ca[n]))^4*(V[n]-VK) + 
                        - gsyn[n]*msyn[n]*(V[n]-Esyn[n]) +
                        + Iapp[n](t)  )
        du[2+n]=(-m[n]+m_inf(V[n],n))*1/tau_m(V[n],n)
        du[4+n]=(-h[n]+h_inf(V[n],n))*1/tau_h(V[n],n)
        du[6+n]=(-mH[n]+mH_inf(V[n],n))*1/tau_mH(V[n],n)
        du[8+n]=(-mt[n]+mt_inf(V[n],n))*1/tau_mt(V[n],n)
        du[10+n]=(-ht[n]+ht_inf(V[n],n))*1/tau_ht(V[n],n)
        du[12+n]=(-mA[n]+mA_inf(V[n],n))*1/tau_mA(V[n],n)
        du[14+n]=(-hA[n]+hA_inf(V[n],n))*1/tau_hA(V[n],n)
        du[16+n]=(-mKD[n]+mKD_inf(V[n],n))*1/tau_mKD(V[n],n)
        du[18+n]=(-mL[n]+mL_inf(V[n],n))*1/tau_mL(V[n],n)
        du[20+n]=(-0.01*Ca[n] - 0.1(gL[n]*mL[n]*(V[n]-VCa) + 0.0*ICa_pump(Ca[n])))/4
        du[22+n]=(-msyn[n]+msyn_inf(V[mod(2*n,3)],n))/tausyn[n]
    end

    # OBSERVER
    dψ = zeros(2,9) 
    dθ = zeros(2,9)
    dP = zeros(9,9,2)

    if mode == "output"
        Vᵢ = V  # ORIGINALLY V
    else
        Vᵢ = V̂
    end
    for n = 1:2
        p = n+2
        du[24+n]=1/C[n]*(- gNâ[n]*m̂[n]^3*ĥ[n]*(Vᵢ[n]-VNa) +
                         - gĤ[n]*mĤ[n]*(Vᵢ[n]-VH) + 
                         - gT̂[n]*mt̂[n]^2*ht̂[n]*(Vᵢ[n]-VCa) +
                         - gÂ[n]*mÂ[n]^4*hÂ[n]*(Vᵢ[n]-VK) + 
                         - gKD̂[n]*mKD̂[n]^4*(Vᵢ[n]-VK) +
                         - gL̂[n]*mL̂[n]*(Vᵢ[n]-VCa) +
                         - gKCâ[n]*(Câ[n]/(15.0+Câ[n]))^4*(Vᵢ[n]-VK) + 
                         - gLeak̂[n]*(Vᵢ[n]-Vleak) +
                         - gsyn̂[n]*msyn̂[n]*(Vᵢ[n]-Esyn[n]) +
                         + Iapp[n](t)) + 
                         + γ*(V[n]-V̂[n]) +
                         + β*ψ[n,:]'*P[n]*ψ[n,:]*(V[n]-V̂[n])
        du[26+n]=(-m̂[n]+m_inf(Vᵢ[n],p))*1/tau_m(Vᵢ[n],p)
        du[28+n]=(-ĥ[n]+h_inf(Vᵢ[n],p))*1/tau_h(Vᵢ[n],p)
        du[30+n]=(-mĤ[n]+mH_inf(Vᵢ[n],p))*1/tau_mH(Vᵢ[n],p)
        du[32+n]=(-mt̂[n]+mt_inf(Vᵢ[n],p))*1/tau_mt(Vᵢ[n],p)
        du[34+n]=(-ht̂[n]+ht_inf(Vᵢ[n],p))*1/tau_ht(Vᵢ[n],p)
        du[36+n]=(-mÂ[n]+mA_inf(Vᵢ[n],p))*1/tau_mA(Vᵢ[n],p)
        du[38+n]=(-hÂ[n]+hA_inf(Vᵢ[n],p))*1/tau_hA(Vᵢ[n],p)
        du[40+n]=(-mKD̂[n]+mKD_inf(Vᵢ[n],p))*1/tau_mKD(Vᵢ[n],p)
        du[42+n]=(-mL̂[n]+mL_inf(Vᵢ[n],p))*1/tau_mL(Vᵢ[n],p)
        du[44+n]=(-0.01*Câ[n] - 0.1(gL[n]*mL̂[n]*(Vᵢ[n]-VCa) + 0.0*ICa_pump(Câ[n])))/4
        du[46+n]=(-msyn̂[n]+msyn_inf(Vᵢ[mod(2*n,3)],p))/tausyn[n]

        ϕ = -1/C[n]*[m̂[n]^3*ĥ[n]*(Vᵢ[n]-VNa)
                    mĤ[n]*(Vᵢ[n]-VH) 
                    mt̂[n]^2*ht̂[n]*(Vᵢ[n]-VCa)
                    mÂ[n]^4*hÂ[n]*(Vᵢ[n]-VK)
                    mKD̂[n]^4*(Vᵢ[n]-VK)
                    mL̂[n]*(Vᵢ[n]-VCa)
                    (Câ[n]/(15.0+Câ[n]))^4*(Vᵢ[n]-VK)
                    (Vᵢ[n]-Vleak)
                    msyn̂[n]*(Vᵢ[n]-Esyn[n])]

        dψ[n,:] = -β*ψ[n,:] + ϕ
        dθ[n,:] = β*P[n]*ψ[n,:]*(V[n]-V̂[n])
        dP[:,:,n] = α₁*P[n] - α₂*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        dP[:,:,n] = (dP[:,:,n]+dP[:,:,n]')/2
    end
    du[49:66]   = reshape(dθ,18)
    du[67:147]  = reshape(dP[:,:,1],81)
    du[148:228] = reshape(dP[:,:,2],81)
    du[229:246] = reshape(dψ,18)
end

function HCO_observer_besancon(du,u,p,t)
    Iapp=p[1] 
    gT=p[2]
    gKD=p[3]
    gH=p[4]
    gNa=p[5]
    gA=p[6]
    gLeak=p[7]
    gL=p[8]
    gKCa=p[9]
    C=p[10]
    tausyn=p[11]
    gsyn=p[12]
    Esyn=p[13]
    α₁,α₂,β,γ = p[14]

    V = u[1:2];         V̂ = u[25:26];
    m = u[3:4];         m̂ = u[27:28];
    h = u[5:6];         ĥ = u[29:30];
    mH = u[7:8];        mĤ = u[31:32];
    mt = u[9:10];       mt̂ = u[33:34];
    ht = u[11:12];      ht̂ = u[35:36];
    mA = u[13:14];      mÂ = u[37:38];
    hA = u[15:16];      hÂ = u[39:40];
    mKD = u[17:18];     mKD̂ = u[41:42];
    mL = u[19:20];      mL̂ = u[43:44];
    Ca = u[21:22];      Câ = u[45:46];
    msyn = u[23:24];    msyn̂ = u[47:48];

    gNâ = u[49:50]
    gĤ = u[51:52]
    gT̂ = u[53:54]
    gÂ = u[55:56]
    gKD̂ = u[57:58]
    gL̂ = u[59:60]
    gKCâ = u[61:62]
    gLeak̂ = u[63:64]

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*m[n]^3*h[n]*(V[n]-VNa) +
                        - gH[n]*mH[n]*(V[n]-VH) + 
                        - gT[n]*mt[n]^2*ht[n]*(V[n]-VCa) +
                        - gA[n]*mA[n]^4*hA[n]*(V[n]-VK) + 
                        - gKD[n]*mKD[n]^4*(V[n]-VK) +
                        - gL[n]*mL[n]*(V[n]-VCa) +
                        - gLeak[n]*(V[n]-Vleak) +
                        - gKCa[n]*(Ca[n]/(15.0+Ca[n]))^4*(V[n]-VK) + 
                        - gsyn[n]*msyn[n]*(V[n]-Esyn[n]) +
                        + Iapp[n](t)  )
        du[2+n]=(-m[n]+m_inf(V[n],n))*1/tau_m(V[n],n)
        du[4+n]=(-h[n]+h_inf(V[n],n))*1/tau_h(V[n],n)
        du[6+n]=(-mH[n]+mH_inf(V[n],n))*1/tau_mH(V[n],n)
        du[8+n]=(-mt[n]+mt_inf(V[n],n))*1/tau_mt(V[n],n)
        du[10+n]=(-ht[n]+ht_inf(V[n],n))*1/tau_ht(V[n],n)
        du[12+n]=(-mA[n]+mA_inf(V[n],n))*1/tau_mA(V[n],n)
        du[14+n]=(-hA[n]+hA_inf(V[n],n))*1/tau_hA(V[n],n)
        du[16+n]=(-mKD[n]+mKD_inf(V[n],n))*1/tau_mKD(V[n],n)
        du[18+n]=(-mL[n]+mL_inf(V[n],n))*1/tau_mL(V[n],n)
        du[20+n]=(-0.01*Ca[n] - 0.1(gL[n]*mL[n]*(V[n]-VCa) + 0.0*ICa_pump(Ca[n])))/4
        du[22+n]=(-msyn[n]+msyn_inf(V[mod(2*n,3)]))/tausyn[n]
    end

    # OBSERVER
    dθ = zeros(2,8)

    for n = 1:2
        p = n+2
        Vᵢ = V̂  # ORIGINALLY V
        du[24+n]=1/C[n]*(- gNâ[n]*m̂[n]^3*ĥ[n]*(Vᵢ[n]-VNa) +
                         - gĤ[n]*mĤ[n]*(Vᵢ[n]-VH) + 
                         - gT̂[n]*mt̂[n]^2*ht̂[n]*(Vᵢ[n]-VCa) +
                         - gÂ[n]*mÂ[n]^4*hÂ[n]*(Vᵢ[n]-VK) + 
                         - gKD̂[n]*mKD̂[n]^4*(Vᵢ[n]-VK) +
                         - gL̂[n]*mL̂[n]*(Vᵢ[n]-VCa) +
                         - gKCâ[n]*(Câ[n]/(15.0+Câ[n]))^4*(Vᵢ[n]-VK) + 
                         - gLeak̂[n]*(Vᵢ[n]-Vleak) +
                         - gsyn[n]*msyn̂[n]*(Vᵢ[n]-Esyn[n]) +
                         + γ*(Vᵢ[n]-V̂[n]) +
                         + Iapp[n](t))
        du[26+n]=(-m̂[n]+m_inf(V[n],p))*1/tau_m(V[n],p)
        du[28+n]=(-ĥ[n]+h_inf(V[n],p))*1/tau_h(V[n],p)
        du[30+n]=(-mĤ[n]+mH_inf(V[n],p))*1/tau_mH(V[n],p)
        du[32+n]=(-mt̂[n]+mt_inf(V[n],p))*1/tau_mt(V[n],p)
        du[34+n]=(-ht̂[n]+ht_inf(V[n],p))*1/tau_ht(V[n],p)
        du[36+n]=(-mÂ[n]+mA_inf(V[n],p))*1/tau_mA(V[n],p)
        du[38+n]=(-hÂ[n]+hA_inf(V[n],p))*1/tau_hA(V[n],p)
        du[40+n]=(-mKD̂[n]+mKD_inf(V[n],p))*1/tau_mKD(V[n],p)
        du[42+n]=(-mL̂[n]+mL_inf(V[n],p))*1/tau_mL(V[n],p)
        du[44+n]=(-0.01*Câ[n] - 0.1(gL[n]*mL̂[n]*(V[n]-VCa) + 0.0*ICa_pump(Câ[n])))/4
        du[46+n]=(-msyn̂[n]+msyn_inf(V[mod(2*n,3)]))/tausyn[n]

        ϕ = -1/C[n]*[m̂[n]^3*ĥ[n]*(V[n]-VNa)
                    mĤ[n]*(V[n]-VH) 
                    mt̂[n]^2*ht̂[n]*(V[n]-VCa)
                    mÂ[n]^4*hÂ[n]*(V[n]-VK)
                    mKD̂[n]^4*(V[n]-VK)
                    mL̂[n]*(V[n]-VCa)
                    (Câ[n]/(15.0+Câ[n]))^4*(V[n]-VK)
                    (V[n]-Vleak)]

        dθ[n,:] = β*ϕ*(V[n]-V̂[n])
    end
    du[49:64]   = reshape(dθ,16)
end