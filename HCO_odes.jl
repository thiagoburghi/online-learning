# The HCO model was adapted directly from the model in 
# "A positive feedback at the cellular level promotes robustness and modulation at the circuit level"
# By Dethier et al, 2015.

function HCO_ode(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    kfA = p[8]
    krA = p[9]
    type = p[10]
    
    V = u[1:2]
    mNa = u[3:4]
    hNa = u[5:6]
    mK = u[7:8]
    mCa = u[9:10]
    hCa = u[11:12]
    mSyn = u[13:14]

    # Choose true model or adaptive observer internal dynamics
    if type == "nominal"
        addind = 0
    else
        addind = 2
    end

    # TRUE SYSTEM
    for n = 1:2
        p = n + addind
        du[0+n]=1/C[n]*(- gNa[n]*mNa[n]^3*hNa[n]*(V[n]-ENa) + 
                        - gK[n]*mK[n]^4*(V[n]-EK) +
                        - gCa[n](t)*mCa[n]^3*hCa[n]*(V[n]-ECa) +
                        - gL[n]*(V[n]-EL) +  
                        - gSyn[n]*mSyn[n]*(V[n]-ESyn) +
                        + Iapp[n](t)  )
        du[2+n]=(-mNa[n]+gamNa_m(V[n],p))*1/tauNa_m(V[n],p)
        du[4+n]=(-hNa[n]+gamNa_h(V[n],p))*1/tauNa_h(V[n],p)
        du[6+n]=(-mK[n]+gamK_m(V[n],p))*1/tauK_m(V[n],p)
        du[8+n]=(-mCa[n]+gamCa_m(V[n],p))*1/tauCa_m(V[n],p)
        du[10+n]=(-hCa[n]+gamCa_h(V[n],p))*1/tauCa_h(V[n],p)
        du[12+n] = kfA[n]*gamSyn(V[mod(2*n,3)],p)*(1-mSyn[n])-krA[n]*mSyn[n] 
    end
end

function HCO_observer_ode(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    kfA = p[8]
    krA = p[9]
    α₁,α₂,β,γ = p[10]

    V = u[1:2];          V̂ = u[15:16]
    mNa = u[3:4];        m̂Na = u[17:18]
    hNa = u[5:6];        ĥNa = u[19:20]
    mK = u[7:8];         m̂K = u[21:22]
    mCa = u[9:10];       m̂Ca = u[23:34]
    hCa = u[11:12];      ĥCa = u[25:26]
    mSyn = u[13:14];     m̂Syn = u[27:28]

    ĝNa = u[29:30]
    ĝK = u[31:32]
    ĝCa = u[33:34]
    ĝL = u[35:36]
    ĝSyn = u[37:38]
    P1 = reshape(u[39:63],5,5);    P1 = (P1+P1')/2
    P2 = reshape(u[64:88],5,5);   P2 = (P2+P2')/2
    P = (P1,P2)
    ψ = reshape(u[89:98],2,5)

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*mNa[n]^3*hNa[n]*(V[n]-ENa) + 
                        - gK[n]*mK[n]^4*(V[n]-EK) +
                        - gCa[n](t)*mCa[n]^3*hCa[n]*(V[n]-ECa) +
                        - gL[n]*(V[n]-EL) +  
                        - gSyn[n]*mSyn[n]*(V[n]-ESyn) +
                        + Iapp[n](t)  )
        du[2+n]=(-mNa[n]+gamNa_m(V[n],n))*1/tauNa_m(V[n],n)
        du[4+n]=(-hNa[n]+gamNa_h(V[n],n))*1/tauNa_h(V[n],n)
        du[6+n]=(-mK[n]+gamK_m(V[n],n))*1/tauK_m(V[n],n)
        du[8+n]=(-mCa[n]+gamCa_m(V[n],n))*1/tauCa_m(V[n],n)
        du[10+n]=(-hCa[n]+gamCa_h(V[n],n))*1/tauCa_h(V[n],n)
        du[12+n] = kfA[n]*gamSyn(V[mod(2*n,3)],n)*(1-mSyn[n])-krA[n]*mSyn[n] 
    end

    # OBSERVER
    dψ = zeros(2,5) 
    dθ = zeros(2,5)
    dP = zeros(5,5,2)

    for n = 1:2
        p = n+2
        du[14+n]=1/C[n]*(- ĝNa[n]*m̂Na[n]^3*ĥNa[n]*(V[n]-ENa) + 
                         - ĝK[n]*m̂K[n]^4*(V[n]-EK) +
                         - ĝCa[n]*m̂Ca[n]^3*ĥCa[n]*(V[n]-ECa) +
                         - ĝL[n]*(V[n]-EL) +  
                         - ĝSyn[n]*m̂Syn[n]*(V[n]-ESyn) + 
                         + Iapp[n](t)) + 
                         + γ*(V[n]-V̂[n]) +
                         + β*ψ[n,:]'*P[n]*ψ[n,:]*(V[n]-V̂[n])
        du[16+n]=(-m̂Na[n]+gamNa_m(V[n],p))*1/tauNa_m(V[n],p)
        du[18+n]=(-ĥNa[n]+gamNa_h(V[n],p))*1/tauNa_h(V[n],p)
        du[20+n]=(-m̂K[n]+gamK_m(V[n],p))*1/tauK_m(V[n],p)
        du[22+n]=(-m̂Ca[n]+gamCa_m(V[n],p))*1/tauCa_m(V[n],p)
        du[24+n]=(-ĥCa[n]+gamCa_h(V[n],p))*1/tauCa_h(V[n],p)
        du[26+n] = kfA[n]*gamSyn(V[mod(2*n,3)],p)*(1-m̂Syn[n])-krA[n]*m̂Syn[n] 

        ϕ = -1/C[n]*[m̂Na[n]^3*ĥNa[n]*(V[n]-ENa),
                     m̂K[n]^4*(V[n]-EK),
                     m̂Ca[n]^3*ĥCa[n]*(V[n]-ECa),
                     (V[n]-EL),
                     m̂Syn[n]*(V[n]-ESyn)]

        dψ[n,:] = -β*ψ[n,:] + ϕ
        dθ[n,:] = β*P[n]*ψ[n,:]*(V[n]-V̂[n])
        dP[:,:,n] = α₁*P[n] - α₂*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        dP[:,:,n] = (dP[:,:,n]+dP[:,:,n]')/2
    end
    du[29:38]   = reshape(dθ,10)
    du[39:63]  = reshape(dP[:,:,1],25)
    du[64:88] = reshape(dP[:,:,2],25)
    du[89:98] = reshape(dψ,10)
end