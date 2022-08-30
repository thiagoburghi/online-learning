function ς(x,lb,ub)
    return min(ub,max(lb,x))
end

function HCO_observer_ode(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    a = p[8]
    b = p[9]
    α,β,γ = p[10]
    η,η̂ = p[11]
    outinj = p[12]

    V = u[1:2];          V̂ = u[15:16]
    mNa = u[3:4];        m̂Na = u[17:18]
    hNa = u[5:6];        ĥNa = u[19:20]
    mK = u[7:8];         m̂K = u[21:22]
    mCa = u[9:10];       m̂Ca = u[23:34]
    hCa = u[11:12];      ĥCa = u[25:26]
    mSyn = u[13:14];     m̂Syn = u[27:28]

    ĝNa = u[29:30];      ĝNaₛ = ς.(ĝNa,0,200)
    ĝK = u[31:32];       ĝKₛ = ς.(ĝK,0,200)       
    ĝCa = u[33:34];      ĝCaₛ = ς.(ĝCa,0,200)
    ĝL = u[35:36];       ĝLₛ = ς.(ĝL,0,200)
    ĝSyn = u[37:38];     ĝSynₛ = ς.(ĝSyn,0,200)
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
        du[2+n]=(-mNa[n]+σNa_m(V[n],η[n]))*1/τNa_m(V[n],η[n])
        du[4+n]=(-hNa[n]+σNa_h(V[n],η[n]))*1/τNa_h(V[n],η[n])
        du[6+n]=(-mK[n]+σK_m(V[n],η[n]))*1/τK_m(V[n],η[n])
        du[8+n]=(-mCa[n]+σCa_m(V[n],η[n]))*1/τCa_m(V[n],η[n])
        du[10+n]=(-hCa[n]+σCa_h(V[n],η[n]))*1/τCa_h(V[n],η[n])
        du[12+n] = a[n]*σSyn(V[mod(2*n,3)],η[n])*(1-mSyn[n])-b[n]*mSyn[n] 
    end

    # OBSERVER
    dψ = zeros(2,5) 
    dθ = zeros(2,5)
    dP = zeros(5,5,2)

    Y = V

    for n = 1:2

        if outinj
            Yᵥ = V 
            ∂₁v̂̇ = 0.0
        else
            Yᵥ = V̂
            ∂₁v̂̇ = -ĝNaₛ[1]*m̂Na[n]^3*ĥNa[n] - ĝKₛ[n]*m̂K[n]^4 - ĝCaₛ[n]*m̂Ca[n]^3*ĥCa[n] - ĝLₛ[n] - ĝSynₛ[n]*m̂Syn[n]  
        end

        du[14+n]=1/C[n]*(- ĝNa[n]*m̂Na[n]^3*ĥNa[n]*(Yᵥ[n]-ENa) + 
                         - ĝK[n]*m̂K[n]^4*(Yᵥ[n]-EK) +
                         - ĝCa[n]*m̂Ca[n]^3*ĥCa[n]*(Yᵥ[n]-ECa) +
                         - ĝL[n]*(Yᵥ[n]-EL) +  
                         - ĝSyn[n]*m̂Syn[n]*(Yᵥ[n]-ESyn) + 
                         + Iapp[n](t)) + 
                         + γ*(Y[n]-V̂[n]) +
                         + γ*ψ[n,:]'*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        du[16+n]=(-m̂Na[n]+σNa_m(Y[n],η̂[n]))*1/τNa_m(Y[n],η̂[n])
        du[18+n]=(-ĥNa[n]+σNa_h(Y[n],η̂[n]))*1/τNa_h(Y[n],η̂[n])
        du[20+n]=(-m̂K[n]+σK_m(Y[n],η̂[n]))*1/τK_m(Y[n],η̂[n])
        du[22+n]=(-m̂Ca[n]+σCa_m(Y[n],η̂[n]))*1/τCa_m(Y[n],η̂[n])
        du[24+n]=(-ĥCa[n]+σCa_h(Y[n],η̂[n]))*1/τCa_h(Y[n],η̂[n])
        du[26+n] = a[n]*σSyn(Y[mod(2*n,3)],η̂[n])*(1-m̂Syn[n])-b[n]*m̂Syn[n] 

        ϕ = -1/C[n]*[m̂Na[n]^3*ĥNa[n]*(Yᵥ[n]-ENa),
                     m̂K[n]^4*(Yᵥ[n]-EK),
                     m̂Ca[n]^3*ĥCa[n]*(Yᵥ[n]-ECa),
                     (Yᵥ[n]-EL),
                     m̂Syn[n]*(Yᵥ[n]-ESyn)]
           
        du[28+n:2:28+10]   = γ*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        dP = α*P[n] + β*Matrix(I, 5, 5) - γ*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        dP = (dP+dP')/2        
        du[39+(n-1)*25:39+(n-1)*25+24] = reshape(dP,25)
        du[88+n:2:88+10] = (-γ+∂₁v̂̇)*ψ[n,:] + ϕ

        # dψ[n,:] = -γ*ψ[n,:] + ϕ
        # dθ[n,:] = γ*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        # dP[:,:,n] = α*P[n] + β*Matrix(I, 5, 5) - γ*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        # dP[:,:,n] = (dP[:,:,n]+dP[:,:,n]')/2
    end
    # du[29:38]   = reshape(dθ,10)
    # du[39:63]  = reshape(dP[:,:,1],25)
    # du[64:88] = reshape(dP[:,:,2],25)
    # du[89:98] = reshape(dψ,10)
end

function HCO_ode(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    a = p[8]
    b = p[9]
    η,η̂ = p[11]

    V = u[1:2];
    mNa = u[3:4];
    hNa = u[5:6];
    mK = u[7:8];
    mCa = u[9:10];
    hCa = u[11:12];
    mSyn = u[13:14];

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*mNa[n]^3*hNa[n]*(V[n]-ENa) + 
                        - gK[n]*mK[n]^4*(V[n]-EK) +
                        - gCa[n](t)*mCa[n]^3*hCa[n]*(V[n]-ECa) +
                        - gL[n]*(V[n]-EL) +  
                        - gSyn[n]*mSyn[n]*(V[n]-ESyn) +
                        + Iapp[n](t)  )
        du[2+n]=(-mNa[n]+σNa_m(V[n],η[n]))*1/τNa_m(V[n],η[n])
        du[4+n]=(-hNa[n]+σNa_h(V[n],η[n]))*1/τNa_h(V[n],η[n])
        du[6+n]=(-mK[n]+σK_m(V[n],η[n]))*1/τK_m(V[n],η[n])
        du[8+n]=(-mCa[n]+σCa_m(V[n],η[n]))*1/τCa_m(V[n],η[n])
        du[10+n]=(-hCa[n]+σCa_h(V[n],η[n]))*1/τCa_h(V[n],η[n])
        du[12+n] = a[n]*σSyn(V[mod(2*n,3)],η[n])*(1-mSyn[n])-b[n]*mSyn[n] 
    end
end

function f_SDE_no_observer!(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    a = p[8]
    b = p[9]
    η,η̂ = p[11]

    V = u[1:2];
    mNa = u[3:4];
    hNa = u[5:6];
    mK = u[7:8];
    mCa = u[9:10];
    hCa = u[11:12];
    mSyn = u[13:14];

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*mNa[n]^3*hNa[n]*(V[n]-ENa) + 
                        - gK[n]*mK[n]^4*(V[n]-EK) +
                        - gCa[n](t)*mCa[n]^3*hCa[n]*(V[n]-ECa) +
                        - gL[n]*(V[n]-EL) +  
                        - gSyn[n]*mSyn[n]*(V[n]-ESyn) +
                        + Iapp[n](t)  )
        du[2+n]=(-mNa[n]+σNa_m(V[n],η[n]))*1/τNa_m(V[n],η[n])
        du[4+n]=(-hNa[n]+σNa_h(V[n],η[n]))*1/τNa_h(V[n],η[n])
        du[6+n]=(-mK[n]+σK_m(V[n],η[n]))*1/τK_m(V[n],η[n])
        du[8+n]=(-mCa[n]+σCa_m(V[n],η[n]))*1/τCa_m(V[n],η[n])
        du[10+n]=(-hCa[n]+σCa_h(V[n],η[n]))*1/τCa_h(V[n],η[n])
        du[12+n] = a[n]*σSyn(V[mod(2*n,3)],η[n])*(1-mSyn[n])-b[n]*mSyn[n] 
    end
end

function g_SDE_no_observer!(du,u,p,t)
    du[:] = zeros(14,1)
end

function f_SDE!(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    gSyn = p[7]
    a = p[8]
    b = p[9]
    α,β,γ = p[10]
    η,η̂ = p[11]
    outinj = p[12]
    σnoise,dt = p[13]

    V = u[1:2];          V̂ = u[15:16]
    mNa = u[3:4];        m̂Na = u[17:18]
    hNa = u[5:6];        ĥNa = u[19:20]
    mK = u[7:8];         m̂K = u[21:22]
    mCa = u[9:10];       m̂Ca = u[23:34]
    hCa = u[11:12];      ĥCa = u[25:26]
    mSyn = u[13:14];     m̂Syn = u[27:28]

    ĝNa = u[29:30];      ĝNaₛ = ς.(ĝNa,0,200)
    ĝK = u[31:32];       ĝKₛ = ς.(ĝK,0,200)       
    ĝCa = u[33:34];      ĝCaₛ = ς.(ĝCa,0,200)
    ĝL = u[35:36];       ĝLₛ = ς.(ĝL,0,200)
    ĝSyn = u[37:38];     ĝSynₛ = ς.(ĝSyn,0,200)
    P1 = reshape(u[39:63],5,5);    P1 = (P1+P1')/2
    P2 = reshape(u[64:88],5,5);   P2 = (P2+P2')/2
    P = (P1,P2)
    ψ = reshape(u[89:98],2,5)
    noise = u[99:100]

    # TRUE SYSTEM
    for n = 1:2
        du[0+n]=1/C[n]*(- gNa[n]*mNa[n]^3*hNa[n]*(V[n]-ENa) + 
                        - gK[n]*mK[n]^4*(V[n]-EK) +
                        - gCa[n](t)*mCa[n]^3*hCa[n]*(V[n]-ECa) +
                        - gL[n]*(V[n]-EL) +  
                        - gSyn[n]*mSyn[n]*(V[n]-ESyn) +
                        + Iapp[n](t)  )
        du[2+n]=(-mNa[n]+σNa_m(V[n],η[n]))*1/τNa_m(V[n],η[n])
        du[4+n]=(-hNa[n]+σNa_h(V[n],η[n]))*1/τNa_h(V[n],η[n])
        du[6+n]=(-mK[n]+σK_m(V[n],η[n]))*1/τK_m(V[n],η[n])
        du[8+n]=(-mCa[n]+σCa_m(V[n],η[n]))*1/τCa_m(V[n],η[n])
        du[10+n]=(-hCa[n]+σCa_h(V[n],η[n]))*1/τCa_h(V[n],η[n])
        du[12+n] = a[n]*σSyn(V[mod(2*n,3)],η[n])*(1-mSyn[n])-b[n]*mSyn[n] 
    end

    # OBSERVER
    Y  = V + noise
    if outinj
        Yᵥ = Y
    else
        Yᵥ = V̂
    end
    # dψ = zeros(2,5) 
    # dθ = zeros(2,5)
    # dP = zeros(5,5,2)

    for n = 1:2
        if outinj
            ∂₁v̂̇ = 0.0
        else
            ∂₁v̂̇ = -ĝNaₛ[n]*m̂Na[n]^3*ĥNa[n] - ĝKₛ[n]*m̂K[n]^4 - ĝCaₛ[n]*m̂Ca[n]^3*ĥCa[n] - ĝLₛ[n] - ĝSynₛ[n]*m̂Syn[n]  
        end
        du[14+n]=1/C[n]*(- ĝNa[n]*m̂Na[n]^3*ĥNa[n]*(Yᵥ[n]-ENa) + 
                         - ĝK[n]*m̂K[n]^4*(Yᵥ[n]-EK) +
                         - ĝCa[n]*m̂Ca[n]^3*ĥCa[n]*(Yᵥ[n]-ECa) +
                         - ĝL[n]*(Yᵥ[n]-EL) +  
                         - ĝSyn[n]*m̂Syn[n]*(Yᵥ[n]-ESyn) + 
                         + Iapp[n](t)) + 
                         + γ*(Y[n]-V̂[n]) +
                         + γ*ψ[n,:]'*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        du[16+n]=(-m̂Na[n]+σNa_m(Y[n],η̂[n]))*1/τNa_m(Y[n],η̂[n])
        du[18+n]=(-ĥNa[n]+σNa_h(Y[n],η̂[n]))*1/τNa_h(Y[n],η̂[n])
        du[20+n]=(-m̂K[n]+σK_m(Y[n],η̂[n]))*1/τK_m(Y[n],η̂[n])
        du[22+n]=(-m̂Ca[n]+σCa_m(Y[n],η̂[n]))*1/τCa_m(Y[n],η̂[n])
        du[24+n]=(-ĥCa[n]+σCa_h(Y[n],η̂[n]))*1/τCa_h(Y[n],η̂[n])
        du[26+n] = a[n]*σSyn(Y[mod(2*n,3)],η̂[n])*(1-m̂Syn[n])-b[n]*m̂Syn[n] 

        ϕ = -1/C[n]*[m̂Na[n]^3*ĥNa[n]*(Yᵥ[n]-ENa),
                     m̂K[n]^4*(Yᵥ[n]-EK),
                     m̂Ca[n]^3*ĥCa[n]*(Yᵥ[n]-ECa),
                     (Yᵥ[n]-EL),
                     m̂Syn[n]*(Yᵥ[n]-ESyn)]
           
        du[28+n:2:28+10]   = γ*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        dP = α*P[n] + β*Matrix(I, 5, 5) - γ*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        dP = (dP+dP')/2        
        du[39+(n-1)*25:39+(n-1)*25+24] = reshape(dP,25)
        du[88+n:2:88+10] = (-γ+∂₁v̂̇)*ψ[n,:] + ϕ
        # dψ[n,:] = -γ*ψ[n,:] + ϕ
        # dθ[n,:] = γ*P[n]*ψ[n,:]*(Y[n]-V̂[n])
        # dP[:,:,n] = α*P[n] + β*Matrix(I, 5, 5) - γ*(P[n]*ψ[n,:])*(P[n]*ψ[n,:])'
        # dP[:,:,n] = (dP[:,:,n]+dP[:,:,n]')/2
    end
    du[99:100] = -1/dt*noise;
    # du[29:38]   = reshape(dθ,10)
    # du[39:63]  = reshape(dP[:,:,1],25)
    # du[64:88] = reshape(dP[:,:,2],25)
    # du[89:98] = reshape(dψ,10)
end

function g_SDE!(du,u,p,t)
    σnoise,dt = p[13]

    # Johnson-Nyquist noise
    gnoise = σnoise/sqrt(dt);

    du[:] = [zeros(88+10,1);gnoise;gnoise]
end