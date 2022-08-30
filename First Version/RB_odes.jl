function RB_ode!(du,u,p,t)
    Iapp=p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    
    V = u[1]
    mNa = u[2]
    hNa = u[3]
    mK = u[4]
    mCa = u[5]
    hCa = u[6]

    # TRUE SYSTEM
    du[1]=1/C*(- gNa*mNa^3*hNa*(V-ENa) + 
                    - gK*mK^4*(V-EK) +
                    - gCa*mCa^3*hCa*(V-ECa) +
                    - gL*(V-EL) +  
                    + Iapp(t)  )
    du[2]=(-mNa+σNa_m(V))*1/τNa_m(V)
    du[3]=(-hNa+σNa_h(V))*1/τNa_h(V)
    du[4]=(-mK+σK_m(V))*1/τK_m(V)
    du[5]=(-mCa+σCa_m(V))*1/τCa_m(V)
    du[6]=(-hCa+σCa_h(V))*1/τCa_h(V)
end

function RB_observer_ode!(du,u,p,t)
    Iapp = p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    α₁,α₂,γ = p[7]

    V = u[1];        V̂ = u[7]
    mNa = u[2];      m̂Na = u[8]
    hNa = u[3];      ĥNa = u[9]
    mK = u[4];       m̂K = u[10]
    mCa = u[5];      m̂Ca = u[11]
    hCa = u[6];      ĥCa = u[12]

    ĝ = u[13:16]
    P = reshape(u[17:32],4,4);    P = (P+P')/2
    Ψ = u[33:36]

    # TRUE SYSTEM
    du[1]=1/C*(- gNa*mNa^3*hNa*(V-ENa) + 
                    - gK*mK^4*(V-EK) +
                    - gCa*mCa^3*hCa*(V-ECa) +
                    - gL*(V-EL) +   
                    + Iapp(t)  )
    du[2]=(-mNa+σNa_m(V))*1/τNa_m(V)
    du[3]=(-hNa+σNa_h(V))*1/τNa_h(V)
    du[4]=(-mK+σK_m(V))*1/τK_m(V)
    du[5]=(-mCa+σCa_m(V))*1/τCa_m(V)
    du[6]=(-hCa+σCa_h(V))*1/τCa_h(V)

    # ADAPTIVE OBSERVER
    dP = zeros(4,4)

    Φ = -1/C*[m̂Na^3*ĥNa*(V-ENa),
              m̂K^4*(V-EK),
              m̂Ca^3*ĥCa*(V-ECa),
              (V-EL)]
    du[7]= dot(Φ,ĝ) + 1/C*Iapp(t) + γ*(V-V̂) + γ*Ψ'*P*Ψ*(V-V̂)
    du[8]=(-m̂Na+σNa_m(V))*1/τNa_m(V)
    du[9]=(-ĥNa+σNa_h(V))*1/τNa_h(V)
    du[10]=(-m̂K+σK_m(V))*1/τK_m(V)
    du[11]=(-m̂Ca+σCa_m(V))*1/τCa_m(V)
    du[12]=(-ĥCa+σCa_h(V))*1/τCa_h(V)

    dΨ = -γ*Ψ + Φ 
    dθ = γ*P*Ψ*(V-V̂)
    dP = α₁*P - α₂*(P*Ψ)*(P*Ψ)';        
    dP = (dP+dP')/2

    du[13:16] = dθ[:]
    du[17:32] = reshape(dP,1,16)
    du[33:36] = dΨ[:]
end

function RB_controlled_ode!(du,u,p,t)
    Īapp = p[1] 
    C = p[2]
    gNa = p[3]
    gK = p[4]
    gCa = p[5]
    gL = p[6]
    α₁,α₂,γ = p[7]
    ḡCa = p[8]

    V = u[1];        V̂ = u[7]
    mNa = u[2];      m̂Na = u[8]
    hNa = u[3];      ĥNa = u[9]
    mK = u[4];       m̂K = u[10]
    mCa = u[5];      m̂Ca = u[11]
    hCa = u[6];      ĥCa = u[12]

    ĝ = u[13:16]
    ĝCa = ĝ[3]
    P = reshape(u[17:32],4,4);    P = (P+P')/2
    Ψ = u[33:36]

    # CONTROL INPUT
    Iapp = -(ḡCa(t)-ĝCa)*m̂Ca^3*ĥCa*(V-ECa) + Īapp(t)

    # TRUE SYSTEM
    du[1]=1/C*(- gNa*mNa^3*hNa*(V-ENa) + 
                    - gK*mK^4*(V-EK) +
                    - gCa*mCa^3*hCa*(V-ECa) +
                    - gL*(V-EL) +   
                    + Iapp  )
    du[2]=(-mNa+σNa_m(V))*1/τNa_m(V)
    du[3]=(-hNa+σNa_h(V))*1/τNa_h(V)
    du[4]=(-mK+σK_m(V))*1/τK_m(V)
    du[5]=(-mCa+σCa_m(V))*1/τCa_m(V)
    du[6]=(-hCa+σCa_h(V))*1/τCa_h(V)

    # ADAPTIVE OBSERVER
    dP = zeros(4,4)

    Φ = -1/C*[m̂Na^3*ĥNa*(V-ENa),
              m̂K^4*(V-EK),
              m̂Ca^3*ĥCa*(V-ECa),
              (V-EL)]
    du[7]= dot(Φ,ĝ) + 1/C*Iapp + γ*(V-V̂) + γ*Ψ'*P*Ψ*(V-V̂)
    du[8]=(-m̂Na+σNa_m(V))*1/τNa_m(V)
    du[9]=(-ĥNa+σNa_h(V))*1/τNa_h(V)
    du[10]=(-m̂K+σK_m(V))*1/τK_m(V)
    du[11]=(-m̂Ca+σCa_m(V))*1/τCa_m(V)
    du[12]=(-ĥCa+σCa_h(V))*1/τCa_h(V)

    dΨ = -γ*Ψ + Φ 
    dθ = γ*P*Ψ*(V-V̂)
    dP = α₁*P - α₂*(P*Ψ)*(P*Ψ)';        
    dP = (dP+dP')/2

    du[13:16] = dθ[:]
    du[17:32] = reshape(dP,1,16)
    du[33:36] = dΨ[:]
end