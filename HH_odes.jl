function ς(x,lb,ub)
    return min(ub,max(lb,x))
end

function τₘ(v)
    Vmax = -38.;
    std = 30.;
    Camp = 0.46;
    Cbase = 0.04;
    return Cbase + Camp*exp(-(v-Vmax).^2/std^2);
end

function τₕ(v)
    Vmax = -67.;
    std = 20.;
    Camp = 7.4;
    Cbase = 1.2;
    return Cbase + Camp*exp(-(v-Vmax).^2/std^2);
end

function τₙ(v)
    Vmax = -79.;
    std = 50.;
    Camp = 4.7;
    Cbase = 1.1;
    return Cbase + Camp*exp(-(v-Vmax).^2/std^2);  
end

function σ(v,r,k)
    return 1 ./ (1+exp(-(v-r)/k));
end

function ∂ᵨσ(v,r,k)
    return -exp(-(v-r)/k) ./ (k*(1+exp(-(v-r)/k))^2)
end
 
function HH_ODE!(dz,z,p,t)
    Iapp        =   p[1]
    θ           =   p[2]
    η           =   p[3]
    (ENa,EK,EL) =   p[4]
    α,β,γ       =   p[5]
    (km,kh,kn)  =   p[6]

    # True system
    v = z[1]
    w = z[2:4]
    m,h,n = w

    ϕ₁ =   [-m^3*h*(v-ENa) ...
            -n^4*(v-EK) ... 
            -(v-EL) ...
            Iapp(t)];

    A = [-1/τₘ(v) 0 0;
          0 -1/τₕ(v) 0;
          0 0 -1/τₙ(v)]

    b = [σ(v,η[1],km)/τₘ(v);
         σ(v,η[2],kh)/τₕ(v);
         σ(v,η[3],kn)/τₙ(v)]

    dv = dot(ϕ₁,θ)
    dw = A*w + b

    # Adaptive observer
    v̂ = z[5]
    ŵ = z[6:8]
    m̂,ĥ,n̂ = ŵ

    θ̂ = z[9:12]
    η̂ = z[13:15]

    P = reshape(z[15+1:15+49],7,7);    
    P = (P+P')/2
    ψ₁ = z[15+49+1:15+49+7]
    ψ₂ = z[15+49+7+1:15+49+7+3]
    Ψ₂ = Diagonal(ψ₂)

    ϕ̂₁ = [-m̂^3*ĥ*(v-ENa) ...
          -n̂^4*(v-EK) ... 
          -(v-EL) ...
          Iapp(t)];

    b̂ = [σ(v,η̂[1],km)/τₘ(v);
         σ(v,η̂[2],kh)/τₕ(v);
         σ(v,η̂[3],kn)/τₙ(v)]

    ∂ᵨb̂ = [∂ᵨσ(v,η̂[1],km)/τₘ(v) 0 0; 
           0 ∂ᵨσ(v,η̂[2],kh)/τₕ(v) 0; 
           0 0 ∂ᵨσ(v,η̂[3],kn)/τₙ(v)]

    ∂₂v̂̇ =  [-θ̂[1]*3*m̂^2*ĥ*(v-ENa) -θ̂[1]*m̂^3*(v-ENa) -θ̂[2]*4*n̂^3*(v-EK)]

    dv̂ = dot(ϕ̂₁,θ̂) + γ*(1+ψ₁'*P*ψ₁)*(v-v̂)
    dŵ = A*ŵ + b̂ + γ*[zeros(4,3);Ψ₂]'*P*ψ₁*(v-v̂)

    dθ̂ = γ*P*ψ₁*(v-v̂);
    dψ₁ = -γ*ψ₁ + [zeros(4,3);Ψ₂]*∂₂v̂̇' + [ϕ̂₁;zeros(3,1)]; 
    dΨ₂ = Ψ₂*A + ∂ᵨb̂;
    dψ₂ = diag(dΨ₂)

    dP = α*P + β*Matrix(I, 7, 7) - γ*((P*ψ₁)*(P*ψ₁)');
    dP = (dP+dP')/2;

    dz[:] = [dv;dw;dv̂;dŵ;dθ̂;dP[:];dψ₁;dψ₂]';
end

function f_SDE!(dz,z,p,t)
    ### Recover parameters
    Iapp             =   p[1]
    θ                =   p[2]
    η                =   p[3]
    ENa,EK,EL        =   p[4]
    α,β,γ            =   p[5]
    km,kh,kn         =   p[6]
    σnoise,dt =   p[7]    # Coloured white noise time constant
    outinj           =   p[8]            

    ### Recover states
    # True system
    v = z[1]
    w = z[2:4];
    m,h,n = w
    # Observer
    v̂ = z[5]
    ŵ = z[6:8];                     
    θ̂ = z[9:12]
    η̂ = z[13:15]
    P = reshape(z[15+1:15+49],7,7); 
    ψ₁ = z[15+49+1:15+49+7]
    ψ₂ = z[15+49+7+1:15+49+7+3];    
    m̂,ĥ,n̂ = ŵ
    P = (P+P')/2
    Ψ₂ = Diagonal(ψ₂)
    # Johnson-Nyquist noise
    noise = z[end]

    ### Define vector fields
    # True system
    ϕ₁ =   [-m^3*h*(v-ENa) ...
            -n^4*(v-EK) ... 
            -(v-EL) ...
            Iapp(t)];

    A = [-1/τₘ(v) 0 0;
          0 -1/τₕ(v) 0;
          0 0 -1/τₙ(v)]

    b = [σ(v,η[1],km)/τₘ(v);
         σ(v,η[2],kh)/τₕ(v);
         σ(v,η[3],kn)/τₙ(v)]

    dv = dot(ϕ₁,θ)
    dw = A*w + b
    
    # Observer
    y  = v + noise
    # m̂,ĥ,n̂ = ς.([m̂,ĥ,n̂],0,1)
    θ̂ₛ = ς.(θ̂,0,200)
    if outinj
        yᵥ = y
        ∂₁v̂̇ = 0.0
    else
        yᵥ = v̂
        ∂₁v̂̇ = -θ̂ₛ[1]*m̂^3*ĥ - θ̂ₛ[2]*n̂^4 - θ̂ₛ[3]
    end

    ϕ̂₁ = [-m̂^3*ĥ*(yᵥ-ENa) ...
          -n̂^4*(yᵥ-EK) ... 
          -(yᵥ-EL) ...
          Iapp(t)];

    b̂ = [σ(y,η̂[1],km)/τₘ(y);
         σ(y,η̂[2],kh)/τₕ(y);
         σ(y,η̂[3],kn)/τₙ(y)]

    ∂ᵨb̂ = [∂ᵨσ(y,η̂[1],km)/τₘ(y) 0 0; 
           0 ∂ᵨσ(y,η̂[2],kh)/τₕ(y) 0; 
           0 0 ∂ᵨσ(y,η̂[3],kn)/τₙ(y)]

    ∂₂v̂̇ = [-θ̂ₛ[1]*3*m̂^2*ĥ*(yᵥ-ENa) -θ̂ₛ[1]*m̂^3*(yᵥ-ENa) -θ̂ₛ[2]*4*n̂^3*(yᵥ-EK)]

    dv̂ = dot(ϕ̂₁,θ̂) + γ*(1+ψ₁'*P*ψ₁)*(y-v̂)
    dŵ = A*ŵ + b̂ + γ*[zeros(4,3);Ψ₂]'*P*ψ₁*(y-v̂)
    
    dθ̂ = γ*P*ψ₁*(y-v̂);

    dψ₁ = (-γ+∂₁v̂̇)*ψ₁ + [zeros(4,3);Ψ₂]*∂₂v̂̇' + [ϕ̂₁;zeros(3,1)]; 
    
    dΨ₂ = Ψ₂*A + ∂ᵨb̂;
    dψ₂ = diag(dΨ₂)
    
    dP = α*P + β*Matrix(I, 7, 7) - γ*((P*ψ₁)*(P*ψ₁)');
    dP = (dP+dP')/2;
    dP = dP[:]

    # Johnson-Nyquist noise
    dnoise = -1/dt*noise;

    dz[:] = [dv;dw;dv̂;dŵ;dθ̂;dP;dψ₁;dψ₂;dnoise];
end

function g_SDE!(dz,z,p,t)
    σnoise,dt = p[7]

    # Johnson-Nyquist noise
    gnoise = σnoise/sqrt(dt);

    dz[:] = [zeros(15+49+7+3,1);gnoise]
end