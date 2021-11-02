# Sodium activation
function gating_m(v)
    Vhalf = -40.;
    k = 9.;              #15 in izhikevich
    Vmax = -38.;
    std = 30.;
    Camp = 0.46;
    Cbase = 0.04;
    τ = Cbase + Camp*exp(-(v-Vmax).^2/std^2);
    σ = 1 ./ (1+exp(-(v-Vhalf)/k));
    return τ, σ
end 

# Sodium inactivation
function gating_h(v)
    Vhalf = -62.;
    k = -7.;
    Vmax = -67.;
    std = 20.;
    Camp = 7.4;
    Cbase = 1.2;
    τ = Cbase + Camp*exp(-(v-Vmax).^2/std^2);
    σ = 1 ./ (1+exp(-(v-Vhalf)/k));
    return τ, σ
end

# Potassium activation
function gating_n(v)
    Vhalf = -53.;
    k = 15.;
    Vmax = -79.;
    std = 50.;
    Camp = 4.7;
    Cbase = 1.1;
    τ = Cbase + Camp*exp(-(v-Vmax).^2/std^2);
    σ = 1 ./ (1+exp(-(v-Vhalf)/k));
    return τ, σ
end

function HH_ode!(dz,z,p,t)
    Iapp =          p[1]
    c =             p[2]
    (gNa,gK,gL) =   p[3]
    (ENa,EK,EL) =   p[4]

    v = z[1]
    m = z[2]
    h = z[3]
    n = z[4]

    (τm,σm) = gating_m(v);
    (τh,σh) = gating_h(v);
    (τn,σn) = gating_n(v);

    g = [gNa; gK; gL];
    phi = [-m^3*h*(v-ENa);-n^4*(v-EK);-(v-EL)];

    dv = 1/c * (dot(phi,g) + Iapp(t));
    dm = 1/τm*(-m + σm);
    dh = 1/τh*(-h + σh);
    dn = 1/τn*(-n + σn);

    dz[:] = [dv,dm,dh,dn]
end

function HH_observer!(dz,z,p,t)
    Iapp =          p[1]
    c =             p[2]
    (gNa,gK,gL) =   p[3]
    (ENa,EK,EL) =   p[4]

    # True system
    v = z[1]
    m = z[2]
    h = z[3]
    n = z[4]

    (τm,σm) = gating_m(v);
    (τh,σh) = gating_h(v);
    (τn,σn) = gating_n(v);

    θ = 1/c*[gNa gK gL gNa*ENa gK*EK gL*EL 1]
    ϕ = [-m^3*h*v ...
         -n^4*v ... 
         -v ...
         m^3*h ...
         n^4 ...
         1 ...
         Iapp(t)];

    dv = dot(ϕ,θ)
    dm = 1/τm*(-m + σm);
    dh = 1/τh*(-h + σh);
    dn = 1/τn*(-n + σn);

    # Adaptive observer
    v̂ = z[5]
    m̂ = z[6]
    ĥ = z[7]
    n̂ = z[8]
    θ̂ = z[9:15]
    P = reshape(z[15+1:15+49],7,7);    
    P = (P+P')/2
    Ψ = z[15+49+1:15+49+7]

    (τm̂,σm̂) = gating_m(v);
    (τĥ,σĥ) = gating_h(v);
    (τn̂,σn̂) = gating_n(v);

    ϕ̂ = [-m̂^3*ĥ*v ...
         -n̂^4*v ... 
         -v ...
         m̂^3*ĥ ...
         n̂^4 ...
         1 ...
         Iapp(t)];

    dv̂ = dot(ϕ̂,θ̂) + γ*(1+Ψ'*P*Ψ)*(v-v̂)
    dm̂ = 1/τm̂*(-m̂ + σm̂);
    dĥ = 1/τĥ*(-ĥ + σĥ);
    dn̂ = 1/τn̂*(-n̂ + σn̂);

    dθ̂ = γ*P*Ψ*(v-v̂);
    dΨ = -γ*Ψ + ϕ̂; 
    dP = α*P - ((P*Ψ)*(P*Ψ)');
    dP = (dP+dP')/2;

    # Besancon 2000 observer
    N = 15+49+7;

    ṽ = z[N+1]
    m̃ = z[N+2]
    h̃ = z[N+3]
    ñ = z[N+4]
    θ̃ = z[N+5:N+11]

    ϕ̃ = [-m̃^3*h̃*v ...
         -ñ^4*v ... 
         -v ...
         m̃^3*h̃ ...
         ñ^4 ...
         1 ...
         Iapp(t)];
            
    (τm̃,σm̃) = gating_m(v);
    (τh̃,σh̃) = gating_h(v);
    (τñ,σñ) = gating_n(v);
    
    dṽ = dot(ϕ̃,θ̃) + γ*(v-ṽ);
    dm̃ = 1/τm̃*(-m̃ + σm̃);
    dh̃ = 1/τh̃*(-h̃ + σh̃);
    dñ = 1/τñ*(-ñ + σñ);
    dθ̃ = γ*ϕ̃*(v-ṽ);

    dz[:] = [dv;dm;dh;dn;dv̂;dm̂;dĥ;dn̂;dθ̂;dP[:];dΨ;dṽ;dm̃;dh̃;dñ;dθ̃]';

end

### REMOVE!
# function HH_loss_ode!(dz,z,p,t)
#     # The only parameter here is ĝNa
#     # This makes it simple to compute the gradient of the loss function
#     # Fixed parameters
#     c =             1.
#     (gNa,gK,gL) =   (120.,36.,0.3)
#     (ENa,EK,EL) =   (55.,-77.,-54.4)

#     # Variable parameter
#     ĝNa = p

#     v = z[1]
#     m = z[2]
#     h = z[3]
#     n = z[4]
#     v̂ = z[5]
#     m̂ = z[6]
#     ĥ = z[7]
#     n̂ = z[8]

#     # True model
#     (τm,σm) = gating_m(v);
#     (τh,σh) = gating_h(v);
#     (τn,σn) = gating_n(v);

#     μ = [gNa; gK; gL];
#     phi = [-m^3*h*(v-ENa);-n^4*(v-EK);-(v-EL)];

#     dv = 1/c * (dot(phi,g) + Iapp(t));
#     dm = 1/τm*(-m + σm);
#     dh = 1/τh*(-h + σh);
#     dn = 1/τn*(-n + σn);

#     # Estimator model
#     (τm̂,σm̂) = gating_m(v̂);
#     (τĥ,σĥ) = gating_h(v̂);
#     (τn̂,σn̂) = gating_n(v̂);

#     ĝ = [ĝNa; gK; gL];
#     pĥi = [-m̂^3*ĥ*(v̂-ENa);-n̂^4*(v̂-EK);-(v̂-EL)];
#     dv̂ = 1/c * (dot(pĥi,ĝ) + Iapp(t));
#     dm̂ = 1/τm̂*(-m̂ + σm̂);
#     dĥ = 1/τĥ*(-ĥ + σĥ);
#     dn̂ = 1/τn̂*(-n̂ + σn̂);

#     dz[:] = [dv,dm,dh,dn,dv̂,dm̂,dĥ,dn̂]
# end

