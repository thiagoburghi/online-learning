const ENa = 50.
const EK = -80.
const ECa = 120.
const EL = -49.
const ESyn = -80.

σ(v,r,k) = 1/(1+exp(-(v-r)/k))
τ(v,c1,c2,r,k) = c1 + c2*exp(-(v-r)^2/k^2)

struct internal_parameters
    rNa_m::Float64
    kNa_m::Float64
    c1Na_m::Float64
    c2Na_m::Float64
    rτNa_m::Float64
    kτNa_m::Float64
    rNa_h::Float64
    kNa_h::Float64
    c1Na_h::Float64
    c2Na_h::Float64
    rτNa_h::Float64
    kτNa_h::Float64
    rK_m::Float64
    kK_m::Float64
    c1K_m::Float64
    c2K_m::Float64
    rτK_m::Float64
    kτK_m::Float64
    rCa_m::Float64        
    kCa_m::Float64  
    c1Ca_m::Float64  
    c2Ca_m::Float64  
    rτCa_m::Float64  
    kτCa_m::Float64  
    rCa_h::Float64  
    kCa_h::Float64  
    c1Ca_h::Float64
    c2Ca_h::Float64  
    rτCa_h::Float64  
    kτCa_h::Float64  
    rSyn::Float64  
    kSyn::Float64  
end

function internal_parameters(σ,rng)                         #Dethier paper      #Drion paper
    internal_parameters(-35.5*(1+σ*2*(0.5-rand(rng))),      #rNa_m      -35.5
                        5.29*(1+σ*2*(0.5-rand(rng))),       #kNa_m      5.29
                        0.0606*(1+σ*2*(0.5-rand(rng))),     #c1Na_m     1.32
                        42.307*(1+σ*2*(0.5-rand(rng))),     #c2Na_m     -1.26
                        -387.917*(1+σ*2*(0.5-rand(rng))),   #rτNa_m     -120
                        133.782*(1+σ*2*(0.5-rand(rng))),    #kτNa_m     25
                        -48.9*(1+σ*2*(0.5-rand(rng))),      #rNa_h      -48.9
                        -5.18*(1+σ*2*(0.5-rand(rng))),      #kNa_h      -5.18
                        1.5*(1+σ*2*(0.5-rand(rng))),        #c1Na_h     1.5
                        1.0*(1+σ*2*(0.5-rand(rng))),        #c2Na_h     1.0    
                        -62.9*(1+σ*2*(0.5-rand(rng))),      #rτNa_h     -62.9
                        10*(1+σ*2*(0.5-rand(rng))),         #kτNa_h     10 
                        -12.3*(1+σ*2*(0.5-rand(rng))),      #rK_m       -12.3
                        11.8*(1+σ*2*(0.5-rand(rng))),       #kK_m       11.8
                        0.799*(1+σ*2*(0.5-rand(rng))),      #c1K_m      7.2
                        5.850*(1+σ*2*(0.5-rand(rng))),      #c2K_m      -6.4
                        -76.624*(1+σ*2*(0.5-rand(rng))),    #rτK_m      -28.3
                        61.420*(1+σ*2*(0.5-rand(rng))),     #kτK_m      19.2
                        -67.1*(1+σ*2*(0.5-rand(rng))),      #rCa_m      -67.1   !!  -21.7
                        7.2*(1+σ*2*(0.5-rand(rng))),        #kCa_m      7.2     !
                        1.007*(1+σ*2*(0.5-rand(rng))),      #c1Ca_m     43.4    !!   cut in half
                        39.023*(1+σ*2*(0.5-rand(rng))),     #c2Ca_m     -42.6   !!   cut in half 
                        -117.586*(1+σ*2*(0.5-rand(rng))),   #rτCa_m     -68.1   !
                        62.873*(1+σ*2*(0.5-rand(rng))),     #kτCa_m     20.5    !
                        -82.1*(1+σ*2*(0.5-rand(rng))),      #rCa_h      -82.1   !!  -32.1
                        -5.5*(1+σ*2*(0.5-rand(rng))),       #kCa_h      -5.5    !
                        40.492*(1+σ*2*(0.5-rand(rng))),     #c1Ca_h     140     !!  105
                        86.018*(1+σ*2*(0.5-rand(rng))),     #c2Ca_h     -100    !!  -89.8 
                        -92.479*(1+σ*2*(0.5-rand(rng))),    #rτCa_h     -55     !
                        -50.237*(1+σ*2*(0.5-rand(rng))),    #kτCa_h     16.9    !     
                        -45*(1+σ*2*(0.5-rand(rng))),        #rSyn       -45     !!  -25
                        2*(1+σ*2*(0.5-rand(rng))))          #kSyn       2       !!  5
end

### Na-current
# m activation
σNa_m(V,η::internal_parameters) = σ(V,η.rNa_m,η.kNa_m)
# m time constant
τNa_m(V,η::internal_parameters) = τ(V,η.c1Na_m,η.c2Na_m,η.rτNa_m,η.kτNa_m)
# h activation
σNa_h(V,η::internal_parameters) = σ(V,η.rNa_h,η.kNa_h)
# h time constant
τNa_h(V,η::internal_parameters) = τ(V,η.c1Na_h,η.c2Na_h,η.rτNa_h,η.kτNa_h) 

### K-current
# mK activation
σK_m(V,η::internal_parameters) = σ(V,η.rK_m,η.kK_m)
# mK time constant
τK_m(V,η::internal_parameters) = τ(V,η.c1K_m,η.c2K_m,η.rτK_m,η.kτK_m)

### Ca (T) current
# mCa activation
σCa_m(V,η::internal_parameters)= σ(V,η.rCa_m,η.kCa_m)
# mCa time constant
τCa_m(V,η::internal_parameters) = τ(V,η.c1Ca_m,η.c2Ca_m,η.rτCa_m,η.kτCa_m)
# hCa activation
σCa_h(V,η::internal_parameters) = σ(V,η.rCa_h,η.kCa_h)
# hCa time constant
τCa_h(V,η::internal_parameters) = τ(V,η.c1Ca_h,η.c2Ca_h,η.rτCa_h,η.kτCa_h)

### Synaptic current
σSyn(V,η::internal_parameters) = σ(V,η.rSyn,η.kSyn)

# SS initial conditions
v_ss =  [-57.6844334480353
-74.06392322770036]

 w_ss =  [0.01481590167511684       
 0.0006853489256503658     
 0.8480836258584368        
 0.9922224212148678        
 0.019665897631623797      
 0.005567067762855896      
 0.7584993338262086        
 0.29451964455198876       
 0.04315740291951285
 0.2359656020849083
 0.0012551175666444348
 0.022645290336750494]

# using LsqFit
# v = -80:0.1:-20
# y = τNa_h.(v,Ref(η[1]))
# #τNa_m, τNa_h, τK_m, τCa_m, τCa_h
# model(v, p) = p[1] .+ p[2]*exp.(-(v .- p[3]).^2/p[4]^2)
# p0 = [0.1,1.0,-40.0,20.0]
# fit = curve_fit(model, v, y, p0)
# param = fit.param

# plot(v,model(v,param))
# plot!(v,y)