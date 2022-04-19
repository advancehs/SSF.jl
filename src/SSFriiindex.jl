
function getRII(::Type{SMAT}, PorC::Int64, pos::NamedTuple, coef::Array{Float64, 1}, 
    Ỹ::Vector, X̃::Matrix, z, Q::Matrix, u, v,  WA::Matrix,  idt::Matrix{Any})
  
    x_pre = X̃*coef[pos.begx : pos.endx]
    z_pre =   coef[pos.begz]  # mu, a scalar
    q_pre = Q*coef[pos.begq : pos.endq]
    u_pre =   coef[pos.begu] # log_σᵤ²
    v_pre =   coef[pos.begv] # log_σᵥ²
    ρ=coef[pos.begρ]
    τ=coef[pos.begτ]
  
    ϵ̃   = PorC*(Ỹ - x_pre)  
    μ   = z_pre
    h   = exp.(q_pre) 
    σᵤ² = exp(u_pre)
    σᵥ² = exp(v_pre)
    
    nofobs = length(Ỹ)
    nofyear = size(idt,1)
    
    RII = zeros(nofobs, 1)
      # bc = zeros(nofobs, 1)
      @inbounds for i in 1:nofyear
        W = WA
        @views N = idt[i,2]
        @views Mᵨ = I(N) + ρ*W
        @views Mₜ = I(N) + τ*W     
        @views Π  = Mᵨ*Mᵨ'*(σᵥ²) 
        @views yr = idt[i,1]
        @views h̃  = Mₜ*h[yr]  
        @views     RII[yr] = @.  h[yr] / h̃ 
  
        end # for_i  
      
        return RII
    end
  
  
  