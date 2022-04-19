
function LL_T(::Type{SMAT}, Ỹ::Union{Vector,Matrix}, X̃::Matrix, z, Q::Matrix, u, v, ρ, τ, WA::Matrix, 
  PorC::Int64, nobs::Int64, po::NamedTuple, rho, idt::Matrix{Any}, ::Nothing)

  β  = rho[1:po.endx]
  δ1 = rho[po.begz]
  ω  = rho[po.begq : po.endq]
  δ2 = rho[po.begu]  
  γ  = rho[po.begv]  # may use rho[po.begw : po.endw][1]
  ρ  = rho[po.begρ]
  τ  = rho[po.begτ]

  μ   = δ1
  σᵤ² = exp(δ2)     
  σᵥ² = exp(γ)   
  h   = exp.(Q*ω) 
  ϵ̃  = PorC*(Ỹ - X̃ * β)

  nofyear = size(idt,1)


  @floop begin
    lik = zero(eltype(ϵ̃))

    @inbounds for i in 1:nofyear
        W = WA
        @views N = idt[i,2]
        @views Mᵨ = I(N) + ρ*W
        @views Mₜ = I(N) + τ*W     

        @views Π  = Mᵨ*Mᵨ'*(σᵥ²) 


        @views yr = idt[i,1]

        @views h̃  = Mₜ*h[yr]  
        @views σₛₛ² = 1.0/(h̃'*inv(Π)*h̃  + 1/σᵤ²)   
        @views μₛₛ  = (μ/σᵤ² - ϵ̃[yr]'*inv(Π)*h̃) * σₛₛ²   
        @views es2 = -0.5*(ϵ̃[yr]'*inv(Π)*ϵ̃[yr] ) 

        @views KK   = -0.5*(idt[i,2])*log(2π) - 0.5*log(det(Π))   

        lik += (KK + es2 
                + 0.5*((μₛₛ^2)/σₛₛ² - μ^2/σᵤ²) 
                + 0.5*log(σₛₛ²) + normlogcdf(μₛₛ/sqrt(σₛₛ²)) 
                - 0.5*δ2 - normlogcdf(μ/sqrt(σᵤ²)) )
    end # for i=1:nofyear
  end

  # if ρ==1 #|| τ<0
  #    lik-=Inf
  # end

  return -lik
end

function LL_T(::Type{SMAH}, Ỹ::Union{Vector,Matrix}, X̃::Matrix, z, Q::Matrix, u, v, ρ, τ, WA::Matrix, 
PorC::Int64, nobs::Int64, po::NamedTuple, rho, idt::Matrix{Any}, ::Nothing)

β  = rho[1:po.endx]
# δ1 = rho[po.begz]
ω  = rho[po.begq : po.endq]
δ2 = rho[po.begu]  
γ  = rho[po.begv]  # may use rho[po.begw : po.endw][1]
ρ  = rho[po.begρ]
τ  = rho[po.begτ]

μ   = 0
σᵤ² = exp(δ2)     
σᵥ² = exp(γ)   
h   = exp.(Q*ω) 
ϵ̃  = PorC*(Ỹ - X̃ * β)

nofyear = size(idt,1)


@floop begin
  lik = zero(eltype(ϵ̃))

  @inbounds for i in 1:nofyear
      W = WA
      @views N = idt[i,2]
      @views Mᵨ = I(N) + ρ*W
      @views Mₜ = I(N) + τ*W     

      @views Π  = Mᵨ*Mᵨ'*(σᵥ²) 


      @views yr = idt[i,1]

      @views h̃  = Mₜ*h[yr]  
      @views σₛₛ² = 1.0/(h̃'*inv(Π)*h̃  + 1/σᵤ²)   
      @views μₛₛ  = (μ/σᵤ² - ϵ̃[yr]'*inv(Π)*h̃) * σₛₛ²   
      @views es2 = -0.5*(ϵ̃[yr]'*inv(Π)*ϵ̃[yr] ) 

      @views KK   = -0.5*(idt[i,2])*log(2π) - 0.5*log(det(Π))   

      lik += (KK + es2 
              + 0.5*((μₛₛ^2)/σₛₛ² - μ^2/σᵤ²) 
              + 0.5*log(σₛₛ²) + normlogcdf(μₛₛ/sqrt(σₛₛ²)) 
              - 0.5*δ2 - normlogcdf(μ/sqrt(σᵤ²)) )
  end # for i=1:nofyear
end

# if ρ==1 #|| τ<0
#    lik-=Inf
# end

return -lik
end


  function LL_T(::Type{SART}, Ỹ::Union{Vector,Matrix}, X̃::Matrix, z, Q::Matrix, u, v, ρ, τ, WA::Matrix, 
    PorC::Int64, nobs::Int64, po::NamedTuple, rho, idt::Matrix{Any}, ::Nothing)

    β  = rho[1:po.endx]
    δ1 = rho[po.begz]
    ω  = rho[po.begq : po.endq]
    δ2 = rho[po.begu]  
    γ  = rho[po.begv]  # may use rho[po.begw : po.endw][1]
    ρ  = rho[po.begρ]
    τ  = rho[po.begτ]

    μ   = δ1
    σᵤ² = exp(δ2)     
    σᵥ² = exp(γ)   
    h   = exp.(Q*ω) 
    ϵ̃  = PorC*(Ỹ - X̃ * β)

    nofyear = size(idt,1)


    @floop begin
      lik = zero(eltype(ϵ̃))

      @inbounds for i in 1:nofyear
          W = WA
          @views N = idt[i,2]
          @views Mᵨ = inv(I(N) - ρ*W)
          @views Mₜ = inv(I(N) - τ*W)     

          @views Π  = Mᵨ*Mᵨ'*(σᵥ²) 
  

          @views yr = idt[i,1]
  
          @views h̃  = Mₜ*h[yr]  
          @views σₛₛ² = 1.0/(h̃'*inv(Π)*h̃  + 1/σᵤ²)   
          @views μₛₛ  = (μ/σᵤ² - ϵ̃[yr]'*inv(Π)*h̃) * σₛₛ²   
          @views es2 = -0.5*(ϵ̃[yr]'*inv(Π)*ϵ̃[yr] ) 

          @views KK   = -0.5*(idt[i,2])*log(2π) - 0.5*log(det(Π))   
  
          lik += (KK + es2 
                  + 0.5*((μₛₛ^2)/σₛₛ² - μ^2/σᵤ²) 
                  + 0.5*log(σₛₛ²) + normlogcdf(μₛₛ/sqrt(σₛₛ²)) 
                  - 0.5*δ2 - normlogcdf(μ/sqrt(σᵤ²)) )
      end # for i=1:nofyear
    end

    # if ρ==1 #|| τ<0
    #    lik-=Inf
    # end

    return -lik
end

function LL_T(::Type{SARH}, Ỹ::Union{Vector,Matrix}, X̃::Matrix, z, Q::Matrix, u, v, ρ, τ, WA::Matrix, 
  PorC::Int64, nobs::Int64, po::NamedTuple, rho, idt::Matrix{Any}, ::Nothing)

  β  = rho[1:po.endx]
  # δ1 = rho[po.begz]
  ω  = rho[po.begq : po.endq]
  δ2 = rho[po.begu]  
  γ  = rho[po.begv]  # may use rho[po.begw : po.endw][1]
  ρ  = rho[po.begρ]
  τ  = rho[po.begτ]

  μ   = 0
  σᵤ² = exp(δ2)     
  σᵥ² = exp(γ)   
  h   = exp.(Q*ω) 
  ϵ̃  = PorC*(Ỹ - X̃ * β)

  nofyear = size(idt,1)


  @floop begin
    lik = zero(eltype(ϵ̃))

    @inbounds for i in 1:nofyear
        W = WA
        @views N = idt[i,2]
        @views Mᵨ = inv(I(N) - ρ*W)
        @views Mₜ = inv(I(N) - τ*W)     

        @views Π  = Mᵨ*Mᵨ'*(σᵥ²) 


        @views yr = idt[i,1]

        @views h̃  = Mₜ*h[yr]  
        @views σₛₛ² = 1.0/(h̃'*inv(Π)*h̃  + 1/σᵤ²)   
        @views μₛₛ  = (μ/σᵤ² - ϵ̃[yr]'*inv(Π)*h̃) * σₛₛ²   
        @views es2 = -0.5*(ϵ̃[yr]'*inv(Π)*ϵ̃[yr] ) 

        @views KK   = -0.5*(idt[i,2])*log(2π) - 0.5*log(det(Π))   

        lik += (KK + es2 
                + 0.5*((μₛₛ^2)/σₛₛ² - μ^2/σᵤ²) 
                + 0.5*log(σₛₛ²) + normlogcdf(μₛₛ/sqrt(σₛₛ²)) 
                - 0.5*δ2 - normlogcdf(μ/sqrt(σᵤ²)) )
    end # for i=1:nofyear
  end

  # if ρ==1 #|| τ<0
  #    lik-=Inf
  # end

  return -lik
end



