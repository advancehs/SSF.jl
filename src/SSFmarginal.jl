function nonConsDataFrame(D::DataFrame, M::Matrix)
  # Given a DataFrame containing the marginal effects 
  # of a set of exogenous determinants $(x1, x2, ..., xn)$
  # on E(u), it return the DataFrame where the marginal 
  # effect of constant $x$s are removed.

  # D: the marginal effect DataFrame; 
  # M: the matrix of (x1, .., xn) where the marginal 
  #    efect is calculated from.

 counter = 0      
 for w in collect(names(D),)
      counter += 1
      if length(unique(M[:, counter])) == 1 # is a constant
          select!(D, Not(Symbol(w)))
      end
 end 
 return D
end


function addDataFrame(Main::DataFrame, A::DataFrame)
  # Combine two DataFrame with unions of columns.
  # For same-name columns, the values are added together.

for k in collect(names(A),) # deal with the wvar
          if k ∈ names(Main)
               Main[:, Symbol(k)] = Main[:, Symbol(k)] + A[:, Symbol(k)]
          else 
               insertcols!(Main, Symbol(k) => A[:, Symbol(k)])
          end
 end 
 return Main  
end 




#? --------------------------------------------------------------------
#? - panel FE SWA, truncated normal, marginal effect function -
#? -------------------------------------------------------------------- 


function marg_swat( # PorC::Int64, 
  cc, Mₜ, pos::NamedTuple, coef::Array{Float64, 1},
  Zmarg, Qmarg, Umarg)


  z_pre = @. Zmarg*coef[pos.begz : pos.endz]  # mu, a scalar
  q_pre = @. Qmarg*coef[pos.begq : pos.endq]
  u_pre = @. Umarg*coef[pos.begu : pos.endu] # log_σᵤ²
  ρ=coef[pos.begρ]
  τ=coef[pos.begτ]

    μ   = z_pre
    h   = exp.(q_pre) 
    σᵤ = exp.(0.5*u_pre)
    h̃ = Mₜ*h
    μₛ   = @. h̃ * μ
    σᵤₛ = @. h̃* h̃ * σᵤ
    Λ  = @. μₛ/σᵤₛ 
 
    uncondU = @.σᵤₛ* (Λ + normpdf(Λ) / normcdf(Λ)) # kx1
    uncondU = uncondU[cc][1]
end   


#? -- panel FE Wang and Ho, truncated normal, , get marginal effect ----

function get_marg(::Type{SMAT}, # PorC::Int64, 
  pos::NamedTuple, num::NamedTuple, coef::Array{Float64, 1}, 
  Z::Matrix, Q::Matrix, U::Matrix,WA::Matrix,idt::Matrix{Any})

     #* Note that Y and X are within-transformed by `getvar`, 
     #* but Q, W, V are still in the original level.


  z_pre =   coef[pos.begz]  # mu, a scalar
  q_pre = Q*coef[pos.begq : pos.endq]
  u_pre =   coef[pos.begu] # log_σᵤ²

  ρ=coef[pos.begρ]
  τ=coef[pos.begτ]


  μ   = z_pre
  h   = exp.(q_pre) 
  σᵤ² = exp(u_pre)

  

  nofyear = size(idt,1)

  mm_z = Array{Float64}(undef, num.nofz, num.nofobs)  
  mm_q = Array{Float64}(undef, num.nofq, num.nofobs)                  
  mm_u = Array{Float64}(undef, num.nofu, num.nofobs)   
  W = WA  
  @inbounds for i in 1:nofyear
              @views N_ind = idt[i,2]
              @views Mᵨ = I(N_ind) + ρ*W
              @views Mₜ = I(N_ind) + τ*W     

              @views yr = idt[i,1]  # yr对应的索引
              @views h̃  = Mₜ*h[yr]  
                     
    @inbounds for index in yr 
                  cc = findall(x->x==index,yr)
                  h̃_c = h̃[cc,:]
                  @views marg = ForwardDiff.gradient(marg -> marg_swat(cc,Mₜ, pos, coef, 
                                                          marg[1 : num.nofz*N_ind],
                                                          marg[num.nofz*N_ind+1 : num.nofz*N_ind+num.nofq*N_ind],
                                                          marg[num.nofz*N_ind+num.nofq*N_ind+1 : num.nofmarg*N_ind]),
                                                          vcat( Z[yr,:], Q[yr,:], U[yr,:]) );                            

                mm_z[:,index] = marg[1 : num.nofz*N_ind][cc,:]
                mm_q[:,index] = marg[num.nofz*N_ind+1 : num.nofz*N_ind+num.nofq*N_ind][cc,:]
                mm_u[:,index] = marg[num.nofz*N_ind+num.nofq*N_ind+1 : end][cc,:]

               end  # for i in 1:num.nofobs
            end  # for yr in 1:nofyear  
  margeff = DataFrame(mm_z', _dicM[:μ]) # the base set
  mm_q = DataFrame(mm_q', _dicM[:hscale])
  mm_u = DataFrame(mm_u', _dicM[:σᵤ²])

  #* purge off the constant var's marginal effect from the DataFrame
  margeff = nonConsDataFrame(margeff, Z)
  mm_q = nonConsDataFrame(mm_q, Q)
  mm_u = nonConsDataFrame(mm_u, U)

  #* if same var in different equations, add up the marg eff
  margeff = addDataFrame(margeff, mm_q)
  margeff = addDataFrame(margeff, mm_u)

   #* prepare info for printing
   margMean = (; zip(Symbol.(names(margeff)) , round.(mean.(eachcol(margeff)); digits=5))...)

   #* modify variable names to indicate marginal effects
   newname = Symbol.(fill("marg_", (size(margeff,2), 1)) .* names(margeff))
   margeff = rename!(margeff, vec(newname))

  return  margeff, margMean
end  
