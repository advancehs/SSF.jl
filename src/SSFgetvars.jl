
#* ---- check multi-collinearity in the data --------

function checkCollinear(modelid, xvar, zvar, qvar, wvar, vvar)                  

    # check multicollinearity in the dataset (matrix)
 
     for i in 1:5
 
         if i == 1 
             themat = xvar
             eqname = " `frontier` " 
             eqsym  = :frontier
         elseif i == 2
             themat = zvar
             eqname = " `μ` "
             eqsym  = :μ
         elseif i == 3
             themat = qvar

                 eqname = " `hscale` "    
                 eqsym  = :hscale
                 
         elseif i == 4
             themat = wvar        
             eqname = " `σᵤ²` "   
             eqsym  =  :σᵤ²
         elseif i == 5
             themat = vvar            
             eqname = " `σᵥ²` "
             eqsym  = :σᵥ²
         end
 
 
         if length(themat) > 0 && size(themat, 2) > 1  # check only if have more than 1 var
 
             _, pivots = rref_with_pivots(themat)  # pivots has unique,
                                                   # non-collinear columns;
                                                   # require `using RowEchelon`
                                                 
             if length(pivots) != size(themat, 2) # number of columns are different, meaning collinear is dropped
                 theColl = filter(x -> x ∉ pivots, 1:size(themat, 2)) # positions of problematic columns
                 bb = ""
                 for i in theColl  # get problematic var's name from the equation
                     bb = bb*String(_dicM[eqsym][i])*String(", ")
                 end
                 printstyled(bb, "appear to have perfect collinearity with other variables in the equation", eqname, ". Try dropping this/these variable(s) from the equation.\n"; color = :red)
                 throw("Multi-collinearity in the variables.")
             end  # if length(pivots)       
         end # if size(themat, 2) > 1
     end # for i in 1:5
   end
 



########################################################
####               check constant                   ####
########################################################


macro requireConst(arg)
    return arg
end

    # ---- check overall constant --------------

function checkConst(mat::Array, S::Symbol, requireConst::Int)
    for i in 1:size(mat,2)
        aa = mat[:,i]
        # varstd = sqrt(sum((aa .- sum(aa)/length(aa)).^2)/length(aa)) ; then varstd <= 1e-8
        bb = unique(aa)
        if requireConst == 1 
            (length(bb) == 1) || throw("The variable in the $(S) function has to be a constant for this model.")
        else
            (length(bb)  > 1) || throw("The $(S) function cannot have a constant for this model.")
        end
    end
end



function get_rowIDT(tvar) 

    # Create a matrix of panel information.

    T  = length(unique(tvar)) # T=number of panels
    id = Array{Any}(undef,T,2)
    id[:,1] = unique(tvar)    # list of year with no repetition
    @inbounds for i = 1:T
        @views id[i,2]= sum(tvar.== id[i,1])    # number of periods for the panel
    end
    @views Nᵢ = id[:,2]

    rowID =  Vector{Vector}()  
    @inbounds for i=1:T
        @views cc = findall(x-> x == id[i,1], tvar) # row index of i'th year
        push!(rowID, cc) # put the year info in the vector; faster than using UnitRange
    end    
  
    rowIDT = hcat(rowID, Nᵢ) 

    return rowIDT # (Nx2): col_1 is panel's row info; col_2 is panel's number of periods
end











  #?--------- panel FE SWA, truncated normal ----------------
  
  function getvar(::Type{SMAT}, dat::DataFrame)
  
    yvar = dat[:, _dicM[:depvar]]   # still a DataFrame
    xvar = dat[:, _dicM[:frontier]]
    zvar = dat[:, _dicM[:μ]]
    qvar = dat[:, _dicM[:hscale]]  
    uvar = dat[:, _dicM[:σᵤ²]]
    vvar = dat[:, _dicM[:σᵥ²]]
    W    =   _dicM[:W] 
    ρ    =   _dicM[:ρ] 
    τ    =   _dicM[:τ] 

    tvar = dat[:, _dicM[:timevar]]
    ivar = dat[:, _dicM[:idvar]] 

    #* --- model info printout --------- 
    modelinfo1 = "true fixed effect model of Wang and Ho (2010 JE), normal and truncated-normal"
    modelinfo2 = begin
     """
     * In the case of type(cost), "- ũᵢₜ" below is changed to "+ ũᵢₜ".
  
     $(_dicM[:depvar][1]) = frontier($(_dicM[:frontier])) + ṽᵢₜ - ũᵢₜ,
     
     where ṽᵢₜ = vᵢₜ + ρ * Wᵢ * vₜ
           vᵢₜ ∼ N(0, σᵥ²),
                σᵥ² = exp(log_σᵥ²) 
                     = exp($(_dicM[:σᵥ²]));
           ũᵢₜ = uᵢₜ + τ * Wᵢ * uₜ
                uᵢₜ = hscaleᵢₜ * uₜˢ,
                  hscaleᵢₜ = exp($(_dicM[:hscale])),
                  uₜˢ ∼ N⁺(μ, σᵤ²),
                     μ = $(_dicM[:μ])
                     σᵤ² = exp(log_σᵤ²) 
                         = exp($(_dicM[:σᵤ²]));
     """
    end
  
    #* --- retrieve and generate important parameters -----
  
    #*   number of obs and number of variables
    nofx = nofz = nofq = nofu = nofv = 0  # to make a complete list
  
    nofobs  = nrow(dat)    
    nofx    = size(xvar,2)  # nofx: number of x vars
    nofz    = size(zvar,2)
    nofq    = size(qvar,2)  
    nofu    = size(uvar,2)
    nofv    = size(vvar,2)
    nofρ    = length([ρ])
    nofτ    = length([τ])
    nofpara = nofx + nofz + nofq + nofu + nofv + nofρ + nofτ
  
    nofvar = (nofobs=nofobs, nofx=nofx, nofz=nofz, nofq=nofq,
              nofu=nofu, nofv=nofv, nofρ=nofρ,nofτ=nofτ, nofpara=nofpara, nofmarg = nofz+nofq+nofu)
  
    #* positions of the variables/parameters
    begx=endx=begz=endz=begq=endq=begu=endu=begv=endv=begρ=endρ=begτ=endτ=0
  
    begx = 1
    endx = nofx
    begz = endx + 1
    endz = begz + nofz-1
    begq = endz + 1
    endq = begq + nofq-1
    begu = endq + 1
    endu = begu + nofu-1
    begv = endu + 1
    endv = begv + nofv-1
    begρ = endv + 1
    endρ = begρ + nofρ-1
    begτ = endρ + 1
    endτ = begτ + nofτ-1

    posvec = (begx=begx, endx=endx, begz=begz, endz=endz,
              begq=begq, endq=endq, begu=begu, endu=endu,
              begv=begv, endv=endv,begρ=begρ, endρ=endρ,
              begτ=begτ, endτ=endτ)
  
    #* create equation names and mark positions for making tables
    eqvec = (frontier = begx + 1, 
                    μ = begz + 1,
           log_hscale = begq + 1,
              log_σᵤ² = begu + 1,
              log_σᵥ² = begv + 1,
                    ρ = begρ + 1,
                    τ = begτ + 1)
  
    #* create equation names and mark positions 
    eqvec2 = (coeff_frontier = (begx:endx), 
                     coeff_μ = (begz:endz),
            coeff_log_hscale = (begq:endq),
               coeff_log_σᵤ² = (begu:endu),
               coeff_log_σᵥ² = (begv:endv),
                     coeff_ρ = (begρ:endρ),
                     coeff_τ = (begτ:endτ))             
  
    #* retrieve variable names for making tables
    xnames  = names(xvar)
    znames  = names(zvar)
    qnames  = names(qvar)
    unames  = names(uvar)
    vnames  = names(vvar)
    ρnames  = "ρ"
    τnames  = "τ"

    varlist = vcat(" ", xnames, znames, qnames, unames, vnames,ρnames,τnames)
   
    #* Converting the dataframe to matrix in order to do computation
    yvar  = convert(Array{Float64}, Matrix(yvar))
    xvar  = convert(Array{Float64}, Matrix(xvar))
    zvar  = convert(Array{Float64}, Matrix(zvar))
    qvar  = convert(Array{Float64}, Matrix(qvar))
    uvar  = convert(Array{Float64}, Matrix(uvar))
    vvar  = convert(Array{Float64}, Matrix(vvar))
  
    tvar  = convert(Array{Float64}, Matrix(tvar))
    ivar  = convert(Array{Float64}, Matrix(ivar))
  
  
  
    # * various functions can and cannot contain a constant, check! ---- *#
    checkConst(xvar, :frontier, @requireConst(0))
    checkConst(zvar, :μ,        @requireConst(1))
    checkConst(qvar, :hscale,   @requireConst(0)) 
    checkConst(uvar, :σᵤ²,      @requireConst(1))
    checkConst(vvar, :σᵥ²,      @requireConst(1))
  
    #* matrix  transformation
  
  
    rowIDT = get_rowIDT(vec(tvar))   # rowIDT (Nx2): col_1 is panel's row info; col_2 is panel's number of periods
   
    yxdata = hcat(yvar, xvar) 
         D = zeros(nofobs, 1+nofx)      # pre-allocate the transformed dataset
         T = length(unique(vec(tvar))) # N=number of panels
  
  
    for i=1:T
        @views D[rowIDT[i,1], :] = (yxdata[rowIDT[i,1], :]) # INVtrM[i] * yxdata[rowIDT[i,1], :] # transform the data
    end    
    
    ỹ = (D[:, 1]) 
    x̃ = (D[:, 2:end]) 
  
    return modelinfo1, modelinfo2, posvec, nofvar, eqvec, eqvec2, ỹ, x̃, zvar, qvar, uvar, vvar, ρ, τ, W,    rowIDT, varlist
  
  
  end
  

