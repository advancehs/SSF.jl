########################################################
####                ssfmodel_spec()                  ####
########################################################

"""
ssfmodel_spec(<keyword arguments>)

Provide specifications of the stochastic frontier model, including the type of model
and names of variables or matrix used in estimating the model. Two ways to
specify: Method 1: Use DataFrame as input, and Method 2: use matrix as input.

# Method 1 (DataFrame input)
Variables come from a DataFrame, and column names of the Dataframe are used in
the variable input. With this method,
equations are identified by macros but not functions (e.g., `@depvar()` but
not `depvar()`).

## Arguments of Method 1
- `sfdist(::Vararg)`: the distribution assumption of the one-sided stochastic
  variable (aka inefficiency term) of the model;
  possible choices include `truncated` (or `trun`, `t`), `half` (or `h`).
- `sftype(::Vararg)`: whether the model is a `production` (or `prod`) frontier
  or a `cost` frontier.
- `ssfpanel(::Vararg)`: the type of panel model. Choices include `SMA`
  (true fixed effect model of  Orea L, Álvarez I C 2019 JE), `SAR` (true fixed
  model of Orea L, Álvarez I C 2019 JE).
- `depvar(::Vararg)`: the dependent variable from a DataFrame.
- `frontier(::Vararg)`: a list of variables, separated by commas, in the frontier function.
- `μ(::Vararg)` : a list of variable, separated by comma,
  in the linear function of μ. (`sftype(trun)` only).
- `σᵥ²(::Vararg)` or `@sigma_v_2(::Vararg)`: a list of variable, separated by comma, in the σᵥ²
  equation.
- `σᵤ²(::Vararg)` or `@sigma_u_2(::Vararg)`: a list of variable, separated by comma, in the σᵤ²
  equation.
- `timevar(::Vararg)`: the variable containing the time period information.
  Panel data model only.
- `idvar(::Vararg)`: the variable identifying each individual. Panel data
  model only.
- message::Bool: Whether printing (=true) or not (=false, the default) the
  confirmation message "A dictionary from ssfmodel_spec() is generated."
  on the screen after `ssfmodel_spec()` is successfully executed.  


# Examples
```julia-repl

ssfmodel_spec(ssfpanel(SMA), ssftype(prod), ssfdist(half),
             timevar(yr), idvar(id),
             depvar(y), 
             frontier(x1, x2, _cons), 
             μ(_cons),
             σᵤ²(_cons),
             σᵥ²(_cons),
             hscale(z1),
             ρ(0.1),
             τ(0.1),
             weight(W), 
             message = false);
```
"""

function ssfmodel_spec(arg::Vararg; message::Bool=false) 

    global _dicM 
           _dicM = Dict{Symbol, Any}()  # nullify and initiate new dictionaries when a new model is specified

         #* -- creates default values ---

        for k in ( :panel, :timevar, :idvar, :dist, :type, :depvar, :frontier, :μ, :hscale, :σᵤ², :σᵥ², :ρ, :τ, :W, :hasDF, :transfer, :misc) # take out :data, :η, :λ, :τ, in this revision
            _dicM[k] = nothing
        end
     
        #* -- replace the defaults with user's values ---          

        for d in :($(arg))
            _dicM[d[1]] = d[2]
        end 

        #* ==== Method 2, matrix input (such as in simulations), create a DataFrame

           _dicM[:hasDF]    = true
           _dicM[:transfer] = false


           if typeof(_dicM[:depvar]) != Array{Symbol,1} # not a DataFrame

            _dicM[:hasDF] = false 
 
            isa(_dicM[:depvar][1], Vector) || isa(_dicM[:depvar][1], Matrix) || throw("
            `depvar()` has to be a Vector or Matrix (e.g., Array{Float64, 1} or Array{Float64, 2}). 
            Check with `isa(your_thing, Matrix)` or `isa(your_thing, Vector)`. 
            Try `convert()`, `reshape()`, `Matrix()`, or something similar.")
 
            comDF = _dicM[:depvar][1]  # create the first data column of comDF
            varname = [:depvar]
            _dicM[:depvar] = [:depvar]
 
            for k in (:timevar, :idvar, :frontier, :μ, :hscale, :σᵤ², :σᵥ²) 
                if _dicM[k] !== nothing # if not nothing, must be Array
                   isa(_dicM[k], Vector) || isa(_dicM[k][1], Vector) || isa(_dicM[k][1], Matrix) || throw("
                      `k` has to be a Vector or Matrix (e.g., Array{Float64, 1} or Array{Float64, 2}). 
                      Check with `isa(your_thing, Matrix)` or `isa(your_thing, Vector)`. 
                      To convert, try `convert()`, `reshape()`, `Matrix()`, or something similar.")
    
                   (isa(_dicM[k], Vector)  && length(_dicM[k][1]) == 1) ?  _dicM[k] = [_dicM[k]] : nothing # ugly fix for pure vector input
    
                   @views comDF = hcat(comDF, _dicM[k][1]) # combine the data
                   aa = Symbol[]
                   for i in 1:size(_dicM[k][1], 2)
                       push!(aa, Symbol(String(k)*"_var$(i)")) # create name for the data
                   end                  
                   varname = vcat(varname, aa) # combine the dataname
                   _dicM[k] = aa
                end 
            end # for k in (...)
 
 
            comDF = DataFrame(comDF, varname)
            _dicM[:sdf] = comDF
 
         end # if typeof(...)
 
      

        #* -- check the model identifier and the model type ---

        (_dicM[:dist] !== nothing) || throw("You need to specify dist().")
        (_dicM[:type] !== nothing) || throw("You need to specify type().")

        #* --- get the model identifier -------

        s = uppercase(String(_dicM[:dist][1])[1:1])
        
        global tagD
        if (_dicM[:panel] == [:SMA]) && (s=="T")
            tagD = Dict{Symbol, Type{SMAT}}()
            tagD[:modelid] = SMAT 
        elseif (_dicM[:panel] == [:SMA]) && (s == "H")  # panel and dist=half
            tagD = Dict{Symbol, Type{SMAH}}()
            tagD[:modelid] = SMAH 
        elseif (_dicM[:panel] == [:SMA])
            throw("The panel SMA model can only have `sfdist(trun)` or `sfdist(half)`.")
        elseif (_dicM[:panel] == [:SAR]) && (s=="T")
            tagD = Dict{Symbol, Type{SART}}()
            tagD[:modelid] = SART 
        elseif (_dicM[:panel] == [:SAR]) && (s=="H")
            tagD = Dict{Symbol, Type{SARH}}()
            tagD[:modelid] = SARH 
        elseif (_dicM[:panel] == [:SAR]) 
            throw("The panel SAR model can only have `sfdist(half)` or `sfdist(half)`.")
        end

        #* ---- check if the model has the correct syntax ---

        # SFrontiers.checksyn(tagD[:modelid])

        #* ----- make return ----------- 
        if message 
          printstyled("A dictionary from ssfmodel_spec() is generated.\n"; color = :green)  
        end  
        return _dicM # for debugging purpose

end  # end of ssfmodel_spec()



########################################################
####                ssfmodel_init()                  ####
########################################################
"""
    ssfmodel_init(<keyword arguments>)

Provide initial values for the stochastic frontier model estimation. The
values could be a vector or scalars. It creates a global dictionary `_dicINI`. Optional.

# Arguments
- `frontier(::Union{Vector, Real})`: initial values of parameters in
  the `frontier()` function
- `μ(::Union{Vector, Real})` : initial values of parameters in the `μ` function
- `hscale(::Union{Vector, Real})`: initial values of parameters in the `hscale()` function
- `σᵤ²(::Union{Vector, Real})` or `sigma_u_2(::Union{Vector, Real})`: initial values of parameters in the
   `σᵤ²` function
- `σᵥ²(::Union{Vector, Real})` or `sigma_v_2(::Union{Vector, Real})`: initial values of parameters in the
   `σᵥ²` function    
-  `ρ(::Real)`: initial values of parameters in the `ρ()` function
-  `τ(::Real)`: initial values of parameters in the `τ()` function
- message::Bool: Whether printing (=true) or not (=false, the default) the
  confirmation message "A dictionary from ssfmodel_init() is generated."
  on the screen after `ssfmodel_init()` is successfully executed.

# Remarks
- Equations do not have to follow specific orders.
- `ssfmodel_init(...)` is optional but is highly recommended. If it is not
  specified or is specified as an empty set, default values are used.
- It is not necessary to specify a complete set of equations. A partial list 
  or even empty lists are acceptable. Default values will be substituted for the
  missing equations.
- The generated `_dicINI` is inheritable in the sense that an exiting
  `_dicINI` (from the previous run of the same or a different model, for
  example) will be used if the current model does not have its own
  `ssfmodel_init(...)`. This design has advantages in a simulations study where
  `ssfmodel_init(...)` needs to be specified only once.

# Examples
```julia-repl
b_ini = ones(2)*0.2
ssfmodel_init( # frontier(bb),             # may skip and use default
             μ(b_ini), 
             hscale(0.1),              
             σᵤ²([-0.1, -0.1]),    # may use a vector
             σᵥ²(-0.1) ,
             ρ(-0.1),
             τ(-0.1) )                   
        
```
"""

function ssfmodel_init(arg::Vararg; message::Bool =false) # create a dictionary of inital vectors

   global _dicINI
          _dicINI = Dict{Symbol, Any}()

    for d in :($(arg))
        _dicINI[d[1]] = d[2]
    end        

    #* If has the key, creates the alias key with the same value.
    
    !(haskey(_dicINI, :μ))      || (_dicINI[:eqz] = _dicINI[:μ])
    !(haskey(_dicINI, :hscale)) || (_dicINI[:eqq] = _dicINI[:hscale]) 
    !(haskey(_dicINI, :σᵤ²))    || (_dicINI[:equ] = _dicINI[:σᵤ²])
    !(haskey(_dicINI, :σᵥ²))    || (_dicINI[:eqv] = _dicINI[:σᵥ²])
    !(haskey(_dicINI, :ρ))    || (_dicINI[:eqρ] = _dicINI[:ρ])
    !(haskey(_dicINI, :τ))    || (_dicINI[:eqτ] = _dicINI[:τ])
    if message 
      printstyled("A dictionary from ssfmodel_init() is generated.\n"; color = :green) 
    end
    return _dicINI # for debugging purpose
end    

########################################################
####                ssfmodel_opt()                   ####
########################################################
"""
    ssfmodel_opt(<keyword arguments>)

Provide options to the optimization algorithms for the maiximum likelihood
estimation. It creates a global dictionary `_dicOPT`. Optional. The `Optim`
package is used for the optimization, and a subset of
`Optim`'s keywords are directly accessible from this API. 

# Arguments
- `warmstart_solver(algorithm)`: The algorithm used in the first-stage ("warmstart")
  optimization process, which serves the purpose of improving upon the initial
  values for the second-stage ("main") estimation. The default is
  `NelderMead()`. Others include `SimulatedAnnealing()`, `SAMIN()`, `ParticleSwarm()`,
  `ConjugateGradient()`, `GradientDescent()`, `BFGS()`, `LBFGS()`,
  `Newton()`, `NewtonTrustRegion()`, and `IPNewton()`. See
  http://julianlsolvers.github.io/Optim.jl/stable/ for details.
  Non-gradient based algorithms are recommended for the warmstart solver. 
- `warmstart_maxIT(::Int64)`: The iteration limit for the warmstart. Default
  is 100.
- `main_solver(algorithm)`: The algorithm used in the main opimization process.
  The default is `Newton()`. Others include `SimulatedAnnealing()`, `SAMIN()`, `ParticleSwarm()`,
  `ConjugateGradient()`, `GradientDescent()`, `BFGS()`, `LBFGS()`,
  `NewtonTrustRegion()`, and `IPNewton()`. See
  http://julianlsolvers.github.io/Optim.jl/stable/ for details.
- `main_maxIT(::Int64)`: The iteration limit for the main estimation. Default
  is 2000.
- `tolerance(::Float64)`: The convergence criterion ("tolerance") based on the
  absolute value of gradients. Default is 1.0e-8. For non-gradient algorithms,
  it controls the main convergence tolerance, which is solver specific. 
  See `Optim`'s `g_tol` option for more information.
- `verbose(::Bool)`: Print on screen (`true`, the default) the information of
  the model and the optimization results.
- `banner(::Bool)`: Print on screen (`true`, the default) a banner to serve as
  a visual indicator of the start of the estimation.
- `ineff_index(::Bool)`: Whether to compute the Jondrow et al. (1982)
  inefficiency index and the Battese and Coelli (1988) efficiency index. The
  defauis `true`.
- `marginal(::Bool)`: Whether to compute the marginal effects of the exogenous
  determinants of inefficiency (if any).
- `rii(::Bool)`: Whether to compute the Idiosyncratic and correlated inefficiency.
- `table_format()`: The format to print the coefficient tables on the screen:
  `text` (default), `html`, or `latex`. A wrapper of `PrettyTables.jl`'s
  `backend` option.
- message::Bool: Whether printing (=true) or not (=false, the default) the
  confirmation message "A dictionary from ssfmodel_opt() is generated."
  on the screen after `ssfmodel_opt()` is successfully executed.


# Remarks
- `ssfmodel_opt(...)` is optional. It can be omitted entirely, or specifying
  only a partial list of the keywords.
- If any of the keywords are missing, default values are used.
- If warmstart is not needed, you need to give empty keyword values to
  warmstart related keys. E.g., either `warmstart_solver()` or
  `warmstart_maxIT()`, or both. Omitting the keyword entirely (i.e., not
  writing down `warmstart_solver` or `warmstart_maxIT`) will not skip the
  warmstart, but will reinstate the default. 
- Users do not need to provide gradient or Hessian functions even if 
  gradient-based optimization algorithms are used. The package uses automatic
  differentiation (https://en.wikipedia.org/wiki/Automatic_differentiation) to 
  compute the derivatives. It is not numerical finite differentiation. It is
  fast and as accurate as the symbolic differentiation.
- The `_dicOPT` is inheritable in the sense that an exiting `_dicOPT` (from
  the previous run of the same or a different model, for example) will be used
  if the current model does not have its own `ssfmodel_opt(...)`. This design
  has advantages in simulation studies where `ssfmodel_opt(...)` needs to be
  specified only once.

# Examples
```julia-repl
ssfmodel_opt(warmstart_solver(NelderMead()),   
            warmstart_maxIT(200),
            main_solver(Newton()), 
            main_maxIT(2000), 
            tolerance(1e-8),
            verbose(true),
            banner(true),
            ineff_index(true),
            marginal(true),
            rii(true),
            table_format(html),
            message = false)
```
"""
function ssfmodel_opt(arg::Vararg; message::Bool=false) # create a dictionary of maximization options

    global _dicOPT
           _dicOPT = Dict{Symbol, Any}()

    #* -- creates the default ---

    _dicOPT[:warmstart_solver] = :(NelderMead())
    _dicOPT[:warmstart_maxIT]  =  100
    _dicOPT[:main_solver]      = :(Newton())
    _dicOPT[:main_maxIT]       =  2000
    _dicOPT[:tolerance]        =  1.0e-8
    _dicOPT[:verbose]          =  true
    _dicOPT[:banner]           =  true
    _dicOPT[:ineff_index]      =  true
    _dicOPT[:marginal]         =  true
    _dicOPT[:rii]              =  true
    _dicOPT[:table_format]     = :(text)

    #* -- replace the defaults with the user's value ---

    for d in :($(arg))
        _dicOPT[d[1]] = d[2]
    end    
    
    #* ---- error checking --

    if (_dicOPT[:main_solver] === nothing) || (_dicOPT[:main_maxIT] === nothing) || (_dicOPT[:tolerance] === nothing)
         throw("You cannot give empty keyword values to `main_solver()`, `main_maxIT()`, or `tolerance()`. If you want to use the default, you may do so by dropping (not emptying keyword values) the keywords.")
    end
    
    if message 
      printstyled("A dictionary from sfmodel_opt() is generated.\n"; color = :green)  
    end  
    return _dicOPT # for debugging purpose

end    # end of ssfmodel_opt()




########################################################
###                 ssfmodel_fit()                   ####
########################################################
"""
    ssfmodel_fit(<keyword arguments>)

Maximum likelihood estimation of the stochastic frontier model specified 
in `ssfmodel_spec(...)`. Estimate the model parameters, calculate Jondrow et al. 
(1982) inefficiency index and Battese and Coelli (1988) efficiency index, 
compute marginal effects of inefficiency determinants (if any).
Return a dictionary with results.

# Arguments
- `useData(::DataFrame)`: The DataFrame used with the Method 1 of
  `ssfmodel_spec(...)`. 

# Remarks
- Use `Optim.jl` to carry out the estimation.
- Users do not need to provide gradient or Hessian functions even if 
  gradient-based optimization algorithms are used. The package uses automatic
  differentiation (https://en.wikipedia.org/wiki/Automatic_differentiation) to 
  compute the derivatives. AD is not numerical finite differentiation. AD is
  fast and as accurate as the symbolic differentiation.

# Examples
```julia-repl
ssfmodel_fit(useData(df))    # Method 1
ssfmodel_fit()               # Method 2
```
"""

function ssfmodel_fit()
  # For Method 2 of `ssfmodel_spec()`.

 !(_dicM[:hasDF]) || throw("Need to specify DataFrame in `ssfmodel_fit()`.")
 
 _dicM[:transfer] = true

#  println("ss",_dicM[:sdf])
 ssfmodel_fit(_dicM[:sdf])

end 

function ssfmodel_fit(sfdat::DataFrame) #, D1::Dict = _dicM, D2::Dict = _dicINI, D3::Dict = _dicOPT)

    (_dicM[:hasDF] || _dicM[:transfer])  || throw("You provided matrix in `ssfmodel_spec()` so you cannot specify a DataFrame in `ssfmodel_fit()`. Leave it blank.")

   #* for simulation, add a flag
   redflag::Bool = 0

  #* ###### Check if the OPT dictionary exists #####

      @isdefined(_dicINI) || ssfmodel_init()  # if not exist, create one with default values
      @isdefined(_dicOPT) || ssfmodel_opt()  
      

     if _dicOPT[:banner] 
        printstyled("\n###------------------------------------###\n"; color=:yellow)
        printstyled("###  Estimating SF models using Julia  ###\n"; color=:yellow)
        printstyled("###------------------------------------###\n\n"; color=:yellow)
      end  


  #* ##### Get variables from dataset #######
    
     # pos: (begx, endx, begz, endz, ...); variables' positions in the parameter vector.
     # num: (nofobs, nofx, ..., nofpara); number of variables in each equation
     # eqvec: ("frontier"=2, "μ"=6,...); named tuple of equation names and equation position in the table
     # eqvec2: (xeq=(1,3), zeq=(4,5),...); named tuple of equation and parameter positions, for ssfmodel_predict
     # varlist: ("x1", "x2",...); variable names for making table

     (minfo1, minfo2, pos, num, eqvec, eqvec2, yvar, xvar, zvar, qvar, uvar, 
      vvar,   ρ, τ,  W,      rowIDT, varlist) = getvar(tagD[:modelid], sfdat)

  #* ### print preliminary information ########

    if _dicOPT[:verbose] 

      printstyled("*********************************\n "; color=:cyan)
      printstyled("      Model Specification:\n"; color=:cyan); 
      printstyled("*********************************\n"; color=:cyan)

      print("Model type: "); printstyled(minfo1; color=:yellow); println();println()
      printstyled(minfo2; color=:yellow); println()
    end

  #* ##### Get the type parameter #######

     _porc::Int64 = 1     

     if (_dicM[:type] == [:cost]) 
         _porc = -1
     end

  #* ########## Process initial value dictionary  #####
     #* --- Get OLS results and other auxiliary values. --- #

     β0     = xvar \ yvar;  # OLS estiamte, uses a pivoted QR factorization;
     resid  = yvar - xvar*β0
     sse    = sum((resid).^2)  
     ssd    = sqrt(sse/(size(resid,1)-1)) # sample standard deviation; σ² = (1/(N-1))* Σ ϵ^2
     ll_ols = sum(normlogpdf.(0, ssd, resid)) # ols log-likelihood
     sk_ols = sum((resid).^3) / ((ssd^3)*(size(resid,1))) # skewnewss of ols residuals

     #* --- Create the dictionary -----------

    #  if (:all_init in keys(_dicINI))
    #      sf_init = _dicINI[:all_init]
    #  else
         #*  Create ini vectors from user's values; if none, use the default.--- #      
         β_ini  = get(_dicINI, :frontier, β0) 
         δ1_ini = get(_dicINI, :eqz, ones(num.nofz) * 0.1) 
         ω_ini  = get(_dicINI, :eqq, ones(num.nofq) * 0.1) 
         δ2_ini = get(_dicINI, :equ, ones(num.nofu) * 0.1)
         γ_ini  = get(_dicINI, :eqv, ones(num.nofv) * 0.1) 
         ρ_ini  = get(_dicINI, :eqρ, ones(num.nofρ) * 0.1)
         τ_ini  = get(_dicINI, :eqτ, ones(num.nofτ) * 0.1)
         #*  Make it Array{Float64,1}; otherwise Array{Float64,2}. ---#     
         #*       Could also use sf_init[:,1]. *#
         sf_init = vcat(β_ini, δ1_ini, ω_ini, δ2_ini, γ_ini,ρ_ini,τ_ini)  
         sf_init = vec(sf_init)   
    #  end # if :all_init


  #* ############ Misc.  ################     
     # --- check if the number of initial values is correct 
        (length(sf_init) == num.nofpara) ||  throw("The number of initial values does not match the number of parameters to be estimated. Make sure the number of init values in ssfmodel_init() matches the number of variabls in ssfmodel_spec().") 

     # --- Make sure there is no numerical issue arising from int vs. Float64.
        sf_init = convert(Array{Float64,1}, sf_init) 

  #* ############# process optimization dictionary  #######

         if (_dicOPT[:warmstart_solver] === nothing) || (_dicOPT[:warmstart_maxIT] === nothing)
             do_warmstart_search = 0
         else 
             do_warmstart_search = 1
             sf_ini_algo  = eval(_dicOPT[:warmstart_solver])  # warmstart search algorithms
             sf_ini_maxit = _dicOPT[:warmstart_maxIT]         # warmstart search iter limit
         end    
             
     # ---- main maximization algorithm -----
         sf_algo  = eval(_dicOPT[:main_solver])    # main algorithm
         sf_maxit = _dicOPT[:main_maxIT] 
         sf_tol   = _dicOPT[:tolerance] 
         sf_table = _dicOPT[:table_format]

  #* ########  Start the Estimation  ##########

    #* ----- Define the problem's Hessian -----#


     _Hessian = TwiceDifferentiable(rho -> LL_T(tagD[:modelid], 
                           yvar, xvar, zvar, qvar, uvar, vvar,  ρ, τ, W, 
                           _porc, num.nofobs, pos, rho,
                                   rowIDT, _dicM[:misc]),
                     sf_init;               
                    autodiff = :finite); 


    #* ---- Make placeholders for dictionary recording purposes *#

    sf_init_1st_dic  = 0
    sf_init_2nd_dic  = 0
    sf_ini_algo_dic  = nothing
    sf_ini_maxit_dic = 0
    sf_total_iter    = 0

    _run = 1  # a counter; use the -if- instead of -for- to avoid using global variables

    if (do_warmstart_search == 1) && (_run == 1)  
 
        if _dicOPT[:verbose] 
            printstyled("The warmstart run...\n\n"; color = :green)
        end
 
        sf_init_1st_dic  = copy(sf_init) # for dict recording
        sf_ini_algo_dic  = sf_ini_algo
        sf_ini_maxit_dic = copy(sf_ini_maxit)

        # @time  
               mfun = optimize(_Hessian, 
                               sf_init,         # initial values  
                               sf_ini_algo,                   
                               Optim.Options(g_tol = sf_tol,
                                             iterations  = sf_ini_maxit, 
                                             store_trace = true,
                                             show_trace  = false))


        sf_total_iter += Optim.iterations(mfun) # for later use

        sf_init = Optim.minimizer(mfun)  # save as initials for the next run
        _run    = 2                      # modify the flag

        if _dicOPT[:verbose] 
            println()
            print("$mfun \n")
            print("The warmstart results are:\n"); printstyled(Optim.minimizer(mfun); color=:yellow); println("\n")
        end

   end  # if  (do_warmstart_search == 1) && (_run == 1)  

   if (do_warmstart_search == 0 ) || (_run == 2) # either no warmstart run and go straight here, or the 2nd run

       sf_init_2nd_dic = copy(sf_init) # for dict recording 

       if _dicOPT[:verbose] 
           println()
           printstyled("Starting the optimization run...\n\n" ; color = :green)
       end 
       
       # @time 
              mfun = optimize(_Hessian, 
                              sf_init,       # initial values  
                              sf_algo,       # different from search run
                              Optim.Options(g_tol = sf_tol,
                                            iterations  = sf_maxit, # different from search run
                                            store_trace = true,
                                            show_trace  = false))
       sf_total_iter += Optim.iterations(mfun)

       if _dicOPT[:verbose] 
             println()
             print("$mfun \n")  
             print("The resulting coefficient vector is:\n"); printstyled(Optim.minimizer(mfun); color=:yellow); println("\n")
       end 


       if isnan(Optim.g_residual(mfun)) || (Optim.g_residual(mfun) > 0.1) 
            redflag = 1
            printstyled("Note that the estimation may not have converged properly. The gradients are problematic (too large, > 0.1, or others).\n\n", color = :red)
       end 


       if Optim.iteration_limit_reached(mfun) 
             redflag = 1
             printstyled("Caution: The number of iterations reached the limit.\n\n"; color= :red)  
       end  
 
   end     # if (do_warmstart_search == 0 )....


  #* ###### Post-estimation process ############### 

      _coevec            = Optim.minimizer(mfun)  # coef. vec.
      numerical_hessian  = hessian!(_Hessian, _coevec)  # Hessain

     #* ------ Check if the matrix is invertible. ----

     var_cov_matrix = try
                         pinv(numerical_hessian)
                      catch err 
                         redflag = 1
                         checkCollinear(tagD[:modelid], xvar, zvar, qvar, uvar, vvar) # check if it is b/c of multi-collinearity in the data         
                         throw("The Hessian matrix is not invertible, indicating the model does not converge properly. The estimation is abort.")
                      end 
                      
          #* In some cases the matrix is invertible but the resulting diagonal
          #*    elements are negative. Check.

          if !all( diag(var_cov_matrix) .> 0 ) # not all are positive
               redflag = 1
               printstyled("Some of the diagonal elements of the var-cov matrix are non-positive, indicating problems in the convergence. The estimation is abort.\n\n"; color = :red)
               checkCollinear(tagD[:modelid], xvar, zvar, qvar, uvar, vvar) # check if it is b/c of multi-collinearity in the data
          end              

     #* ------- JLMS and BC index -------------------

     if _dicOPT[:ineff_index] 
        @views (_jlms,_bc) = jlmsbc(tagD[:modelid], _porc, pos, _coevec, 
                                     yvar, xvar, zvar, qvar, uvar, vvar,  W,       rowIDT)
        _jlmsM = mean(_jlms) # sum(_jlms)/length(_jlms)
        _bcM   = mean(_bc)   # sum(_bc)/length(_bc)
     else
        _jlms  = nothing
        _bc    = nothing
        _jlmsM = nothing
        _bcM   = nothing
     end 


     #* ---- marginal effect on E(u) -------------- 
 
     if _dicOPT[:marginal] 
        margeff, margMinfo = get_marg(tagD[:modelid], pos, num, _coevec, zvar, qvar, uvar,W,       rowIDT)
     else
        margeff, margMinfo = nothing, ()
     end                         

     #* ---- marginal effect on E(u) -------------- 
 
     if _dicOPT[:rii] 
      RII_ratio = getRII(tagD[:modelid], _porc, pos, _coevec,  yvar, xvar, zvar, qvar, uvar, vvar,  W,       rowIDT)
     else
        RII_ratio  = nothing 
     end                         



     #* ------- Make Table ------------------

     stddev  = sqrt.(diag(var_cov_matrix)) # standard error
     t_stats = _coevec ./ stddev          # t statistics
     p_value = zeros(num.nofpara)   # p values
     ci_low  = zeros(num.nofpara) # confidence interval
     ci_upp  = zeros(num.nofpara) 
     tt      = cquantile(Normal(0,1), 0.025)

     for i = 1:num.nofpara 
         @views p_value[i,1] = pvalue(TDist(num.nofobs - num.nofpara), t_stats[i,1]; tail=:both)
         @views ci_low[i,1] = _coevec[i,1] - tt*stddev[i,1]
         @views ci_upp[i,1] = _coevec[i,1] + tt*stddev[i,1]
     end  

       #* Build the table columns *#

    table = zeros(num.nofpara, 7)  # 7 columns in the table
    table[:,2] = _coevec   # estiamted coefficients
    table[:,3] = stddev    # std deviation
    table[:,4] = t_stats   # t statistic
    table[:,5] = p_value   # p value
    table[:,6] = ci_low
    table[:,7] = ci_upp
    table      = [" " "Coef." "Std. Err." "z" "P>|z|" "95%CI_l" "95%CI_u"; table]  # add to top of the table

       #*  creating a column of function names 

    table[:, 1] .= ""
    for i in 1:length(eqvec)
        @views j = eqvec[i]
        @views table[j,1] = keys(eqvec)[i]
    end

       #*  Add the column of variable names

    table = hcat(varlist, table)                      # combine the variable names column (hcat, horizontal concatenate; see also vcat)
    table[:,1], table[:,2] = table[:,2], table[:,1]   # swap the first name column and the function column
    table[1,2] = "Var."

     # * ------ Print Results ----------- *#

     if _dicOPT[:verbose] 

         printstyled("*********************************\n "; color=:cyan)
         printstyled("      Estimation Results:\n"; color=:cyan); 
         printstyled("*********************************\n"; color=:cyan)

         print("Model type: "); printstyled(minfo1; color=:yellow); println()
         print("Number of observations: "); printstyled(num.nofobs; color=:yellow); println()
         print("Number of total iterations: "); printstyled(sf_total_iter; color=:yellow); println()
         if Optim.converged(mfun) 
             print("Converged successfully: "); printstyled(Optim.converged(mfun); color=:yellow); println()
         elseif Optim.converged(mfun) == false
             print("Converged successfully: "); printstyled(Optim.converged(mfun); color=:red); println()
             redflag = 1
         end         
         print("Log-likelihood value: "); printstyled(round(-1*Optim.minimum(mfun); digits=5); color=:yellow); println()
         println()
     
         pretty_table(table[2:end,:],    # could print the whole table as is, but this prettier
                      header=["", "Var.", "Coef.", "Std.Err.", "z", "P>|z|", 
                              "95%CI_l", "95%CI_u"],
                      formatters = ft_printf("%5.4f", 3:8),
                      compact_printing = true,
                      backend = Val(sf_table))
         println()


         # *----- Auxiliary Table, log parameters to original scales --------

         auxtable = Array{Any}(undef,2,3)
         rn = 0 # row index

         if size(uvar,2) == 1 # single variable
            varstd = sqrt(sum((uvar .- sum(uvar)/length(uvar)).^2)/length(uvar)) 
            if varstd  <= 1e-7 # constant, assuming =1 anyway
                rn += 1
                auxtable[rn, 1] = :σᵤ²
                auxtable[rn, 2] = exp(_coevec[pos.begu])
                auxtable[rn, 3] = exp(_coevec[pos.begu])*stddev[pos.begu]
            end
         end

         if size(vvar,2) == 1 # single variable
            varstd = sqrt(sum((vvar .- sum(vvar)/length(vvar)).^2)/length(vvar)) 
            if varstd  <= 1e-7 # constant
                rn += 1
                auxtable[rn, 1] = :σᵥ²
                auxtable[rn, 2] = exp(_coevec[pos.begv])
                auxtable[rn, 3] = exp(_coevec[pos.begv])*stddev[pos.begv]
            end
         end

         if rn >= 1  # table is non-empty
             println("Convert the constant log-parameter to its original scale, e.g., σ² = exp(log_σ²):")   
             pretty_table(auxtable[1:rn,:],
                          header=["", "Coef.", "Std.Err."],
                          formatters = ft_printf("%5.4f", 2:3),
                          compact_printing = true,
                          backend = Val(sf_table))

             print("\nTable format: "); printstyled("$(sf_table)"; color=:yellow); println(". Use ssfmodel_opt() to choose between text, html, and latex.")
             println()
         end

         printstyled("***** Additional Information *********\n"; color=:cyan)
 

         print("* OLS (frontier-only) log-likelihood: "); printstyled(round(ll_ols; digits=5); color=:yellow); println("")
         print("* Skewness of OLS residuals: "); printstyled(round(sk_ols; digits=5); color=:yellow); println("")
         if _dicOPT[:ineff_index] 
            print("* The sample mean of the JLMS inefficiency index: "); printstyled(round(_jlmsM; digits=5); color=:yellow); println("")
            print("* The sample mean of the BC efficiency index: "); printstyled(round(_bcM; digits=5); color=:yellow); println("\n")
         end
        #  if length(margMinfo) >= 1
        #     print("* The sample mean of inefficiency determinants' marginal effects on E(u): " ); printstyled(margMinfo; color=:yellow); println("")
        #     println("* Marginal effects of the inefficiency determinants at the observational level are saved in the return. See the follows.\n")
        #  end

         println("* Use `name.list` to see saved results (keys and values) where `name` is the return specified in `name = ssfmodel_fit(..)`. Values may be retrieved using the keys. For instance:")
         println("   ** `name.loglikelihood`: the log-likelihood value of the model;")
         println("   ** `name.jlms`: Jondrow et al. (1982) inefficiency index;")
         println("   ** `name.bc`: Battese and Coelli (1988) efficiency index;")
         println("   ** `name.marginal`: a DataFrame with variables' (if any) marginal effects on E(u).")
         println("* Use `keys(name)` to see available keys.")

         printstyled("**************************************\n\n\n"; color=:cyan)

     end  # if_verbose

  #* ########### create a dictionary and make a tuple for return ########### *#
     
      _dicRES = OrderedDict{Symbol, Any}()     
      _dicRES[:converged]          = Optim.converged(mfun)
      _dicRES[:iter_limit_reached] = Optim.iteration_limit_reached(mfun)
      _dicRES[:_______________] = "___________________"  #33
      _dicRES[:n_observations]  = num.nofobs
      _dicRES[:loglikelihood]   = -Optim.minimum(mfun)
      _dicRES[:table]           = [table][1]
      _dicRES[:coeff]           = _coevec
      _dicRES[:std_err]         = stddev
      _dicRES[:var_cov_mat]     = [var_cov_matrix][1]
      _dicRES[:jlms]            = _jlms
      _dicRES[:bc]              = _bc
      _dicRES[:OLS_loglikelihood] = ll_ols
      _dicRES[:OLS_resid_skew]    = sk_ols
      _dicRES[:marginal]      = margeff
      _dicRES[:marginal_mean] = margMinfo
      _dicRES[:rii] = RII_ratio
      _dicRES[:_____________] = "___________________"  #31      
      _dicRES[:model]         = minfo1      
#     _dicRES[:data]          = "$sfdat"
      _dicRES[:depvar]        = _dicM[:depvar]
      _dicRES[:frontier]      = _dicM[:frontier]
      _dicRES[:μ]             = _dicM[:μ]
      _dicRES[:hscale]        = _dicM[:hscale]        

      _dicRES[:σᵤ²]           = _dicM[:σᵤ²]
      _dicRES[:σᵥ²]           = _dicM[:σᵥ²]

      _dicRES[:log_σᵤ²]       = _dicM[:σᵤ²] 
      _dicRES[:log_σᵥ²]       = _dicM[:σᵥ²]
      _dicRES[:type]          = _dicM[:type]
      _dicRES[:dist]          = _dicM[:dist]
      _dicRES[:PorC]          = _porc
      _dicRES[:timevar]       = _dicM[:timevar]  
      _dicRES[:idvar]         = _dicM[:idvar] # for bootstrap marginal effect
      _dicRES[:table_format]  = _dicOPT[:table_format]
      _dicRES[:modelid]       = tagD[:modelid]
      _dicRES[:verbose]       = _dicOPT[:verbose]
      _dicRES[:hasDF]         = _dicM[:hasDF]
      _dicRES[:transfer]      = _dicM[:transfer]

    for i in 1:length(eqvec2)
        _dicRES[keys(eqvec2)[i]] = _coevec[eqvec2[i]]
    end

      _dicRES[:________________]  = "___________________" #34
      _dicRES[:Hessian]           = [numerical_hessian][1]
      _dicRES[:gradient_norm]     = Optim.g_residual(mfun)
    # _dicRES[:trace]             = Optim.trace(mfun)     # comment out because not very informative and size could be large
      _dicRES[:actual_iterations] = Optim.iterations(mfun)
      _dicRES[:______________] = "______________________" #32
      _dicRES[:warmstart_solver] = sf_ini_algo_dic
      _dicRES[:warmstart_ini]    = sf_init_1st_dic
      _dicRES[:warmstart_maxIT]  = sf_ini_maxit_dic
      _dicRES[:main_solver]      = sf_algo
      _dicRES[:main_ini]         = sf_init_2nd_dic
      _dicRES[:main_maxIT]       = sf_maxit
      _dicRES[:tolerance]        = sf_tol
      _dicRES[:eqpo]             = eqvec2

      _dicRES[:redflag]          = redflag

     #* ----- Delete optional keys that have value nothing, 

         for k in (:μ, :hscale,  ) 
             if _dicRES[k] === nothing 
                delete!(_dicRES, k)
             end 
         end

     #* ----- Create a NamedTuple from the dic as the final output; 
     #* -----     put the dic in the tuple.

         _ntRES = NamedTuple{Tuple(keys(_dicRES))}(values(_dicRES))
         _ntRES = (; _ntRES..., list    = _dicRES)

     #* ---- Create a gloal dictionary for sf_predict ---- 

        global _eqncoe 
        _eqncoe = Dict{Symbol, Vector}()  # nullify and initiate new dictionaries when a new model is specified

        for i in 1:length(eqvec2)
            _eqncoe[keys(eqvec)[i]]  = _coevec[eqvec2[i]] # for sf_predict
        end

  #* ############  make returns  ############ *#

      return _ntRES

end # ssfmodel_fit

