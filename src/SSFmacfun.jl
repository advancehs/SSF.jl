function misc(arg::Any)
    return (:misc, arg)
end

#? ---- functions for ssfmodel_spec() 
function weight(arg::Matrix)    # for smodel_init(b0)
    return (:W,  arg)
  end

#--------------------------------------------------------------------------------

function ssfpanel(arg::Vararg)
    (length(arg)==1) || throw("`ssfpanel()` in ssfmodel_spec should have one single string.")
    
    if !(arg[1] ∈ (SMA,SAR))
         throw("The keyword of `ssfpanel` in ssfmodel_spec is specified incorrectly.")
    end 

    return (:panel, [Symbol(arg[1])])
 end


 #--------------------------------------------------------------------------
"""
ssfdist(arg::Vararg)

An argument in `ssfmodel_sepc()`. Specify the distribution assumption of the one-sided stochastic
variable (aka inefficiency term) of the model. Possible choices include `truncated` (or `trun`, `t`), `half` (or `h`).

See the help on `ssfmodel_spec()` for more information.

# Examples
```julia-repl
ssfmodel_spec(sfdist(t), ...)
ssfmodel_spec(sfdist(h), ...)
```
"""    
function ssfdist(arg::Vararg)

   (length(arg)==1) || throw("`sfdist()` in ssfmodel_spec should have one and only one word.")

    if arg[1] ∈ (Trun, trun, t)
       ha = :t
    elseif arg[1] ∈ (Half, half, h)
       ha = :h
    else   
       throw("The inefficiency distribution in `sfdist()` of ssfmodel_spec() is specified incorrectly.")             
    end     

   return (:dist, [ha])
end  

#--------------------------------------------------------------------------
"""
 ssftype(arg::Vararg)

An argument in `ssfmodel_sepc()`. Specify whether the model is a `production` (or `prod`) frontier
or a `cost` frontier.

See the help on `ssfmodel_spec()` for more information.

# Examples
```julia-repl
ssftype(production)
ssftype(cost)
```
"""    
function ssftype(arg::Vararg) 
   (length(arg)==1) || throw("`ssftype()` in ssfmodel_spec should be either production, prod, or cost.")    
   if (arg[1]!=production) && (arg[1]!=prod) && (arg[1]!=cost)
        throw("`ssftype()` in ssfmodel_spec should be either production, prod, or cost.")    
   end
   return (:type, [Symbol(arg[1])])
end  
#--------------------------------------------------------------------------

"""
depvar(arg::Vararg)
 
An argument in `ssfmodel_sepc()`. Specify the dependent variable using a
column name from a DataFrame. 
 
See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×5 DataFrame
│ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
│     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
│ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
│ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
│ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
│ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
⋮
│ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
│ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
│ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
│ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
│ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │

ssfmodel_spec(depvar(yvar), ...)
```
"""
 function depvar(arg::Vararg)
     return (:depvar, collect(:($(arg))))
 end

 # ------------------------------------------
 """
 frontier(arg::Vararg)
  
 An argument in `ssfmodel_sepc()`. Specify the variables in the `frontier`
 function using column names from a DataFrame. 
  
 See the help on `ssfmodel_spec()` for more information.
  
 # Examples
 ```julia-repl
 julia> df
 100×5 DataFrame
 │ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
 │     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
 ├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
 │ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
 │ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
 │ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
 │ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
 ⋮
 │ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
 │ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
 │ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
 │ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
 │ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │
 
 ssfmodel_spec(frontier(xvar1, xvar2, _cons), ...)
 ```
 """
  function frontier(arg::Vararg)   # for ssfmodel_spec()          
      return (:frontier, collect(:($(arg))))
  end 
  
 #---------------------------------------------------
 
  """
  μ(arg::Vararg)
   
  An argument in `ssfmodel_sepc()`. Specify the variables in the `μ`
  function using column names from a DataFrame. 
   
  See the help on `ssfmodel_spec()` for more information.
   
  # Examples
  ```julia-repl
  julia> df
  100×5 DataFrame
  │ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
  │     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
  ├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
  │ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
  │ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
  │ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
  │ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
  ⋮
  │ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
  │ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
  │ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
  │ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
  │ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │
  
  ssfmodel_spec( μ(zvar, _cons), ...)
  ```
  """  
 function μ(arg::Vararg)
     return (:μ, collect(:($(arg))))
 end
 
# -------------------
 
"""
σᵤ²(arg::Vararg)
 
An argument in `ssfmodel_sepc()`. Specify the variables in the `σᵤ²`
function using column names from a DataFrame. 
 
See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×5 DataFrame
│ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
│     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
│ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
│ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
│ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
│ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
⋮
│ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
│ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
│ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
│ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
│ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │

ssfmodel_spec( σᵤ²(zvar, _cons), ...)
```
"""  
function σᵤ²(arg::Vararg)
   return (:σᵤ², collect(:($(arg))))
end


# -------------------
 
"""
@σᵥ²(arg::Vararg)
 
An argument in `ssfmodel_sepc()`. Specify the variables in the `σᵥ²`
function using column names from a DataFrame. 
 
See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×5 DataFrame
│ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
│     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
│ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
│ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
│ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
│ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
⋮
│ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
│ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
│ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
│ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
│ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │

ssfmodel_spec( σᵥ²(zvar, _cons), ...)
```
"""    
function σᵥ²(arg::Vararg)
   return (:σᵥ², collect(:($(arg))))
end


"""
hscale(arg::Vararg)
 
An argument in `ssfmodel_sepc()`. Specify the variables in the `hscale`
function using column names from a DataFrame. 
 
See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×5 DataFrame
│ Row │ yvar  │ xvar1     │ xvar2     │ zvar      │ _cons   │
│     │ Int64 │ Float64   │ Float64   │ Float64   │ Float64 │
├─────┼───────┼───────────┼───────────┼───────────┼─────────┤
│ 1   │ 1     │ 0.0306449 │ 0.452148  │ 0.808817  │ 1.0     │
│ 2   │ 2     │ 0.460691  │ 0.296092  │ 0.454545  │ 1.0     │
│ 3   │ 3     │ 0.897503  │ 0.376972  │ 0.907454  │ 1.0     │
│ 4   │ 4     │ 0.682894  │ 0.776861  │ 0.161721  │ 1.0     │
⋮
│ 96  │ 96    │ 0.329647  │ 0.0914057 │ 0.825032  │ 1.0     │
│ 97  │ 97    │ 0.0781165 │ 0.338999  │ 0.761652  │ 1.0     │
│ 98  │ 98    │ 0.41394   │ 0.0063118 │ 0.295372  │ 1.0     │
│ 99  │ 99    │ 0.516381  │ 0.285415  │ 1.91995   │ 1.0     │
│ 100 │ 100   │ 0.944     │ 0.702226  │ -0.539848 │ 1.0     │

ssfmodel_spec( hscale(zvar, _cons), ...)
```
"""   
function hscale(arg::Vararg)
   return (:hscale, collect(:($(arg))))
end

#-----------------------------------
"""
timevar(arg::Vararg)
 
An argument in `ssfmodel_sepc()` for panel data models. Specify the column
name of a DataFrame that contain the time information of the panel data. 

See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×4 DataFrame
│ Row │ year  │ firm  │ yvar  │ xvar1     │
│     │ Int64 │ Int64 │ Int64 │ Float64   │
├─────┼───────┼───────┼───────┼───────────┤
│ 1   │ 2019  │ 1     │ 1     │ 0.77645   │
│ 2   │ 2020  │ 1     │ 2     │ 0.0782388 │
│ 3   │ 2021  │ 1     │ 3     │ 0.222884  │
│ 4   │ 2022  │ 1     │ 4     │ 0.762864  │
⋮
│ 96  │ 2022  │ 24    │ 96    │ 0.590184  │
│ 97  │ 2019  │ 25    │ 97    │ 0.364425  │
│ 98  │ 2020  │ 25    │ 98    │ 0.639463  │
│ 99  │ 2021  │ 25    │ 99    │ 0.500526  │
│ 100 │ 2022  │ 25    │ 100   │ 0.239137  │

ssfmodel_spec( timvar(year), ...)
```
"""   
function timevar(arg::Vararg)
   return (:timevar, collect(:($(arg))))
end

#----------------------------------

"""
idvar(arg::Vararg)
 
An argument in `ssfmodel_sepc()` for panel data models. Specify the column
name of a DataFrame that contain the individual's id information of the panel data. 

See the help on `ssfmodel_spec()` for more information.
 
# Examples
```julia-repl
julia> df
100×4 DataFrame
│ Row │ year  │ firm  │ yvar  │ xvar1     │
│     │ Int64 │ Int64 │ Int64 │ Float64   │
├─────┼───────┼───────┼───────┼───────────┤
│ 1   │ 2019  │ 1     │ 1     │ 0.77645   │
│ 2   │ 2020  │ 1     │ 2     │ 0.0782388 │
│ 3   │ 2021  │ 1     │ 3     │ 0.222884  │
│ 4   │ 2022  │ 1     │ 4     │ 0.762864  │
⋮
│ 96  │ 2022  │ 24    │ 96    │ 0.590184  │
│ 97  │ 2019  │ 25    │ 97    │ 0.364425  │
│ 98  │ 2020  │ 25    │ 98    │ 0.639463  │
│ 99  │ 2021  │ 25    │ 99    │ 0.500526  │
│ 100 │ 2022  │ 25    │ 100   │ 0.239137  │

ssfmodel_spec( idvar(firm), ...)
```
"""   
function idvar(arg::Vararg)
    return (:idvar, collect(:($(arg))))
end




#? ---- functions for ssfmodel_init() 
"""
frontier(arg::Vector)
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `frontier` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init(frontier([0.1, 0.2, 0.5]), ...)
b0 = ones(3)*0.1
ssfmodel_init( frontier(b0), ...)
```
"""
function frontier(arg::Vector)    # for smodel_init(b0)
    return (:frontier,  arg)
end

#-------------------------------------------
"""
 μ(arg::Vector)
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `μ` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init( μ([0.1, 0.2, 0.5]), ...)

b0 = ones(3)*0.1
ssfmodel_init( μ(b0), ...)
```
"""
 function μ(arg::Vector)
     return (:μ, arg)
 end


#-------------------------------------

 """
 hscale(arg::Vector)
  
 An argument in `sfmodel_init()`. Specify the initial values for coefficients in
 the `hscale` function.
  
 See the help on `sfmodel_init()` for more information.
  
 # Examples
 ```julia-repl
 sfmodel_init( hscale(0.1, 0.2, 0.5), ...)
 
 b0 = ones(3)*0.1
 sfmodel_init( hscale(b0), ...)
 ```
 """
 function hscale(arg::Vector)
     return (:hscale, arg)
 end



 #-------------------------------------
 """
 σᵤ²(arg::Vector)
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `σᵤ²` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init( σᵤ²(0.1, 0.2, 0.5), ...)

b0 = ones(3)*0.1
ssfmodel_init( σᵤ²(b0), ...)
```
"""
 function σᵤ²(arg::Vector)
     return (:σᵤ², arg)
 end
 

 #---------------------------------------------
 """
 σᵥ²(arg::Vector)
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `σᵥ²` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init( σᵥ²([0.1, 0.2, 0.5]), ...)

b0 = ones(3)*0.1
ssfmodel_init( σᵥ²(b0), ...)
```
"""
 function σᵥ²(arg::Vector)
     return (:σᵥ², arg)
 end
 
 #-----------------------------
 """
 ρ(arg::Real) 
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `ρ` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init( ρ(0.1), ...)

```
"""
 function ρ(arg::Real)
     return (:ρ, arg)
 end
 
 #-----------------------------
 """
 τ(arg::Real) 
 
An argument in `ssfmodel_init()`. Specify the initial values for coefficients in
the `τ` function.
 
See the help on `ssfmodel_init()` for more information.
 
# Examples
```julia-repl
ssfmodel_init( τ(0.1), ...)

```
"""
 function τ(arg::Real)
     return (:τ, arg)
 end


 #? ----  functions for ssfmodel_opt --------
 

 """
 warmstart_solver(arg)
  
 An argument in `ssfmodel_opt()`. Specify the algorithm used in the first-stage ("warmstart")
   optimization process.
 
 The default is `NelderMead()`. Others include `SimulatedAnnealing()`, `SAMIN()`, `ParticleSwarm()`,
   `ConjugateGradient()`, `GradientDescent()`, `BFGS()`, `LBFGS()`,
   `Newton()`, `NewtonTrustRegion()`, and `IPNewton()`. See
   http://julianlsolvers.github.io/Optim.jl/stable/ for details.
   Non-gradient based algorithms are recommended for the warmstart solver. 
 
 See the help on `ssfmodel_opt()` for more information.
  
 # Examples
 ```julia-repl
 ssfmodel_opt( warmstart_solver(NelderMead()), ...)
 ```
 """
  function warmstart_solver(arg=nothing)
      return (:warmstart_solver, arg)
  end
  

#----------------------------------

"""
warmstart_maxIT(arg)
 
An argument in `ssfmodel_opt()`. Specify the iteration limit for the warmstart. Default
is 100.

See the help on `ssfmodel_opt()` for more information.
 
# Examples
```julia-repl
ssfmodel_opt( warmstart_maxIT(400), ...)
```
"""
 function warmstart_maxIT(arg=nothing)
     return (:warmstart_maxIT, arg)
 end

#--------------------------------

    """
    main_solver(arg)
    
    An argument in `ssfmodel_opt()`. Specify the algorithm used in the 2nd-stage ("main")
    optimization process.

    The default is `Newton()`. Others include `SimulatedAnnealing()`, `SAMIN()`, `ParticleSwarm()`,
    `ConjugateGradient()`, `GradientDescent()`, `BFGS()`, `LBFGS()`,
    `NewtonTrustRegion()`, and `IPNewton()`. See
    http://julianlsolvers.github.io/Optim.jl/stable/ for details.


    See the help on `ssfmodel_opt()` for more information.
    
    # Examples
    ```julia-repl
    ssfmodel_opt( main_solver(Newton()), ...)
    ```
    """
function main_solver(arg=nothing)
    return (:main_solver, arg)
end

#--------------------------------
    """
    main_maxIT(arg)
    
    An argument in `ssfmodel_opt()`. Specify the iteration limit for the main solver. Default
    is 2000.
    
    See the help on `ssfmodel_opt()` for more information.
    
    # Examples
    ```julia-repl
    ssfmodel_opt( main_maxIT(2500), ...)
    ```
    """
 function main_maxIT(arg::Any=nothing)
     return (:main_maxIT, arg)
 end

#--------------------------------
    """
    tolerance(arg::Float64)
    
    An argument in `ssfmodel_opt()`. Specify the convergence criterion ("tolerance") based on the
    absolute value of gradients. Default is 1.0e-8. For non-gradient algorithms,
    it controls the main convergence tolerance, which is solver specific. 
    See `Optim`'s `g_tol` option for more information.

    Also see the help on `ssfmodel_opt()` for more information.
    
    # Examples
    ```julia-repl
    ssfmodel_opt( tolerance(1.0e-6), ...)
    ```
    """
 function tolerance(arg::Float64=nothing)
     return (:tolerance, arg)
 end



   # --------------------------
   """
   verbose(arg::Bool)
    
   An argument in `ssfmodel_opt()`. Specify whether to print on screen (`true`,
   the default) the information of the model and the optimization results.
   
   See the help on `ssfmodel_opt()` for more information.
    
   # Examples
   ```julia-repl
   ssfmodel_opt( verbose(false), ...)
   ```
   """ 
    function verbose(arg::Bool=false)
        return (:verbose, arg)
    end
    
   
    """
    banner(arg::Bool)
     
    An argument in `ssfmodel_opt()`. Specify whether to print on screen (`true`,
    the default) a banner to serve as a visual indicator of the start of the 
    estimation.
    
    See the help on `ssfmodel_opt()` for more information.
     
    # Examples
    ```julia-repl
    ssfmodel_opt( banner(false), ...)
    ```
    """ 
    function banner(arg::Bool=true)
       return (:banner, arg)
    end
    
   
    """
    ineff_index(arg::Bool)
     
    An argument in `ssfmodel_opt()`. Specify whether (`true`, the default) to compute the Jondrow et al. (1982)
    inefficiency index and the Battese and Coelli (1988) efficiency index.
    
    See the help on `ssfmodel_opt()` for more information.
     
    # Examples
    ```julia-repl
    ssfmodel_opt( ineff_index(false), ...)
    ```
    """ 
    function ineff_index(arg::Bool=true)
       return (:ineff_index, arg)
    end
   
# -------------------------
   
   """
   marginal(arg::Bool)
     
    An argument in `ssfmodel_opt()`. Specify whether (`true`, the default) to
    compute the marginal effects of the exogenous determinants of inefficiency (if any).
    
    See the help on `ssfmodel_opt()` for more information.
     
    # Examples
    ```julia-repl
    ssfmodel_opt( marginal(false), ...)
    ```
   """ 
    function marginal(arg::Bool=true)
       return (:marginal, arg)
    end
   
     # -------------------------
     """
     rii(arg::Bool)
       
      An argument in `ssfmodel_opt()`. Specify whether (`true`, the default) to
      compute the Idiosyncratic and correlated inefficiency.
      
      See the help on `ssfmodel_opt()` for more information.
       
      # Examples
      ```julia-repl
      ssfmodel_opt( marginal(false), ...)
      ```
     """ 
      function rii(arg::Bool=true)
         return (:rii, arg)
      end
     
       # -------------------------
   
     """
     table_format(arg)
       
      An argument in `ssfmodel_opt()`. Specify the format to print the coefficient
      tables on the screen: `text` (default), `html`, or `latex`. A wrapper of `PrettyTables.jl`'s
      `backend` option.
   
      See the help on `ssfmodel_opt()` for more information.
       
      # Examples
      ```julia-repl
      ssfmodel_opt( table_format(html), ...)
      ```
     """ 
    function table_format(arg=nothing)
       
       if !(arg ∈ (text, html, latex ))
           throw("The keyword of `table_format` in `sfmodel_opt()` is specified incorrectly. Only allow `text`, `html`, or `latex`. Got `$(arg)` instead.")
      end 
   
       return (:table_format, Symbol(arg))
   end
   





  #? ------ functions for ssfmodel_fit ----------

   """
 useData(D::DataFrame)
   
  An argument in `ssfmodel_fit()`. Specify the name of the DataFrame that
  contains the estimation data.

  See the help on `ssfmodel_fit()` for more information.
   
# Examples
```julia-repl
  julia> mydf
  100×4 DataFrame
  │ Row │ year  │ firm  │ yvar  │ xvar1     │
  │     │ Int64 │ Int64 │ Int64 │ Float64   │
  ├─────┼───────┼───────┼───────┼───────────┤
  │ 1   │ 2019  │ 1     │ 1     │ 0.77645   │
  │ 2   │ 2020  │ 1     │ 2     │ 0.0782388 │
  │ 3   │ 2021  │ 1     │ 3     │ 0.222884  │
  │ 4   │ 2022  │ 1     │ 4     │ 0.762864  │
  ⋮
  │ 96  │ 2022  │ 24    │ 96    │ 0.590184  │
  │ 97  │ 2019  │ 25    │ 97    │ 0.364425  │
  │ 98  │ 2020  │ 25    │ 98    │ 0.639463  │
  │ 99  │ 2021  │ 25    │ 99    │ 0.500526  │
  │ 100 │ 2022  │ 25    │ 100   │ 0.239137  │

  ssfmodel_fit(useData(mydf), ...)
```
""" 
 function useData(D::DataFrame)  # doesn't work using macro (perhaps because of data), so...
     return D
 end
 