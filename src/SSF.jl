module SSF


export ssfmodel_spec, ssfmodel_init, ssfmodel_opt, 
       ssfmodel_fit, #sfmodel_predict, drawTrun,
       #genTrun, sf_demean, sfmodel_boot_marginal,
       # likelihood functions 
        LL_T, 
       # macros for _spec; 
        depvar, timevar, idvar,
        ssfdist, ssftype,
         
        ssfpanel,
        
       # functions for sfmodel_init 
        frontier, 
        μ,  σᵤ²,  σᵥ²,
        ρ,  τ,
        misc,
        hscale,
        weight,
       # functions for sfmodel_opt
        warmstart_solver, warmstart_maxIT,
        main_solver, main_maxIT, tolerance, verbose, banner,
        ineff_index, marginal, table_format, rii,
       # functions for sfmodel_fit
         useData,
       # functions for JLMS and BC index
        jlmsbc,
       # the table for regular and mixed Chi-square test
        sfmodel_MixTable, sfmodel_ChiSquareTable,
       # struct 
        Sfmodeltype,
        Trun, trun, truncated, t, 
        Half, half, h,
        
        production, cost, #* not prod, b/c conflict with Base
        text, html, latex,
        SART, SARH, SMAT, SMAH, get_marg, #* export for testing purpose
        SMA, SAR,
      # Optim's algorithms  
        NelderMead, SimulatedAnnealing, SAMIN, ParticleSwarm,
        ConjugateGradient, GradientDescent, BFGS, LBFGS,
        Newton, NewtonTrustRegion, IPNewton


        
using Optim
using DataFrames
using NLSolversBase              # for hessian!
using StatsFuns                  # for normlogpdf(), normlogcdf()
using Statistics                 #
using HypothesisTests            # for pvalue()
using LinearAlgebra              # extract diagnol and Matrix(I,...)
using Distributions              # for TDist, Normal
using DataStructures             # for OrderedDict
using PrettyTables               # making tables 
using ForwardDiff                # for marginal effect
using QuadGK                     # for TFE_CSW2014 model
using RowEchelon                 # for checkCollinear, check multi-collinearity
using FLoops                     # multithreading
using KahanSummation             # for time decay model, true random effect model
using Random           


#############################
##   Define Model Types    ##
#############################


abstract type SSfmodeltype end
struct Trun      <: SSfmodeltype end
struct truncated <: SSfmodeltype end
struct trun      <: SSfmodeltype end
struct t         <: SSfmodeltype end 
struct Half      <: SSfmodeltype end
struct half      <: SSfmodeltype end
struct h         <: SSfmodeltype end 


struct SMAH      <: SSfmodeltype end
struct SMAT      <: SSfmodeltype end
struct SARH      <: SSfmodeltype end
struct SART      <: SSfmodeltype end
abstract type PanelModel end
  struct SMA  <: PanelModel end
  struct SAR <: PanelModel end


  abstract type PorC end
  struct production <: PorC end
  struct cost <: PorC end

  abstract type TableFormat end
  struct text   <: TableFormat end
  struct html   <: TableFormat end
  struct latex  <: TableFormat end


################################################
##    include other files; order important    ##
################################################


include("SSFmacfun.jl")
include("SSFloglikefunc.jl")
# include("SFcheck.jl")
include("SSFgetvars.jl")
include("SSFindex.jl")
include("SSFriiindex.jl")
# include("SSFpredict.jl")
include("SSFmarginal.jl")
include("SSFmainfun.jl")




end # module
