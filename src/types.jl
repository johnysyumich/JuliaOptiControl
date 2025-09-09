# Author: Siyuan Yu and Congkai Shen
# Position: Graduate Student in University of Michigan
# Date Created: 02/28/2024

"""
    JuliaOptimalControl.jl Type Definitions

This module defines the core data structures for solving optimal control problems (OCPs)
using direct methods. The types support various numerical integration schemes and
both fixed and variable time horizon problems.

Key Components:
- State and control variable definitions
- Solver configuration and settings
- Problem formulation with constraints and dynamics
- Solution results storage
"""

using JuMP
using Parameters
using DataFrames
using Ipopt

export  States,
        Control,
        Solver,
        OCPFormulation,
        OCPSetting,
        OCPBound,
        OCPParameter,
        OCPResults,
        OCP


_Ipopt_MPC_defaults = (
    "mu_strategy" => "adaptive",              # Barrier parameter update strategy
    "max_iter" => 100,                       # Maximum solver iterations for real-time
    "tol" => 6e-3,                          # Optimality tolerance (relaxed for speed)
    "dual_inf_tol" => 2.,                   # Dual infeasibility tolerance
    "constr_viol_tol" => 3e-1,              # Constraint violation tolerance
    "compl_inf_tol" => 3e-1,                # Complementarity tolerance
    "acceptable_tol" => 2e-2,               # Acceptable optimality tolerance
    "acceptable_constr_viol_tol" => 0.02,   # Acceptable constraint violation
    "acceptable_dual_inf_tol" => 1e10,      # Acceptable dual infeasibility
    "acceptable_compl_inf_tol" => 0.02,     # Acceptable complementarity
    "warm_start_init_point" => "yes",       # Enable warm starting for MPC
    "fixed_variable_treatment" => "relax_bounds", # Handle fixed variables
    "max_cpu_time" => 0.1,                  # Maximum solve time for real-time MPC
    "print_level" => 0                      # Suppress solver output
    )

_Ipopt_defaults = (
    "mu_strategy" => "adaptive",              # Barrier parameter update strategy
    "max_iter" => 1000,                      # Maximum solver iterations (higher for accuracy)
    "warm_start_init_point" => "no",         # Disable warm starting by default
    "fixed_variable_treatment" => "relax_bounds", # Handle fixed variables
    "max_cpu_time" => 10.0,                  # Maximum solve time (longer for offline problems)
    "print_level" => 0                       # Suppress solver output
    )

"""
State variable configuration structure.
"""
@with_kw mutable struct States
    num::Int                                    = 0                                 # Number of States
    name::Vector{Symbol}                        = Vector{Symbol}[]                  # Names of States Variable
    pts::Int                                    = 0                                 # Number of points in States Variable
end

"""
Control variable configuration structure.
"""
@with_kw mutable struct Control
    num::Int                                    = 0                                 # Number of Control
    name::Vector{Symbol}                        = Vector{Symbol}[]                  # Names of Control Variable
    pts::Int                                    = 0                                 # Number of points in Control Variable
end

"""
NLP solver configuration structure.
"""
@with_kw mutable struct Solver
    name::Symbol                                = :Ipopt                            # Default solver name
    settings::Tuple                             = _Ipopt_defaults                   # Default settings
end

"""
Mathematical problem formulation structure.
"""
@with_kw mutable struct OCPFormulation{ T <: Number }
    tfDV::Bool                                  = false                             # Determines whether tf is a design variable
    Np::Int64                                   = 0                                 # Number of discretization points
    IntegrationScheme::Vector{Symbol}           = Vector{Symbol}()                  # Integration scheme
    tf::Any                                     = Any                               # Total time, type depends on tfDV
    tw::Vector{Float64}                         = Vector{Float64}()                 # Time weights
    TInt::Vector{Any}                           = Vector{Any}[]                     # Time intervals
    mdl::JuMP.Model                             = JuMP.Model()                      # JuMP model
    dx::Vector{Any}                             = Vector{Any}()                     # Dynamics function handles
    cons::Vector{Any}                           = Vector{Any}()                     # Constraints function handles
    expr::Vector{Any}                           = Vector{Any}()                     # Expression function handles
    params::Matrix{Any}                         = Matrix{Any}(undef, 0, 0)          # Parameter placeholder
end

"""
Problem configuration settings structure.
"""
@with_kw mutable struct OCPSetting{ T <: Number }
    states::States                              = States()                          # States structure in OCPSetting
    control::Control                            = Control()                         # Control structure in OCPSetting
    solver::Solver                              = Solver()                          # Solver structure in OCPSetting
    InternalLogging::Bool                       = true                              # Bool for logging data internally
    X0slack::Bool                               = false                             # Boolean for using tolerance on initial states
    XFslack::Bool                               = false                             # Boolean for using tolerance on final states
end

"""
Bounds and boundary conditions structure.
"""
@with_kw mutable struct OCPBound{ T <: Number }
    tfMin::Float64                              = 0.0                               # tf minimum
    tfMax::Float64                              = 100                               # tf maximum
    X0::Vector{T}                               = Vector{T}[]                       # Initial Condition
    X0_tol::Vector{T}                           = Vector{T}[]                       # Initial Condition tolerance
    XF::Vector{T}                               = Vector{T}[]                       # Final Condition
    XF_tol::Vector{T}                           = Vector{T}[]                       # Final Condition tolerance
    XL::Vector{T}                               = Vector{T}[]                       # Lower bound of states variable
    XU::Vector{T}                               = Vector{T}[]                       # Upper bound of states variable
    CL::Vector{T}                               = Vector{T}[]                       # Lower bound of control variable
    CU::Vector{T}                               = Vector{T}[]                       # Upper bound of control variable
end

"""
JuMP variables and parameters structure.
"""
@with_kw mutable struct OCPParameter{ T <: Number }
    x::Matrix{Any}                              = Matrix{Any}(undef,0,0)            # State variables matrix
    u::Matrix{VariableRef}                      = Matrix{VariableRef}(undef,0,0)    # Control variables matrix
    params::Matrix{VariableRef}                 = Matrix{VariableRef}(undef, 0, 0)  # Parameters for dynamics
    tV::Any                                     = Any                               # Time vector
    xvar::Matrix{VariableRef}                   = Matrix{VariableRef}(undef,0, 0)   # Intermediate variables for RK methods
    δx::Matrix{Any}                             = Matrix{Any}(undef,0,0)            # State derivatives
end

"""
Optimization results and solution data structure.
"""
@with_kw mutable struct OCPResults{ T <: Number }
    X::Matrix{Float64}                          = Matrix{Float64}(undef,0,0)        # State trajectory
    U::Matrix{Float64}                          = Matrix{Float64}(undef,0,0)        # Control trajectory
    Tst::Vector{Float64}                        = Vector{Float64}()                 # Time points
    dt::Vector{Float64}                         = Vector{Float64}()                 # Time intervals
    Status::Symbol                              = :InFeasible                       # Solution status
    IterNum::Int64                              = 0                                 # Solver iterations
    EvalNum::Int64                              = 0                                 # Function evaluations
    TerminalStatus::MOI.TerminationStatusCode   = MOI.OTHER_ERROR                   # Termination status
    TSolve::Float64                             = 0.0                               # Solve time
    Objval::Float64                             = 0.0                               # Objective value
    Dfs::Vector{DataFrame}                      = Vector{DataFrame}()               # DataFrame results
end

"""
Main optimal control problem container.
"""
@with_kw mutable struct OCP{ T <: Number }
    s::OCPSetting{T}                            = OCPSetting{T}()        # Settings
    b::OCPBound{T}                              = OCPBound{T}()          # Bounds
    f::OCPFormulation{T}                        = OCPFormulation{T}()    # Formulation
    p::OCPParameter{T}                          = OCPParameter{T}()      # Parameters
    r::OCPResults{T}                            = OCPResults{T}()        # Results
end
