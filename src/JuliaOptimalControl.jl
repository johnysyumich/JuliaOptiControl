isdefined(Base, :__precompile__) && __precompile__()

"""
    JuliaOptimalControl

A Julia package for solving optimal control problems using direct methods.

This package provides a high-level interface for formulating and solving
optimal control problems (OCPs) through discretization and nonlinear programming.

# Key Features
- Multiple integration schemes (Euler, Runge-Kutta 2/3/4, Trapezoidal)
- Variable and fixed time horizon problems
- Path constraints and boundary conditions
- Warm starting for MPC applications
- JuMP-based formulation with Ipopt solver

# Basic Usage
```julia
using JuliaOptimalControl

# Define problem structure
ocp = defineOCP(numStates=2, numControls=1, ...)

# Configure problem formulation
formulation = ConfigurePredefined(ocp,
    tf=10.0, Np=51, IntegrationScheme=:RK4,
    dx=(x,u,p) -> [x[2]; u[1]])

# Build and solve
OCPdef!(ocp, formulation)
OptSolve!(ocp)

# Extract results
plot(ocp.r.Tst, ocp.r.X)
```

# Mathematical Background
Solves optimal control problems of the form:
    min  ∫₀ᵗᶠ L(x,u,p,t) dt + φ(x(tf))
    s.t. ẋ = f(x,u,p,t)
         g(x,u,p,t) ≥ 0
         x(0) = x₀, x(tf) ∈ Sf
         XL ≤ x ≤ XU, CL ≤ u ≤ CU

Through direct transcription into nonlinear programming problems.
"""
module JuliaOptimalControl

using JuMP
using Parameters
using DataFrames
using Ipopt

include("types.jl")

include("utils.jl")

export  OptSolve!,
        ExprIntegral,
        RetrieveSolveStatus,
        GetOptimizeValue!,
        ResultsToDataFrame,
        CreateEmptyFormulation,
        DeleteElement,
        WarmStart,
        UpdateX0!

include("setup.jl")

export  defineOCP,
        defineStates!,
        defineControls!,
        defineTolerance!,
        ValidateScheme,
        ConfigurePredefined,
        CheckOCPFormulation,
        defineSolver!,
        CalXvar,
        OCPdef!

end