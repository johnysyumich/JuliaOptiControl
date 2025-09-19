# Comprehensive Tutorial for JuliaOptimalControl.jl

This is a detailed tutorial for the JuliaOptimalControl.jl package, covering installation, API reference, and comprehensive examples.

## 1. Installation

### Prerequisites
- Julia 1.10.2 or higher
- Package manager access

### Installing JuliaOptimalControl.jl

**Quick Installation** (Recommended):

```julia
using Pkg

# Install directly from GitHub - this automatically installs ALL dependencies!
Pkg.add("https://github.com/johnysyumich/JuliaOptiControl")

# Verify installation
using JuliaOptimalControl
```

That's it! Julia automatically installs all required dependencies (JuMP, Ipopt, Parameters, DataFrames, etc.) from the `Project.toml` file.

**Alternative installation methods:**

```julia
# Install specific branch/version
Pkg.add(url="https://github.com/johnysyumich/JuliaOptiControl", branch="main")

# Development mode for local development
Pkg.develop(path="/path/to/JuliaOptimalControl")

# Local installation if you have the source code
Pkg.activate(".")
Pkg.instantiate()
```

### Import Required Packages

```julia
using JuliaOptimalControl  # Our optimal control package
using JuMP                 # Mathematical optimization
using Plots                # Visualization
using DataFrames           # Data manipulation
using LinearAlgebra        # Linear algebra operations
```

## 2. Package Overview

JuliaOptimalControl.jl is a high-level package for solving optimal control problems using direct methods. It transforms continuous optimal control problems into nonlinear programming problems through discretization.

### Mathematical Foundation

The package solves optimal control problems of the form:

$$\min \int_{0}^{t_f} L(x(t), u(t), p, t) \, dt + \phi(x(t_f))$$

subject to:
- **Dynamics**: $\dot{x}(t) = f(x(t), u(t), p, t)$
- **Path constraints**: $g(x(t), u(t), p, t) \geq 0$
- **Boundary conditions**: $x(0) = x_0$, $x(t_f) \in S_f$
- **Variable bounds**: $x_L \leq x(t) \leq x_U$, $u_L \leq u(t) \leq u_U$

Where:
- $x(t)$ are the state variables
- $u(t)$ are the control inputs
- $p$ are parameters
- $L(\cdot)$ is the cost integrand
- $f(\cdot)$ defines the system dynamics
- $g(\cdot)$ defines path constraints

## 3. API Reference

### 3.1 Core Data Structures

#### OCP Structure
The main container for optimal control problems with five components:

```julia
mutable struct OCP{T <: Number}
    s::OCPSetting{T}      # Problem settings
    b::OCPBound{T}        # Bounds and boundary conditions
    f::OCPFormulation{T}  # Mathematical formulation
    p::OCPParameter{T}    # Optimization variables
    r::OCPResults{T}      # Solution results
end
```

### 3.2 Problem Definition Functions

#### defineOCP()
**Purpose**: Initialize an optimal control problem structure

**Signature**:
```julia
defineOCP(;
    numStates::Int,     # Number of state variables
    numControls::Int,   # Number of control variables
    X0,                 # Initial state values [numStates]
    XF,                 # Final state values [numStates] (NaN = free)
    XL,                 # State lower bounds [numStates]
    XU,                 # State upper bounds [numStates]
    CL,                 # Control lower bounds [numControls]
    CU                  # Control upper bounds [numControls]
) -> OCP
```

**Inputs**:
- `numStates::Int`: Number of state variables (required, > 0)
- `numControls::Int`: Number of control inputs (required, > 0)
- `X0::Vector`: Initial conditions. Use `NaN` for free initial states
- `XF::Vector`: Final conditions. Use `NaN` for free final states
- `XL::Vector`: State variable lower bounds. Use `-Inf` for unbounded
- `XU::Vector`: State variable upper bounds. Use `Inf` for unbounded
- `CL::Vector`: Control variable lower bounds. Use `-Inf` for unbounded
- `CU::Vector`: Control variable upper bounds. Use `Inf` for unbounded

**Outputs**:
- Returns an `OCP` structure ready for configuration

#### defineStates!() and defineControls!()
**Purpose**: Assign symbolic names to variables for clarity

**Signatures**:
```julia
defineStates!(ocp::OCP, states::Vector{Symbol})
defineControls!(ocp::OCP, controls::Vector{Symbol})
```

**Inputs**:
- `ocp::OCP`: The problem structure
- `states/controls::Vector{Symbol}`: Variable names (length must match `numStates`/`numControls`)

**Example**:
```julia
defineStates!(ocp, [:position, :velocity, :angle])
defineControls!(ocp, [:force, :torque])
```

#### defineTolerance!()
**Purpose**: Enable tolerance-based boundary conditions instead of exact equality

**Signature**:
```julia
defineTolerance!(ocp::OCP;
    X0_tol::Vector = fill(NaN, numStates),
    XF_tol::Vector = fill(NaN, numStates)
)
```

**Inputs**:
- `X0_tol::Vector`: Tolerance for initial conditions. Changes constraint to: `X0[i] - tol[i] ≤ x(0)[i] ≤ X0[i] + tol[i]`
- `XF_tol::Vector`: Tolerance for final conditions. Changes constraint to: `XF[i] - tol[i] ≤ x(tf)[i] ≤ XF[i] + tol[i]`

### 3.3 Problem Configuration

#### ConfigurePredefined()
**Purpose**: Configure the mathematical formulation and discretization

**Signature**:
```julia
ConfigurePredefined(ocp::OCP;
    # Required arguments
    Np::Int,                    # Number of discretization points
    dx::Function,               # Dynamics function f(x,u,p)

    # Time configuration (choose one)
    tf::Real,                   # Fixed final time
    tfDV::Bool,                 # Variable final time flag

    # Integration scheme
    IntegrationScheme::Symbol,  # :RK1, :RK2, :RK3, :RK4, :trapezoidal, :bkwEuler

    # Optional arguments
    cons::Function,             # Path constraints g(x,u,p) ≥ 0
    expr::Function,             # Cost integrand L(x,u,p)
    params::Matrix              # Parameter values
) -> OCPFormulation
```

**Inputs**:

**Required**:
- `Np::Int`: Number of discretization points (typically 50-200)
- `dx::Function`: Dynamics function with signature `f(x, u, p) -> Vector` returning state derivatives

**Time Configuration (choose one)**:
- `tf::Real`: Fixed final time horizon
- `tfDV::Bool`: Set to `true` for variable final time optimization

**Integration Schemes**:
- `:RK1`: Forward Euler (1st order, fastest)
- `:RK2`: Runge-Kutta 2nd order
- `:RK3`: Runge-Kutta 3rd order
- `:RK4`: Runge-Kutta 4th order (most accurate)
- `:trapezoidal`: Trapezoidal rule
- `:bkwEuler`: Backward Euler (implicit, stable)

**Optional**:
- `cons::Function`: Path constraints with signature `g(x, u, p) -> Vector` (must return ≥ 0)
- `expr::Function`: Cost integrand with signature `L(x, u, p) -> Scalar`
- `params::Matrix`: Parameter values, size `(Np × num_params)` or `(1 × num_params)` for broadcasting

**Outputs**:
- Returns an `OCPFormulation` structure containing the discretized problem

#### OCPdef!()
**Purpose**: Build the complete optimization problem in JuMP

**Signature**:
```julia
OCPdef!(ocp::OCP, formulation::OCPFormulation)
```

**Function**: Transforms the mathematical formulation into a solvable nonlinear programming problem by:
- Creating JuMP decision variables for states and controls
- Applying bounds and boundary conditions
- Generating dynamics constraints using the specified integration scheme
- Adding path constraints
- Setting up the solver

### 3.4 Optimization and Results

#### OptSolve!()
**Purpose**: Solve the optimal control problem

**Signature**:
```julia
OptSolve!(ocp::OCP)
```

**Function**:
- Calls the nonlinear solver (Ipopt by default)
- Extracts solution status and timing information
- Populates the results structure

**Outputs** (stored in `ocp.r`):
- `X::Matrix{Float64}`: State trajectory, size `(Np × numStates)`
- `U::Matrix{Float64}`: Control trajectory, size `(Np × numControls)`
- `Tst::Vector{Float64}`: Time points, length `Np`
- `dt::Vector{Float64}`: Time intervals, length `Np-1`
- `Status::Symbol`: Solution status (`:Optimal`, `:UserLimit`, `:Infeasible`)
- `TSolve::Float64`: Solve time in seconds
- `Objval::Float64`: Optimal objective value
- `IterNum::Int64`: Number of solver iterations
- `Dfs::Vector{DataFrame}`: Results in DataFrame format

### 3.5 Utility Functions

---

#### ExprIntegral()

**Purpose**: Compute numerical integration of cost expressions

**Signature**:
```julia
ExprIntegral(ocp::OCP) -> JuMP.AffExpr
```

**Functionality**:
• Integrates cost expressions using the same integration scheme as dynamics
• Automatically applies proper quadrature rules based on chosen integration method
• Returns a JuMP expression suitable for use in objective functions

---

#### WarmStart()

**Purpose**: Initialize optimization variables with previous solution for faster convergence

**Signature**:
```julia
WarmStart(ocp::OCP)
```

**Functionality**:
• Sets initial guess for state and control variables using previous solution
• Enables faster solver convergence for iterative applications like MPC
• Automatically checks if solver supports warm starting

---

#### UpdateX0!()

**Purpose**: Update initial conditions and apply warm start (for MPC)

**Signature**:
```julia
UpdateX0!(ocp::OCP, X0::Vector)
```

**Functionality**:
• Updates stored initial conditions
• Forces new equality constraints for initial states
• Automatically applies warm starting with shifted previous solution

---

#### ResultsToDataFrame()

**Purpose**: Convert results to DataFrame format for easy analysis

**Signature**:
```julia
ResultsToDataFrame(ocp::OCP) -> DataFrame
```

**Functionality**:
• Creates structured data table with named columns
• Includes time vector and all state/control trajectories
• Uses symbolic names defined in `defineStates!()` and `defineControls!()`

**Outputs**: DataFrame with columns:
- `:t`: Time points
- State variable columns (using names from `defineStates!()`)
- Control variable columns (using names from `defineControls!()`)

## 4. Examples

The following examples demonstrate the use of JuliaOptimalControl.jl for various optimal control problems. Complete implementations with detailed explanations are available in the [`examples/`](examples/) directory.

### Available Examples

#### 4.1 CartPole Swing-Up Control
**Location**: [`examples/cartpole/`](examples/cartpole/)

Classic underactuated control problem demonstrating energy-based strategies for swinging up an inverted pendulum.

**Features**:
- Nonlinear underactuated dynamics
- Energy transfer strategies
- Large-angle pendulum motion
- Trajectory visualization and animation

**Key Concepts**: Underactuated systems, energy control, nonlinear optimization

---

#### 4.2 Vehicle Obstacle Avoidance
**Location**: [`examples/ObstacleAvoidance/`](examples/ObstacleAvoidance/)

Optimal trajectory planning for vehicles navigating around static obstacles using bicycle vehicle dynamics.

**Features**:
- 7-state bicycle vehicle model
- Circular obstacle avoidance constraints
- Multi-objective cost function
- Real-time capable optimization

**Key Concepts**: Path planning, collision avoidance, autonomous vehicles

---

#### 4.3 Rocket Landing Control
**Location**: [`examples/rocket_landing/`](examples/rocket_landing/)

Optimal guidance for powered rocket landing with fuel-optimal trajectories and thrust vector control.

**Key Concepts**: Aerospace applications, fuel optimization, landing guidance

---


### Running Examples

Each example directory contains:
- **`README.md`**: Detailed problem description and mathematical formulation
- **`main.jl`**: Main execution script
- **Supporting files**: Dynamics models, parameters, and utility functions

To run any example:

```julia
# Navigate to the desired example directory
cd("examples/cartpole")  # or vehicleOpt, racing, etc.

# Execute the example
include("main.jl")
```
