"""
Create and initialize an optimal control problem structure.
"""
function defineOCP(;
    numStates::Int = 0,
    numControls::Int = 0,
    X0 = fill(NaN,numStates),
    XF = fill(NaN,numStates),
    XL = fill(NaN,numStates),
    XU = fill(NaN,numStates),
    CL = fill(NaN,numControls),
    CU = fill(NaN,numControls))::OCP
    ocp = OCP{Float64}()

    # Input validation
    if numStates < 0 error("States number must be larger than 0") end
    if numControls < 0 error("Control number must be larger than 0") end
    if length(X0) != numStates error("Number of states do not match X0") end
    if length(XF) != numStates error("Number of states do not match XF") end
    if length(XL) != numStates error("Number of states do not match XL") end
    if length(XU) != numStates error("Number of states do not match XU") end
    if length(CL) != numControls error("Number of controls do not match CL") end
    if length(CU) != numControls error("Number of controls do not match CU") end

    # Initialize problem dimensions
    ocp.s.states.num    = numStates
    ocp.s.control.num   = numControls

    # Set boundary conditions and bounds
    ocp.b.X0            = X0        # Initial states
    ocp.b.X0_tol        = fill(NaN, numStates)  # No tolerance by default
    ocp.b.XF_tol        = fill(NaN, numStates)  # No tolerance by default
    ocp.b.XL            = XL        # State lower bounds
    ocp.b.XF            = XF        # Final states
    ocp.b.XU            = XU        # State upper bounds
    ocp.b.CL            = CL        # Control lower bounds
    ocp.b.CU            = CU        # Control upper bounds


    return ocp
end

"""
Assign symbolic names to state variables.
"""
function defineStates!(ocp::OCP, states::Vector{Symbol})
    if length(states) != ocp.s.states.num error("Wrong number of state variables") end
    ocp.s.states.name = states
    return nothing
end

"""
Assign symbolic names to control variables.
"""
function defineControls!(ocp::OCP, controls::Vector{Symbol})
    if length(controls) != ocp.s.control.num error("Wrong number of control variables") end
    ocp.s.control.name = controls
    return nothing
end
"""
Enable tolerance-based boundary conditions.
"""
function defineTolerance!(ocp::OCP; X0_tol=fill(NaN, ocp.s.states.num), XF_tol=fill(NaN, ocp.s.states.num))
    # Enable initial state tolerance if any values are provided (not all NaN)
    if !all(isnan.(X0_tol))
        ocp.s.X0slack = true        # Switch from equality to inequality constraints
        ocp.b.X0_tol = X0_tol       # Store tolerance values
    end

    # Enable final state tolerance if any values are provided (not all NaN)
    if !all(isnan.(XF_tol))
        ocp.s.XFslack = true        # Switch from equality to inequality constraints
        ocp.b.XF_tol = XF_tol       # Store tolerance values
    end
    return nothing
end


function ValidateScheme(scheme)
    num = maximum(size(scheme))
    for i in 1:num
        # Check if each scheme is in the list of supported integration methods
        if scheme[i] ∉ [:bkwEuler, :trapezoidal, :RK1, :RK2, :RK3, :RK4]
            return false, i  # Return failure and index of invalid scheme
        end
    end
    return true  # All schemes are valid
end


"""
Configure a predefined optimal control formulation.

┌─────────────────────────────────────────────────────────────────────────────┐
│                           KEYWORD ARGUMENTS                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│ REQUIRED ARGUMENTS:                                                         │
│   Np::Int               │ Number of discretization points                   │
│   dx::Function          │ Dynamics function f(x,u,p) → state derivatives    │
│                                                                             │
│ TIME CONFIGURATION (choose one):                                            │
│   tf::Real              │ Fixed final time horizon                          │
│   tfDV::Bool            │ Set to true for variable final time               │
│                                                                             │
│ INTEGRATION SCHEMES:                                                        │
│   IntegrationScheme::Symbol │ Integration method for all intervals:         │
│     :RK1                │ Forward Euler (1st order)                         │
│     :RK2                │ Runge-Kutta 2nd order                             │
│     :RK3                │ Runge-Kutta 3rd order                             │
│     :RK4                │ Runge-Kutta 4th order                             │
│     :trapezoidal        │ Trapezoidal rule                                  │
│     :bkwEuler           │ Backward Euler                                    │
│                                                                             │
│ OPTIONAL ARGUMENTS:                                                         │
│   cons::Function        │ Path constraint function g(x,u,p) ≥ 0             │
│   expr::Function        │ Cost integrand function L(x,u,p)                  │
│   params::Matrix        │ Parameter values for dynamics/constraints         │
│                         │   Size: (Np x num_params) or (1 x num_params)     │
└─────────────────────────────────────────────────────────────────────────────┘

Time Configuration Rules:
• If tfDV=true: tf becomes optimization variable with bounds [tfMin, tfMax]
• If tfDV=false: Must provide fixed tf value
• Time intervals: Δt = tf x (1/(Np-1)) for uniform spacing

Parameter Broadcasting:
• If params provided: Available as p argument in dynamics/constraints
• If no params: Use Nonparam placeholder in function calls
• Broadcasting: (1xn) parameters repeated for all time points
"""
function ConfigurePredefined(ocp::OCP; kwargs...)::OCPFormulation
    # Create empty formulation structure
    OCPForm = OCPFormulation{Float64}()
    OCPForm.mdl = JuMP.Model()     # Initialize JuMP optimization model
    set_silent(OCPForm.mdl)        # Suppress solver output by default
    kw = Dict(kwargs)              # Convert kwargs to dictionary for easier access
    # Configure final time as design variable or fixed parameter
    if haskey(kw, :tfDV)
        OCPForm.tfDV = get(kw, :tfDV, 0)
    end
    if haskey(kw, :tf)
        # Fixed time horizon provided - override tfDV setting
        if OCPForm.tfDV == true
            @warn "Conflict intructions to put tf as design variable, change to fixed horizon"
            OCPForm.tfDV = false
        end
        OCPForm.tf = get(kw, :tf, 0)
        # Validate time horizon bounds
        if OCPForm.tf <= ocp.b.tfMin || OCPForm.tf >= ocp.b.tfMax
            error("Please make sure tf ∈ [$(ocp.b.tfMin), $(ocp.b.tfMax)]")
        end
    else
        # No fixed time provided - must use tf as design variable
        if OCPForm.tfDV == false
            error("Either treat tf as design variable or use fixed time horizon")
        end
    end

    # Configure number of discretization points
    if haskey(kw, :Np)
        Np = get(kw, :Np, 0)
        NpType = typeof(Np)
        if NpType != Int64
            if NpType <: Real
                @warn "Round Np to nearest integer"
                OCPForm.Np = Int(round(Np))  # Convert to integer
            else
                error("Wrong type of number of points, should be Int64")
            end
        else
            OCPForm.Np = Np
        end
        # Validate positive integer
        if OCPForm.Np <= 0
            error("Make sure Np is an integer and is larger than 0")
        end
    else
        error("No number of point input")
    end

    # Configure integration scheme for all intervals
    if haskey(kw, :IntegrationScheme)
        IntegrationScheme_single = get(kw, :IntegrationScheme, 0)
        if IntegrationScheme_single ∈ [:bkwEuler, :trapezoidal, :RK1, :RK2, :RK3, :RK4]
            # Apply same scheme to all intervals (Np-1 intervals total)
            OCPForm.IntegrationScheme = vec(repeat([IntegrationScheme_single],OCPForm.Np-1))
        else
            @warn "$IntegrationScheme_single is not implemented, use default (bkwEuler)"
        end
    end

    # Set up time discretization
    if OCPForm.tfDV == false
        # Fixed time horizon: compute time intervals directly
        OCPForm.tw = 1 / (OCPForm.Np - 1) .* ones(OCPForm.Np - 1)  # Equal time weights
        OCPForm.TInt = OCPForm.tf .*  OCPForm.tw                   # Δt = tf * weight
    elseif OCPForm.tfDV == true
        # Variable time horizon: tf is optimization variable
        OCPForm.tw = 1 / (OCPForm.Np - 1) .* ones(OCPForm.Np - 1)  # Equal time weights
        OCPForm.tf = @variable(OCPForm.mdl, ocp.b.tfMin <= tf <= ocp.b.tfMax)  # tf as JuMP variable
        OCPForm.TInt = @expression(OCPForm.mdl, [idx = 1:OCPForm.Np - 1], OCPForm.tf * OCPForm.tw[idx])  # Δt expressions
    end

    if !haskey(kw, :dx)
        error("No dynamics here")
    else
        OCPForm.dx = Vector{Any}(nothing, OCPForm.Np)
        OCPForm.dx[1:OCPForm.Np] .= get(kw, :dx, 0)
    end
    if haskey(kw, :params)
        params = get(kw, :params, 0)
        if size(params, 1) == OCPForm.Np
            OCPForm.params = Matrix{Any}(undef, size(params, 1), size(params, 2))
            OCPForm.params[1:size(params, 1), 1:size(params, 2)] = params
        else
            if size(params, 2) == 1 && size(params, 1) != 1
                OCPForm.params = repeat(params', OCPForm.Np)
            elseif size(params, 1) == 1 && size(params, 2) != 1
                OCPForm.params = repeat(params, OCPForm.Np)
            elseif size(params, 1) == 1 && size(params, 2) == 1
                OCPForm.params = repeat(params, OCPForm.Np)
            else
                error("Size of param ($(size(params, 1)) * $(size(params, 2)) is not expected)")
            end
        end
    end
    if size(OCPForm.params, 1) == 0 || size(OCPForm.params, 2) == 0
        OCPForm.params = Matrix{Any}(undef, OCPForm.Np, 1)
    end

    OCPForm.cons = Vector{Any}(nothing, OCPForm.Np)
    if haskey(kw, :cons)
        OCPForm.cons .= get(kw, :cons, 0)
    end

    OCPForm.expr = Vector{Any}(nothing, OCPForm.Np)

    if haskey(kw, :expr)
        OCPForm.expr .= get(kw, :expr, 0)
    end
    return OCPForm
end

function CheckOCPFormulation(ocp::OCP, OCPForm::OCPFormulation)
    tw = OCPForm.tw
    # Verify time weights sum to 1 (normalized)
    if abs(1 - sum(tw)) > 1e-4
        error("Wrong weights of tf")
    end

    # Basic size consistency checks
    if length(OCPForm.tw) != length(OCPForm.TInt) error("Size of tw and TInt do not match") end
    if (OCPForm.Np < 0) || (typeof(OCPForm.Np) != Int64) error("Wrong input of Np") end

    # Validate tf type based on whether it's a design variable
    if OCPForm.tfDV == true && typeof(OCPForm.tf) != JuMP.VariableRef error("Wrong type of tf since tfDV is $(OCPForm.tfDV)") end
    if OCPForm.tfDV == false
        if !(typeof(OCPForm.tf) <: Real)
            error("Wrong type of tf since tfDV is $(OCPForm.tfDV)")
        else
            # Check fixed time horizon is within bounds
            if OCPForm.tf < ocp.b.tfMin || OCPForm.tf > ocp.b.tfMax
                error("Please make sure tf ∈ [$(ocp.b.tfMin), $(ocp.b.tfMax)]")
            end
        end
    end

    for i in eachindex(OCPForm.TInt)
        if OCPForm.tfDV == true 
            if OCPForm.TInt[i] != OCPForm.tw[i] * OCPForm.tf
                error("Weight and TInt do not correspond")
            end
        else
            if abs(OCPForm.TInt[i] - OCPForm.tw[i] * OCPForm.tf) > 1e-4
                error("Weight and TInt do not correspond")
            end
        end
    end
    ValidationInfo = ValidateScheme(OCPForm.IntegrationScheme)
    if length(ValidationInfo) == 2 # Midpoint Collocation?
        error("$(OCPForm.IntegrationScheme[ValidationInfo[2]]) is not implemented")
    end

    if size(OCPForm.params, 1) != OCPForm.Np
        error("Size [$(size(OCPForm.params, 1)), $(size(OCPForm.params, 2))] do not match number of points $(OCPForm.Np) ")
    end

    return nothing
end


function defineSolver!(OCPForm::OCPFormulation, SolverName::Symbol, Options::Tuple)
    if solver_name(OCPForm.mdl) == "No optimizer attached."
        if SolverName == :Ipopt
            set_optimizer(OCPForm.mdl, Ipopt.Optimizer)
            set_attributes(OCPForm.mdl, Options...)
        else
            error("Solver $SolverName is not implemented")
        end
    else
        @warn "$(solver_name(OCPForm.mdl)) is set, can not change to $SolverName"
    end
    return nothing
end


"""
Calculate number of intermediate variables needed for Runge-Kutta schemes.
"""
function CalXvar(Scheme)
    count = 0
    for i in 1:maximum(size(Scheme))
        if Scheme[i] == :RK2
            count += 1
        elseif Scheme[i] == :RK3
            count += 2
        elseif Scheme[i] == :RK4
            count += 3
        end
    end
    return Int32(ceil(count))
end

"""
Build the complete discretized optimal control problem in JuMP.
"""
function OCPdef!(ocp::OCP, OCPForm::OCPFormulation)
    ocp.f = OCPForm
    CheckOCPFormulation(ocp, OCPForm)

    # Configure solver if none attached
    if solver_name(OCPForm.mdl) == "No optimizer attached."
        defineSolver!(OCPForm, ocp.s.solver.name, ocp.s.solver.settings)
    end

    # Create decision variables with bounds
    ocp.s.states.pts = ocp.s.control.pts = OCPForm.Np  # Both use same discretization
    # State variables: x[time_point, state_index] with bounds XL ≤ x ≤ XU
    ocp.p.x = @variable(OCPForm.mdl, ocp.b.XL[i] <= x[j in 1:ocp.s.states.pts, i in 1:ocp.s.states.num] <= ocp.b.XU[i])
    # Control variables: u[time_point, control_index] with bounds CL ≤ u ≤ CU
    ocp.p.u = @variable(OCPForm.mdl, ocp.b.CL[i] <= u[j in 1:ocp.s.control.pts, i in 1:ocp.s.control.num] <= ocp.b.CU[i])
    # Store variables in model for easy access
    OCPForm.mdl[:x] = ocp.p.x
    OCPForm.mdl[:u] = ocp.p.u

    # Apply initial state constraints
    if !ocp.s.X0slack
        # Fixed initial states: x(0) = X0
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.X0[st])
                fix(OCPForm.mdl[:x][1, st], ocp.b.X0[st]; force = true)  # Fix state variable
            end
        end
    else
        # Tolerance-based initial states: X0 - tol ≤ x(0) ≤ X0 + tol
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.X0[st])
                if !isnan(ocp.b.X0_tol[st])
                    # Add inequality constraint with tolerance
                    @constraint(OCPForm.mdl, ocp.b.X0[st] - ocp.b.X0_tol[st] <= OCPForm.mdl[:x][1, st] <= ocp.b.X0[st] + ocp.b.X0_tol[st])
                else
                    # No tolerance specified, use equality
                    @constraint(OCPForm.mdl, OCPForm.mdl[:x][1, st] == ocp.b.X0[st])
                end
            end
        end
    end

    # Apply final state constraints
    if !ocp.s.XFslack
        # Fixed final states: x(tf) = XF
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.XF[st])
                fix(OCPForm.mdl[:x][end, st], ocp.b.XF[st]; force = true)  # Fix final state
            end
        end
    else
        # Tolerance-based final states: XF - tol ≤ x(tf) ≤ XF + tol
        for st = 1:1:ocp.s.states.num
            if !isnan(ocp.b.XF[st])
                if !isnan(ocp.b.XF_tol[st])
                    # Add inequality constraint with tolerance
                    @constraint(OCPForm.mdl, ocp.b.XF[st] - ocp.b.XF_tol[st] <= OCPForm.mdl[:x][end, st] <= ocp.b.XF[st] + ocp.b.XF_tol[st])
                else
                    # No tolerance specified, use equality
                    @constraint(OCPForm.mdl, OCPForm.mdl[:x][end, st] == ocp.b.XF[st])
                end
            end
        end
    end
    
    # Parameters
    Nonparam = Any
    if isassigned(OCPForm.params, 1)
        ocp.p.params = @variable(OCPForm.mdl, Params[i = 1:size(OCPForm.params, 1), j = 1:size(OCPForm.params, 2)] in Parameter(OCPForm.params[i, j]) )
    end





    # Dynamical Constraints
    δx = Matrix{Any}(undef, ocp.s.states.pts, ocp.s.states.num)

    ColumnCounter = 1

    numXvar = CalXvar(OCPForm.IntegrationScheme)
    if numXvar > 0
        ocp.p.xvar = @variable(OCPForm.mdl, xvar[j in 1:numXvar, i in 1:ocp.s.states.num])
    end



    # Generate dynamics constraints for each time interval
    for j in 1:ocp.s.states.pts - 1  # Loop through all intervals (Np-1 total)
        # Get parameters for this time point
        if isassigned(OCPForm.params, j)
            param = ocp.p.params[j, :]  # Use provided parameters
        else
            param = Nonparam            # No parameters
        end

        # Apply integration scheme based on method chosen for this interval
        if OCPForm.IntegrationScheme[j] == :RK2
            # Runge-Kutta 2nd order: x_{k+1} = x_k + (k1 + k2)/2 * Δt
            xk1 = xvar[ColumnCounter, :]  # Intermediate state at t + Δt/2
            ColumnCounter += 1
            # k1 = f(x_k, u_k) - slope at beginning of interval
            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            # k2 = f(x_k + k1*Δt, u_k) - slope at end of interval
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            # Constraint: xk1 = x_k + k1 * Δt (intermediate point)
            @constraint(OCPForm.mdl, [i = 1:ocp.s.states.num],  xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j])
            # Average slope for RK2
            δx[j, :] = @expression(OCPForm.mdl, k1 ./ 2 .+ k2 ./ 2)
            # Integration constraint: x_{k+1} = x_k + δx * Δt
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])
        elseif OCPForm.IntegrationScheme[j] == :RK3
            xk1 = xvar[ColumnCounter,:]
            xk2 = xvar[ColumnCounter+1,:]
            ColumnCounter += 2

            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            k3 = @expression(OCPForm.mdl, OCPForm.dx[j](xk2, ocp.p.u[j, :], param))
            
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk2[i] - ocp.p.x[j, i] == -k1[i] * OCPForm.TInt[j] + 2*k2[i] * OCPForm.TInt[j])

            δx[j, :] = @expression(OCPForm.mdl, (k1 .+  4*k2 .+ k3) ./ 6 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :RK4
            xk1 = xvar[ColumnCounter,:]
            xk2 = xvar[ColumnCounter+1,:]
            xk3 = xvar[ColumnCounter+2,:]
            ColumnCounter += 3

            k1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            k2 = @expression(OCPForm.mdl, OCPForm.dx[j](xk1, ocp.p.u[j, :], param))
            k3 = @expression(OCPForm.mdl, OCPForm.dx[j](xk2, ocp.p.u[j, :], param))
            k4 = @expression(OCPForm.mdl, OCPForm.dx[j](xk3, ocp.p.u[j, :], param))
            
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk1[i] - ocp.p.x[j, i] == k1[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk2[i] - ocp.p.x[j, i] == k2[i] * OCPForm.TInt[j]/2)
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], xk3[i] - ocp.p.x[j, i] == k3[i] * OCPForm.TInt[j])

            δx[j, :] = @expression(OCPForm.mdl, (k1 .+  2*k2 .+ 2*k3 .+ k4) ./ 6 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :RK1
            δx[j, :] = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        elseif OCPForm.IntegrationScheme[j] == :bkwEuler
            δx[j, :] = @expression(OCPForm.mdl, OCPForm.dx[j + 1](ocp.p.x[j + 1, :], ocp.p.u[j + 1, :], param))
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] == ocp.p.x[j, i] + δx[j, i] * OCPForm.TInt[j] )

        elseif OCPForm.IntegrationScheme[j] == :trapezoidal
            δx1 = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :], param))
            δx2 = @expression(OCPForm.mdl, OCPForm.dx[j + 1](ocp.p.x[j + 1, :], ocp.p.u[j + 1, :], param))
            δx[j, :] = @expression(OCPForm.mdl, δx1 ./ 2 .+ δx2 ./ 2 )
            @constraint(OCPForm.mdl, [i=1:ocp.s.states.num], ocp.p.x[j + 1, i] - ocp.p.x[j, i] == δx[j, i] * OCPForm.TInt[j])

        end

    end

    ocp.p.δx = δx
  
    
    # Inner States Constraints
    for j in 1:ocp.s.states.pts 
        if isassigned(OCPForm.params, j)
            param = ocp.p.params[j, :]
        else
            param = Nonparam
        end

        if !isnothing(OCPForm.cons[j])
            @constraints(OCPForm.mdl, begin OCPForm.cons[j](OCPForm.mdl[:x][j, :], OCPForm.mdl[:u][j, :], param) >= 0 end)
        end
    end

    ocp.p.tV = [0.0; cumsum(OCPForm.TInt)]
    ocp.f = OCPForm
    return nothing

end

