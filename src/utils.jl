"""
Solve the optimal control problem and extract results.
"""
function OptSolve!(ocp::OCP)
    JuMP.optimize!(ocp.f.mdl)
    ocp.r.TerminalStatus = termination_status(ocp.f.mdl)
    ocp.r.Status = RetrieveSolveStatus(ocp.r.TerminalStatus)
    ocp.r.Objval = objective_value(ocp.f.mdl)
    ocp.r.TSolve = solve_time(ocp.f.mdl)
    ocp.r.IterNum = barrier_iterations(ocp.f.mdl)
    ocp.r.EvalNum = ocp.r.EvalNum + 1
    GetOptimizeValue!(ocp)
    return nothing
end



"""
Compute numerical integration of cost expressions over time horizon.
"""
function ExprIntegral(ocp::OCP)
    cost = @expression(ocp.f.mdl, 0)  # Initialize integral to zero
    Nonparam = Any
    for j in 2:ocp.f.Np  # Start from j=2 (skip initial point)
        if !isnothing(ocp.f.expr[j])  # Check if cost expression exists at this point
            # Get parameters for this time point
            if isassigned(ocp.f.params, j)
                param = ocp.p.params[j, :]
            else
                param = Nonparam  # No parameters
            end

            # Apply integration rule based on scheme used
            if ocp.f.IntegrationScheme[j - 1] ∈ [:RK1, :RK2, :RK3, :RK4]
                # Runge-Kutta methods: use left endpoint evaluation
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j - 1, :], ocp.p.u[j - 1, :], param) * ocp.f.TInt[j - 1])
            elseif ocp.f.IntegrationScheme[j - 1] == :trapezoidal
                # Trapezoidal rule: use right endpoint with averaged control
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j, :], (ocp.p.u[j - 1, :] + ocp.p.u[j, :]) / 2, param) * ocp.f.TInt[j - 1])
            elseif ocp.f.IntegrationScheme[j - 1] == :bkwEuler
                # Backward Euler: use right endpoint evaluation
                cost = @expression(ocp.f.mdl, cost + ocp.f.expr[j](ocp.p.x[j, :], ocp.p.u[j, :], param) * ocp.f.TInt[j - 1])
            end
        end
    end
    return cost
end




"""
Convert MOI termination status to simplified status symbol.
"""
function RetrieveSolveStatus(status::MOI.TerminationStatusCode)
    SolvingStatus = [:Optimal, :UserLimit, :Infeasible]
    OptimalList = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
    LimitList = [MOI.ITERATION_LIMIT, MOI.TIME_LIMIT, MOI.NODE_LIMIT, MOI.SOLUTION_LIMIT, MOI.MEMORY_LIMIT, MOI.OBJECTIVE_LIMIT, MOI.NORM_LIMIT, MOI.OTHER_LIMIT ]

    if status ∈ OptimalList
        return SolvingStatus[1]
    elseif status ∈ LimitList
        return SolvingStatus[2]
    else
        return SolvingStatus[3]
    end
end

"""
Extract optimal variable values and store in results structure.
"""
function GetOptimizeValue!(ocp::OCP)
    if ocp.r.Status ∈ [:Optimal, :UserLimit]
        ocp.r.X = value.(ocp.p.x)
        ocp.r.U = value.(ocp.p.u)
        if ocp.f.tfDV
            ocp.r.Tst = value.(ocp.p.tV)
            ocp.r.dt = value.(ocp.f.TInt)
        else
            ocp.r.Tst = ocp.p.tV
            ocp.r.dt = ocp.f.TInt
        end
        dfs = ResultsToDataFrame(ocp)
        if ocp.s.InternalLogging
            push!(ocp.r.Dfs, dfs)
        else
            ocp.r.Dfs = Vector{DataFrame}()
            push!(ocp.r.Dfs, dfs)
        end
    else
        @warn "Solution is infeasible"
    end
end

"""
Convert optimization results to DataFrame format.
"""
function ResultsToDataFrame(ocp::OCP)
    dfs = DataFrame()
    dfs[!, :t] = ocp.r.Tst
    for st in 1:ocp.s.states.num
        dfs[!, ocp.s.states.name[st]] = ocp.r.X[:,st]
    end
    for ctr in 1:ocp.s.control.num
        dfs[!, ocp.s.control.name[ctr]] = ocp.r.U[:,ctr]
    end
    return dfs
end


"""
Create an empty OCPFormulation structure.
"""
function CreateEmptyFormulation()::OCPFormulation
    return OCPFormulation{Float64}()
end

"""
Remove constraint or expression at specified index.
"""
function DeleteElement(ConsExpr, index)
    # Check if the vector can hold nothing
    if eltype(ConsExpr) >: Nothing || eltype(ConsExpr) == Any
        ConsExpr[index] = nothing
    else
        # For vectors that can't hold nothing, remove the element
        deleteat!(ConsExpr, index)
    end
end


"""
Set warm start values using previous solution for faster convergence.
"""
function WarmStart(ocp::OCP)
    flag = false
    # Check if solver supports warm starting
    try get_attribute(ocp.f.mdl, "warm_start_init_point")
      if get_attribute(ocp.f.mdl, "warm_start_init_point") == "yes"
        flag = true  # Warm start enabled
      else
        flag = false  # Warm start disabled
      end
    catch
      flag = false  # Solver doesn't support warm start attribute
    end

    if flag == true && !isempty(ocp.r.X) && !isempty(ocp.r.U) && size(ocp.r.X, 1) > 1
      # Initialize state variables: [X0; previous_trajectory[2:end]]
      set_start_value.(ocp.p.x, [ocp.b.X0'; ocp.r.X[2:end, :]])
      # Initialize control variables: shift previous solution + extend last value
      set_start_value.(ocp.p.u, [ocp.r.U[2:end, :]; ocp.r.U[end, :]'])
    end
    return nothing
end
  

"""
Update initial conditions and apply warm start for MPC.
"""
function UpdateX0!(ocp::OCP, X0)
    ocp.b.X0 = X0  # Update stored initial conditions
    # Fix all initial state variables to new values
    for i = 1:1:ocp.s.states.num
        fix(ocp.p.x[1, i], ocp.b.X0[i]; force = true)  # Force update constraint
    end
    WarmStart(ocp)  # Apply warm start with previous solution
end
