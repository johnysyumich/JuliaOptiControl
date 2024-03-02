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

function RetrieveSolveStatus(status::MOI.TerminationStatusCode)
    SolvingStatus = [:Optimal, :UserLimit, :InFeasible]
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