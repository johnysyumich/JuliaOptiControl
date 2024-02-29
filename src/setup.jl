include("types.jl")

function defineOCP(;
    numStates::Int = 0,
    numControls::Int = 0,
    X0 = fill(NaN,numStates),
    XF = fill(NaN,numStates),
    XL = fill(NaN,numStates),
    XU = fill(NaN,numStates),
    CL = fill(NaN,numControls),
    CU = fill(NaN,numControls))::OCP
    ocp = OCP()
    if numStates < 0 error("States number must be larger than 0") end
    if numControls < 0 error("Control number must be larger than 0") end
    if length(X0) != numStates error("Number of states do not match X0") end
    if length(XF) != numStates error("Number of states do not match XF") end
    if length(XL) != numStates error("Number of states do not match XL") end
    if length(XU) != numStates error("Number of states do not match XU") end
    if length(CL) != numControls error("Number of controls do not match CL") end
    if length(CU) != numControls error("Number of controls do not match CU") end

    ocp.s.states.num    = numStates
    ocp.s.control.num   = numControls
    ocp.b.X0            = X0
    ocp.b.XL            = XL
    ocp.b.XU            = XU
    ocp.b.CL            = CL
    ocp.b.CU            = CU


    return ocp
end

function defineStates!(ocp::OCP, states::Vector{Symbol})
    if length(states) != ocp.s.states.num error("Wrong number of state variable") end
    ocp.s.states.name = states
    return nothing
end

function defineControls!(ocp::OCP, controls::Vector{Symbol})
    if length(controls) != ocp.s.control.num error("Wrong number of state variable") end
    ocp.s.control.name = controls
    return nothing
end


function ConfigurePredefined(ocp::OCP; kwargs...)
    # all inputs here: tfDV, tf, Np, Integration Scheme
    # Integration scheme: bkwEuler, trapezoidal
    # Trajectory methods: SingleShooting Collocation
    OCPForm = OCPFormulation{Float64}()
    OCPForm.mdl = JuMP.Model()
    kw = Dict(kwargs)
    if haskey(kw, :tfDV)
        OCPForm.tfDV = get(kw, :tfDV, 0)
    end
    if haskey(kw, :tf)
        if OCPForm.tfDV == true
            @warn "Conflict intructions to put tf as design variable, change to fixed horizon"
            OCPForm.tfDV = false
        end
        OCPForm.tf = get(kw, :tf, 0)
        if OCPForm.tf <= ocp.b.tfMin || OCPForm.tf >= ocp.b.tfMax
            error("Please make sure tf ∈ [$(ocp.b.tfMin), $(ocp.b.tfMax)]")
        end
    else
        if OCPForm.tfDV == false
            error("Either treat tf as design variable or use fixed time horizon")
        end
    end


    if haskey(kw, :TrajMethod)
        TrajectoryMethod = get(kw, :TrajMethod, 0)
        if TrajectoryMethod ∈ [:SingleShooting, :Collocation]
            ocp.s.TrajMethod = TrajectoryMethod
        end
    end

    if haskey(kw, :IntegrationScheme)
        IntegrationScheme = get(kw, :IntegrationScheme, 0)
        if IntegrationScheme ∈ [:bkwEuler, :trapezoidal]
            OCPForm.IntegrationScheme = IntegrationScheme
        else
            @warn "$IntegrationScheme is not implemented, use default (bkwEuler)"
        end
    end

    if haskey(kw, :Np)
        Np = get(kw, :Np, 0)
        NpType = typeof(Np)
        if NpType != Int64 
            if NpType <: Real
                @warn "Round Np to nearest integer"
                OCPForm.Np = Int(round(a))
            else
                error("Wrong type of number of points, should be Int64")
            end
        else
            OCPForm.Np = Np
        end
        if OCPForm.Np <= 0
            error("Make sure Np is an integer and is larger than 0")
        end
    else
        error("No number of point input")
    end
    if OCPForm.tfDV == false && OCPForm.IntegrationScheme ∈ [:bkwEuler, :trapezoidal]
        OCPForm.TInt = OCPForm.tf / (OCPForm.Np - 1) .* ones(OCPForm.Np - 1)
    elseif OCPForm.tfDV == true && OCPForm.IntegrationScheme ∈ [:bkwEuler, :trapezoidal]
        OCPForm.tf = @variable(OCPForm.mdl, ocp.b.tfMin <= tf <= ocp.b.tfMax)
        OCPForm.TInt = @expression(OCPForm.mdl, [idx = 1:OCPForm.Np - 1], OCPForm.tf / (OCPForm.Np - 1))
    end
    return OCPForm
end