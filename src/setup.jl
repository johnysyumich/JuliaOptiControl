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


function Configure!(ocp::OCP)