# Simple test for Rocket Landing
using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra

# Simple rocket dynamics
function rocket_dynamics(states, controls, params)
    h, v, m = states
    T = controls[1]

    g = 9.81
    ve = 2500.0

    h_dot = v
    v_dot = T/m - g
    m_dot = -T/ve

    return [h_dot, v_dot, m_dot]
end

# Simple cost function
function rocket_cost(states, controls, params)
    h, v, m = states
    T = controls[1]

    return T  # Just minimize fuel (thrust)
end

println("🚀 Testing Rocket Landing setup...")

# Problem setup - more conservative
X0 = [500.0, -20.0, 1200.0]     # Start: 500m high, falling slowly, with fuel
XF = [NaN, NaN, NaN]             # End: free final state
XL = [0.0, -100.0, 800.0]       # Allow more mass consumption
XU = [1000.0, 50.0, 1200.0]
CL = [0.0]                       # No negative thrust
CU = [15000.0]                   # Reduced max thrust

ocp = defineOCP(numStates=3, numControls=1, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU)
defineStates!(ocp, [:altitude, :velocity, :mass])
defineControls!(ocp, [:thrust])

println("✅ OCP defined!")

# Configure solver
formulation = ConfigurePredefined(ocp,
    tf=20.0,
    Np=21,
    IntegrationScheme=:RK4,
    dx=rocket_dynamics,
    expr=rocket_cost
)

OCPdef!(ocp, formulation)

# Add simple final constraints for landing
h = ocp.p.x[:, 1]
v = ocp.p.x[:, 2]

# Final altitude should be low
@constraint(ocp.f.mdl, h[end] <= 10.0)
# Final velocity should be reasonable
@constraint(ocp.f.mdl, v[end] >= -10.0)

obj = ExprIntegral(ocp)
# Add terminal cost for landing accuracy
terminal_cost = 100 * h[end]^2 + 10 * (v[end] + 2)^2
@objective(ocp.f.mdl, Min, obj + terminal_cost)

println("✅ Problem built!")

# Solve
println("🚀 Solving...")
@time OptSolve!(ocp)

println("✅ Solution status: $(ocp.r.Status)")

if ocp.r.Status == :Optimal
    println("📊 Objective: $(round(ocp.r.Objval, digits=4))")

    # Simple plot
    p = plot(layout=(2,2), size=(800,600))

    plot!(p[1], ocp.r.Tst, ocp.r.X[:, 1], title="Altitude", label="h [m]")
    plot!(p[2], ocp.r.Tst, ocp.r.X[:, 2], title="Velocity", label="v [m/s]")
    plot!(p[3], ocp.r.Tst, ocp.r.X[:, 3], title="Mass", label="m [kg]")
    plot!(p[4], ocp.r.Tst, ocp.r.U[:, 1]./1000, title="Thrust", label="T [kN]")

    display(p)
    savefig(p, "rocket_test.png")
    println("💾 Plot saved as rocket_test.png")

    # Landing analysis
    println("🎯 Landing results:")
    println("   Final altitude: $(round(ocp.r.X[end, 1], digits=2)) m")
    println("   Final velocity: $(round(ocp.r.X[end, 2], digits=2)) m/s")
    println("   Fuel used: $(round(X0[3] - ocp.r.X[end, 3], digits=1)) kg")
else
    println("❌ Solution failed!")
end