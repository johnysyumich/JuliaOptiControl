# Simple rocket landing with fixed mass
using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra

# Simple rocket dynamics (fixed mass)
function rocket_dynamics_fixed(states, controls, params)
    h, v = states
    T = controls[1]

    g = 9.81
    m = 1000.0  # Fixed mass

    h_dot = v
    v_dot = T/m - g

    return [h_dot, v_dot]
end

# Simple cost function
function rocket_cost(states, controls, params)
    h, v = states
    T = controls[1]

    return 0.001 * T  # Minimize fuel (thrust)
end

println("🚀 Testing Fixed-Mass Rocket Landing...")

# Problem setup
X0 = [200.0, -10.0]             # Start: 200m high, falling
XF = [NaN, NaN]                 # End: free final state
XL = [0.0, -50.0]               # Bounds
XU = [500.0, 20.0]
CL = [0.0]                      # No negative thrust
CU = [15000.0]                  # Max thrust

ocp = defineOCP(numStates=2, numControls=1, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU)
defineStates!(ocp, [:altitude, :velocity])
defineControls!(ocp, [:thrust])

println("✅ OCP defined!")

# Configure solver
formulation = ConfigurePredefined(ocp,
    tf=15.0,
    Np=21,
    IntegrationScheme=:RK4,
    dx=rocket_dynamics_fixed,
    expr=rocket_cost
)

OCPdef!(ocp, formulation)

# Add simple final constraints for landing
h = ocp.p.x[:, 1]
v = ocp.p.x[:, 2]

# Final altitude should be low
@constraint(ocp.f.mdl, h[end] <= 5.0)
# Final velocity should be small downward
@constraint(ocp.f.mdl, v[end] >= -5.0)
@constraint(ocp.f.mdl, v[end] <= 1.0)

obj = ExprIntegral(ocp)
# Add terminal cost for landing accuracy
terminal_cost = 500 * h[end]^2 + 100 * (v[end] + 1)^2
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
    hline!(p[1], [0], color=:red, linestyle=:dash, alpha=0.5)

    plot!(p[2], ocp.r.Tst, ocp.r.X[:, 2], title="Velocity", label="v [m/s]")
    hline!(p[2], [0], color=:black, linestyle=:dash, alpha=0.5)

    plot!(p[3], ocp.r.Tst, ocp.r.U[:, 1]./1000, title="Thrust", label="T [kN]")

    # Phase portrait
    plot!(p[4], ocp.r.X[:, 2], ocp.r.X[:, 1], title="Phase Portrait",
          xlabel="Velocity [m/s]", ylabel="Altitude [m]", label="Trajectory")
    scatter!(p[4], [ocp.r.X[1, 2]], [ocp.r.X[1, 1]], color=:green, markersize=6, label="Start")
    scatter!(p[4], [ocp.r.X[end, 2]], [ocp.r.X[end, 1]], color=:red, markersize=6, label="Land")

    display(p)
    savefig(p, "rocket_fixed_mass.png")
    println("💾 Plot saved as rocket_fixed_mass.png")

    # Landing analysis
    println("🎯 Landing results:")
    println("   Final altitude: $(round(ocp.r.X[end, 1], digits=2)) m")
    println("   Final velocity: $(round(ocp.r.X[end, 2], digits=2)) m/s")
    println("   Landing time: $(round(ocp.r.Tst[end], digits=1)) s")
else
    println("❌ Solution failed!")
end