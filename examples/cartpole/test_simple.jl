# Simple test for CartPole dynamics
using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra

# Simple cart-pole dynamics
function cartpole_dynamics(states, controls, params)
    x, x_dot, theta, theta_dot = states
    F = controls[1]

    # System parameters
    M = 1.0    # Cart mass
    m = 0.1    # Pole mass
    L = 0.5    # Pole half-length
    g = 9.81   # Gravity

    # Simplified dynamics
    sin_theta = sin(theta)
    cos_theta = cos(theta)

    denominator = M + m * sin_theta^2

    x_ddot = (F + m * L * theta_dot^2 * sin_theta - m * g * sin_theta * cos_theta) / denominator
    theta_ddot = (-F * cos_theta - m * L * theta_dot^2 * sin_theta * cos_theta + (M + m) * g * sin_theta) / (L * denominator)

    return [x_dot, x_ddot, theta_dot, theta_ddot]
end

# Simple cost function
function cartpole_cost(states, controls, params)
    x, x_dot, theta, theta_dot = states
    F = controls[1]

    return x^2 + 100*theta^2 + 0.1*x_dot^2 + theta_dot^2 + 0.01*F^2
end

println("🎯 Testing CartPole setup...")

# Problem setup
X0 = [0.0, 0.0, π, 0.0]      # Start: pole hanging down
XF = [0.0, 0.0, 0.0, 0.0]    # End: pole upright
XL = [-3.0, -5.0, -3π, -10.0]
XU = [3.0, 5.0, 3π, 10.0]
CL = [-20.0]
CU = [20.0]

ocp = defineOCP(numStates=4, numControls=1, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU)
defineStates!(ocp, [:cart_pos, :cart_vel, :pole_angle, :pole_vel])
defineControls!(ocp, [:force])

println("✅ OCP defined successfully!")

# Configure solver
formulation = ConfigurePredefined(ocp,
    tf=3.0,
    Np=21,  # Smaller for faster testing
    IntegrationScheme=:RK4,
    dx=cartpole_dynamics,
    expr=cartpole_cost
)

println("✅ Formulation configured!")

# Build problem
OCPdef!(ocp, formulation)
obj = ExprIntegral(ocp)
@objective(ocp.f.mdl, Min, obj)

println("✅ Problem built!")

# Solve
println("🚀 Solving...")
@time OptSolve!(ocp)

println("✅ Solution status: $(ocp.r.Status)")

if ocp.r.Status == :Optimal
    println("📊 Objective: $(round(ocp.r.Objval, digits=4))")

    # Simple plot
    p = plot(layout=(2,2), size=(800,600))

    plot!(p[1], ocp.r.Tst, ocp.r.X[:, 1], title="Cart Position", label="x [m]")
    plot!(p[2], ocp.r.Tst, rad2deg.(ocp.r.X[:, 3]), title="Pole Angle", label="θ [°]")
    plot!(p[3], ocp.r.Tst, ocp.r.X[:, 2], title="Cart Velocity", label="ẋ [m/s]")
    plot!(p[4], ocp.r.Tst, ocp.r.U[:, 1], title="Control Force", label="F [N]")

    display(p)
    savefig(p, "cartpole_test.png")
    println("💾 Plot saved as cartpole_test.png")
else
    println("❌ Solution failed!")
end