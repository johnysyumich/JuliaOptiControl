"""
CartPole Swing-Up Control Example

This example demonstrates optimal control of a cart-pole system using JuliaOptimalControl.
The goal is to swing up the pole from hanging down to the upright position while keeping
the cart near the center.

Problem Setup:
- Initial state: Pole hanging down (θ = π), cart at center
- Final state: Pole upright (θ = 0), cart at center, minimal velocities
- Constraints: Cart position bounds, control force limits

This is a classic benchmark problem in optimal control that demonstrates:
- Nonlinear dynamics
- Unstable equilibrium control
- Energy-based swing-up strategy
"""

using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra

# Include the dynamics model
include("cartpole_dynamics.jl")

"""
    setup_cartpole_problem()

Set up the cart-pole optimal control problem with state bounds, initial/final conditions,
and dynamics.
"""
function setup_cartpole_problem()
    println("Setting up CartPole swing-up problem...")

    # Problem dimensions
    num_states = 4     # [x, ẋ, θ, θ̇]
    num_controls = 1   # [F]

    # Initial conditions: pole hanging down, cart at center
    X0 = [0.0,    # Cart position at center
          0.0,    # Cart velocity = 0
          π,      # Pole hanging down (π radians from vertical)
          0.0]    # Pole angular velocity = 0

    # Final conditions: pole approximately upright, cart near center (relaxed)
    XF = [NaN,    # Cart position (free)
          NaN,    # Cart velocity (free)
          0.0,    # Pole upright (0 radians from vertical)
          NaN]    # Pole angular velocity (free)

    # State bounds
    XL = [-3.0,   # Cart position lower bound [m]
          -5.0,   # Cart velocity lower bound [m/s]
          -3π,    # Pole angle lower bound [rad]
          -10.0]  # Pole angular velocity lower bound [rad/s]

    XU = [3.0,    # Cart position upper bound [m]
          5.0,    # Cart velocity upper bound [m/s]
          3π,     # Pole angle upper bound [rad]
          10.0]   # Pole angular velocity upper bound [rad/s]

    # Control bounds (force applied to cart)
    CL = [-20.0]  # Maximum force in negative direction [N]
    CU = [20.0]   # Maximum force in positive direction [N]

    # Create the OCP
    ocp = defineOCP(
        numStates=num_states,
        numControls=num_controls,
        X0=X0,
        XF=XF,
        XL=XL,
        XU=XU,
        CL=CL,
        CU=CU
    )

    # Define meaningful state and control names
    defineStates!(ocp, [:cart_pos, :cart_vel, :pole_angle, :pole_vel])
    defineControls!(ocp, [:force])

    println("CartPole problem setup complete!")
    return ocp
end

"""
    configure_solver(ocp; time_horizon=3.0, num_points=51)

Configure the optimal control formulation with dynamics, cost function, and solver settings.
"""
function configure_solver(ocp; time_horizon=3.0, num_points=51)
    println("Configuring solver...")

    # Configure the OCP formulation
    formulation = ConfigurePredefined(ocp,
        tf=time_horizon,              # Time horizon [s]
        Np=num_points,                # Number of discretization points
        IntegrationScheme=:RK4,       # 4th-order Runge-Kutta integration
        dx=cartpole_dynamics,         # Dynamics function
        expr=cartpole_cost            # Cost function
    )

    # Build the optimization problem
    OCPdef!(ocp, formulation)

    # Set up the objective function (minimize integral of stage cost)
    obj = ExprIntegral(ocp)
    @objective(ocp.f.mdl, Min, obj)

    println("Solver configuration complete!")
    return formulation
end

"""
    solve_cartpole(ocp)

Solve the cart-pole optimal control problem and return the solution.
"""
function solve_cartpole(ocp)
    println("Solving CartPole problem...")

    # Solve the optimization problem
    @time OptSolve!(ocp)

    # Check solution status
    if ocp.r.Status == :Optimal
        println("Solution found! Status: $(ocp.r.Status)")
        println("Objective value: $(round(ocp.r.Objval, digits=4))")
        println("Solve time: $(round(ocp.r.TSolve, digits=3)) seconds")
        println("Iterations: $(ocp.r.IterNum)")
    else
        println("Optimization failed with status: $(ocp.r.Status)")
    end

    return ocp.r.Status == :Optimal
end

"""
    visualize_results(ocp; save_plots=true)

Create comprehensive visualizations of the cart-pole solution including:
- State trajectories over time
- Control input over time
- Phase portraits
- Animation frames for creating GIFs
"""
function visualize_results(ocp; save_plots=true)
    println("Creating visualizations...")

    # Extract solution data
    time = ocp.r.Tst
    states = ocp.r.X
    controls = ocp.r.U

    # Extract individual trajectories
    cart_pos = states[:, 1]
    cart_vel = states[:, 2]
    pole_angle = states[:, 3]
    pole_vel = states[:, 4]
    force = controls[:, 1]

    # Create comprehensive plot layout
    p = plot(layout=(2, 3), size=(1200, 800))

    # Plot 1: Cart position vs time
    plot!(p[1], time, cart_pos,
          title="Cart Position",
          xlabel="Time [s]",
          ylabel="Position [m]",
          linewidth=2,
          color=:blue,
          legend=false)

    # Plot 2: Pole angle vs time (convert to degrees for readability)
    plot!(p[2], time, rad2deg.(pole_angle),
          title="Pole Angle",
          xlabel="Time [s]",
          ylabel="Angle [°]",
          linewidth=2,
          color=:red,
          legend=false)
    hline!(p[2], [0], color=:black, linestyle=:dash, alpha=0.5)

    # Plot 3: Velocities vs time
    plot!(p[3], time, cart_vel,
          title="Velocities",
          xlabel="Time [s]",
          ylabel="Velocity",
          linewidth=2,
          color=:blue,
          label="Cart [m/s]")
    plot!(p[3], time, pole_vel,
          linewidth=2,
          color=:orange,
          label="Pole [rad/s]")

    # Plot 4: Control force vs time
    plot!(p[4], time, force,
          title="Control Force",
          xlabel="Time [s]",
          ylabel="Force [N]",
          linewidth=2,
          color=:green,
          legend=false)
    hline!(p[4], [0], color=:black, linestyle=:dash, alpha=0.5)

    # Plot 5: Phase portrait (position vs velocity)
    plot!(p[5], cart_pos, cart_vel,
          title="Cart Phase Portrait",
          xlabel="Position [m]",
          ylabel="Velocity [m/s]",
          linewidth=2,
          color=:blue,
          legend=false)
    scatter!(p[5], [cart_pos[1]], [cart_vel[1]], color=:green, markersize=6, label="Start")
    scatter!(p[5], [cart_pos[end]], [cart_vel[end]], color=:red, markersize=6, label="End")

    # Plot 6: Pole phase portrait
    plot!(p[6], rad2deg.(pole_angle), pole_vel,
          title="Pole Phase Portrait",
          xlabel="Angle [°]",
          ylabel="Angular Velocity [rad/s]",
          linewidth=2,
          color=:red,
          legend=false)
    scatter!(p[6], [rad2deg(pole_angle[1])], [pole_vel[1]], color=:green, markersize=6)
    scatter!(p[6], [rad2deg(pole_angle[end])], [pole_vel[end]], color=:red, markersize=6)

    # Add overall title
    plot!(p, plot_title="CartPole Swing-Up Optimal Control Solution")

    # Display the plot
    display(p)

    # Save the plot if requested
    if save_plots
        savefig(p, "figures/cartpole_solution.png")
        println("Saved comprehensive plot as 'figures/cartpole_solution.png'")
    end

    return p
end

"""
    create_animation_frames(ocp; sample_frames=10)

Create animation frames showing the cart-pole motion for visualization.
"""
function create_animation_frames(ocp; sample_frames=10)
    println("Creating animation frames...")

    # Extract solution data
    time = ocp.r.Tst
    states = ocp.r.X

    # Sample frames evenly from the solution
    n_total = length(time)
    frame_indices = round.(Int, range(1, n_total, length=sample_frames))

    # Animation parameters
    cart_width = 0.3
    cart_height = 0.2
    pole_length = POLE_LENGTH

    frames = []

    for (frame_num, idx) in enumerate(frame_indices)
        # Current state
        x_cart = states[idx, 1]
        theta = states[idx, 3]

        # Calculate pole end position
        x_pole = x_cart + pole_length * sin(theta)
        y_pole = cart_height/2 + pole_length * cos(theta)

        # Create frame
        p = plot(xlims=(-2.5, 2.5), ylims=(-0.3, 1.2),
                aspect_ratio=:equal,
                title="CartPole Animation - Frame $frame_num/$sample_frames",
                xlabel="Position [m]",
                ylabel="Height [m]",
                legend=false,
                grid=true)

        # Draw ground
        hline!([0], color=:black, linewidth=3)

        # Draw cart
        cart_x = [x_cart - cart_width/2, x_cart + cart_width/2,
                  x_cart + cart_width/2, x_cart - cart_width/2, x_cart - cart_width/2]
        cart_y = [0, 0, cart_height, cart_height, 0]
        plot!(cart_x, cart_y, color=:blue, linewidth=2, fill=true, fillalpha=0.3)

        # Draw pole
        plot!([x_cart, x_pole], [cart_height/2, y_pole], color=:red, linewidth=4)

        # Draw pole mass as circle
        scatter!([x_pole], [y_pole], color=:red, markersize=8)

        # Add time annotation
        annotate!(1.8, 1.0, text("t = $(round(time[idx], digits=2)) s", 12))

        push!(frames, p)
    end

    println("Created $sample_frames animation frames")
    return frames
end

"""
    main()

Main function to run the complete cart-pole example.
"""
function main()
    println("CartPole Swing-Up Optimal Control Example")
    println("=" ^ 50)

    # Set up the problem
    ocp = setup_cartpole_problem()

    # Configure solver
    formulation = configure_solver(ocp, time_horizon=5.0, num_points=31)

    # Solve the problem
    success = solve_cartpole(ocp)

    if success
        # Visualize results
        visualize_results(ocp, save_plots=true)

        # Create animation frames
        frames = create_animation_frames(ocp, sample_frames=6)

        # Save a few key frames
        display(frames[1])              # Initial frame
        display(frames[length(frames)÷2]) # Middle frame
        display(frames[end])            # Final frame

        println("\nCartPole example completed successfully!")
        println("Summary:")
        println("   - Initial: Pole hanging down (180°)")
        println("   - Final: Pole upright (0°)")
        println("   - Time horizon: $(formulation.tf) seconds")
        println("   - Discretization points: $(formulation.Np)")
        println("   - Integration scheme: $(formulation.IntegrationScheme[1])")

        return ocp, frames
    else
        println("CartPole example failed to find solution")
        return nothing, nothing
    end
end

