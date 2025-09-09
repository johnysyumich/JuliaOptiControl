"""
Rocket Landing Optimal Control Example

This example demonstrates optimal control of a rocket performing a vertical landing
maneuver, similar to SpaceX Falcon 9 landings. The goal is to land the rocket softly
while minimizing fuel consumption.

Problem Setup:
- Initial state: High altitude, some downward velocity, full fuel
- Final state: Ground level, soft landing velocity, fuel remaining
- Objective: Minimize fuel consumption
- Constraints: Thrust limits, no ground collision

This demonstrates:
- Fuel-optimal control
- Inequality path constraints
- Terminal constraints
- Mass-varying dynamics
"""

using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra
using Statistics

# Include the dynamics model
include("rocket_dynamics.jl")

"""
    setup_rocket_problem()

Set up the rocket landing optimal control problem with state bounds, initial/final
conditions, and dynamics.
"""
function setup_rocket_problem()
    println("Setting up Rocket Landing problem...")

    # Problem dimensions (simplified: fixed mass)
    num_states = 2     # [h, v] - altitude and velocity
    num_controls = 1   # [T] - thrust

    # Initial conditions: high altitude, falling
    initial_altitude = 500.0   # Start at 500m altitude (more reasonable)
    initial_velocity = -30.0   # Falling at 30 m/s

    X0 = [initial_altitude,   # Altitude [m]
          initial_velocity]   # Velocity [m/s] (negative = downward)

    # Final conditions: near ground, soft landing
    XF = [NaN,    # Final altitude (optimized, but should be small)
          NaN]    # Final velocity (optimized, but should be small)

    # State bounds
    XL = [0.0,            # Altitude ≥ 0 (no ground penetration)
          -100.0]         # Velocity lower bound [m/s]

    XU = [1000.0,         # Altitude upper bound [m]
          50.0]           # Velocity upper bound [m/s]

    # Control bounds (thrust force)
    CL = [0.0]           # Minimum thrust (can't suck in exhaust)
    CU = [MAX_THRUST]    # Maximum thrust [N]

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
    defineStates!(ocp, [:altitude, :velocity])
    defineControls!(ocp, [:thrust])

    println("Rocket problem setup complete!")
    return ocp
end

"""
    configure_solver(ocp; time_horizon=30.0, num_points=51)

Configure the optimal control formulation with dynamics, cost function, and solver settings.
"""
function configure_solver(ocp; time_horizon=30.0, num_points=51)
    println("Configuring solver...")

    # Configure the OCP formulation
    formulation = ConfigurePredefined(ocp,
        tf=time_horizon,              # Time horizon [s]
        Np=num_points,                # Number of discretization points
        IntegrationScheme=:RK4,       # 4th-order Runge-Kutta integration
        dx=rocket_dynamics,           # Dynamics function
        expr=rocket_cost              # Stage cost function
    )

    # Build the optimization problem
    OCPdef!(ocp, formulation)

    # Set up the objective function (minimize integral of stage cost + terminal cost)
    stage_cost = ExprIntegral(ocp)

    # Add terminal cost for landing accuracy
    h_final = ocp.p.x[end, 1]  # Final altitude
    v_final = ocp.p.x[end, 2]  # Final velocity

    # Terminal cost: penalize deviation from perfect landing
    terminal_cost = 500.0 * h_final^2 + 100.0 * (v_final + 1.0)^2

    # Total objective
    @objective(ocp.f.mdl, Min, stage_cost + terminal_cost)

    println("Solver configuration complete!")
    return formulation
end

"""
    add_path_constraints(ocp)

Add path constraints to ensure safe landing trajectory.
"""
function add_path_constraints(ocp)
    println("Adding safety constraints...")

    # Extract state variables
    h = ocp.p.x[:, 1]  # Altitude
    v = ocp.p.x[:, 2]  # Velocity
    T = ocp.p.u[:, 1]  # Thrust

    # Constraint 1: Final altitude should be low (near ground)
    @constraint(ocp.f.mdl, h[end] <= 5.0)

    # Constraint 2: Final landing velocity should be reasonable
    @constraint(ocp.f.mdl, v[end] >= -8.0)
    @constraint(ocp.f.mdl, v[end] <= 2.0)

    println("Safety constraints added!")
end

"""
    solve_rocket(ocp)

Solve the rocket landing optimal control problem and return the solution.
"""
function solve_rocket(ocp)
    println("Solving Rocket Landing problem...")

    # Add path constraints
    add_path_constraints(ocp)

    # Solve the optimization problem
    @time OptSolve!(ocp)

    # Check solution status
    if ocp.r.Status == :Optimal
        println("Solution found! Status: $(ocp.r.Status)")
        println("Objective value: $(round(ocp.r.Objval, digits=4))")
        println("Solve time: $(round(ocp.r.TSolve, digits=3)) seconds")
        println("Iterations: $(ocp.r.IterNum)")

        # Landing analysis
        final_altitude = ocp.r.X[end, 1]
        final_velocity = ocp.r.X[end, 2]
        max_thrust = maximum(ocp.r.U[:, 1])
        avg_thrust = mean(ocp.r.U[:, 1])

        println("Landing analysis:")
        println("   - Final altitude: $(round(final_altitude, digits=2)) m")
        println("   - Final velocity: $(round(final_velocity, digits=2)) m/s")
        println("   - Max thrust: $(round(max_thrust/1000, digits=1)) kN")
        println("   - Avg thrust: $(round(avg_thrust/1000, digits=1)) kN")

    else
        println("Optimization failed with status: $(ocp.r.Status)")
    end

    return ocp.r.Status == :Optimal
end

"""
    visualize_results(ocp; save_plots=true)

Create comprehensive visualizations of the rocket landing solution.
"""
function visualize_results(ocp; save_plots=true)
    println("Creating visualizations...")

    # Extract solution data
    time = ocp.r.Tst
    states = ocp.r.X
    controls = ocp.r.U

    # Extract individual trajectories
    altitude = states[:, 1]
    velocity = states[:, 2]
    thrust = controls[:, 1]

    # Calculate derived quantities
    avg_mass = DRY_MASS + FUEL_MASS/2  # Use average mass
    thrust_to_weight = thrust ./ (avg_mass * GRAVITY)  # Thrust-to-weight ratio

    # Create comprehensive plot layout
    p = plot(layout=(2, 3), size=(1200, 800))

    # Plot 1: Altitude vs time
    plot!(p[1], time, altitude,
          title="Altitude Profile",
          xlabel="Time [s]",
          ylabel="Altitude [m]",
          linewidth=2,
          color=:blue,
          legend=false)
    hline!(p[1], [0], color=:red, linestyle=:dash, alpha=0.5, label="Ground")

    # Plot 2: Velocity vs time
    plot!(p[2], time, velocity,
          title="Velocity Profile",
          xlabel="Time [s]",
          ylabel="Velocity [m/s]",
          linewidth=2,
          color=:red,
          legend=false)
    hline!(p[2], [0], color=:black, linestyle=:dash, alpha=0.5)

    # Plot 3: Thrust vs time
    plot!(p[3], time, thrust./1000,  # Convert to kN
          title="Thrust Profile",
          xlabel="Time [s]",
          ylabel="Thrust [kN]",
          linewidth=2,
          color=:green,
          legend=false)

    # Plot 4: Thrust-to-weight ratio
    plot!(p[4], time, thrust_to_weight,
          title="Thrust-to-Weight Ratio",
          xlabel="Time [s]",
          ylabel="T/W [-]",
          linewidth=2,
          color=:orange,
          legend=false)
    hline!(p[4], [1], color=:black, linestyle=:dash, alpha=0.5)

    # Plot 5: Altitude vs velocity (phase portrait)
    plot!(p[5], velocity, altitude,
          title="Trajectory Phase Portrait",
          xlabel="Velocity [m/s]",
          ylabel="Altitude [m]",
          linewidth=2,
          color=:purple,
          legend=false)
    scatter!(p[5], [velocity[1]], [altitude[1]], color=:green, markersize=6, label="Start")
    scatter!(p[5], [velocity[end]], [altitude[end]], color=:red, markersize=6, label="Land")

    # Plot 6: Acceleration profile
    acceleration = thrust ./ avg_mass .- GRAVITY
    plot!(p[6], time, acceleration,
          title="Acceleration Profile",
          xlabel="Time [s]",
          ylabel="Acceleration [m/s²]",
          linewidth=2,
          color=:brown,
          legend=false)
    hline!(p[6], [0], color=:black, linestyle=:dash, alpha=0.5)

    # Add overall title
    plot!(p, plot_title="Rocket Landing Optimal Control Solution")

    # Display the plot
    display(p)

    # Save the plot if requested
    if save_plots
        savefig(p, "figures/rocket_landing_solution.png")
        println("Saved comprehensive plot as 'figures/rocket_landing_solution.png'")
    end

    return p
end

"""
    create_landing_animation(ocp; sample_frames=10)

Create animation frames showing the rocket landing sequence.
"""
function create_landing_animation(ocp; sample_frames=10)
    println("Creating landing animation frames...")

    # Extract solution data
    time = ocp.r.Tst
    states = ocp.r.X
    controls = ocp.r.U

    # Sample frames evenly from the solution
    n_total = length(time)
    frame_indices = round.(Int, range(1, n_total, length=sample_frames))

    frames = []

    for (frame_num, idx) in enumerate(frame_indices)
        # Current state
        h = states[idx, 1]      # Altitude
        v = states[idx, 2]      # Velocity
        T = controls[idx, 1]    # Thrust
        avg_mass = DRY_MASS + FUEL_MASS/2  # Average mass

        # Normalize thrust for visualization (0 to 1)
        thrust_viz = T / MAX_THRUST

        # Create frame
        max_alt = maximum(states[:, 1])
        p = plot(xlims=(-50, 50), ylims=(0, max_alt * 1.1),
                aspect_ratio=:equal,
                title="Rocket Landing - Frame $frame_num/$sample_frames",
                xlabel="Position [m]",
                ylabel="Altitude [m]",
                legend=false,
                grid=true)

        # Draw ground
        hline!([0], color=:brown, linewidth=5, label="Ground")

        # Draw rocket body (rectangle)
        rocket_width = 4
        rocket_height = 15
        rocket_x = [-rocket_width/2, rocket_width/2, rocket_width/2, -rocket_width/2, -rocket_width/2]
        rocket_y = [h, h, h + rocket_height, h + rocket_height, h]
        plot!(rocket_x, rocket_y, color=:silver, linewidth=2, fill=true, fillalpha=0.7)

        # Draw thrust plume (if thrusting)
        if T > 1000  # Only show visible thrust
            plume_length = thrust_viz * 20  # Scale plume with thrust
            plume_width = rocket_width * 0.8
            plume_x = [-plume_width/2, 0, plume_width/2, -plume_width/2]
            plume_y = [h, h - plume_length, h, h]
            plot!(plume_x, plume_y, color=:orange, linewidth=2, fill=true, fillalpha=0.8)
        end

        # Add velocity arrow
        if abs(v) > 1.0  # Only show significant velocity
            arrow_scale = min(abs(v) / 10.0, 10.0)  # Scale arrow
            arrow_end_y = h + rocket_height/2 + (v > 0 ? arrow_scale : -arrow_scale)
            plot!([0, 0], [h + rocket_height/2, arrow_end_y],
                  color=:red, linewidth=3, arrow=true)
        end

        # Add information text
        annotate!(30, max_alt * 0.9, text("t = $(round(time[idx], digits=1)) s", 10))
        annotate!(30, max_alt * 0.85, text("h = $(round(h, digits=1)) m", 10))
        annotate!(30, max_alt * 0.8, text("v = $(round(v, digits=1)) m/s", 10))
        annotate!(30, max_alt * 0.75, text("T = $(round(T/1000, digits=1)) kN", 10))
        annotate!(30, max_alt * 0.7, text("T/W = $(round(T/(avg_mass*GRAVITY), digits=1))", 10))

        push!(frames, p)
    end

    println("Created $sample_frames animation frames")
    return frames
end

"""
    main()

Main function to run the complete rocket landing example.
"""
function main()
    println("Rocket Landing Optimal Control Example")
    println("=" ^ 50)

    # Set up the problem
    ocp = setup_rocket_problem()

    # Configure solver
    formulation = configure_solver(ocp, time_horizon=25.0, num_points=41)

    # Solve the problem
    success = solve_rocket(ocp)

    if success
        # Visualize results
        visualize_results(ocp, save_plots=true)

        # Create animation frames
        frames = create_landing_animation(ocp, sample_frames=8)

        # Display key frames
        display(frames[1])              # Initial frame
        display(frames[length(frames)÷2]) # Middle frame
        display(frames[end])            # Landing frame

        println("\nRocket Landing example completed successfully!")
        println("Summary:")
        println("   - Landing maneuver completed")
        println("   - Fuel-optimal trajectory computed")
        println("   - Time horizon: $(formulation.tf) seconds")
        println("   - Discretization points: $(formulation.Np)")
        println("   - Integration scheme: $(formulation.IntegrationScheme[1])")

        return ocp, frames
    else
        println("Rocket Landing example failed to find solution")
        return nothing, nothing
    end
end

