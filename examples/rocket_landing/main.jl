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
using Interpolations
cd(@__DIR__)  # Set working directory to script location
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

    # Create comprehensive plot layout with enhanced styling
    p = plot(layout=(2, 3), size=(1400, 900),
             background_color=:white, foreground_color=:black)

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
    create_landing_animation(ocp; dt=0.1)

Create animation frames showing the rocket landing sequence with smooth interpolation.
"""
function create_landing_animation(ocp; dt::Float64=0.1)
    println("Creating landing animation with Δt = $dt s...")

    # Extract solution data
    time = ocp.r.Tst
    states = ocp.r.X
    controls = ocp.r.U
    @assert length(time) == size(states, 1) == size(controls, 1) "time, states, and controls length mismatch"

    # Create uniform time grid for smooth animation
    new_time = collect(time[1]:dt:time[end])

    # Interpolate states and controls on uniform grid
    nstate = size(states, 2)
    ncontrol = size(controls, 2)
    new_states = Array{Float64}(undef, length(new_time), nstate)
    new_controls = Array{Float64}(undef, length(new_time), ncontrol)

    # Interpolate each dimension independently
    for j in 1:nstate
        itp = Interpolations.linear_interpolation(time, states[:, j])
        new_states[:, j] = itp.(new_time)
    end
    for j in 1:ncontrol
        itp = Interpolations.linear_interpolation(time, controls[:, j])
        new_controls[:, j] = itp.(new_time)
    end

    # Rocket visualization parameters
    rocket_width = 15.0
    rocket_height = 30.0
    max_alt = maximum(states[:, 1])

    # Prepare trajectory trail
    trail_h = new_states[:, 1]
    trail_positions = zeros(length(new_time))  # x-position for trail (rocket is vertical)

    frames = Vector{Plots.Plot{Plots.GRBackend}}()

    for (frame_num, t) in enumerate(new_time)
        # Current state
        h = new_states[frame_num, 1]      # Altitude
        v = new_states[frame_num, 2]      # Velocity
        T = new_controls[frame_num, 1]    # Thrust

        # Normalize thrust for visualization effects
        thrust_viz = T / MAX_THRUST
        thrust_normalized = max(0.0, min(1.0, thrust_viz))

        # Create simple frame
        p = plot(xlims=(-60, 60), ylims=(-10, max_alt * 1.05),
                aspect_ratio=:equal,
                title="Rocket Landing Simulation",
                xlabel="Horizontal Position [m]",
                ylabel="Altitude [m]",
                legend=false,
                grid=true,
                background_color=:white,
                xticks=[-60, 0, 60])

        # Draw ground
        hline!([0], color=:brown, linewidth=3)

        # Simple rocket body
        rocket_x_center = 0.0

        # Main body (simple rectangle)
        body_x = [rocket_x_center - rocket_width/2, rocket_x_center + rocket_width/2,
                  rocket_x_center + rocket_width/2, rocket_x_center - rocket_width/2,
                  rocket_x_center - rocket_width/2]
        body_y = [h, h, h + rocket_height, h + rocket_height, h]
        plot!(body_x, body_y, color=:gray, linewidth=2, fill=true, fillalpha=0.7)

        # Simple nose cone
        nose_x = [rocket_x_center - rocket_width/3, rocket_x_center, rocket_x_center + rocket_width/3,
                  rocket_x_center - rocket_width/3]
        nose_y = [h + rocket_height, h + rocket_height + 6, h + rocket_height, h + rocket_height]
        plot!(nose_x, nose_y, color=:darkgray, linewidth=2, fill=true, fillalpha=0.7)

        # Side annotations (completely outside the main plot area)
        annotate!(200, max_alt * 0.9, text("t = $(round(t, digits=2)) s", 10, :black))
        annotate!(200, max_alt * 0.85, text("h = $(round(h, digits=1)) m", 10, :black))
        annotate!(200, max_alt * 0.8, text("v = $(round(v, digits=1)) m/s", 10, :black))

        push!(frames, p)
    end

    println("Created $(length(new_time)) animation frames.")
    return frames
end

function save_frames_as_gif(frames; fps::Int=30, path::AbstractString="figures/rocket_landing_animation.gif")
    mkpath(dirname(path))
    anim = Animation()
    for f in frames
        frame(anim, f)
    end
    gif(anim, path, fps=fps)
    println("Saved GIF to '$path' with $(length(frames)) frames at $fps fps.")
    return path
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

        # Create animation frames with smooth interpolation
        frames = create_landing_animation(ocp, dt=0.2)  # 0.05s for very smooth animation
        save_frames_as_gif(frames; fps=30, path="figures/rocket_landing_animation.gif")

        # Display key frames for preview
        display(frames[1])              # Initial frame
        display(frames[length(frames)÷3]) # Early descent
        display(frames[2*length(frames)÷3]) # Mid descent
        display(frames[end])            # Landing frame

        println("\nRocket Landing example completed successfully!")
        println("Summary:")
        println("   - Landing maneuver completed")
        println("   - Fuel-optimal trajectory computed")
        println("   - Smooth animation with $(length(frames)) frames created")
        println("   - Time horizon: $(formulation.tf) seconds")
        println("   - Discretization points: $(formulation.Np)")
        println("   - Integration scheme: $(formulation.IntegrationScheme[1])")

        return ocp, frames
    else
        println("Rocket Landing example failed to find solution")
        return nothing, nothing
    end
end

main()

nothing

