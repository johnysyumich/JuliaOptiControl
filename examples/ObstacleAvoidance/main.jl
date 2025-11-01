"""
Vehicle Obstacle Avoidance Optimal Control Example

This example demonstrates optimal control of a vehicle navigating around a circular
obstacle using JuliaOptimalControl. The goal is to reach a target while avoiding
collision and minimizing control effort.

Problem Setup:
- Initial state: Vehicle starting behind obstacle
- Obstacle: Fixed circular obstacle in path
- Objective: Reach target while avoiding collision and minimizing control effort
- Constraints: Collision avoidance, vehicle dynamics limits

This demonstrates:
- Static obstacle avoidance
- Circular obstacle constraints
- Vehicle dynamics (bicycle model)
- Path optimization with geometric constraints
"""

using JuliaOptimalControl
using Plots
using JuMP
using LinearAlgebra
using Statistics
using Interpolations
using DataFrames
cd(@__DIR__)  # Set working directory to script location

# Include dynamics and parameters
include("bicycleModel.jl")
include("parameters.jl")

# Vehicle parameters
const VEHICLE_LENGTH = 4.5  # Vehicle length [m]
const VEHICLE_WIDTH = 2.0   # Vehicle width [m]

"""
    setup_obs_avoid_problem()

Set up the vehicle obstacle avoidance optimal control problem with bicycle model
dynamics and circular obstacle constraints.
"""
function setup_obs_avoid_problem()
    println("Setting up Vehicle Obstacle Avoidance problem...")

    # Problem dimensions: 7 states, 2 controls
    # States: [x, y, v, r, ψ, ux, δf]
    # Controls: [ax, dδf] - longitudinal acceleration, steering rate

    # State bounds
    XL = [-40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12]
    XU = [300, 20, 3, pi/5, pi/2, 15.0, pi/12]

    # Control bounds
    CL = [-2.6, -0.1]  # [accel_min, steering_rate_min]
    CU = [2.6, 0.1]    # [accel_max, steering_rate_max]

    # Initial conditions: vehicle starting position
    X0 = [-10.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0]

    # Final conditions: free (reach target area)
    XF = [NaN, NaN, NaN, NaN, NaN, NaN, NaN]

    # Create optimal control problem
    ocp = defineOCP(
        numStates=7,
        numControls=2,
        X0=X0,
        XF=XF,
        XL=XL,
        XU=XU,
        CL=CL,
        CU=CU
    )

    # Define meaningful state and control names
    defineStates!(ocp, [:x, :y, :v, :r, :ψ, :ux, :δf])
    defineControls!(ocp, [:ax, :dδf])

    println("Vehicle Obstacle Avoidance problem setup complete!")
    return ocp
end

"""
    configure_obs_avoid_formulation(ocp)

Configure the discretization and dynamics for the vehicle optimization problem.
"""
function configure_obs_avoid_formulation(ocp)
    println("Configuring solver...")

    # Configure problem formulation
    formulation = ConfigurePredefined(ocp,
        Np = 81,                          # 81 discretization points
        tfDV = false,                     # Fixed final time
        tf = 8.0,                         # 8 second time horizon
        IntegrationScheme = :RK1,         # Forward Euler
        dx = bicycleModel_expr,           # Bicycle model dynamics
        expr = bicycle_cost               # Base cost function
    )

    println("Solver configuration complete!")
    return formulation
end

"""
    add_obstacle_constraints!(ocp, formulation)

Add circular obstacle avoidance constraints and custom objective function.
"""
function add_obstacle_constraints!(ocp, formulation)
    println("Adding obstacle avoidance constraints...")

    # Build the optimization problem
    OCPdef!(ocp, formulation)

    # Obstacle parameters
    block_center_x_info = 30.0
    block_center_y_info = 2.0
    block_radius_info = 6.0

    # Add obstacle parameters as JuMP parameters
    ocp.f.mdl[:block_center_x] = @variable(ocp.f.mdl, block_center_x in Parameter(block_center_x_info))
    ocp.f.mdl[:block_center_y] = @variable(ocp.f.mdl, block_center_y in Parameter(block_center_y_info))
    ocp.f.mdl[:block_radius] = @variable(ocp.f.mdl, block_radius in Parameter(block_radius_info))

    # Set parameter values
    set_parameter_value.(ocp.f.mdl[:block_center_x], block_center_x_info)
    set_parameter_value.(ocp.f.mdl[:block_center_y], block_center_y_info)
    set_parameter_value.(ocp.f.mdl[:block_radius], block_radius_info)

    # Get state and control variables
    xpos = ocp.p.x[:, 1]  # X position
    y = ocp.p.x[:, 2]     # Y position
    ux = ocp.p.x[:, 6]    # Longitudinal velocity
    δf = ocp.p.x[:, 7]    # Front wheel angle
    ax = ocp.p.u[:, 1]    # Longitudinal acceleration
    dδf = ocp.p.u[:, 2]   # Steering rate

    # Circular obstacle avoidance constraint
    block_center_x = ocp.f.mdl[:block_center_x]
    block_center_y = ocp.f.mdl[:block_center_y]
    block_radius = ocp.f.mdl[:block_radius]

    obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts],
        block_radius^2 <= ((xpos[i] - block_center_x)^2 + (y[i] - block_center_y)^2))

    # Custom objective function
    obj = @expression(ocp.f.mdl,
        sum((0.05 * y[j]^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 +
             0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1]
            for j in 2:ocp.f.Np))

    @objective(ocp.f.mdl, Min, obj)

    println("Obstacle avoidance constraints added!")

    # Return obstacle info for visualization
    obstacle_info = (x=block_center_x_info, y=block_center_y_info, r=block_radius_info)
    return obstacle_info
end

"""
    solve_obs_avoid(ocp)

Solve the vehicle obstacle avoidance optimal control problem.
"""
function solve_obs_avoid(ocp)
    println("Solving Vehicle Obstacle Avoidance problem...")

    # Solve the optimization problem
    @time OptSolve!(ocp)

    # Check solution status
    if ocp.r.Status == :Optimal
        println("Solution found! Status: $(ocp.r.Status)")
        println("Objective value: $(round(ocp.r.Objval, digits=4))")
        println("Solve time: $(round(ocp.r.TSolve, digits=3)) seconds")
        println("Iterations: $(ocp.r.IterNum)")
        return true
    else
        println("Failed to find optimal solution. Status: $(ocp.r.Status)")
        return false
    end
end

"""
    circleShape(h, k, r)

Generate points for drawing a circle with center (h,k) and radius r.
"""
function circleShape(h, k, r)
    theta = LinRange(0, 2*pi, 500)
    return h .+ r*sin.(theta), k .+ r*cos.(theta)
end

"""
    visualize_obs_avoid_results(ocp, obstacle_info; save_plots=true)

Create comprehensive visualization of the vehicle obstacle avoidance results.
"""
function visualize_obs_avoid_results(ocp, obstacle_info; save_plots=true)
    println("Creating visualizations...")

    # Create comprehensive plot
    p = plot(layout=(2,2), size=(1000, 800))

    # Plot 1: Vehicle trajectory with obstacle
    plot!(p[1], ocp.r.X[:, 1], ocp.r.X[:, 2],
          title="Vehicle Path Planning with Obstacle",
          xlabel="X Position [m]",
          ylabel="Y Position [m]",
          linewidth=3, color=:blue, label="Vehicle Path")

    scatter!(p[1], [ocp.r.X[1,1]], [ocp.r.X[1,2]],
             color=:green, markersize=8, label="Start")
    scatter!(p[1], [ocp.r.X[end,1]], [ocp.r.X[end,2]],
             color=:red, markersize=8, label="End")

    # Add circular obstacle
    circle_x, circle_y = circleShape(obstacle_info.x, obstacle_info.y, obstacle_info.r)
    plot!(p[1], circle_x, circle_y,
          color=:red, linewidth=2, fill=true, fillalpha=0.3, label="Obstacle")

    # Plot 2: Velocity and heading
    plot!(p[2], ocp.r.Tst, ocp.r.X[:, 6],
          title="Vehicle Heading and Velocity",
          xlabel="Time [s]",
          ylabel="Longitudinal Velocity [m/s]",
          linewidth=2, color=:blue, label="Velocity ux")

    plot!(p[2], ocp.r.Tst, rad2deg.(ocp.r.X[:, 5]),
          ylabel="Heading [deg]",
          linewidth=2, color=:red, label="Heading ψ")

    # Plot 3: Steering control
    plot!(p[3], ocp.r.Tst, rad2deg.(ocp.r.X[:, 7]),
          title="Steering Control",
          xlabel="Time [s]",
          ylabel="Front Wheel Angle [deg]",
          linewidth=2, color=:green, label="δf")

    # Plot 4: Control inputs
    plot!(p[4], ocp.r.Tst, ocp.r.U[:, 1],
          title="Control Inputs",
          xlabel="Time [s]",
          ylabel="Acceleration [m/s²]",
          linewidth=2, color=:red, label="Longitudinal ax")

    plot!(p[4], ocp.r.Tst, rad2deg.(ocp.r.U[:, 2]),
          ylabel="Steering Rate [deg/s]",
          linewidth=2, color=:blue, label="Steering Rate dδf")

    plot!(p, plot_title="Vehicle Obstacle Avoidance Optimal Control")

    if save_plots
        mkpath("figures")
        savefig(p, "figures/obs_avoid_solution.png")
        println("Saved comprehensive plot as 'figures/obs_avoid_solution.png'")
    end

    return p
end

"""
    create_obs_avoid_animation_frames(ocp, obstacle_info; dt=0.1)

Create animation frames showing the vehicle navigating around the obstacle.
"""
function create_obs_avoid_animation_frames(ocp, obstacle_info; dt=0.05)
    println("Creating animation frames with Δt = $dt s ...")

    # Extract solution
    time = ocp.r.Tst
    states = ocp.r.X
    @assert length(time) == size(states, 1) "time and states length mismatch"

    # Uniform time grid for animation
    new_time = collect(time[1]:dt:time[end])

    # Interpolate vehicle states
    nstate = size(states, 2)
    new_states = Array{Float64}(undef, length(new_time), nstate)
    for j in 1:nstate
        itp = Interpolations.linear_interpolation(time, states[:, j])
        new_states[:, j] = itp.(new_time)
    end

    # Animation parameters
    vehicle_length = VEHICLE_LENGTH
    vehicle_width = VEHICLE_WIDTH

    frames = Vector{Plots.Plot{Plots.GRBackend}}()

    for (frame_num, t) in enumerate(new_time)
        # Vehicle position and orientation
        x_veh = new_states[frame_num, 1]
        y_veh = new_states[frame_num, 2]
        psi = new_states[frame_num, 5]  # Vehicle heading
        ux = new_states[frame_num, 6]   # Longitudinal velocity

        # Create frame
        p = plot(xlims=(-15, 80), ylims=(-15, 20),
                 aspect_ratio=:equal,
                 title="Vehicle Obstacle Avoidance - Frame $frame_num/$(length(new_time))",
                 xlabel="X Position [m]", ylabel="Y Position [m]",
                 legend=:topright, grid=true)

        # Vehicle trajectory (past)
        past_idx = findall(time .<= t)
        if length(past_idx) > 1
            plot!(p, states[past_idx, 1], states[past_idx, 2],
                  color=:blue, linewidth=2, alpha=0.6, label="Vehicle Path")
        end

        # Draw circular obstacle
        circle_x, circle_y = circleShape(obstacle_info.x, obstacle_info.y, obstacle_info.r)
        plot!(p, circle_x, circle_y,
              color=:red, linewidth=2, fill=true, fillalpha=0.4, label="")

        # Draw vehicle (as oriented rectangle)
        cos_psi, sin_psi = cos(psi), sin(psi)
        veh_corners_x = [
            x_veh + cos_psi * vehicle_length/2 - sin_psi * vehicle_width/2,
            x_veh + cos_psi * vehicle_length/2 + sin_psi * vehicle_width/2,
            x_veh - cos_psi * vehicle_length/2 + sin_psi * vehicle_width/2,
            x_veh - cos_psi * vehicle_length/2 - sin_psi * vehicle_width/2,
            x_veh + cos_psi * vehicle_length/2 - sin_psi * vehicle_width/2
        ]
        veh_corners_y = [
            y_veh + sin_psi * vehicle_length/2 + cos_psi * vehicle_width/2,
            y_veh + sin_psi * vehicle_length/2 - cos_psi * vehicle_width/2,
            y_veh - sin_psi * vehicle_length/2 - cos_psi * vehicle_width/2,
            y_veh - sin_psi * vehicle_length/2 + cos_psi * vehicle_width/2,
            y_veh + sin_psi * vehicle_length/2 + cos_psi * vehicle_width/2
        ]
        plot!(p, veh_corners_x, veh_corners_y,
              color=:blue, linewidth=2, fill=true, fillalpha=0.3, label="")

        # Add velocity vector
        v_scale = 3.0
        v_x = ux * cos_psi * v_scale
        v_y = ux * sin_psi * v_scale
        plot!(p, [x_veh, x_veh + v_x], [y_veh, y_veh + v_y],
              arrow=true, color=:green, linewidth=2, label="")

        # Add start and target markers
        scatter!(p, [states[1,1]], [states[1,2]],
                color=:green, markersize=8, label="")
        scatter!(p, [states[end,1]], [states[end,2]],
                color=:magenta, markersize=8, label="")

        annotate!(p, x_veh + 5, y_veh + 8,
                 text("t = $(round(t, digits=2)) s\\nSpeed = $(round(ux, digits=1)) m / s", 10))

        push!(frames, p)
    end

    println("Created $(length(new_time)) frames.")
    return frames
end

"""
    save_obs_avoid_frames_as_gif(frames; fps=10, path="figures/obs_avoid_animation.gif")

Save animation frames as GIF file.
"""
function save_obs_avoid_frames_as_gif(frames; fps=10, path="figures/obs_avoid_animation.gif")
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

Main function to run the complete vehicle obstacle avoidance example.
"""
function main()
    println("Vehicle Obstacle Avoidance Optimal Control Example")
    println("=" ^ 50)

    # Set up the problem
    ocp = setup_obs_avoid_problem()
    formulation = configure_obs_avoid_formulation(ocp)
    obstacle_info = add_obstacle_constraints!(ocp, formulation)

    # Solve the problem
    success = solve_obs_avoid(ocp)

    if success
        # Visualize results
        visualize_obs_avoid_results(ocp, obstacle_info, save_plots=true)

        # Create animation frames
        frames = create_obs_avoid_animation_frames(ocp, obstacle_info, dt=0.05)

        # Display key frames
        display(frames[1])              # Initial frame
        display(frames[length(frames)÷2]) # Middle frame
        display(frames[end])            # Final frame

        # Save as GIF
        save_obs_avoid_frames_as_gif(frames, fps=20, path="figures/obs_avoid_animation.gif")

        println("\nVehicle Obstacle Avoidance example completed successfully!")
        println("Summary:")
        println("   - Obstacle successfully avoided")
        println("   - Target reached with optimal path")
        println("   - Time horizon: $(formulation.tf) seconds")
        println("   - Discretization points: $(formulation.Np)")
        println("   - Integration scheme: $(formulation.IntegrationScheme[1])")

        return ocp, frames
    else
        println("Vehicle Obstacle Avoidance example failed to find solution")
        return nothing, nothing
    end
end

# Run the example
main()

nothing