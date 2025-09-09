# Racing Vehicle Optimal Control Example

This example demonstrates high-performance racing vehicle trajectory optimization using JuliaOptimalControl. The problem involves computing the fastest possible path around the Thunderhill West raceway while respecting track boundaries, vehicle dynamics limits, and safety constraints.

## Problem Description

The racing optimization problem seeks to find the minimum-time trajectory for a high-performance vehicle navigating a race track. This involves sophisticated vehicle dynamics, tire models, track constraints, and advanced optimization techniques used in professional motorsports.

### System Model

The vehicle uses an advanced 8-state dynamic model capturing essential racing vehicle behavior:

**States (8):**
- `x`: Global X position [m]
- `y`: Global Y position [m]
- `v`: Lateral velocity [m/s]
- `r`: Yaw rate [rad/s]
- `ψ`: Heading angle [rad]
- `ux`: Longitudinal velocity [m/s]
- `sa`: Steering angle [rad]
- `ax`: Longitudinal acceleration [m/s²]

**Controls (2):**
- `sr`: Steering rate [rad/s]
- `jx`: Longitudinal jerk [m/s³]

### Vehicle Dynamics

The model incorporates realistic racing vehicle dynamics:

```julia
# Position kinematics
ẋ = ux*cos(ψ) - v*sin(ψ)
ẏ = ux*sin(ψ) + v*cos(ψ)

# Lateral dynamics with tire forces
v̇ = (Fyf + Fyr)/M - ux*r

# Yaw dynamics
ṙ = (lf*Fyf - lr*Fyr)/Izz

# Longitudinal dynamics
u̇x = ax

# Control derivatives
ṡa = sr
ȧx = jx
```

Where tire forces `Fyf, Fyr` are computed using advanced tire models accounting for:
- Normal load transfer during acceleration/braking
- Lateral load transfer during cornering
- Combined slip conditions
- Friction circle constraints

## Track Modeling

The Thunderhill West raceway is modeled using:

### Track Geometry
- **Center line**: Splined reference trajectory
- **Boundaries**: Left and right track limits
- **Curvature**: Variable radius corners and straights
- **Elevation**: Track banking and elevation changes

### Safety Constraints

**Drivability constraints** using smoothed barrier functions:
```julia
# Track boundary constraint (log-sum-exp smoothing)
h(x,y) = (1/ρ)*log(∑exp(ρ*(-d_track(x,y)))) ≥ -safety_margin
```

**Control barrier functions** for collision avoidance:
```julia
# CBF constraint ensuring safe trajectory evolution
ḣ ≥ -α*h + compensation_term
```

Where the compensation term accounts for system dynamics uncertainty.

## Cost Function

The racing optimization minimizes lap time while penalizing aggressive control:

```julia
cost = w₁*obstácle_cost + w₂*progress_cost + w₃*sr² + w₄*sa² + w₅*ax² + w₆*jx² + w₇*curvature² + w₈*lateral_velocity²
```

This balances:
- **Progress maximization**: Advancing around the track
- **Safety**: Maintaining distance from track boundaries
- **Control smoothness**: Avoiding jerky inputs
- **Vehicle stability**: Limiting lateral velocity and curvature

## Advanced Features

### Model Predictive Control (MPC)
The example implements real-time MPC with:
- **Prediction horizon**: 4 seconds
- **Discretization**: 25 time points
- **Update rate**: 10 Hz for real-time capability
- **Warm starting**: Using previous solution for faster convergence

### Robust Optimization
- **Uncertainty compensation**: Accounts for model uncertainty
- **Safety margins**: Conservative constraints for racing conditions
- **Adaptive parameters**: Adjusts based on vehicle state

## Usage

Run the racing optimization example:

```julia
# Navigate to the racing example directory
cd("examples/racing")

# Ensure track data is available
# Note: Requires "processed_ThunderHill_West.mat" file

# Run the main script
include("main.jl")
```

The example will:
1. Load Thunderhill West track data
2. Set up advanced vehicle dynamics model
3. Configure MPC optimization with safety constraints
4. Run real-time simulation with trajectory updates
5. Generate racing line visualization

## Results

The optimal solution demonstrates professional racing techniques:

![Racing Solution](figures/racing_solution.png)

The racing line shows:
- **Corner entry**: Late braking with trail-braking technique
- **Apex targeting**: Geometric racing line through corners
- **Corner exit**: Early throttle application for acceleration
- **Straight sections**: Maximum speed with stability

### Code Structure

```julia
# Main MPC execution script
include("main.jl")

# Advanced vehicle dynamics
include("VehicleModel.jl")

# Track processing utilities
include("racing_utils.jl")

# Key functions:
DefineOCP()                    # Set up racing optimization problem
LC500Model()                   # High-fidelity vehicle model
FindCenterLine()               # Track geometry processing
GetLengthParams()              # Progress parameterization
OptSolve!()                    # Real-time optimization
```

### Performance Characteristics

- **Lap time optimization**: Finds theoretically fastest racing line
- **Real-time capability**: 10 Hz MPC updates for closed-loop control
- **Safety assurance**: CBF constraints prevent track boundary violations
- **Vehicle limits**: Respects tire friction and stability constraints

## Physical Interpretation

This represents the state-of-the-art in autonomous racing and high-performance vehicle control:

### Applications
- **Autonomous racing**: Formula Student Driverless, Roborace
- **Performance optimization**: Track day and racing driver aids
- **Vehicle testing**: Automated testing of vehicle dynamics
- **Motorsports**: Racing line analysis and driver training

### Key Racing Insights

1. **Geometric vs. dynamic optimization**: Trading curvature for time
2. **Tire management**: Optimal friction circle utilization
3. **Energy management**: Balancing speed with control stability
4. **Safety margins**: Racing at the limit with mathematical guarantees

## Numerical Properties

- **Time horizon**: 4 seconds (receding horizon)
- **Discretization**: 25 time points
- **Integration**: Backward Euler for stability
- **Solver**: Ipopt with custom tolerances for real-time performance
- **Update rate**: 10 Hz MPC implementation
- **Convergence**: Typically 5-15 iterations per solve

## Extensions

Enhance the racing example:

1. **Multi-vehicle racing**: Overtaking and defensive strategies
2. **Weather conditions**: Wet track dynamics and reduced friction
3. **Tire degradation**: Long-stint optimization with changing grip
4. **Fuel strategy**: Endurance racing with fuel consumption
5. **3D tracks**: Elevation changes and banking effects
6. **Real vehicle integration**: Hardware-in-the-loop testing

This example represents the cutting edge of autonomous racing technology, combining advanced optimal control theory with practical real-time implementation for high-performance applications.