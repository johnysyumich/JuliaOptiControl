# Moving Obstacle Avoidance Example

This example demonstrates optimal trajectory planning for a vehicle avoiding a moving obstacle using JuliaOptimalControl. The problem involves computing a collision-free path from a starting position to a target destination while avoiding a dynamic obstacle that moves with known trajectory.

## Problem Description

The vehicle must navigate from an initial position to a target location while avoiding collision with a moving obstacle. The obstacle follows a predetermined path with known position and velocity at each time instant. The optimization finds the optimal vehicle trajectory that minimizes travel time and control effort while maintaining safe separation.

### System Model

The vehicle is modeled using a 3DOF bicycle model that captures the essential vehicle dynamics:

**States (8):**
- `x`: X position [m]
- `y`: Y position [m]
- `v`: Lateral velocity [m/s]
- `r`: Yaw rate [rad/s]
- `ψ`: Heading angle [rad]
- `sa`: Steering angle [rad]
- `ux`: Longitudinal velocity [m/s]
- `ax`: Longitudinal acceleration [m/s²]

**Controls (2):**
- `sr`: Steering rate [rad/s]
- `jx`: Longitudinal jerk [m/s³]

### Vehicle Dynamics

The dynamics are based on a bicycle model with longitudinal and lateral coupling:

```julia
# Position dynamics
ẋ = ux*cos(ψ) - v*sin(ψ)
ẏ = ux*sin(ψ) + v*cos(ψ)

# Velocity dynamics
v̇ = (Fyf + Fyr)/m - ux*r
u̇x = ax

# Yaw dynamics
ṙ = (lf*Fyf - lr*Fyr)/Iz
ψ̇ = r

# Steering dynamics
ṡa = sr
ȧx = jx
```

Where tire forces are computed using the Pacejka tire model for realistic vehicle behavior.

## Obstacle Avoidance

The moving obstacle is modeled as an ellipse with time-varying position:

```julia
# Obstacle constraint (ellipsoidal safety region)
1 ≤ ((x - x_obs(t))²/a² + (y - y_obs(t))²/b²)
```

Where:
- `x_obs(t), y_obs(t)`: Obstacle position at time t
- `a, b`: Safety margin dimensions around obstacle

## Cost Function

The objective balances multiple competing goals:

```julia
cost = ∫(w₁*(x - x_ref)² + w₂*v² + w₃*sr² + w₄*jx²)dt + tf + w₅*(y_final - y_target)² + w₆*(x_final - x_target)²
```

This encourages:
- Following a reference lateral position (lane keeping)
- Minimizing lateral velocity for stability
- Smooth steering inputs
- Smooth longitudinal control
- Fast arrival time (minimize tf)
- Accurate final positioning

## Usage

Run the moving obstacle avoidance example:

```julia
# Navigate to the movingobs example directory
cd("examples/movingobs")

# Run the main script
include("main.jl")
```

The example will:
1. Set up the vehicle model and obstacle trajectory
2. Configure ellipsoidal constraint for collision avoidance
3. Solve the trajectory optimization problem
4. Generate visualization showing vehicle and obstacle paths

## Results

The optimal solution demonstrates intelligent collision avoidance behavior:

![Moving Obstacle Solution](figures/movingobs_solution.png)

The trajectory shows:
- **Initial phase**: Vehicle maintains reference trajectory
- **Avoidance phase**: Vehicle deviates laterally to avoid collision
- **Recovery phase**: Vehicle returns to target trajectory after obstacle passes

### Code Structure

```julia
# Main execution script
include("main.jl")

# Vehicle dynamics model
include("ThreeDOF_Bicycle.jl")

# System parameters
include("parameter.jl")

# Key components:
defineOCP()                    # Set up optimization problem
ConfigurePredefined()          # Configure discretization
@constraint()                  # Add obstacle avoidance constraint
@objective()                   # Define cost function
OptSolve!()                   # Solve optimization
```

### Solution Characteristics

- **Solution time**: ~2-5 seconds depending on complexity
- **Collision avoidance**: Maintains safe ellipsoidal separation
- **Control smoothness**: RK2 integration ensures realistic vehicle motion
- **Adaptability**: Automatically adjusts path timing based on obstacle motion

## Physical Interpretation

This problem represents fundamental challenges in autonomous vehicle navigation:

### Applications
- **Autonomous vehicles**: Highway lane change and obstacle avoidance
- **Robot navigation**: Mobile robot path planning in dynamic environments
- **UAV control**: Drone collision avoidance in shared airspace
- **Marine vehicles**: Ship navigation with moving traffic

### Key Engineering Insights

1. **Predictive avoidance**: Using known obstacle trajectory for optimal planning
2. **Safety margins**: Ellipsoidal constraints provide robust collision avoidance
3. **Multi-objective optimization**: Balancing safety, comfort, and efficiency
4. **Real-time capability**: Fast enough for model predictive control applications

## Numerical Properties

- **Time horizon**: Variable (optimized)
- **Discretization**: 25 time points
- **Integration**: 2nd-order Runge-Kutta (RK2)
- **Constraints**: Ellipsoidal obstacle avoidance
- **Solver**: Ipopt nonlinear programming

## Extensions

Try modifying the example:

1. **Multiple obstacles**: Add more moving obstacles with different trajectories
2. **Uncertainty**: Add robust constraints for uncertain obstacle motion
3. **Traffic rules**: Include lane boundaries and traffic regulations
4. **Vehicle limits**: Add more detailed tire friction and stability constraints
5. **Comfort constraints**: Limit lateral acceleration for passenger comfort
6. **Emergency scenarios**: High-speed obstacle avoidance maneuvers

This example showcases the power of trajectory optimization for safety-critical autonomous vehicle applications where collision avoidance is paramount.