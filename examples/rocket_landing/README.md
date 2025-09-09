# Rocket Landing Optimal Control Example

This example demonstrates optimal control of a rocket performing a vertical landing maneuver, inspired by SpaceX Falcon 9 landings. The goal is to land the rocket softly on the ground while minimizing fuel consumption.

## Problem Description

### System Dynamics
The rocket system consists of:
- **Rocket Body**: Variable mass due to fuel consumption
- **Propulsion**: Single engine providing vertical thrust
- **Environment**: Constant gravitational field

### States
- `h`: Altitude above ground [m]
- `v`: Vertical velocity [m/s] (positive = upward, negative = downward)
- `m`: Total mass [kg] (decreases as fuel burns)

### Control
- `T`: Thrust force [N], bounded to [0, 25000] N

### Objective
Land the rocket safely while minimizing fuel consumption:
- Start from high altitude with downward velocity
- End at ground level with soft landing velocity
- Minimize thrust usage (proportional to fuel consumption)

## Mathematical Model

The rocket dynamics in vertical flight are:

**Altitude equation:**
```
ḣ = v
```

**Velocity equation (Newton's second law):**
```
mv̇ = T - mg
v̇ = T/m - g
```

**Mass equation (rocket equation):**
```
ṁ = -T/ve
```

Where:
- `g = 9.81 m/s²`: gravitational acceleration
- `ve = 2500 m/s`: effective exhaust velocity (typical for LOX/RP-1)

## Cost Function

The optimization minimizes:
- **Primary objective**: Fuel consumption (∝ thrust integral)
- **Terminal cost**: Landing accuracy
  - Altitude deviation from ground: `1000 × h_final²`
  - Velocity deviation from soft landing: `100 × (v_final + 1)²`

## Physical Parameters

- **Dry mass**: 1000 kg (empty rocket structure)
- **Initial fuel**: 500 kg
- **Maximum thrust**: 25,000 N (thrust-to-weight ratio ≈ 1.7)
- **Exhaust velocity**: 2500 m/s (specific impulse ≈ 255 s)

## Constraints

### State Constraints
- `h ≥ 0`: No ground penetration
- `m ≥ dry_mass`: Can't burn rocket structure
- Velocity and altitude bounds for numerical stability

### Path Constraints
- `T ≥ 0`: No negative thrust (engine can't suck exhaust)
- `T ≤ T_max`: Thrust limit
- Monotonic mass decrease: `ṁ ≤ 0`

### Terminal Constraints
- `h_final = 0`: Land on ground
- `-5 ≤ v_final ≤ 0`: Soft landing velocity range

## Files

- `rocket_dynamics.jl`: Contains the dynamics model and cost functions
- `main.jl`: Main script to set up and solve the optimal control problem
- `README.md`: This documentation file

## Usage

Run the rocket landing optimization example:

```julia
# Navigate to the rocket landing example directory
cd("examples/rocket_landing")

# Run the complete example
include("main.jl")

# Execute the optimization
ocp, frames = main()
```

## Results

The optimal solution demonstrates the classic rocket landing strategy:

![Rocket Landing Solution](figures/rocket_landing_solution.png)

The trajectory shows:
1. **Initial Phase**: High thrust to decelerate rapidly
2. **Coasting Phase**: Minimal thrust to conserve fuel
3. **Final Approach**: Precise thrust control for soft landing

### Code Structure

```julia
# Main execution script
include("main.jl")

# System dynamics and cost functions
include("rocket_dynamics.jl")

# Key functions:
setup_rocket_problem()       # Define states, controls, bounds
configure_solver()           # Set up discretization and integration
solve_rocket()              # Solve the optimization problem
visualize_results()         # Create comprehensive plots
create_landing_animation()  # Generate animation frames
main()                      # Execute complete example
```

### Typical Solution Characteristics

- **Landing time**: ~20-25 seconds
- **Fuel efficiency**: 60-80% fuel consumption
- **Thrust profile**: High initial thrust, then throttled approach
- **Velocity profile**: Rapid deceleration followed by controlled descent

### Visualizations

The example generates comprehensive plots saved to `figures/`:
- **Altitude profile**: Height vs. time trajectory
- **Velocity profile**: Speed vs. time (shows deceleration)
- **Thrust profile**: Engine throttle vs. time
- **Thrust-to-weight ratio**: T/W over time
- **Phase portrait**: Altitude vs. velocity trajectory
- **Animation frames**: Landing sequence visualization

## Physical Interpretation

This problem represents fundamental challenges in aerospace engineering:

### Applications
- **Rocket landing systems** (SpaceX Falcon 9, Blue Origin New Shepard)
- **Lunar/planetary landers** (Apollo LM, Chang'e series)
- **Spacecraft orbital maneuvers**
- **Launch vehicle trajectory optimization**

### Key Engineering Insights

1. **Fuel-optimality vs. Time-optimality**: Trade-off between fast landing and fuel conservation
2. **Thrust throttling**: Deep throttling capability crucial for precise landing
3. **Gravity losses**: Hovering wastes fuel; optimal trajectory minimizes hover time
4. **Terminal guidance**: Final approach requires high precision control

## Numerical Properties

- **Time horizon**: 25 seconds (optimized)
- **Discretization**: 41 time points
- **Integration**: 4th-order Runge-Kutta (RK4)
- **Solver**: Ipopt with automatic differentiation
- **Typical solve time**: < 2 seconds
- **Convergence**: Usually 20-60 iterations

## Extensions and Variations

Try modifying the example:

1. **3D Landing**: Add horizontal position and lateral thrust
2. **Atmospheric effects**: Include drag forces
3. **Landing pad constraints**: Target specific landing zone
4. **Engine dynamics**: Add thrust buildup/shutdown delays
5. **Uncertainty**: Robust control with wind disturbances
6. **Multi-stage**: Optimize entire launch-to-landing trajectory
7. **Landing legs**: Model gear deployment and ground contact
8. **Fuel slosh**: Add fuel dynamics effects

## Advanced Topics

### Optimal Control Theory
- This is a **Bolza problem** (integral + terminal cost)
- Demonstrates **singular arcs** in thrust control
- Shows **bang-bang-singular** control structure

### Rocket Equation
- Classic **Tsiolkovsky rocket equation** in action
- Exponential relationship between velocity change and fuel consumption
- Importance of specific impulse in mission design

### Landing Guidance
- Real rockets use **convex optimization** for real-time guidance
- This example shows the underlying optimal control structure
- Practical systems add robustness and computational constraints

## References

1. Açıkmeşe, B., & Ploen, S. R. (2007). *Convex programming approach to powered descent guidance for mars landing*. Journal of guidance, control, and dynamics, 30(5), 1353-1366.

2. Blackmore, L., Açıkmeşe, B., & Scharf, D. P. (2010). *Minimum-landing-error powered-descent guidance for Mars landing using convex optimization*. Journal of guidance, control, and dynamics, 33(4), 1161-1171.

3. Mease, K. D., Chen, D. T., Teufel, P., & Schönenberger, H. (2002). *Reduced-order optimal control of the powered descent phase of Mars landing*. Journal of guidance, control, and dynamics, 25(2), 389-397.

4. Sagliano, M. (2019). *Pseudospectral convex optimization for powered descent and landing*. Journal of Guidance, Control, and Dynamics, 42(7), 1500-1514.