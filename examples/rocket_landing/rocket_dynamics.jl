"""
Rocket Landing Dynamics

This file contains the dynamics model for a vertical rocket landing problem.
The rocket must navigate from an initial altitude and velocity to a soft landing
while minimizing fuel consumption and avoiding ground collision.

States:
- h:  Altitude [m] (height above ground)
- v:  Vertical velocity [m/s] (positive = upward)
- m:  Mass [kg] (decreases as fuel burns)

Controls:
- T:  Thrust force [N] (positive = upward)

Physical Parameters:
- g:  Gravitational acceleration [m/s²]
- ve: Effective exhaust velocity [m/s]
"""

# Physical parameters for the rocket system
const GRAVITY = 9.81           # Gravitational acceleration [m/s²]
const EXHAUST_VELOCITY = 2500.0  # Effective exhaust velocity [m/s] (typical for LOX/RP-1)
const DRY_MASS = 1000.0       # Dry mass (empty rocket) [kg]
const FUEL_MASS = 500.0       # Initial fuel mass [kg]
const MAX_THRUST = 25000.0    # Maximum thrust [N]

"""
    rocket_dynamics(states, controls, params)

Compute the rocket landing system dynamics (simplified fixed-mass version).

The rocket dynamics in vertical flight are:
- Altitude rate: ḣ = v
- Velocity rate: v̇ = T/m - g (with fixed mass)

# Arguments
- states: [h, v] - current state vector (altitude, velocity)
- controls: [T] - control input (thrust force)
- params: unused parameter vector

# Returns
- [ḣ, v̇] - state derivatives
"""
function rocket_dynamics(states, controls, params)
    # Extract states
    h, v = states

    # Extract control input
    T = controls[1]

    # System parameters
    g = GRAVITY
    m = DRY_MASS + FUEL_MASS/2  # Use average mass for simplicity

    # State derivatives
    h_dot = v                    # Altitude rate = velocity
    v_dot = T/m - g             # Velocity rate = thrust/mass - gravity

    # Return state derivatives: [ḣ, v̇]
    return [h_dot, v_dot]
end

"""
    rocket_cost(states, controls, params)

Cost function for rocket landing optimization (simplified version).

Minimizes fuel consumption (proportional to thrust).

# Arguments
- states: [h, v] - current state vector
- controls: [T] - control input
- params: unused parameter vector

# Returns
- cost: scalar cost value
"""
function rocket_cost(states, controls, params)
    # Extract states and controls
    h, v = states
    T = controls[1]

    # Cost weights
    fuel_penalty = 0.001    # Fuel consumption penalty (proportional to thrust)
    control_penalty = 1e-8  # Control smoothness penalty

    # Primary objective: minimize fuel consumption (proportional to thrust)
    fuel_cost = fuel_penalty * T

    # Secondary objective: smooth control
    control_effort = control_penalty * T^2

    # Total cost
    cost = fuel_cost + control_effort

    return cost
end

"""
    rocket_final_cost(states, controls, params)

Terminal cost function for rocket landing.

Heavily penalizes final state deviations from desired landing conditions:
- Final altitude should be 0 (on ground)
- Final velocity should be small and negative (soft touchdown)
- Remaining fuel is not penalized (can be saved)

# Arguments
- states: [h, v, m] - final state vector
- controls: [T] - final control input
- params: unused parameter vector

# Returns
- cost: scalar terminal cost value
"""
function rocket_final_cost(states, controls, params)
    # Extract final states
    h, v, m = states

    # Terminal cost weights
    altitude_penalty = 1000.0    # Heavy penalty for not reaching ground
    velocity_penalty = 500.0     # Heavy penalty for hard landing

    # Desired final conditions
    target_altitude = 0.0        # Land on ground
    target_velocity = -1.0       # Soft landing velocity (slightly downward)

    # Terminal cost
    altitude_cost = altitude_penalty * (h - target_altitude)^2
    velocity_cost = velocity_penalty * (v - target_velocity)^2

    # Total terminal cost
    terminal_cost = altitude_cost + velocity_cost

    return terminal_cost
end