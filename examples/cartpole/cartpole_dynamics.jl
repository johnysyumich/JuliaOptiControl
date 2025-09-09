"""
CartPole Dynamics

This file contains the dynamics model for the classic cart-pole system.
The cart-pole is a benchmark problem in control theory where we must balance
a pole on a cart by applying horizontal forces to the cart.

States:
- x:  Cart position [m]
- ẋ:  Cart velocity [m/s]
- θ:  Pole angle from vertical [rad] (0 = upright)
- θ̇:  Pole angular velocity [rad/s]

Controls:
- F:  Horizontal force applied to cart [N]

Physical Parameters:
- M:  Cart mass [kg]
- m:  Pole mass [kg]
- L:  Pole half-length [m]
- g:  Gravitational acceleration [m/s²]
"""

# Physical parameters for the cart-pole system
const CART_MASS = 1.0      # Cart mass [kg]
const POLE_MASS = 0.1      # Pole mass [kg]
const POLE_LENGTH = 0.5    # Pole half-length [m]
const GRAVITY = 9.81       # Gravitational acceleration [m/s²]

"""
    cartpole_dynamics(states, controls, params)

Compute the cart-pole system dynamics.

The cart-pole dynamics are derived from Lagrangian mechanics:
- Cart equation: (M + m)ẍ + mL(θ̈cos(θ) - θ̇²sin(θ)) = F
- Pole equation: mLẍcos(θ) + mL²θ̈ = mgLsin(θ)

# Arguments
- states: [x, ẋ, θ, θ̇] - current state vector
- controls: [F] - control input (horizontal force)
- params: unused parameter vector

# Returns
- [ẋ, ẍ, θ̇, θ̈] - state derivatives
"""
function cartpole_dynamics(states, controls, params)
    # Extract states
    x, x_dot, theta, theta_dot = states

    # Extract control input
    F = controls[1]

    # System parameters
    M = CART_MASS
    m = POLE_MASS
    L = POLE_LENGTH
    g = GRAVITY

    # Trigonometric functions
    sin_theta = sin(theta)
    cos_theta = cos(theta)

    # Common terms for efficiency
    denominator = M + m * sin_theta^2

    # Cart acceleration (from solving the coupled equations)
    x_ddot = (F + m * L * theta_dot^2 * sin_theta - m * g * sin_theta * cos_theta) / denominator

    # Pole angular acceleration
    theta_ddot = (-F * cos_theta - m * L * theta_dot^2 * sin_theta * cos_theta +
                  (M + m) * g * sin_theta) / (L * denominator)

    # Return state derivatives: [ẋ, ẍ, θ̇, θ̈]
    return [x_dot, x_ddot, theta_dot, theta_ddot]
end

"""
    cartpole_cost(states, controls, params)

Cost function for cart-pole balancing.

Penalizes:
- Pole angle deviation from upright (θ = 0)
- Cart position deviation from center (x = 0)
- High velocities
- Large control efforts

# Arguments
- states: [x, ẋ, θ, θ̇] - current state vector
- controls: [F] - control input
- params: unused parameter vector

# Returns
- cost: scalar cost value
"""
function cartpole_cost(states, controls, params)
    # Extract states and controls
    x, x_dot, theta, theta_dot = states
    F = controls[1]

    # Cost weights
    Q_x = 1.0      # Cart position penalty
    Q_x_dot = 0.1  # Cart velocity penalty
    Q_theta = 100.0  # Pole angle penalty (high priority!)
    Q_theta_dot = 1.0  # Pole angular velocity penalty
    R = 0.01       # Control effort penalty

    # Quadratic cost function
    cost = Q_x * x^2 +
           Q_x_dot * x_dot^2 +
           Q_theta * theta^2 +
           Q_theta_dot * theta_dot^2 +
           R * F^2

    return cost
end