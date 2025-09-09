using NLOCPSolver
using Plots
using JuMP
include("bicycleModel.jl")
include("parameters.jl")

XL = [-40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12]
XU = [300, 20, 3, pi/5, pi/2, 15.0, pi/12]
CL = [-2.6, -0.1]
CU = [2.6, 0.1]
X0 = [-10.0, 0, 0, 0, 0, 10.0, 0]
XF = [NaN, NaN, NaN, NaN, NaN, NaN, NaN]
ocp = defineOCP(numStates=7, numControls=2, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU);
defineStates!(ocp, [:x, :y, :v, :r, :ψ, :ux, :δf])
defineControls!(ocp, [:ax, :dδf])
OCPForm = ConfigurePredefined(ocp; (:Np => 81), (:tfDV => false), (:tf => 8), (:IntegrationScheme => :RK1), (:dx => bicycleModel_expr), (:expr => bicycle_cost))
# user_options = ()


OCPdef!(ocp, OCPForm)
block_center_x_info = 30;
block_center_y_info = 2;
block_radius_info = 6;
ocp.f.mdl[:block_center_x] = @variable(ocp.f.mdl, block_center_x in Parameter(block_center_x_info));
ocp.f.mdl[:block_center_y] = @variable(ocp.f.mdl, block_center_y in Parameter(block_center_y_info));
ocp.f.mdl[:block_radius] = @variable(ocp.f.mdl, block_radius in Parameter(block_radius_info));
xpos = ocp.p.x[:, 1];
y = ocp.p.x[:, 2]; dδf = ocp.p.u[:, 2]; ax = ocp.p.u[:, 1]; ux = ocp.p.x[:, 6]; δf = ocp.p.x[:, 7];
# obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts], 36 <= ((xpos[i] - 30).^2 + (y[i] - 2).^2));
set_parameter_value.(ocp.f.mdl[:block_center_x], block_center_x_info);
set_parameter_value.(ocp.f.mdl[:block_center_y], block_center_y_info);
set_parameter_value.(ocp.f.mdl[:block_radius], block_radius_info);
obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts],  block_radius^2 <= ((xpos[i] - block_center_x).^2 + (y[i] - block_center_y).^2));
# obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts-1], (xpos[i+1] - block_center_x).^2 + (y[i+1] - block_center_y).^2 - (xpos[i] - block_center_x).^2 - (y[i] - block_center_y).^2 >= -1 * ((xpos[i] - block_center_x).^2 + (y[i] - block_center_y).^2 - block_radius^2))
obj = @expression(ocp.f.mdl, sum((0.05 * (y[j])^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 + 0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
@objective(ocp.f.mdl, Min, obj)
@time OptSolve!(ocp)
function circleShape(h,k,r)
    theta = LinRange(0,2*pi,500)
    h .+ r*sin.(theta), k .+ r*cos.(theta)
end
# Create comprehensive visualization
p = plot(layout=(2,2), size=(1000, 800))

# Plot 1: Vehicle trajectory with obstacle
plot!(p[1], ocp.r.X[:, 1], ocp.r.X[:, 2],
      title="Vehicle Path Planning with Obstacle",
      xlabel="X Position [m]",
      ylabel="Y Position [m]",
      linewidth=3, color=:blue, label="Vehicle Path")
plot!(p[1], circleShape(30,2,6), seriestype = [:shape,], lw=2, c=:red,
      linecolor=:black, fillalpha=0.3, aspect_ratio=1, label="Obstacle")
scatter!(p[1], [ocp.r.X[1,1]], [ocp.r.X[1,2]], color=:green, markersize=8, label="Start")
scatter!(p[1], [ocp.r.X[end,1]], [ocp.r.X[end,2]], color=:red, markersize=8, label="End")

# Plot 2: Longitudinal velocity
plot!(p[2], ocp.r.Tst, ocp.r.X[:, 6],
      title="Longitudinal Velocity",
      xlabel="Time [s]",
      ylabel="Velocity [m/s]",
      linewidth=2, color=:green, legend=false)

# Plot 3: Steering angle
plot!(p[3], ocp.r.Tst, rad2deg.(ocp.r.X[:, 7]),
      title="Steering Angle",
      xlabel="Time [s]",
      ylabel="Steering Angle [deg]",
      linewidth=2, color=:orange, legend=false)

# Plot 4: Control inputs
plot!(p[4], ocp.r.Tst, ocp.r.U[:, 1],
      title="Longitudinal Acceleration",
      xlabel="Time [s]",
      ylabel="Acceleration [m/s²]",
      linewidth=2, color=:red, label="Longitudinal")
plot!(p[4], ocp.r.Tst, rad2deg.(ocp.r.U[:, 2]),
      linewidth=2, color=:blue, label="Steering Rate [deg/s]")

plot!(p, plot_title="Vehicle Obstacle Avoidance Optimal Control")
display(p)

# Save the plot
savefig(p, "figures/vehicleOpt_solution.png")
println("Saved vehicle optimization plot as 'figures/vehicleOpt_solution.png'")
