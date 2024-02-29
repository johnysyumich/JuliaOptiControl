include("src/setup.jl")
XL = [-2.6, 0, NaN, NaN, -6.283185307179586, -0.62, 10, -2.6]
XU = [1.801, 300, NaN, NaN, 6.283185307179586, 0.62, 20,  2.6]
CL=[-2.1, -6.2]
CU=[1.9, 6.2]
ocp = defineOCP(numStates=8,numControls=2,X0=[1.8, 0., 0.0, 0.0, pi/2, 0.0, 15., 0.],XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);
defineStates!(ocp, [:x,:y,:v,:r,:psi,:sa,:ux,:ax])
defineControls!(ocp, [:sr, :jx])
ocpf = ConfigurePredefined(ocp; (:Np=>11), (:tfDV => true), (:IntegrationScheme=>:trapezoidal))
