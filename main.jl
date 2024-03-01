include("src/setup.jl")
include("parameter.jl")
include("ThreeDOF_Bicycle.jl")
XL = [-2.6, 0, NaN, NaN, -6.283185307179586, -0.62, 10, -2.6]
XU = [1.801, 300, NaN, NaN, 6.283185307179586, 0.62, 20,  2.6]
CL=[-2.1, -6.2]
CU=[1.9, 6.2]
ocp = defineOCP(numStates=8,numControls=2,X0=[1.8, 0., 0.0, 0.0, pi/2, 0.0, 15., 0.],XF=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN], XL=XL, XU=XU, CL=CL,CU=CU);
defineStates!(ocp, [:x,:y,:v,:r,:psi,:sa,:ux,:ax])
defineControls!(ocp, [:sr, :jx])
OCPForm = ConfigurePredefined(ocp; (:Np=>11), (:tfDV => false), (:tf=>3), (:dynamics => ThreeDOFBicycle_expr))

CheckOCPFormulation(ocp, OCPForm)
## TODO Configure the model setting
DefineSolver!(OCPForm, ocp.s.solver.name, ocp.s.solver.settings)
## Currently only collocation.
## TODO Write the initial condition constraints
ocp.s.states.pts = ocp.s.control.pts = OCPForm.Np
ocp.p.x = @variable(OCPForm.mdl, ocp.b.XL[i] <= sts[j in 1:ocp.s.states.pts, i in 1:ocp.s.states.num] <= ocp.b.XU[i])
ocp.p.u = @variable(OCPForm.mdl, ocp.b.CL[i] <= ctr[j in 1:ocp.s.control.pts, i in 1:ocp.s.control.num] <= ocp.b.CU[i])

# Dynamical Constraints
if OCPForm.IntegrationScheme == :bkwEuler
    δx = Matrix{Any}(undef, ocp.s.states.pts, ocp.s.states.num)
    for j in 1:ocp.s.states.pts
        δx[j, :] = @expression(OCPForm.mdl, OCPForm.dx[j](ocp.p.x[j, :], ocp.p.u[j, :]))
    end
    dyncon = Matrix{Any}(undef, ocp.s.states.pts - 1, ocp.s.states.num)
    for k in 1:ocp.s.states.num
        for i in 1:ocp.s.states.pts - 1
            dyncon[i, k] = @constraint(OCPForm.mdl, ocp.p.x[i + 1, k]- ocp.p.x[i, k] == δx[i + 1, k] * OCPForm.TInt[i])
        end
    end
end

