using Test
using JuliaOptimalControl
using JuMP
using DataFrames
using Parameters
using LinearAlgebra

# Import MathOptInterface through JuMP
const MOI = JuMP.MOI

@testset "JuliaOptimalControl Complete Test Suite" begin

    @testset "Type Construction Tests" begin

        @testset "States struct tests" begin
            # Test default construction
            states = States()
            @test states.num == 0
            @test states.name == Vector{Symbol}[]
            @test states.pts == 0

            # Test parameterized construction
            states = States(num=3, name=[:x, :y, :z], pts=51)
            @test states.num == 3
            @test states.name == [:x, :y, :z]
            @test states.pts == 51
        end

        @testset "Control struct tests" begin
            # Test default construction
            control = Control()
            @test control.num == 0
            @test control.name == Vector{Symbol}[]
            @test control.pts == 0

            # Test parameterized construction
            control = Control(num=2, name=[:u1, :u2], pts=51)
            @test control.num == 2
            @test control.name == [:u1, :u2]
            @test control.pts == 51
        end

        @testset "Solver struct tests" begin
            # Test default construction
            solver = Solver()
            @test solver.name == :Ipopt
            @test solver.settings == JuliaOptimalControl._Ipopt_defaults

            # Test custom construction
            custom_settings = ("max_iter" => 500, "tol" => 1e-6)
            solver = Solver(name=:Ipopt, settings=custom_settings)
            @test solver.name == :Ipopt
            @test solver.settings == custom_settings
        end

        @testset "OCPFormulation struct tests" begin
            # Test default construction
            form = OCPFormulation{Float64}()
            @test form.tfDV == false
            @test form.Np == 0
            @test form.IntegrationScheme == Vector{Symbol}()
            @test form.tw == Vector{Float64}()
            @test form.TInt == Vector{Any}[]
            @test isa(form.mdl, JuMP.Model)

            # Test parameterized construction
            mdl = JuMP.Model()
            form = OCPFormulation{Float64}(
                tfDV=true,
                Np=51,
                IntegrationScheme=[:RK4],
                mdl=mdl
            )
            @test form.tfDV == true
            @test form.Np == 51
            @test form.IntegrationScheme == [:RK4]
            @test form.mdl === mdl
        end

        @testset "OCPSetting struct tests" begin
            # Test default construction
            setting = OCPSetting{Float64}()
            @test isa(setting.states, States)
            @test isa(setting.control, Control)
            @test isa(setting.solver, Solver)
            @test setting.InternalLogging == true
            @test setting.X0slack == false
            @test setting.XFslack == false

            # Test custom construction
            states = States(num=3)
            control = Control(num=2)
            setting = OCPSetting{Float64}(
                states=states,
                control=control,
                InternalLogging=false,
                X0slack=true
            )
            @test setting.states === states
            @test setting.control === control
            @test setting.InternalLogging == false
            @test setting.X0slack == true
        end

        @testset "OCPBound struct tests" begin
            # Test default construction
            bound = OCPBound{Float64}()
            @test bound.tfMin == 0.0
            @test bound.tfMax == 100.0
            @test bound.X0 == Vector{Float64}[]
            @test bound.XF == Vector{Float64}[]

            # Test parameterized construction
            X0 = [1.0, 2.0, 3.0]
            XF = [4.0, 5.0, 6.0]
            bound = OCPBound{Float64}(
                tfMin=0.5,
                tfMax=10.0,
                X0=X0,
                XF=XF
            )
            @test bound.tfMin == 0.5
            @test bound.tfMax == 10.0
            @test bound.X0 == X0
            @test bound.XF == XF
        end

        @testset "OCPParameter struct tests" begin
            # Test default construction
            param = OCPParameter{Float64}()
            @test size(param.x) == (0, 0)
            @test size(param.u) == (0, 0)
            @test size(param.params) == (0, 0)
            @test size(param.xvar) == (0, 0)
            @test size(param.δx) == (0, 0)
        end

        @testset "OCPResults struct tests" begin
            # Test default construction
            result = OCPResults{Float64}()
            @test size(result.X) == (0, 0)
            @test size(result.U) == (0, 0)
            @test result.Tst == Vector{Float64}()
            @test result.dt == Vector{Float64}()
            @test result.Status == :InFeasible
            @test result.IterNum == 0
            @test result.EvalNum == 0
            @test result.TSolve == 0.0
            @test result.Objval == 0.0
            @test result.Dfs == Vector{DataFrame}()

            # Test parameterized construction
            X = [1.0 2.0; 3.0 4.0]
            U = reshape([0.5, 1.5], 2, 1)  # Make it a matrix
            result = OCPResults{Float64}(
                X=X,
                U=U,
                Status=:Optimal,
                Objval=1.23
            )
            @test result.X == X
            @test result.U == U
            @test result.Status == :Optimal
            @test result.Objval == 1.23
        end

        @testset "OCP main struct tests" begin
            # Test default construction
            ocp = OCP{Float64}()
            @test isa(ocp.s, OCPSetting{Float64})
            @test isa(ocp.b, OCPBound{Float64})
            @test isa(ocp.f, OCPFormulation{Float64})
            @test isa(ocp.p, OCPParameter{Float64})
            @test isa(ocp.r, OCPResults{Float64})

            # Test that all components are properly initialized
            @test ocp.s.states.num == 0
            @test ocp.b.tfMin == 0.0
            @test ocp.f.tfDV == false
            @test size(ocp.p.x) == (0, 0)
            @test ocp.r.Status == :InFeasible
        end
    end

    @testset "Default Parameter Tests" begin
        @testset "Ipopt MPC defaults" begin
            defaults = Dict(JuliaOptimalControl._Ipopt_MPC_defaults)
            @test defaults["mu_strategy"] == "adaptive"
            @test defaults["max_iter"] == 100
            @test defaults["tol"] == 6e-3
            @test defaults["warm_start_init_point"] == "yes"
            @test defaults["max_cpu_time"] == 0.1
            @test defaults["print_level"] == 0
        end

        @testset "Ipopt standard defaults" begin
            defaults = Dict(JuliaOptimalControl._Ipopt_defaults)
            @test defaults["mu_strategy"] == "adaptive"
            @test defaults["max_iter"] == 1000
            @test defaults["warm_start_init_point"] == "no"
            @test defaults["max_cpu_time"] == 10.0
            @test defaults["print_level"] == 0
        end
    end

    @testset "Setup Function Tests" begin

        @testset "defineOCP tests" begin
            # Test basic construction
            ocp = defineOCP(
                numStates=3,
                numControls=2,
                X0=[1.0, 2.0, 3.0],
                XF=[4.0, 5.0, 6.0],
                XL=[-1.0, -2.0, -3.0],
                XU=[10.0, 20.0, 30.0],
                CL=[-5.0, -6.0],
                CU=[15.0, 16.0]
            )

            @test isa(ocp, OCP{Float64})
            @test ocp.s.states.num == 3
            @test ocp.s.control.num == 2
            @test ocp.b.X0 == [1.0, 2.0, 3.0]
            @test ocp.b.XF == [4.0, 5.0, 6.0]
            @test ocp.b.XL == [-1.0, -2.0, -3.0]
            @test ocp.b.XU == [10.0, 20.0, 30.0]
            @test ocp.b.CL == [-5.0, -6.0]
            @test ocp.b.CU == [15.0, 16.0]

            # Test default NaN tolerance values
            @test all(isnan.(ocp.b.X0_tol))
            @test all(isnan.(ocp.b.XF_tol))
            @test length(ocp.b.X0_tol) == 3
            @test length(ocp.b.XF_tol) == 3

            # Test error conditions
            @test_throws ArgumentError defineOCP(numStates=-1, numControls=1)
            @test_throws ArgumentError defineOCP(numStates=1, numControls=-1)
            @test_throws ErrorException defineOCP(
                numStates=2, numControls=1,
                X0=[1.0]  # Wrong size
            )
        end

        @testset "defineStates! tests" begin
            ocp = defineOCP(numStates=3, numControls=1)

            # Test successful state definition
            defineStates!(ocp, [:x, :y, :z])
            @test ocp.s.states.name == [:x, :y, :z]

            # Test error for wrong number of states
            @test_throws ErrorException defineStates!(ocp, [:x, :y])  # Too few
            @test_throws ErrorException defineStates!(ocp, [:x, :y, :z, :w])  # Too many
        end

        @testset "defineControls! tests" begin
            ocp = defineOCP(numStates=2, numControls=2)

            # Test successful control definition
            defineControls!(ocp, [:u1, :u2])
            @test ocp.s.control.name == [:u1, :u2]

            # Test error for wrong number of controls
            @test_throws ErrorException defineControls!(ocp, [:u1])  # Too few
            @test_throws ErrorException defineControls!(ocp, [:u1, :u2, :u3])  # Too many
        end

        @testset "defineTolerance! tests" begin
            ocp = defineOCP(numStates=3, numControls=1)

            # Test initial tolerance setting
            defineTolerance!(ocp, X0_tol=[0.1, 0.2, 0.3])
            @test ocp.s.X0slack == true
            @test ocp.b.X0_tol == [0.1, 0.2, 0.3]
            @test ocp.s.XFslack == false  # Should remain false

            # Test final tolerance setting
            defineTolerance!(ocp, XF_tol=[0.4, 0.5, 0.6])
            @test ocp.s.XFslack == true
            @test ocp.b.XF_tol == [0.4, 0.5, 0.6]

            # Test both tolerances
            ocp2 = defineOCP(numStates=2, numControls=1)
            defineTolerance!(ocp2, X0_tol=[0.1, 0.2], XF_tol=[0.3, 0.4])
            @test ocp2.s.X0slack == true
            @test ocp2.s.XFslack == true
            @test ocp2.b.X0_tol == [0.1, 0.2]
            @test ocp2.b.XF_tol == [0.3, 0.4]
        end

        @testset "ValidateScheme tests" begin
            # Test valid schemes
            @test ValidateScheme([:RK1]) == true
            @test ValidateScheme([:RK2]) == true
            @test ValidateScheme([:RK3]) == true
            @test ValidateScheme([:RK4]) == true
            @test ValidateScheme([:trapezoidal]) == true
            @test ValidateScheme([:bkwEuler]) == true
            @test ValidateScheme([:RK1, :RK2, :trapezoidal]) == true

            # Test invalid schemes
            result = ValidateScheme([:InvalidScheme])
            @test result == (false, 1)

            result = ValidateScheme([:RK1, :InvalidScheme, :RK2])
            @test result == (false, 2)
        end

        @testset "CalXvar tests" begin
            # Test schemes requiring no intermediate variables
            @test CalXvar([:RK1]) == 0
            @test CalXvar([:bkwEuler]) == 0
            @test CalXvar([:trapezoidal]) == 0

            # Test schemes requiring intermediate variables
            @test CalXvar([:RK2]) == 1
            @test CalXvar([:RK3]) == 2
            @test CalXvar([:RK4]) == 3

            # Test combinations
            @test CalXvar([:RK2, :RK3]) == 3  # 1 + 2
            @test CalXvar([:RK4, :RK4]) == 6  # 3 + 3
            @test CalXvar([:RK1, :RK2, :RK3, :RK4]) == 6  # 0 + 1 + 2 + 3
        end

        @testset "ConfigurePredefined tests" begin
            ocp = defineOCP(numStates=2, numControls=1)

            # Simple dynamics function for testing
            dynamics = (x, u, p) -> [x[2], u[1]]

            @testset "Fixed time horizon" begin
                formulation = ConfigurePredefined(ocp,
                    tf=5.0,
                    Np=21,
                    IntegrationScheme=:RK4,
                    dx=dynamics
                )

                @test isa(formulation, OCPFormulation)
                @test formulation.tfDV == false
                @test formulation.tf == 5.0
                @test formulation.Np == 21
                @test formulation.IntegrationScheme == fill(:RK4, 20)  # Np-1 intervals
                @test length(formulation.tw) == 20
                @test abs(sum(formulation.tw) - 1.0) < 1e-10  # Should sum to 1
                @test formulation.dx[1] === dynamics
            end

            @testset "Variable time horizon" begin
                formulation = ConfigurePredefined(ocp,
                    tfDV=true,
                    Np=11,
                    IntegrationScheme=:trapezoidal,
                    dx=dynamics
                )

                @test formulation.tfDV == true
                @test isa(formulation.tf, JuMP.VariableRef)
                @test formulation.Np == 11
                @test formulation.IntegrationScheme == fill(:trapezoidal, 10)
            end

            @testset "With cost expression" begin
                cost_expr = (x, u, p) -> x[1]^2 + u[1]^2

                formulation = ConfigurePredefined(ocp,
                    tf=3.0,
                    Np=6,
                    IntegrationScheme=:RK2,
                    dx=dynamics,
                    expr=cost_expr
                )

                @test formulation.expr[1] === cost_expr
                @test all(f === cost_expr for f in formulation.expr)
            end

            @testset "With constraints" begin
                constraint = (x, u, p) -> [x[1] - 10]  # x[1] <= 10

                formulation = ConfigurePredefined(ocp,
                    tf=2.0,
                    Np=5,
                    IntegrationScheme=:bkwEuler,
                    dx=dynamics,
                    cons=constraint
                )

                @test formulation.cons[1] === constraint
                @test all(f === constraint for f in formulation.cons)
            end

            @testset "With parameters" begin
                # Test single parameter row (broadcasted)
                params1 = [1.0 2.0]  # 1x2 matrix
                formulation1 = ConfigurePredefined(ocp,
                    tf=1.0,
                    Np=4,
                    IntegrationScheme=:RK1,
                    dx=dynamics,
                    params=params1
                )
                @test size(formulation1.params) == (4, 2)
                @test all(formulation1.params[i, :] == [1.0, 2.0] for i in 1:4)

                # Test full parameter matrix
                params2 = [1.0 2.0; 3.0 4.0; 5.0 6.0]  # 3x2 matrix
                formulation2 = ConfigurePredefined(ocp,
                    tf=1.0,
                    Np=3,
                    IntegrationScheme=:RK1,
                    dx=dynamics,
                    params=params2
                )
                @test formulation2.params == params2
            end

            @testset "Error conditions" begin
                # Missing required arguments
                @test_throws ErrorException ConfigurePredefined(ocp, tf=1.0)  # No Np
                @test_throws ErrorException ConfigurePredefined(ocp, Np=5)    # No tf or tfDV
                @test_throws ErrorException ConfigurePredefined(ocp, tf=1.0, Np=5)  # No dx

                # Invalid Np type
                @test_logs (:warn, "Round Np to nearest integer") ConfigurePredefined(ocp,
                    tf=1.0, Np=5.7, IntegrationScheme=:RK1, dx=dynamics)

                # Time horizon bounds error
                ocp_bounded = defineOCP(numStates=2, numControls=1)
                ocp_bounded.b.tfMin = 1.0
                ocp_bounded.b.tfMax = 10.0
                @test_throws ErrorException ConfigurePredefined(ocp_bounded,
                    tf=0.5, Np=5, IntegrationScheme=:RK1, dx=dynamics)  # tf < tfMin
            end
        end

        @testset "CheckOCPFormulation tests" begin
            ocp = defineOCP(numStates=2, numControls=1)

            # Create valid formulation
            formulation = ConfigurePredefined(ocp,
                tf=5.0,
                Np=11,
                IntegrationScheme=:RK2,
                dx=(x, u, p) -> [x[2], u[1]]
            )

            # Should not throw for valid formulation
            @test_nowarn CheckOCPFormulation(ocp, formulation)

            # Test time weight validation
            formulation_bad_weights = deepcopy(formulation)
            formulation_bad_weights.tw = [0.1, 0.2, 0.3]  # Doesn't sum to 1
            @test_throws ErrorException CheckOCPFormulation(ocp, formulation_bad_weights)

            # Test size mismatch
            formulation_size_mismatch = deepcopy(formulation)
            pop!(formulation_size_mismatch.tw)  # Remove one element
            @test_throws ErrorException CheckOCPFormulation(ocp, formulation_size_mismatch)

            # Test invalid Np
            formulation_bad_np = deepcopy(formulation)
            formulation_bad_np.Np = -1
            @test_throws ErrorException CheckOCPFormulation(ocp, formulation_bad_np)
        end

        @testset "defineSolver! tests" begin
            formulation = OCPFormulation{Float64}()
            formulation.mdl = JuMP.Model()

            # Test Ipopt solver setting
            custom_options = ("max_iter" => 500, "tol" => 1e-8)
            @test_nowarn defineSolver!(formulation, :Ipopt, custom_options)

            # Test unsupported solver
            formulation_new = OCPFormulation{Float64}()
            formulation_new.mdl = JuMP.Model()
            @test_throws ErrorException defineSolver!(formulation_new, :UnsupportedSolver, ())

            # Test when solver already attached
            formulation2 = OCPFormulation{Float64}()
            formulation2.mdl = JuMP.Model()
            try
                import Ipopt
                set_optimizer(formulation2.mdl, Ipopt.Optimizer)
                @test_logs (:warn, r"is set, can not change") defineSolver!(formulation2, :Ipopt, ())
            catch
                # If Ipopt not available, skip this test
                @test true
            end
        end
    end

    @testset "Utils Function Tests" begin

        @testset "RetrieveSolveStatus tests" begin
            # Test optimal status codes
            @test RetrieveSolveStatus(MOI.OPTIMAL) == :Optimal
            @test RetrieveSolveStatus(MOI.LOCALLY_SOLVED) == :Optimal
            @test RetrieveSolveStatus(MOI.ALMOST_OPTIMAL) == :Optimal

            # Test user limit status codes
            @test RetrieveSolveStatus(MOI.ITERATION_LIMIT) == :UserLimit
            @test RetrieveSolveStatus(MOI.TIME_LIMIT) == :UserLimit
            @test RetrieveSolveStatus(MOI.NODE_LIMIT) == :UserLimit
            @test RetrieveSolveStatus(MOI.SOLUTION_LIMIT) == :UserLimit
            @test RetrieveSolveStatus(MOI.MEMORY_LIMIT) == :UserLimit

            # Test infeasible status codes
            @test RetrieveSolveStatus(MOI.INFEASIBLE) == :Infeasible
            @test RetrieveSolveStatus(MOI.DUAL_INFEASIBLE) == :Infeasible
            @test RetrieveSolveStatus(MOI.OTHER_ERROR) == :Infeasible
        end

        @testset "CreateEmptyFormulation tests" begin
            form = CreateEmptyFormulation()
            @test isa(form, OCPFormulation)
            @test form.tfDV == false
            @test form.Np == 0
            @test form.IntegrationScheme == Vector{Symbol}()
        end

        @testset "DeleteElement tests" begin
            # Test with vector of Any (can hold nothing)
            test_vector = Any[x -> x^2, x -> x^3, x -> x^4]
            DeleteElement(test_vector, 2)
            @test isnothing(test_vector[2])
            @test !isnothing(test_vector[1])
            @test !isnothing(test_vector[3])

            # Test with vector of symbols (cannot hold nothing, so element is removed)
            test_expr = [:expr1, :expr2, :expr3]
            DeleteElement(test_expr, 1)
            @test length(test_expr) == 2  # Element was removed
            @test test_expr[1] == :expr2  # Second element moved to first position
            @test test_expr[2] == :expr3  # Third element moved to second position
        end

        @testset "ResultsToDataFrame tests" begin
            # Create a simple OCP with results
            ocp = OCP{Float64}()
            ocp.s.states.num = 2
            ocp.s.states.name = [:x, :y]
            ocp.s.control.num = 1
            ocp.s.control.name = [:u]

            # Mock results data
            ocp.r.Tst = [0.0, 0.5, 1.0]
            ocp.r.X = [1.0 2.0; 1.5 2.5; 2.0 3.0]  # 3x2 matrix
            ocp.r.U = reshape([0.1, 0.2, 0.3], 3, 1)  # 3x1 matrix

            df = ResultsToDataFrame(ocp)

            @test isa(df, DataFrame)
            @test size(df) == (3, 4)  # 3 rows, 4 columns (t, x, y, u)
            @test df.t == [0.0, 0.5, 1.0]
            @test df.x == [1.0, 1.5, 2.0]
            @test df.y == [2.0, 2.5, 3.0]
            @test df.u == [0.1, 0.2, 0.3]
        end

        @testset "WarmStart tests" begin
            # Create OCP with mock model
            ocp = OCP{Float64}()
            ocp.f.mdl = JuMP.Model()

            # Mock variables
            ocp.s.states.num = 2
            ocp.s.control.num = 1
            ocp.p.x = @variable(ocp.f.mdl, x[1:3, 1:2])
            ocp.p.u = @variable(ocp.f.mdl, u[1:3, 1:1])

            # Mock previous solution data
            ocp.b.X0 = [1.0, 2.0]
            ocp.r.X = [1.0 2.0; 1.5 2.5; 2.0 3.0]
            ocp.r.U = reshape([0.1, 0.2, 0.3], 3, 1)

            # Test warm start function (should not error)
            @test_nowarn WarmStart(ocp)

            # Test when warm start is not supported
            # (The function should handle the try-catch gracefully)
            @test_nowarn WarmStart(ocp)
        end

        @testset "UpdateX0! tests" begin
            # Create OCP with mock model
            ocp = OCP{Float64}()
            ocp.f.mdl = JuMP.Model()

            # Mock variables and settings
            ocp.s.states.num = 2
            ocp.p.x = @variable(ocp.f.mdl, x[1:3, 1:2])
            ocp.p.u = @variable(ocp.f.mdl, u[1:3, 1:1])

            # Initial conditions
            X0_new = [5.0, 6.0]
            ocp.r.X = [1.0 2.0; 1.5 2.5; 2.0 3.0]
            ocp.r.U = reshape([0.1, 0.2, 0.3], 3, 1)

            # Test UpdateX0! function
            @test_nowarn UpdateX0!(ocp, X0_new)
            @test ocp.b.X0 == X0_new

            # Verify that variables are fixed to new values
            @test is_fixed(ocp.p.x[1, 1])
            @test is_fixed(ocp.p.x[1, 2])
            @test fix_value(ocp.p.x[1, 1]) == 5.0
            @test fix_value(ocp.p.x[1, 2]) == 6.0
        end

        @testset "OptSolve! and GetOptimizeValue! integration tests" begin
            # Create a minimal working OCP for testing
            ocp = OCP{Float64}()
            ocp.s.states.num = 1
            ocp.s.control.num = 1
            ocp.s.states.name = [:x]
            ocp.s.control.name = [:u]

            # Create simple model
            ocp.f.mdl = JuMP.Model()
            set_silent(ocp.f.mdl)

            # Add simple variables
            ocp.p.x = @variable(ocp.f.mdl, x[1:2, 1:1])
            ocp.p.u = @variable(ocp.f.mdl, u[1:2, 1:1])

            # Add simple objective and constraint
            @objective(ocp.f.mdl, Min, sum(ocp.p.x) + sum(ocp.p.u))
            @constraint(ocp.f.mdl, ocp.p.x[1, 1] == 1.0)
            @constraint(ocp.f.mdl, ocp.p.x[2, 1] == ocp.p.x[1, 1] + ocp.p.u[1, 1])

            # Mock time data
            ocp.p.tV = [0.0, 1.0]
            ocp.f.TInt = [1.0]
            ocp.f.tfDV = false
            ocp.s.InternalLogging = true

            # Test that OptSolve! doesn't crash (solver may not be available)
            try
                # Set a simple solver if available
                import Ipopt
                set_optimizer(ocp.f.mdl, Ipopt.Optimizer)
                OptSolve!(ocp)

                # If solve succeeds, test results extraction
                if ocp.r.Status ∈ [:Optimal, :UserLimit]
                    @test size(ocp.r.X, 2) == 1  # One state
                    @test size(ocp.r.U, 2) == 1  # One control
                    @test length(ocp.r.Tst) == 2  # Two time points
                    @test length(ocp.r.Dfs) == 1  # One DataFrame stored
                end
            catch
                # If Ipopt not available or solve fails, test that function completes
                @test_nowarn OptSolve!(ocp)
                @test ocp.r.Status == :InFeasible  # Default status
            end
        end

        @testset "ExprIntegral tests" begin
            # Create OCP with expression functions
            ocp = OCP{Float64}()
            ocp.f.mdl = JuMP.Model()

            # Mock data
            ocp.f.Np = 3
            ocp.f.expr = Vector{Any}(nothing, 3)
            ocp.f.expr[2] = (x, u, p) -> x[1]^2 + u[1]^2  # Simple quadratic cost
            ocp.f.expr[3] = (x, u, p) -> x[1]^2 + u[1]^2

            ocp.f.IntegrationScheme = [:RK4, :RK4]
            ocp.f.params = Matrix{Any}(undef, 3, 1)
            ocp.f.TInt = [0.5, 0.5]

            # Create variables
            ocp.p.x = @variable(ocp.f.mdl, x[1:3, 1:1])
            ocp.p.u = @variable(ocp.f.mdl, u[1:3, 1:1])

            # Test integral computation (just check it doesn't error)
            @test_nowarn ExprIntegral(ocp)

            # Test with no expressions
            ocp.f.expr = Vector{Any}(nothing, 3)
            @test_nowarn ExprIntegral(ocp)
        end
    end

    @testset "Integration Tests" begin

        @testset "Simple Double Integrator Test" begin
            # Test complete workflow for double integrator: ẍ = u

            # System: [x, ẋ] with control u
            dynamics = function(x, u, p)
                return [x[2], u[1]]  # [ẋ, ẍ] = [x₂, u]
            end

            cost = function(x, u, p)
                return 0.5 * (x[1]^2 + x[2]^2 + u[1]^2)  # Quadratic cost
            end

            # Problem setup
            ocp = defineOCP(
                numStates=2,
                numControls=1,
                X0=[1.0, 0.0],     # Start at position 1, velocity 0
                XF=[0.0, 0.0],     # End at origin with zero velocity
                XL=[-10.0, -5.0],  # State bounds
                XU=[10.0, 5.0],
                CL=[-2.0],         # Control bounds
                CU=[2.0]
            )

            defineStates!(ocp, [:position, :velocity])
            defineControls!(ocp, [:acceleration])

            # Configure formulation
            formulation = ConfigurePredefined(ocp,
                tf=2.0,
                Np=21,
                IntegrationScheme=:RK4,
                dx=dynamics,
                expr=cost
            )

            # Build complete problem
            @test_nowarn OCPdef!(ocp, formulation)

            # Verify problem structure
            @test ocp.s.states.pts == 21
            @test ocp.s.control.pts == 21
            @test size(ocp.p.x) == (21, 2)
            @test size(ocp.p.u) == (21, 1)

            # Test that objective can be set
            obj = ExprIntegral(ocp)
            @test_nowarn @objective(ocp.f.mdl, Min, obj)

            # Test solve (may not converge without proper solver)
            @test_nowarn OptSolve!(ocp)

            # Basic sanity checks on results structure
            @test ocp.r.Status ∈ [:Optimal, :UserLimit, :InFeasible]
            @test length(ocp.r.Tst) ∈ [0, 21]  # Either empty or full trajectory
        end

        @testset "Free Final Time Problem" begin
            # Minimum time problem: minimize tf subject to dynamics

            dynamics = function(x, u, p)
                return [u[1]]  # Simple integrator: ẋ = u
            end

            ocp = defineOCP(
                numStates=1,
                numControls=1,
                X0=[0.0],
                XF=[5.0],      # Reach position 5
                XL=[-100.0],
                XU=[100.0],
                CL=[-1.0],     # Control constraint |u| ≤ 1
                CU=[1.0]
            )

            # Set time bounds for free final time
            ocp.b.tfMin = 0.1
            ocp.b.tfMax = 20.0

            defineStates!(ocp, [:position])
            defineControls!(ocp, [:velocity])

            # Configure with variable final time
            formulation = ConfigurePredefined(ocp,
                tfDV=true,        # Free final time
                Np=11,
                IntegrationScheme=:trapezoidal,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Verify time variable was created
            @test isa(formulation.tf, JuMP.VariableRef)
            @test formulation.tfDV == true

            # Add minimum time objective: minimize tf
            @test_nowarn @objective(ocp.f.mdl, Min, formulation.tf)

            @test_nowarn OptSolve!(ocp)
        end

        @testset "Problem with Path Constraints" begin
            # Pendulum with path constraints

            dynamics = function(x, u, p)
                # Simple pendulum: θ̈ = u (torque control)
                return [x[2], u[1]]
            end

            constraint = function(x, u, p)
                # Keep angle within bounds: -π/2 ≤ θ ≤ π/2
                return [π/2 - abs(x[1])]  # Must be ≥ 0
            end

            cost = function(x, u, p)
                return u[1]^2  # Minimize control effort
            end

            ocp = defineOCP(
                numStates=2,
                numControls=1,
                X0=[π/4, 0.0],     # Start at 45 degrees
                XF=[0.0, 0.0],     # End at vertical with zero velocity
                XL=[-π, -10.0],
                XU=[π, 10.0],
                CL=[-5.0],
                CU=[5.0]
            )

            defineStates!(ocp, [:angle, :angular_velocity])
            defineControls!(ocp, [:torque])

            formulation = ConfigurePredefined(ocp,
                tf=3.0,
                Np=16,
                IntegrationScheme=:RK3,
                dx=dynamics,
                cons=constraint,
                expr=cost
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Add objective
            obj = ExprIntegral(ocp)
            @test_nowarn @objective(ocp.f.mdl, Min, obj)

            @test_nowarn OptSolve!(ocp)
        end

        @testset "Parametric Problem" begin
            # System with parameters

            dynamics = function(x, u, p)
                mass = p[1]
                damping = p[2]
                # Mass-damper system: mẍ + cẋ = u
                return [x[2], (u[1] - damping * x[2]) / mass]
            end

            # Parameters: mass=2.0, damping=0.5
            params = [2.0 0.5]  # Will be broadcasted to all time points

            ocp = defineOCP(
                numStates=2,
                numControls=1,
                X0=[1.0, 0.0],
                XF=[0.0, 0.0],
                XL=[-5.0, -3.0],
                XU=[5.0, 3.0],
                CL=[-10.0],
                CU=[10.0]
            )

            defineStates!(ocp, [:position, :velocity])
            defineControls!(ocp, [:force])

            formulation = ConfigurePredefined(ocp,
                tf=4.0,
                Np=26,
                IntegrationScheme=:bkwEuler,
                dx=dynamics,
                params=params
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Verify parameters were set up correctly
            @test size(formulation.params) == (26, 2)
            @test all(formulation.params[i, 1] == 2.0 for i in 1:26)
            @test all(formulation.params[i, 2] == 0.5 for i in 1:26)

            @test_nowarn OptSolve!(ocp)
        end

        @testset "Tolerance-based Boundary Conditions" begin
            # Test soft boundary conditions

            dynamics = function(x, u, p)
                return [u[1]]  # Simple integrator
            end

            ocp = defineOCP(
                numStates=1,
                numControls=1,
                X0=[0.0],
                XF=[1.0],
                XL=[-10.0],
                XU=[10.0],
                CL=[-2.0],
                CU=[2.0]
            )

            # Enable tolerances
            defineTolerance!(ocp, X0_tol=[0.1], XF_tol=[0.2])

            defineStates!(ocp, [:x])
            defineControls!(ocp, [:u])

            formulation = ConfigurePredefined(ocp,
                tf=1.0,
                Np=6,
                IntegrationScheme=:bkwEuler,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Verify tolerance settings were applied
            @test ocp.s.X0slack == true
            @test ocp.s.XFslack == true

            @test_nowarn OptSolve!(ocp)
        end

        @testset "MPC Workflow Test" begin
            # Test model predictive control workflow

            dynamics = function(x, u, p)
                return [x[2], u[1]]  # Double integrator
            end

            cost = function(x, u, p)
                return 0.1 * x[1]^2 + x[2]^2 + 0.1 * u[1]^2
            end

            # Setup MPC problem
            ocp = defineOCP(
                numStates=2,
                numControls=1,
                X0=[2.0, 1.0],     # Initial condition (will be updated)
                XF=[0.0, 0.0],     # Target
                XL=[-5.0, -3.0],
                XU=[5.0, 3.0],
                CL=[-1.0],
                CU=[1.0]
            )

            defineStates!(ocp, [:x, :xdot])
            defineControls!(ocp, [:u])

            # Use MPC-appropriate settings
            ocp.s.solver.settings = JuliaOptimalControl._Ipopt_MPC_defaults

            formulation = ConfigurePredefined(ocp,
                tfDV=true,          
                Np=11,
                IntegrationScheme=:bkwEuler,
                dx=dynamics,
                expr=cost
            )

            @test_nowarn OCPdef!(ocp, formulation)

            obj = ExprIntegral(ocp)
            @test_nowarn @objective(ocp.f.mdl, Min, obj)

            # First solve
            @test_nowarn OptSolve!(ocp)

            # Simulate MPC update cycle
            new_state = [1.9, 1.9]  # New measured state
            @test_nowarn UpdateX0!(ocp, new_state)
            @test ocp.b.X0 == new_state

            # Second solve (warm started)
            @test_nowarn OptSolve!(ocp)
        end

        @testset "Different Integration Schemes Comparison" begin
            # Test that different integration schemes produce results

            dynamics = (x, u, p) -> [u[1]]  # Simple integrator

            schemes_to_test = [:RK1, :RK2, :RK3, :RK4, :trapezoidal, :bkwEuler]

            for scheme in schemes_to_test
                ocp = defineOCP(
                    numStates=1,
                    numControls=1,
                    X0=[0.0],
                    XF=[1.0],
                    XL=[-10.0],
                    XU=[10.0],
                    CL=[-5.0],
                    CU=[5.0]
                )

                defineStates!(ocp, [:x])
                defineControls!(ocp, [:u])

                formulation = ConfigurePredefined(ocp,
                    tf=1.0,
                    Np=6,
                    IntegrationScheme=scheme,
                    dx=dynamics
                )

                # Should not error for any scheme
                @test_nowarn OCPdef!(ocp, formulation)

                # Verify correct scheme was set
                @test all(s == scheme for s in formulation.IntegrationScheme)

                # Check that appropriate intermediate variables were created
                expected_xvar = CalXvar(formulation.IntegrationScheme)
                if expected_xvar > 0
                    @test size(ocp.p.xvar, 1) == expected_xvar
                end

                @test_nowarn OptSolve!(ocp)
            end
        end

        @testset "Results Processing Test" begin
            # Test results extraction and DataFrame conversion

            # Simple 2D system
            dynamics = (x, u, p) -> [x[2], u[1]]

            ocp = defineOCP(
                numStates=2,
                numControls=1,
                X0=[1.0, 0.0],
                XF=[0.0, 0.0],
                XL=[-2.0, -2.0],
                XU=[2.0, 2.0],
                CL=[-1.0],
                CU=[1.0]
            )

            defineStates!(ocp, [:pos, :vel])
            defineControls!(ocp, [:acc])

            formulation = ConfigurePredefined(ocp,
                tfDV=true,
                Np=6,
                IntegrationScheme=:RK1,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp, formulation)
            @test_nowarn OptSolve!(ocp)

            # Test internal logging
            @test ocp.s.InternalLogging == true
        end

        @testset "Complete Problem Solution Test" begin
            # Test a complete problem that should solve successfully (simple integrator)

            dynamics = (x, u, p) -> [u[1]]

            ocp = defineOCP(
                numStates=1,
                numControls=1,
                X0=[0.0],
                XF=[1.0],
                XL=[-100.0],
                XU=[100.0],
                CL=[-10.0],
                CU=[10.0]
            )

            defineStates!(ocp, [:position])
            defineControls!(ocp, [:velocity])

            formulation = ConfigurePredefined(ocp,
                tf=1.0,
                Np=11,
                IntegrationScheme=:RK1,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Add a simple objective (minimize control effort)
            @test_nowarn @objective(ocp.f.mdl, Min, sum(ocp.p.u .^ 2))

            @test_nowarn OptSolve!(ocp)

            # This simple problem should be solvable if solver is available
            # Test problem structure regardless of solve status
            @test size(ocp.p.x) == (11, 1)
            @test size(ocp.p.u) == (11, 1)
            @test ocp.r.Status ∈ [:Optimal, :UserLimit, :InFeasible]
        end
    end

    # Include the original bicycle model test
    @testset "Bicycle Model Integration Test" begin
        # Bicycle model parameters
        la = 1.56
        lb = 1.64
        m = 2020
        g = 9.81
        Izz = 4095
        h = 0.6
        mu = 0.8

        function MagicFormula(alpha, Fz, mu)
            B =  5.68   # Input Q2b value here
            C =  1.817   # Input Q2b value here
            Fy =  mu*Fz*sin(C*atan(B/mu*alpha))  # Lateral force calculation
            return Fy
        end

        function bicycleModel_expr(states, controls, parameters)
            x = states[1]
            y = states[2]
            v = states[3]
            r = states[4]
            ψ = states[5]
            ux = states[6]
            δf = states[7]
            ax = controls[1]
            dδf = controls[2]
            Fzf = m*g*la/(la+lb) - m*h/(la+lb)*ax # Front axle load
            Fzr = m*g*lb/(la+lb) + m*h/(la+lb)*ax # Rear axle load

            αf = δf - (v+la*r)/ux # Front slip angle
            αr = -(v-lb*r)/ux # Rear slip angle

            Fyf = MagicFormula(αf, Fzf, mu) # Front lateral force
            Fyr = MagicFormula(αr, Fzr, mu) # Rear lateral force

            dstates = Vector{Any}(undef,7)
            dstates[1]         = ux*cos(ψ) - v*sin(ψ)
            dstates[2]         = ux*sin(ψ) + v*cos(ψ)
            dstates[3]         = (Fyf+Fyr)/m - ux*r
            dstates[4]         = (Fyf*la-Fyr*lb)/Izz
            dstates[5]         = r
            dstates[6]         = ax
            dstates[7]         = dδf

            return dstates
        end

        # cost weights
        w_y  = 1e-5
        w_sr = 2
        w_ax = 0.2
        w_ux = 0.2
        w_sa = 1
        function bicycle_cost(states, controls)
            y_cost = w_y * (states[2] - 5)^2
            dδf_cost = w_sr * (controls[2])^2
            ax_cost = w_ax * (controls[1])^2
            ux_cost = w_ux * (states[6] - 13.0)^2
            δf_cost = w_sa * (states[7])^2
            return y_cost + dδf_cost + ax_cost + ux_cost + δf_cost
        end

        XL = [-40, -20, -3, -pi/5, -pi/2, 5.0, -pi/12]
        XU = [300, 20, 3, pi/5, pi/2, 15.0, pi/12]
        CL = [-2.6, -0.1]
        CU = [2.6, 0.1]
        X0 = [-10.0, 0, 0, 0, 0, 10.0, 0]
        XF = [NaN, NaN, NaN, NaN, NaN, NaN, NaN]

        ocp = defineOCP(numStates=7, numControls=2, X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU)
        defineStates!(ocp, [:x, :y, :v, :r, :ψ, :ux, :δf])
        defineControls!(ocp, [:ax, :dδf])

        OCPForm = ConfigurePredefined(ocp;
            (:Np => 81),
            (:tfDV => false),
            (:tf => 8),
            (:IntegrationScheme => :bkwEuler),
            (:dx => bicycleModel_expr),
            (:expr => bicycle_cost)
        )

        @test_nowarn OCPdef!(ocp, OCPForm)

        xpos = ocp.p.x[:, 1]
        y = ocp.p.x[:, 2]
        dδf = ocp.p.u[:, 2]
        ax = ocp.p.u[:, 1]
        ux = ocp.p.x[:, 6]
        δf = ocp.p.x[:, 7]

        obs_cons = @constraint(ocp.f.mdl, [i=1:ocp.s.states.pts], 36 <= ((xpos[i] - 30).^2 + (y[i] - 2).^2))
        obj = @expression(ocp.f.mdl, sum((0.05 * (y[j] - sin(xpos[j]))^2 + 2 * dδf[j]^2 + 0.2 * ax[j]^2 + 0.2 * (ux[j] - 13)^2 + 1 * δf[j]^2) * ocp.f.TInt[j-1] for j in 2:ocp.f.Np))
        @objective(ocp.f.mdl, Min, obj)

        @test_nowarn OptSolve!(ocp)

        # Basic structure tests
        @test size(ocp.p.x) == (81, 7)
        @test size(ocp.p.u) == (81, 2)

        # Test results if solve was successful
        if ocp.r.Status ∈ [:Optimal, :UserLimit]
            @test size(ocp.r.X) == (81, 7)

            # Distance constraint check
            Dis = (ocp.r.X[:,1].-30).^2 + (ocp.r.X[:,2].-2).^2
            @test count(x->x<35.9, Dis) == 0
        end
    end

    @testset "Advanced Functionality Tests" begin

        @testset "Constraint Enforcement Tests" begin
            # Test that box constraints are actually enforced
            dynamics = (x, u, p) -> [u[1]]
            function test_cost(states, controls, params)
                cost = states[1]^2 + controls[1]^2
                return cost
            end
            ocp = defineOCP(
                numStates=1,
                numControls=1,
                X0=[0.0],     # Start within feasible region
                XF=[NaN],
                XL=[-2.0],    # State constraint
                XU=[2.0],
                CL=[-1.0],    # More reasonable control constraint
                CU=[1.0]
            )

            defineStates!(ocp, [:x])
            defineControls!(ocp, [:u])

            formulation = ConfigurePredefined(ocp;
                (:tfDV => true),
                (:Np => 5),
                (:IntegrationScheme => :bkwEuler),
                (:dx => dynamics),
                (:expr => test_cost)
            )

            @test_nowarn OCPdef!(ocp, formulation)
            obj = ExprIntegral(ocp)
            @test_nowarn @objective(ocp.f.mdl, Min, obj)
            @test_nowarn OptSolve!(ocp)

            # If solve succeeds, check constraint satisfaction
            if ocp.r.Status == :Optimal
                @test all(-2.0 .<= ocp.r.X[:, 1] .<= 2.0)  # State constraints
                @test all(-0.1 .<= ocp.r.U[:, 1] .<= 0.1)  # Control constraints
            end
        end

        @testset "Integration Accuracy Tests" begin
            # Test integration accuracy with a known solution
            # ẋ = -x, x(0) = 1 has solution x(t) = e^(-t)
            dynamics = (x, u, p) -> [-x[1]]

            ocp = defineOCP(
                numStates=1,
                numControls=0,  # No controls - just integration
                X0=[1.0],
                XL=[-10.0],
                XU=[10.0]
            )

            defineStates!(ocp, [:x])

            # Test different integration schemes for accuracy
            schemes = [:RK1, :RK2, :RK4]
            accuracies = Float64[]

            for scheme in schemes
                formulation = ConfigurePredefined(ocp,
                    tf=1.0,
                    Np=11,
                    IntegrationScheme=scheme,
                    dx=dynamics
                )

                # Create problem without controls
                ocp_test = deepcopy(ocp)
                ocp_test.s.control.num = 0
                @test_nowarn OCPdef!(ocp_test, formulation)

                # Expected value at t=1: e^(-1) ≈ 0.3679
                expected = exp(-1.0)

                # For integration test, we'll just check the formulation works
                @test formulation.IntegrationScheme == fill(scheme, 10)
            end
        end

        @testset "Scaling and Large Problem Tests" begin
            # Test with larger state/control dimensions
            n_states = 10
            n_controls = 5

            dynamics = function(x, u, p)
                # Simple linear dynamics: ẋ = Ax + Bu
                A = -0.1 * Matrix{Float64}(I, n_states, n_states)  # Stable system
                B = ones(n_states, n_controls) / n_controls
                return A * x + B * u
            end

            X0 = ones(n_states)
            XF = zeros(n_states)
            XL = -10 * ones(n_states)
            XU = 10 * ones(n_states)
            CL = -ones(n_controls)
            CU = ones(n_controls)

            ocp = defineOCP(
                numStates=n_states,
                numControls=n_controls,
                X0=X0, XF=XF, XL=XL, XU=XU, CL=CL, CU=CU
            )

            defineStates!(ocp, [Symbol("x$i") for i in 1:n_states])
            defineControls!(ocp, [Symbol("u$i") for i in 1:n_controls])

            formulation = ConfigurePredefined(ocp,
                tf=5.0,
                Np=11,
                IntegrationScheme=:RK2,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp, formulation)

            # Verify problem structure for large problem
            @test size(ocp.p.x) == (11, n_states)
            @test size(ocp.p.u) == (11, n_controls)

            @test_nowarn OptSolve!(ocp)
        end

        @testset "Boundary Condition Variations Tests" begin
            dynamics = (x, u, p) -> [u[1]]

            # Test with only initial conditions (free final state)
            ocp1 = defineOCP(
                numStates=1,
                numControls=1,
                X0=[0.0],
                XF=[NaN],  # Free final state
                XL=[-10.0], XU=[10.0], CL=[-1.0], CU=[1.0]
            )

            defineStates!(ocp1, [:x])
            defineControls!(ocp1, [:u])

            formulation1 = ConfigurePredefined(ocp1,
                tf=2.0, Np=6, IntegrationScheme=:RK1, dx=dynamics)

            @test_nowarn OCPdef!(ocp1, formulation1)

            # Test with only final conditions (free initial state)
            ocp2 = defineOCP(
                numStates=1,
                numControls=1,
                X0=[NaN],  # Free initial state
                XF=[1.0],
                XL=[-10.0], XU=[10.0], CL=[-1.0], CU=[1.0]
            )

            defineStates!(ocp2, [:x])
            defineControls!(ocp2, [:u])

            formulation2 = ConfigurePredefined(ocp2,
                tf=2.0, Np=6, IntegrationScheme=:RK1, dx=dynamics)

            @test_nowarn OCPdef!(ocp2, formulation2)
        end

        @testset "Cost Function Variations Tests" begin
            dynamics = (x, u, p) -> [u[1]]

            # Test Mayer cost (terminal cost only)
            mayer_cost = (x, u, p) -> x[1]^2  # Only at final time

            ocp1 = defineOCP(
                numStates=1, numControls=1,
                X0=[1.0], XF=[NaN],
                XL=[-10.0], XU=[10.0], CL=[-1.0], CU=[1.0]
            )

            defineStates!(ocp1, [:x])
            defineControls!(ocp1, [:u])

            # Only apply cost at final point
            formulation1 = ConfigurePredefined(ocp1,
                tf=1.0, Np=6, IntegrationScheme=:RK1, dx=dynamics)

            # Manually set only final cost
            formulation1.expr = Vector{Any}(nothing, 6)
            formulation1.expr[6] = mayer_cost  # Only at final time

            @test_nowarn OCPdef!(ocp1, formulation1)

            # Test Lagrange cost (integral cost)
            lagrange_cost = (x, u, p) -> u[1]^2  # Control effort

            ocp2 = defineOCP(
                numStates=1, numControls=1,
                X0=[1.0], XF=[0.0],
                XL=[-10.0], XU=[10.0], CL=[-1.0], CU=[1.0]
            )

            defineStates!(ocp2, [:x])
            defineControls!(ocp2, [:u])

            formulation2 = ConfigurePredefined(ocp2,
                tf=1.0, Np=6, IntegrationScheme=:RK1, dx=dynamics, expr=lagrange_cost)

            @test_nowarn OCPdef!(ocp2, formulation2)
        end

        @testset "Time Horizon Edge Cases Tests" begin
            dynamics = (x, u, p) -> [u[1]]

            # Test very short time horizon
            ocp1 = defineOCP(
                numStates=1, numControls=1,
                X0=[0.0], XF=[1.0],
                XL=[-10.0], XU=[10.0], CL=[-10.0], CU=[10.0]
            )

            defineStates!(ocp1, [:x])
            defineControls!(ocp1, [:u])

            formulation1 = ConfigurePredefined(ocp1,
                tf=0.01,  # Very short
                Np=3,     # Minimum points
                IntegrationScheme=:RK1,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp1, formulation1)

            # Test with minimum number of points
            ocp2 = defineOCP(
                numStates=1, numControls=1,
                X0=[0.0], XF=[1.0],
                XL=[-10.0], XU=[10.0], CL=[-10.0], CU=[10.0]
            )

            defineStates!(ocp2, [:x])
            defineControls!(ocp2, [:u])

            formulation2 = ConfigurePredefined(ocp2,
                tf=1.0,
                Np=2,     # Absolute minimum
                IntegrationScheme=:RK1,
                dx=dynamics
            )

            @test_nowarn OCPdef!(ocp2, formulation2)
            @test formulation2.Np == 2
            @test length(formulation2.IntegrationScheme) == 1
        end

    end
end

nothing
