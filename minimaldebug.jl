using JuMP
model = Model()
nx = @variable(model, states[1:8])
expr = @expression(model, nx[1] + (nx[7] * cos(nx[5])) - ((nx[3] + 1.5521 * nx[4]) * sin(nx[5])) * 0.5       )
x = expr
obj = @expression(model, x + 1.8)
