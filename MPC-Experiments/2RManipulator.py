import casadi as ca

# System dynamics
def dynamics(x, u):
    dx = u - x**2
    return dx

# Discretize using Runge-Kutta
def rk4_step(x, u, dt):
    k1 = dynamics(x, u)
    k2 = dynamics(x + dt/2 * k1, u)
    k3 = dynamics(x + dt/2 * k2, u)
    k4 = dynamics(x + dt * k3, u)
    x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    return x_next

# NMPC setup
N = 20  # prediction horizon
dt = 0.1  # time step
x0 = ca.MX(0)  # initial condition
x = x0
u = ca.MX.sym('u')  # control input

# Objective and constraints
Q = 1
R = 0.1
cost = 0
g = []  # constraints
u_lb = -1
u_ub = 1
for i in range(N):
    cost += Q*(x**2) + R*(u**2)
    x = rk4_step(x, u, dt)
    g.append(x)
    g.append(u)

# Define optimization variables
OPT_variables = ca.vertcat(u)
n_vars = OPT_variables.size()[0]

# Define constraints
g = ca.vertcat(*g)
n_constraints = g.size()[0]

# Set up the optimization problem
nlp = {'x':OPT_variables, 'f':cost, 'g':g}

# Set up solver
opts = {"ipopt.print_level": 0, "ipopt.tol": 1e-8, "print_time": 0}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

# Solve NMPC
u0 = [0]*N
lbx = [u_lb]*N
ubx = [u_ub]*N
lbg = [-ca.inf]*n_constraints
ubg = [ca.inf]*n_constraints
res = solver(x0=u0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)

print(res['x'])  # optimized control inputs over the prediction horizon

