import numpy as np
import cvxpy as cp

import matplotlib.pyplot as plt

# Define system dynamics: x_next = A*x + B*u
A = np.array([[1]])
B = np.array([[1]])

# Define system constraints
umax = 1.0
umin = -1.0

# MPC setup
N = 5  # Prediction horizon
target = np.array([[10]])  # Target position

# Initialize system state
x0 = np.array([[0]])

# MPC loop
num_iterations = 20
x_traj = [x0]
u_traj = []

for _ in range(num_iterations):
    x = cp.Variable((1, N+1))
    u = cp.Variable((1, N))

    # Cost function
    cost = 0
    constraints = []
    
    for k in range(N):
        cost += cp.sum_squares(x[:, k] - target) + cp.sum_squares(u[:, k])
        constraints += [x[:, k+1] == A @ x[:, k] + B @ u[:, k],
                        umin <= u[:, k], u[:, k] <= umax]

    # Initial condition constraint
    constraints += [x[:, 0] == x0.flatten()]

    # Solve the optimization problem
    problem = cp.Problem(cp.Minimize(cost), constraints)
    problem.solve()

    # Apply the first control input to the system
    u0 = u[:, 0].value
    u_traj.append(u0)
    x0 = A @ x0 + B * u0
    x_traj.append(x0)

# Display results
print("State trajectory:", x_traj)
print("Control trajectory:", u_traj)

# Convert trajectories to arrays for plotting
x_values = np.array(x_traj).flatten()
u_values = np.array(u_traj).flatten()

# Plot state trajectory
plt.figure()
plt.plot(x_values, '-o', label='State (x)')
plt.axhline(y=target[0,0], color='r', linestyle='--', label='Target')
plt.title('State Trajectory')
plt.xlabel('Time step')
plt.ylabel('State')
plt.legend()
plt.grid(True)

# Plot control trajectory
plt.figure()
plt.plot(u_values, '-o', label='Control (u)')
plt.title('Control Trajectory')
plt.xlabel('Time step')
plt.ylabel('Control Input')
plt.grid(True)

plt.show()