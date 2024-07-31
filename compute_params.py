import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.optimize import minimize
from matplotlib.pyplot import cm
import sys


# Define constants
d1 = float(sys.argv[1])
d2 = float(sys.argv[2])
epsilon = float(sys.argv[3])


### x is epsilon 1
def fff(x, alpha):

    p = 1/(1+np.exp(x))

    s1 = p*(1-p)/(1-2*p)**2

    s2 = alpha**2 * d1 + (1-alpha)**2 * d2

    s3 = (2*(1-p)**2/(1-2*p)**2)/(epsilon - x )**2

    s4 = alpha**2 + (1-alpha)**2

    return s1*s2 + s3*s4

# Define the objective function to minimize
def objective(x_alpha):
    x, alpha = x_alpha
    return fff(x, alpha)

# Define bounds for x and alpha
bounds = [(epsilon*0.05, epsilon*0.95), (0.05, 0.95)]

# Initial guess
x0 = np.array([0.5*epsilon, 0.5])

# Minimize the objective function
result = minimize(objective, x0, bounds=bounds)

# Get the coordinates of the minimum point
min_x = result.x[0]
min_alpha = result.x[1]
min_fff = result.fun

# Print the result
# print("Global minimum value of fff:", result.fun)
# print("Corresponding x value:", result.x[0])
# print("Corresponding alpha value:", result.x[1])
print(min_x)
print(min_alpha)



# Generate values for x and alpha
x_values = np.linspace(epsilon*0.3, epsilon*0.7, 20)
alpha_values = np.linspace(0, 1, 20)

# Create meshgrid for x and alpha
X, Alpha = np.meshgrid(x_values, alpha_values)

# Compute the values of fff for each combination of x and alpha
Z = fff(X, Alpha)

# # Plot the 3D surface
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Alpha, Z, cmap=cm.Blues)
# surf = ax.plot_surface(X, Alpha, Z, cmap=cm.Greys, alpha=1, zorder=0)
ax.set_xlabel(r'$\epsilon_1$')  # Label epsilon1 in math mode
ax.set_ylabel(r'$\alpha$')  # Label alpha in math mode
ax.set_zlabel('variance')

start = np.ceil(epsilon * 0.3 * 5) / 5  # Round up epsilon * 0.3 to the nearest multiple of 0.2
end = np.floor(epsilon * 0.7 * 5) / 5  # Round down epsilon * 0.7 to the nearest multiple of 0.2
x_ticks = np.arange(start, end + 0.01, 0.2)  # Add 0.01 to end to ensure the endpoint is included
ax.set_xticks(x_ticks)  # Set the x ticks using the generated values
ax.view_init(elev=6, azim=-41)
# plt.savefig("fig.pdf")
# plt.savefig("sigmod-submission/curves/fig.pdf")
# plt.savefig("sigmod-submission/curves/{}-{}-{}_fig.pdf".format(d1, d2, epsilon))

plt.close()


# Plot the curves
plt.style.use('default')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2
# for the varying plots 
plt.rcParams["figure.figsize"] = (6,5)
plt.rcParams['pdf.fonttype'] = 42

Z_alpha0 = fff(x_values, 0)
Z_alpha1 = fff(x_values, 1)
plt.style.use('default')
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2
# for the varying plots 
plt.rcParams["figure.figsize"] = (6,5)
plt.rcParams['pdf.fonttype'] = 42


# Plot the curve for alpha = 0
plt.plot(x_values, Z_alpha0, color='blue', label=r'$\alpha = 0$', linewidth=2.0)
# Plot the curve for alpha = 1
plt.plot(x_values, Z_alpha1, color='red', label=r'$\alpha = 1$', linewidth=2.0)

# Plot the curve for alpha = 0.5
alpha_05 = 0.5
Z_alpha05 = fff(x_values, alpha_05)
plt.plot(x_values, Z_alpha05, color='green', label=r'$\alpha = 0.5$', linewidth=2.0)

# Plot the horizontal line at the value of min_fff
plt.axhline(y=min_fff, color='gray', linestyle='--', label='Global minimum', linewidth=2.0)

# Add labels and legend
plt.xlabel(r'$\epsilon_1$', fontsize=20)  # Label epsilon1 in math mode
plt.xticks(x_ticks, fontsize=20)
plt.yticks(fontsize=20)
# plt.ylabel('L2 loss', fontsize=20)
plt.title(r"L2 loss when $d_u= {}$, $d_w = {}$, and $\epsilon = {}$".format(int(d1), int(d2), epsilon), 
    fontsize=15)
plt.subplots_adjust(left=0.15, bottom=0.15)
plt.legend(fontsize=15)
plt.savefig("sigmod-submission/curves/{}-{}-{}_same_plot.pdf".format(d1, d2, epsilon))

# plt.plot(x_values, Z_alpha0, color='blue', label=r'$\alpha = 0$')
# plt.xlabel(r'$\epsilon_1$')  # Label epsilon1 in math mode
# plt.ylabel('variance')
# plt.title('Curve for alpha = 0')
# plt.legend()
# plt.savefig("curve0.pdf")
# plt.close()

# plt.plot(x_values, Z_alpha1, color='blue', label=r'$\alpha = 1$')
# plt.xlabel(r'$\epsilon_1$')  # Label epsilon1 in math mode
# plt.ylabel('variance')
# plt.title('Curve for alpha = 1')
# plt.legend()
# plt.savefig("curve1.pdf")
# plt.close()

