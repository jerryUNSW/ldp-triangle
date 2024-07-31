import numpy as np
import sympy as sp
from scipy.optimize import root_scalar
import sys

# Define the function
def f(x, d, epsilon):
    p = 1 / (1 + np.exp(x))
    factor__ = (1 - p) / (1 - 2 * p)**2
    return factor__ * (p*d/4 + (1-p)/ (epsilon - x)**2 )

if __name__ == "__main__":
    # Check if correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python script_name.py d epsilon")
        sys.exit(1)

    # Get values of d and epsilon from command line arguments
    d = float(sys.argv[1])
    epsilon = float(sys.argv[2])

    # Define x as a symbolic variable
    x = sp.Symbol('x')

    # Define the function
    p = 1 / (1 + sp.exp(x))
    factor__ = (1 - p) / (1 - 2 * p)**2

    avg_f_expr = factor__ * (p * d / 4 + (1 - p) / (epsilon - x)**2)

    f_prime_expr = sp.diff(avg_f_expr, x)

    # Convert symbolic expression to a function
    f_prime_func = sp.lambdify(x, f_prime_expr)

    def f_prime_numeric(x):
        return f_prime_func(x)

    # Use a root-finding algorithm to find the roots of f'
    search_range = [0.01 * epsilon, 0.99 * epsilon]
    root_result = root_scalar(f_prime_numeric, bracket=search_range)

    print(root_result.root)
