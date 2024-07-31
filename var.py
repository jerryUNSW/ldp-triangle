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
    # if len(sys.argv) != 3:
    #     print("Usage: python script_name.py d epsilon")
    #     sys.exit(1)

    # # Get values of d and epsilon from command line arguments
    # d = float(sys.argv[1])
    # epsilon = float(sys.argv[2])

    # Define x as a symbolic variable
    x = sp.Symbol('x')
    d = sp.Symbol('d')
    epsilon = sp.Symbol('epsilon')

    # Define the function
    p = 1 / (1 + sp.exp(x))

    factor__ = (1 - p) / (1 - 2 * p)**2

    avg_f_expr = factor__ * (p * d / 4 + (1 - p) / (epsilon - x)**2)

    # print("avg_f_expr = ", avg_f_expr)
    
    # Expand and simplify the expression
    simplified_expr = sp.simplify(avg_f_expr)

    print("simplified_expr = ", simplified_expr)

    expanded_expr = simplified_expr.subs(sp.sinh(x / 2)**2, (sp.exp(x) - 2 + sp.exp(-x)) / 4)

    print('expanded_expr = ', sp.simplify(expanded_expr))


    f_prime_expr = sp.diff(avg_f_expr, x)
    # print(f_prime_expr)
    # Factor the derivative
    factored_expr = sp.factor(f_prime_expr)
    # print("fprime = ", factored_expr)

    # Extract only the numerator
    numerator_expr, denominator_expr = factored_expr.as_numer_denom()

    ff = numerator_expr/sp.exp(x)

    # print("ff = ", ff)

    ff_prime_expr = sp.diff(ff, x)    

    # print("ff' = ", ff_prime_expr)

    # # Convert symbolic expression to a function
    # f_prime_func = sp.lambdify((x, d, epsilon), f_prime_expr)

    # # Define values
    # d_value = 3
    # epsilon_value = 2


    # def f_prime_numeric(x):
    #     return f_prime_func(x, d_value, epsilon_value)

 
    # x_values = np.arange(0.4 * epsilon_value, 0.6* epsilon_value, 0.01*epsilon_value)

    # # Print values
    # for x_val in x_values:
    #     result = f_prime_numeric(x_val, d_value, epsilon_value)
    #     print(f"f'({x_val}) =", result)


    # def find_root(func, target, left, right, tol=1e-6, max_iter=1000):
    #     """
    #     Find a root of a monotonously increasing function using binary search.

    #     Parameters:
    #         func (callable): The function for which to find the root.
    #         target (float): The target value that the function should be close to.
    #         left (float): The left boundary of the search interval.
    #         right (float): The right boundary of the search interval.
    #         tol (float): The tolerance for the root.
    #         max_iter (int): Maximum number of iterations.

    #     Returns:
    #         float: The approximate root found.
    #     """
    #     for _ in range(max_iter):
    #         mid = (left + right) / 2
    #         if func(mid) < target:
    #             left = mid
    #         else:
    #             right = mid
    #         if abs(func(mid) - target) < tol:
    #             return mid
    #     return None

    # # Usage
    # root = find_root(f_prime_numeric, 0, 0.05 * epsilon_value, 0.95 * epsilon_value, tol = 0.01)
    # print("Numeric root found:", root)

    # # Use a root-finding algorithm to find the roots of f'
    # search_range = [0.01 * epsilon, 0.99 * epsilon]
    # root_result = root_scalar(f_prime_numeric, bracket=search_range)

    # print(root_result.root)
