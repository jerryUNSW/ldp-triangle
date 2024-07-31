import sympy as sp

# Define the symbolic variables
t, p = sp.symbols('t p')

# Define the MGF of a Bernoulli random variable X
MGF_X = (1 - p) + p * sp.exp(t)

# Define the number of Bernoulli trials (d)
# this is the degree
d = sp.Symbol('d')

# Define the MGF of the sum Y = sum_{i=1}^d X_i
MGF_Y = MGF_X**d

# Simplify the expression (if needed)
MGF_Y_simplified = sp.simplify(MGF_Y)

# Print the MGF of Y in a readable format
print("MGF of Y (sum of d Bernoulli variables):")
print(MGF_Y_simplified)

# Print the MGF of Y in a LaTeX-like format
# print("LaTeX-like format:")
# sp.pretty_print(MGF_Y_simplified, use_unicode=True)



# need to compute E(\tilde(f)^2) based on 
# variance and MGF. 
