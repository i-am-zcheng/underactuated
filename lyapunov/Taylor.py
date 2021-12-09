from sympy import *
import math

def Taylor_polynomial_sympy(function_expression, variable_list, evaluation_point, degree):
    """
    Mathematical formulation reference:
    https://math.libretexts.org/Bookshelves/Calculus/Supplemental_Modules_(Calculus)/Multivariable_Calculus/3%3A_Topics_in_Partial_Derivatives/Taylor__Polynomials_of_Functions_of_Two_Variables
    :param function_expression: Sympy expression of the function
    :param variable_list: list. All variables to be approximated (to be "Taylorized")
    :param evaluation_point: list. Coordinates, where the function will be expressed
    :param degree: int. Total degree of the Taylor polynomial
    :return: Returns a Sympy expression of the Taylor series up to a given degree, of a given multivariate expression, approximated as a multivariate polynomial evaluated at the evaluation_point
    """
    from sympy import factorial, Matrix, prod
    import itertools

    n_var = len(variable_list)
    point_coordinates = [(i, j) for i, j in (zip(variable_list, evaluation_point))]  # list of tuples with variables and their evaluation_point coordinates, to later perform substitution

    deriv_orders = list(itertools.product(range(degree + 1), repeat=n_var))  # list with exponentials of the partial derivatives
    deriv_orders = [deriv_orders[i] for i in range(len(deriv_orders)) if sum(deriv_orders[i]) <= degree]  # Discarding some higher-order terms
    n_terms = len(deriv_orders)
    deriv_orders_as_input = [list(sum(list(zip(variable_list, deriv_orders[i])), ())) for i in range(n_terms)]  # Individual degree of each partial derivative, of each term

    polynomial = 0
    for i in range(n_terms):
        partial_derivatives_at_point = function_expression.diff(*deriv_orders_as_input[i]).subs(point_coordinates)  # e.g. df/(dx*dy**2)
        denominator = prod([factorial(j) for j in deriv_orders[i]])  # e.g. (1! * 2!)
        distances_powered = prod([(Matrix(variable_list) - Matrix(evaluation_point))[j] ** deriv_orders[i][j] for j in range(n_var)])  # e.g. (x-x0)*(y-y0)**2
        polynomial += partial_derivatives_at_point / denominator * distances_powered
    return polynomial

import numpy as np

x0,x1,x2,x3 = symbols('x0 x1 x2 x3')

m = 1
M = 5
L = 2
g = -10
d = 1
u = 0
Sx = sin(x2)
Cx = cos(x2)
D = m*L*L*(M+m*(1-Cx**2))

f1 = x1
f2 = (1/D)*(-(m**2)*(L**2)*g*Cx*Sx + m*(L**2)*(m*L*(x3**2)*Sx - d*x1)) + m*L*L*(1/D)*u
f3 = x3
f4 = (1/D)*((m+M)*m*g*L*Sx - m*L*Cx*(m*L*(x3**2)*Sx - d*x1)) - m*L*Cx*(1/D)*u


x=[x0,x1,x2,x3]
mu=[0, 0, pi, 0]
mygee = Taylor_polynomial_sympy(f4,x,mu,3)
print(mygee)

f = 7*x1*(x2 - pi)**2/100 - x1/10 + 6*x2 - x3**2*(x2 - pi)/5 - 11*(x2 - pi)**3/5 - 6*pi

print(f.subs([(x2, pi)]))
