from sympy import *


def monomToIndicies(monom):
    """
    Given a monomial, return the indices of the corresponding term.
    Example: The term x_2x_3 has the monom == (0, 0, 1, 1, ...)
             Returns 2, 3
    """
    # A monom may only have one term, means its a diagonal
    if sum(monom) == 1:
        for i in range(len(monom)):
            if monom[i] == 1:
                return i, i


    # Otherwise, there must only be 2 terms
    elif sum(monom) == 2:
        # Find the first non-zero index
        for i in range(len(monom)):
            # x^2 == x for binary variables, return the diagonal
            if monom[i] == 2:
                return i, i

            elif monom[i] == 1:
                # Found the first index
                break

        # Find the second non-zero index
        for j in range(i + 1, len(monom)):
            if monom[j] == 1:
                return i, j

    else:
        return None, None


def makeCoeffMatrix(poly, n):
    """
    Given a polynomial, return the n by n matrix of its coefficients.
    """
    coeffs = poly.coeffs()
    monoms = poly.monoms()

    # Sanity check, must be equal lengths
    if len(coeffs) != len(monoms):
        print("Can't make matrix")
        return None

    # Initialize nxn matrix
    m = [
        [0]*n for _ in range(n)
    ]

    for i in range(len(monoms)):
        x, y = monomToIndicies(monoms[i])
        if x is not None:
            m[x][y] += coeffs[i]

    return Matrix(m)
