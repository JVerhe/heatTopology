import numpy as np

# Adjoint method according to ZhangEtAl2007
def adjoint(T, v, corners, p=1):
    """
    Returns a column vector of the derrivative of the const function with respect to the fraction of metal for each element.
    Args:
        T (np.ndarray): The global temperature vector during a certain iteration
        v (np.ndarray): Global vector with fractions of metal in each element
        corners(dict): Dictionary with the element idx as key, the node indices in a list as value
        p (float): p > 1 is an float that pushes the fractions of metal to either 0 or 1. Higher values push harder.

    Returns:
        gradJv (np.ndarray): A column vector of size(n) containing the derivative of the cost function with respect to the fraction of metal in each of the n elements.

    """
    elements = T.size
    gradJv = np.zeros((elements,1))
    km = 65
    kp = 0.2
    K0 = np.matrix([[2/3, -1/6, -1/3, -1/6],[-1/6, 2/3, -1/6, -1/3],[-1/3, -1/6, 2/3, -1/6],[-1/6, -1/3, -1/6, 2/3]])

    for el in range(elements):
        elementVertices = corners.get(el)
        Te = np.take(T, elementVertices)
        gradJv[el] = - 0.5 * p * (v[el] ** (p - 1)) * (km - kp) * Te.T @ K0 @ Te

    return gradJv
