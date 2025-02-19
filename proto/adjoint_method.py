import numpy as np

n = 9 # total nodes
e = 4 # total elements
K = np.random.rand(n,n)
F = np.random.rand(n,1)

km = 65
kp = 0.2
p = 1
K0 = np.matrix([[2/3, -1/6, -1/3, -1/6],[-1/6, 2/3, -1/6, -1/3],[-1/3, -1/6, 2/3, -1/6],[-1/6, -1/3, -1/6, 2/3]])
v = np.zeros((9,1))

cornerDict = {
    1: [1,2,4,5],
    2: [2,3,5,6],
    3: [4,5,7,8],
    4: [5,6,8,9],
}

def adjoint(K,F,v):
    gradJv = np.zeros((e,1))

    T = np.linalg.solve(K,F)

    for i in range(e):
        idx = i + 1
        corners = cornerDict.get(idx)
        Te = np.array([corners]).T
        gradJv[i] = - 0.5 * p * v[i] ** (p - 1) * (km - kp) * Te.T @ K0 @ Te

    return gradJv
    
grad = adjoint(K,F,v)
print(grad)

