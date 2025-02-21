import sys
sys.path.insert(1, "../proto/")
from adjointMethod import adjoint
import numpy as np

T = 273 + np.random.rand(9,1) * 20
corners = [
    [3,4,0,1],
    [4,5,1,2],
    [6,7,3,4],
    [7,8,4,5],
]
v = np.random.rand(4,1)

gradJv = adjoint(T, v, corners)

assert(gradJv.size == v.size)
assert(gradJv.size == len(corners))

for i in np.nditer(gradJv):
    assert not (np.isinf(i))
    assert not (i is None)