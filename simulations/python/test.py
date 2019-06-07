from comp_methods import *
import numpy as np


x    = np.zeros(4)
b    = np.array([0, 1/3, 1/3, 1/6])
vec  = np.array([1,-3,6,-3,-3,6,-3,-3,3])
rows = np.array([0,1,1,1,2,2,2,3,3], dtype="uint16")
cols = np.array([0,0,1,2,1,2,3,2,3], dtype="uint16")
A    = EvanMatrix(vec, rows, cols, 4)

A.print()
diag = A.get_diag()

print(diag)

method = "jacobi"

linear_solve(A, x, b, method, False, 1000, 1e-16, 1e-7);



