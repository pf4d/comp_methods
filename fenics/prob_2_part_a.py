from fenics import *

mesh = IntervalMesh(3,0,1)
Q    = FunctionSpace(mesh, "CG", 1)
k    = Constant(1.0)
f    = Constant(1.0)
g    = Constant(1.0)
u_0  = Constant(0.0)
u    = TrialFunction(Q)
v    = TestFunction(Q)

a    = k * u.dx(0) * v.dx(0) * dx
l    = f * v * dx + g * v * ds

def left(x, on_boundary):  return x[0] == 0 and on_boundary
bc   = DirichletBC(Q, u_0, left)

u    = Function(Q)
solve(a == l, u, bc)

print(u.vector().get_local()[::-1])
