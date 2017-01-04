
# This script implements FEM solution for a pure steady advective problem with non-zero forcing term and Dirichlet 
# boundary conditions: 
#
#                          transport grad u = f    
#                          u = u_0      on the inflow boundary 
#
# The forcing term is computed directly by knowing the exact solution and the transport field.
# The transport field is a 2D divergence-free vector field with corresponding inflow boundary easily identified.   
#
# The script solves the steady pure advective equation by a SUPG stabilization method. 
# A convergence analysis and error estimation in the L2 norm and in the infinity norm is carried out in order to check the # implementation. 
#
# Created by: Alessandro Fermi - 2 March 2016
#        

from dolfin import *
import numpy
 
# Load mesh with nx subintervals
nx = 120
mesh = RectangleMesh(0.0, 0.0, 10.0, 10.0, nx, nx, 'left')  #mesh on a 2D rectangle 
h = CellSize(mesh)
h1 = 10.0/nx
print h1

# Create Lagrange finite elements Function Spaces of specific degree
degree = 1
Q = FunctionSpace(mesh, "Lagrange", degree)
#V = Q * Q
V = VectorFunctionSpace(mesh, "Lagrange", degree, dim=2)

#define the transport field as an element in the vector function space V
transport = Function(V)
V_assigners = [FunctionAssigner(V.sub(i), Q) for i in range(V.num_sub_spaces())]   #specify the subspaces of V
tra_x_aux = Expression('x[1]')
tra_x = interpolate(tra_x_aux, Q)
tra_y_aux = Expression('x[0]')
tra_y = interpolate(tra_y_aux,Q)
comps = [tra_x, tra_y]
[V_assigners[i].assign(transport.sub(i), comp) for i,comp in enumerate(comps)]

#define the inflow boundary and Dirichlet boundary conditions  
def inflowboundary(x):
	if (near(x[0], 0.0) | near(x[1], 0.0)):
		return True
	else:
		return False

u_0 = Expression('(1/pi)*sin(pi*x[0])*cos(pi*x[1])')   
bc = DirichletBC(Q, u_0, inflowboundary)

# Test and trial functions
u = TrialFunction(Q)
v = TestFunction(Q)

# Define the forcing term
f = Expression('x[1]*cos(pi*x[0])*cos(pi*x[1]) - x[0]*sin(pi*x[0])*sin(pi*x[1])')

#define exact solution and corresponding finite element approximation 
u_exact = Expression('(1/pi)*sin(pi*x[0])*cos(pi*x[1])')
u_exact = interpolate(u_exact, Q)
#u_exact = project(u_exact, Q)

# Galerkin variational problem
F = v*inner(transport, nabla_grad(u))*dx - f*v*dx 

# Residual for the SUPG stabilization term
r = inner(transport, nabla_grad(u)) - f

# Add SUPG stabilisation terms
vnorm = sqrt(dot(transport, transport))
delta = h/(2.0*vnorm)
aux = inner(transport, nabla_grad(v))
F += delta*aux*r*dx
#F += h/(2.0*vnorm)*dot(transport, grad(v))*r*dx

# Create bilinear and linear forms
a = lhs(F)
L = rhs(F) 

# Define functions for the solution and for the mean transport field
u_fin = Function(Q)
W = Function(V)
	
# Assemble matrix and right hand side, apply boundary conditions and solve
A = assemble(a)
b = assemble(L)
#print len(b.array())
#print b.array()
#print A.array().shape
#print A.array()
	
bc.apply(A, b)
U_fin = u_fin.vector()
solve(A, U_fin, b)
#problem = LinearVariationalProblem(a, L, u_fin, bc)
#solver = LinearVariationalSolver(problem)
#solver.solve()

#If u_fin is a streamfunction, compute derivatives of u_fin and the component of the underlying velocity field
grad_u_fin = project(grad(u_fin), V)
(dx_u_fin, dy_u_fin) = grad_u_fin.split(deepcopy=True)
dx_u_fin = project(dx_u_fin, Q)
dy_u_fin = project(dy_u_fin, Q)
dx_u_fin_array = dx_u_fin.vector().array()
dy_u_fin_array = dy_u_fin.vector().array()
dx_u_fin.vector()[:] = - dx_u_fin_array
comps = [dy_u_fin, dx_u_fin]
[V_assigners[i].assign(W.sub(i), comp) for i,comp in enumerate(comps)]
	
#File output  
#file_StreamFunction = File("StreamFunction.xml")
#file_MeanTransport = File("MeanTransport.xml")
#file_StreamFunction << u_fin
#file_MeanTransport << W
	
#convergence analysis and error estimation in the L2 norm and in the infinity norm
Q_exact = FunctionSpace(mesh, 'Lagrange', degree = degree+4)
u_ex = interpolate(u_exact, Q_exact)
u_sol = interpolate(u_fin, Q_exact)
max_err = abs(u_sol.vector().array() - u_ex.vector().array()).max()
e_Qe = Function(Q_exact)
e_Qe.vector()[:] = u_ex.vector().array() - u_sol.vector().array()
error = e_Qe**2*dx
E = sqrt(assemble(error))
print '\nError estimate in the L2 norm:', E
print '\nError estimate in the infinity norm (at nodal values): ', max_err
	
#plot the results
plot(mesh)
plot(transport, title = 'Transport field')
plot(u_ex, title='Exact solution')
plot(u_fin, title= 'First solution ever')
#plot(W , title='First mean transport wind with FEM ever')

# Hold plot
interactive()
