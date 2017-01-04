
# This script implements FEM solution for a pure advective problem with non zero forcing term parametrized by time.
# In order to check the implementation a gaussian is chosen and is transported with a constant transport field in the  
# computational domain. 
# The unknown function u is interpreted as the streamfunction of a 2D divergence-free velocity field W, that 
# is computed by taking the proper directional derivative of u.
# 
# The script solves a family of pure advective equations parametrized by time, using a SUPG stabilization method. More 
# precisely, there is no time derivative of the unknown function in the advection model and the transport field is constant 
# in time, but the forcing term and Dirichlet boundary conditions are time-varying. 
# Therefore the exact solution happens to depend on time. 
# In other words we are solving a NON-recursive set of stationary problems uncoupled with each other.
#
# REFERENCE: [1] Burman "Consistent SUPG-method for transient transport problems:Stability and convergence", 
#                Comput. Methods Appl. Mech. Engrg. 199, 2010 
#
# Created by: Alessandro Fermi - 8 March 2016

from dolfin import *
import numpy
 
# Load mesh with nx subintervals
nx = 45
mesh = RectangleMesh(0.0, 0.0, 1.0, 1.0, nx, nx, 'left')  #mesh on a 2D rectangle where the Uji network is located
h = CellSize(mesh)
h1 = 1.0/nx
print h1

# Create finite elements Function Spaces of specific degree
degree = 1
Q = FunctionSpace(mesh, "Lagrange", degree)
#V = Q * Q
V = VectorFunctionSpace(mesh, "Lagrange", degree, dim=2)
transport = Function(V)
V_assigners = [FunctionAssigner(V.sub(i), Q) for i in range(V.num_sub_spaces())]  #component of vector function space V 

# definition of the transport field 
tra_x_aux = Expression('1.0')     
tra_x = interpolate(tra_x_aux, Q)
tra_y_aux = Expression('1.0')		
tra_y = interpolate(tra_y_aux,Q)
comps = [tra_x, tra_y]
[V_assigners[i].assign(transport.sub(i), comp) for i,comp in enumerate(comps)] 

#definition of inflow boundary for the transport field in this example
def inflowboundary(x):
	if near(x[0], 0.0) or near(x[1], 0.0):
		return True
	else:
		return False 

#time parameters 
T = 1.2
dt = 0.02
t = 0.0

#exact solution - gaussian function and corresponding finite element approximation  
u0 = Expression('exp(-(30*pow(x[0] - 0.3 - t,2) + 30*pow(x[1] - 0.3 - t,2)))', t=t)
u_exact = interpolate(u0, Q)  
plot(u_exact, title='Exact solution at t = 0')
interactive()

#Dirichlet boundary conditions
bc = DirichletBC(Q, u0, inflowboundary)

# Test and trial functions
u = TrialFunction(Q)
v = TestFunction(Q)

# forcing term time-dependent
f  = Expression('- exp(-(30*pow(x[0]-0.3-t,2) + 30*pow(x[1]-0.3-t,2)))*60*(x[0]-0.3-t) - exp(-(30*pow(x[0]-0.3-t,2) + 30*pow(x[1]-0.3-t,2)))*60*(x[1]-0.3-t)', t=t)

# Galerkin variational problem
#F = dt*v*inner(transport, nabla_grad(u))*dx - f*v*dx 
F = v*inner(transport, nabla_grad(u))*dx - f*v*dx 

# Residual for the SUPG stabilization term
#r = dt*inner(transport, nabla_grad(u)) - f
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

#assemble stiffness matrix - time-independent because the transport field is constant in time
A = assemble(a)

# Define functions for the solution and for the mean transport field
u_fin = Function(Q)
W = Function(V)

while (t <= T):
	
	#time-dependent exact solution update
	u0.t = t
	u_exact = interpolate(u0, Q)
	f.t = t
	print t	
	
	# Assemble right hand side - this depends on time - and solve the algebraic system
	b = assemble(L)
	bc.apply(A, b)
	solve(A, u_fin.vector(), b)

	#Being u_fin the streamfunction, compute derivatives as the component of the velocity field of the flow
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
	#file_StreamFunction << (u_fin, t)
	#file_MeanTransport << (W, t)
	
	#convergence analysis and error estimation in the L2 norm and in the infinity norm
	Q_exact = FunctionSpace(mesh, 'Lagrange', degree = degree+2)
	u_ex = interpolate(u_exact, Q_exact)
	u_sol = interpolate(u_fin, Q_exact)
	max_err = abs(u_sol.vector().array() - u_ex.vector().array()).max()
	e_Qe = Function(Q_exact)
	e_Qe.vector()[:] = u_ex.vector().array() - u_sol.vector().array()
	error = e_Qe**2*dx
	E = sqrt(assemble(error))
	print '\nError estimate in the L2 norm:', E
	print '\nError estimate in the infinity norm (at nodal values): ', max_err
	
	#convergence analysis in the SUPG stabilized L2norm - see [1]  	
	grad_u0 = project(nabla_grad(u_exact), V)
	grad_usol = project(grad(u_fin), V)
	error1 = (u_fin - u_exact)**2*dx
	E1 = assemble(error1)
	error2 = delta**2*inner(transport, grad_usol - grad_u0)**2*dx
	E2 = assemble(error2)
	E_tot = sqrt(E1 + E2)
	print '\nError estimate in the STABILIZED L2 norm:', E_tot
	
	#plot the results
	plot(mesh, title='Computational mesh')
	plot(transport, title = 'Transport field')
	plot(u_fin, title= 'FEM streamfunction solution', rescale=False)
#	plot(u_fin, title= 'FEM streamfunction solution')
	plot(W , title='Divergence-free velocity field described by the solution', rescale=False)
#	plot(W , title='Divergence-free velocity field described by the solution')

	t +=dt
	# Hold plot
	#interactive()
