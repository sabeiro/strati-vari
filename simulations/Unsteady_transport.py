
# This script implements FEM solution for a transient advective problem with forcing term equal to zero.
# The transport field is constant in time and divergence-free.
# In this example the inflow boundary is empty, allowing therefore to have empty Dirichlet boundary conditions.
# 
# The exact solution is a gaussian (therefore the solution is smooth) transported around the computational domain given 
# by a 2D full disk. The initial condition is this gaussian computed for t = 0 
#
# The script solves the unsteady pure advective equation by using spatially a SUPG stabilization method. 
# The time discretization can be chosen between a Backward Euler finite difference scheme and a Crank-Nicolson finite 
# difference scheme. 
# 
# We have included a convergence analysis in the L2 norm, the infinity norm and the STABILIZED L2 norm, in order to check 
# the implementation. 
# 
# REFERENCE: [1] Burman "Consistent SUPG-method for transient transport problems:Stability and convergence", 
#                Comput. Methods Appl. Mech. Engrg. 199, 2010 
#
# Created by: Alessandro Fermi - 1st March 2016
# Modified: 8th March 2016    

from dolfin import *
import numpy

# Load mesh with nx subintervals
mesh = CircleMesh(Point(0,0), 1.0, 0.08)
h = CellSize(mesh)

print mesh.hmax()
#mesh.init(1)
#print mesh.num_faces()

# Create finite elements FunctionSpaces of specific degree
degree = 2
Q = FunctionSpace(mesh, "Lagrange", degree)
#V = Q * Q
V = VectorFunctionSpace(mesh, "Lagrange", degree, dim=2)

#define transport field and assign to the nodes the chosen value
transport = Function(V)
V_assigners = [FunctionAssigner(V.sub(i), Q) for i in range(V.num_sub_spaces())]
tra_x_aux = Expression('- x[1]')
tra_x = interpolate(tra_x_aux, Q)
tra_y_aux = Expression('x[0]')
tra_y = interpolate(tra_y_aux,Q)
comps = [tra_x, tra_y]
[V_assigners[i].assign(transport.sub(i), comp) for i,comp in enumerate(comps)]

#define inflow boundary and Dirichlet boundary conditions - in this case the inflow boundary is empty
def inflowboundary(x, on_boundary):
	return on_boundary

u_boundary = Constant(0.0)
bc = DirichletBC(Q, u_boundary, inflowboundary)

# Test and trial functions
u = TrialFunction(Q)
v = TestFunction(Q)

#parameters - theta = 1 for Backward Euler and theta = 0.5 for Crank-Nicolson
theta = 0.5
T = 7.0
dt = 0.05
t = 0 
r1 = 0.42426

#forcing term 
f = Constant(0.0)

#exact solution and initial condition - smooth initial data interpolated on the finite element space Q
u0 = Expression('exp(-(10*pow(x[0] - r1*cos(t + pi/4),2) + 10*pow(x[1] - r1*sin(t + pi/4),2)))', r1=r1, t=0)
u1 = interpolate(u0, Q)

#Linear combination of discrete solutions for Backward Euler and Crank-Nicolson time discretization
#u_mid = 0.5*(u1 + u)
u_mid = theta*u + (1.0 - theta)*u1

# Residual
r = u-u1 + dt*(inner(transport, grad(u_mid)) - f)

# Galerkin variational problem
F = v*(u - u1)*dx + dt*(v*inner(transport, nabla_grad(u_mid))*dx - f*v*dx) 

# Add SUPG stabilisation terms
vnorm = sqrt(dot(transport, transport))
delta = h/(2.0*vnorm)  #stabilization parameter
aux = inner(transport, nabla_grad(v))
F += delta*aux*r*dx
#F += h/(2.0*vnorm)*dot(transport, grad(v))*r*dx

# Create bilinear and linear forms
a = lhs(F)
L = rhs(F) 

# Assemble matrix - the stiffness matrix is time-independent and therefore computed once and for all time levels
A = assemble(a)

# Define functions for the solution and for the mean transport field
u_fin = Function(Q)
W = Function(V)

#set initial condition
u_fin = u1

while (t <= T):
	#assemble right hand side for each time level and solve the algebraic system	
	b = assemble(L)
	u0.t = t
	bc.apply(A, b)
	solve(A, u_fin.vector(), b)

	#If u_fin is a streamfunction, compute the underlying velocity field
#	grad_u_fin = project(grad(u_fin), V)
#	(dx_u_fin, dy_u_fin) = grad_u_fin.split(deepcopy=True)
#	dx_u_fin = project(dx_u_fin, Q)
#	dy_u_fin = project(dy_u_fin, Q)
#	dx_u_fin_array = dx_u_fin.vector().array()
#	dy_u_fin_array = dy_u_fin.vector().array()
#	dx_u_fin.vector()[:] = - dx_u_fin_array
#	comps = [dy_u_fin, dx_u_fin]
#	[V_assigners[i].assign(W.sub(i), comp) for i,comp in enumerate(comps)]
	
	#File output
	#file_StreamFunction = File("StreamFunction.xml")
	#file_MeanTransport = File("MeanTransport.xml")
	#file_StreamFunction << (u_fin, t)
	#file_MeanTransport << (W, t)
	
	#increment the time level
	t +=dt

	#copy solution from previous interval
	u1.assign(u_fin)

	#convergence analysis and error estimation in the L2 norm and in the infinity norm
	Q_exact = FunctionSpace(mesh, 'Lagrange', degree = degree+4)
	u_ex = interpolate(u0, Q_exact)
	u_sol = interpolate(u_fin, Q_exact)
	max_err = abs(u_sol.vector().array() - u_ex.vector().array()).max()
	e_Qe = Function(Q_exact)
	e_Qe.vector()[:] = u_ex.vector().array() - u_sol.vector().array()
	error = e_Qe**2*dx
	E = sqrt(assemble(error))
	print '\nError estimate in the L2 norm:', E
	print 'Error estimate in the infinity norm (at nodal values): ', max_err
	
	#convergence analysis following theory developed in [1]
	u_exact = interpolate(u0, Q)	
	grad_u0 = project(nabla_grad(u_exact), V)
	grad_usol = project(grad(u_fin), V)
	error1 = (u_fin - u_exact)**2*dx
	E1 = assemble(error1)
	error2 = delta**2*inner(transport, grad_usol - grad_u0)**2*dx
	E2 = assemble(error2)
	E_tot = sqrt(E1 + E2)
	print '\nError estimate in the STABILIZED L2 norm:', E_tot
	
	#plot results
	plot(mesh)
	plot(transport, title = 'Transport field')
#	plot(u_ex, title='Exact solution')
	plot(u_fin, title= 'Solution of unsteady transport equation')
#	plot(W , title='Divergence-free velocity field described by the solution')

	# Hold plot
#	interactive()
