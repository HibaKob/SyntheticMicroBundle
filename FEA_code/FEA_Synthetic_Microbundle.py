# =============================================================================
# NeoHookean Material Model:
#   Nearly Incompressible Transversely Isotropic Active Contraction: 
#        Active Strain Formulation
# =============================================================================

from dolfin import *
from fenics import *
import numpy as np
import os
import matplotlib.pyplot as plt

#################### User-defined functions ####################

# Generate time series activation
def timeseries_act(var_seed,frequency,frames,amplitude):
	rndm = np.random.default_rng(seed=var_seed)
	x = np.linspace(0,2*np.pi,frames)
	y = -amplitude*np.cos((frequency)*x)
	y[y < 0.008] = 0
	n_noise_pts = len(y)
	noise = rndm.normal(0,0.001,n_noise_pts)
	# Correlation coefficient of Random Walk noise
	r = 0.95
	noise_r = np.zeros(n_noise_pts)
	noise_r[0] = noise[0]
	for ii in range(n_noise_pts-1):
		noise_r[ii+1] = r*noise_r[ii] + (1 - r**2)**(1/2)*noise[ii+1]
	y_n = y + noise_r 
	return y_n

# Generate spatially heterogeneous activation: Active and passive regions
def hetero_activation(x,y,incl_centers, incl_min_rad):
    d_min = incl_min_rad # radius of passive inclusion
    d_min_ss = d_min**2
    alpha = -80 # parameter of smooth min function
    beta = 0.9
    num = 0 # numerator of smooth min function
    denum = 0 # denominator of smooth min function
    for point in incl_centers:
        dist = (point[0]-x)**2 + (point[1]-y)**2
        num += dist*exp(alpha*dist)
        denum += exp(alpha*dist)
    d = num/denum 
    return conditional(lt(d, d_min), 1, d_min_ss/(beta*d_min_ss + (1-beta)*d**2))

# Get vector variables at slice
def get_disp(v, lst_x_y, z_slice):
	d_x = []
	d_y = []
	d_z = []
	for ii in range(len(lst_x_y)):
		x_coord = lst_x_y[ii][0]
		y_coord = lst_x_y[ii][1]
		v_x = v([x_coord,y_coord,z_slice])[0]
		d_x.append(v_x)
		v_y = v([x_coord,y_coord,z_slice])[1]
		d_y.append(v_y)
		v_z = v([x_coord,y_coord,z_slice])[2]
		d_z.append(v_z)
	return d_x, d_y, d_z

# Get tensor variables at slice (output all strains)
def get_strain(t, lst_x_y, z_slice):
	S_xx = []
	S_xy = []
	S_xz = []
	S_yy = []
	S_yz = []
	S_zz = []
	for ii in range(len(lst_x_y)):
		x_coord = lst_x_y[ii][0]
		y_coord = lst_x_y[ii][1]
		t_xx = t([x_coord,y_coord,z_slice])[0]
		t_xy = t([x_coord,y_coord,z_slice])[1]
		t_xz = t([x_coord,y_coord,z_slice])[2]
		t_yy = t([x_coord,y_coord,z_slice])[4]
		t_yz = t([x_coord,y_coord,z_slice])[5]
		t_zz = t([x_coord,y_coord,z_slice])[8]
		S_xx.append(t_xx) 
		S_xy.append(t_xy)
		S_xz.append(t_xz)
		S_yy.append(t_yy)
		S_yz.append(t_yz) 
		S_zz.append(t_zz)		
	return S_xx, S_xy, S_xz, S_yy, S_yz, S_zz

# Get tissue coordinates
def get_Coord(subdomain, subdomain_ID):
	tissue_cells = SubsetIterator(subdomain, subdomain_ID)
	midp = []
	for cell in tissue_cells:
		tiss_mid= cell.midpoint()
		if tiss_mid[2] > -0.02 - DOLFIN_EPS:
			mid_xy = np.array([tiss_mid[0],tiss_mid[1]])
			midp.append(mid_xy)
	return midp
################################################################

# Define mesh
# XDMF Method
# Import mesh and define function space
mesh_file = XDMFFile(MPI.comm_world, "MicroTug_RefDim_3D.xdmf")
mesh = Mesh()
mesh_file.read(mesh);

mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
mesh_file.read(mvc, "regions")
sub = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Change domain IDs from [1 1 1 ... 2 2 2] to [0 0 0 .... 1 1 1]
# get element indices for a given domain id:
domain1_ids = sub.where_equal(1)
domain2_ids = sub.where_equal(2)

# first create a marker and flush it with zeros everywhere:
marker = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
# create an integer array of the size equal to number of elements in mesh
indarray = np.zeros(np.size(marker.array()),'int64')
indarray[domain1_ids] = 0
indarray[domain2_ids] = 1
marker.set_values(indarray)

dx = Measure('dx', domain=mesh, subdomain_data=marker)

#File('subdomains3D_RefDim_MicroTug.pvd') << marker

flag_quad = 2
# Compiler settings / optimization options 
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = flag_quad

ffc_options = {"optimize": True, \
                "eliminate_zeros": True, \
                "precompute_basis_const": True, \
                "precompute_ip_const": True}

# MPI communicator
comm = MPI.comm_world
rank = comm.Get_rank()
size = comm.Get_size()

parameters['allow_extrapolation'] = True

# Define function space to solve the problem (Taylor-Hood elements)
P2 = VectorElement("Lagrange", mesh.ufl_cell(), flag_quad)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

S = FunctionSpace(mesh,'CG', 1)
V = VectorFunctionSpace(mesh,'CG', 2)
T1= TensorFunctionSpace(mesh,'DG', 1) 

V0 = FunctionSpace(mesh,'DG',0)     

# Define traction on the boundary and body forces
T  = Constant((0.0, 0.0, 0.0)) 
B  = Constant((0.0, 0.0, 0.0))

# Define finite element problem
u_ = Function(W)                      
du_ = TrialFunction(W)
v_ = TestFunction(W)                  

# Displacement and hydrostatic_pressure
(u, p) = split(u_) 
(v, q) = split(v_)

# Mark boundaries
left_pillar = CompiledSubDomain('on_boundary && x[0]>0.205+tol && x[0]<0.595+tol && x[1]>0.185+tol && x[1]<0.615+tol && near(x[2],-0.995, tol)',tol=1E-5) 
right_pillar = CompiledSubDomain('on_boundary && x[0]>1.405+tol && x[0]<1.795+tol && x[1]>0.185+tol && x[1]<0.615+tol && near(x[2],-0.995, tol)',tol=1E-5)
 
# Boundary conditions 
lftBC = DirichletBC(W.sub(0), Constant((0.0,0.0,0.0)), left_pillar) 
rgtBC = DirichletBC(W.sub(0), Constant((0.0,0.0,0.0)), right_pillar)  
bcs = [lftBC,rgtBC]

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all (0)
left_pillar.mark(boundaries ,1)
right_pillar.mark(boundaries ,2)
ds = ds(domain=mesh , subdomain_data=boundaries)

# Kinematics
dim = mesh.geometry().dim()
I = Identity(dim)           # Identity tensor
F = I + grad(u)             # Deformation gradient
F = variable(F)
J = det(F)
F_bar = variable(pow(J,-1/dim)*F)

# Define variable fiber direction with respect to tissue depth
f0 = Expression(("cos((alpha_0*(0.4-abs(x[2])) + alpha_1*abs(x[2]))/0.4)", "sin((alpha_0*(0.4-abs(x[2])) + alpha_1*abs(x[2]))/0.4)", "0.0"), alpha_0 = 9.33*pi/180, alpha_1 = 15.33*pi/180, degree=1, cell=mesh.ufl_cell()) # Tissue depth = 0.4
f0f0 = outer(f0, f0)

# Define activation as an updatable expression
activation = Expression("(act_value)", act_value = 0.0, degree = 0)   
mgamma = 1 - activation # Comment this line when implementing heterogeneous activation 

''' Define spatially heterogeneous activation (uncomment below code to implement this option) '''
# XX = SpatialCoordinate(mesh)
# centers = [[1,0.4]]
# min_r = 0.1
# activation_hetero = activation - activation*hetero_activation(XX[0], XX[1],centers,min_r)
#mgamma = 1 - activation_hetero

# Deformation gradient: Multiplicative split F = Fe*Fa 
Fa = mgamma * f0f0 + pow(mgamma, -1.0 / float(dim - 1)) * (I - f0f0) 
Fe = variable(F * inv(Fa))
Fe_bar= variable(F_bar * inv(Fa))

# Right Cauchy-Green tensor
C = variable(F.T*F)                 
C_bar = variable(F_bar.T*F_bar) 
Ce = variable(Fe.T*Fe) 
Ce_bar = variable(Fe_bar.T*Fe_bar)        

# First invariant of deformation tensor
I_1 = tr(C)
I_1_e = tr(Ce)
I_1_bar = tr(C_bar)
I_1_e_bar = tr(Ce_bar)

# Tissue material parameters
mu = 0.0127 #MPa (dimensions in mm) # Comment this line for adding inclusions 
kappa = 10.0 #MPa

''' Define a stiffer/less stiff circular inclusion (uncomment below code to implement this option)'''
#mu = Expression('((pow((x[0]-inc_centrX),2)+pow((x[1]-inc_centrY),2))<inc_r*inc_r) ? factor*0.0127 : 0.0127', inc_centrX=1,inc_centrY=0.4,inc_r=0.1,factor=5,degree=0, cell=mesh.ufl_cell())
#mu = project(mu,V0)

#Pillars
E = 1.6 #MPa
nu = 0.47
mu_p = E/(2*(1+nu))
lmbda = E*nu/((1+nu)*(1-2*nu))

psi_Pillars = 1/2*mu_p*( inner(F,F) - 3 - 2*ln(det(F)) ) + 1/2*lmbda*(1/2*(det(F)**2 - 1) - ln(det(F)))

psi_vol = (kappa/2)*pow(ln(J),2) # Doesn't make a differnce since not used in solving
psi = 1/2*mu*(I_1_e_bar - dim) - p*(J-1)

# First Piola-Kirchhoff stress
Piola_first = diff(psi,F)
# First Piola-Kirchhoff stress (Pillars)
Piola_first_Pillars = diff(psi_Pillars,F)

stress_work = inner(grad(v),Piola_first)*dx(0) + inner(grad(v),Piola_first_Pillars)*dx(1)
ext_work = dot(B, v)*dx('everywhere') + dot(T, v)*ds 
incomp_const = q*((ln(J)/J) - (1/kappa)*p)*dx

G_all = stress_work - ext_work + incomp_const

dG = derivative(G_all, u_, du_)

# Green-Lagrange strain
EE_s = 0.5*(C - I)

# Second Piola-Kirchhoff stress
Piola_secnd = inv(F)*Piola_first  

# Cauchy stress
sig = F*Piola_secnd*F.T*((1/det(F))*I) #in terms of second Piola

uu = Function(V, name="Displacement")
pp = Function(S, name="Pressure")
jj = Function(S, name="J")
Cauchy = Function(T1, name="Cauchy Stress")
fPiola = Function(T1, name="First Piola")
E_strain = Function(T1, name="Green Lagrange Strain")

Coord_tiss = get_Coord(marker,0)

comm.barrier()
Coord_tiss_all = np.vstack(comm.allgather(Coord_tiss))
comm.barrier()

if rank == 0:
	np.savetxt('Tissue_Slice_Coordinates.txt',np.array(Coord_tiss_all),fmt='%0.5e')

# Define timeseries activation values
max_act = 0.1 
activation_values = timeseries_act(100,8,200,max_act)

file_results = XDMFFile("Neo_Homog_Act{0}_VFA.xdmf".format(max_act))
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

# Create folder to save results
save_folder = "FEA_Results_Homog_MaxAct{0}_VFA/".format(max_act)
if not os.path.exists(save_folder):
    os.makedirs(save_folder, exist_ok=True)

z_coord = 0.0

step = 0
# Stop simulation once all activation values are passed
for act in activation_values:
	# Update activation
	activation.act_value = act
    
	# Solve
	solve(G_all == 0, u_, bcs, J=dG, form_compiler_parameters=ffc_options)
	(u, p) = u_.split(True) # Split the mixed solution using deepcopy 
	
	disp_x, disp_y, disp_z = get_disp(u,Coord_tiss,z_coord)
	E_xx, E_xy, E_xz, E_yy, E_yz, E_zz = get_strain(project(EE_s,T1),Coord_tiss,z_coord)
	
	comm.barrier()
	disp_x_all = np.hstack(comm.allgather(disp_x))
	disp_y_all = np.hstack(comm.allgather(disp_y))
	disp_z_all = np.hstack(comm.allgather(disp_z))
	
	disp_all = np.array([disp_x_all,disp_y_all,disp_z_all])
		
	E_xx_all = np.hstack(comm.allgather(E_xx))
	E_xy_all = np.hstack(comm.allgather(E_xy))	
	E_xz_all = np.hstack(comm.allgather(E_xz))
	E_yy_all = np.hstack(comm.allgather(E_yy))
	E_yz_all = np.hstack(comm.allgather(E_yz))
	E_zz_all = np.hstack(comm.allgather(E_zz))
	
	E_all = np.array([E_xx_all,E_xy_all,E_xz_all,E_yy_all,E_yz_all,E_zz_all])	
	comm.barrier()

	if rank == 0:
		np.savetxt(save_folder + 'disp_all_Step{0}.txt'.format(step),np.array(disp_all),fmt='%0.6e')
		np.savetxt(save_folder + 'Estrain_all_Step{0}.txt'.format(step),np.array(E_all),fmt='%0.6e')

	uu.assign(u)
	pp.assign(p)
	jj.assign(project(J,S))
	fPiola.assign(project(Piola_first, T1))
	Cauchy.assign(project(sig,T1))
	E_strain.assign(project(EE_s,T1)) 
	
	file_results.write(uu, step)
	file_results.write(pp, step)
	file_results.write(jj, step)
	file_results.write(fPiola, step)
	file_results.write(Cauchy, step)
	file_results.write(E_strain, step)

	step +=1
