from lmm2d.geometry import TriMesh
from lmm2d.metamaterial import Composite
from lmm2d.writer import write_lam


import numpy as np
from math import pi, sqrt, cos, sin

msh_file = 'meta-chain-w13-r17-lc02' 
tri_mesh = TriMesh(msh_file) # units = 1 mm

# CGS units

thickness = 0.318 # cm approx. 1/8"
Y_mod = 1.68e7 # 0.1 Pa (polyurethane PU 40A, 1.68 MPa)
rho_e = 1.25 # 1.05-1.26 g/cm^3 (polyurethane PU 40A)
rho_m = 7.5 

rho_factor = 1e4 # increase density to decrease time step

# resize and recenter
L_scale = 0.1 # mesh units to cm
tri_mesh.V_into_xyz = L_scale * tri_mesh.V_into_xyz
tri_mesh.V_into_xyz = tri_mesh.V_into_xyz - np.mean(tri_mesh.V_into_xyz, axis=0)

##############################################################################
# Magnetize

# CGS to SI units
c = 3e8 # speed of light (m/s)
eps = 1e-7 # J 
sig = 1e-2 # m 
eps0 = 8.8541878128e-12 # F/m
mu0 = 1/(eps0 * c**2)

M_scale =  c * sqrt(4*pi*eps0 * eps/sig**3) # A/m = 1 emu/cm^3
B_scale =  1/c * sqrt( 1/(4*pi*eps0) * eps/sig**3) # T = 1 Gauss

# magnetization
B_r = 1.25 # T (McMaster-Carr Neodymium round magnets 5862K202)
M_r = B_r/mu0 / M_scale

diam_m = 0.31 # cm
m_dip = M_r * (pi * diam_m**2 * thickness / 4)

"""
Specify the 4 dipole moments per unit cell for the order parameter Q.
"""
# NOTE: Q is not unique

# # Q = -1
# filename = msh_file + '-Qneg1'
# m_dip_per_block = 5*[None]
# m_dip_per_block[1] = m_dip * np.array([0, 1, 0])
# m_dip_per_block[2] = m_dip * np.array([0, -1, 0])
# m_dip_per_block[3] = m_dip * np.array([0, -1, 0])
# m_dip_per_block[4] = m_dip * np.array([0, 1, 0])

# # Q = -0.5
# filename = msh_file + '-Qneg5'
# m_dip_per_block = 5*[None]
# theta = 3*pi/4
# for k in range(4):
#     m_dip_per_block[k+1] = m_dip * np.array([cos(theta + k*pi/2), sin(theta + k*pi/2), 0])

# # Q = 0
# filename = msh_file + '-Q0'
# m_dip_per_block = 5*[None]
# for k in range(4):
#     m_dip_per_block[k+1] = m_dip * np.array([0, 0, 0])

# # Q = 0.5
# filename = msh_file + '-Q5'
# m_dip_per_block = 5*[None]
# theta = 3*pi/4
# for k in range(2):
#     m = m_dip * np.array([cos(theta - k*pi/2), sin(theta - k*pi/2), 0])
#     m_dip_per_block[k+1] = m
#     m_dip_per_block[k+3] = m

# # Alternative Q = 0.5 configuration with stronger magnetic interaction
# filename = msh_file + '-Q5-strong'
# m_dip_per_block = 5*[None]
# theta = pi/4
# for k in range(2):
#     m = m_dip * np.array([cos(theta + k*pi/2), sin(theta + k*pi/2), 0])
#     m_dip_per_block[k+1] = m
#     m_dip_per_block[k+3] = m

# Q = 1
filename = msh_file + '-Q1'
m_dip_per_block = 5*[None]
for k in range(2):
    m = m_dip * np.array([0, 1, 0])
    m_dip_per_block[k+1] = m
    m_dip_per_block[k+3] = m

##############################################################################
# Create the metamaterial

Y_mod_per_block = [Y_mod, None, None, None, None]

"""
We want to pre-stretch the elastomer since the magnet is larger than the hole. This will give an initial configuration closer to the experiments.
To do this, we replace the magnet atoms with a large spherical dipole. To avoid double-counting the mass, we set the density of the part to be very small,
but not zero as LAMMPS would complain. Type 3 atoms can be deleted or excluded from the time-integration. Note that the moment of inertia will be incorrect, 
but this is not important for quasi-static simulations. If we wish to keep the mass on atoms, we still need to redistribute the atoms to get the correct MOI. 
Type 4 atoms are connected to the elastomer, so they should not be deleted. The pre-stretching is done by a repulsive 'soft' potential between the dipole 
(type 5) and type 4 atoms. After this, their positions are locked by 'fix rigid'.

NOTE: The magnet and pins are treated as rigid parts bonded to the elastomer, so this may not be realistic especially at large strains.
"""
rho_per_block = [rho_e *rho_factor, 1e-3, 1e-3, 1e-3, 1e-3] 

# the last 4 elements in the list are the end segments: [lower, pin, upper, pin]
meta = Composite(tri_mesh, thickness, 10*rho_per_block + [rho_e *rho_factor, 1., rho_e *rho_factor, 1.], 
                 10*Y_mod_per_block + [Y_mod, None, Y_mod, None], 10*m_dip_per_block + [None, None, None, None])

##############################################################################
# Define elastic boundary atoms	

E_phys = meta.mesh.get_cells_type("line")	
E_phys_ID = meta.mesh.get_cell_data("gmsh:physical", "line")	
boundary_E = np.where(E_phys_ID == 7) # Check the Gmsh file for the physical line ID
boundary_V = np.unique(E_phys[boundary_E])	
meta.V_into_atom_type[boundary_V] = np.where(meta.V_into_atom_type[boundary_V] == 1, 2, meta.V_into_atom_type[boundary_V])	

##############################################################################
# Modify diameter and density; recalculate mass
"""
Correct sequence:
1. set_diam_bVs : convert boundary atoms to spheres for contact model
2. compute_rho_from_mesh : update point masses to densities for boundary atoms
3. set_diam_per_type : overwrite diameter of dipoles (default is shortest bond length)
4. set_rho_per_type : overwrite density of dipoles
"""
meta.set_diam_bVs()
meta.compute_rho_from_mesh()
meta.set_diam_per_type(5, diam_m)
meta.set_rho_per_type(5, rho_m * (pi * diam_m**2 * thickness / 4)/(pi * diam_m**3 / 6) *rho_factor)

##############################################################################	
# Write LAMMPS data file
write_lam([meta], filename=filename)

##############################################################################	
# Visualize
import matplotlib.pyplot as plt

x, y = meta.V_into_xyz[:,0], meta.V_into_xyz[:,1]
plt.scatter(x, y)

# Boundary atoms
xb, yb = meta.V_into_xyz[boundary_V,0], meta.V_into_xyz[boundary_V,1]
plt.scatter(xb, yb, c='tab:orange')

# Dipoles
m_id = np.where(meta.V_into_atom_type == 5)[0]
m_x = meta.V_into_q_mu[m_id, 1]
m_y = meta.V_into_q_mu[m_id, 2]

plt.quiver(x[m_id], y[m_id], m_x, m_y)

plt.axis('equal')
plt.show()