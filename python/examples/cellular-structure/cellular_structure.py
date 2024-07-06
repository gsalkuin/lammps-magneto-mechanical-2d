from lmm2d.geometry import TriMesh
from lmm2d.metamaterial import Composite
from lmm2d.writer import write_lam

import numpy as np
from math import pi, sqrt

msh_file = 'cell-lc025'
tri_mesh = TriMesh(msh_file) # units = 1 mm

# CGS units

# solid properties
thickness = 1.5 # cm
Y_mod = 1.4e7 # 0.1 Pa
rho_e = 1.25 # (guess)
rho_m = 7.5 # (guess)

# increase rho to increase time step
rho_factor = 1e6

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

diam_m = 0.25

m_dip = 5 * 0.014 / (M_scale*sig**3) 
m_dip_per_block = [None] + 9 * [m_dip * np.array([0., 0., 1.])]

##############################################################################
# Create the metamaterial

Y_mod_per_block = [Y_mod] + 9*[None]
rho_per_block = [rho_e * rho_factor] + 9*[rho_m * rho_factor] 

meta = Composite(tri_mesh, thickness, rho_per_block, Y_mod_per_block, m_dip_per_block)

##############################################################################
# Define elastic boundary atoms	

E_phys = meta.mesh.get_cells_type("line")	
E_phys_ID = meta.mesh.get_cell_data("gmsh:physical", "line")	
boundary_E = np.where(E_phys_ID == 1)	
boundary_V = np.unique(E_phys[boundary_E])	
meta.V_into_atom_type[boundary_V] = np.where(meta.V_into_atom_type[boundary_V] == 1, 2, meta.V_into_atom_type[boundary_V])	

##############################################################################
# Modify diameter and density; recalculate mass
"""
Correct sequence:
1. set_diam_bVs : convert boundary atoms to spheres for contact model
2. compute_rho_from_mesh : update point masses to densities for boundary atoms
"""
meta.set_diam_bVs()
meta.compute_rho_from_mesh()

##############################################################################	
# Write LAMMPS data file
write_lam([meta], filename=msh_file)

##############################################################################	
# Visualize
import matplotlib.pyplot as plt

x, y = meta.V_into_xyz[:,0], meta.V_into_xyz[:,1]

m_id = np.where(meta.V_into_atom_type == 5)[0]
m_x = meta.V_into_q_mu[m_id, 1]
m_y = meta.V_into_q_mu[m_id, 2]

plt.scatter(x, y)
plt.quiver(x[m_id], y[m_id], m_x, m_y, norm=None)

plt.axis('equal')
plt.show()
