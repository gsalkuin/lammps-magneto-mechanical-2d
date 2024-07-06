from lmm2d.geometry import TriMesh
from lmm2d.metamaterial import Composite
from lmm2d.writer import write_lam

import numpy as np
from math import pi, sqrt

msh_file = 'square-achiral-compact-3x3'
tri_mesh = TriMesh(msh_file) # units = mm

# CGS units

# solid properties
thickness = 0.15 # cm
Y_mod = 3e7 # 0.1 Pa (PLA, 3 MPa)
rho_e = 1.25 # g/cm^3
rho_m = 7.5 # g/cm^3

# resize and recenter
L_scale = 0.1 # mesh units to cm
tri_mesh.V_into_xyz = L_scale * tri_mesh.V_into_xyz
tri_mesh.V_into_xyz = tri_mesh.V_into_xyz - np.mean(tri_mesh.V_into_xyz, axis=0)

##############################################################################
# Magnetize

# magnetization
M_r = 9.45e5 # A/m 

# CGS to SI units
c = 3e8 # speed of light (m/s)
eps = 1e-7 # J 
sig = 1e-2 # m 
eps0 = 8.8541878128e-12 # F/m
mu0 = 1/(eps0 * c**2)

M_scale =  c * sqrt(4*pi*eps0 * eps/sig**3) # A/m = 1 emu/cm^3
B_scale =  1/c * sqrt( 1/(4*pi*eps0) * eps/sig**3) # T = 1 Gauss

m_vol = 0.3 * 0.2 * 0.1 # cm^3
m_dip = M_r/M_scale * m_vol

# order depends on the mesh cell blocks
m_dip_per_block = 16*[m_dip*np.array([1, 0, 0])] + [None]
for i in (1,3,4,6,9,11,12,14):
    m_dip_per_block[i] = m_dip*np.array([-1, 0, 0])

##############################################################################
# Create the metamaterial

Y_mod_per_block = 16*[None] + [Y_mod]
rho_per_block = 16*[rho_m] + [rho_e]

meta = Composite(tri_mesh, thickness, rho_per_block, Y_mod_per_block, m_dip_per_block)

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
