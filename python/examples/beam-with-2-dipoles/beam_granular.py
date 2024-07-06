from lmm2d.geometry import TriMesh
from lmm2d.metamaterial import Composite
from lmm2d.writer import write_lam

import numpy as np

msh_file = 'beam-w2dip'
tri_mesh = TriMesh(msh_file) # units = 1 cm

# CGS units

thickness = 1. # cm
Y_mod = 1e7 # 0.1 Pa --> 1e7 = 1MPa
rho_e = 1.25 # g/cm^3
rho_m = 7.5 # g/cm^3

# recenter
tri_mesh.V_into_xyz = tri_mesh.V_into_xyz - np.mean(tri_mesh.V_into_xyz, axis=0)

##############################################################################
# NOTE: Magnetization is set in the LAMMPS input file.

m_dip_per_block = [None, np.array([0., 0., 0.]), np.array([0., 0., 0.]), None, None]

##############################################################################
# Create the metamaterial

Y_mod_per_block = [Y_mod, None, None, None, None]
rho_per_block = [rho_e, rho_m, rho_m, rho_m, rho_m] 

beam = Composite(tri_mesh, thickness, rho_per_block, Y_mod_per_block, m_dip_per_block)

##############################################################################
# Define boundary atoms; modify diameter and density; recalculate mass

# Here, we use type 4 atoms for contact
beam.set_diam_bVs(atom_type=4)
beam.compute_rho_from_mesh()

##############################################################################
# Set rigid ID since there are two adjacent rigid (non-bonded) molecules
beam.set_rigid_ID([1], 0) # non-rigid is always 0
beam.set_rigid_ID([2,4], 1)
beam.set_rigid_ID([3,5], 2)

##############################################################################	
# Write LAMMPS data file
filename = msh_file + '-granular-' + str(Y_mod/1e7) + 'MPa'
write_lam([beam], filename=filename)

##############################################################################	
# Visualize
import matplotlib.pyplot as plt

x, y = beam.V_into_xyz[:,0], beam.V_into_xyz[:,1]
plt.scatter(x, y)

# Boundary atoms
xb, yb = beam.V_into_xyz[beam.V_into_atom_type==4,0], beam.V_into_xyz[beam.V_into_atom_type==4,1]
plt.scatter(xb, yb, c='tab:orange')

# Dipoles
m_id = np.where(beam.V_into_atom_type == 5)[0]
m_x = beam.V_into_q_mu[m_id, 1]
m_y = beam.V_into_q_mu[m_id, 2]

plt.quiver(x[m_id], y[m_id], m_x, m_y)

plt.axis('equal')
plt.show()