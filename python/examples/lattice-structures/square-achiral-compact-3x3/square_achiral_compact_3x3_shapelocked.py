from lmm2d.geometry import TriMesh
from lmm2d.metamaterial import Composite
from lmm2d.writer import write_lam

import numpy as np
import os, glob, re

msh_file = 'square-achiral-compact-3x3'
tri_mesh = TriMesh(msh_file) # units = mm

# CGS units

# solid properties
thickness = 0.15 # cm
Y_mod = 2.4e10 # 0.1 Pa (PLA at room T, 2.4 GPa)
rho_e = 1.25 # g/cm^3
rho_m = 7.5 # g/cm^3

# modulus is too high, so we increase the density to decrease the time step
rho_factor = 1000

##############################################################################
# Set dipole moment to None since magnet interactions are negligible
# NOTE: This removes type 5
m_dip_per_block = None

##############################################################################
# Create the metamaterial

Y_mod_per_block = 16*[None] + [Y_mod]
rho_per_block = 16*[rho_m *rho_factor] + [rho_e *rho_factor]

# Get the directory of the executed script
script_dir = os.path.dirname(os.path.abspath(__file__))

pattern = re.compile(r'B(\d+)mT')

# Loop over all .xyz files in the curved-xyz directory
for file_path in glob.glob(os.path.join(script_dir, 'curved-xyz/*.xyz')):
    if os.path.isfile(file_path):
      match = pattern.search(file_path)
      if match:
         B = match.group(1) # only used for naming the file

         filename = msh_file + '-B{}mT-curved'.format(B) 
         curved_txyz = np.loadtxt(file_path, skiprows=2)
         tri_mesh.V_into_xyz = curved_txyz[:,1:]

         # recenter
         tri_mesh.V_into_xyz = tri_mesh.V_into_xyz - np.mean(tri_mesh.V_into_xyz, axis=0)

         meta = Composite(tri_mesh, thickness, rho_per_block, Y_mod_per_block, m_dip_per_block)
         write_lam([meta], filename=filename)
