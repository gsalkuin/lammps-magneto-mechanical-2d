"""
Generate a mesh for the cross-section of the cellular solid in Tipton, C. R., Han, E., & Mullin, T. (2012). Magneto-elastic buckling of a soft cellular solid. Soft Matter, 8(26), 6880-6883.
"""

filename = "cell-lc025"
write_path = filename + ".msh"

from math import pi, cos, sin
import gmsh
import sys

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()
gmsh.model.add(filename)

# Set the kernel to OpenCASCADE
gmsh.option.setNumber("General.Terminal", 1)
# Set the tolerance for combining points
gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
# Set the target mesh size globally
gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
gmsh.option.setNumber("Mesh.Algorithm", 5)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)

# Parameters (in Fig. 2c)
L = 10 # length of the small square
w = 0.5
R = (L-w)/2 # hole radius 
r = 1.25 # magnet radius

#################################################################################################
# Create the elastomer

ela_l = gmsh.model.occ.addRectangle(0, 0, 0, 4*L, 4*L)

# Create the holes
holes = []
for i in range(4):
    for j in range(4):
        h = gmsh.model.occ.addDisk(L/2+i*L, L/2+j*L, 0, R, R)
        holes.append(h)

# Remove the holes from the elastomer
gmsh.model.occ.cut([(2, ela_l)], [(2, hole) for hole in holes], removeObject=True)

# Create smaller region
ela_s = gmsh.model.occ.addRectangle(L/2, L/2, 0, 3*L, 3*L)

# Create the cellular structure
cell = gmsh.model.occ.intersect([(2, ela_l)], [(2, ela_s)], removeObject=True)[0]

gmsh.model.occ.synchronize()

# Add the magnets
mags = []
for i in range(3):
    for j in range(3):
        mag = gmsh.model.occ.addDisk(L+i*L, L+j*L, 0, r, r)
        mags.append(mag)

# Separate the magnets from the cell
gmsh.model.occ.cut(cell, [(2, mag) for mag in mags], removeObject=True, removeTool=False)

#################################################################################################
# Meshing

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)


# #################################################################################################
# # Physical groups

cell_b = [abs(tag) for [dim, tag] in gmsh.model.getBoundary(cell)]

# Boundary curve of the elastic region
gmsh.model.addPhysicalGroup(1, [abs(curve) for curve in cell_b], name = "elastic boundary")

gmsh.model.addPhysicalGroup(2, [cell[0][1]], name = "Elastic")
gmsh.model.addPhysicalGroup(2, mags, name = "Magnetic")

# #################################################################################################
# # Write to file

gmsh.write(write_path)

# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
# command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# We can use this to clear all the model data:
gmsh.clear()

gmsh.finalize()

