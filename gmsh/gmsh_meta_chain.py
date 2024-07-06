"""
Generate a mesh for the metamaterial chain in Liang, X., Fu, H., & Crosby, A. J. (2022). 
Phase-transforming metamaterial with magnetic interactions. Proceedings of the National Academy of Sciences, 119(1), e2118161119.
"""

filename = "meta-chain-w13-r17-lc02"
write_path = filename + ".msh"

from math import pi
import gmsh
import sys

gmsh.initialize()
gmsh.model.add(filename)

# Set the kernel to OpenCASCADE
gmsh.option.setNumber("General.Terminal", 1)
# Set the tolerance for combining points
gmsh.option.setNumber("Geometry.Tolerance", 1e-6)

# Set the target mesh size globally
lc = 0.2
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

# Parameters in Fig. 1c (mm units)
L = 6 # length of the small square
# w = L/6
w = 0.13*L
r = 0.17 # aspect ratio
b = (L-w)/(1+r)
a = b*r
R = 1.45 # hole radius < magnet radius 1.55

#################################################################################################
# Create the unit cell (square with 4 magnets)
"""
Tags: dim=1
1: center ellpise hole
2: half-ellipse (long axis)
    3: right
    4: left
5: half-ellipse (short axis)
    6: up
    7: down
quarter-ellipse (4 corners)
    8: upper right
    9: upper left
    10: lower left
    11: lower right
"""

# central ellipse hole
ctr_ell = gmsh.model.occ.addEllipse(0, 0, 0, b, a)
ctr_ell_c = gmsh.model.occ.addCurveLoop([ctr_ell])

# right and left half-ellipses (along long axis), for some reason addEllipse can't have R_x < R_y.
half_ell_long = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=2, angle2=pi)

half_ell_r = gmsh.model.occ.copy([(1, 2)])
gmsh.model.occ.rotate(half_ell_r, 0, 0, 0, 0, 0, 1, pi/2)
gmsh.model.occ.translate(half_ell_r, L, 0, 0)

half_ell_l = gmsh.model.occ.copy([(1, 2)])
gmsh.model.occ.rotate(half_ell_l, 0, 0, 0, 0, 0, 1, -pi/2)
gmsh.model.occ.translate(half_ell_l, -L, 0, 0)

# right and left half-ellipses (along short axis)
half_ell_short = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=5, angle1=pi/2, angle2=3*pi/2)

half_ell_u = gmsh.model.occ.copy([(1, 5)])
gmsh.model.occ.rotate(half_ell_u, 0, 0, 0, 0, 0, 1, pi/2)
gmsh.model.occ.translate(half_ell_u, 0, L, 0)

half_ell_d = gmsh.model.occ.copy([(1, 5)])
gmsh.model.occ.rotate(half_ell_d, 0, 0, 0, 0, 0, 1, -pi/2)
gmsh.model.occ.translate(half_ell_d, 0, -L, 0)

qtr_ell_ur = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=8, angle1=pi, angle2=3*pi/2)
gmsh.model.occ.translate([(1, qtr_ell_ur)], L, L, 0)

qtr_ell_ul = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=9, angle1=3*pi/2, angle2=2*pi)
gmsh.model.occ.translate([(1, qtr_ell_ul)], -L, L, 0)

qtr_ell_ll = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=10, angle1=0, angle2=pi/2)
gmsh.model.occ.translate([(1, qtr_ell_ll)], -L, -L, 0)

qtr_ell_lr = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=11, angle1=pi/2, angle2=pi)
gmsh.model.occ.translate([(1, qtr_ell_lr)], L, -L, 0)

# Delete the template ellipses
gmsh.model.occ.remove([(1, half_ell_long), (1, half_ell_short)])
gmsh.model.occ.removeAllDuplicates()

# Connect the boundary ellipses (visualize the model first), start ccw from +x axis
b1 = gmsh.model.occ.addLine(4, 15)
b2 = gmsh.model.occ.addLine(14, 11)
b3 = gmsh.model.occ.addLine(10, 17)
b4 = gmsh.model.occ.addLine(16, 7)

b5 = gmsh.model.occ.addLine(6, 19)
b6 = gmsh.model.occ.addLine(18, 13)
b7 = gmsh.model.occ.addLine(12, 21)
b8 = gmsh.model.occ.addLine(20, 5)

# NOTE: elliptic arcs are ccw, so we add a negative sign.
b_loop = [3, b1, -8, b2, -6, b3, -9, b4, -4, b5, -10, b6, -7, b7, -11, b8]
b_loop_tag = gmsh.model.occ.addCurveLoop(b_loop)

# Create the magnet holes
m1 = gmsh.model.occ.addCircle(L/2, L/2, 0, R)
m2 = gmsh.model.occ.addCircle(-L/2, L/2, 0, R)
m3 = gmsh.model.occ.addCircle(-L/2, -L/2, 0, R)
m4 = gmsh.model.occ.addCircle(L/2, -L/2, 0, R)

m_c1 = gmsh.model.occ.addCurveLoop([m1])
m_c2 = gmsh.model.occ.addCurveLoop([m2])
m_c3 = gmsh.model.occ.addCurveLoop([m3])
m_c4 = gmsh.model.occ.addCurveLoop([m4])

gmsh.model.occ.removeAllDuplicates()

# Define the surfaces
ela_s = gmsh.model.occ.addPlaneSurface([b_loop_tag, ctr_ell_c, m_c1, m_c2, m_c3, m_c4])
m1_s = gmsh.model.occ.addPlaneSurface([m_c1])
m2_s = gmsh.model.occ.addPlaneSurface([m_c2])
m3_s = gmsh.model.occ.addPlaneSurface([m_c3])
m4_s = gmsh.model.occ.addPlaneSurface([m_c4])

gmsh.model.occ.removeAllDuplicates()

#################################################################################################
# Meshing

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)

#################################################################################################
# Rotate the unit cell by 90 ccw

points = gmsh.model.getEntities(0)
surfs = gmsh.model.getEntities(2)

gmsh.model.occ.rotate(points, 0, 0, 0, 0, 0, 1, pi/2) # prevents stray points
gmsh.model.occ.rotate(surfs, 0, 0, 0, 0, 0, 1, pi/2)

gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

#################################################################################################
# Duplicate the unit cell

for k in range(1, 10):
    temp = gmsh.model.occ.copy(gmsh.model.getEntities())
    gmsh.model.occ.translate(temp, 0, k*2*L, 0)

gmsh.model.occ.removeAllDuplicates()

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)

#################################################################################################
# Create the end segments

lo_half_ell = gmsh.model.occ.addEllipse(0, -L, 0, b, a, tag=1001, angle1=pi, angle2=2*pi)

lo_qtr_ell_r = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=1002, angle1=pi/2, angle2=pi)
gmsh.model.occ.rotate([(1, lo_qtr_ell_r)], 0, 0, 0, 0, 0, 1, pi/2) # rot-90
gmsh.model.occ.translate([(1, lo_qtr_ell_r)], L, -L, 0)

lo_qtr_ell_l = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=1003, angle1=0, angle2=pi/2)
gmsh.model.occ.rotate([(1, lo_qtr_ell_l)], 0, 0, 0, 0, 0, 1, -pi/2)
gmsh.model.occ.translate([(1, lo_qtr_ell_l)], -L, -L, 0)

gmsh.model.occ.synchronize()

# Visualize first to get tags

lo_l1 = gmsh.model.occ.addLine(19, 6)
lo_l2 = gmsh.model.occ.addLine(7, 16)
lo_l3 = gmsh.model.occ.addLine(374, 373)

lo_cl = gmsh.model.occ.addCurveLoop([lo_l1, lo_half_ell, lo_l2, lo_qtr_ell_l, lo_l3, lo_qtr_ell_r])
lo_s = gmsh.model.occ.addPlaneSurface([lo_cl])

lo_r = gmsh.model.occ.addRectangle(-L, -2*L-(a+w), 0, 2*L, 2*(a+w))

lo_s = gmsh.model.occ.fuse([(2, lo_s)], [(2, lo_r)], removeObject=True, removeTool=True)[0]

lo_h = gmsh.model.occ.addDisk(0, -2*L + 2*w, 0, 0.75*R, 0.75*R)

lo_end_s = gmsh.model.occ.cut(lo_s, [(2, lo_h)], removeObject=True, removeTool=False)[0]

gmsh.model.occ.synchronize()

hi_half_ell = gmsh.model.occ.addEllipse(0, 19*L, 0, b, a, tag=2001, angle1=0, angle2=pi)

hi_qtr_ell_r = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=2002, angle1=0, angle2=pi/2)
gmsh.model.occ.rotate([(1, hi_qtr_ell_r)], 0, 0, 0, 0, 0, 1, pi/2)
gmsh.model.occ.translate([(1, hi_qtr_ell_r)], L, 19*L, 0)

hi_qtr_ell_l = gmsh.model.occ.addEllipse(0, 0, 0, b, a, tag=2003, angle1=pi/2, angle2=pi)
gmsh.model.occ.rotate([(1, hi_qtr_ell_l)], 0, 0, 0, 0, 0, 1, -pi/2)
gmsh.model.occ.translate([(1, hi_qtr_ell_l)], -L, 19*L, 0)

gmsh.model.occ.synchronize()

# Visualize first
hi_l1 = gmsh.model.occ.addLine(368, 357)
hi_l2 = gmsh.model.occ.addLine(358, 359)
hi_l3 = gmsh.model.occ.addLine(391, 394)

hi_cl = gmsh.model.occ.addCurveLoop([hi_l1, hi_half_ell, hi_l2, hi_qtr_ell_r, hi_l3, hi_qtr_ell_l])
hi_s = gmsh.model.occ.addPlaneSurface([hi_cl])

hi_r = gmsh.model.occ.addRectangle(-L, 20*L-(a+w), 0, 2*L, 2*(a+w))

hi_s = gmsh.model.occ.fuse([(2, hi_s)], [(2, hi_r)], removeObject=True, removeTool=True)[0]

hi_h = gmsh.model.occ.addDisk(0, 20*L - 2*w, 0, 0.75*R, 0.75*R)

# Instead of removing the pin hole surfaces and introducing new 'pin' particles in LAMMPS, we apply fixes directly to the surfaces.
hi_end_s = gmsh.model.occ.cut(hi_s, [(2, hi_h)], removeObject=True, removeTool=False)[0]

gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

#################################################################################################
# Physical groups

n_surfs = len(gmsh.model.getEntities(2))
surfs = list(range(1,n_surfs+1))

ela = gmsh.model.addPhysicalGroup(2, surfs[0:50:5]+[lo_end_s[0][1]]+[hi_end_s[0][1]], name = "Elastic")
gmsh.model.addPhysicalGroup(2, surfs[1:50:5], name = "Magnet 2")
gmsh.model.addPhysicalGroup(2, surfs[2:50:5], name = "Magnet 3")
gmsh.model.addPhysicalGroup(2, surfs[3:50:5], name = "Magnet 4")
gmsh.model.addPhysicalGroup(2, surfs[4:50:5], name = "Magnet 1")
gmsh.model.addPhysicalGroup(2, [lo_h, hi_h], name = "Pins")

# Get surface tags of elastic group
ela_s_tags = gmsh.model.getEntitiesForPhysicalGroup(2, ela)

ela_boundary = []
for tag in ela_s_tags:
    b_tags = [abs(tg) for [dim, tg] in gmsh.model.getBoundary([(2, tag)])]
    ela_boundary.extend(b_tags)

# Remove the shared edges from stacking
from collections import Counter
ela_boundary = [tag for tag, count in Counter(ela_boundary).items() if count == 1]

ela_b_curves = gmsh.model.addPhysicalGroup(1, [curve for curve in ela_boundary], name = "elastic boundary")

#################################################################################################
# Re-meshing

gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.mesh.generate(2)

#################################################################################################
# Write to file

gmsh.write(write_path)

# To visualize the model we can run the graphical user interface with
# `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
# command line arguments:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.clear()
gmsh.finalize()

