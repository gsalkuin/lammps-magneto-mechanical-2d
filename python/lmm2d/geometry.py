"""
Read the .msh file created using Gmsh
"""
import os
import numpy as np
import meshio
import igl

class TriMesh:
    def __init__(self, dot_msh):
        """
        Parameters:
        dot_msh: <string> name of the .msh file without the extension
        NOTE:
        1) If using Gmsh v.4 ASCII format, check 'save all elements' but not 'save parametric coordinates'.
        2) If using physical groups, do not save all elements. Define physical groups for everything you want saved.
        """
        self.mesh = self.read_msh(dot_msh) # meshio.Mesh object

        self.V_into_xyz = self.mesh.points
        self.F_into_V = self.mesh.get_cells_type("triangle")
        self.E_into_V = igl.edge_flaps(self.F_into_V)[0]

        etol = 1e-2
        if np.min(self.E_lengths(self.E_into_V)) < etol*np.max(self.E_lengths(self.E_into_V)):
            raise Exception("Edge length too short. This may affect stability in simulations.")
        
        self.nV = self.V_into_xyz.shape[0]        
        self.nF = self.F_into_V.shape[0]
        self.nE = self.E_into_V.shape[0]

        self.F_into_V_per_block = [cell_block.data for cell_block in self.mesh.cells if cell_block.type == "triangle"]
        self.E_into_V_per_block = [igl.edge_flaps(F_into_V)[0] for F_into_V in self.F_into_V_per_block]
                
        self.n_cell_blocks = len(self.F_into_V_per_block)
        self.nF_per_block = [F_into_V.shape[0] for F_into_V in self.F_into_V_per_block]

        # returns ordered vertices of only the longest boundary loop (oriented ccw), e.g. no holes
        self.b_loop_Vs_per_block = [igl.boundary_loop(F_into_V) for F_into_V in self.F_into_V_per_block]

        # returns a list of lists containing at index i the adjacent vertices of vertex i
        self.V_adj_list = igl.adjacency_list(self.F_into_V) 

    def read_msh(self, dot_msh):
        cwd = os.getcwd()
        msh_file = os.path.join(cwd, dot_msh + '.msh')
        mesh = meshio.read(msh_file)
        return mesh

    def F_areas(self):
        f_areas = 0.5 * igl.doublearea(self.V_into_xyz, self.F_into_V)
        return f_areas
    
    def F_areas_per_block(self):
        return [0.5 * igl.doublearea(self.V_into_xyz, F_into_V) for F_into_V in self.F_into_V_per_block]

    def E_lengths(self, E_into_V):
        V1, V2 = E_into_V[:, 0], E_into_V[:, 1]
        P1, P2 = self.V_into_xyz[V1, :], self.V_into_xyz[V2, :] 
        E_lengths = np.linalg.norm(P1 - P2, axis=1)
        return E_lengths[:, None]
        
    def distance(self, V1, V2):
        P1, P2 = self.V_into_xyz[V1], self.V_into_xyz[V2]
        return np.linalg.norm(P1 - P2)

    # returns the shortest edge length connected to each vertex
    def V_shortest_edge(self):
        V_sE = np.zeros(self.nV)
        for V, adj_Vs in enumerate(self.V_adj_list):
            V_sE[V] = min([self.distance(V, adj_V) for adj_V in adj_Vs])
        return V_sE
    
