"""
Part and Composite classes 
"""
import numpy as np

# Creates a part from a cell block of a TriMesh object.
class _Part:
    def __init__(self, tri_mesh, cell_block_index, thickness, rho, Y_mod=None, m_vec=None):
        """
        tri_mesh: TriMesh object
        cell_block_index: <int> cell block # in tri_mesh
        thickness: <float>
        rho: <float>
        Y_mod: <float> or None
        m_vec: [<float>, <float>, <float>] or None
        """
        
        self.tri_mesh = tri_mesh

        self.magnet_flag = False
        
        # check if elastic or rigid/magnetic
        if Y_mod:
            if m_vec is not None or np.any(m_vec == 0.): 
                raise ValueError("Magnet cannot be elastic.")
            self.rigid_flag = False
            self.bond_flag = True
            self.Y_mod = Y_mod
        else:
            self.rigid_flag = True
            self.bond_flag = False

            if m_vec is not None:
                if m_vec.size != 3: raise ValueError("m_vec dimensions must be 3.")
                self.magnet_flag = True
                self.m_vec = m_vec

        self.thickness = thickness
        self.rho = rho
        
        self.F_into_V = tri_mesh.F_into_V_per_block[cell_block_index]
        self.F_areas = tri_mesh.F_areas_per_block()[cell_block_index]

        self.b_loop_Vs = tri_mesh.b_loop_Vs_per_block[cell_block_index]

        # Bonds
        if self.bond_flag:    
            self.E_into_V = tri_mesh.E_into_V_per_block[cell_block_index] 

            # Get boundary and interior edges
            E_is_boundary = np.all(np.isin(self.E_into_V, self.b_loop_Vs), axis=1)
            self.boundary_E_into_V = self.E_into_V[E_is_boundary]
            self.interior_E_into_V = self.E_into_V[~E_is_boundary]

            # Get interior bond and bond_coeffs list
            self.interior_bonds = self.interior_E_into_V
            self.interior_bond_coeffs = self._assemble_bond_coeffs(self.interior_bonds, self.Y_mod)

            # bonds on the boundary will only have half the stiffness since they will be double counted.
            self.boundary_bonds = self.boundary_E_into_V
            self.boundary_bond_coeffs = self._assemble_bond_coeffs(self.boundary_bonds, self.Y_mod / 2)

            self.bonds = np.concatenate((self.interior_bonds, self.boundary_bonds), axis=0)
            self.bond_coeffs = np.concatenate((self.interior_bond_coeffs, self.boundary_bond_coeffs), axis=0)

    # bond style harmonic; should contain everything except ID
    def _assemble_bond_coeffs(self, E_into_V, Y_mod):
        S_mod = Y_mod * self.thickness
        bond_coeff = np.sqrt(3) / 4 * S_mod * np.ones(E_into_V.shape[0])[:, None]
        return np.concatenate((bond_coeff, self.tri_mesh.E_lengths(E_into_V)), axis=1)
    

#######################################################################################################################
# Assembles all the parts into one composite object.
class Composite:
    def __init__(self, tri_mesh, thickness_per_cell_block, rho_per_cell_block, Y_mod_per_cell_block=None, m_vec_per_cell_block=None):
        """
        tri_mesh: TriMesh object
        thickness_per_cell_block: list of <float>
        rho_per_cell_block: list of <float>
        Y_mod_per_cell_block: list of either <float> or None
        m_vec_per_cell_block: list of either [<float>, <float>, <float>] or None
        """
        self.tri_mesh = tri_mesh
        self.mesh = tri_mesh.mesh # meshio Mesh object
        self.n_parts = tri_mesh.n_cell_blocks
        
        # consistency checks
        if Y_mod_per_cell_block is None: Y_mod_per_cell_block = [None] * self.n_parts
        elif isinstance(Y_mod_per_cell_block, float): Y_mod_per_cell_block = [Y_mod_per_cell_block] * self.n_parts        
        if len(Y_mod_per_cell_block) != self.n_parts: raise ValueError("Number of parts does not match number of moduli.")

        if m_vec_per_cell_block is None: m_vec_per_cell_block = [None] * self.n_parts
        if len(m_vec_per_cell_block) != self.n_parts: raise ValueError("Number of parts does not match number of dipoles.")

        if isinstance(thickness_per_cell_block, float): thickness_per_cell_block = [thickness_per_cell_block] * self.n_parts            
        if len(thickness_per_cell_block) != self.n_parts: raise ValueError("Number of parts does not match number of thicknesses.")
        
        if isinstance(rho_per_cell_block, float): rho_per_cell_block = [rho_per_cell_block] * self.n_parts
        if len(rho_per_cell_block) != self.n_parts: raise ValueError("Number of parts does not match number of densities.")

        # instantiate the parts
        self.parts = []
        for index in range(self.n_parts):
            self.parts.append(_Part(tri_mesh, index, thickness_per_cell_block[index], rho_per_cell_block[index], Y_mod_per_cell_block[index], m_vec_per_cell_block[index]))

        # Bonds
        try:
            self.bonds = np.concatenate([part.bonds for part in self.parts if part.bond_flag], axis=0)
            self.bond_coeffs = np.concatenate([part.bond_coeffs for part in self.parts if part.bond_flag], axis=0)
        except:
            self.bonds = None
            self.bond_coeffs = None

	    # Atoms
        self.V_into_xyz = tri_mesh.V_into_xyz
        self.n_atoms = tri_mesh.nV

        # initialize charges and dipoles, types, etc.
        self.V_into_q_mu = np.zeros((self.n_atoms, 4), float) # q, mux, muy, muz       
        self.V_into_atom_type = np.ones(self.n_atoms, int)
        self.atom_types = [1, 2]
        self.n_atom_types = len(self.atom_types)
        self.rigid_flag = False
        self.magnet_flag = False

        self.V_into_diam = np.zeros((self.n_atoms), float) # 0 = point mass
        self.V_into_rho = np.zeros(self.n_atoms, float) # mass or density, depending on diam
        
        ####################################################################################################################
        # Assign the part IDs = molecule IDs in LAMMPS
        self.V_into_part_ID = np.zeros(self.n_atoms, int)
        
        # for now, we just assign a part-ID in sequence
        self.part_IDs = list(range(1, self.n_parts+1)) 

        # rigid part-IDs take precedence on rigid-elastic boundary atoms otherwise the boundary bonds will be incorrect
        # for atoms on rigid-rigid boundaries, the larger part-ID is kept
        self.rigid_flags = np.array([part.rigid_flag for part in self.parts], dtype = int)
        for index in np.argsort(self.rigid_flags):
            self.set_part_ID(self.parts[index], self.part_IDs[index])

        # checks for unassigned/stray atoms indicating a problem with the geometry/mesh
        if self.V_into_part_ID.any == 0: 
            print("WARNING: Atoms with mol-ID = 0 found! Please check the .lam file for stray atoms.")
        
        ####################################################################################################################
        # Rigidify
        self._rigidify() # may add new atom types

        # rigid IDs (unnecessary if we use part-ID to distinguish rigid objects)
        self.V_into_rigid_ID = None # by default, rigid_ID = part_ID

        ####################################################################################################################
        # create dipoles
        self.n_magnets = self._magnetize()
        
        # set the mass/density
        self.compute_rho_from_mesh()

        
    #################################################################################################################### 

    def set_part_ID(self, part, ID):
        for _, vertex in np.ndenumerate(part.F_into_V):
            self.V_into_part_ID[vertex] = ID

    def set_rigid_ID(self, list_of_part_IDs, rigid_ID):
        """
        set a custom ID for each part; creates a <filename>-rigid-ID.txt. 
        This can be used to define an atomfile-style variable in LAMMPS, e.g.:
            ...
            variable rigidID atomfile atomfile.txt
            fix ID group-ID rigid custom v_rigidID
            ...
        """
        if self.V_into_rigid_ID is None:
            self.V_into_rigid_ID = np.copy(self.V_into_part_ID) 
            # set all elastic atoms to 0
            self.V_into_rigid_ID = np.where(np.isin(self.V_into_atom_type, [1,2]), 0, self.V_into_rigid_ID)
 
        self.V_into_rigid_ID = np.where(np.isin(self.V_into_part_ID, list_of_part_IDs), rigid_ID, self.V_into_rigid_ID)
        return

    def _rigidify(self):
        """
        WARNING: This will overwrite the atom types of rigid parts.
        """
        for part in self.parts:
            if part.rigid_flag:
                self.rigid_flag = True
                vertices = part.F_into_V.flatten()
                self.V_into_atom_type[vertices] = 3
                self.V_into_atom_type[part.b_loop_Vs] = 4      

        if self.rigid_flag:
            self.atom_types.append(3) # bulk, no charge -> delete in LAMMPS
            self.atom_types.append(4) # boundary, no charge if dipole is used

        self.n_atom_types = len(self.atom_types)
        return   

    # dipole model
    def _magnetize(self):
        n_magnets = 0

        # set the diameter to be very small to not significantly affect MOI
        # we only need the dipole to not be a point mass, no need to update rho (= type 3)
        mdiam = np.min(self.bond_coeffs[:,1]) # shortest bond
        for part in self.parts:
            if part.magnet_flag:
                self.magnet_flag = True
                n_magnets += 1
                # find the atom closest to the center of the part
                Vs = np.unique(part.F_into_V.flatten())
                x_ctr = np.mean(self.V_into_xyz[Vs], axis=0)
                V_ctr = Vs[np.argmin(np.linalg.norm(self.V_into_xyz[Vs]-x_ctr, axis=1))]
                assert self.V_into_atom_type[V_ctr] == 3, "The dipole must be an interior atom of the magnet."
                self.V_into_atom_type[V_ctr] = 5 # dipole
                self.V_into_q_mu[V_ctr] = np.array([0.0, *part.m_vec])
                self.V_into_diam[V_ctr] = mdiam

        if self.magnet_flag:
            self.atom_types.append(5)
            print("{} atoms turned into dipoles with diameter {}.".format(n_magnets, mdiam))

        return n_magnets

    def set_rho_per_type(self, type, rho):
        """
        Override the density computed using compute_rho_from_mesh() for the specified atom type.
        type: <int> 1-5
        rho: <float>
        """
        self.V_into_rho[np.where(np.asarray(self.V_into_atom_type, int) == type)] = rho
        return

    def set_diam_per_type(self, type, diam):
        """
        Set the diameter of the given atom type.
        type: <int> 1-5
        diam: <float>
        """
        self.V_into_diam[np.where(np.asarray(self.V_into_atom_type, int) == type)] = diam
        return
    
    def set_diam_bVs(self, atom_type=2):
        """
        Set the diameter of each boundary atom with the specified atom type based on the distance of its two neighboring boundary atoms.
        Boundary must be a disjoint union of simply connected closed curves.
        NOTES: 
        1) only works if each boundary vertex has exactly two adjacent boundary vertices; 
           make sure non-neighboring b vertices are not part of the same triangle element
        2) use special_bonds to exclude interactions between neighboring boundary atoms (may not work with granular)
        """
        boundary_Vs = np.where(self.V_into_atom_type == atom_type)[0]
        bV_adj_list = [self.tri_mesh.V_adj_list[V] for V in boundary_Vs]

        # get only the adjacent vertices on the boundary
        bV_adj_bV = [] # final shape must be (n_bV, 2)
        for adj_list in bV_adj_list:
            adj_V = np.intersect1d(adj_list, boundary_Vs, assume_unique=True)
            if adj_V.size != 2: raise Exception("Boundary vertex should have exactly two adjacent boundary vertices.")
            bV_adj_bV.append(adj_V)
        
        bV_adj_bV = np.asarray(bV_adj_bV).reshape(-1, 2)

        # get the edge indices
        E_into_V = np.sort(self.tri_mesh.E_into_V, axis=1)
        bE1 = np.zeros((bV_adj_bV.shape), int)
        bE1[:, 0] = boundary_Vs
        bE2 = bE1.copy()
        bE1[:, 1] = bV_adj_bV[:, 0]
        bE2[:, 1] = bV_adj_bV[:, 1]

        # get the edge indices in E_into_V for each boundary vertex
        bV_adj_E1 = np.argwhere((np.expand_dims(np.sort(bE1, axis=1), 1) == E_into_V).all(axis=2))[:, 1]
        bV_adj_E2 = np.argwhere((np.expand_dims(np.sort(bE2, axis=1), 1) == E_into_V).all(axis=2))[:, 1]

        # get the edge lengths
        E_lengths = self.tri_mesh.E_lengths(self.tri_mesh.E_into_V)
        E1_lengths = E_lengths[bV_adj_E1]
        E2_lengths = E_lengths[bV_adj_E2]

        # radii can overlap if bonded atoms are excluded from contact
        self.V_into_diam[boundary_Vs] = np.maximum(E1_lengths, E2_lengths).flatten()
        return
    
    def compute_rho_from_mesh(self):
        # compute the point masses from the mesh
        if np.any(self.V_into_rho != 0.):
            print("WARNING: Masses have previously been set and will be reset.")
            self.V_into_rho = np.zeros(self.n_atoms, float)

        for part in self.parts:
            for face, verts in enumerate(part.F_into_V):
                self.V_into_rho[verts] += part.rho * part.thickness * part.F_areas[face] / 3.
        
        # convert to density if not a point mass
        sph_vol = np.pi * self.V_into_diam**3 / 6.
        self.V_into_rho = np.where(self.V_into_diam == 0., self.V_into_rho, self.V_into_rho/sph_vol)
        return


