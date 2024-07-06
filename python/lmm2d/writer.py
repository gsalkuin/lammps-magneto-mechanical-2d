"""
Writes the LAMMPS data file.
"""
import numpy as np

def write_lam(list_of_composites, filename):
    
    n_atoms = [] # list of nV for each composite, used to shift the global atom index
    atom_types = []
    
    V_into_xyz = []
    V_into_mol_ID = []
    V_into_atom_type = []
    V_into_q_mu = []
    V_into_rho = []
    V_into_diam = []
    V_into_rigid_ID = []
    
    bond_flag = False
    bonds = []
    bond_coeffs = []

    mol_ID_offset = 0
    for composite in list_of_composites:
        V_into_xyz.append(composite.V_into_xyz)

        V_into_mol_ID.append(composite.V_into_part_ID + mol_ID_offset * np.ones(composite.n_atoms, int))
        mol_ID_offset += composite.n_parts

        V_into_atom_type.append(composite.V_into_atom_type)
        atom_types.append(composite.atom_types)

        V_into_q_mu.append(composite.V_into_q_mu)
        V_into_rho.append(composite.V_into_rho)
        V_into_diam.append(composite.V_into_diam)

        if composite.bonds is not None:
            bond_flag = True
            bonds.append(composite.bonds + np.sum(n_atoms)) # offset the local atom index for global bond topology
            bond_coeffs.append(composite.bond_coeffs)

        n_rigid = 0
        if composite.V_into_rigid_ID is not None:
            V_into_rigid_ID.append(composite.V_into_rigid_ID + n_rigid)
            n_rigid += np.max(composite.V_into_rigid_ID) # not necessarily the number of rigid bodies
            
        n_atoms.append(composite.n_atoms) 
        
    N_atoms = np.sum(n_atoms)

    atom_types = np.concatenate(atom_types, axis=0)
    n_atom_types = np.unique(atom_types).size

    V_into_xyz = np.concatenate(V_into_xyz, axis=0)
    V_into_mol_ID = np.concatenate(V_into_mol_ID, axis=0)[:, None]
    V_into_atom_type = np.concatenate(V_into_atom_type, axis=0)[:, None]
    V_into_q_mu = np.concatenate(V_into_q_mu, axis=0)
    V_into_rho = np.concatenate(V_into_rho, axis=0)[:, None]
    V_into_diam = np.concatenate(V_into_diam, axis=0)[:, None]

    if V_into_rigid_ID != []:
        V_into_rigid_ID = np.concatenate(V_into_rigid_ID, axis=0)[:, None]
        V_into_rigid_ID[V_into_atom_type == 1] = 0 # elastic nodes

        # atomfile variable style 
        with open(filename + '-rigid-ID.txt', 'w') as file:    
            file.write(str(N_atoms) + '\n')
            data = np.concatenate((np.arange(1, N_atoms+1)[:, None], V_into_rigid_ID), axis=1)
            template = '%i %f \n' * N_atoms
            file.write(template % tuple(data.ravel()))
            print('Generated the rigid-ID atomfile ' + filename + '-rigid-ID.txt')
    
    if bond_flag:
        bonds = np.concatenate(bonds, axis=0)
        bond_coeffs = np.concatenate(bond_coeffs, axis=0)
        N_bonds = bonds.shape[0] 

    # data file
    with open(filename + '.lam', 'w') as file:
        file.write('LAMMPS DATA FILE \n\n')

        file.write(str(N_atoms) + ' atoms \n')
        if bond_flag:
            file.write(str(N_bonds) + ' bonds \n\n')

        file.write(str(n_atom_types) + ' atom types \n')
        if bond_flag:
            file.write(str(N_bonds) + ' bond types \n\n')
      
        # create box
        x, y =  V_into_xyz[:,0], V_into_xyz[:,1]
        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        x_range, y_range = x_max - x_min, y_max - y_min
        z_range = np.min(bond_coeffs[:,1])

        file.write('%f %f xlo xhi \n'%(x_min - 0.1*x_range, x_max + 0.1*x_range))
        file.write('%f %f ylo yhi \n'%(y_min - 0.1*y_range, y_max + 0.1*y_range))
        file.write('%f %f zlo zhi \n\n'%(-z_range/2, z_range/2))

        if bond_flag:
            file.write('Bond Coeffs # harmonic \n\n')
            ids = 1 + np.arange(N_bonds)[:, None]
            bond_coeffs = np.concatenate((ids, bond_coeffs), axis=1)
            template = '%i %g %f \n' * N_bonds
            file.write(template % tuple(bond_coeffs.ravel()))
            file.write('\n')

        file.write('Atoms # hybrid dipole molecular sphere \n\n')
        ids = 1 + np.arange(N_atoms)[:, None]
        atoms = np.concatenate((ids, V_into_atom_type, V_into_xyz, V_into_q_mu, V_into_mol_ID, V_into_diam, V_into_rho), axis=1)
        template = '%i %i %f %f %f %i %g %g %g %i %g %g \n' * N_atoms # atom-id atom-type x y z q mux muy muz mol-ID diam rho
        file.write(template % tuple(atoms.ravel()))
        file.write('\n')

        if bond_flag:
            file.write('Bonds \n\n')
            bonds = bonds + 1 # re-index to 1-based
            ids = 1 + np.arange(N_bonds)[:, None]
            bonds = np.concatenate((ids, ids, bonds), axis=1)
            template = '%i %i %i %i \n' * N_bonds
            file.write(template % tuple(bonds.ravel()))
            file.write('\n')

    print('Generated the LAMMPS file ' + filename + '.lam')
    return

