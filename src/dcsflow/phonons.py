########## Phonons module ##########

import os
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure
from phonopy.cui.create_force_sets import create_FORCE_SETS

def displace_structure(calc, supercell):
    """ Run Phonopy to create displaced supercell structure
    Arg: 
    calc (str): Calculator mode (Options: VASP, DFTB+, ...)
    supercell (list): Supercell size. 
    #######################################
    Return:
    Displaced structure files named POSCAR-001 ... 
    """
    # Read the relaxed structure (CONTCAR)
    if calc == 'vasp':
        path = 'CONTCAR' 
    if calc == 'dftb+': 
        path = 'geo_end.gen' 
    
    unitcell, _ = read_crystal_structure(path, interface_mode=calc)
    
    phonon = Phonopy(unitcell,
                     supercell_matrix=[[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]],
                     primitive_matrix=[[0, 0.5, 0.5],
                                       [0.5, 0, 0.5],
                                       [0.5, 0.5, 0]])
    phonon.generate_displacements()
    supercells = phonon.supercells_with_displacements

def get_force_sets(calc):
    """ Obtain FORCE_SETS file for specific calculator from Phonopy.
    Docs:
    FORCE_SETS file format
    1st line: number of atoms.
    2nd line: number of supercell displacement (Each supercell only involves one displaced atom).
    In each set:
    First, the atom number in supercell is written. 
    Second, the atomic displacement in Cartesian coordinates is written. 
    Last, atomic forces in Cartesian coordinates are successively written. 
    Args:
    calc (str): Calculator mode (Options: VASP, DFTB+, ...)
    #########################################
    Return:
    FORCE_SETS file
    """
    dirlist = sorted([x.name for x in os.scandir() if x.is_dir()])
    
    if calc == 'vasp':
        filepath = [x + '/vasprun.xml' for x in dirlist]

    if calc == 'dftb':
        filepath = [x + '/results.tag' for x in dirlist]
        
    create_FORCE_SETS(interface_mode=calc,
                      force_filenames=filepath,
                      disp_filename='phonopy_disp.yaml')

def run_mesh(mesh):
    """ Run Phonopy to create mesh.yaml file which saves normal mode frequencies and eigenvectors.
    mesh (list): Mesh size for Phonopy eigenvector problems. 
    """
    phonon = phonopy.load('phonopy_disp.yaml')
    force_constant = phonon.produce_force_constants()
    print(" Creating phonopy_params.yaml file which saves force constants. ")
    phonon.save(settings={'force_constants': True})
    print(" Done creating phonopy_params.yaml. ")
    
    phonon.run_mesh(mesh,with_eigenvectors=True)
    print(" Creating mesh.yaml file which saves normal mode frequencies. ")
    phonon.write_yaml_mesh() 
    print(" Done creating mesh.yaml. ")