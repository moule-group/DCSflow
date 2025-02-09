########## Phonons module ##########

import os
import numpy as np
import phonopy
import subprocess
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure
from phonopy.interface.calculator import write_supercells_with_displacements

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
    if calc == 'dftbp': 
        path = 'geo_end.gen' 
    
    unitcell, _ = read_crystal_structure(path, interface_mode=calc)
    
    supercell_matrix = np.diag(supercell)
    phonon = Phonopy(unitcell,
                     supercell_matrix=supercell_matrix,
                     )
    phonon.generate_displacements()
    cells_with_disps = phonon.supercells_with_displacements

    write_supercells_with_displacements(
        calc,
        phonon.supercell,
        cells_with_disps,
        optional_structure_info=None, # Default in Phonopy is None
        additional_info=None, # Default in Phonopy is None
    )
    phonon.save('phonopy_disp.yaml')

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
    root, dirs, files = next(os.walk(os.getcwd()))
    dirs = [d for d in dirs if len(d) <= 3]
    max_dir = max(dirs)
    if calc == 'vasp':
        subprocess.run(f"phonopy -f {{001..{max_dir}}}/vasprun.xml", shell=True)

    if calc == 'dftbp':
        subprocess.run(f"phonopy -f {{001..{max_dir}}}/results.tag --dftb+", shell=True)

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