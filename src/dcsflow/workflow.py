# This script is the workflow for INS spectrum analysis - DCS-flow

import ase.io
import os
import argparse
import glob
import subprocess
import shutil
import dcsflow.vasp as va
import dcsflow.dftb as db
import dcsflow.phonons as ph
import dcsflow.oclimax as ocl
import dcsflow.restart_dftbmd as re
import dcsflow.restart_collector as rc

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    Args:
    path: The current directory (use os.getcwd())
    #########################################
    Return:
    file[0]: The structure file in the path
    """
    file = glob.glob(path + "/*.cif") + glob.glob(path + "/POSCAR") 
    
    return file[0]

# 1st flow: DFT-FD
def dftfd(local, gpu, hpc, kpts, disp, fmax, mesh, supercell):
    """ Run DFT simulation with VASP, Phonopy and Oclimax to get INS spectrum
    local (boolean): If True, run VASP on local machine, otherwise run on HPC 
    gpu (boolean): If True, run VASP with GPU, otherwise run on CPU 
    hpc (list): srun setting for VASP simulation, Defaults: srun -n {8} -c {32} -G {8} --cpu-bind=cores --gpu-bind=none vasp_std'
    kpts (list): number of KPOINTS. (for Phonon simulation, Defaults to [1,1,1], need user input for relaxation.)
    disp (boolean): If True, set ivdw=12, which is DFT-D3 method with Becke-Johnson damping function (Defaults to True)
    fmax: Maximum allowed force for convergence between atoms (Defaults to 1e-2)
    mesh (list): Mesh size for Phonopy eigenvector problems. (Defaults to [8,8,8])
    supercell: Supercell size. (Defaults to [2,2,2])
    """
    main_path = os.getcwd()
    if not os.path.exists(os.path.join(main_path, "1-relax" ,"CONTCAR")): # Check whether the relaxation simulation is finished. 
        print(" Writing VASP input files! ")
        os.makedirs(os.path.join(main_path, "1-relax"), exist_ok=True)
        os.chdir(os.path.join(main_path, "1-relax")) # Change directory to 1-relax
        atom = ase.io.read(getGeometry(main_path))
        va.dft(disp, fmax, kpts, atoms=atom, mode=1)
        if local:
            subprocess.run('vasp_std', shell=True) # Run VASP on local machine
        if gpu:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} -G {hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std', shell=True) # Run VASP on hpc gpu nodes
        else:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} --cpu-bind=cores vasp_std', shell=True) # Run VASP on hpc cpu nodes
        os.chdir(main_path)
            
    else:
        print(" The relaxation is already finished, please continue simulating phonons!!! ")

    ph_kpts = [1,1,1] # The kpoints for phonon simulation
    if not os.path.exists(os.path.join(main_path, "2-phonons", "mesh.yaml")): # Check whether the phonon simulation is finished. 
        print(" Using Phonopy to create displaced structures! ")
        src = os.path.join(main_path, "1-relax', 'CONTCAR") # source file
        dst = os.path.join(main_path, "2-phonons") # destination
        os.makedirs(dst, exist_ok=True)
        os.chdir(dst)
        shutil.copy(src, dst+'POSCAR')
        
        ph.displace_structure(calc='vasp', supercell=[2,2,2]) # Create displaced POSCAR files by Phonopy
        root, dirs, files = next(os.walk(dst))
        
        for file in files:
            folder_name = file.split('-')[1]  # Extracts '009' from 'POSCAR-009'
            folder_path = os.path.join(dst, folder_name)
            os.makedirs(folder_path, exist_ok=True) # Create the folder if it doesn't exist
    
            source = os.path.join(dst, file) # Define the source and destination file paths
            new_file = "POSCAR"  # Rename the file to 'POSCAR'
            destination = os.path.join(folder_path, new_file)
            shutil.move(source, destination) # Move the file into the corresponding folder
                
            root2, dirs2, files2 = next(os.walk(dst)) # Loop through every sub-folder in 2-phonons to write VASP input files and run script
            for d in dirs2:
                os.chdir(d)
                atom = ase.io.read('POSCAR')
                va.dft(disp, fmax, kpts=ph_kpts, atoms=atom, mode=2)
                if local:
                    subprocess.run('vasp_std', shell=True) # Run VASP on local machine
                if gpu:
                    subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} -G {hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std', shell=True) # Run VASP on hpc gpu nodes
                else:
                    subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} --cpu-bind=cores vasp_std', shell=True) # Run VASP on hpc cpu nodes
                
                os.chdir('..')
        
        ph.get_force_sets(calc='vasp')
        ph.run_mesh(mesh)

    else:
        print(" The phonon simulation is already finished!!! ")
    if not os.path.exists(os.path.join(main_path, "3-oclimax")):
        ocl.run_oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=1)

# 2nd flow: DFTB-FD
def dftbfd(kpts, disp, fmax, mesh, supercell):
    """ Run DFTB+, Phonopy and Oclimax to get INS spectrum
    kpts (list): number of KPOINTS. (for Phonon simulation, Defaults to [1,1,1], need user input for relaxation.)
    disp (boolean): If True, which is DFT-D3 method with Becke-Johnson damping function (Defaults to True)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-2)
    mesh (list): Mesh size for Phonopy eigenvector problems. (Defaults to [8,8,8])
    supercell: Supercell size. (Defaults to [2,2,2])
    """
    main_path = os.getcwd()
    if not os.path.exists(os.path.join(main_path, "1-relax", "geo_end.gen")): # Check whether the relaxation simulation is finished. 
        print(" Writing DFTB+ input files! ")
        os.makedirs(os.path.join(main_path, "1-relax"), exist_ok=True)
        os.chdir(os.path.join(main_path, "1-relax")) # Change directory to 1-relax
        atom = ase.io.read(getGeometry(main_path))
        db.relax(disp,fmax,kpts,atoms=atom)
        subprocess.run('ulimit -s unlimited', shell=True)
        subprocess.run('dftb+ > relax.out', shell=True)

    else:
        print(" The Relaxation is finished, skip to the step 2! ")

    if not os.path.exists(os.path.join(main_path, "2-phonons" ,"mesh.yaml")): # Check whether the phonon simulation is finished. 
        print(" Using Phonopy to create displaced structures! ")
        src = os.path.join(main_path, '1-relax', 'geo_end.gen') # source file
        dst = os.path.join(main_path, '2-phonons') # destination
        os.makedirs(dst, exist_ok=True)
        os.chdir(dst)
        shutil.copy(src, dst)
        
        ph.displace_structure(calc='dftb+', supercell=[2,2,2]) # Create displaced geo.genS files
        root, dirs, files = next(os.walk(dst))
        
        for file in files:
            folder_name = file.split('-')[1]  # Extracts '009' from 'geo.genS-009'
            folder_path = os.path.join(dst, folder_name)
            os.makedirs(folder_path, exist_ok=True) # Create the folder if it doesn't exist
    
            source = os.path.join(dst, file) # Define the source and destination file paths
            new_file = "geo_end.gen"  # Rename the file to 'POSCAR'
            destination = os.path.join(folder_path, new_file)
            shutil.move(source, destination) # Move the file into the corresponding folder
        
            root2, dirs2, files2 = next(os.walk(dst))
            for d in dirs2:
                os.chdir(d)
                atom = ase.io.read("geo_end.gen")
                db.force(disp, kpts=ph_kpts, atoms=atom)
                subprocess.run('dftb+ > force.out', shell=True)
                os.chdir('..')
        ph.get_force_sets(calc='dftb+')
        ph.run_mesh(mesh)
    
    else:
        print(" The phonon simulation is already done!!! ")
    
    if not os.path.exists(os.path.join(main_path, "3-oclimax")):
        ocl.run_oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=1)

# 3rd flow: DFT-MD

def dftmd(local, gpu, hpc, kpts, disp, max, nsw1, nsw2, steps, supercell, tebeg, teend):
    """ Run DFT-MD simulation with VASP and Oclimax to get INS spectrum
    local (boolean): If True, run VASP on local machine, otherwise run on HPC 
    gpu (boolean): If True, run VASP with GPU, otherwise run on CPU 
    kpts (list): number of KPOINTS.
    disp (boolean): If True, set ivdw=12, which is DFT-D3 method with Becke-Johnson damping function (Defaults to True)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-2)
    nsw1 (int): Number of NVT MD steps, Defaults to 4000 (1ps)
    nsw2 (int): Number of NVE MD steps, Defaults to 4000 (1ps)
    steps (int): Multiple of NVE-MD steps; (i.e. 4000*steps is the total steps for NVE-MD) Defaults to 10
    supercell: Supercell size. (Defaults to [2,2,2])
    tebeg (float): The initial temperature for MD simulation, Defaults to 150K
    teend (float): The final temperature for MD simulation, Defaults to 150K. this is useful when training MLFF using NPT-MD
    """
    main_path = os.getcwd()
    if not os.path.exists(os.path.join(main_path, "1-relax", "CONTCAR")): # Check whether the relaxation simulation is finished. 
        print(" Writing VASP input files! ")
        os.makedirs(os.path.join(main_path, "1-relax"), exist_ok=True)
        os.chdir(os.path.join(main_path, "1-relax")) # Change directory to 1-relax
        atom = ase.io.read(getGeometry(main_path))
        va.dft(disp, fmax, kpts, atoms=atom, mode='relax')
        if local:
            subprocess.run('vasp_std', shell=True) # Run VASP on local machine
        if gpu:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} -G {hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std', shell=True) # Run VASP on hpc gpu nodes
        else:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} --cpu-bind=cores vasp_std', shell=True) # Run VASP on hpc cpu nodes
        
        os.chdir(main_path)
            
    else:
        print(" The relaxation is already finished, Continue to NVT-MD !!! ")

    md_kpts = [1,1,1] # The kpoints for MD simulation
    if not os.path.exists(os.path.join(main_path, "2-nvtmd", "vasprun.xml")):
        print(" Writing VASP input files! ")
        os.makedirs(os.path.join(main_path, "2-nvtmd"), exist_ok=True)
        atom = ase.io.read(os.path.join(main_path,'1-relax','CONTCAR'))
        atom *= supercell
        va.md(disp,nsw1,tebeg,teend,kpts=md_kpts,atoms=atom,ensemble=1)
        if local:
            subprocess.run('vasp_std', shell=True) # Run VASP on local machine
        if gpu:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} -G {hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std', shell=True) # Run VASP on hpc gpu nodes
        else:
            subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} --cpu-bind=cores vasp_std', shell=True) # Run VASP on hpc cpu nodes

    else:
        print(" The thermalization (NVT-MD) is already finished, Continue to NVE-MD !!! ")

    if not os.path.exists(os.path.join(main_path, "3-nvemd")):
        print(" Writing VASP input files! ")
        os.makedirs(os.path.join(main_path, "3-nvemd"), exist_ok=True)
        atom = ase.io.read(os.path.join(main_path,'2-nvtmd','CONTCAR'))
        os.makedirs(os.path.join(main_path, "3-nvemd", "1"), exist_ok=True)
        va.md(disp,nsw2,tebeg,teend,kpts=md_kpts,atoms=atom,ensemble=2)
        subprocess.run('vasp_std', shell=True)
        os.chdir('..')
        for i in range(1,steps):
            atom = ase.io.read('3-nvemd/'+str(i)+'/CONTCAR')
            os.makedirs(os.path.join(main_path, "3-nvemd", str(i+1)), exist_ok=True)
            os.chdir(os.path.join(main_path, "3-nvemd", str(i+1)))
            va.md(disp,nsw2,temp,kpts=md_kpts,atoms=atom,ensemble=2)
            if local:
                subprocess.run('vasp_std', shell=True) # Run VASP on local machine
            if gpu:
                subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} -G {hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std', shell=True) # Run VASP on hpc gpu nodes
            else:
                subprocess.run(f'srun -n {hpc[0]} -c {hpc[1]} --cpu-bind=cores vasp_std', shell=True) # Run VASP on hpc cpu nodes
            os.chdir('..')

    else:
        print(" The (NVE-MD) is already finished, Continue to Oclimax !!! ")
    
    if not os.path.exists(os.path.join(main_path, "4-oclimax")):
        ocl.run_oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=2)

# 4th flow: DFTB-MD
def dftbmd(kpts,disp,fmax,nsw1,nsw2,steps,supercell,temp):
    """ Run DFTB-MD simulation with DFTB+ and Oclimax to get INS spectrum
    kpts (list): number of KPOINTS.
    disp (boolean): If True, set DFT-D3 method with Becke-Johnson damping function (Defaults to True)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-2)
    nsw1 (int): Number of NVT MD steps, Defaults to 4000 (1ps)
    nsw2 (int): Number of NVE MD steps, Defaults to 4000 (1ps)
    steps (int): Multiple of NVE-MD steps; (i.e. 4000*steps is the total steps for NVE-MD) Defaults to 10
    supercell: Supercell size. (Defaults to [2,2,2])
    temp (float): The temperature for MD simulation, Defaults to 150K
    """
    main_path = os.getcwd()
    if not os.path.exists(os.path.join(main_path, "1-relax", "geo_end.gen")): # Check whether the relaxation simulation is finished. 
        print(" Writing DFTB+ input files! ")
        os.makedirs(os.path.join(main_path, "1-relax"), exist_ok=True)
        os.chdir(os.path.join(main_path, "1-relax")) # Change directory to 1-relax
        atom = ase.io.read(getGeometry(main_path))
        db.relax(disp,fmax,kpts,atoms=atom)
        subprocess.run('ulimit -s unlimited', shell=True)
        subprocess.run('dftb+ > relax.out', shell=True)

    else:
        print(" The Relaxation is finished, skip to the step 2! ")

    md_kpts = [1,1,1]
    if not os.path.exists(os.path.join(main_path, "2-nvtmd", "nvt.out")): # Check whether the nvtmd simulation is finished.
        os.makedirs(os.path.join(main_path, "2-nvtmd"), exist_ok=True)
        atom = ase.io.read(os.path.join(main_path, '1-relax', 'geo_end.gen'))
        atom *= supercell
        db.md(disp,nsw1,temp,kpts=md_kpts,atoms=atom,ensemble=1)
        subprocess.run('dftb+ > md.out', shell=True)
                              
    else:
        print(" The NVT-MD simulation is already finished!!! ")
                
    if not os.path.exists(os.path.join(main_path, "3-nvemd", "f{steps}")):
        print(" Writing DFTB+ input files! ")
        os.makedirs(os.path.join(main_path, "3-nvemd"), exist_ok=True) # Here we will create subfolders for NVE-MD simulation.
        atom = ase.io.read(os.path.join(main_path,'2-nvtmd','geo_end.gen'))
        os.makedirs(os.path.join(main_path, "3-nvemd", "1"), exist_ok=True)
        db.md(disp,nsw2,temp,kpts=md_kpts,atoms=atom,ensemble=2)
        subprocess.run('dftb+ > nve.out', shell=True)
        re.make_files(max_iter,extra_files,output_dir,
                    self_copy,write_over,force_restart,
                    restart_from)
        os.chdir('..')
        for i in range(1,steps):
            src = os.path.join(main_path, "3-nvemd", str(i), 'restart') # source folder
            dst = os.path.join(main_path, "3-nvemd", str(i+1)) # destination folder
            shutil.copytree(src, dst)
            os.chdir(os.path.join(main_path, "3-nvemd", str(i+1)))
            subprocess.run('dftb+ > md.out', shell=True)
            re.make_files(max_iter,extra_files,output_dir,
                        self_copy,write_over,force_restart,
                        restart_from)
            os.chdir('..')

        # Have to collect the trajectory files and put it together.
        rc.collect(extra_files,
                output_dirname,
                input_dirname,
                properties,
                lattice,
                add_mode,
                consequtive
                )

    else:
        print(" The NVE-MD simulation is already finished!!! ")
    
    if not os.path.exists(os.path.join(main_path, "4-oclimax")):
        ocl.run_oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=3)