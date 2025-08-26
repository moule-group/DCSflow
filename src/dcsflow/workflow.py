# This script is the workflow for INS spectrum analysis - DCS-flow

import ase.io
import time
import os
import argparse
import glob
import subprocess
import shutil
import sys
import string
import re
import dcsflow.vasp as va
import dcsflow.dftb as db
import dcsflow.phonons as ph
import dcsflow.oclimax as ocl
import dcsflow.restart_dftbmd as rd
import dcsflow.restart_collector as rc
import dcsflow.utils as ut

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    Args:
    path: The current directory (use os.getcwd())
    #########################################
    Return:
    file[0]: The structure file in the path
    """
    file = glob.glob(path + "/*.cif") + glob.glob(path + "/POSCAR") + glob.glob(path + "/geo_end.gen")
    
    if len(file) == 0:
        raise FileNotFoundError

    return file[0]

# 1st flow: DFT-FD 
class Dftfd:
    """ Run DFT simulation with VASP, Phonopy and Oclimax to get INS spectrum
        local (boolean): If True, run VASP on local machine, otherwise run on HPC 
        Args:
        gpu (boolean): If True, run VASP with GPU, otherwise run on CPU 
        hpc (list): srun setting for VASP simulation, Defaults: srun -n {8} -c {32} -G {8} --cpu-bind=cores --gpu-bind=none vasp_std'
        kpts (list): number of KPOINTS. (for Phonon simulation, Defaults to [1,1,1], need user input for relaxation.)
        disp (boolean): If True, set ivdw=12, which is DFT-D3 method with Becke-Johnson damping function (Defaults to True)
        fmax: Maximum allowed force for convergence between atoms (Defaults to 1e-2)
        mesh (list): Mesh size for Phonopy eigenvector problems. (Defaults to [8,8,8])
        supercell: Supercell size. (Defaults to [2,2,2])
    """
    def __init__(self, account, local, gpu, hpc, time, kpts, disp, fmax, mesh, supercell):
        self.account = account
        self.local = local
        self.gpu = gpu
        self.hpc = hpc 
        self.time = time
        self.kpts = kpts   
        self.disp = disp
        self.fmax = fmax
        self.mesh = mesh
        self.supercell = supercell
        self.main_path = os.getcwd()
    
    def relax(self):
        """ Optimization using VASP DFT
        """
        if not os.path.exists(os.path.join(self.main_path, "1-relax" ,"CONTCAR")): # Check whether the relaxation simulation is finished. 
            print(" Relaxation! ")
            os.makedirs(os.path.join(self.main_path, "1-relax"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "1-relax")) # Change directory to 1-relax 
            try:
                atom = ase.io.read(getGeometry(self.main_path))
            except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.") 
                sys.exit(1)  # Exit the script with an error
            va.dft(self.disp, self.fmax, self.kpts, atoms=atom, mode=1)

            if self.local:
                subprocess.run('vasp_std > relax.out', shell=True) # Run VASP on local machine

            if self.gpu:
                slurm_script = ["#!/bin/bash\n",
                f"#SBATCH -A {self.account}\n", 
                "#SBATCH -C gpu\n",
                "#SBATCH -q regular\n",
                "#SBATCH -N 2\n",
                f"#SBATCH -t {self.time}\n",
                "#SBATCH -J vasp_relax\n",
                "#SBATCH -o vasp_relax-%j.out\n",
                "#SBATCH -e vasp_relax-%j.err\n",
                "\n",
                "module load vasp/6.4.3-gpu\n",
                "\n",
                "export OMP_NUM_THREADS=1\n",
                "export OMP_PLACES=threads\n",
                "export OMP_PROC_BIND=spread\n",
                "\n",
                f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std"
                ]
                slurm_filename = "vasp_relax.slurm"
                with open(slurm_filename, "w") as f:
                    f.writelines(slurm_script)
                subprocess.run(f"sbatch {slurm_filename}", shell=True)
                print("VASP relaxation job submitted!")

            else:
                slurm_script = ["#!/bin/bash\n",
                f"#SBATCH -A {self.account}\n", 
                "#SBATCH -C cpu\n",
                "#SBATCH -q regular\n",
                "#SBATCH -N 2\n",
                f"#SBATCH -t {self.time}\n",
                "#SBATCH -J vasp_phonon\n",
                "#SBATCH -o vasp_phonon-%j.out\n",
                "#SBATCH -e vasp_phonon-%j.err\n",
                "\n",
                "module load vasp/6.4.3-cpu\n",
                "\n",
                "export OMP_NUM_THREADS=2\n",
                "export OMP_PLACES=threads\n",
                "export OMP_PROC_BIND=spread\n",
                "\n",
                f"srun -n 128 -c 4 --cpu-bind=cores vasp_std" 
                ]
                slurm_filename = "vasp_relax.slurm"
                with open(slurm_filename, "w") as f:
                    f.writelines(slurm_script)
                subprocess.run(f"sbatch {slurm_filename}", shell=True)
                print("VASP relaxation job submitted!") # Run VASP on hpc cpu nodes
            os.chdir(self.main_path)
            
        else:
            print(" The relaxation is already finished, please continue simulating phonons!!! ")

    def phonons(self):
        """ Phonopy simulation
        """
        ph_kpts = [1,1,1] # The kpoints for phonon simulation
        if not os.path.exists(os.path.join(self.main_path, "2-phonons", "SPOSCAR")): # Check whether the displaced structure is created. 
            print(" Using Phonopy to create displaced structures! ")
            src = os.path.join(self.main_path, "1-relax", "CONTCAR") # source file
            dst = os.path.join(self.main_path, "2-phonons") # destination
            os.makedirs(dst, exist_ok=True)
            os.chdir(dst)
            if not os.path.exists(src):
                ut.print_error(f"Structure file '{src}' is missing. Exiting.")
                sys.exit(1) # Exit the function if the source file is missing
            shutil.copy(src, dst)
        
            ph.displace_structure(calc='vasp', supercell=self.supercell) # Create displaced POSCAR files by Phonopy
            print(" Displaced structures are created! ")
            root, dirs, files = next(os.walk(dst))
        
            for file in files:
                if '-' not in file: # Ensures file name is in correct format
                    continue
                folder_name = file.split('-')[1]  # Extracts number, e.g. '009' from 'POSCAR-009'
                folder_path = os.path.join(dst, folder_name)
                os.makedirs(folder_path, exist_ok=True) # Create the folder if it doesn't exist
    
                source = os.path.join(dst, file) # Define the source and destination file paths
                new_file = "POSCAR"  # Rename the file to 'POSCAR'
                destination = os.path.join(folder_path, new_file)
                shutil.move(source, destination) # Move the file into the corresponding folder

                os.chdir(folder_path)
                atom = ase.io.read(getGeometry(folder_path))
                va.dft(self.disp, self.fmax, kpts=ph_kpts, atoms=atom, mode=2)

        print(" Single point energy calculation using VASP for each displaced structure! ")    
        root2, dirs2, files2 = next(os.walk(dst)) # Loop through every sub-folder in 2-phonons to write VASP input files and run script

        # Write directories to phonon_dirs.txt for SLURM job array
        phonon_dir_list = os.path.join(self.main_path, 'phonon_dirs.txt')
        with open(phonon_dir_list, 'w') as f:
            for d in dirs2:
                f.write(f"{d}\n")
                
        # Generate SLURM job array script
        node_type = "gpu" if self.gpu else "cpu"
        vasp_module = f"vasp/6.4.3-{node_type}"
        omp_threads = 1 if self.gpu else 2
        cpu_bind = "--cpu-bind=cores"
        gpu_bind = "--gpu-bind=none" if self.gpu else ""
        srun_line = (
            f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} {cpu_bind} {gpu_bind} vasp_std"
            if self.gpu else
            "srun -n 128 -c 4 --cpu-bind=cores vasp_std"
        )

        slurm_script = [
            "#!/bin/bash\n",
            f"#SBATCH -A {self.account}\n",
            f"#SBATCH -C {node_type}\n",
            "#SBATCH -q regular\n",
            "#SBATCH -N 2\n",
            f"#SBATCH -t {self.time}\n",
            "#SBATCH -J vasp_ph_array\n",
            "#SBATCH -o vasp_ph_array-%A_%a.out\n",
            "#SBATCH -e vasp_ph_array-%A_%a.err\n",
            f"#SBATCH --array=0-{len(dirs2)-1}\n",
            "\n",
            f"module load {vasp_module}\n",
            f"export OMP_NUM_THREADS={omp_threads}\n",
            "export OMP_PLACES=threads\n",
            "export OMP_PROC_BIND=spread\n",
            "\n",
            f"DIR=$(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" {phonon_dir_list})\n",
            f"cd {os.path.join(self.main_path, '2-phonons')}/$DIR || exit 1\n",
            f"{srun_line.strip()}\n"
        ]

        # Submit SLURM job array
        slurm_filename = os.path.join(self.main_path, "vasp_phonon_array.slurm")
        with open(slurm_filename, "w") as f:
            f.writelines(slurm_script)

        subprocess.run(f"sbatch {slurm_filename}", shell=True)
        print("VASP phonon job array submitted!")
         
    def oclimax(self):
        """ oclimax simulation
        """
        if not os.path.exists(os.path.join(self.main_path, "2-phonons", "FORCE_SETS")):
            os.chdir(os.path.join(self.main_path, '2-phonons')) # Phonopy mesh simulation
            ph.get_force_sets(calc='vasp')
            ph.run_mesh(self.mesh)
            os.chdir(self.main_path)

        if not os.path.exists(os.path.join(self.main_path, "3-oclimax", "ocl.out")):
            ocl.oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=1)

# 2nd flow: DFTB-FD

class Dftbfd:
    """ Run DFTB+, Phonopy and Oclimax to get INS spectrum
    kpts (list): number of KPOINTS. (for Phonon simulation, Defaults to [1,1,1], need user input for relaxation.)
    disp (boolean): If True, which is DFT-D3 method with Becke-Johnson damping function (Defaults to True)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-2)
    mesh (list): Mesh size for Phonopy eigenvector problems. (Defaults to [8,8,8])
    supercell: Supercell size. (Defaults to [2,2,2])
    """
    def __init__(self, kpts, disp, fmax, mesh, supercell):     
        self.kpts = kpts
        self.disp = disp
        self.fmax = fmax
        self.mesh = mesh
        self.supercell = supercell
        self.main_path = os.getcwd()

    def relax(self):
        """ Relaxation using DFTB+
        """
        if not os.path.exists(os.path.join(self.main_path, "1-relax", "geo_end.gen")): # Check whether the relaxation simulation is finished. 
            print(" Relaxation! ")
            os.makedirs(os.path.join(self.main_path, "1-relax"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "1-relax")) # Change directory to 1-relax
            try:
                atom = ase.io.read(getGeometry(self.main_path))
            except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.") 
                sys.exit(1)  # Exit the script with an error
            db.relax(self.disp,self.fmax,self.kpts,atoms=atom)
            subprocess.run('ulimit -s unlimited', shell=True)
            subprocess.run('dftb+ > relax.out', shell=True)

        else:
            print(" The Relaxation is finished, skip to phonon simulation! ")

    def phonons(self):
        """ Phonopy simulation
        """
        ph_kpts = [1,1,1] # The kpoints for phonon simulation
        if not os.path.exists(os.path.join(self.main_path, "2-phonons", "mesh.yaml")): # Check whether the phonon simulation is finished. 
            print(" Using Phonopy to create displaced structures! ")
            src = os.path.join(self.main_path, '1-relax', 'geo_end.gen') # source file
            dst = os.path.join(self.main_path, '2-phonons') # destination
            os.makedirs(dst, exist_ok=True)
            os.chdir(dst)
            if not os.path.exists(src):
                ut.print_error(f"Structure file '{src}' is missing. Exiting.")
                sys.exit(1) # Exit the function if the source file is missing
            shutil.copy(src, dst)
        
            ph.displace_structure(calc='dftbp', supercell=self.supercell) # Create displaced geo.genS files
            print(" Displaced structures are created! ")
            root, dirs, files = next(os.walk(dst))
        
            for file in files:
                if '-' not in file:
                    continue
                folder_name = file.split('-')[1]  # Extracts '009' from 'geo.genS-009'
                folder_path = os.path.join(dst, folder_name)
                os.makedirs(folder_path, exist_ok=True) # Create the folder if it doesn't exist
    
                source = os.path.join(dst, file) # Define the source and destination file paths
                new_file = "geo_end.gen"  # Rename the file to 'geo_end.gen'
                destination = os.path.join(folder_path, new_file)
                shutil.move(source, destination) # Move the file into the corresponding folder

            print(" Single point energy calculation using DFTB+ for each displaced structure! ")
            root2, dirs2, files2 = next(os.walk(dst))
            for d in dirs2:
                if len(d) <= 3:
                    geo = os.path.join(self.main_path, '2-phonons', d)
                    os.chdir(geo)
                    if not os.path.exists(geo+'results.tag'):
                        atom = ase.io.read(getGeometry(geo))
                        db.force(self.disp, kpts=ph_kpts, atoms=atom)
                        subprocess.run('dftb+ > force.out', shell=True)
                    os.chdir('..')
        
        ph.get_force_sets(calc='dftbp')
        ph.run_mesh(self.mesh)

    def oclimax(self):
        """ oclimax simulation
        """
        if not os.path.exists(os.path.join(self.main_path, "3-oclimax", "ocl.out")):
            ocl.oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=1)

# 3rd flow: DFT-MD

class Dftmd:
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

    def __init__(self, account, local, gpu, hpc, time, kpts, disp, fmax, nsw1, nsw2, steps, supercell, tebeg, teend):
        self.account = account
        self.local = local
        self.gpu = gpu
        self.hpc = hpc
        self.time = time
        self.kpts = kpts
        self.disp = disp
        self.fmax = fmax
        self.nsw1 = nsw1
        self.nsw2 = nsw2
        self.steps = steps
        self.supercell = supercell
        self.tebeg = tebeg
        self.teend = teend
        self.main_path = os.getcwd()
    
    def relax(self):
        """ Relaxation using VASP DFT
        """
        if not os.path.exists(os.path.join(self.main_path, "1-relax", "CONTCAR")): # Check whether the relaxation simulation is finished. 
            print(" Relaxation! ")
            os.makedirs(os.path.join(self.main_path, "1-relax"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "1-relax")) # Change directory to 1-relax
            try:
                atom = ase.io.read(getGeometry(self.main_path))
            except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.") 
                sys.exit(1)  # Exit the script with an error
            va.dft(self.disp, self.fmax, self.kpts, atoms=atom, mode='relax')

            if local:
                subprocess.run('vasp_std > relax.out', shell=True) # Run VASP on local machine

            if self.gpu:
                slurm_script = ["#!/bin/bash\n",
                f"#SBATCH -A {self.account}\n", 
                "#SBATCH -C gpu\n",
                "#SBATCH -q regular\n",
                "#SBATCH -N 2\n",
                f"#SBATCH -t {self.time}\n",
                "#SBATCH -J vasp_relax\n",
                "#SBATCH -o vasp_relax-%j.out\n",
                "#SBATCH -e vasp_relax-%j.err\n",
                "\n",
                "module load vasp/6.4.3-gpu\n",
                "\n",
                "export OMP_NUM_THREADS=1\n",
                "export OMP_PLACES=threads\n",
                "export OMP_PROC_BIND=spread\n",
                "\n",
                f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std" 
                ]
                slurm_filename = "vasp_relax.slurm"
                with open(slurm_filename, "w") as f:
                    f.writelines(slurm_script)
                subprocess.run(f"sbatch {slurm_filename}", shell=True)
                print("VASP relaxation job submitted!")
            
            else:
                subprocess.run(f'srun -n {self.hpc[0]} -c {self.hpc[1]} --cpu-bind=cores vasp_st > relax.out', shell=True) # Run VASP on hpc cpu nodes
        
            os.chdir(self.main_path)
            
        else:
            print(" The relaxation is already finished, Continue to NVT-MD !!! ")

    def nvtmd(self):
        """ NVT-MD simulation
        """
        md_kpts = [1,1,1] # The kpoints for MD simulation
        if not os.path.exists(os.path.join(self.main_path, "2-nvtmd", "vasprun.xml")):
            print(" NVT-MD simulation! ")
            os.makedirs(os.path.join(self.main_path, "2-nvtmd"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "2-nvtmd"))
            atom = ase.io.read(os.path.join(self.main_path,'1-relax','CONTCAR'))
            atom *= self.supercell
            va.md(self.disp,self.nsw1,self.tebeg,self.teend,kpts=md_kpts,atoms=atom,ensemble=1)

            if local:
                subprocess.run('vasp_std > nvtmd.out', shell=True) # Run VASP on local machine

            if self.gpu:
                slurm_script = ["#!/bin/bash\n",
                f"#SBATCH -A {self.account}\n",
                "#SBATCH -C gpu\n",
                "#SBATCH -q regular\n",
                "#SBATCH -N 2\n",
                f"#SBATCH -t {self.time}\n",
                "#SBATCH -J vasp_nvtmd\n",
                "#SBATCH -o vasp_nvtmd-%j.out\n",
                "#SBATCH -e vasp_nvtmd-%j.err\n",
                "\n",
                "module load vasp/6.4.3-gpu\n",
                "\n",
                "export OMP_NUM_THREADS=1\n",
                "export OMP_PLACES=threads\n",
                "export OMP_PROC_BIND=spread\n",
                "\n",
                f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std"
                ]
                slurm_filename = "vasp_nvtmd.slurm"
                with open(slurm_filename, "w") as f:
                    f.writelines(slurm_script)
                subprocess.run(f"sbatch {slurm_filename}", shell=True)
                print("VASP NVTMD job submitted!")

            else:
                subprocess.run(f'srun -n {self.hpc[0]} -c {self.hpc[1]} --cpu-bind=cores vasp_std > nvtmd.out', shell=True) # Run VASP on hpc cpu nodes

        else:
            print(" The thermalization (NVT-MD) is already finished, Continue to NVE-MD !!! ")

    def nvemd(self):
        """ NVE-MD simulation
        """
        md_kpts = [1,1,1]
        if not os.path.exists(os.path.join(self.main_path, "3-nvemd", f"{self.steps}", "vasprun.xml")):
            print(" NVE-MD simulation! ")
            os.makedirs(os.path.join(self.main_path, "3-nvemd"), exist_ok=True)
            if not os.path.exists(os.path.join(self.main_path, "3-nvemd", "1", "vasprun.xml")):
                atom = ase.io.read(os.path.join(self.main_path,'2-nvtmd','CONTCAR'))
                os.makedirs(os.path.join(self.main_path, "3-nvemd", "1"), exist_ok=True)
                os.chdir(os.path.join(self.main_path, "3-nvemd", "1"))
                va.md(self.disp,self.nsw2,self.tebeg,self.teend,kpts=md_kpts,atoms=atom,ensemble=2)
                
                if local:
                    subprocess.run('vasp_std > nvemd.out', shell=True) # Run VASP on local machine

                else:
                    if self.gpu:
                        slurm_script = ["#!/bin/bash\n",
                        f"#SBATCH -A {self.account}\n",
                        "#SBATCH -C gpu\n",
                        "#SBATCH -q regular\n",
                        "#SBATCH -N 2\n",
                        f"#SBATCH -t {self.time}\n",
                        "#SBATCH -J vasp_nvemd\n",
                        "#SBATCH -o vasp_nvemd-%j.out\n",
                        "#SBATCH -e vasp_nvemd-%j.err\n",
                        "\n",
                        "module load vasp/6.4.3-gpu\n",
                        "\n",
                        "export OMP_NUM_THREADS=1\n",
                        "export OMP_PLACES=threads\n",
                        "export OMP_PROC_BIND=spread\n",
                        "\n",
                        f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std" 
                        ]
                        slurm_filename = "vasp_nvemd.slurm"
                        with open(slurm_filename, "w") as f:
                            f.writelines(slurm_script)
                        subprocess.run(f"sbatch {slurm_filename}", shell=True)
                        print(" Folder 1 VASP NVEMD simulation submitted! ")

                    else: # CPU
                        subprocess.run(f'srun -n {self.hpc[0]} -c {self.hpc[1]} --cpu-bind=cores vasp_std > nvemd.out', shell=True) # Run VASP on hpc cpu nodes

                os.chdir('..')

            for i in range(1,self.steps):
                atom = ase.io.read('3-nvemd/'+str(i)+'/CONTCAR')
                os.makedirs(os.path.join(self.main_path, "3-nvemd", str(i+1)), exist_ok=True)
                os.chdir(os.path.join(self.main_path, "3-nvemd", str(i+1)))
                if not os.path.exists(os.path.join(self.main_path, "3-nvemd", f"{i}", "vasprun.xml")):
                    print(f" Running folder {i} MD simulation! ")
                    va.md(self.disp,self.nsw2,self.temp,kpts=md_kpts,atoms=atom,ensemble=2)

                    if local:
                        subprocess.run('vasp_std > nvemd.out', shell=True) # Run VASP on local machine

                    else:
                        if self.gpu:
                            slurm_script = ["#!/bin/bash\n",
                            f"#SBATCH -A {self.account}\n",
                            "#SBATCH -C gpu\n",
                            "#SBATCH -q regular\n",
                            "#SBATCH -N 2\n",
                            f"#SBATCH -t {self.time}\n",
                            "#SBATCH -J vasp_nvemd\n",
                            "#SBATCH -o vasp_nvemd-%j.out\n",
                            "#SBATCH -e vasp_nvemd-%j.err\n",
                            "\n",
                            "module load vasp/6.4.3-gpu\n",
                            "\n",
                            "export OMP_NUM_THREADS=1\n",
                            "export OMP_PLACES=threads\n",
                            "export OMP_PROC_BIND=spread\n",
                            "\n",
                            f"srun -n {self.hpc[0]} -c {self.hpc[1]} -G {self.hpc[2]} --cpu-bind=cores --gpu-bind=none vasp_std" 
                            ]
                            slurm_filename = "vasp_nvemd.slurm"
                            with open(slurm_filename, "w") as f:
                                f.writelines(slurm_script)
                            subprocess.run(f"sbatch {slurm_filename}", shell=True)
                            print(f" Folder {i} VASP NVEMD job submitted!")

                        else:
                            subprocess.run(f'srun -n {self.hpc[0]} -c {self.hpc[1]} --cpu-bind=cores vasp_std > nvemd.out', shell=True) # Run VASP on hpc cpu nodes
                os.chdir('..')

        else:
            print(" The (NVE-MD) is already finished, Continue to Oclimax !!! ")
    
    def oclimax(self):
        """ oclimax simulation
        """
        if not os.path.exists(os.path.join(self.main_path, "4-oclimax", "ocl.out")):
            ocl.oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=2)

# 4th flow: DFTB-MD

class Dftbmd:
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

    def __init__(self, kpts, disp, fmax, nsw1, nsw2, steps, supercell, temp):
        self.kpts = kpts    
        self.disp = disp
        self.fmax = fmax
        self.nsw1 = nsw1
        self.nsw2 = nsw2
        self.steps = steps
        self.supercell = supercell
        self.temp = temp
        self.main_path = os.getcwd()
    
    def relax(self):
        """ Relaxation using DFTB+
        """
        if not os.path.exists(os.path.join(self.main_path, "1-relax", "geo_end.gen")): # Check whether the relaxation simulation is finished. 
            print(" Relaxation! ")
            os.makedirs(os.path.join(self.main_path, "1-relax"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "1-relax")) # Change directory to 1-relax
            try:
                atom = ase.io.read(getGeometry(self.main_path))
            except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.")
                sys.exit(1)  # Exit the script with an error
            db.relax(self.disp,self.fmax,self.kpts,atoms=atom)
            subprocess.run('ulimit -s unlimited', shell=True)
            subprocess.run('dftb+ > relax.out', shell=True)

        else:
            print(" The Relaxation is finished, skip to the phonon simulation! ")
    
    def nvtmd(self):
        """ NVT-MD simulation
        """
        md_kpts = [1,1,1]
        if not os.path.exists(os.path.join(self.main_path, "2-nvtmd", "nvtmd.out")): # Check whether the nvtmd simulation is finished.
            print(" NVT-MD simulation! ")
            os.makedirs(os.path.join(self.main_path, "2-nvtmd"), exist_ok=True)
            os.chdir(os.path.join(self.main_path, "2-nvtmd"))
            atom = ase.io.read(os.path.join(self.main_path, '1-relax', 'geo_end.gen'))
            atom *= self.supercell
            db.md(self.disp,self.nsw1,self.temp,kpts=md_kpts,atoms=atom,ensemble=1)
            subprocess.run('dftb+ > nvtmd.out', shell=True)
            rd.make_files(max_iter=0,extra_files=[],output_dir=None,
                        self_copy=False,write_over=False,force_restart=False,
                        restart_from=-1)
                              
        else:
            print(" The NVT-MD simulation is already finished!!! ")

    def nvemd(self):
        """ NVE-MD simulation
        """ 
        md_kpts = [1,1,1]
        letters = list(string.ascii_lowercase) # Naming the folders with letters
        if not os.path.exists(os.path.join(self.main_path, "3-nvemd", f"{letters[self.steps]}", "nvemd.out")): # Check whether the nvemd simulation is finished.
            print(" NVE-MD simulation! ")
            os.makedirs(os.path.join(self.main_path, "3-nvemd"), exist_ok=True) # Here we will create subfolders for NVE-MD simulation.
            if not os.path.exists(os.path.join(self.main_path, "3-nvemd", "a", "nvemd.out")): # Check whether the "a" simulation is finished.
                src = os.path.join(self.main_path, "2-nvtmd", "restart") # source folder
                dst = os.path.join(self.main_path, "3-nvemd", "a") # destination folder
                shutil.copytree(src, dst)
                os.chdir(os.path.join(self.main_path, "3-nvemd", "a"))
                with open("dftb_in.hsd", "r") as f1:
                    data = f1.read()
            
                data = re.sub(
                        r"Thermostat\s*\{[^{}]*\{[\s\S]*?\}[\s\S]*?\}",  # Match entire "Thermostat { ... }" with nested braces
                        "Thermostat {None{}",  # Replace with an empty Thermostat block, this is for NVE-MD simulation
                        data,
                        flags=re.MULTILINE)
                data = re.sub(
                        r"Steps\s*=\s*\d+",  # Match "Steps = <number>"
                        f"Steps = {self.nsw2}",  # Replace with the desired number
                        data)
                with open("dftb_in.hsd", "w") as f2:
                    f2.write(data)

                subprocess.run('dftb+ > nvemd.out', shell=True)
                rd.make_files(max_iter=0,extra_files=[],output_dir=None,
                            self_copy=False,write_over=False,force_restart=False,
                            restart_from=-1)
                os.chdir('..')

            for i in range(0,self.steps-1): # We separate the NVE-MD simulation into multiple folders
                if not os.path.exists(os.path.join(self.main_path, "3-nvemd", letters[i+1], "nvemd.out")):
                    print(f" Running folder {letters[i+1]} MD simulation! ")
                    src_ = os.path.join(self.main_path, "3-nvemd", letters[i], 'restart') # source folder
                    dst_ = os.path.join(self.main_path, "3-nvemd", letters[i+1]) # destination folder
                    shutil.copytree(src_, dst_)
                    os.chdir(os.path.join(self.main_path, "3-nvemd", letters[i+1]))
                    subprocess.run('dftb+ > nvemd.out', shell=True)
                    rd.make_files(max_iter=0,extra_files=[],output_dir=None,
                            self_copy=False,write_over=False,force_restart=False,
                            restart_from=-1)
                    os.chdir('..')

            # Have to collect the trajectory files and put it together.
            rc.collect(steps=self.steps)

        else:
            print(" The NVE-MD simulation is already finished!!! ")
    
    def oclimax(self):
        """ oclimax simulation
        """
        if not os.path.exists(os.path.join(self.main_path, "4-oclimax", "ocl.out")):
            ocl.oclimax(dt=1.0,params=None,task=0,e_unit=0,mode=3)