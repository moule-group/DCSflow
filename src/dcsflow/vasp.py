#### VASP DFT simulation with DCS-Flow ####
# The author: Ray
# Starting date: 01/06/2024 

import os
import argparse
import shutil
import glob
from ase.io import read
import subprocess
from ase.calculators.vasp import Vasp
import time

#def runfile(name):
#    """ Create run file for using HPC (NERSC)
#    Args:
#    name: the name of the run file
    ##########################
#    Return: run.{name} file in the folder for HPC job submission
#    """
#    with open(f'run.{name}', 'w') as f:
#        cmds = ["#!/bin/bash\n",
#                f"#SBATCH -J {name}\n",
#                "#SBATCH -C gpu\n",
#                "#SBATCH -A m2734\n",
#                "#SBATCH -q regular\n",
#                "#SBATCH -t 04:00:00\n",
#                "#SBATCH -N 2\n",
#               "#SBATCH -G 8\n",
#               "#SBATCH --exclusive\n",
#                "#SBATCH -o %x-%j.out\n",
#                "#SBATCH -e %x-%j.err\n",
#                "\n",
#                "module load vasp/6.4.1-gpu\n", # Load VASP module, might have to change the version
#                "\n",
#                "export OMP_NUM_THREADS=1\n",
#                "export OMP_PLACES=threads\n",
#                "export OMP_PROC_BIND=spread\n",
#                "\n",
#                "srun -n 8 -c 32 --cpu-bind=cores --gpu-bind=none -G 8 vasp_std"]
#        f.writelines(cmds)
    
def dft(disp, fmax, kpts, atoms, mode):
    """ Create VASP DFT input files
    Args:
    disp (int): Set ivdw in INCAR, 12: DFT-D3 method with Becke-Johnson damping function 
    fmax (float): Maximum allowed force for convergence between atoms 
    kpts (list): number of kpts for DFT simulation. 
    atoms (str): ASE atoms object.
    mode (int): 1. relax 2. phonons
    """
    if mode == 1:
        calc = Vasp(prec="Accurate",
                    xc="pbe",
                    encut=520.000000,
                    ismear=0,
                    sigma=0.050000,
                    ediff=1.00e-06,
                    ediffg=(-1*fmax),
                    nsw=1000,
                    ivdw=disp,
                    isif=3,
                    ibrion=2,
                    kpts=kpts,
                    ncore=32,
                    lreal=False,
                    lcharg=False,
                    lwave=False,
                    nwrite=1,
                    gamma=True)
             
    if mode == 2:
        calc = Vasp(prec="Accurate",
                    xc="pbe",
                    encut=520.000000,
                    ismear=0,
                    sigma=0.050000,
                    ediff=1.00e-06,
                    ediffg=(-1*fmax),
                    nsw=0,
                    ivdw=disp,
                    isif=3,
                    ibrion=-1,
                    kpts=kpts,
                    ncore=32,
                    lreal=False,
                    lcharg=False,
                    lwave=False,
                    nwrite=1,
                    gamma=True)

    atoms.set_calculator(calc) 
    calc.write_input(atoms)
    
def md(disp, nsw, tebeg, teend, kpts, atoms, ensemble):
    """ Create VASP MD input files
    Args:
    disp (int): Set ivdw in INCAR, 12: DFT-D3 method with Becke-Johnson damping function 
    nsw (int): Number of MD steps
    tebeg (float): Beginning temp (K)
    teend (float): Ending temp (K)
    kpts (list): number of KPOINTS
    atoms (str): ASE atoms object
    ensemble (int): 1: NVT; 2:NVE 3:NPT
    """ 
    if ensemble == 1: # Nose-Hoover NVT
        calc = Vasp(prec='Accurate',
                    xc='pbe',
                    ismear=0,
                    sigma=0.050000,
                    encut=520.000000,
                    isym=0,
                    ibrion=0, # MD simulation
                    isif=2, # fix cell volume
                    nsw=nsw, # number of steps
                    ivdw=disp,
                    ediff=1.00e-06,
                    ediffg=-1.00e-2,
                    mdalgo=2, # Nose Hoover
                    potim=0.25, # timestpe size
                    smass=0, # the NH mass controls frequency of temperature oscillations
                    tebeg=tebeg,
                    kpts=kpts,
                    ncore=32,
                    lreal=False,
                    lcharg=False,
                    lwave=False,
                    nwrite=1,
                    gamma=True)
             
        atoms.set_calculator(calc) 
        calc.write_input(atoms)
        
    if ensemble == 2: # NVE
        calc = Vasp(prec='Accurate',
                    xc='pbe',
                    ismear=0,
                    sigma=0.050000,
                    encut=520.000000,
                    isym=0,
                    ibrion=0, # MD simulation
                    isif=2, # fix cell volume
                    nsw=nsw, # number of steps
                    ivdw=disp,
                    ediff=1.00e-06,
                    ediffg=-1.00e-2,
                    mdalgo=1, # Andersen
                    andersen_prob=0, # for NVE the additional setting 
                    potim=0.25, # timestep size
                    tebeg=tebeg,
                    kpts=kpts,
                    ncore=32,
                    lreal=False,
                    lcharg=False,
                    lwave=False,
                    nwrite=1,
                    gamma=True)
             
        atoms.set_calculator(calc) 
        calc.write_input(atoms)
        
    if ensemble == 3: # NPT (specific for MLFF training)
        calc = Vasp(prec="Accurate",
                    xc="pbe",
                    encut=520.000000,
                    isym=0,
                    ismear=0,
                    sigma=0.050000,
                    ediff=1.00e-06,
                    ediffg=-1.00e-2,
                    kpts=kpts,
                    nsw=nsw,
                    ivdw=disp,
                    ibrion=0,
                    potim=0.25,
                    tebeg=tebeg,
                    teend=teend,
                    mdalgo=3,
                    isif=3,
                    ncore=32,
                    lreal=False,
                    lcharg=False,
                    lwave=False,
                    nwrite=1,
                    gamma=True)

        atoms.set_calculator(calc) 
        calc.write_input(atoms)
