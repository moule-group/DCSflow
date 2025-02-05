###################### DFTB tool by ASE ###########################
# Create and modified by Ray started in May 2024 

# Updation 7/20/24: Adding DFT-D3 dispersion correction 
# Updation 7/23/24: 1. Adding Options_WriteChargesAsText='Yes', this can be used for partial charges in MD sim.; 2. Adding parser to set temperature value
# Updation 10/21/24 Adding DFTB-MD function
# Updation 11/5/24 Adding DFTB-MD dispersion function

###################################################################
import os
import argparse
import glob
from ase.io import read, write
import subprocess
import shutil
from ase.calculators.dftb import Dftb

# Constants

k_to_au = 0.316681534524639E-05 # 1 Kelvin
fs_to_au = 0.413413733365614E+02 # 1 fs
wn_to_au = 0.725163330219952E-06 # 1 cm-1
        
#def runfile(name):
   # """ Create run file for using NERSC HPC, this is written for GPU node.
   # Args:
   # name: The name of the run file.
    ##########################
   # Return: run.sh file in the folder
   # """
    
  #  with open(f'run_{name}', 'w') as f:
  #      cmds = ["#!/bin/bash\n",
   #             f"#SBATCH -J {name}\n",
   #             "#SBATCH -C gpu\n",
   #             "#SBATCH -A m2734\n",
   #             "#SBATCH -q regular\n",
   #             "#SBATCH -t 04:00:00\n",
   #             "#SBATCH -N 1\n",
   #             "#SBATCH -G 4\n",
   #             "#SBATCH --exclusive\n",
   #             "#SBATCH -o %x-%j.out\n",
   #             "#SBATCH -e %x-%j.err\n",
   #             "\n",
   #             "ulimit -s unlimited # for DFTB\n"
   #            "\n",
   #             "export OMP_NUM_THREADS=16\n",
   #             "export OMP_PLACES=threads\n",
   #             "export OMP_PROC_BIND=spread\n",
   #             "\n",
    #            f"srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=none dftb+ > {name}.out"]
   #     f.writelines(cmds)

def md(disp,nsw,temp,kpts,atoms,ensemble):
    """ Run MD using DFTB calculatior 
    Args:
    disp (boolean): Dispersion correction using DFT-D3  
    nsw (int): Number of MD steps
    temp (float): Temperature 
    kpts (list): Number of kpts.
    atoms (str): ASE atoms object.
    ensemble (int): 1:NVT; 2:NVE
    """
    if not disp:
        calc_nvt = Dftb(label='nvt', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='NoseHoover',
                        Driver_Thermostat_Temperature=temp*k_to_au,  # 1000K = 0.0031668 au. Here is 150K
                        Driver_Thermostat_CouplingStrength=3000*wn_to_au, # 3000 cm-1
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kpts,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )

        calc_nve = Dftb(label='nve', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='None',
                        Driver_Thermostat_empty='',
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kpts,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )

    else: 
        calc_nvt = Dftb(label='nvt', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='NoseHoover',
                        Driver_Thermostat_Temperature=temp*k_to_au,  # 1000K = 0.0031668 au. Here is 150K
                        Driver_Thermostat_CouplingStrength=3000*wn_to_au, # 3000 cm-1
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kpts,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Hamiltonian_Dispersion_='DftD3',
                        Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                        Hamiltonian_Dispersion_Damping_a1='0.5719',
                        Hamiltonian_Dispersion_Damping_a2='3.6017',
                        Hamiltonian_Dispersion_s6='1.0',
                        Hamiltonian_Dispersion_s8='0.5883',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )
        
        calc_nve = Dftb(label='nve', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='None',
                        Driver_Thermostat_empty='',
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kpts,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Hamiltonian_Dispersion_='DftD3',
                        Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                        Hamiltonian_Dispersion_Damping_a1='0.5719',
                        Hamiltonian_Dispersion_Damping_a2='3.6017',
                        Hamiltonian_Dispersion_s6='1.0',
                        Hamiltonian_Dispersion_s8='0.5883',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )
        
    if ensemble == 1:
        atoms.calc = calc_nvt
        calc_nvt.write_input(atoms) 
     
    elif ensemble == 2:
        atoms.calc = calc_nve
        calc_nve.write_input(atoms) 

def relax(disp,fmax,kpts,atoms):
    """ Relax structure by using DFTB
    Documentation: 
    Driver: determines how the geometry should be changed 
    during the calculation.
    Driver_Moveatoms="1:-1": Move all atoms in the system
    Hamiltonian_MaxAngularMomentum_C='p':
    the value of the highest orbital angular momentum each element
    ##############################################
    Args:
    kpts (list): number of kpoints for relaxation. (Defaults to [1,1,1], it should be increased.)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-3)
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    atoms (str): ASE atoms object.
    """
    if not disp:
        calc = Dftb(label='relax', 
                Driver_='ConjugateGradient', # relax algorithm (This command will be change in later version)
                Driver_MovedAtoms='1:-1',
                Driver_MaxForceComponent=fmax, # converge if forces on all atoms < fmax
                Driver_MaxSteps=1000,
                Driver_LatticeOpt='Yes',
                Driver_Isotropic='Yes',
                kpts=kpts,
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=300), # Set this temp = 300K
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Options_='',
                Options_WriteChargesAsText='Yes'
                )
    
    else:
        calc = Dftb(label='relax', 
                Driver_='ConjugateGradient', # relax algorithm (This command will be change in later version)
                Driver_MovedAtoms='1:-1',
                Driver_MaxForceComponent=fmax, # converge if forces on all atoms < fmax
                Driver_MaxSteps=1000,
                Driver_LatticeOpt='Yes',
                Driver_Isotropic='Yes',
                kpts=kpts,
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=300), # Set this temp = 300K
                Hamiltonian_MaxAngularMomentum_='', 
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Hamiltonian_Dispersion_='DftD3',
                Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                Hamiltonian_Dispersion_Damping_a1='0.5719',
                Hamiltonian_Dispersion_Damping_a2='3.6017',
                Hamiltonian_Dispersion_s6='1.0',
                Hamiltonian_Dispersion_s8='0.5883',
                Options_='',
                Options_WriteChargesAsText='Yes'
                )

    #atoms.set_calculator(calc) 
    #calc.calculate(atoms)
    #with open('dftb_in.hsd', 'w') as f:
    #calc.write_dftb_in(f) 
    calc.write_input(atoms) 

def force(disp, kpts, atoms):
    """ Single point energy calculation and extract force by using DFTB,
    this step is for getting phonons with phonopy.
    Args:   
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    kpts (list): kpoints list for single point energy simulation. (Defaults to [1,1,1])
    atoms (str): ASE atoms object.
    """
    if not disp:
        calc = Dftb(label='Force', 
                kpts=kpts,
                Analysis_='',
                Analysis_CalculateForces='Yes', # This command changed in dftb+ 24.1
                Options_="",
                Options_WriteResultsTag='Yes',
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=300),
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p')
    else:
        calc = Dftb(label='Force', 
                kpts=kpts,
                Analysis_='',
                Analysis_CalculateForces='Yes', # This command changed in dftb+ 24.1
                Options_='',
                Options_WriteResultsTag='Yes',
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=300),
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Hamiltonian_Dispersion_='DftD3',
                Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                Hamiltonian_Dispersion_Damping_a1='0.5719',
                Hamiltonian_Dispersion_Damping_a2='3.6017',
                Hamiltonian_Dispersion_s6='1.0',
                Hamiltonian_Dispersion_s8='0.5883')
        
    calc.write_input(atoms) 
    #os.system("dftb+ 1>> forces.out 2>> forces.err") # Run this in order to get result.tag

        