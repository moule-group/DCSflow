# Davis Computational Spectroscopy Workflow 

DCS-flow is a tool to simulate inelastic neutron scattering (INS) spectrum, it is still under development.


## Installation

1. DCSflow

```

```

2. DFTB+: Also need download Slater-Koster files from DFTB website https://dftb.org/parameters/download/all-sk-files .

```
    mamba install 'dftbplus=*=nompi_*'
```


3. Oclimax: Please refer to the page: https://sites.google.com/site/ornliceman/download .

4. Setting up environmental variables in configuration file (.bashrc)

```
export DFTB_PREFIX="your_path/slako/mio/mio-1-1/" # SK files for DFTB
export PATH="your_path/oclimax/bin":$PATH # Oclimax path
export VASP_PP_PATH="your_path" # VASP PP path
```


## Description

There are 2 methods simulating phonon spectrum and convert into inelastic neutron scattering spectrum

## Documentation

1. 2 Methods for simulating INS spectrum. One is finite displacement method (FD) and the other is molecular dynamics method (MD). The theory of these 2 methods are written in next section. Currently, there are 2 calculators in DCS-flow, one is VASP and the other is DFTB+.

2. 

    (a) 1st flow: DFT-FD: Relaxation --> Phonon --> Oclimax

    (b) 2nd flow: DFTB-FD: Relaxation --> Phonon --> Oclimax

    (c) 3rd flow: DFT-MD: Relaxation --> Nvt-md --> Nve-md --> Oclimax

    (d) 4th flow: DFTB-MD: Relaxation --> Nvt-md --> Nve-md --> Oclimax

3. Commands table

|  Command |Explanation| DFT-FD  |  DFTB-FD |  DFT-MD  | DFTB-MD  |   
|----------|----------|----------|----------|----------|----------|
|    -w    | workflow (int) |   1     |     2    |     3    |     4    | 
|    -r    | relaxation (bootlean) |   True/**False**  |    True/**False** |  True/**False**  |   True/**False**  |
|    -p    | phonons (bootlean), apply when workflow = 1 or 2 | True/**False**    |    True/**False** |   N/A    |   N/A    | 
|    -nvt  | NVT-MD (bootlean), apply when workflow = 3 or 4 |   N/A    |    N/A   |   True/**False**  |   True/**False**  |
|    -nve  | NVE-MD (bootlean), apply when workflow = 3 or 4|   N/A    |    N/A   |   True/**False**  |   True/**False**  |
|    -o    | oclimax (bootlean)|   True/**False**  |    True/**False** |   True/**False**  |   True/**False**  |
|    -l    | Running VASP in local machine, apply when workflow = 1 or 3 (bootlean)|   True/**False**     |    N/A   |     True/**False**    |   N/A    |  
|    -g    | Running VASP on HPC GPU mode, apply when workflow = 1 or 3 (bootlean). If Flase, run on CPU mode |  **True**/False   |   N/A   |  **True**/False |  N/A   | 
|    -a    | If -g, set your account # in slurm script (str) | your account #      |    N/A   |  your account #  |    N/A   |
|    -t    | If -g, set time limit on HPC | **01:00:00** | N/A | **01:00:00** | N/A |  
|    -H    | srun setting for VASP simulation, ex: srun -n 8 -c 32 -G 8 --cpu-bind=cores --gpu-bind=none vasp_std"| **8 32 8** |  N/A  | **8 32 8** |  N/A  |  
|    -f    | Converge if forces on all atoms < fmax| Float  **1e-2**   |  Float **1e-2**   |  Float **1e-2**   | Float **1e-2**   |
|    -t1   | Initial temperature for MD simulation (K) |   N/A    |   N/A    |   Float **150.0**   |    Float **150.0**   |
|    -t2   | Final temperature for MD simulation (K), currently not under usage|   N/A    |    N/A   |  Float **150.0**   |    N/A   |
|    -k    | Kpoints for simulation. We recommend increase K-grid in relaxation | **1 1 1**  | **1 1 1**  | **1 1 1**  | **1 1 1**  |
|    -s    | Supercell size | **2 2 2**  | **2 2 2**  | **2 2 2**  | **2 2 2**  |
|    -d    | Apply dispersion correction, refer to vasp docs. For DFTB, only support DFT-D3 (Becke-Johnson) if -d > 0| int **12**  | int **12** | int **12** | int **12** |
|    -m    | The mesh grid running phonopy mesh | **8 8 8**  | **8 8 8**  |   N/A   |   N/A   |
|    -n1   | Number of NVT MD steps, apply when workflow = 3 or 4|  N/A | N/A  | int **4000**   | int **4000** |
|    -n2   | Number of NVE MD steps for each loop, apply when workflow = 3 or 4  |  N/A  |  N/A  | int **4000**  | int **4000** |
|    -step | Multiple of NVE-MD steps; (i.e. 4000*steps is the total steps for NVE-MD) |  N/A  | N/A  | int **10**  | int **10** |

4. Usage:

Specify workflow using -w, and then specify the steps you are running (relaxation: -r; phonon: -p; NVTMD: -nvt; NVEMD: -nve; Oclimax:-o).

Example: DFTB-FD relaxation
``` 
dcsflow -w 2 -r
```

## Theory

1. Finite Displacement Method (Phonopy)

2. Molecular Dynamics Method