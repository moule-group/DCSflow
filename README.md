# Davis Computational Spectroscopy Workflow 

DCS-flow is a tool to simulate inelastic neutron scattering (INS) spectrum, it is still under development.


## Installation

1. DCSflow

```

```

2. DFTB+: 

```
    mamba install 'dftbplus=*=nompi_*'
```
Download Slater-Koster files from DFTB website https://dftb.org/parameters/download/all-sk-files .

For more detail, please refer to 

3. Oclimax: Please refer to the page: https://sites.google.com/site/ornliceman/download .

## Description

There are 2 methods simulating phonon spectrum and convert into inelastic neutron scattering spectrum

## Documentation

1. 2 Methods for simulating INS spectrum. One is finite displacement method (FD) and the other is molecular dynamics method (MD). The theory of these 2 methods are written in next section. 

2. There are 2 calculators in DCS-flow, one is VASP and the other is DFTB+.

3. Commands table

|  Command |  DFT-FD  |  DFTB-FD |  DFT-MD  | DFTB-MD  |   
|----------|----------|----------|----------|----------|
|    -w    |    1     |     2    |     3    |     4    | 
|    -r    |   False  |    False |   False  |   False  |
|    -p    | False    |    False |   N/A    |   N/A    | 
|    -nvt  |   N/A    |    N/A   |   False  |   False  |
|    -nve  |   N/A    |    N/A   |   False  |   False  |
|    -o    |   False  |    False |   False  |   False  |
|    -l    |    f     |    N/A   |     f    |   N/A    |  
|    -g    |    t     |    N/A   |     t    |    N/A   | 
|    -H    | [8,32,8] |    N/A   | [8,32,8] |    N/A   |  
|    -f    |   1e-2   |   1e-2   |   1e-2   |   1e-2   |
|    -t1   |   N/A    |   N/A    |    150   |    150   |
|    -t2   |   N/A    |    N/A   |    150   |    N/A   |
|    -k    | [1,1,1]  | [1,1,1]  | [1,1,1]  | [1,1,1]  |
|    -s    | [2,2,2]  | [2,2,2]  | [2,2,2]  | [2,2,2]  |
|    -d    |    12    |    12    |    12    |    12    |
|    -m    | [8,8,8]  | [8,8,8]  |   N/A    |    N/A   |
|    -n1   |    N/A   |   N/A    |   4000   |   4000   |
|    -n2   |    N/A   |   N/A    |   4000   |   4000   |
|    -step |    N/A   |   N/A    |    10    |    10    |

4. Usage:

Simulation using dftb-fd method for relaxation, must specify -w --workflow to 2 and set -r.
``` 
dcsflow -w 2 -r
```

## Theory

1. Finite Displacement Method (Phonopy)

2. Molecular Dynamics Method