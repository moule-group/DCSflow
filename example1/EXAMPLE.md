# Example 1 of DCSflow: workflow 1

This example will show you the how to use workflow 1 (DFT-FD)

## Usage

1. HPC (NERSC as example here)

Use interactive jobs, if VASP 6, we recommend using GPU.

GPU
```
salloc --nodes 2 --qos interactive --time 01:00:00 --constraint gpu --gpus 8 --account mxxxx
```

CPU 
```
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account mxxxx
```

Load vasp module 
```
module load vasp/6.4.2-gpu 
```


