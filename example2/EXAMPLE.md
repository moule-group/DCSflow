# Example of DCSflow: workflow 2

This example will show you the how to use workflow 2 (DFTB-FD)

## Usage

1. HPC (NERSC as example here)

Use interactive jobs

GPU
```
salloc --nodes 2 --qos interactive --time 01:00:00 --constraint gpu --gpus 8 --account mxxxx
```

CPU 
```
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account mxxxx
```

1. Relaxation, need specify account number 

```
dcsflow -w 2 -r 
```

2. Phonon
```
dcsflow -w 2 -p 
```

3. Oclimax
```
dcsflow -w 2 -o
```
