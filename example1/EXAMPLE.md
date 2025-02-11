# Example of DCSflow: workflow 1

This example will show you the how to use workflow 1 (DFT-FD)

## Usage

We will show the example running on HPC (NERSC as example here)

1. Relaxation, need specify account number 

```
dcsflow -w 1 -r -a mxxxx
```

2. Phonon
```
dcsflow -w 1 -p 
```

3. Oclimax
```
dcsflow -w 1 -o
```




