
# Example of DCSflow: workflow 4

This example will show you the how to use workflow 4 (DFTB-MD)

## Usage

We will show the example running on HPC (NERSC as example here)

1. Relaxation, need specify account number 

```
dcsflow -w 4 -r 
```

2. NVT-MD, this step is to thermalize the system.
```
dcsflow -w 4 -nvt 
```

3. NVE-MD, produce trajectory
```
dcsflow -w 4 -nve
```

4. Oclimax
```
dcsflow -w 4 -o
```
