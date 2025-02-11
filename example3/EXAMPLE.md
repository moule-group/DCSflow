
# Example of DCSflow: workflow 3

This example will show you the how to use workflow 3 (DFT-MD)

## Usage

We will show the example running on HPC (NERSC as example here)

1. Relaxation, need specify account number 

```
dcsflow -w 3 -r -a mxxxx
```

2. NVT-MD, this step is to thermalize the system.
```
dcsflow -w 3 -nvt 
```

3. NVE-MD, produce trajectory. Default: 4000 steps in each loop (-n2) and it will loop 10 times (-step)
```
dcsflow -w 3 -nve
```

4. Oclimax
```
dcsflow -w 3 -o
```

