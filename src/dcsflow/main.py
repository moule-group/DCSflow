import dcsflow.workflow as workflow
import argparse

def main():
    """ The main function to run the DCS-Flow 
    """
    print(" ----------------------------------------------------------------- ")
    print(" Starting Davis Computation Spectroscopy Flow (DCS-Flow)  !!! ")
    print(" Made by Moule Group at UC Davis (version: 0.1.0) ")
    print(" #########  ##########   #######                ")
    print(" $       $  $            $        ")
    print(" $      $   $            #######           ")
    print(" $     $    $                  $  ")
    print(" ######     ##########   #######  ")
    print(" ----------------------------------------------------------------- ")

    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-w", "--workflow", type=int, default=1, help="Please choose the workflow you want to run (type the number)") # Add an argument: workflow
    parser.add_argument("-l", "--local", action='store_true', default=False, help="VASP simulation running in local desktops!") # Add an argument: locak
    parser.add_argument("-g", "--gpu", action='store_true', default=True, help="VASP 6 simulation running in HPC GPU nodes!") # Add an argument: gpu
    parser.add_argument("-H", "--hpc", type=int, nargs=3, default=[8,32,8], help="srun setting for VASP simulation, ex: srun -n {8} -c {32} -G {8} --cpu-bind=cores --gpu-bind=none vasp_std") # Add an argument: hpc
    parser.add_argument("-f", "--force", type=float, default=1e-2, help="Converge if forces on all atoms < fmax; Defaults to 1e-2 (DFTB can be set to 1e-3!") # Add an argument: force
    parser.add_argument("-t1", "--temp1", type=float, default=150.0, help="Initial temperature for MD simulation; Defaults to 150K! ") # Add an argument: temp1
    parser.add_argument("-t2", "--temp2", type=float, default=150.0, help="Final temperature for MD simulation; Defaults to 150K! (Only useful when training MLFF using NPT)") # Add an argument: temp2
    parser.add_argument("-k", "--kpts", type=int, nargs=3, default=[1,1,1], help="Kpoints for relaxation. Enter numbers separated by spaces") # Add an argument: kpts
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=[2,2,2], help="Supercell size, Defaults to (2,2,2). Enter numbers separated by spaces:") # Add an argument: supercell   
    parser.add_argument("-d", "--dispersion", type=int, default=12, help="Apply DFT-D3 dispersion?; Defaults to 12 (DFT-D3)!") # Add an argument: dispersion
    parser.add_argument("-m", "--mesh",type=int, nargs=3, default=[8,8,8], help="The mesh grid running phonopy mesh, Defaults to [8,8,8]") # Add an argument: mesh
    parser.add_argument("-n1", "--nsw1", type=int, default=4000, help="Number of NVT MD steps") # Add an argument: nsw1
    parser.add_argument("-n2", "--nsw2", type=int, default=4000, help="Number of NVE MD steps") # Add an argument: nsw2
    parser.add_argument("-step", "--steps", type=int, default=10, help="Multiple of NVE-MD steps; (i.e. 4000*steps is the total steps for NVE-MD) Defaults to 10!") # Add an argument: steps
    args = parser.parse_args() # Parse the argument

    if args.workflow == 1:
        print(" ----------------------------------------------------------------- ")
        print(" 1st flow: DFT-FD: VASP DFT --> Phonopy --> Oclimax ")
        print(" ----------------------------------------------------------------- ")
        workflow.dftfd(local=args.local, gpu=args.gpu, hpc=args.hpc, kpts=args.kpts, disp=args.dispersion, fmax=args.force, mesh=args.mesh, supercell=args.supercell)

    if args.workflow == 2:
        print(" ----------------------------------------------------------------- ")
        print(" 2nd flow: DFTB-FD: DFTB+ --> Phonopy --> Oclimax ")
        print(" ----------------------------------------------------------------- ")
        workflow.dftbfd(kpts=args.kpts, disp=args.dispersion, fmax=args.force, mesh=args.mesh, supercell=args.supercell)
        
    if args.workflow == 3: 
        print(" ----------------------------------------------------------------- ")
        print(" 3rd flow: DFT-MD: VASP DFT-MD --> Oclimax ")  
        print(" ----------------------------------------------------------------- ")
        workflow.dftmd(local=args.local, gpu=args.gpu, hpc=args.hpc, disp=args.dispersion, fmax=args.force, kpts=args.kpts, nsw1=args.nsw1, nsw2=args.nsw2, supercell=args.supercell, tebeg=args.temp1, teend=args.temp2)

    if args.workflow == 4:
        print(" ----------------------------------------------------------------- ")
        print(" 4th flow: DFTB-MD: DFTB+ DFTB-MD --> Oclimax ")
        print(" ----------------------------------------------------------------- ")
        workflow.dftbmd(disp=args.dispersion, fmax=args.force, kpts=args.kpts, nsw1=args.nsw1, nsw2=args.nsw2, steps=args.steps, supercell=args.supercell, temp=args.temp1)