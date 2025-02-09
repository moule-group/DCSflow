import os
import argparse
import itertools
import shutil

main_path = os.getcwd()

def load_iter_range(filename):
    """ Load iteration range from file iter_range.txt
    """
    try:
        with open(filename, "r") as f:
            first = f.readline().strip()
            iter_from =  int(first) if first else None
            second = f.readline().strip()
            iter_until =  int(second) if second else None
            iter_range = {"from": iter_from, "until": iter_until}
    except:
      iter_range = {"from": None, "until": None}
    
    return iter_range

def load_frames(frame_filename, iter_from=0):
    """ Read frame file to give frame list which have modified xyz lines
    frame_filename (str): Frame file name
    iter_from (int): Iteration number to start counting from
    """
    frames = []
    with open(frame_filename, "r") as f:
        lines = []    
        header = f.readline() 
        try:
            natoms = int(header) # find the number of atoms
        except ValueError as e:
            raise ValueError('Expected xyz header but got: {}'.format(e))

        lines.append(header)
        comment = f.readline()
        iter_number = int(comment[comment.find('iter:') + 5:].split()[0]) + 0
        lines.append(f"iter:{iter_number} \n")
        for i in range(natoms):
            lines.append(f.readline())
        frames.append(lines)

def get_iter_from_frame(frame):
    """ Get MD iteration number from given xyz frame index
    Args:
    frame (list): XYZ frame
    """ 
    comment = frame[1]
    iter_number = int(comment[comment.find('iter:') + 5:].split()[0])
    return iter_number

# Main function
def collect(steps):
    """ Main function to collect DFTB-MD trajectory 
    Args:
    steps (int): Number of folders to collect
    """
    current_iter = 0
    collected_frames = []

    root, dirs, files = next(os.walk(os.path.join(main_path, "3-nvemd"))) # Walk through all folders in 3-nvemd
    i = 1
    for d in dirs:
        if len(d) <= 2:
            os.chdir(os.path.join(main_path, "3-nvemd", d))
            iter_from = load_iter_range("iter_range.txt")["from"] # Append frames
            frames = load_frames(os.path.join(main_path, "3-nvemd", d, "geo_end.xyz"), iter_from=iter_from)
            if i < steps:
                collected_frames += frames[:-1] # Because the last frame is the same as the first frame of the next iteration
            else: 
                collected_frames += frames
            i += 1
            current_iter = get_iter_from_frame(frames[-1])
            os.chdir('..')
    
    with open('geo_end.xyz', "w") as f:
        f.writelines(list(itertools.chain.from_iterable(collected_frames)))
