#!/usr/bin/env python

import mdtraj as md
import numpy as np
import os
import argparse

from utils import preprocess
from utils import recursive_fitting

MAX_WOKERS = 4

def main():
    parser = argparse.ArgumentParser(description='fitting trajectory.')
    parser.add_argument('-t', '--trajectory', required=True, help='trajectory file (.trr)')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-o', '--out', required=True, help='output file path (.trr)')
    args = parser.parse_args()


    ### read file ###
    trj_mdtraj = md.load_trr(args.trajectory, top=args.topology)
    n_frames = trj_mdtraj.n_frames
    print(f'Trajectory Info ({n_frames} frames, {trj_mdtraj.n_atoms} atoms)')


    ### preprocess ###
    trj_mdtraj, atomlist, wlist = preprocess(trj_mdtraj)


    ### fitting ###
    trj_array = recursive_fitting(trj_mdtraj.xyz, wlist, MAX_WOKERS)


    ### ndarray to trr ###
    topology = trj_mdtraj.topology
    trj_mdtraj = md.Trajectory(trj_array, topology)


    ### save ###
    trj_mdtraj.save_trr(args.out)


if __name__=='__main__':
    main()