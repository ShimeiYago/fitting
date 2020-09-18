#!/usr/bin/env python

import mdtraj as md
import numpy as np
import os
import argparse

from utils import preprocess
from utils import recursive_fitting

def main():
    parser = argparse.ArgumentParser(description='fitting trajectory.')
    parser.add_argument('-t', '--trajectory', required=True, help='trajectory file (.trr)')
    parser.add_argument('--npy', required=True, help='.npy')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-r', '--recursive', action='store_true', default=False, help='do fitting 2 times')
    parser.add_argument('-o', '--out', required=True, help='output file path (.trr or .npy)')
    parser.add_argument('-w', '--max_wokers', default=2, type=int, help='max_wokers of multi-process')
    args = parser.parse_args()


    ### read file ###
    trj_mdtraj = md.load_trr(args.trajectory, top=args.topology)
    if args.npy:
        trj_mdtraj.xyz = np.load(args.npy)
    n_frames = trj_mdtraj.n_frames
    print(f'Trajectory Info ({n_frames} frames, {trj_mdtraj.n_atoms} atoms)')


    ### preprocess ###
    trj_mdtraj, atomlist, wlist = preprocess(trj_mdtraj)


    ### fitting ###
    trj_array = recursive_fitting(trj_mdtraj.xyz, wlist, args.max_wokers, args.recursive)


    ### ndarray to trr ###
    topology = trj_mdtraj.topology
    trj_mdtraj = md.Trajectory(trj_array, topology)


    ### save ###
    ext = os.path.splitext(args.out)[1]
    if ext == ".trr":
        trj_mdtraj.save_trr(args.out)
    elif ext == ".npy":
        np.save(args.out, trj_mdtraj.xyz)


if __name__=='__main__':
    main()