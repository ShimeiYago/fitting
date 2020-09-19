#!/usr/bin/env python

import mdtraj as md
import numpy as np
import pandas as pd
import os
import argparse

from utils import preprocess
from utils import recursive_fitting

MAINCHAIN_TRR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mainchain.trr")

def main():
    parser = argparse.ArgumentParser(description='fitting trajectory.')
    parser.add_argument('-t', '--trr', default=MAINCHAIN_TRR, help='trajectory file (.trr)')
    parser.add_argument('--trj', required=True, help='.npy')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-r', '--recursive', action='store_true', default=False, help='do fitting 2 times')
    parser.add_argument('-o', '--out', required=True, help='output file path (filename.trr or filename.npy or filename)')
    parser.add_argument('-w', '--max_wokers', default=2, type=int, help='max_wokers of multi-process')
    args = parser.parse_args()


    ### read file ###
    trj_mdtraj = md.load_trr(args.trr, top=args.topology)

    if args.trj:
        trj_mdtraj.xyz = load_trj(args.trj)

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
    else:
        trj_mdtraj.save_trr(args.out + '.trr')
        np.save(args.out + '.npy', trj_mdtraj.xyz)


def load_trj(fp):
    ext = os.path.splitext(fp)[1]
    if ext == '.npy':
        return np.load(fp)

    elif ext == '.xvg':
        data = np.loadtxt(fp, comments='@', delimiter='\t', skiprows=14)[:, 1:]
        data = data.reshape(data.shape[0], -1, 3)
        return data


if __name__=='__main__':
    main()