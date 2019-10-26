#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import os
import sys
import argparse
from lib.weightdict import atomicWeightsDecimal as wdict


def main():
    parser = argparse.ArgumentParser(description='fit center of gravity to origin coordinates')
    parser.add_argument('-t', '--trajectory', required=True, help='trajectory file (.trr)')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-o', '--out', default="./", help='output dir or file path(.npz)')
    args = parser.parse_args()

    if os.path.isdir(args.out):
        outpath = os.path.join(args.out, "preprocessed.npz")
    
    else:
        outpath = args.out
        

    ### read file ###
    trj = md.load_trr(args.trajectory, top=args.topology)
    n_frames = trj.n_frames

    print(f'Trajectory Info ({n_frames} frames, {trj.n_atoms} atoms)')


    ### make weight list ###
    wlist = make_weightlist(trj)


    ### centering ###
    trj.center_coordinates(mass_weighted=True) # centering


    ### save ###
    np.savez(outpath, trj=trj.xyz, wlist=wlist)


def make_weightlist(trj, weight_key='standard'):
    atoms_dict = trj.topology.to_dataframe()[0].element

    wlist = [float(wdict[atom][weight_key]) for atom in atoms_dict]
    return wlist


if __name__ == '__main__':
    main()