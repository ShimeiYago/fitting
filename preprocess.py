#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import os
import argparse
from dicts.weightdict import atomicWeightsDecimal as wdict


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
    atomlist, wlist = make_weightlist(trj.topology)


    ### centering ###
    trj.center_coordinates(mass_weighted=True) # centering


    ### save ###
    np.savez(outpath, trj=trj.xyz, atomlist=atomlist, wlist=wlist)


def make_weightlist(top, weight_key='standard'):
    atomlist = [atom for atom in top.to_dataframe()[0].name]
    wlist = [float(wdict[atom][weight_key]) for atom in top.to_dataframe()[0].element]
    return atomlist, wlist


if __name__ == '__main__':
    main()