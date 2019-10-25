#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='convert npz to trr')
    parser.add_argument('-n', '--npz', required=True, help='trajectory and weight-list (.npz)')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-o', '--outprefix', default="./fitted", help='output file prefix')
    args = parser.parse_args()


    ### read file ###
    npz = np.load(args.npz)
    trj = npz['trj']    

    topology = md.load(args.topology).topology

    ### save ###
    outpath = f"{args.outprefix}.trr"
    trj_mdtraj = md.Trajectory(trj, topology)
    trj_mdtraj.save_trr(outpath)


if __name__ == '__main__':
    main()