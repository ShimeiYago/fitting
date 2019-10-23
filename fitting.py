#!/usr/bin/env python3

import numpy as np
import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description='fitting trajectory by super-impose')
    parser.add_argument('-n', '--npz', required=True, help='trajectory and weight-list (.npz)')
    parser.add_argument('-o', '--outprefix', default="./fitted", help='output file prefix')
    parser.add_argument('-i', action='store_true', help='fit to init structure. if not called, fit to mean structure.')
    args = parser.parse_args()


    ### read file ###
    npz = np.load(args.npz)
    trj = npz['trj']    
    n_frames = trj.shape[0]
    print(f'fitting the trajectory ({n_frames} frames)')

    wlist = npz['wlist']


    ### decide reference_structure ###
    if args.i:
        reference_structure = trj[0]
    
    else:
        reference_structure = trj.mean(axis=0) # mean structure


    ### fitting ###
    fitted_trj = np.empty_like(trj)
    for i in range(n_frames):
        fitted_structure = super_impose(trj[i], reference_structure, wlist)
        fitted_trj[i] = fitted_structure

        if i % 10 == 0 or i+1 == n_frames:
            progress_frames = i+1
            progress_percentage = int(progress_frames/n_frames * 100)
            print(f"\rprogress: {progress_percentage}% ({progress_frames} frames)", end="")
    
    print(f"\ncompleted!")


    ### save ###
    outpath = f"{args.outprefix}.npz"
    np.savez(outpath, trj=fitted_trj, wlist=wlist)



def make_weightlist(trj, weight_key='standard'):
    atoms_dict = trj.topology.to_dataframe()[0].element

    wlist = [float(wdict[atom][weight_key]) for atom in atoms_dict]
    return wlist


def super_impose(target_structure:np.ndarray, reference_structure:np.ndarray, wlist):
    ### matrix U ###
    U = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            U[i][j] = sum( [wlist[n]*target_structure[n,i]*reference_structure[n,j] for n in range(target_structure.shape[0])] )

    ### matrix OMEGA ###
    OMEGA = np.empty((6,6))
    for i in range(3):
        for j in range(3):
            OMEGA[i][j] = 0
            OMEGA[i+3][j+3] = 0
            OMEGA[i+3][j] = U[i][j]
            OMEGA[i][j+3] = U.T[i][j]


    ### resolve Eigenvalue problem ###
    eig_val, eig_vec =np.linalg.eig(OMEGA)

    omegas = eig_vec.T

    ### split eig_vec to h and k ###
    hks = np.array([omegas[i] for i in range(6) if eig_val[i] > 0])
    hks = hks * np.sqrt(2) # root2

    k = np.empty((3,3))
    h = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            k[i][j] = hks[j][i]
            h[i][j] = hks[j][i+3]
    

    ### rotation matrix R ###
    R = np.empty((3,3))
    for i in range(3):
        for j in range(3):
            R[i][j] = sum( [k[i][a]*h[j][a] for a in range(len(h))] )

    
    ### target_trj to fitted_trj ###
    fitted_structure = np.empty_like(target_structure)
    for i in range(3):
        for n in range(fitted_structure.shape[0]):
            fitted_structure[n,i] = sum( [R[i,j]*target_structure[n,j] for j in range(3)] )
    
    
    return fitted_structure


if __name__ == '__main__':
    main()