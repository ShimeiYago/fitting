#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import os
import sys
import argparse
from weightdict import atomicWeightsDecimal as wdict


def main():
    parser = argparse.ArgumentParser(description='fitting trajectory by super-impose')
    parser.add_argument('-t', '--trajectory', required=True, help='trajectory file (.trr)')
    parser.add_argument('-p', '--topology', required=True, help='topology file (.gro, .pdb)')
    parser.add_argument('-o', '--output', default="fitted.trr", help='output file path (default: output/fitted.trr)')
    args = parser.parse_args()


    ### read file ###
    trj = md.load_trr(args.trajectory, top=args.topology)
    n_frames = trj.n_frames

    ### make weight list ###
    wlist = make_weightlist(trj)

    ### centering ###
    trj.center_coordinates(mass_weighted=True)


    ### fitting ###
    print(f'fitting the trajectory ({n_frames} frames, {trj.n_atoms} atoms)')

    trj_array = trj.xyz
    fitted_trj_array = np.empty_like(trj_array)

    for i in range(n_frames):
        fitted_structure = super_impose(trj_array[i], trj_array[0], wlist)
        fitted_trj_array[i] = fitted_structure

        if i % 10 == 0 or i+1 == n_frames:
            progress_frames = i+1
            progress_percentage = int(progress_frames/n_frames * 100)
            print(f"\rprogress: {progress_percentage}% ({progress_frames} frames)", end="")
    
    print(f"\ncompleted!")

    ### save ###
    fitted_trj = md.Trajectory(fitted_trj_array, trj.topology)
    fitted_trj.save_trr(args.output)



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