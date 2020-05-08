import mdtraj as md
import numpy as np
import os
import argparse
from .weightdict import atomicWeightsDecimal as wdict


def preprocess(trj:md.core.trajectory.Trajectory):

    ### make weight list ###
    atomlist, wlist = make_weightlist(trj.topology)


    ### centering ###
    trj.center_coordinates(mass_weighted=True) # centering


    ### save ###
    return trj, atomlist, wlist


def make_weightlist(top, weight_key='standard'):
    atomlist = [atom for atom in top.to_dataframe()[0].name]
    wlist = [float(wdict[atom][weight_key]) for atom in top.to_dataframe()[0].element]
    return atomlist, wlist

