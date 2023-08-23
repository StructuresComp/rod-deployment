#!/usr/bin/env python3
import numpy as np
import math
from planner.utils import visualizePatternAndTraj, computeSegLen, computeSegLen,\
                          computeTangent, computeCurvature, convert2Matrix, computePatternTangent
import time
from scipy.spatial.transform import Rotation as R
import os
import subprocess
from multiprocessing import Pool
import matplotlib.pyplot as plt
import sys

import torch
from planner.deployment_planner import DeploymentPlanner
import time


#############################################################
# RunBASim function will simply execute a single command
def RunBASim(cmdline):
    BASimProcess = subprocess.Popen(cmdline, shell=True)
    BASimProcess.communicate()
#############################################################


def computeOptPathwithSim(ls, Lgb, ks, pose, pattern):
    rot0 = R.from_quat(pose[3:])
    rot0 = rot0.as_matrix()
    m0 = np.matrix([[0, 0, 1], [1, 0, 0], [0, 1 , 0]], dtype = 'float32')
    segL = computeSegLen(pattern)
    Kappa = computeCurvature(pattern, segL)

    Tangent = computeTangent(pattern, Kappa, segL)
    # compute curvature
    z_axis = np.array([0, 0, 1.0], dtype="float32")
    # get the inputs
    ls_ = (ls - np.cumsum(segL)).reshape(segL.shape)/Lgb
    _kappa = np.absolute(Kappa * Lgb)
    _ks = np.full(_kappa.shape, ks)
    _ls = np.zeros((pattern.shape[0], 1), dtype='float32')
    _ls[0] = ls/Lgb
    _ls[1:] = ls_
    _ls = _ls.squeeze()
    _kappa = _kappa.squeeze()
    _ks = _ks.squeeze()

    Pose = np.zeros((pattern.shape[0], 7), dtype="float32")
    # intialized guess
    kap0 = 0
    if _ls[0] > 10:
        ratio = 0.1
    elif _ls[0] < 1:
        ratio = 0.9
    else:
        ratio = -0.06*(_ls[0]-1) + 0.9

    if _ls[0] < 12:
        xTop = _ls[0] * ratio
        zTop = _ls[0] * (1-ratio)
    else:
        xTop = 3
        zTop = _ls[0] - xTop

    yTop = 0
    for i in range(pattern.shape[0]):
        ls = _ls[i]
        ks = _ks[i]
        kap = _kappa[i]
        # run sim
        cmd = "./simOpt option.txt -- hangLength {:} -- normKs {:}\
        -- kappaI {:} -- kappaT {:} -- xTop {:} -- yTop {:} -- zTop {:}".format(ls, ks, kap0, kap, xTop, yTop, zTop)
        # print(cmd)
        result = subprocess.run(cmd, capture_output=True, text=True, shell = True)

        # read the contents
        res = np.loadtxt("datafiles/simDER.txt")
        # get the last
        res = res.reshape((-1, 16))
        res = res[-1, :]
        Pose[i, :3] = res[:3] * Lgb
        # define rotation
        tmpP = pattern[i, :]
        sign = 1
        if Kappa[i] < 0:
            sign = -1
        tmpT = Tangent[i, :]
        tmpY_axis = np.cross(z_axis, tmpT)
        tmpRot = np.concatenate((tmpT.reshape(3, 1), tmpY_axis.reshape(3,1), z_axis.reshape(3,1)), axis = 1)
        rot = res[-9:]
        rot = rot.reshape((3, 3))
        rot = rot @ m0
        if Kappa[i] > 0:
            rot[0, 0] = -rot[0, 0]
            # rot[2, 0] = -rot[2, 0]
            rot[1, 1] = -rot[1, 1]
            rot[2, 1] = -rot[2, 1]

            # rot[:, 0] = np.cross(rot[:, 1], rot[:, 2])
        rot[:, 0] = -rot[:, 0]
        rot[:, 1] = -rot[:, 1]

        tmpRot = tmpRot @ rot
        if i == 0:
            tmpRot0 = tmpRot
            rot0 = np.matrix([[1, 0, 0],
                              [0, -1, 0],
                              [0, 0,  -1]], dtype = 'float32')
        tmpdcm = tmpRot @ tmpRot0.transpose() @ rot0
        quat = R.from_matrix(tmpdcm)
        Pose[i, 3:] = quat.as_quat()
        Pose[i, :3] = tmpP + Pose[i, 0] * tmpT - sign * Pose[i, 1] * tmpY_axis + Pose[i, 2] * z_axis

        kap0 = kap
        xTop, yTop, zTop = res[:3]
        if i%10 == 0:
            print(f"{i+1}/{pattern.shape[0]} is done")

    return Pose

if __name__ == "__main__":
    script_name = sys.argv[0]
    isInitial = sys.argv[1]
    solver = sys.argv[2]
    fileName = sys.argv[3]

    # define parameters (one example)
    h = 1.6e-3 # rod radius
    rho = 1180 # volume density
    Lgb = 1.8e-2 # gravito-bending length
    ls = 0.875 # initial suspended length
    g = 10

    E = Lgb**3 * 8 * rho * g/h**2

    ks = E*math.pi*h**2
    kb = E*math.pi*h**4/4.0

    normks = Lgb**2 * ks/kb
    # load pattern
    pattern = np.loadtxt(f"patterns/{fileName}")

    pose = [0, 0, ls, 1, 0, 0, 0]
    print(isInitial)
    if isInitial.lower() == "true":
        dp = DeploymentPlanner(Lgb, normks)
        Pose = dp.intuitiveTraj(ls, pose, pattern)
    else:
        if solver == "NN":
            ## NN solver
            dp = DeploymentPlanner(Lgb, normks)
            print(f"Begin planning for {fileName}")
            start_time = time.time()
            Pose = dp.optimalPath(ls, pose, pattern)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed time: {elapsed_time} seconds")
        else:
            ## Numeric Solver
            print(f"Begin planning for {fileName}")
            start_time = time.time()
            Pose = computeOptPathwithSim(ls, Lgb, normks, pose, pattern)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed time for {fileName}: {elapsed_time} seconds")

    # visulize result
    visualizePatternAndTraj(pattern, Pose)

    np.savetxt(f"pattern.txt", pattern, "%.5f")
    np.savetxt(f"traj_Sim.txt", Pose, "%.5f")
