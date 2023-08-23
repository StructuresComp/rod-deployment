import numpy as np
import pickle
import torch
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from copy import deepcopy
from matplotlib.cm import ScalarMappable
from scipy.spatial.transform import Rotation as R
# from shapely.geometry import Polygon
from planner.trainedSolver import PoseInterpolator
from planner.utils import float_tensor, get_numpy, ON_GPU, computeSegLen, computeTangent, computeCurvature, convert2Matrix, computePatternTangent

class DeploymentPlanner:
    def __init__(self, Lgb, ks):
        self._nn = PoseInterpolator()
        if ON_GPU:
            self._nn = self._nn.cuda()
            self._device = 'cuda'
        else:
            self._device = 'cpu'

        self.Lgb = Lgb
        self.ks = ks


    def getLocalPose(self, _ls, _ks, _kappa):
        _std_data = np.concatenate((_ls, _ks, _kappa), axis=1)
        _std_data = self._nn.preprocess_input(_std_data)[1]
        # data_tensor = float_tensor(_std_data)
        out = get_numpy(self._nn(_std_data))
        return out

    def intuitiveTraj(self, ls, pose, pattern):
        # current suspended length,  current pose, and current pattern
        # define the initial rotation
        rot0 = R.from_quat(pose[3:])
        rot0 = rot0.as_matrix()
        segL = computeSegLen(pattern)
        Tangent = computePatternTangent(pattern)
        # determine the manipulating direction
        z_axis = np.array([0, 0, 1.0], dtype="float32")
        Pose = np.zeros((pattern.shape[0], 7), dtype="float32")
        # define the first point
        Pose[0, :3] = pose[:3]
        Pose[0, 3] = ls
        tmpdcm = rot0
        quat = R.from_matrix(tmpdcm)
        quat = quat.as_quat()
        Pose[0, 3:] = quat

        for i in range(pattern.shape[0]-1):
            tmpP = pattern[i, :]
            tmpT = Tangent[i, :]
            tmpY_axis = np.cross(z_axis, tmpT)
            tmpRot = np.concatenate((tmpT.reshape(3, 1), tmpY_axis.reshape(3,1), z_axis.reshape(3,1)), axis = 1)
            ls -= segL[i]
            Pose[i+1, :3] = tmpP + tmpT * segL[i] + ls * z_axis
            if i == 0:
                tmpRot0 = tmpRot

            tmpdcm = tmpRot @ tmpRot0.transpose() @ rot0
            quat = R.from_matrix(tmpdcm)
            quat = quat.as_quat()
            Pose[i+1, 3:] = quat

        return Pose

    def optimalPath(self, ls, pose, pattern):
        rot0 = R.from_quat(pose[3:])
        rot0 = rot0.as_matrix()
        m0 = np.matrix([[0, 0, 1], [1, 0, 0], [0, 1 , 0]], dtype = 'float32')
        segL = computeSegLen(pattern)
        Kappa = computeCurvature(pattern, segL)

        Tangent = computeTangent(pattern, Kappa, segL)
        # compute curvature
        z_axis = np.array([0, 0, 1.0], dtype="float32")
        # get the inputs
        ls_ = (ls - np.cumsum(segL)).reshape(segL.shape)/self.Lgb
        _kappa = np.absolute(Kappa * self.Lgb)
        _ks = np.full(_kappa.shape, self.ks)
        _ls = np.zeros((pattern.shape[0], 1), dtype='float32')
        _ls[0] = ls/self.Lgb
        _ls[1:] = ls_
        _out = self.getLocalPose( _ls, _ks, _kappa)
        _pose = _out[:,:3] * self.Lgb
        _axisA = _out[:,3:]
        rotTensor = convert2Matrix(_axisA, m0, Kappa)
        # convert local 2 global
        Pose = np.zeros((pattern.shape[0], 7), dtype="float32")
        # define the first position
        for i in range(pattern.shape[0]):
            tmpP = pattern[i, :]
            sign = 1
            if Kappa[i] < 0:
                sign = -1
            tmpT = Tangent[i, :]
            tmpY_axis = np.cross(z_axis, tmpT)
            tmpRot = np.concatenate((tmpT.reshape(3, 1), tmpY_axis.reshape(3,1), z_axis.reshape(3,1)), axis = 1)
            rot = rotTensor[i, : ,: ].reshape(3,3)
            tmpRot = tmpRot @ rot
            if i == 0:
                tmpRot0 = tmpRot
                # the rotation matrix of initial pose
                rot0 = np.matrix([[1, 0, 0],
                                  [0, -1, 0],
                                  [0, 0, -1]], dtype='float32')
            tmpdcm = tmpRot @ tmpRot0.transpose() @ rot0
            quat = R.from_matrix(tmpdcm)
            Pose[i, 3:] = quat.as_quat()
            # define the coordinate
            Pose[i, :3] = tmpP + _pose[i, 0] * tmpT - sign * _pose[i, 1] * tmpY_axis + _pose[i, 2] * z_axis

        return Pose
