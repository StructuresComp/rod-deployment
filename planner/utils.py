import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.spatial.transform import Rotation as R
import torch
import copy
from scipy import interpolate

global ON_GPU
if torch.cuda.is_available():
    float_tensor = torch.cuda.FloatTensor
    get_numpy = lambda x: x.detach().cpu().numpy()
    ON_GPU = True
else:
    float_tensor = torch.FloatTensor
    get_numpy = lambda x: x.data.numpy()
    ON_GPU = False


def computeSegLen(pattern):
    tangent = np.diff(pattern, axis = 0)
    segL = np.sum(tangent**2, 1)**0.5
    # define voronoiLen
    segL = segL.reshape((segL.shape[0], 1))
    return segL

def computePatternTangent(pattern):
    tangent = np.diff(pattern, axis = 0)
    segL = np.sum(tangent**2, 1)**0.5
    # compute tangent
    tangent = tangent/segL[:, None]
    # allocate tangent for voronoi distribution
    Tangent = tangent.reshape((tangent.shape[0], 3))
    return Tangent

def computeCurvature(pattern, Tangent, segL):
    Kappa = np.zeros((Tangent.shape[0], 1), dtype='float32')
    for i in range(Tangent.shape[0]-1):
        t1 = Tangent[i, :]
        t2 = Tangent[i+1,:]
        kappa = 2 * np.cross(t1, t2)/(np.dot(t1, t2) + 1)
        if kappa[2] < 0:
            Kappa[i+1] = - np.linalg.norm(kappa)/segL[i]
        else:
            Kappa[i+1] = np.linalg.norm(kappa)/segL[i]

    Kappa[0] = Kappa[1]
    return Kappa

def computeCurvature(pattern, segL):
    Kappa = np.zeros((pattern.shape[0], 1), dtype='float32')
    Tangent = computePatternTangent(pattern)
    for i in range(Tangent.shape[0]-1):
        t1 = Tangent[i, :]
        t2 = Tangent[i+1,:]
        kappa = 2 * np.cross(t1, t2)/(np.dot(t1, t2) + 1)
        if kappa[2] < 0:
            Kappa[i+1] = - np.linalg.norm(kappa)/segL[i]
        else:
            Kappa[i+1] = np.linalg.norm(kappa)/segL[i]

    Kappa[0] = Kappa[1]
    Kappa[-1] = Kappa[-2]
    return Kappa


def computeTangent(pattern, Kappa, segL):
    tangent = computePatternTangent(pattern)
    Tangent = np.zeros((tangent.shape[0]+1, 3), dtype="float32")
    Tangent[0, :] = tangent[0, :]
    Tangent[1:-1,:] = (tangent[:-1,:] + tangent[1:, :])/2
    # compute the last tangent
    kapn1 = Kappa[-1] * segL[-1]
    phi = 2* np.arctan(abs(kapn1)/2.0)[0]
    if kapn1 < 0:
        phi = -phi
    rot = np.matrix([[np.cos(phi), -np.sin(phi), 0],
                     [np.sin(phi), np.cos(phi), 0],
                     [0, 0, 1]], dtype = 'float32')
    t1 = tangent[-1,:] @ rot.transpose()
    Tangent[-1, :] = (t1 + tangent[-1, :])/2.0
    Tangent = Tangent/(np.sum(Tangent**2, 1)**0.5)[:, None]

    return Tangent


def convert2Matrix(axisA, m0, Kappa):
    Rot = np.zeros((axisA.shape[0], 3, 3), dtype='float32')
    for i in range(axisA.shape[0]):
        e = axisA[i, :]
        theta = np.linalg.norm(e)
        e = e/theta
        e_x = np.matrix([[0, -e[2], e[1]],[e[2], 0, -e[0]],[-e[1], e[0], 0]], dtype='float32')
        rot = np.identity(3) + np.sin(theta) * e_x + (1 - np.cos(theta)) * (e_x @ e_x)
        rot = rot @ m0
        # update rot
        if Kappa[i] > 0:
            rot[0, 0] = -rot[0, 0]
            rot[1, 1] = -rot[1, 1]
            rot[2, 1] = -rot[2, 1]

        rot[:,0] = -rot[:,0]
        rot[:,1] = -rot[:,1]
        Rot[i, : ] = rot.reshape(1, 3, 3)
    return Rot


def visualizePatternAndTraj(pattern, Traj):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot(pattern[:,0], pattern[:,1], pattern[:,2], 'red')
    ax.plot(Traj[:,0], Traj[:,1], Traj[:,2], 'blue')

    # add quiver for the end-effector frame
    _m1 = np.zeros((Traj.shape[0], 3))
    _m2 = np.zeros(_m1.shape)
    _t = np.zeros(_m1.shape)

    for i in range(Traj.shape[0]):
        quat = Traj[i, 3:]
        rot = R.from_quat(quat)
        rot = rot.as_matrix()
        _m1[i, :] = rot[:, 0].reshape(1, 3)
        _m2[i, :] = rot[:, 1].reshape(1, 3)
        _t[i, :] = rot[:, 2].reshape(1, 3)

    ind = range(0, Traj.shape[0], 2)

    ax.quiver(Traj[ind,0], Traj[ind, 1], Traj[ind,2], _m1[ind, 0], _m1[ind,1], _m1[ind,2], length = 0.05, color='blue')
    ax.quiver(Traj[ind,0], Traj[ind, 1], Traj[ind,2], _m2[ind, 0], _m2[ind,1], _m2[ind,2], length = 0.05, color='red')
    ax.quiver(Traj[ind,0], Traj[ind, 1], Traj[ind,2], _t[ind, 0], _t[ind,1], _t[ind,2], length = 0.05, color="green")
    # ax.set_aspect("equal")



    X = np.array([-(Traj[:,1].max() + 0.1)/2.0, (Traj[:,1].max() + 0.1)/2.0])
    Y = np.array([-0.1, Traj[:,1].max()+0.1])
    Z = np.array([-0.1, Traj[:,2].max()+0.1])
    # ax.plot(X, Y, Z, 'w')
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.set_xlabel("x")
    ax.set_ylabel("y")



    plt.show()

def preProcessPattern(pattern):
    tangent = np.diff(pattern, axis = 0)
    segL = np.sum(tangent**2, 1)**0.5
    segL = np.cumsum(segL)
    segL = np.insert(segL, 0, 0)
    factor0 = segL[-1]
    segL /= factor0
    pattern /= factor0
    x = pattern[:, 0]
    y = pattern[:, 1]
    tcx = interpolate.splrep(segL, x, s = 0.001)
    tcy = interpolate.splrep(segL, y, s = 0.001)


    sl = np.linspace(0, 1, pattern.shape[0])
    sx = np.expand_dims(interpolate.splev(sl, tcx), 1)
    sy = np.expand_dims(interpolate.splev(sl, tcy), 1)


    node = np.concatenate((sx, sy, np.zeros((sx.shape[0], 1))), axis = 1)

    tangent = np.diff(node, axis = 0)
    segL = np.sum(tangent**2, 1)**0.5
    segL = np.cumsum(segL)
    factor = segL[-1]
    node = node/factor *  factor0

    return node
