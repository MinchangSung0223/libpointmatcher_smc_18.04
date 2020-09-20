from ctypes import *
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np
import open3d as o3d
import copy
import string
from ctypes import Structure,c_int,c_float,POINTER,cast,pointer,byref,CDLL
from scipy.spatial.transform import Rotation as R
from math import cos, sin, pi
def draw_registration_result(source, target, transformation):
    source_temp = copy.deepcopy(source)
    target_temp = copy.deepcopy(target)
    source_temp.paint_uniform_color([1, 0.706, 0])
    target_temp.paint_uniform_color([0, 0.651, 0.929])
    source_temp.transform(transformation)
    o3d.visualization.draw_geometries([source_temp, target_temp])
def RPY2T(r,p,y,t):
    T = np.eye(4,dtype=float)
    Rz = np.array([[cos(y), -sin(y),0],[sin(y),cos(y),0],[0,0,1]],dtype=float)
    Rx = np.array([[1,0,0],[0,cos(p),-sin(p)],[0,sin(p),cos(p)]],dtype=float)
    Ry = np.array([[cos(r), 0,sin(r)],[0,1,0],[-sin(r),0,cos(r)]],dtype=float)
    print(Rz)
    print(Rx)
    print(Ry)
    R = Rz*Ry*Rx
    T[0:3,0:3] = R
    T[0:3,3] = np.transpose(t)
    return T
lib = cdll.LoadLibrary("examples/icp_simple.so")  
test_icp= lib.Test_icp
test_icp.restype = ndpointer(dtype=ctypes.c_float, shape=(16))#ctypes.POINTE

class Row(Structure):
    _fields_ = [('cols_count', c_int), 
                ('cols', POINTER(c_float))]
    def __init__(self,cols):
        self.cols_count = cols
        # Allocate an array of character pointers
        pc = (c_float * cols)()
        self.cols = cast(pc,POINTER(c_float))            

class Unit(Structure):
    _fields_ = [('rows_count', c_int),
                ('rows',POINTER(Row))]
    def __init__(self,rows,cols):
        self.rows_count = rows
        # Allocate an array of Row structures.
        # This does NOT call __init__.
        pr = (Row * rows)()
        # Call init manually with the column size.
        for r in pr:
            r.__init__(cols)
        self.rows = cast(pr,POINTER(Row))

def changeUnit(u, pcd1,pcd2):
    points1 = np.asarray(pcd1.points).T
    points2 = np.asarray(pcd2.points).T
    r1 = np.asarray(pcd1.points).T.shape[1]
    c1 = np.asarray(pcd1.points).T.shape[0]
    r2 = np.asarray(pcd2.points).T.shape[1]
    c2 = np.asarray(pcd2.points).T.shape[0]
    print(r1,c1,r2,c2)
    for i in range(0,r1+r2):
       for j in range(c1):
           if i<r1:
              u.rows[i].cols[j] = points1[j,i]
           else:
              u.rows[i].cols[j] = points2[j,i-r1]
    return u
def dataSend(u,r1,r2):
    T=test_icp(byref(u),r1,r2)
    T = np.reshape(T,(4,4))
    return T.T

threshold = 0.02
target = o3d.io.read_point_cloud('examples/nivea_final.ply')
r = np.array(R.from_quat([0, 0, np.sin(np.pi/2), np.cos(np.pi)]).as_matrix())

T = np.eye(4)
T[:3, :3] = r
T[0, 3] = 0.01
T[1, 3] = 0.01
print(T)
source = copy.deepcopy(target).transform(T)
draw_registration_result(source,target,np.eye(4))
reg_p2l = o3d.registration.registration_icp(
        source, target, threshold, np.eye(4),
        o3d.registration.TransformationEstimationPointToPoint())
print(reg_p2l)
print("Transformation is:")
print(reg_p2l.transformation)
draw_registration_result(source,target, reg_p2l.transformation)

unit = Unit(np.asarray(target.points).shape[0]+np.asarray(source.points).shape[0],3)
u = changeUnit(unit,source,target)
T = dataSend(u,np.asarray(source.points).shape[0],np.asarray(target.points).shape[0])
print(T)
draw_registration_result(source,target, T)

