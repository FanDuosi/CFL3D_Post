
from fortio import FortranFile
import numpy as np


def read_mesh(filename):
    f = FortranFile(filename)
    
    # read nblock
    nblock = f.read_record('i4')
    # read ni,nj,nk
    ndim = f.read_record('i4')
    
    # read coordinate
    tmp = f.read_record('f4',shape=(4,ndim[0],ndim[1],ndim[2]))
    
    f.close()
    
    return tmp[0,:,:,:],tmp[1,:,:,:],tmp[2,:,:,:]


def read_solution(filename):
    f = FortranFile(filename)
    
    # read nblock
    nblock = f.read_record('i4')
    
    # read ni,nj,nk
    ndim = f.read_record('i4')
    
    # read ma, alpha, Re,time
    ma,alpha,re,time = f.read_record('i4')
    
    # read q data (rho, rhou, rhov, rhow, rhoe)
    q = f.read_record('f4',shape=(5,ndim[0],ndim[1],ndim[2]))
    
    f.close()
    
    return q

def con2prim(q):
    gamma = 1.4
    
    w = np.zeros_like(q)
    
    # density
    w[0,:,:,:] = q[0,:,:,:]
    
    # velocity
    w[1,:,:,:] = q[1,:,:,:]/q[0,:,:,:]
    w[2,:,:,:] = q[2,:,:,:]/q[0,:,:,:]
    w[3,:,:,:] = q[3,:,:,:]/q[0,:,:,:]
    
    # pressure
    w[4,:,:,:] = (gamma - 1.) * (q[4,:,:,:] - 
                                 0.5 * q[0,:,:,:] *
                                 (w[1,:,:,:]**2 + w[2,:,:,:]**2 + w[3,:,:,:]**2))
    
    return w