from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_matrix
import time
import scipy.io as sio
from netCDF4 import Dataset


def solve(method,A,b,lev):
        N = b.size
   #     h = 1
    
        if method == "inbuilt":
        #
                u= spsolve(A,b[lev,:])
        
        u = np.reshape(u, [150,165])
    #

 #  u = u_tmp.copy()
        
 #   midpt = (n+1)/2
 #   lapl_midpt = (u[midpt+1,midpt] + u[midpt-1,midpt] - 4*u[midpt,midpt] + u[midpt,midpt + 1] + u[midpt,midpt -1])/h**2
 #   err = 2. - lapl_midpt
 #   print("Laplacian at mid point = {}, Error = {}".format(lapl_midpt, err)) 
        
        return u



if __name__ == '__main__':
        start=time.clock()
        n = 150
        N = 150*165
        h = 36000
        A=sio.mmread('A1.mtx')
        A = A/h/h
        A.tocsr()
        
        f = Dataset('mawa_u2012_24h_fff.nc')
        ff = f.variables['ff']
        fa= np.array(ff)
        b = np.reshape(fa,[34,150*165])

        adssv1 = np.zeros((34,150,165))
        for lev in range(0,34):
                adssv1[lev,:,:] = solve("inbuilt",A,-b,lev)
        
        da=Dataset("mawa_adssv1","w",format="NETCDF4")
        da.createDimension("i",150)
        da.createDimension("j",165)
        da.createDimension("k",34)
        da.createVariable("adssv1","f8",("k","i","j")) 
        da.variables["adssv1"][:]=adssv1
        da.close()













