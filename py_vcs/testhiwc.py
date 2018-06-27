##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

get_ipython().run_line_magic('pylab', 'inline')
# For interactive plotting
# %matplotlib notebook
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import numdifftools as nd
import pandas as pd
import xarray as xr
import glob, re, os
import scipy
import pickle
import datetime
import netCDF4

## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger


#%%

## Load data and functions
with open('tmpcay.p','rb') as f:
    psd,nml,_dump_t,bin_div,_dump_bin_mid,_dump_bin_diff = pickle.load(f)
    bin_mid = (bin_div[1:]+bin_div[:-1])/2
    bin_diff = np.diff(bin_div)

moments = np.array([0,2,3])
mobs = np.empty_like(moments)
for szmoment in range(len(mobs)):
    mobs[szmoment] = np.sum( psd*bin_diff*bin_mid**(moments[szmoment]) )
    
### IGF routine part2 ellipsoid generator
sa2ds = minimum(5.13*bin_mid**2,63e3) * 1280; # 2DS in um
sapip = minimum(3*bin_mid**2/0.658, 260e3) * 6400; # PIP in um
sa = np.where(bin_mid<1e3, sa2ds, sapip)/1e12 # m^2


#%%

"""
Adjustable parameters:
    bin_div (bin_mid, bin_diff) -> 2ds+pip
    sa -> 2ds+pip
    moments -> np.array([0,2,3])
"""
# Define functions
f_mygamma = lambda nml, x: 10**nml[0]*x**nml[1]*np.exp(-nml[2]*x)
def f_geneigvector( psd, bin_div=bin_div, moments=moments ):
    ## Calculate the counts
    global bin_mid, bin_diff
    count = f_count(psd,bin_mid,bin_diff,200*5)
    dchi = f_delta_chisquare( psd, count, moments, bin_mid, bin_diff, mobs )

    jac = np.matrix(nd.Jacobian(f_sum_chisquare,step=np.array([1,1,1e-3])*1e-3)( nml, moments, bin_mid, bin_diff, mobs ))
    hes = np.matrix(nd.Hessian(f_sum_chisquare,step=np.array([1,1,1e-3])*1e-3)( nml, moments, bin_mid, bin_diff, mobs ))

    # HV = Vd = VD
    d, V = scipy.linalg.eigh(hes)
    # D = np.asmatrix(np.diag(d))
    V = np.asmatrix(V)
    centroid = nml
    return centroid, V, d, dchi # sync1

def f_plot_E(args):
    centroid, V, d, dchi = args # sync1
    rx, ry, rz = np.sqrt(2*dchi/d)
    u, v = np.mgrid[0:2*pi:20j, -pi/2:pi/2:10j]

    x = rx*cos(u)*cos(v)
    y = ry*sin(u)*cos(v)
    z = rz*sin(v)

    E = np.stack([x,y,z],axis=-1)
    # BUG found Jun 10, 2017, use V.T instead of V
    E = np.dot(E,np.array(V).T) + centroid
    # Move the xyz dimension from the last to the first so that it can be assigned.
    E = np.rollaxis(E, axis = -1)
    
    x, y, z = E
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X=x, Y=y, Z=z, color='b', alpha = 0.75)
    plt.show()

# Plot families of PSD based on nml
from scipy.stats import chi2
def f_plot_Ef(nmls, prob=0.95,angles=None):
    """Plot families of PSD based on nml.
    nmls: needs to be shape [3,time]
    prob: confidential interval
    angles: final plot view angles (elevation, azimuth)"""
    
    # Getting rid of the outliers based on prob first
    centroid = np.mean(nmls,axis=1)
    dist = np.sqrt(np.sum((nmls.T-centroid).T**2,axis=0))
    ind = dist < np.percentile(dist, 100*prob)
    nmls = nmls[:,ind]
    
    hes = np.cov(nmls)
    # HV = Vd = VD
    d, V = scipy.linalg.eigh(hes)
    # D = np.asmatrix(np.diag(d))
    V = np.asmatrix(V)
    centroid = np.mean(nmls,axis=1)
    # return centroid, V, d, dchi # sync1

    dchi = chi2.isf(1-prob,3) # inverse survival function, 99 percent, degree of freedom 3
    rx, ry, rz = np.sqrt(dchi*d) # This is different from the Hessian matrix, see
    # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/#comment-190
    u, v = np.mgrid[0:2*pi:20j, -pi/2:pi/2:10j]

    x = rx*cos(u)*cos(v)
    y = ry*sin(u)*cos(v)
    z = rz*sin(v)

    E = np.stack([x,y,z],axis=-1)
    # BUG found Jun 10, 2017, use V.T instead of V
    E = np.dot(E,np.array(V).T) + centroid
    # Move the xyz dimension from the last to the first so that it can be assigned.
    E = np.rollaxis(E, axis = -1)

    x, y, z = E
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if angles is not None:
        ax.view_init(*angles)

    ax.plot_surface(X=x, Y=y, Z=z, color='g', alpha = 0.5)
    ax.set_xlabel('log10(N0) (1/L/um)')
    ax.set_ylabel('$\mu$')
    ax.set_zlabel('$\lambda$ (1/um)')

    x,y,z = nmls
    ax.scatter(x, y, z, color='b', alpha = 0.75)
    plt.show()
    return

def f_count( psd, bin_mid, bin_diff, tasdt ):
    global sa
    return psd*sa*tasdt*bin_diff*1e3 # convert psd in /L to /m3

def f_delta_chisquare( psd, count, moments, bin_mid, bin_diff, mobs ):
    count[count==0] = nan
    psd[psd==0] = nan
    return nansum( np.array([ nansum( (bin_mid**m *bin_diff*psd/mob)**2/count )
                          for m,mob in zip(moments,mobs)]))

def f_sum_chisquare( nml, moments, bin_mid, bin_diff, mobs ):
    psd = f_mygamma(nml, bin_mid)
    ### We may drop this condition if the code works fine.
#     if any(np.isnan(psd)) or any(np.isinf(psd)):
#         return np.inf
    mfit = np.empty(mobs.shape)
    for szmoment in range(len(mfit)):
        mfit[szmoment] = np.sum( psd*bin_diff*bin_mid**(moments[szmoment]) )
    return np.sum( (mfit-mobs)**2/mobs/mfit )

## For median mass diameter calculation
def f_mmd(bin_div, msd):
    # bin_diff = np.diff(bin_div)
    cmsd = np.concatenate( (np.array([0]),np.cumsum(msd*bin_diff)) )
    if cmsd[-1]<=0:
        return np.NaN
    cmsd /= cmsd[-1]
    indtmp = np.where(np.diff(cmsd>0.5)==1)[0]
    x1,x2,y1,y2 = bin_div[indtmp], bin_div[indtmp+1], cmsd[indtmp], cmsd[indtmp+1]
    mmd = (x2-x1)/(y2-y1)*(0.5-y1)+x1
    return mmd

%matplotlib notebook
# %matplotlib inline
# f_plot_E(f_geneigvector(psd))
# nmltest = np.random.randn(3,100)
nmltest = (np.dot(np.random.randn(3,3),np.random.randn(3,500)).T*np.array([2,4,70])+2).T
f_plot_Ef(nmltest)
plt.show()## Needed only for initializing data
test = xr.open_dataset('tmp/cayenne/cayenne_sync_bulk.nc',group='/saffire')
ikp = xr.open_dataset('tmp/cayenne/cayenne_sync_bulk.nc',group='/ikp')['XKBZR5s']
tmppsd = xr.open_dataset('tmp/cayenne/cayenne_sync_gz.nc',group='/lamp')
lampproc = xr.open_dataset('tmp/cayenne/cayenne_sync_bulk.nc',group='/lampproc')
bin_div = tmppsd.attrs['bin_div']
bin_diff = np.diff(bin_div)
bin_mid = (bin_div[:-1]+bin_div[1:])/2.

with open('tmpcay.p','wb') as f:
    pickle.dump([psd,nml,t,bin_div,bin_mid,bin_diff],f)
# psd = tmppsd.psddmax[ikp>2,:][2,:]# Code for KDE, reference:
# http://stackoverflow.com/a/6658307/5426033
import numpy as np
import scipy.ndimage as ndi

data = np.random.rand(30000,2)           ## create random dataset
bin_div = np.linspace(0,1,257)
gridpts,null,null = np.histogram2d(data[:,0],data[:,1],bins=[bin_div,bin_div])

img2 = ndi.gaussian_filter(gridpts, (10,10))
plt.pcolor(img2)
plt.colorbar()import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pi = np.pi
sin = np.sin
cos = np.cos
# minimum volume enclosing ellipsoid
def mvee(points, tol = 0.001):
    """
    Finds the ellipse equation in "center form"
    (x-c).T * A * (x-c) = 1
    """
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T
    err = tol+1.0
    u = np.ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = np.dot(np.dot(Q, np.diag(u)), Q.T)
        M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
    c = np.dot(u,points)        
    A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points)
               - np.multiply.outer(c,c))/d
    return A, c