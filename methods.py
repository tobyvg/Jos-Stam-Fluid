import numpy as np
import matplotlib.pyplot as plt
from project_utilities import *
import time
init_mpl(150,mat_settings = True)
from IPython.display import clear_output
import pygame
from numba import prange

@numba.njit()
def set_bnd(N,b,x):
    if b == 0:
        for i in prange(N+2):
            x[0,i] = x[1,i]
            x[i,0] = x[i,1]
            x[N+1,i] = x[N,i]
            x[i,N+1] = x[i,N]
    elif b == 1:
        for i in prange(N+2):
            x[0,i] = -x[1,i]
            x[N+1,i] = -x[N,i]
            x[i,0] = x[i,1]
            x[i,N+1] = x[i,N]
    elif b == 2:
        for i in prange(N+2):
            x[0,i] = x[1,i]
            x[N+1,i] = x[N,i]
            x[i,0] = -x[i,1]
            x[i,N+1] = -x[i,N]


    x[0,0] = 1/2*(x[1,0]+x[0,1])
    x[0,N+1] = 1/2*(x[1,N+1]+x[0,N])
    x[N+1,0] = 1/2*(x[N,0]+x[N+1,1])
    x[N+1,N+1] = 1/2*(x[N,N+1]+x[N+1,N])

@numba.njit()
def add_source(x,s,dt):
    x += dt*s

@numba.njit()
def diffuse(N,b,x,x0,diff,dt):
    a = dt*diff*N**2
    for k in range(20):
        for i in range(1,N+1):
            for j in range(1,N+1):
                x[i,j] = (x0[i,j] + a*(x[i-1,j]+x[i+1,j]+x[i,j-1]+x[i,j+1]))/(1+4*a)
    set_bnd(N,b,x)

@numba.njit()
def advect(N,b,d,d0,u,v,dt):
    dt0 = N*dt
    for i in prange(1,N+1):
        for j in prange(1,N+1):
            x = i-dt0*u[i,j]
            y = j-dt0*v[i,j]
            if x < 0.5:
                x = 0.5
            elif x > N + 0.5:
                x = N + 0.5
            if y < 0.5:
                y = 0.5
            elif y > N + 0.5:
                y = N + 0.5
            i0 = int(np.floor(x))
            i1 = i0+1
            j0 = int(np.floor(y))
            j1 = j0+1
            s1 = x -i0
            s0 = 1 - s1
            t1 = y-j0
            t0= 1-t1
            d[i,j] = s0*(t0*d0[i0,j0]+t1*d0[i0,j1]) + s1*(t0*d0[i1,j0] + t1*d0[i1,j1])
    set_bnd(N,b,d)


@numba.njit()
def dens_step(N,x,x0,u,v,diff,dt,s):
    add_source(x,s,dt)
    diffuse(N,0,x0,x,diff,dt)
    advect(N,0,x,x0,u,v,dt)

@numba.njit()
def project(N,u,v,p,div):
    h = 1/N
    for i in prange(1,N+1):
        for j in prange(1,N+1):
            div[i,j] = -0.5*h*(u[i+1,j]-u[i-1,j] + v[i,j+1]-v[i,j-1])
            p[i,j] = 0
    set_bnd(N,0,div)
    set_bnd(N,0,p)
    for k in range(20):
        for i in range(1,N+1):
            for j in range(1,N+1):
                p[i,j] = (div[i,j]+p[i-1,j]+p[i+1,j]+p[i,j-1]+p[i,j+1])/4
    set_bnd(N,0,p)

    for i in prange(1,N+1):
        for j in prange(1,N+1):
            u[i,j] -= 0.5*(p[i+1,j]-p[i-1,j])/h
            v[i,j] -= 0.5*(p[i,j+1]-p[i,j-1])/h
    set_bnd(N,1,u)
    set_bnd(N,2,v)

@numba.njit()
def vel_step(N,u,v,u0,v0,visc,dt,su,sv):
    add_source(u,su,dt)
    add_source(v,sv,dt)
    diffuse(N,1,u0,u,visc,dt)
    diffuse(N,1,v0,v,visc,dt)
    project(N,u0,v0,u,v)
    advect(N,1,u,u0,u0,v0,dt)
    advect(N,2,v,v0,u0,v0,dt)
    project(N,u,v,u0,v0)

import matplotlib as mpl
from matplotlib import cm

class MplColorHelper:
  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)
