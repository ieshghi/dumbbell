from numpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def get_traj(path):
    q = loadtxt(path)
    x1 = q[:,0]
    x2 = q[:,1]
    y1 = vpot(x1)
    y2 = vpot(x2)

    return x1[::10],y1[::10],x2[::10],y2[::10]

def parabpot(x):
    u=1
    z=((x%10+10)%10)
    if z<9:
        out=u/81*(z-9)**2
    else:
        out=u*(z-9)**2
    return out
vpot = vectorize(parabpot)

def drawpot(xmin,xmax):
    xvals = linspace(xmin,xmax,1000)
    yvals = vpot(xvals)

    return xvals,yvals

