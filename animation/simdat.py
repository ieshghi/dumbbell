from numpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def get_traj(path):
    q = loadtxt(path)
    x1 = q[:,0]
    x2 = q[:,1]
    y1 = vpot(x1)
    y2 = vpot(x2)

    return x1[::100],y1[::100],x2[::100],y2[::100]

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

