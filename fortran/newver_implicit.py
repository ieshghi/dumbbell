from numpy import *
import mesh
import scipy.sparse as spr
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
def runsim(nx,ny,nt):
    p,g,pars = init_stuff(nx,ny)
    psi_i = pinit(p[:,0],p[:,1])
    psi = psi_i.copy()
    np = psi.size
    psi_hist = zeros((np,nt//1))
    j = 0
#    psi_hist[:,0] = psi
    error = 1
    ls1 = lhs(p,g,pars,ifcrank = 1)
    ls2 = lhs(p,g,pars,ifcrank = 2)
#    ls = lhs(p,g,pars,ifcrank = 0)
#    ls = preplhs(p,g,pars,ls)
    ls1 = preplhs(p,g,pars,ls1)
    i = 0
    while error > 1e-10:
#    for i in range(nt):
        psi_new = imp_tstep(p,g,pars,ls1,psi,ls2)
        error = max(abs(psi_new-psi))
        print(i,sum(abs(psi_new)),sum(abs(psi_new))-sum(abs(psi)),error)
        psi = psi_new.copy()
        i = i+1
    return psi_hist,psi_i,psi,p

def pinit(xg,yg):
    p = exp(-(4*(xg)**2+(yg-5)**2))
    return p/sum(p)

def init_stuff(nx,ny):
    constvec = consts()
    ymax = constvec[4]+constvec[5]
    k = constvec[0]
    tbar = constvec[2]
    kmult = constvec[3]
    xspr = sqrt(2*tbar/k)
    xmax = kmult*xspr

    p,g,dx,dy = mesh.gencylmesh(nx,ny,xmax,ymax)
    dt = find_dt(p[:,0],p[:,1],constvec,dx,dy)
    pars = hstack((constvec,xmax,ymax,dx,dy,dt))
    return p,g,pars

def find_dt(xgrid,ygrid,pars,dx,dy):
    axa = getbig(abs(ax(xgrid,ygrid,pars)))
    aya = getbig(abs(ay(xgrid,ygrid,pars)))
    b = barray(pars)

    maxinv = 4*(b[0,0]/(dx**2)+b[1,1]/(dy**2))+1/(2*dx*dy)*(dx*aya+dy*axa+4*b[0,1])
    dt = (maxinv)**(-1)
    print(dx,dy,dt)
    return dt

def imp_tstep(p,g,pars,ls,psi,ls2 = 0):
    rs = rhs(p,g,pars,psi,ls2)
    pn = spsolve(ls,rs)
    return pn

def rhs(p,g,pars,psi,ls2):
    if sum(ls2**2) == 0:
        rhs = psi
    else:
        rhs = ls2.dot(psi)

    ind = logical_or(g[:,1] == -1,g[:,3] == -1)
    rhs[ind] = 0#(b[0,1]*(psi[g[ind,0]]-psi[g[ind,2]])/(2*dy) + psi[ind]*dxax1[ind])/(b[0,0]-ax1[ind])
    return rhs

def preplhs(p,g,pars,ls_s):
    ls = ls_s.toarray()
    dx = pars[-3]
    dy = pars[-2]
    dt = pars[-1]
    xmax = pars[-5]
    ymax = pars[-6]
    n = p[:,0].size
    ax1 = ax(p[:,0],p[:,1],pars)
    ay1 = ay(p[:,0],p[:,1],pars)
    div = diva(p[:,0],p[:,1],pars)
    b = barray(pars)
 
    for i in range(n):
        if g[i,1] == -1:
            ls[i,:] = 0
            ls[i,i] = -3/(2*dx)*(b[0,0])
            ls[i,g[i,3]] = 2/dx*(b[0,0])
            ls[i,g[g[i,3],3]] = -1/(2*dx)*(b[0,0])
                       
            ls[i,i] = ls[i,i] - ax1[i]
            ls[i,g[i,0]] = b[0,1]/(2*dy)
            ls[i,g[i,2]] = -b[0,1]/(2*dy)
        elif g[i,3] == -1:
            ls[i,:] = 0
            ls[i,i] = 3/(2*dx)*(b[0,0])
            ls[i,g[i,1]] = -2/dx*(b[0,0])
            ls[i,g[g[i,1],1]] = 1/(2*dx)*(b[0,0])

            ls[i,i] = ls[i,i] - ax1[i]
            ls[i,g[i,0]] = b[0,1]/(2*dy)
            ls[i,g[i,2]] = -b[0,1]/(2*dy)
    return spr.csr_matrix(ls)

def lhs(p,g,pars,ifcrank = 0):
    dx = pars[-3]
    dy = pars[-2]
    dt = pars[-1]
    xmax = pars[-5]
    ymax = pars[-6]
    n = p[:,0].size
    ax1 = ax(p[:,0],p[:,1],pars)
    ay1 = ay(p[:,0],p[:,1],pars)
    div = diva(p[:,0],p[:,1],pars)
    b = barray(pars)
    ls = zeros((n,n))

    if ifcrank == 1:
        dtc = dt/2
    elif ifcrank == 2:
        dtc = -dt/2
    else:
        dtc = dt

    d = zeros((n))
    d = d + 1 + dtc*div + 2*dtc*b[0,0]/(dx**2) + 2*dtc*b[1,1]/(dy**2)
    fill_diagonal(ls,d)
    for i in range(n):
        ls[i,g[i,0]] = dtc*ay1[i]/(2*dy)-dtc*b[1,1]/(dy**2)
        ls[i,g[i,2]] = -dtc*ay1[i]/(2*dy)-dtc*b[1,1]/(dy**2)
        ls[i,g[i,3]] = dtc*ax1[i]/(2*dx)-dtc*b[0,0]/(dx**2)
        ls[i,g[i,1]] = -dtc*ax1[i]/(2*dx)-dtc*b[0,0]/(dx**2)
        ls[i,g[g[i,3],0]] = dtc*b[0,1]/(2*dx*dy)
        ls[i,g[g[i,1],2]] = dtc*b[0,1]/(2*dx*dy)
        ls[i,g[g[i,3],2]] = -dtc*b[0,1]/(2*dx*dy)
        ls[i,g[g[i,1],0]] = -dtc*b[0,1]/(2*dx*dy)
    sls = spr.csr_matrix(ls)
    return sls
                    
def barray(pars):
    tbar = pars[2]
    tau = pars[1]
    return 2*tbar*array([[4,tau/2],[tau/2,1]])

def consts():
    k = 1
    tau = 0
    tbar = 1
    kmult = 3
    lsmall = 1
    lbig = 9
    umax = 1
    
    vec = array([k,tau,tbar,kmult,lsmall,lbig,umax])
    return vec

def diva(x,y,pars):
    return -2*pars[0] + upp(y-x/2) + upp(y+x/2)

def dxax(x,y,pars):
    return 2*pars[0]*ones(x.shape) + 0.5*(upp(y-x/2)+upp(y+x/2))
def dyay(x,y,pars):
    return 0.5*(upp(y-x/2)+upp(y+x/2))

def upp(x):
    return 0

def ax(x,y,pars):
    k = pars[0]
    um = uprime(y-x/2,pars)
    up = uprime(y+x/2,pars)

    return -2*k*x+up-um

def ay(x,y,pars):
    k = pars[0]
    um = uprime(y-x/2,pars)
    up = uprime(y+x/2,pars)

    return 0.5*(up+um)

def uprime(x,pars):
    lsmall = pars[4]
    lbig = pars[5]
    l = lsmall + lbig
    umax = pars[6]
    
    uprime = zeros(x.shape)
    ind = (x%l)<lsmall
    uprime[ind] = umax/lsmall
    uprime[~ind] = -umax/lbig

    return uprime

def getbig(array):
    b = reshape(array,array.size)
    sort(b)
    return b[-1]

def specder(f,xmin,xmax,n,ax = 0):
    if ax == 1:
        f = f.T
    fp = fft.fftshift(fft.fft(f))
    k = fft.fftshift(fft.fftfreq(n))
    k = 2*pi*k/(xmax-xmin)
    fp = fp*1j*k
    fpout = n*real(fft.ifft(fft.ifftshift(fp)))
    
    if ax == 1:
        return fpout.T
    else:
        return fpout

def spec2der(f,xmin,xmax,n,ax = 0):
    return specder(specder(f,xmin,xmax,n,ax),xmin,xmax,n,ax)
