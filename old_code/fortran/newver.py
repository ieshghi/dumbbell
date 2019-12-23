from numpy import *

def runsim(nx,ny,nt):
    xg,yg,pars = init_stuff(nx,ny)
    pi = pinit(xg,yg)
    p = pi.copy()
    p_hist = zeros((nx,ny,nt//1))
    j = 0
    p_hist[:,:,0] = p
    for i in range(nt):
        p = p + ftcs(p,pars,xg,yg)
        print(sum(abs(p)))
        if i%1==0:
            p_hist[:,:,j] = p
            j = j+1
    return p_hist,pi,p

def pinit(xg,yg):
    p = exp(-(4*(xg)**2+(yg-5)**2))
    return p/sum(p)

def init_stuff(nx,ny):
    constvec = consts()
    xgrid,ygrid,dx,dy,xmax,ymax = gengrid(nx,ny,constvec)
    dt = find_dt(xgrid,ygrid,constvec,dx,dy)
    pars = hstack((constvec,xmax,ymax,dx,dy,dt))
    return xgrid,ygrid,pars

def gengrid_fft(nx,ny,pars):
    ymax = pars[4]+pars[5]
    k = pars[0]
    tbar = pars[2]
    kmult = pars[3]

    xspr = sqrt(2*tbar/k)
    xmax = kmult*xspr

    xv = linspace(-xmax,xmax,nx,endpoint = False) # This might be problematic for FFT...
    yv = linspace(0,ymax,ny,endpoint = False)
    xg,yg = meshgrid(xv,yv)
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]

    return xg,yg,dx,dy,xmax,ymax


def gengrid(nx,ny,pars):
    ymax = pars[4]+pars[5]
    k = pars[0]
    tbar = pars[2]
    kmult = pars[3]

    xspr = sqrt(2*tbar/k)
    xmax = kmult*xspr

    xv = linspace(-xmax,xmax,nx)
    yv = linspace(0,ymax,ny)
    xg,yg = meshgrid(xv,yv)
    dx = xv[1]-xv[0]
    dy = yv[1]-yv[0]

    return xg,yg,dx,dy,xmax,ymax

def find_dt(xgrid,ygrid,pars,dx,dy):
    axa = getbig(abs(ax(xgrid,ygrid,pars)))
    aya = getbig(abs(ay(xgrid,ygrid,pars)))
    b = barray(pars)

    maxinv = 4*(b[0,0]/(dx**2)+b[1,1]/(dy**2))+1/(2*dx*dy)*(dx*aya+dy*axa+4*b[0,1])
    dt = 1/10*(maxinv)**(-1)
    print(dx,dy,dt)
    return dt

def ftcs_fft(p,pars,xgrid,ygrid): #assume periodic BC in y-directions, take spectral ders.
    dx = pars[-3]
    dy = pars[-2]
    dt = pars[-1]
    xmax = pars[-5]
    ymax = pars[-6]
    nx = size(p[:,1])
    ny = size(p[1,:])
    axa = ax(xgrid,ygrid,pars)
    aya = ay(xgrid,ygrid,pars)
    pax = p*axa
    pay = p*aya
    ymin = 0
    xmin = -xmax
    b = barray(pars)
    pnew = zeros(p.shape)
    pnew = -(specder(pay,ymin,ymax,ny,ax=1) + specder(pax,xmin,xmax,nx,ax=0))
    pnew = pnew + b[0,0]*spec2der(p,xmin,xmax,nx,ax=0) + b[1,1]*spec2der(p,ymin,ymax,ny,ax=1)
    pnew = pnew + 2*b[0,1]*specder(specder(p,xmin,xmax,nx,ax=0),ymin,ymax,ny,ax=1)
    return dt*pnew

def ftcs(p,pars,xgrid,ygrid):
    dx = pars[-3]
    dy = pars[-2]
    dt = pars[-1]
    nx = size(p[:,1])
    ny = size(p[1,:])
    axa = ax(xgrid,ygrid,pars)
    aya = ay(xgrid,ygrid,pars)
    b = barray(pars)
    p = init_ghost(p,axa,aya,b,dx,dy)
    pax = p*axa
    pay = p*aya
    pnew = zeros(p.shape)
    pnew = -1./2*((roll(pay,1,1)-roll(pay,-1,1))/dy + (roll(pax,1,0)-roll(pax,-1,0))/dx)
    pnew = pnew + b[0,0]*(2*p-roll(p,1,0)-roll(p,-1,0))/(dx**2)
    pnew = pnew + b[1,1]*(2*p-roll(p,1,1)-roll(p,-1,1))/(dy**2)
    pnew = pnew + b[0,1]*(roll(roll(p,1,0),1,1) + roll(roll(p,-1,0),-1,1) - roll(roll(p,-1,0),1,1) - roll(roll(p,1,0),-1,1))/(2*dx*dy)
    return dt*pnew

def init_ghost(p,ax,ay,b,dx,dy,type = 0):
    nx = size(p[:,1])
    ny = size(p[1,:])
    dyp = (roll(p,1,1)-roll(p,-1,1))/(2*dy)
    pg = p.copy()
    b = -b
    if type == 0:
        pg[:,0] = pg[:,-2] #periodic boundaries along y
        pg[:,-1] = pg[:,1]
    pg[0,:] = (ax[2,:]+b[0,0])/(ax[0,:]+b[0,0])*p[2,:] + (2*dx*b[0,1])/(ax[0,:]+b[0,0])*dyp[1,:]
    pg[-1,:] = (ax[-3,:]+b[0,0])/(ax[-1,:]+b[0,0])*p[-3,:] - (2*dx*b[0,1])/(ax[-1,:]+b[0,0])*dyp[-2,:]
    
    return pg
                    
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
    return b[0]

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
