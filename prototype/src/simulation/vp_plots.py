from numpy import *
from pylab import *
import tables
import string
rundir='/Users/sonnen/Codes/selalib/prototype/keen1d02/'

thdiag = loadtxt(rundir+'thdiag.dat')
efield = loadtxt(rundir+'exdiag.dat')
rho = loadtxt(rundir+'rhodiag.dat')
param_file = open(rundir+'param_out.dat')
param = param_file.readline().split()
case = param[0]
xmin = string.atof(param[1])
xmax = string.atof(param[2])
ncx = string.atoi(param[3])
dx = (xmax-xmin) / ncx
vmin = string.atof(param[4])
vmax = string.atof(param[5])
ncv = string.atoi(param[6])
dt = string.atof(param[7])
nbiter = string.atoi(param[8])
freqdiag = string.atoi(param[9])
if param[10]== 0:
    is_delta_f = True
else:
    is_delta_f = False
params_drive_file = open(rundir+'param_out_drive.dat')
params_drive = params_drive_file.readline().split()
omega_dr = string.atof(params_drive[8])
edr = loadtxt(rundir+"eappdiag.dat")
x=linspace(xmin,xmax,ncx+1)
v=linspace(vmin,vmax,ncv+1)
X,V=meshgrid(x,v)
# Maxwellian needs to be added to f in delta_f case as delta_f is stored
if is_delta_f:
    fmaxw = exp(-0.5*V*V)/sqrt(2*pi)
else:
    fmaxw = zeros_like(X)
    
def read_h5(it):
    """reads hdf5 file and returns dist_func array"""
    h5file=tables.openFile(rundir+'dist_func'+string.zfill(it,4)+'-dist_func.h5','r')
    obj = h5file.getNode("/dist_func")
    a = obj.read()
    h5file.close()
    return a

def plot_f(it,ivmin=0,ivmax=-1,fmax=None):
    clf()
    f=read_h5(it)+fmaxw
    fplot = f[ivmin:ivmax,:]
    if fmax == None:
        fmax = fplot.max()
    ratio = (xmax-xmin) / (v[ivmax]-v[ivmin])
    imshow(fplot,origin='lower',extent=(xmin,xmax,v[ivmin],v[ivmax]),
           aspect=ratio,vmax=fmax)
    colorbar()

def movie_f(ivmin=0,ivmax=-1,fmax=None,vphase=0.):
    clf()
    f=read_h5(1)+fmaxw
    fplot = f[ivmin:ivmax,:]
    ratio = (xmax-xmin) / (v[ivmax]-v[ivmin])
    if fmax == None:
        fmax = fplot.max()
    im=imshow(fplot,origin='lower',extent=(xmin,xmax,v[ivmin],v[ivmax]),
           aspect=ratio,vmax=fmax)
    colorbar()
    for it in range(2,nbiter/freqdiag):
        f = read_h5(it)+fmaxw
        f1 = f[ivmin:ivmax,:]
        # phase shift
        if vphase <> 0.0:
            fplot=zeros((f1.shape[0],ncx),dtype=float)
            for i in range(ncx):
                xshift = mod(i*dx + (vphase*freqdiag*(it-1))*dt, xmax-xmin)
                k = int(xshift/dx)
                alpha = xshift - k*dx
                fplot[:,i] = (1.-alpha)*f1[:,k] + alpha*f1[:,k+1]
        else:
            fplot = f1[:,:-1]
        #title("f at time "+str((it-1)*dt))
        im.set_data(fplot)
        draw()

def plot_fhatv(it):
    clf()
    f=read_h5(it)
    f_hatv=sum(f,axis=1)
    plot(v,f_hatv)

def plot_thdiag():
    """Plot time history monitoring data: average, L1 norm, kinetic entropy, L2 norm."""
    fig=figure()
    fig.subplots_adjust(hspace=0.4)
    subplot(221)
    plot(thdiag[:,0], thdiag[:,1])
    title('average')
    subplot(222)
    plot(thdiag[:,0], thdiag[:,2])
    title('L1 norm')
    subplot(223)
    plot(thdiag[:,0], thdiag[:,3])
    title('momentum')
    subplot(224)
    plot(thdiag[:,0], thdiag[:,4])
    title('L2 norm')

def plot_energy():
    fig=figure()
    fig.subplots_adjust(hspace=0.4)
    subplot(311)
    plot(thdiag[:,0], thdiag[:,5])
    title('kinetic energy')
    subplot(312)
    plot(thdiag[:,0], thdiag[:,6])
    title('potential energy')
    subplot(313)
    plot(thdiag[:,0], thdiag[:,7])
    title('total energy')

def plot_E_PF():
    adr = loadtxt(rundir+"adrdiag.dat")
    plot(adr[:,0],adr[:,1])

def plot_avg(x, u, tstart=0):
    """plot the average over time (starting at tstart) of function u"""
    uavg = sum(u[tstart:,:],axis=0)
    plot(x,uavg)
    title("average")
    draw()

def plot1d(u,it=0, name="u"):
    """Plots 1d function u versus x at time t[it]."""
    plot(x,u[it])
    title(name + ' at time '+ str(it*dt))
    
def movie1d(u,name='movie',imin=0, imax=nbiter):
    umin = u.min()
    umax = u.max()
    for i in range(imin,imax):        
        plot(x,u[i])
        #ylim(umin,umax)
        draw()


def phase(a):
    """Compute phase of array of complex number"""
    b=arctan2(a.imag,a.real)
    c=empty_like(b)
    arg = 0.
    c[0] = b[0]
    #c[1] = b[1]
    for i in range(1,a.size-1):
        # When arctan has a large jump increment phase by +- 2*pi
        if b[i]-b[i-1] > pi:
            arg = arg - 2*pi
        elif b[i-1]-b[i] > pi:
            arg = arg - 2*pi
        c[i] = b[i] + arg
    c[-1] = b[-1] + arg
    return c

def plot_phase_modes_time(u,f=phase):
    uk = rfft(u)/u.shape[1]
    t=arange(u.shape[0]) * freqdiag * dt
    clf()
    hold(True)
    plot(t,f(uk[:,1])/t,'r',label='k=1')
    plot(t,f(uk[:,2])/(2*t),'b',label='k=2')
    plot(t,f(uk[:,3])/(3*t),'orange',label='k=3')
    plot(t,f(uk[:,4])/(4*t),'c',label='k=4')
    plot(t,f(uk[:,5])/(5*t),'k',label='k=5')
    plot(t,f(uk[:,6])/(6*t),'m',label='k=6')
    plot(t,f(uk[:,7])/(7*t),'y',label='k=7')
    title('rho_k vs time') 
    legend()
    hold(False)

def plot_modes_time(u,f=abs):
    uk = rfft(u[:,:-1])/(ncx+1)
    t=arange(u.shape[0]) * freqdiag * dt
    clf()
    hold(True)
    plot(t,f(uk[:,1]),'r',label='k=1')
    plot(t,f(uk[:,2]),'b',label='k=2')
    plot(t,f(uk[:,3]),'orange',label='k=3')
    plot(t,f(uk[:,4]),'c',label='k=4')
    plot(t,f(uk[:,5]),'k',label='k=5')
    plot(t,f(uk[:,6]),'m',label='k=6')
    plot(t,f(uk[:,7]),'y',label='k=7')
    title('rho_k vs time') 
    legend()
    hold(False) 
        

def plot_modes_frequency(u,istart=0,iend=-1,f=abs,step=1):
    """Plot modes versus frequency using the butterfly transform in time
    of the Fourier modes starting at time t[itstart].
    The butterfly of a function u defined in [t1,t2] consists in extending
    it to [t1-t2,t2] by parity. It is then considered periodic on this interval."""
    t=arange(u[::step,:].shape[0]) * freqdiag * dt * step
    nx = ncx+1
    # compute frequency range
    freq = arange(t[istart:iend].size+1)*pi/(t[iend]-t[istart])/omega_dr
    # compute Fourier transform of u
    uk = rfft(u[::step,:-1])/nx
    # compute the butterfly transform of the first 7 modes of uk (excluding 0)
    s = uk.shape
    uk_butterfly = empty((8,2*s[0]-1))
    for i in range(8):
        uk_butterfly[i,s[0]-1:] = uk[:,i] 
        uk_butterfly[i,:s[0]-1] = uk[:0:-1,i]
    # Fourier transform of uk_butterfly
    ukom = fft(uk_butterfly)/t.size
    # plot positive half of spectrum
    clf()
    hold(True)
    plot(freq,f(ukom[1,istart:s[0]]),'r',label='k=1')
    plot(freq,f(ukom[2,istart:s[0]]),'b',label='k=2')
    plot(freq,f(ukom[3,istart:s[0]]),'orange',label='k=3')
    plot(freq,f(ukom[4,istart:s[0]]),'c',label='k=4')
    plot(freq,f(ukom[5,istart:s[0]]),'k',label='k=5')
    plot(freq,f(ukom[6,istart:s[0]]),'m',label='k=6')
    plot(freq,f(ukom[7,istart:s[0]]),'y',label='k=7')
    title('rho_k vs frequency') 
    legend()
    hold(False) 

def print_conserved():
    """Prints errors in conserved quantities"""
    print 'mass:', (thdiag[:,1].max()-thdiag[:,1].min())/thdiag[:,1].max()
    print 'L1 norm:',(thdiag[:,2].max()-thdiag[:,2].min())/thdiag[:,2].max()
    print 'Total momentum:',(thdiag[:,3].max()-thdiag[:,3].min())
    print 'Total energy:',(thdiag[:,7].max()-thdiag[:,7].min())/thdiag[:,7].max()
    
