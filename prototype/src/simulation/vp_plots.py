from pylab import *

rundir='run_k=0.26/'

x=loadtxt(rundir+'x.txt')
v=loadtxt(rundir+'v.txt')
v_df = loadtxt(rundir+'deltaf_movie_v.txt')
t=loadtxt(rundir+'rho_movie_t.txt')
rho=loadtxt(rundir+'rho_movie.txt')
#rhok = loadtxt(rundir+'rho_from_E_movie.txt')
thdiag = loadtxt(rundir+'time_history.txt')
phi=loadtxt(rundir+'phi_movie.txt')
efield = loadtxt(rundir+'efield_movie.txt')
edr = loadtxt(rundir+"Edr.txt")
deltaf=loadtxt(rundir+'deltaf_movie.txt')
nx=x.size
nv=v.size
deltaf.shape = (deltaf.shape[0]/nx,nx,deltaf.shape[1])
fhat0v = sum(deltaf,axis=1)

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
    title('kinetic entropy')
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
    adr = loadtxt(rundir+"adr.txt")
    plot(adr[:,0],adr[:,1])

def plot_avg(x, u, tstart=0):
    """plot the average over time (starting at tstart) of function u"""
    uavg = sum(u[tstart:,:],axis=0)
    plot(x,uavg)
    title("average")
    draw()
    
def movie_rho(maxit=-1):
    t=loadtxt(rundir+'rho_movie_t.txt')
    rho=loadtxt(rundir+'rho_movie.txt')
    umin = rho.min()
    umax = rho.max()
    if maxit < 0:
        imax = t.size
    else:
        imax = maxit
    for i in range(imax):        
        plot(x,rho[i,:])
        ylim(umin,umax)
        title('rho at time '+ str(t[i]))
        draw()

def plot1d(u,it=0, name="u"):
    """Plots 1d function u versus x at time t[it]."""
    plot(x,u[it])
    title(name + ' at time '+ str(t[it]))
    
def movie1d(x,u,name='movie',imin=0, imax=201):
    umin = u.min()
    umax = u.max()
    for i in range(imin,imax):        
        plot(x,u[i,:])
        ylim(umin,umax)
        draw()

def plot2d(u,it=0, name="u"):
    """Plots 2 function at time t[it]."""
    imshow(u[it,:,:].transpose(), origin='lower', aspect=4.0,
           extent=(x[0],x[-1], v_df[0],v_df[-1]))
    colorbar()
    title(name + ' at time '+ str(t[it]))
    
def movie2d(u):
    hold(False)
    for i in range(u.shape[0]/2):
        clf()
        title('frame='+str(i))
        imshow(u[i,:,:].transpose(), origin='lower', aspect=4.0,
               extent=(x[0],x[-1], v_df[0],v_df[-1]))
        colorbar()
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
    uk = rfft(u)/nx
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
    uk = rfft(u)/nx
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
        

def plot_modes_frequency(u,istart=0,iend=t.size-1,f=abs):
    """Plot modes versus frequency using the butterfly transform in time
    of the Fourier modes starting at time t[itstart].
    The butterfly of a function u defined in [t1,t2] consists in extending
    it to [t1-t2,t2] by parity. It is then considered periodic on this interval."""

    # compute frequency range
    freq = arange(t[istart:iend].size+1)*pi/(t[iend]-t[istart])
    # compute Fourier transform of u
    uk = rfft(u)/nx
    # compute the butterfly transform of the first 7 modes of uk (excluding 0)
    s = uk.shape
    uk_butterfly = empty((8,2*s[0]-1))
    for i in range(8):
        uk_butterfly[i,s[0]-1:] = uk[:,i] 
        uk_butterfly[i,:s[0]-1] = uk[:0:-1,i]
    # Fourier transform of uk_butterfly
    ukom = fft(uk_butterfly)/(iend-istart+1)
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
