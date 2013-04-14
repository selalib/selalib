from matplotlib.pylab import *
import numpy as np
from matplotlib.image import NonUniformImage
import h5py as h5

#--------------------------------------------------
# Load an HDF5 file
#--------------------------------------------------
class loadHDF5():

    def __init__(self,filename):
        fh5       = h5.File(filename,'r');
        var_in_h5 = fh5.keys()
        self.keys = var_in_h5 
        for str_var in var_in_h5:
            exec "%s = fh5['%s']" %(str_var,str_var)
            exec "sh = shape(%s)" %(str_var)
            exec "self.%s = %s[:].reshape(%s)" %(str_var,str_var,sh)
        fh5.close()                


#---------------------------------------------
# Ask the initial and the final times
#---------------------------------------------
def Ask_time(Ttime,ask_diagstep=True,msg=''):
    from utils_func import str2num

    nb_time = np.size(Ttime)
    tinit = raw_input(' init time ' + msg + \
                      ' (between '+ str(Ttime[0]) + \
                      ' and '+ str(Ttime[nb_time-1]) + \
                      ' (default '+ str(Ttime[0]) + ')) ? : ')
    if (tinit !=''):
        t_init = str2num(tinit)
    else:
        t_init = Ttime[0]
    it1 = search_pos(t_init,Ttime)

    tfin = raw_input(' end time ' + msg + \
                     ' (between ' + str(Ttime[it1]) + \
                     ' and ' + str(Ttime[nb_time-1]) + \
                     '  (default '+ str(Ttime[nb_time-1]) + ')) ? : ')
    if (tfin !=''):
        t_fin = str2num(tfin)
    else:
        t_fin = Ttime[nb_time-1]
    it2 = search_pos(t_fin,Ttime)

    if (ask_diagstep):
        if ( it1!=it2 ):
            time_step = raw_input(' --> Diag step (between 1 and ' + \
                                  str(nb_time/3) + ') ? (default 1) : ')
            if (time_step !=''):
                istep = str2num(time_step)
            else:
                istep = 1
        else:
            istep = 1
        return [it1,it2,istep]
    else:
        return [it1,it2]


#--------------------------------------------------
# tw is equal to :
#  True if the variables are keeped in the
#  dictionarie (default)
#  False if the variables are defined as global
#--------------------------------------------------
def DctCorr (dct,tw=True) :
    dctret = {}
    for k in dct.keys() :
        idx=k.find ('\x00')
        if idx != -1:
            k2=k[:idx]
        else:
            k2=k
        dctret[k2] = dct[k]
    if tw :
        glo=globals()
        for k in dctret.keys() :
            glo[k]=dctret[k]
    return dctret


#--------------------------------------------------
# definition of the size of the axes
#--------------------------------------------------
def my_axes (ax,fs=12) :
    xticklabels = get(ax, 'xticklabels');
    yticklabels = get(ax, 'yticklabels');
    setp(xticklabels, 'color', 'k', fontsize=fs);
    setp(yticklabels, 'color', 'k', fontsize=fs);
    grid(True)


#--------------------------------------------------
# definition of the size of the legend
#--------------------------------------------------
def my_legend (fs=12,visibl=True,lo='best') :
    leg    = legend(loc=lo);
    ltext  = leg.get_texts()
    setp(ltext, 'fontsize', fs)
    lframe = leg.get_frame()
    setp(lframe,fc='w',ec='w',visible=visibl)
    return leg


#--------------------------------------------------
# definition of the size of xlabel
#--------------------------------------------------
def my_xlabel (lab="",s=12,tex=0) :
    if (tex == 0):
        xlab = xlabel(r'$'+lab+'$',size=s)
    else:
        xlab = xlabel(lab,size=s)
    return xlab


#--------------------------------------------------
# definition of the size of xlabel
#--------------------------------------------------
def my_ylabel (lab="",s=12,tex=0) :
    if (tex == 0):
        ylab = ylabel(r'$'+lab+'$',size=s)
    else:
        ylab = ylabel(lab,size=s)
    return ylab


#--------------------------------------------------
# definition of the size of title
#--------------------------------------------------
def my_title (lab="",s=12,tex=0) :
    if (tex == 0):
        tit = title(r'$'+lab+'$',size=s)
    else:
        tit = title(lab,size=s)        
    return tit


#--------------------------------------------------
# definition of the size of title
#--------------------------------------------------
def my_graymap():
    graymap = get_cmap('gray')
    graymap._segmentdata = {'blue': ((0.0,1,1), (1.0,0,0)),
                            'green': ((0.0,1,1), (1.0,0,0)),
                            'red': ((0.0,1,1), (1.0,0,0))}
    return graymap


#--------------------------------------------------
# definition of the format of the colorbar
#--------------------------------------------------
def my_colorbar ():
    # 2.041
    #colorbar(format='%.3f');
    # 2.0e+00
    colorbar(format='%1.1e');
    # 2.0e+00
    #colorbar(format='%24.12e');


#--------------------------------------------------
# definition of the close('all')
#--------------------------------------------------
def closa ():
    close('all')


#--------------------------------------------------
# search index position
#--------------------------------------------------
def search_pos(x_value,x_array) :
    x_indx = find(x_array>x_value);
    
    if (size(x_indx) != 0):
        x_indx = x_indx[0]-1;
    else:
        x_indx = size(x_array)-1;
        
    return x_indx


#-------------------------------------------------
# First Derivative
#   Input: F        = function to be derivate
#          dx       = step of the variable for derivative
#          periodic = 1 if F is periodic
#                   = 0 otherwise (by default)
#   Output: dFdx = first derivative of F
#-------------------------------------------------
def Derivee1(F,dx,periodic=0):
    nx = len(F)

    dFdx = 0.*F

    c0 = 2./3.
    dFdx[2:nx-2] = c0/dx * ( F[3:nx-1]-F[1:nx-3]  - (F[4:nx]-F[0:nx-4])/8. )

    c1 = 4./3.
    c2 = 25./12.
    c3 = 5./6.
    if (periodic == 0):
        dFdx[0]    = (-F[4]/4.    + c1*F[3]    - 3.*F[2]     + 4.*F[1]     - c2*F[0])/dx
        dFdx[nx-1] = ( F[nx-5]/4. - c1*F[nx-4] + 3.*F[nx-3]  - 4.*F[nx-2]  + c2*F[nx-1])/dx
        dFdx[1]    = ( F[4]/12.   - F[3]/2.     + F[2]/c0    - c3*F[1]    - F[0]/4.)/dx
        dFdx[nx-2] = ( F[nx-1]/4. + c3*F[nx-2] - F[nx-3]/c0 + F[nx-4]/2. - F[nx-5]/12.)/dx
    else:
        dFdx[0]    = c0/dx * ( F[1]- F[nx-1]   - (F[2]-F[nx-2])/8. )
        dFdx[1]    = c0/dx * ( F[2]- F[0]    - (F[3]-F[nx-1])/8. )
        dFdx[nx-1] = c0/dx * ( F[0]- F[nx-2] - (F[1]-F[nx-3])/8. )
        dFdx[nx-2] = c0/dx * ( F[nx-1]-F[nx-3] - (F[0]-F[nx-4])/8. )

    return dFdx


#-----------------------------------------------------------
# Second Derivative
#   Input: F        = function to be derivate
#          dx       = step of the variable for derivative
#          periodic = 1 if F is periodic
#                   = 0 otherwise (by default)
#   Output: dFdx = second derivative of F
#------------------------------------------------------------
def Derivee2(F,dx,periodic=0):
    nx  = len(F)
    dx2 = dx*dx

    d2Fdx = 0.*F
    d2Fdx[2:nx-2] = (-30.*F[2:nx-2] + 16.*(F[3:nx-1]+F[1:nx-3]) - \
                (F[0:nx-4]+F[4:nx])) / (12.*dx2);

    c1 = 11./12.;
    c2 = 14./3.;
    c3 = 9.5;
    c4 = 26./3.;
    c5 = 35./12.;
    c6 = 5./3.;
    c7 = 11./12.;
    if (periodic == 0):
       d2Fdx[0]    = (c1*F[4] - c2*F[3] + c3*F[2] - c4*F[1] + c5*F[0])/dx2;
       d2Fdx[nx-1] = (c1*F[nx-5] - c2*F[nx-4] + c3*F[nx-3] - c4*F[nx-2] + c5*F[nx-1])/dx2;
       d2Fdx[1]    = (-F[4]/12. + F[3]/3. + F[2]/2. - c6*F[1] + c7*F[0])/dx2;
       d2Fdx[nx-2] = (-F[nx-5]/12. + F[nx-4]/3. + F[nx-3]/2. - c6*F[nx-2] + c7*F[nx-1])/dx2;
    else:
       d2Fdx[0]    = (-30.*F[0] + 16.*(F[1]+F[nx-1]) - (F[nx-2]+F[2])) / (12.*dx2);
       d2Fdx[1]    = (-30.*F[1] + 16.*(F[2]+F[0]) - (F[nx-1]+F[3])) / (12.*dx2);
       d2Fdx[nx-1] = (-30.*F[nx-1] + 16.*(F[0]+F[nx-2]) - (F[nx-3]+F[1])) / (12.*dx2);
       d2Fdx[nx-2] = (-30.*F[nx-2] + 16.*(F[nx-1]+F[nx-3]) - (F[nx-4]+F[0])) / (12.*dx2);
    #end if

    return d2Fdx


#-------------------------------------------------
# Personal FFT2D function
#-------------------------------------------------
def Fourier2D(F0,y0,x0):
    nx0 = len(x0)
    nx  = 2*int(nx0/2)
    hnx = nx/2
    ny0 = len(y0)
    ny  = 2*int(ny0/2)
    hny = ny/2

    x   = x0[0:nx]
    y   = y0[0:ny]
    F   = F0[0:ny,0:nx]

    Lx  = x[nx-1]-x[0]
    dx  = x[1] - x[0]
    dkx = 2.*pi/(Lx+dx)
    kx  = np.zeros(nx)
    temp       = -dkx*r_[1:hnx+1]
    kx[0:hnx]  = temp[::-1]
    kx[hnx:nx] = dkx*r_[0:hnx]
    
    Ly  = y[ny-1]-y[0]
    dy  = y[1] - y[0]
    dky = 2.*pi/(Ly+dy)
    ky  = np.zeros(ny)
    temp       = -dky*r_[1:hny+1]
    ky[0:hny]  = temp[::-1]
    ky[hny:ny] = dky*r_[0:hny]

    TFF = np.zeros(( ny, nx ), dtype=complex)
    AA  = np.zeros(( ny, nx ), dtype=complex)
    var = fft2(F)/(nx*ny)
    
    AA[:,0:hnx]  = var[:,hnx:nx]
    AA[:,hnx:nx] = var[:,0:hnx]
    TFF[0:hny,:]  = AA[hny:ny,:]
    TFF[hny:ny,:] = AA[0:hny,:]
    
    return TFF,kx,ky


#-------------------------------------------------
# Personal FFT1D function
#-------------------------------------------------
def Fourier1D(F0,x0):
    nx0 = len(x0)
    nx  = 2*int(nx0/2)
    hnx = nx/2

    dim=ndim(F0)

    x   = x0[0:nx]
    if (dim == 1):
        F   = F0[0:nx]
    elif (dim == 2):
        F   = F0[:,0:nx]
    else:
        print 'dim not supported'

    Lx  = x[nx-1]-x[0]
    dx  = x[1] - x[0]
    dkx = 2.*pi/(Lx+dx)
    kx  = np.zeros(nx)
    temp       = -dkx*r_[1:hnx+1]
    kx[0:hnx]  = temp[::-1]
    kx[hnx:nx] = dkx*r_[0:hnx]
    
    var = fft(F)/(nx)
    TFF=0.*var
    
    if (dim == 1):
        TFF[0:hnx]  = var[hnx:nx]
        TFF[hnx:nx] = var[0:hnx]
    elif (dim == 2):
        TFF[:,0:hnx]  = var[:,hnx:nx]
        TFF[:,hnx:nx] = var[:,0:hnx]

    return TFF,kx


#-------------------------------------------------
# Personal inverse FFT1D function
#-------------------------------------------------
def iFourier1D(TFF,kx):
    nx = len(kx)
    hnx = nx/2

    dkx = abs(kx[1]-kx[0])
    x = r_[0:nx]*2.*pi/(nx*dkx)

    sTFF = np.zeros(nx, dtype=complex)
    sTFF[0:hnx] = TFF[hnx:nx]
    sTFF[hnx:nx]  = TFF[0:hnx]

    F = real(ifft(nx*sTFF))
    
    return F,x

#-----------------------------------------------------
#   find the minimum and maximum indices in a array 
#   corresponding to the minimal and maximal values
#   asked
#-----------------------------------------------------
def find_min_max(x_array,x_name):
  xinf_str = x_array[0];
  xsup_str = x_array[size(x_array)-1];
  xmin_str = raw_input('minimum ' + x_name + ' (' + `xinf_str` + ' to ' + `xsup_str` + ') [default ' + `xinf_str` + '] ? : ');
  xmax_str = raw_input('maximum ' + x_name + ' (' + `xinf_str` + ' to ' + `xsup_str` + ') [default ' + `xsup_str` + '] ? : ');

  if (xmin_str==''):
    indx_min = 0;
  else:
    xmin     = float(xmin_str);
    indx_min = find(x_array<=xmin);
    indx_min = indx_min[size(indx_min)-1];

  if (xmax_str==''):
    indx_max = size(x_array)-1;
  else:
    xmax     = float(xmax_str);
    indx_max = find(x_array>=xmax);
    indx_max = indx_max[0];

  return indx_min, indx_max;


#-----------------------------------------------------
#   find the path for executable file
#-----------------------------------------------------
def my_execfile(filename):
   import os
   my_path     = os.environ.get('PYTHONPATH')
   directories = my_path.split(os.pathsep);
   ifind       = 0
   for idir in directories:
      filesearch = idir + '/' +filename
      if os.path.isfile(filesearch):
         ifind = 1
         execfile(filesearch);
   if ifind==0:
      print filename + ' not found';


#-----------------------------------------------------
# Print the dictionnary   
#-----------------------------------------------------
def printDict(di, format="%-25s %s"):
    for (key, val) in di.items():
        print format % (str(key)+':', val)


#-----------------------------------------------------
#   construct the file name according to the number
#    ifile
#-----------------------------------------------------
def create_filename(acronym,ifile):
    if ifile < 10:
        file_name = acronym + '_000' + `ifile` + '.mat';
    elif ifile < 100 :
        file_name = acronym + '_00' + `ifile` + '.mat';
    elif ifile < 1000 :
        file_name = acronym + '_0' + `ifile` + '.mat';
    else:
        file_name = acronym + '_' + `ifile` + '.mat';
    return file_name;


#-----------------------------------------------------
# Test for non-uniform images
#-----------------------------------------------------
def nonuniform_imshow(x, y, C, **kwargs):
    """Plot image with nonuniform pixel spacing.
    
    This function is a convenience method for calling image.NonUniformImage.
    """
    ax = plt.gca()
    im = NonUniformImage(ax, **kwargs)
    im.set_data(x, y, C)
    ax.images.append(im)
    return im


#-----------------------------------------------------
# To convert string to integer or float
#-----------------------------------------------------
def str2num(datum):
    try:
        return int(datum)
    except:
        try:
            return float(datum)
        except:
            return datum

#-----------------------------------------------------
# Write in a file
#-----------------------------------------------------
def write_log(file_log,s):
    F = open(file_log,'a')

    if (not os.path.exists(file_log)):
        print str(s)
    else:
        F.write(s+str("\n"))

    F.close()


#-----------------------------------------------------
#
#-----------------------------------------------------
def poloidal_plot(r,t,y,conts=[]):
    xx      = dot(cos(t.reshape(len(t),1)),r.reshape(1,len(r)))
    yy      = dot(sin(-t.reshape(len(t),1)),r.reshape(1,len(r)))
    y_remap = copy(y.reshape((len(t),len(r))))
    figure(figsize=(8,6))
    if len(conts)>0:
        contour(xx,yy,y_remap,conts,colors='w')
    pcolormesh(xx,yy,y_remap,cmap='spectral')
    axis('equal')
    ax=gca()
    ax.set_xticks([])
    ax.set_yticks([])
    axis('tight')
    colorbar()
