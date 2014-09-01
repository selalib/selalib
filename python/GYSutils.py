from matplotlib.pylab import *
import numpy as np
from matplotlib.image import NonUniformImage
import h5py as h5


#--------------------------------------------------
# Load an HDF5 file
#--------------------------------------------------
class loadHDF5():
    """ Load an HDF5 file"""

    def __init__(self, filename):

        fh5       = h5.File(filename, 'r')
        var_in_h5 = fh5.keys()
        self.keys = var_in_h5
        for str_var in var_in_h5:
            #--> replace '%' by '_' in variable names
            str_var_new = str_var.replace('%','_')
            var_tmp     = fh5[str_var]
            sh          = shape(var_tmp)
            self.__dict__[str_var_new] = var_tmp[:].reshape(sh)
        fh5.close()                
#end class loadHDF5


#--------------------------------------------------
# import or reload a module already existent
#  Input:
#    - modulename (str) : name of the module
#--------------------------------------------------
def impall(modulename):
    """ Import or reload a module already existent"""

    exec "exist_module = '%s' in sys.modules.keys()" % (modulename)
    if (exist_module):
        exec "import %s" % (modulename)
        exec "reload(%s)" % (modulename)
    else:
        exec "import %s" % (modulename)
    exec "module_load = %s" % (modulename)
    return module_load
#end def impall


#-------------------------------------------------
#--> Find negative values of 'varname' variable
#     in 'filename' file
#-------------------------------------------------
def Find_NegativValues(filename, varname):
    """
    Find negative values of 'varname' variable
    in 'filename' file
    """

    fh5 = loadHDF5(filename)
    exec "H5f_var = fh5.%s" % (varname)
    findneg_val = (H5f_var < 0.).nonzero()
    nbneg_val   = np.shape(findneg_val)[1]
    return [nbneg_val, findneg_val]
#end def Find_NegativValues


#---------------------------------------------
# Ask the initial and the final times
#---------------------------------------------
def Ask_firstandlast_time(Ttime, ask_diagstep=True, msg=''):
    """ Ask the initial and the final times"""

    nb_time = np.size(Ttime)
    tinit   = raw_input(' init time ' + msg +
                        ' (between ' + str(Ttime[0]) +
                        ' and ' + str(Ttime[nb_time - 1]) +
                        ' (default ' + str(Ttime[0]) + ')) ? : ')
    if (tinit != ''):
        t_init = str2num(tinit)
    else:
        t_init = Ttime[0]
    it1 = search_pos(t_init, Ttime)

    tfin = raw_input(' end time ' + msg +
                     ' (between ' + str(Ttime[it1]) +
                     ' and ' + str(Ttime[nb_time - 1]) +
                     '  (default ' + str(Ttime[nb_time - 1]) + ')) ? : ')
    if (tfin != ''):
        t_fin = str2num(tfin)
    else:
        t_fin = Ttime[nb_time - 1]
    it2 = search_pos(t_fin, Ttime)

    if (ask_diagstep):
        if (it1 != it2):
            time_step = raw_input(' --> Diag step (between 1 and ' +
                                  str(nb_time / 3) + ') ? (default 1) : ')
            if (time_step != ''):
                istep = str2num(time_step)
            else:
                istep = 1
        else:
            istep = 1
        return [it1, it2, istep]
    else:
        return [it1, it2]
#end def Ask_firstandlast_time


#-----------------------------------------------------------------
# Ask for a specific time
#  -> Default value is given in function of 'give_default' value,
#      . give_default = 0 : no default value
#      . give_default = 1 : first time by default
#      . give_default = 2 : last time by default
#-----------------------------------------------------------------
def Ask_time(Ttime, ask_msg='', first_indx_time=0, give_default=0):
    """
    Ask for a specific time
    -> Default value is given in function of 'give_default' value,
    . give_default = 0 : no default value
    . give_default = 1 : first time by default
    . give_default = 2 : last time by default
    """

    nb_time = np.size(Ttime)
    if (ask_msg == ''):
        ask_msg = ' Time'

    if (give_default == 0):
        str_time_ask = raw_input(ask_msg +
                                 ' (between ' + str(Ttime[first_indx_time]) +
                                 ' and ' + str(Ttime[nb_time - 1]) + ') ? : ')
        if (str_time_ask != ''):
            time_ask     = str2num(str_time_ask)
            int_time_ask = search_pos(time_ask, Ttime)
        else:
            int_time_ask = nb_time + 10
    #end if

    if (give_default > 0):
        if (give_default == 1):
            indx_default = 0
        if (give_default == 2):
            indx_default = nb_time - 1
        str_time_ask = raw_input(ask_msg +
                                 ' (between ' + str(Ttime[first_indx_time]) +
                                 ' and ' + str(Ttime[nb_time - 1]) +
                                 ' [default ' + str(Ttime[indx_default]) +
                                 ']) ? : ')
        if (str_time_ask != ''):
            time_ask     = str2num(str_time_ask)
            int_time_ask = search_pos(time_ask, Ttime)
        else:
            if (give_default == 1):
                time_ask = Ttime[first_indx_time]
            if (give_default == 2):
                time_ask = Ttime[nb_time - 1]
            int_time_ask = search_pos(time_ask, Ttime)
            str_time_ask = str(time_ask)
        #end if
    #end if

    return [str_time_ask, int_time_ask]
#end def Ask_time


#--------------------------------------------------
# tw is equal to :
#  True if the variables are keeped in the
#  dictionarie (default)
#  False if the variables are defined as global
#--------------------------------------------------
def DctCorr(dct, tw=True):
    """
    Check if variable is in the dictionarie or global
    tw is equal to :
    True if the variables are keeped in the dictionarie (default)
    False if the variables are defined as global
    """

    dctret = {}
    for k in dct.keys():
        idx = k.find('\x00')
        if idx != -1:
            k2 = k[:idx]
        else:
            k2 = k
        dctret[k2] = dct[k]
    if tw:
        glo = globals()
        for k in dctret.keys():
            glo[k] = dctret[k]
    return dctret
#end def DctCorr


#--------------------------------------------------
# definition of the size of the axes
#--------------------------------------------------
def my_axes(ax, fs=16):
    """ Definition of the size of the axes"""

    xticklabels = get(ax, 'xticklabels')
    yticklabels = get(ax, 'yticklabels')
    setp(xticklabels, 'color', 'k', fontsize=fs)
    setp(yticklabels, 'color', 'k', fontsize=fs)
#end def my_axes


#--------------------------------------------------
# definition of the size of the legend
#--------------------------------------------------
def my_legend(fs=16, visibl=True, lo='best'):
    """ Definition of the size of the legend"""

    leg   = legend(loc=lo)
    ltext = leg.get_texts()
    setp(ltext, 'fontsize', fs)
    lframe = leg.get_frame()
    setp(lframe, fc='w', ec='w', visible=visibl)
    return leg
#end def my_legend


#--------------------------------------------------
# definition of the size of xlabel
#--------------------------------------------------
def my_xlabel(lab="", s=18, tex=0):
    """ Definition of the size of xlabel"""

    if (tex == 0):
        xlab = xlabel(r'$' + lab + '$', size=s)
    else:
        xlab = xlabel(lab, size=s)
    return xlab
#end def my_xlabel


#--------------------------------------------------
# definition of the size of ylabel
#--------------------------------------------------
def my_ylabel(lab="", s=18, tex=0):
    """ Definition of the size of ylabel"""

    if (tex == 0):
        ylab = ylabel(r'$' + lab + '$', size=s)
    else:
        ylab = ylabel(lab, size=s)
    return ylab
#end def my_ylabel


#--------------------------------------------------
# definition of the size of title
#--------------------------------------------------
def my_title(lab="", s=16, tex=0):
    """ Definition of the size of title"""

    if (tex == 0):
        tit = title(r'$' + lab + '$', size=s)
    else:
        tit = title(lab, size=s)
    return tit
#end def my_title


#--------------------------------------------------
# Title time
#--------------------------------------------------
def my_suptitle(lab="", s=18, tex=0):
    """ Title time"""

    if (tex == 0):
        tit = suptitle(r'$' + lab + '$', size=s)
    else:
        tit = suptitle(lab, size=s)
    return tit
#end def title_time


#--------------------------------------------------
# definition of the Gray map
#--------------------------------------------------
def my_graymap():
    """ Definition of the Gray map"""

    graymap = get_cmap('gray')
    graymap._segmentdata = {'blue': ((0.0, 1, 1), (1.0, 0, 0)),
                            'green': ((0.0, 1, 1), (1.0, 0, 0)),
                            'red': ((0.0, 1, 1), (1.0, 0, 0))}
    return graymap
#end def my_graymap


#--------------------------------------------------
# definition of the format of the colorbar
#--------------------------------------------------
def my_colorbar():
    """ Definition of the format of the colorbar"""

    # 2.041
    #colorbar(format='%.3f')
    # 2.0e+00
    colorbar(format='%1.1e')
    # 2.0e+00
    #colorbar(format='%24.12e')
#end def my_colorbar


#--------------------------------------------------
# definition of the close('all')
#--------------------------------------------------
def closa():
    """ Definition of the close('all')"""

    close('all')
#end def closa


#--------------------------------------------------
# search index position
#--------------------------------------------------
def search_pos(x_value, x_array):
    """ Search index position"""

    x_indx = find(x_array > x_value)

    if (size(x_indx) != 0):
        x_indx = x_indx[0] - 1
    else:
        x_indx = size(x_array) - 1

    return x_indx
#end def search_pos


#-------------------------------------------------
# First Derivative
#   Input: F        = function to be derivate
#          dx       = step of the variable for derivative
#          periodic = 1 if F is periodic
#                   = 0 otherwise (by default)
#   Output: dFdx = first derivative of F
#-------------------------------------------------
def Derivee1(F, dx, periodic=0):
    """
    First Derivative
       Input: F        = function to be derivate
              dx       = step of the variable for derivative
              periodic = 1 if F is periodic
       Output: dFdx = first derivative of F
    """

    nx   = len(F)
    dFdx = 0. * F

    c0 = 2. / 3.
    dFdx[2:nx - 2] = c0 / dx * (F[3:nx - 1] - F[1:nx - 3] - (F[4:nx] -
                                                             F[0:nx - 4]) / 8.)

    c1 = 4. / 3.
    c2 = 25. / 12.
    c3 = 5. / 6.
    if (periodic == 0):
        dFdx[0]      = (-F[4] / 4. + c1 * F[3] - 3. * F[2] + 4. * F[1] -
                        c2 * F[0]) / dx
        dFdx[nx - 1] = (F[nx - 5] / 4. - c1 * F[nx - 4] + 3. *
                        F[nx - 3] - 4. * F[nx - 2] + c2 * F[nx - 1]) / dx
        dFdx[1]      = (F[4] / 12. - F[3] / 2. + F[2] / c0 - c3 *
                        F[1] - F[0] / 4.) / dx
        dFdx[nx - 2] = (F[nx - 1] / 4. + c3 * F[nx - 2] - F[nx - 3] /
                        c0 + F[nx - 4] / 2. - F[nx - 5] / 12.) / dx
    else:
        dFdx[0]      = c0 / dx * (F[1] - F[nx - 1] - (F[2] - F[nx - 2]) / 8.)
        dFdx[1]      = c0 / dx * (F[2] - F[0] - (F[3] - F[nx - 1]) / 8.)
        dFdx[nx - 1] = c0 / dx * (F[0] - F[nx - 2] - (F[1] - F[nx - 3]) / 8.)
        dFdx[nx - 2] = c0 / dx * (F[nx - 1] - F[nx - 3] -
                                  (F[0] - F[nx - 4]) / 8.)

    return dFdx
#end def Derivee1


#-----------------------------------------------------------
# Second Derivative
#   Input: F        = function to be derivate
#          dx       = step of the variable for derivative
#          periodic = 1 if F is periodic
#                   = 0 otherwise (by default)
#   Output: dFdx = second derivative of F
#------------------------------------------------------------
def Derivee2(F, dx, periodic=0):
    """
    Second Derivative
      Input: F        = function to be derivate
             dx       = step of the variable for derivative
             periodic = 1 if F is periodic
                      = 0 otherwise (by default)
      Output: dFdx = second derivative of F
    """

    nx  = len(F)
    dx2 = dx * dx

    d2Fdx = 0. * F
    d2Fdx[2:nx - 2] = (-30. * F[2:nx - 2] + 16. * (F[3:nx - 1] + F[1:nx - 3]) -
                       (F[0:nx - 4] + F[4:nx])) / (12. * dx2)

    c1 = 11. / 12.
    c2 = 14. / 3.
    c3 = 9.5
    c4 = 26. / 3.
    c5 = 35. / 12.
    c6 = 5. / 3.
    c7 = 11. / 12.
    if (periodic == 0):
        d2Fdx[0]      = (c1 * F[4] - c2 * F[3] + c3 * F[2] -
                         c4 * F[1] + c5 * F[0]) / dx2
        d2Fdx[nx - 1] = (c1 * F[nx - 5] - c2 * F[nx - 4] + c3 * F[nx - 3] -
                         c4 * F[nx - 2] + c5 * F[nx - 1]) / dx2
        d2Fdx[1]      = (-F[4] / 12. + F[3] / 3. + F[2] / 2. - c6 * F[1] +
                         c7 * F[0]) / dx2
        d2Fdx[nx - 2] = (-F[nx - 5] / 12. + F[nx - 4] / 3. + F[nx - 3] / 2. -
                         c6 * F[nx - 2] + c7 * F[nx - 1]) / dx2
    else:
        d2Fdx[0]      = (-30. * F[0] + 16. * (F[1] + F[nx - 1]) -
                         (F[nx - 2] + F[2])) / (12. * dx2)
        d2Fdx[1]      = (-30. * F[1] + 16. * (F[2] + F[0]) -
                         (F[nx - 1] + F[3])) / (12. * dx2)
        d2Fdx[nx - 1] = (-30. * F[nx - 1] + 16. * (F[0] + F[nx - 2]) -
                         (F[nx - 3] + F[1])) / (12. * dx2)
        d2Fdx[nx - 2] = (-30. * F[nx - 2] + 16. * (F[nx - 1] + F[nx - 3]) -
                         (F[nx - 4] + F[0])) / (12. * dx2)
    #end if

    return d2Fdx
#end def Derivee2


#-------------------------------------------------
# Personal FFT1D function
#-------------------------------------------------
def Fourier1D(F0, x0):
    """ Personal FFT1D function"""

    nx0 = len(x0)
    nx  = 2 * int(nx0 / 2)
    hnx = nx / 2

    dim = np.ndim(F0)
    x   = x0[0:nx]
    if (dim == 1):
        F = F0[0:nx]
    elif (dim == 2):
        F = F0[:, 0:nx]
    else:
        print '\033[0;31m dim not supported !! \033[1;m '

    Lx   = x[nx - 1] - x[0]
    dx   = x[1] - x[0]
    dkx  = 2. * np.pi / (Lx + dx)
    kx   = np.zeros(nx)
    temp = -dkx * np.r_[1:hnx + 1]
    kx[0:hnx]  = temp[::-1]
    kx[hnx:nx] = dkx * np.r_[0:hnx]

    var = np.conjugate(np.fft.fft(np.conjugate(F))) / float(nx)
    TFF = 0. * var

    if (dim == 1):
        TFF[0:hnx]  = var[hnx:nx]
        TFF[hnx:nx] = var[0:hnx]
    elif (dim == 2):
        TFF[:, 0:hnx]  = var[:, hnx:nx]
        TFF[:, hnx:nx] = var[:, 0:hnx]

    return TFF, kx
#end def Fourier1D


#-------------------------------------------------
# Personal FFT2D function
#-------------------------------------------------
def Fourier2D(F0, y0, x0):
    """ Personal FFT2D function"""

    nx0 = len(x0)
    nx  = 2 * int(nx0 / 2)
    hnx = nx / 2
    ny0 = len(y0)
    ny  = 2 * int(ny0 / 2)
    hny = ny / 2

    x = x0[0:nx]
    y = y0[0:ny]
    F = F0[0:ny, 0:nx]

    Lx   = x[nx - 1] - x[0]
    dx   = x[1] - x[0]
    dkx  = 2. * pi / (Lx + dx)
    kx   = np.zeros(nx)
    temp = -dkx * r_[1:hnx + 1]
    kx[0:hnx]  = temp[::-1]
    kx[hnx:nx] = dkx * r_[0:hnx]

    Ly   = y[ny - 1] - y[0]
    dy   = y[1] - y[0]
    dky  = 2. * pi / (Ly + dy)
    ky   = np.zeros(ny)
    temp = -dky * r_[1:hny + 1]
    ky[0:hny]  = temp[::-1]
    ky[hny:ny] = dky * r_[0:hny]

    TFF = np.zeros((ny, nx), dtype=complex)
    AA  = np.zeros((ny, nx), dtype=complex)
    var = np.conjugate(np.fft.fft2(np.conjugate(F))) / float((nx * ny))

    AA[:, 0:hnx]   = var[:, hnx:nx]
    AA[:, hnx:nx]  = var[:, 0:hnx]
    TFF[0:hny, :]  = AA[hny:ny, :]
    TFF[hny:ny, :] = AA[0:hny, :]

    return TFF, kx, ky
#end def Fourier2D


#-------------------------------------------------
# Personal inverse FFT1D function
#-------------------------------------------------
def iFourier1D(TFF, kx):
    """ Personal inverse FFT1D function"""

    nx  = len(kx)
    hnx = nx / 2

    dkx = abs(kx[1] - kx[0])
    x   = r_[0:nx] * 2. * pi / (nx * dkx)

    sTFF         = np.zeros(nx, dtype=complex)
    sTFF[0:hnx]  = TFF[hnx:nx]
    sTFF[hnx:nx] = TFF[0:hnx]

    F = real(np.conjugate(np.fft.ifft(nx * np.conjugate(sTFF))))

    return F, x
#end def iFourier1D


#-----------------------------------------------------
#   find the minimum and maximum indices in a array
#   corresponding to the minimal and maximal values
#   asked
#-----------------------------------------------------
def find_min_max(x_array, x_name):
    """
    Find the minimum and maximum indices in a array
    corresponding to the minimal and maximal values
    asked
    """

    xinf_str = x_array[0]
    xsup_str = x_array[size(x_array) - 1]
    xmin_str = raw_input('minimum ' + x_name +
                         ' (' + repr(xinf_str) + ' to ' +
                         repr(xsup_str) + ') [default ' +
                         repr(xinf_str) + '] ? : ')
    xmax_str = raw_input('maximum ' + x_name +
                         ' (' + repr(xinf_str) + ' to ' + repr(xsup_str) +
                         ') [default ' + repr(xsup_str) + '] ? : ')

    if (xmin_str == ''):
        indx_min = 0
    else:
        xmin     = float(xmin_str)
        indx_min = find(x_array <= xmin)
        indx_min = indx_min[size(indx_min) - 1]

    if (xmax_str == ''):
        indx_max = size(x_array) - 1
    else:
        xmax     = float(xmax_str)
        indx_max = find(x_array >= xmax)
        indx_max = indx_max[0]

        return indx_min, indx_max
#end def find_min_max


#-----------------------------------------------------
#   find the path for executable file
#-----------------------------------------------------
def my_execfile(filename):
    """ Find the path for executable file"""

    import os
    my_path     = os.environ.get('PYTHONPATH')
    directories = my_path.split(os.pathsep)
    ifind = 0
    for idir in directories:
        filesearch = idir + '/' + filename
        if os.path.isfile(filesearch):
            ifind = 1
            execfile(filesearch)
    if ifind == 0:
        print '\033[0;31m filename not found !! \033[1;m '
#end def my_execfile


#-----------------------------------------------------
# Print the dictionnary
#-----------------------------------------------------
def printDict(di, format=" %-25s %s"):
    """ Print the dictionnary"""

    for (key, val) in di.items():
        print format % (str(key) + ':', val)
#end def printDict


#-----------------------------------------------------
#   construct the file name according to the number
#    ifile
#-----------------------------------------------------
def create_filename(acronym, ifile):
    """ Construct the file name according
    to the number ifile"""

    if ifile < 10:
        file_name = acronym + '_000' + repr(ifile) + '.mat'
    elif ifile < 100:
        file_name = acronym + '_00' + repr(ifile) + '.mat'
    elif ifile < 1000:
        file_name = acronym + '_0' + repr(ifile) + '.mat'
    else:
        file_name = acronym + '_' + repr(ifile) + '.mat'
    return file_name
#end def create_filename


#-----------------------------------------------------
# Test for non-uniform images
#-----------------------------------------------------
def nonuniform_imshow(x, y, C, **kwargs):
    """
    Plot image with nonuniform pixel spacing.
    This function is a convenience method for
    calling image.NonUniformImage.
    """
    ax = plt.gca()
    im = NonUniformImage(ax, **kwargs)
    im.set_data(x, y, C)
    ax.images.append(im)
    return im
#end def nonuniform_imshow


#-----------------------------------------------------
# To convert string to integer or float
#-----------------------------------------------------
def str2num(datum):
    """ Convert string to integer or float"""

    try:
        return int(datum)
    except:
        try:
            return float(datum)
        except:
            return datum
#end def str2num


#-----------------------------------------------------
# Write in a file
#-----------------------------------------------------
def write_log(file_log, s):
    """ Write in a file"""

    F = open(file_log, 'a')

    if (not os.path.exists(file_log)):
        print str(s)
    else:
        F.write(s + str("\n"))

    F.close()
#end def write_log


#-----------------------------------------------------
#
#-----------------------------------------------------
def poloidal_plot(r, t, y, conts=[]):
    """ Plolidal plot"""

    xx = dot(cos(t.reshape(len(t), 1)), r.reshape(1, len(r)))
    yy = dot(sin(-t.reshape(len(t), 1)), r.reshape(1, len(r)))
    y_remap = copy(y.reshape((len(t), len(r))))
    figure(figsize=(8, 6))
    if len(conts) > 0:
        contour(xx, yy, y_remap, conts, colors='w')
    pcolormesh(xx, yy, y_remap, cmap='spectral')
    axis('equal')
    ax = gca()
    ax.set_xticks([])
    ax.set_yticks([])
    axis('tight')
    colorbar()
#end def poloidal_plot
