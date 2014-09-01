import sys
import numpy as np

import zealpy_gyrokin_
from zealpy import zealpy
from pyreldisp_params import params, pyreldisp_exception

#-------------------------------------------------------------
# Class 'zealpy_gyrokin'  
#-------------------------------------------------------------
class zealpy_gyrokin(zealpy):
    def __init__(self,Zi,NNr,rmin,Lr,B0):
        super(zealpy_gyrokin,self).__init__()
        self.setup_pymodule(zealpy_gyrokin_)

        required_params={'kpar'    : float,
                         'ktheta'  : float,
                         'Zi'      : int,
                         'NNr'     : int,
                         'dr'      : float,
                         'rmesh'   : np.ndarray,
                         'Ti'      : np.ndarray,
                         'dlogTi'  : np.ndarray, 
                         'Te'      : np.ndarray,
                         'n0'      : np.ndarray,
                         'dlogn0'  : np.ndarray,
                         'ddlogn0' : np.ndarray,
                         'btheta'  : np.ndarray,
                         'bz'      : np.ndarray,
                         'ordre_grandeur' : int,
                         'B0' : float
                       }

        self.params = params(required_params)
        self.params.set_value('ordre_grandeur',1)

        #--> Initialize the charge
        self.params.set_value('Zi',Zi)
        self.params.set_value('B0',B0)

        #--> Initialize the mesh
        [dr,rmesh] = init_mesh(NNr,rmin,Lr)
        self.params.set_value('NNr',NNr)
        self.params.set_value('dr',dr)
        self.params.set_value('rmesh',rmesh)
    #end def __init__


    #-------------------------------------------------------------
    # Initialize the dispersion relation
    #-------------------------------------------------------------
    def init_rel_disp(self,kpar,ktheta):
        #--> Fill the different parameters required in fun_gyrokin.f90
        self.params.set_value('kpar',kpar)
        self.params.set_value('ktheta',ktheta)
        
        #--> Fill fortran variable
        self.fill_fortran_var()
        
        #--> Call the initialization of the dispersion relation
        self.function_input_module.init()
    #end def init_rel_disp
#end class zealpy_gyrokin


#-------------------------------------------------------------
# Fill the variables for the mesh definition:
#   - dr, Lr, Lz and rmesh
#-------------------------------------------------------------
def init_mesh(NNr,rmin,Lr):
    dr    = np.float(Lr)/(NNr-1)
    rmesh = rmin + dr*np.arange(NNr)
    return [dr,rmesh]
#end def init_mesh
