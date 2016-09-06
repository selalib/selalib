import sys
import numpy as np

import zealpy_gyrokin_cleanup_
from zealpy import zealpy
from pyreldisp_params import params, pyreldisp_exception

#-------------------------------------------------------------------------------
# Class 'zealpy_gyrokin_cleanup'
#-------------------------------------------------------------------------------
class zealpy_gyrokin_cleanup(zealpy):
    def __init__(self,Zi,NNr,rmin,Lr,B0):
        super(zealpy_gyrokin_cleanup,self).__init__()
        self.setup_pymodule(zealpy_gyrokin_cleanup_)

        required_params={'kzeta'   : float,
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
                         'B0' : float }

        self.params = params(required_params)
        self.params.set_value('ordre_grandeur',1)

        # Initialize the charge
        self.params.set_value('Zi',Zi)
        self.params.set_value('B0',B0)

        # Initialize the mesh
        [dr,rmesh] = init_mesh(NNr,rmin,Lr)
        self.params.set_value('NNr',NNr)
        self.params.set_value('dr',dr)
        self.params.set_value('rmesh',rmesh)
   

    # Initialize the dispersion relation
    def init_rel_disp(self,kzeta,ktheta):
        # Fill the different parameters required in fun_gyrokin_cleanup.f90
        self.params.set_value('kzeta',kzeta)
        self.params.set_value('ktheta',ktheta)
        
        # Fill fortran variable
        self.fill_fortran_var()
        
        # Call the initialization of the dispersion relation
        self.function_input_module.init()
   

#-------------------------------------------------------------------------------
# Fill the variables for the mesh definition: dr and rmesh
#-------------------------------------------------------------------------------
def init_mesh(NNr,rmin,Lr):
    dr    = np.float(Lr)/(NNr-1)
    rmesh = rmin+dr*np.arange(NNr)
    return [dr,rmesh]

