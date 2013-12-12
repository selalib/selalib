import sys
import numpy as np

import zealpy_BOT_
from zealpy import zealpy
from pyreldisp_params import params, pyreldisp_exception

#-------------------------------------------------------------
# Class 'zealpy_BOT'  
#-------------------------------------------------------------
class zealpy_BOT(zealpy):
    def __init__(self,kmode,a_bump,v0,sigmab,kbT,l0ld):
        super(zealpy_BOT,self).__init__()
        self.setup_pymodule(zealpy_BOT_)

        required_params={
            'kmode'  : float,
            'a_bump' : float,
            'v0'     : float,
            'kbT'    : float,
            'sigmab' : float,
            'l0ld'   : float,
            }

        self.params = params(required_params)
        self.params.set_value('kmode',kmode)
        self.params.set_value('a_bump',a_bump)
        self.params.set_value('v0',v0)
        self.params.set_value('kbT',kbT)
        self.params.set_value('sigmab',sigmab)
        self.params.set_value('l0ld',l0ld)

        #--> Fill fortran variable
        self.fill_fortran_var()        
    #end def __init__
#end class zealpy_BOT
