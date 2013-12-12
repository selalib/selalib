from zealpy import zealpy
import zealpy_landau_ 

#-------------------------------------------------------------
# Class 'zealpy_landau'  
#-------------------------------------------------------------
class zealpy_landau(zealpy):
    def __init__(self,kmode=0.5):
        super(zealpy_landau,self).__init__()
        self.setup_pymodule(zealpy_landau_)

        self.function_input_module.kmode = kmode
    #end def __init__
#end class zealpy_landau

