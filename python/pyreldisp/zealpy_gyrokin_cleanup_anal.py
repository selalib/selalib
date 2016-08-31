import sys
import numpy as np
from zealpy_gyrokin_cleanup import zealpy_gyrokin_cleanup
from pyreldisp_params import params, pyreldisp_exception

#-------------------------------------------------------------
# Class 'zealpy_gyrokin_cleanup_anal'  
#-------------------------------------------------------------
class zealpy_gyrokin_cleanup_anal(zealpy_gyrokin_cleanup):
    def __init__(self,**anal_arg):
        required_params={'Zi'      : int,
                         'NNr'     : int,
                         'rmin'    : float,
                         'Lr'      : float,
                         'invLTi'  : float,
                         'deltarTi': float,
                         'invLTe'  : float,
                         'deltarTe': float,
                         'invLn0'  : float,
                         'deltarn0': float,
                         'iota'    : float,
                         'B0'    : float
                        }
        self.param_anal = params(required_params,anal_arg)
        Zi   = self.param_anal.get_value('Zi')
        NNr  = self.param_anal.get_value('NNr')
        rmin = self.param_anal.get_value('rmin')
        Lr   = self.param_anal.get_value('Lr')
        iota   = self.param_anal.get_value('iota')
        B0   = self.param_anal.get_value('B0')
        super(zealpy_gyrokin_cleanup_anal,self).__init__(Zi,NNr,rmin,Lr,B0)    

        #--> Compute the analytic profiles
        dr    = self.params.get_value('dr')
        rmesh = self.params.get_value('rmesh')
        anorm = rmesh[-1]-rmesh[0]
        rpeak = rmesh[0]+0.5*anorm

        invLTi   = self.param_anal.get_value('invLTi')
        deltarTi = self.param_anal.get_value('deltarTi')
        invLTe   = self.param_anal.get_value('invLTe')
        deltarTe = self.param_anal.get_value('deltarTe')
        invLn0   = self.param_anal.get_value('invLn0')
        deltarn0 = self.param_anal.get_value('deltarn0')

        [Ti,dlogTi]         = compute_Ti_dlogTi_anal(
            rmesh,rpeak,invLTi,deltarTi)
        [Te,dlogTe]         = compute_Te_dlogTe_anal(
            rmesh,rpeak,invLTe,deltarTe)
        [n0,dlogn0,ddlogn0] = compute_n0_dlogn0_anal(
            rmesh,rpeak,invLn0,deltarn0)
        [btheta,bz] = compute_b_anal(
            rmesh,iota)
        self.eta = dlogTi/dlogn0

        Zi = self.param_anal.get_value('Zi')
        self.params.set_value('Zi',Zi)
        self.params.set_value('Ti',Ti)
        self.params.set_value('dlogTi',dlogTi)
        self.params.set_value('Te',Te)
        self.params.set_value('n0',n0)
        self.params.set_value('dlogn0',dlogn0)
        self.params.set_value('ddlogn0',ddlogn0)
        self.params.set_value('btheta',btheta)
        self.params.set_value('bz',bz)
        B0 = self.param_anal.get_value('B0')
        self.params.set_value('B0',B0)
        if (not self.param_anal.check_validity()):
            raise pyreldisp_exception("zealpy_gyrkin_cleanup_anal: Parameters not valid")
#end class zealpy_gyrokin_cleanup_anal    


#-------------------------------------------------
# Compute ion temperature profile Ti and 
#   dlogTi = d/dr(ln Ti(r))
# with
#   1/Ti(r) dTi(r)/dr = -(1/LTi)*cosh^-2(r-rpeak/deltar)   
#     where (1/Ln) = kappa_Ti/R   
#-------------------------------------------------
def compute_Ti_dlogTi_anal(rmesh,rpeak,invLTi,deltarTi):
    Ti     = np.exp(-invLTi*deltarTi * 
                     np.tanh((rmesh-rpeak)/deltarTi))
    dlogTi = -invLTi/np.cosh((rmesh-rpeak)/deltarTi)**2
    return [Ti,dlogTi]
#end def compute_Ti_dlogTi_anal


#-------------------------------------------------
# Compute electron temperature profile Te and 
#   dlogTe = d/dr(ln Te(r))
# with
#   1/Te(r) dTe(r)/dr = -(1/LTe)*cosh^-2(r-rpeak/deltar)   
#     where (1/Ln) = kappa_Te/R   
#-------------------------------------------------
def compute_Te_dlogTe_anal(rmesh,rpeak,invLTe,deltarTe):
    Te     = np.exp(-invLTe*deltarTe * 
                     np.tanh((rmesh-rpeak)/deltarTe))
    dlogTe = -invLTe/np.cosh((rmesh-rpeak)/deltarTe)**2
    return [Te,dlogTe]
#end def compute_Te_dlogTe_anal


#-------------------------------------------------
# Compute electron density profile n0,
#   dlogn0  = d/dr(ln n0(r)) and 
#   ddlogn0 = d/dr (dlogn0)
# with
#   1/n0(r) dn0(r)/dr = -(1/Ln0)*cosh^-2(r-rpeak/deltar)   
#     where (1/Ln) = kappa_n0/R   
#-------------------------------------------------
def compute_n0_dlogn0_anal(rmesh,rpeak,invLn0,deltarn0):
    n0     = np.exp(-invLn0*deltarn0 * 
                     np.tanh((rmesh-rpeak)/deltarn0))
    dlogn0 = -invLn0/np.cosh((rmesh-rpeak)/deltarn0)**2
    ddlogn0 = 2*invLn0/deltarn0*np.sinh((rmesh-rpeak)/deltarn0) / \
        np.cosh((rmesh-rpeak)/deltarn0)**3
    
    return [n0,dlogn0,ddlogn0]
#end def compute_n0_dlogn0_anal

#-------------------------------------------------
# Compute c,btheta,bz
#   c  = iota*r  normally iota*r/R0 
#-------------------------------------------------
def compute_b_anal(rmesh,iota):
  #c = iota*rmesh
  #iota = 0.8/239.8
  btheta = iota*rmesh/(1.0+(iota*rmesh)**2)**(0.5)
  bz = 1.0/(1.0+(iota*rmesh)**2)**(0.5)   
  #bz = 1.
  
  return [btheta,bz]
#end def compute_b_anal
