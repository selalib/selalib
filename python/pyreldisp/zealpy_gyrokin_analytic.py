import sys
import numpy as np
from types            import FunctionType
from zealpy_gyrokin   import zealpy_gyrokin
from pyreldisp_params import params, pyreldisp_exception

#-------------------------------------------------------------------------------
# Class 'zealpy_gyrokin_cleanup_anal'  
#-------------------------------------------------------------------------------
class zealpy_gyrokin_analytic( zealpy_gyrokin ):
    def __init__( self, **kwargs ):
        required_params={'Zi'      : int,
                         'NNr'     : int,
                         'rmin'    : float,
                         'Lr'      : float,
                         'R0'      : float,
                         'invLTi'  : float,
                         'deltarTi': float,
                         'invLTe'  : float,
                         'deltarTe': float,
                         'invLn0'  : float,
                         'deltarn0': float,
                         'B0'      : float,
                         'iota'    : FunctionType }

        self.params_analytic = params( required_params, kwargs )
        Zi   = self.params_analytic.get_value('Zi')
        NNr  = self.params_analytic.get_value('NNr')
        rmin = self.params_analytic.get_value('rmin')
        Lr   = self.params_analytic.get_value('Lr')
        R0   = self.params_analytic.get_value('R0')
        B0   = self.params_analytic.get_value('B0')
        iota = self.params_analytic.get_value('iota')
        super( zealpy_gyrokin_analytic, self ).__init__( Zi, NNr, rmin, Lr, B0 )

        # Compute the analytic profiles
        dr    = self.params.get_value('dr')
        rmesh = self.params.get_value('rmesh')
        anorm = rmesh[-1]-rmesh[0]
        rpeak = rmesh[0]+0.5*anorm

        invLTi   = self.params_analytic.get_value('invLTi')
        deltarTi = self.params_analytic.get_value('deltarTi')
        invLTe   = self.params_analytic.get_value('invLTe')
        deltarTe = self.params_analytic.get_value('deltarTe')
        invLn0   = self.params_analytic.get_value('invLn0')
        deltarn0 = self.params_analytic.get_value('deltarn0')

        [Ti,dlogTi] = compute_Ti_dlogTi_analytic( rmesh, rpeak, invLTi, deltarTi )
        [Te,dlogTe] = compute_Te_dlogTe_analytic( rmesh, rpeak, invLTe, deltarTe )
        [n0,dlogn0,ddlogn0] = compute_n0_dlogn0_analytic(
                              rmesh, rpeak, invLn0, deltarn0)
        [btheta,bz] = compute_b_analytic( rmesh, iota( rmesh ), R0 )

        Zi = self.params_analytic.get_value('Zi')
        self.params.set_value('Zi',Zi)
        self.params.set_value('Ti',Ti)
        self.params.set_value('dlogTi',dlogTi)
        self.params.set_value('Te',Te)
        self.params.set_value('n0',n0)
        self.params.set_value('dlogn0',dlogn0)
        self.params.set_value('ddlogn0',ddlogn0)
        self.params.set_value('btheta',btheta)
        self.params.set_value('bz',bz)
        B0 = self.params_analytic.get_value('B0')
        self.params.set_value('B0',B0)
        if (not self.params_analytic.check_validity()):
           raise pyreldisp_exception("zealpy_gyrokin_analytic: Parameters not valid")

#-------------------------------------------------------------------------------
# Compute ion temperature profile Ti(r) and dlogTi(r)=d/dr(ln Ti(r))
#-------------------------------------------------------------------------------
def compute_Ti_dlogTi_analytic(rmesh,rpeak,invLTi,deltarTi):

    Ti     = np.exp(-invLTi*deltarTi*np.tanh((rmesh-rpeak)/deltarTi))
    dlogTi = -invLTi/np.cosh((rmesh-rpeak)/deltarTi)**2

    return [Ti,dlogTi]

#-------------------------------------------------------------------------------
# Compute electron temperature profile Te(r) and dlogTe(r)=d/dr(ln Te(r))
#-------------------------------------------------------------------------------
def compute_Te_dlogTe_analytic(rmesh,rpeak,invLTe,deltarTe):

    Te     = np.exp(-invLTe*deltarTe*np.tanh((rmesh-rpeak)/deltarTe))
    dlogTe = -invLTe/np.cosh((rmesh-rpeak)/deltarTe)**2

    return [Te,dlogTe]

#-------------------------------------------------------------------------------
# Compute electron density profile n0(r), dlogn0(r)=d/dr(ln n0(r)) and
# ddlogn0(r)=d/dr(dlogn0)
#-------------------------------------------------------------------------------
def compute_n0_dlogn0_analytic(rmesh,rpeak,invLn0,deltarn0):

    n0      = np.exp(-invLn0*deltarn0*np.tanh((rmesh-rpeak)/deltarn0))
    dlogn0  = -invLn0/np.cosh((rmesh-rpeak)/deltarn0)**2
    ddlogn0 = 2.0*invLn0/deltarn0*np.sinh((rmesh-rpeak)/deltarn0) \
              /np.cosh((rmesh-rpeak)/deltarn0)**3

    return [n0,dlogn0,ddlogn0]

#-------------------------------------------------------------------------------
# Compute btheta(r) and bz(r)
#-------------------------------------------------------------------------------
def compute_b_analytic( rmesh, iota, R0 ):

    xi     = iota*rmesh/R0
    btheta = xi /np.sqrt(1.0+xi**2)
    bz     = 1.0/np.sqrt(1.0+xi**2)

    return [btheta,bz]
