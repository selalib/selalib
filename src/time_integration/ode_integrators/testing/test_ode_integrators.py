from __future__ import print_function

import os, sys
import subprocess

#==============================================================================
executable = 'test_ode_integrators'
directory  = os.path.dirname( os.path.realpath( __file__ ) )

vals_ode_selector    = [0]
vals_odeint_selector = [0, 1, 2, 3, 4, 5, 10, 12]
vals_h_factor        = [1.0, 2.0, 4.0, 8.0]

cmd_template = "{:s}/{:s} {:d} {:d} {:f}"

#==============================================================================
for ode_selector in vals_ode_selector:
    for odeint_selector in vals_odeint_selector:
        for h_factor in vals_h_factor:
            cmd = cmd_template.format( directory,
                    executable,
                    ode_selector,
                    odeint_selector,
                    h_factor
                    )
            line = '-' * len( cmd )
            print( line ) 
            print( cmd  )
            print( line )
            sys.stdout.flush()
            subprocess.check_call( cmd, shell=True )

#==============================================================================
# TODO: compute convergence rates
print( "PASSED" )

