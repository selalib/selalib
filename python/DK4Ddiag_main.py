#-------------------------------------------------
#  file : DK4Ddiag_main.py
#  date : 25/01/2014
#
#   Diagnostics for DK4D code
#-------------------------------------------------
import GYSutils
GYSut      = GYSutils.impall('GYSutils')
DK4D_Phi2D = GYSutils.impall('DK4Ddiag_Phi2D')

#-----------------------------------------------
# Diagnostic choices
#-----------------------------------------------
def DK4Ddiag(resu):
    """ Diagnostic choices"""

    choice = '1';
    while choice != '' :
        print ' '
        print ' **** PHI2D DIAGNOSTICS ****'
        print '  *   1 : Phi(x,y)'

        print ' **** F2D DIAGNOSTICS ****'
        print '  *   21 : f2D(x,y)'

        choice = raw_input(' choice ? : ');
        choice = GYSut.str2num(choice)

        if (choice != ''):
            DK4Ddiag_choice(choice,resu)

    #end while

    print 'Bye bye!'
#end def DK4Ddiag


#-----------------------------------------------------------------
#  Execution of DK4D diagnostics according to the user choice
#-----------------------------------------------------------------
def DK4Ddiag_choice(ichoice,resu):
    """
    Execution of GYSELA diagnostics according to the user choice
    """

    # Phi2D diagnostics
    if (ichoice == 1):
        DK4D_Phi2D.GPhi2D_xy(resu)
    #elif (ichoice == 21):
    #    DK4D_f2D_xy(resu)

#end def DK4Ddiag_choice

