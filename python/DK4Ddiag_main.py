#-------------------------------------------------
#  file : DK4Ddiag_main.py
#  date : 25/01/2014
#
#   Diagnostics for DK4D code
#-------------------------------------------------
import GYSutils
GYSut      = GYSutils.impall('GYSutils')
DK4D_Phi2D = GYSutils.impall('DK4Ddiag_Phi2D')
DK4D_f2D   = GYSutils.impall('DK4Ddiag_f2D')
DK4D_CL   = GYSutils.impall('DK4Ddiag_CL')

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
        print ' '

        print ' '
        print ' **** PHI3D DIAGNOSTICS ****'
        print '  *   2 : int Phi(x,y,z)**2 jac dx dy dz'
        print ' '

        print ' **** F2D DIAGNOSTICS ****'
        print '  *   21 : f2D(x,y)'
        print '  *   22 : f2D(z,vpar)'
        print ' '

        print ' **** CONSERVATION LAWS ****'
        print '  *   31 : L1-Norm, L2-Norm and Linf-Norm'
        print '  *   32 : Energy'
        print ' '

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
        DK4D_Phi2D.Phi2D_xy_plot(resu)
    elif (ichoice == 21):
        DK4D_f2D.f2D_xy_plot(resu)
    elif (ichoice == 22):
        DK4D_f2D.f2D_zvpar_plot(resu)
    elif (ichoice == 31):
        DK4D_CL.L1_L2_Linf_norms_plot(resu)
    elif (ichoice == 32):
        DK4D_CL.energies_plot(resu)
    elif (ichoice == 2):
        DK4D_Phi2D.Int_phisquare_plot(resu)

#end def DK4Ddiag_choice

