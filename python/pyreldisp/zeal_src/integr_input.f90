!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Integration_Input_Module
!**DATE WRITTEN   1999/06/01  (year/month/day)
!**AUTHORS
!
!  Peter Kravanja, Marc Van Barel
!
!    Department of Computer Science
!    Katholieke Universiteit Leuven
!    Celestijnenlaan 200 A
!    B-3001 Heverlee
!    BELGIUM
!    E-mail: Peter.Kravanja@na-net.ornl.gov
!            Marc.VanBarel@cs.kuleuven.ac.be
!
!  Omiros Ragos, Michael N. Vrahatis, and Filareti A. Zafiropoulos
!
!    Department of Mathematics
!    University of Patras
!    GR-261.10 Patras
!    GREECE
!    E-mail: {ragos,vrahatis,phikapa}@math.upatras.gr
!
!**DESCRIPTION
!
!  The user can specify the values of the various input parameters 
!  needed for numerical integration. 
!
!**END PROLOGUE Integration_Input_Module 
!-----------------------------------------------------------------------

MODULE Integration_Input_Module

     USE Precision_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC

!-----------------------------------------------------------------------
!  Specify the values of NUMABS and NUMREL.
!  These variables determine the absolute and relative accuracy to which
!  the integrals that calculate the number of zeros are to be evaluated.
!  If NUMREL = 0.0_DP, then only an absolute criterion will be used.
!  If NUMABS = 0.0_DP, then only a relative criterion will be used.
!
!  If NUMABS and NUMREL are both too small, then the numerical
!  integration may be time-consuming. If they are both too large, then
!  the calculated number of zeros may be wrong. The default values of
!  NUMABS and NUMREL are 0.07_DP and 0.0_DP, respectively.
!
     REAL(KIND=DP) :: NUMABS = 0.07_DP
     REAL(KIND=DP) :: NUMREL = 0.0_DP

!-----------------------------------------------------------------------
!  Specify the values of INTABS and INTREL.
!  These variables determine the absolute and relative accuracy to which
!  the integrals that are used to compute approximations for the zeros
!  are to be calculated.
!  If INTREL = 0.0_DP, then only an absolute criterion will be used.
!  If INTABS = 0.0_DP, then only a relative criterion will be used.
!
!  If INTABS and INTREL are both too small, then the numerical
!  integration may be time-consuming. If they are both too large, then
!  the approximations for the zeros may be very inaccurate and Newton's
!  method, which is used to refine these approximations (see NEWTONZ
!  and NEWTONF), may fail. The default values of INTABS and INTREL are
!  0.0_DP and 1.0E-12_DP, respectively.
!
     REAL(KIND=DP) :: INTABS = 0.0_DP
     REAL(KIND=DP) :: INTREL = 1.0E-12_DP

!-----------------------------------------------------------------------
!  Specify the value of EPS_STOP.
!  This variable is used in the stopping criterion that determines the
!  number of mutually distinct zeros. If EPS_STOP is too large, then
!  the computed number of mutually distinct zeros may be too small.
!  If EPS_STOP is too small, then the computed number of mutually
!  distinct zeros may be too large, especially in case the function has
!  many multiple zeros. A recommended value is 1.0E-08_DP.
!
     REAL(KIND=DP) :: EPS_STOP = 1.0E-08_DP

END MODULE Integration_Input_Module

