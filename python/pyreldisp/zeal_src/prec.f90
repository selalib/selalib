!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Precision_Module 
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
!**END PROLOGUE Precision_Module 
!-----------------------------------------------------------------------

MODULE Precision_Module 

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC

!-----------------------------------------------------------------------
!  To which precision are the floating point calculations to be done?
!  The following corresponds to Fortran 77's "double precision".
!
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,70)

!-----------------------------------------------------------------------
!  The following auxiliary constants are used in several parts of the code.
!  They may be used in the definition of the function and its derivative
!  (subroutine FDF).
!
     REAL(KIND=DP), PARAMETER    :: ZERO = 0.0_DP
     REAL(KIND=DP), PARAMETER    :: ONE  = 1.0_DP
     REAL(KIND=DP), PARAMETER    :: TWO  = 2.0_DP
     COMPLEX(KIND=DP), PARAMETER :: I = (0.0_DP,1.0_DP)

END MODULE Precision_Module
