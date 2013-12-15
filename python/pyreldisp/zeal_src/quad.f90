!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Quad_Module
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
!  This module contains certain global variables that are needed 
!  for the calculation of the total number of zeros. 
!
!**END PROLOGUE Quad_Module
!-----------------------------------------------------------------------

MODULE Quad_Module

     USE Precision_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC

!-----------------------------------------------------------------------
!  The total number of zeros is calculated via numerical integration.
!  The value of INTACC is used in tests concerning the estimated 
!  quadrature error. The error should be < 0.5 as we are calculating 
!  an integer.
!
     REAL(KIND=DP), PARAMETER :: INTACC = 0.45_DP

!-----------------------------------------------------------------------
!  The numerical integration is along the edges of the rectangular box
!  whose lower left vertex is given by (X0,Y0) and whose edges have 
!  length DX and DY, respectively.
!
     REAL(KIND=DP) :: X0, Y0, DX, DY


END MODULE Quad_Module
