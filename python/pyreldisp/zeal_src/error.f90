!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Error_Module
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
!  This module contains a global variable that is used to take care
!  of the error handling.
!
!**END PROLOGUE Error_Module
!-----------------------------------------------------------------------

MODULE Error_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC

!-----------------------------------------------------------------------
!  INFO = 0  =>  Improper input parameters
!         1  =>  No errors  
!         2  =>  Calculation of the total number of zeros has failed
!         3  =>  Isolation of the zeros has failed
!         4  =>  Computation of the zeros has failed
!
     INTEGER, SAVE :: INFO = 1


END MODULE Error_Module
