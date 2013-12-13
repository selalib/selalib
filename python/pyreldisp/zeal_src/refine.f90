!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Refine_Module 
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
!  This module contains subroutines that are used to refine the 
!  computed approximations for the zeros. We have implemented the 
!  modified Newton's method, which takes into account the multiplicity 
!  of a zero and converges quadratically,
!
!**END PROLOGUE Refine_Module
!-----------------------------------------------------------------------

MODULE Refine_Module 

     USE Precision_Module
     USE Function_Input_Module
     USE Zeal_Input_Module
      
     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PRIVATE
     PUBLIC :: REFINE, REFINE_EMPTY 


CONTAINS


     SUBROUTINE REFINE(POINT,STEP,ZEROS,FZEROS,MULTIPLICITIES, &
                       NUMBERDISTINCT,REFINE_OK)
!-----------------------------------------------------------------------
!**PURPOSE
!  Each computed approximation for a zero is refined via the modified 
!  Newton's method. Depending on whether this procedure is successful or
!  not, the corresponding element of REFINE_OK is set .TRUE. or .FALSE.
!-----------------------------------------------------------------------
!  Parameters
!
     INTEGER, INTENT(IN)                     :: NUMBERDISTINCT
     INTEGER, DIMENSION(M), INTENT(IN)       :: MULTIPLICITIES
     LOGICAL, DIMENSION(M), INTENT(OUT)      :: REFINE_OK
     REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: POINT, STEP
     COMPLEX(KIND=DP), DIMENSION(M), INTENT(INOUT) :: ZEROS 
     COMPLEX(KIND=DP), DIMENSION(M), INTENT(INOUT) :: FZEROS

!-----------------------------------------------------------------------
!  Local variables
!
     INTEGER :: H, J, L, MULT
     INTEGER, PARAMETER :: NTRIAL = 8
     LOGICAL :: NEWTOK
     REAL(KIND=DP) :: BOXDIAM, PI 
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: DISTANCE 
     COMPLEX(KIND=DP) :: Z, ZZ, F

!-----------------------------------------------------------------------
     INTRINSIC SQRT, ABS, REAL, AIMAG, MIN, EXP

!
! For each approximation ZEROS(J) we find the nearest other
! approximation in the vector ZEROS. This information is used
! in case Newton's method, started at ZEROS(J), fails.
!
     BOXDIAM = SQRT(STEP(1)**2 + STEP(2)**2)
     ALLOCATE(DISTANCE(NUMBERDISTINCT))
     
     DO J = 1, NUMBERDISTINCT
       Z = ZEROS(J)
       DISTANCE(J) = BOXDIAM
       DO L = 1, NUMBERDISTINCT
         IF ( L .NE. J ) THEN
           ZZ = ZEROS(L)
           IF ( ABS(Z-ZZ) .LT. DISTANCE(J) ) DISTANCE(J) = ABS(Z-ZZ) 
         END IF 
       END DO
     END DO
!
!  Now we take into account the distance to the edges of the box.
!
     DO J = 1, NUMBERDISTINCT
       Z = ZEROS(J)
       DISTANCE(J) = MIN(DISTANCE(J),ABS(REAL(Z,DP)-POINT(1)), &
                         ABS(POINT(1)+STEP(1)-REAL(Z,DP)),     &
                         ABS(AIMAG(Z)-POINT(2)),               &
                         ABS(POINT(2)+STEP(2)-AIMAG(Z)))
       DISTANCE(J) = 0.7_DP*DISTANCE(J)
     END DO

     PI = 4.0_DP*ATAN(ONE)
!
!  Newton's method
!
     DO J = 1, NUMBERDISTINCT
       REFINE_OK(J) = .FALSE.
       Z = ZEROS(J)
       MULT = MULTIPLICITIES(J)
       CALL NEWTON(POINT,STEP,MULT,Z,F,NEWTOK)
!
!  Did Newton converge to a zero that we have already found?
!
       DO L = 1, J-1
         IF ( ABS(Z-ZEROS(L)) .LT. NEWTONZ ) NEWTOK = .FALSE.
       END DO

       IF ( NEWTOK ) THEN
         REFINE_OK(J) = .TRUE.
       ELSE
         INTDO: DO L = 1, NTRIAL
           Z = ZEROS(J) + DISTANCE(J)*EXP(L*TWO*PI*I/NTRIAL)
           CALL NEWTON(POINT,STEP,MULT,Z,F,NEWTOK)
           DO H = 1, J-1
             IF ( ABS(Z-ZEROS(H)) .LT. NEWTONZ ) NEWTOK = .FALSE.
           END DO
           IF ( NEWTOK ) THEN
             REFINE_OK(J) = .TRUE.
             EXIT INTDO
           END IF
         END DO INTDO
       END IF

       IF ( REFINE_OK(J) ) THEN
         ZEROS(J) = Z
         FZEROS(J) = F  
       END IF
     END DO

     DEALLOCATE(DISTANCE)

     END SUBROUTINE REFINE
!-----------------------------------------------------------------------


     SUBROUTINE NEWTON(POINT,STEP,MULT,Z,F,NEWTOK)
!-----------------------------------------------------------------------
!**PURPOSE
!  Implementation of the modified Newton's method.
!-----------------------------------------------------------------------
!  Parameters
!
     LOGICAL, INTENT(OUT) :: NEWTOK
     INTEGER, INTENT(IN)  :: MULT
     REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: POINT, STEP
     COMPLEX(KIND=DP), INTENT(INOUT)  :: Z
     COMPLEX(KIND=DP), INTENT(OUT)    :: F

!-----------------------------------------------------------------------
!  Local variables 
!
     INTEGER :: K
     INTEGER, PARAMETER :: MAXIT = 20
     COMPLEX(KIND=DP) :: DF, ZNEW

!-----------------------------------------------------------------------
     INTRINSIC ABS

     NEWTOK = .FALSE.

     CALL FDF(Z,F,DF)
     IF ( ABS(F) .LT. NEWTONF ) THEN
       NEWTOK = .TRUE.
       RETURN
     END IF

     DO K = 1, MAXIT
       ZNEW = Z - MULT*F/DF

       IF ( .NOT. INSIDEBOX(ZNEW,POINT,STEP,0.1_DP) ) RETURN

       CALL FDF(ZNEW,F,DF)

       IF ( ABS(F) .LT. NEWTONF .OR.            &
            ABS(ZNEW - Z)/ABS(Z) .LT. NEWTONZ ) THEN
         Z = ZNEW
         NEWTOK = .TRUE.
         RETURN
       END IF

       Z = ZNEW
     END DO

     END SUBROUTINE NEWTON
!-----------------------------------------------------------------------


     FUNCTION INSIDEBOX(Z,POINT,STEP,ETA)
!-----------------------------------------------------------------------
!**PURPOSE
!  Check if the point Z lies inside the box specified by POINT and STEP.
!  Allow a tolerance ETA, i.e., consider the box that is 2*ETA larger
!  in each direction.  
!-----------------------------------------------------------------------
     LOGICAL INSIDEBOX
     REAL(KIND=DP), INTENT(IN) :: ETA
     REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: POINT, STEP
     COMPLEX(KIND=DP), INTENT(IN)  :: Z
 
     INTRINSIC REAL, AIMAG
 
     INSIDEBOX = REAL(Z,DP) .GE. POINT(1)-ETA*STEP(1)           .AND.   &
                 REAL(Z,DP) .LE. POINT(1)+STEP(1)+ETA*STEP(1)   .AND.   &
                 AIMAG(Z)   .GE. POINT(2)-ETA*STEP(2)           .AND.   &
                 AIMAG(Z)   .LE. POINT(2)+STEP(2)+ETA*STEP(2)

     END FUNCTION INSIDEBOX
!-----------------------------------------------------------------------


     SUBROUTINE REFINE_EMPTY(REFINE_OK)
!-----------------------------------------------------------------------
!**PURPOSE
!  FOR TESTING PURPOSES ONLY!
!-----------------------------------------------------------------------
!  Parameters
!
     LOGICAL, DIMENSION(M), INTENT(OUT)      :: REFINE_OK

!-----------------------------------------------------------------------
     REFINE_OK = .TRUE.

     END SUBROUTINE REFINE_EMPTY


END MODULE Refine_Module 
