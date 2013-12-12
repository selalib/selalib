!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Zeal_Module 
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
!  This modules contains the subroutine ZEAL, which is the main subroutine
!  of the package, and also the subroutine CHECK_INPUT, which is used to
!  check the input parameters specified by the user.
!----------------------------------------------
! 
! modified by Eric Sonnendrucker (University of Strasbourg)
!             sonnen@math.u-strasbg.fr
!  
! for binding with python
!
!**END PROLOGUE Zeal_Module 
!-----------------------------------------------------------------------

MODULE Zout

     USE Precision_Module
     USE Error_Module
     USE Zeal_Input_Module
     USE Split_Module
     USE Zeros_Module
     USE Refine_Module

     IMPLICIT NONE

     INTEGER :: TOTALNUMBER, DISTINCTNUMBER, REFINEDNUMBER, INFORMATION 
     INTEGER, DIMENSION(:), ALLOCATABLE          :: MULTIPLICITIES
     LOGICAL, DIMENSION(:), ALLOCATABLE          :: REFINEMENT_OK
     COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZEROS, FZEROS
     LOGICAL :: FAIL


!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
!     PRIVATE
!     PUBLIC :: ZEAL 


CONTAINS


!     SUBROUTINE ZEAL(TOTALNUMBER,DISTINCTNUMBER,ZEROS,FZEROS,       &
!                     MULTIPLICITIES,REFINEDNUMBER,REFINEMENT_OK)
  SUBROUTINE ZEAL
!-----------------------------------------------------------------------
!**PURPOSE
!
!  Given a rectangular region W in the complex plane and an analytic
!  function f, such that the boundary of W does not contain zeros of f,
!  the package ZEAL calculates all the zeros of f that lie inside W,
!  together with their respective multiplicities.
!
!  ZEAL uses an integral formula to calculate the total number of zeros
!  (counting multiplicities) of f that lie inside W. Then, by using the
!  same procedure, W is subdivided into subregions that contain at most
!  M zeros (counting multiplicities). Approximations for these zeros are
!  calculated via an algorithm that is based on numerical integration
!  along the boundaries of the subregions and, also, on generalized
!  eigenvalue problems . The multiplicities of the zeros are calculated
!  by solving a Vandermonde system. The approximations for the zeros are
!  refined via the modified Newton's method, which takes into account the
!  multiplicity of a zero and converges quadratically.
!
!-----------------------------------------------------------------------
!  Parameters
!
!     INTEGER, INTENT(OUT)                     :: TOTALNUMBER
!     INTEGER, INTENT(OUT)                     :: DISTINCTNUMBER
!     INTEGER, DIMENSION(:), POINTER           :: MULTIPLICITIES 
!     COMPLEX(KIND=DP), DIMENSION(:), POINTER  :: ZEROS, FZEROS
!
!  The approximations for the zeros that have been calculated by the 
!  subroutine APPROXIMATE are refined iteratively by the subroutine REFINE.
!  REFINEDNUMBER is equal to the number of approximations for the zeros 
!  that the procedure REFINE has been able to refine successfully.
!  The boolean array REFINEMENT_OK indicates whether this has been the
!  case or not for a particular zero. 
!
!     INTEGER, INTENT(OUT)                     :: REFINEDNUMBER
!     LOGICAL, DIMENSION(:), POINTER           :: REFINEMENT_OK
!-----------------------------------------------------------------------
!  Local variables
!
     LOGICAL :: FLAG, ICON_4_FLAG
     INTEGER :: NRPERT, NR_ISOLATED_BOXES, J, K, L
     INTEGER :: REFINED_OK_AND_PRINTED
     INTEGER :: PARTIALSUM_NUMBERDISTINCT
     INTEGER :: NUMBER, NUMBER1, NUMBER2
     REAL(KIND=DP) :: EPSMCH, DLAMCH, TOL
     REAL(KIND=DP), DIMENSION(2) :: POINT, STEP, LVPERT, HPERT
     REAL(KIND=DP), DIMENSION(2) :: POINT1, STEP1, POINT2, STEP2
     REAL(KIND=DP), DIMENSION(6) :: INT, ERR, INT1, ERR1, INT2, ERR2

     TYPE :: ISOLATED_BOX
       INTEGER                     :: TOTAL_NUMBER_OF_ZEROS
       INTEGER                     :: DISTINCT_NUMBER_OF_ZEROS
       REAL(KIND=DP), DIMENSION(2) :: POINT, STEP
       TYPE(ISOLATED_BOX), POINTER :: NEXT
     END TYPE 
!
!  The pointer START_LIST_OF_ISOLATED_BOXES points to a list of boxes
!  that contain no more than M zeros. 
!  The number of boxes in this list is given by NR_ISOLATED_BOXES. 
!
     TYPE(ISOLATED_BOX), POINTER :: START_LIST_OF_ISOLATED_BOXES
     TYPE(ISOLATED_BOX), POINTER :: END_OF_LIST_OF_ISOLATED_BOXES
     TYPE(ISOLATED_BOX), POINTER :: PTR_ISOL_BOXES 
!
!  We need a stack to keep track of the boxes that still need to
!  be subdivided.
!
     TYPE :: STACK_BOX
       INTEGER                     :: NUMBER
       REAL(KIND=DP), DIMENSION(2) :: POINT, STEP
       REAL(KIND=DP), DIMENSION(6) :: INT, ERR
       TYPE(STACK_BOX), POINTER    :: PREVIOUS
     END TYPE

     TYPE(STACK_BOX), POINTER  :: STACK, PTR 
!
!  The subroutine APPROXIMATE computes approximations for the zeros 
!  that lie in one of the boxes in the list START_LIST_OF_ISOLATED_BOXES. 
!  The following variables are used:
!
     INTEGER                         :: NDIST_BOX
     LOGICAL, DIMENSION(M)           :: REFINE_OK
     INTEGER, DIMENSION(M)           :: MULT_BOX
     COMPLEX(KIND=DP), DIMENSION(M)  :: ZEROS_BOX, FZEROS_BOX
!
!  On output, the arrays MULTIPLICITIES, ZEROS, FZEROS and REFINEMENT_OK 
!  will be of size (DISTINCTNUMBER). They can only be allocated at the
!  very end, when the value of DISTINCTNUMBER is known. Thus, we define 
!  the following counterparts, which we will allocate size (TOTALNUMBER).
!
     INTEGER, DIMENSION(:), ALLOCATABLE          :: MULTIPLICITIES_TN
     COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: ZEROS_TN, FZEROS_TN
     LOGICAL, DIMENSION(:), ALLOCATABLE          :: REFINEMENT_OK_TN

!-----------------------------------------------------------------------
     INTRINSIC SQRT, REAL, AIMAG
     EXTERNAL DLAMCH

  
     IF ( VERBOSE ) THEN
       IF ( ICON .EQ. 4 ) THEN 
         PRINT 7000, LV, H, M, NR, ICON, FILES
       ELSE
         PRINT 7010, LV, H, M, ICON, FILES
       END IF
     END IF
!
!  Initialize FAIL
!
     FAIL = .FALSE.
!
!  Initialize INFO 
!
     INFO = 1
!
!  Improper input parameters ?
!
     CALL CHECK_INPUT

     IF ( INFO .EQ. 0 ) THEN
       CALL ERROR_EXIT
       RETURN
     END IF
!
!  Calculate the machine precision, set the value of TOL and 
!  perturb the rectangular box specified by the user.
!
     EPSMCH = DLAMCH('P')
     TOL = 10.0_DP*SQRT(EPSMCH)

     POINT(1) = LV(1) - 1.1_DP*TOL
     POINT(2) = LV(2) - 1.3_DP*TOL
     STEP(1)  = H(1) + 2.4_DP*TOL
     STEP(2)  = H(2) + 2.8_DP*TOL

!
!  Call INBOX to compute the total number of zeros.
!
     CALL INBOX(POINT,STEP,INT,ERR,NRPERT,FLAG)
     IF ( INFO .EQ. 2 ) THEN
       IF ( IFAIL .NE. 1 ) PRINT 7020
       CALL ERROR_EXIT
       RETURN
     END IF
!
!  If NRPERT < 0, then the calculation of the total number of zeros
!  has failed due to some problem that occured during the numerical 
!  evaluation of the integrals.
!
     IF ( NRPERT .LT. 0 ) INFO = 2
     IF ( INFO .EQ. 2 ) THEN
       IF ( IFAIL .NE. 1 ) PRINT 7020
       CALL ERROR_EXIT
       RETURN
     END IF
!
     TOTALNUMBER = NRPERT
!
!  Allocate ZEROS_TN, FZEROS_TN, MULTIPLICITIES_TN and REFINEMENT_OK_TN.
!
     ALLOCATE(ZEROS_TN(TOTALNUMBER))
     ALLOCATE(FZEROS_TN(TOTALNUMBER))
     ALLOCATE(MULTIPLICITIES_TN(TOTALNUMBER))
     ALLOCATE(REFINEMENT_OK_TN(TOTALNUMBER))
!
!  Get the new geometry of the box.
!
     LVPERT = POINT
     HPERT = STEP 
!
!  Note: NRPERT is equal to the number of zeros in the previous (perturbed)
!  box. As INBOX may have extended the user's box, NRPERT can be strictly 
!  larger than the number of zeros in the box specified by the user.
!
!  Print a warning message if the geometry of the initial box has been
!  changed by INBOX.
!
     IF ( FLAG .AND. VERBOSE ) PRINT 7030
!
!  The value of NR is changed in case it is larger than NRPERT, the
!  calculated total number of zeros.
!
     IF ( ICON .EQ. 4 .AND. NR .GT. NRPERT ) THEN
       IF ( VERBOSE ) PRINT 7040
       NR = NRPERT
     END IF
!
!  Print the box that we are considering and the value of NRPERT.
!
     IF ( VERBOSE ) THEN
       PRINT 7050, LVPERT(1), LVPERT(2), HPERT(1), HPERT(2)
       PRINT 7060, NRPERT
     END IF    
!
!  If there are no zeros or ICON = 1, then we may stop.
!
     IF ( NRPERT .EQ. 0 .OR. ICON .EQ. 1 ) RETURN
!
!  If the box specified by LVPERT and HPERT contains no more than M 
!  zeros, then we put it in a list of such boxes and, if required,
!  we calculate all the zeros that lie inside this box, together
!  with their corresponding multiplicities. 
!
     ALLOCATE(START_LIST_OF_ISOLATED_BOXES)
     ALLOCATE(END_OF_LIST_OF_ISOLATED_BOXES)
     NR_ISOLATED_BOXES = 0 

     DISTINCTNUMBER = 0
     REFINEDNUMBER = 0

     IF ( NRPERT .LE. M ) THEN                     
       START_LIST_OF_ISOLATED_BOXES%TOTAL_NUMBER_OF_ZEROS = NRPERT
       START_LIST_OF_ISOLATED_BOXES%POINT = LVPERT
       START_LIST_OF_ISOLATED_BOXES%STEP  = HPERT
       NULLIFY(START_LIST_OF_ISOLATED_BOXES%NEXT) 
       END_OF_LIST_OF_ISOLATED_BOXES => START_LIST_OF_ISOLATED_BOXES

       NR_ISOLATED_BOXES = 1
!
!  The approximations for the zeros are computed (in case the user
!  asked for it) via the subroutine APPROXIMATE.
!
       IF ( ICON .GE. 3 ) THEN
         CALL APPROXIMATE(LVPERT,HPERT,NRPERT,ZEROS_BOX,FZEROS_BOX, &
                          MULT_BOX,NDIST_BOX)
!
!  Did the procedure for the computation of the zeros fail?
!
         IF ( INFO .EQ. 4 ) THEN
           IF ( IFAIL .NE. 1 ) PRINT 7110
           CALL ERROR_EXIT
           RETURN
         END IF
!
!  We refine the computed approximations for the zeros. The boolean 
!  array REFINE_OK indicates (for 1,...,NDIST_BOX) if the procedure 
!  has been successful.
!
         CALL REFINE(LVPERT,HPERT,ZEROS_BOX,FZEROS_BOX, &
                     MULT_BOX,NDIST_BOX,REFINE_OK)
!--
!  If you don't want the computed approximations for the zeros to be
!  refined, for example in case you want to check the accuracy obtained
!  by APPROXIMATE, comment the previous call and use the following one.
!
!        CALL REFINE_EMPTY(REFINE_OK)
!--
 
!
!  We update our data structures.
!
         DO J = 1, NDIST_BOX
           ZEROS_TN(DISTINCTNUMBER+J) = ZEROS_BOX(J)
           FZEROS_TN(DISTINCTNUMBER+J) = FZEROS_BOX(J)
           MULTIPLICITIES_TN(DISTINCTNUMBER+J) = MULT_BOX(J)
           REFINEMENT_OK_TN(DISTINCTNUMBER+J) = REFINE_OK(J) 
           IF ( REFINE_OK(J) ) REFINEDNUMBER = REFINEDNUMBER + 1
         END DO
         DISTINCTNUMBER = DISTINCTNUMBER + NDIST_BOX
         START_LIST_OF_ISOLATED_BOXES%DISTINCT_NUMBER_OF_ZEROS = NDIST_BOX
       END IF
     ELSE

!  We split the box specified by LVPERT and HPERT.
!  The variables POINT, STEP, INT and ERR still contain the correct
!  corresponding values. We have to set the value of NUMBER, though.

       NUMBER = TOTALNUMBER
!
!  Nullify the stack that we will use to keep track of the boxes
!  that still need to be subdivided.
!
       NULLIFY(STACK)
       EXTDO2: DO
         CALL SPLITBOX(POINT,STEP,INT,ERR,NUMBER,       &
                       POINT1,STEP1,INT1,ERR1,NUMBER1,  &
                       POINT2,STEP2,INT2,ERR2,NUMBER2)
!
!  The procedure for the isolation of the zeros fails in case the
!  numerical evaluation of the integrals fails. 
!
         IF ( INFO .EQ. 3 ) EXIT EXTDO2
!
!  We process both the subboxes.
!
!  Box #1.
!
         IF ( NUMBER1 .GT. 0 .AND. NUMBER1 .LE. M ) THEN

           IF ( NR_ISOLATED_BOXES .EQ. 0 ) THEN
             START_LIST_OF_ISOLATED_BOXES%TOTAL_NUMBER_OF_ZEROS = NUMBER1 
             START_LIST_OF_ISOLATED_BOXES%POINT = POINT1 
             START_LIST_OF_ISOLATED_BOXES%STEP  = STEP1 
             NULLIFY(START_LIST_OF_ISOLATED_BOXES%NEXT)
             END_OF_LIST_OF_ISOLATED_BOXES => START_LIST_OF_ISOLATED_BOXES
             NR_ISOLATED_BOXES = NR_ISOLATED_BOXES + 1
           ELSE
             ALLOCATE(END_OF_LIST_OF_ISOLATED_BOXES%NEXT)
             END_OF_LIST_OF_ISOLATED_BOXES =>  &
               END_OF_LIST_OF_ISOLATED_BOXES%NEXT
             END_OF_LIST_OF_ISOLATED_BOXES%TOTAL_NUMBER_OF_ZEROS = NUMBER1
             END_OF_LIST_OF_ISOLATED_BOXES%POINT = POINT1
             END_OF_LIST_OF_ISOLATED_BOXES%STEP = STEP1
             NULLIFY(END_OF_LIST_OF_ISOLATED_BOXES%NEXT) 
             NR_ISOLATED_BOXES = NR_ISOLATED_BOXES + 1
           END IF

           IF ( ICON .GE. 3 ) THEN
             CALL APPROXIMATE(POINT1,STEP1,NUMBER1,ZEROS_BOX, &
                              FZEROS_BOX,MULT_BOX,NDIST_BOX)

             IF ( INFO .EQ. 4 ) THEN
               IF ( IFAIL .NE. 1 ) PRINT 7110
               CALL ERROR_EXIT
               RETURN
             END IF

             CALL REFINE(POINT1,STEP1,ZEROS_BOX, &
                         FZEROS_BOX,MULT_BOX,NDIST_BOX,REFINE_OK)

             DO J = 1, NDIST_BOX
               ZEROS_TN(DISTINCTNUMBER+J) = ZEROS_BOX(J)
               FZEROS_TN(DISTINCTNUMBER+J) = FZEROS_BOX(J)
               MULTIPLICITIES_TN(DISTINCTNUMBER+J) = MULT_BOX(J)
               REFINEMENT_OK_TN(DISTINCTNUMBER+J) = REFINE_OK(J)
               IF ( REFINE_OK(J) ) REFINEDNUMBER = REFINEDNUMBER + 1
             END DO
             DISTINCTNUMBER = DISTINCTNUMBER + NDIST_BOX
             END_OF_LIST_OF_ISOLATED_BOXES%DISTINCT_NUMBER_OF_ZEROS = &
               NDIST_BOX
           END IF

         ELSE IF ( NUMBER1 .GT. M ) THEN

           IF ( .NOT. ASSOCIATED(STACK) ) THEN
             ALLOCATE(STACK)
             NULLIFY(STACK%PREVIOUS)
           ELSE
             ALLOCATE(PTR)
             PTR%PREVIOUS => STACK
             STACK => PTR
           END IF
           STACK%NUMBER = NUMBER1
           STACK%POINT = POINT1
           STACK%STEP = STEP1
           STACK%INT = INT1 
           STACK%ERR = ERR1 

         END IF
!
!  Box #2.
!
         IF ( NUMBER2 .GT. 0 .AND. NUMBER2 .LE. M ) THEN

           IF ( NR_ISOLATED_BOXES .EQ. 0 ) THEN
             START_LIST_OF_ISOLATED_BOXES%TOTAL_NUMBER_OF_ZEROS = NUMBER2
             START_LIST_OF_ISOLATED_BOXES%POINT = POINT2
             START_LIST_OF_ISOLATED_BOXES%STEP  = STEP2
             NULLIFY(START_LIST_OF_ISOLATED_BOXES%NEXT)
             END_OF_LIST_OF_ISOLATED_BOXES => START_LIST_OF_ISOLATED_BOXES
             NR_ISOLATED_BOXES = NR_ISOLATED_BOXES + 1
           ELSE
             ALLOCATE(END_OF_LIST_OF_ISOLATED_BOXES%NEXT)
             END_OF_LIST_OF_ISOLATED_BOXES =>  &
               END_OF_LIST_OF_ISOLATED_BOXES%NEXT
             END_OF_LIST_OF_ISOLATED_BOXES%TOTAL_NUMBER_OF_ZEROS = NUMBER2
             END_OF_LIST_OF_ISOLATED_BOXES%POINT = POINT2
             END_OF_LIST_OF_ISOLATED_BOXES%STEP = STEP2
             NULLIFY(END_OF_LIST_OF_ISOLATED_BOXES%NEXT)
             NR_ISOLATED_BOXES = NR_ISOLATED_BOXES + 1
           END IF

           IF ( ICON .GE. 3 ) THEN
             CALL APPROXIMATE(POINT2,STEP2,NUMBER2,ZEROS_BOX, &
                              FZEROS_BOX,MULT_BOX,NDIST_BOX)

             IF ( INFO .EQ. 4 ) THEN
               IF ( IFAIL .NE. 1 ) PRINT 7110
               CALL ERROR_EXIT
               RETURN
             END IF

             CALL REFINE(POINT2,STEP2,ZEROS_BOX, &
                         FZEROS_BOX,MULT_BOX,NDIST_BOX,REFINE_OK)

             DO J = 1, NDIST_BOX
               ZEROS_TN(DISTINCTNUMBER+J) = ZEROS_BOX(J)
               FZEROS_TN(DISTINCTNUMBER+J) = FZEROS_BOX(J)
               MULTIPLICITIES_TN(DISTINCTNUMBER+J) = MULT_BOX(J)
               REFINEMENT_OK_TN(DISTINCTNUMBER+J) = REFINE_OK(J)
               IF ( REFINE_OK(J) ) REFINEDNUMBER = REFINEDNUMBER + 1
             END DO
             DISTINCTNUMBER = DISTINCTNUMBER + NDIST_BOX
             END_OF_LIST_OF_ISOLATED_BOXES%DISTINCT_NUMBER_OF_ZEROS = &
               NDIST_BOX
           END IF

         ELSE IF ( NUMBER2 .GT. M ) THEN

           IF ( .NOT. ASSOCIATED(STACK) ) THEN
             ALLOCATE(STACK)
             NULLIFY(STACK%PREVIOUS)
           ELSE
             ALLOCATE(PTR)
             PTR%PREVIOUS => STACK
             STACK => PTR
           END IF
           STACK%NUMBER = NUMBER2
           STACK%POINT = POINT2
           STACK%STEP = STEP2
           STACK%INT = INT2
           STACK%ERR = ERR2

         END IF
!
!  Take the next box from the stack.
!
         INTDO: DO
           IF ( ASSOCIATED(STACK) ) THEN
             POINT = STACK%POINT
             STEP = STACK%STEP
             INT = STACK%INT
             ERR = STACK%ERR 
             NUMBER = STACK%NUMBER
             STACK => STACK%PREVIOUS

             IF ( SQRT(STEP(1)**2+STEP(2)**2) .LT. NEWTONZ ) THEN
               PRINT 7200, POINT(1), POINT(2), STEP(1), STEP(2), NUMBER, M
               CYCLE INTDO
             END IF           

             CYCLE EXTDO2
           END IF
           EXIT EXTDO2
         END DO INTDO
       END DO EXTDO2
     END IF

!
!  Did the procedure for the isolation of the zeros fail?
!
     IF ( INFO .EQ. 3 ) THEN
       IF ( IFAIL .NE. 1 ) PRINT 7070
       CALL ERROR_EXIT
       RETURN
     END IF
!
!  We print the boxes that contain at most M zeros.
!
     IF ( VERBOSE ) THEN
       PRINT 7080, M, NR_ISOLATED_BOXES
       PRINT 7090
       PTR_ISOL_BOXES => START_LIST_OF_ISOLATED_BOXES
       DO J = 1, NR_ISOLATED_BOXES
         PRINT 7100, J, PTR_ISOL_BOXES%POINT(1), PTR_ISOL_BOXES%POINT(2), &
                        PTR_ISOL_BOXES%STEP(1), PTR_ISOL_BOXES%STEP(2),   &
                        PTR_ISOL_BOXES%TOTAL_NUMBER_OF_ZEROS
         PTR_ISOL_BOXES => PTR_ISOL_BOXES%NEXT
       END DO
     END IF
!
!  If ICON = 2, then we may stop.
!
     IF ( ICON .EQ. 2 ) RETURN
!
!  Set the output parameters.
!
     IF(ALLOCATED(ZEROS)) DEALLOCATE(ZEROS)
     IF(ALLOCATED(FZEROS)) DEALLOCATE(FZEROS)
     IF(ALLOCATED(MULTIPLICITIES)) DEALLOCATE(MULTIPLICITIES)
     IF(ALLOCATED(REFINEMENT_OK)) DEALLOCATE(REFINEMENT_OK)

     ALLOCATE(ZEROS(DISTINCTNUMBER))
     ALLOCATE(FZEROS(DISTINCTNUMBER))
     ALLOCATE(MULTIPLICITIES(DISTINCTNUMBER))
     ALLOCATE(REFINEMENT_OK(DISTINCTNUMBER))
  
     ZEROS(1:DISTINCTNUMBER) = ZEROS_TN(1:DISTINCTNUMBER)
     FZEROS(1:DISTINCTNUMBER) = FZEROS_TN(1:DISTINCTNUMBER)
     MULTIPLICITIES(1:DISTINCTNUMBER) = MULTIPLICITIES_TN(1:DISTINCTNUMBER)
     REFINEMENT_OK(1:DISTINCTNUMBER) = REFINEMENT_OK_TN(1:DISTINCTNUMBER)
!
!  Print the computed approximations for the zeros and related output.
!
     IF ( VERBOSE ) THEN
       IF ( ICON .EQ. 4 ) PRINT 7120, NR
       PRINT 7130

       PTR_ISOL_BOXES => START_LIST_OF_ISOLATED_BOXES
       PARTIALSUM_NUMBERDISTINCT =  &
              PTR_ISOL_BOXES%DISTINCT_NUMBER_OF_ZEROS
       K = 1
       REFINED_OK_AND_PRINTED = 0

       ICON_4_FLAG = .FALSE.
       EXTDO3: DO J = 1, NR_ISOLATED_BOXES
         PRINT 7140, J, PTR_ISOL_BOXES%DISTINCT_NUMBER_OF_ZEROS
         DO L = K, PARTIALSUM_NUMBERDISTINCT
           PRINT 7150, REAL(ZEROS(L),DP),  AIMAG(ZEROS(L)),    &
                       REAL(FZEROS(L),DP), AIMAG(FZEROS(L)),   &
                       MULTIPLICITIES(L)
           IF ( REFINEMENT_OK(L) ) REFINED_OK_AND_PRINTED   &
                                   = REFINED_OK_AND_PRINTED + 1 
           IF ( .NOT. REFINEMENT_OK(L) ) PRINT 7160
           IF ( ICON .EQ. 4 .AND. REFINED_OK_AND_PRINTED .EQ. NR ) THEN
             ICON_4_FLAG = .TRUE.
             EXIT EXTDO3
           END IF
         END DO
         IF ( J .LT. NR_ISOLATED_BOXES ) THEN
           PTR_ISOL_BOXES => PTR_ISOL_BOXES%NEXT
           K = PARTIALSUM_NUMBERDISTINCT + 1
           PARTIALSUM_NUMBERDISTINCT =             &
             PARTIALSUM_NUMBERDISTINCT +           &
             PTR_ISOL_BOXES%DISTINCT_NUMBER_OF_ZEROS
         END IF
       END DO EXTDO3
       IF ( (.NOT. ICON_4_FLAG) .AND.                                   &
            ((ICON .EQ. 3 .AND. DISTINCTNUMBER .NE. REFINEDNUMBER) .OR. &
             (ICON .EQ. 4 .AND. REFINEDNUMBER .LT. NR)) ) PRINT 7170
     END IF
!
!  If (FILES), then we create the files "zeros.dat", "fzeros.dat" 
!  and "mult.dat".
!
     IF ( FILES ) THEN
       OPEN(UNIT=1,FILE="zeros.dat",STATUS="REPLACE",ACTION="READWRITE")
       OPEN(UNIT=2,FILE="fzeros.dat",STATUS="REPLACE",ACTION="READWRITE")
       OPEN(UNIT=3,FILE="mult.dat",STATUS="REPLACE",ACTION="READWRITE")
       WRITE(UNIT=1,FMT=7180) ZEROS(1:DISTINCTNUMBER)
       WRITE(UNIT=2,FMT=7180) FZEROS(1:DISTINCTNUMBER)
       WRITE(UNIT=3,FMT=7190) MULTIPLICITIES(1:DISTINCTNUMBER)
       CLOSE(UNIT=1)
       CLOSE(UNIT=2)
       CLOSE(UNIT=3)
     END IF
!
!  Last statements of ZEAL.
!
     DEALLOCATE(ZEROS_TN)
     DEALLOCATE(FZEROS_TN)
     DEALLOCATE(MULTIPLICITIES_TN)
     DEALLOCATE(REFINEMENT_OK_TN)

     RETURN


7000 FORMAT ( /3X, 'This is ZEAL. Version of June 1999.',           &
             //3X, 'Input:',                                        &
             //3X, 'LV     =   ', G22.15, G23.15,                   &
              /3X, 'H      =   ', G22.15, G23.15,                   &
             //3X, 'M      =   ', I2,                               &
              /3X, 'NR     =   ', I2,                               &
              /3X, 'ICON   =    ', I1,                              &
             //3X, 'FILES  =    ', L1,                              &
             //3X, 63('=')                                          &
             //3X, 'Results:')
7010 FORMAT ( /3X, 'This is ZEAL. Version of June 1999.',           &
             //3X, 'Input:',                                        &
             //3X, 'LV     =   ', G22.15, G23.15,                   &
              /3X, 'H      =   ', G22.15, G23.15,                   &
             //3X, 'M      =   ', I2,                               &
              /3X, 'ICON   =    ', I1,                              &
             //3X, 'FILES  =    ', L1,                              &
             //3X, 63('=')                                          &
             //3X, 'Results:')
7020 FORMAT (/3X, '*** The procedure for the calculation of the ',  &
                  'total number of',                                &
             /3X, '    zeros has failed ***'//)
7030 FORMAT (/3X, 'The initial box has been expanded.',             &
             /3X, 'Zeros close to the boundary of the given box ',  &
                  'have been detected.')
7040 FORMAT (/3X, 'The requested number of zeros exceeds the ',     &
                  'total number of zeros',                          &
             /3X, 'that lie inside the box. ',                      &
                  'The value of NR is changed.')
7050 FORMAT (/3X, 'The following box has been considered:',         &
            //5X, 'LV = ', G22.15, G24.15,                          &
             /5X, 'H  = ', G22.15, G24.15)
7060 FORMAT (/3X, 'Total number of zeros inside this box',9X,'= ',I4)
7070 FORMAT (/3X, '*** The procedure for the isolation of the ',    &
                  'zeros has failed ***'//)
7080 FORMAT (/3X, 'Number of boxes containing at most',I4,          &
                  '  zeros = ',I4)
7090 FORMAT (/3X, 'These boxes are given by:',/3X, 24('-'))
7100 FORMAT (/2X, I4,')  LV = ',G22.15, G24.15,                     &
                   /9X, 'H  = ',G22.15, G24.15,                     &
                   /9X, 'Total number of zeros inside this box = ',I4)
7110 FORMAT (/3X, '*** The procedure for the computation of the ',  &
                  'zeros has failed ***'//)
7120 FORMAT (/3X, 'Requested number of mutually distinct zeros = ',I4)
7130 FORMAT (/3X, 'Final approximations for the zeros and ',        &
                  'verification:',                                  &
             /3X, 51('-')) 
7140 FORMAT (/2X, I4,')  Number of mutually distinct zeros inside ',&
                  'this box',3X,'= ',I4)
7150 FORMAT (/9X, 'z    = (', G22.15, ',', G22.15, ' )',            &
             /9X, 'f(z) = (', G22.15, ',', G22.15, ' )',            &
             /9X, 'multiplicity = ', I4)
7160 FORMAT (/9X, '*** Warning: Refinement procedure has failed ',  &
                  'for this zero ***')
7170 FORMAT (/3X, '*** Warning: Unable to refine the requested ',   &
                  'number of zeros ***')
7180 FORMAT (2X,G22.15,'     ',G22.15)
7190 FORMAT (2X,I5)
7200 FORMAT (/3X, '*** Error message from the subroutine ZEAL ***', &
             /3X, '*** The isolation procedure for the box',        &
             /3X, '      LV = ',G22.15, G24.15,                     &
             /3X, '      H  = ',G22.15, G24.15,                     &
             /3X, '*** has failed. The total number of zeros ',     &
                  'inside this box = ',I5,                          &
             /3X, '*** whereas M = ',I5,                            &
             /3X, '*** This box is too small. You may want to ',    &
                  'increase the value of M',                        &
             /3X, '*** or decrease the value of NEWTONZ.')

     END SUBROUTINE ZEAL


!-----------------------------------------------------------------------
     SUBROUTINE CHECK_INPUT
!-----------------------------------------------------------------------
!**PURPOSE
!  Check the user's input.  
!  
!-----------------------------------------------------------------------
!  Local variables
!
     LOGICAL :: HFLAG, MFLAG, ICONFLAG, NRFLAG, VALREGFLAG, EPSSTOPFLAG
     REAL(KIND=DP) :: EPSMCH, DLAMCH

!-----------------------------------------------------------------------
     INTRINSIC MAX, SQRT
     EXTERNAL DLAMCH 

     HFLAG = .FALSE.
     IF ( H(1) .LE. ZERO .OR. H(2) .LE. ZERO ) THEN
       HFLAG = .TRUE.
       INFO = 0
     END IF

     MFLAG = .FALSE.
     IF ( M .LE. 0 ) THEN
       MFLAG = .TRUE.
       INFO = 0
     END IF

     ICONFLAG = .FALSE.
     IF ( ICON .LT. 1 .OR. ICON .GT. 4 ) THEN
       ICONFLAG = .TRUE.
       INFO = 0
     END IF

     NRFLAG = .FALSE.
     IF ( ICON .EQ. 4 .AND. NR .LE. 0 ) THEN
       NRFLAG = .TRUE.
       INFO = 0
     END IF

     VALREGFLAG = .FALSE.
     IF ( .NOT. VALREG(LV,H) ) THEN
       VALREGFLAG = .TRUE.
       INFO = 0
     END IF

     EPSSTOPFLAG = .FALSE.
     IF ( EPS_STOP .LE. ZERO ) THEN
       EPSSTOPFLAG = .TRUE.
       INFO = 0
     END IF

     IF ( INFO .EQ. 0 ) THEN
       IF ( VERBOSE .AND. IFAIL .NE. 1 ) THEN
         PRINT 8000
         IF ( HFLAG ) PRINT 8010
         IF ( MFLAG ) PRINT 8020
         IF ( ICONFLAG ) PRINT 8030
         IF ( NRFLAG ) PRINT 8040
         IF ( VALREGFLAG ) PRINT 8050 
         IF ( EPSSTOPFLAG ) PRINT 8060
       END IF
     END IF
!
!  Compute the machine precision
!
     EPSMCH = DLAMCH('P')
 !
!  Determine if the input values of NUMABS and NUMREL are proper.
!
     IF ( MAX(NUMABS,NUMREL) .LT. 100.0_DP*EPSMCH  .OR.   &
          MAX(NUMABS,NUMREL) .GT. 0.1_DP                ) THEN
       NUMABS = 0.07_DP
       NUMREL = 0.0_DP
       IF ( VERBOSE .AND. IFAIL .NE. 1 ) PRINT 8070 
     END IF
!
!  Determine if the input values of NEWTONZ and NEWTONF are proper.
!
      IF ( NEWTONZ .GT. 0.001_DP .OR.                  &
           NEWTONZ .LT. SQRT(EPSMCH)/100.0_DP .OR.     &
           NEWTONF .LT. 10.0_DP*EPSMCH               ) THEN
        NEWTONZ = SQRT(EPSMCH)/100.0_DP
        NEWTONF = 10.0_DP*EPSMCH
        IF ( VERBOSE .AND. IFAIL .NE. 1 ) PRINT 8080, NEWTONZ, NEWTONF
      END IF

8000 FORMAT (/3X, '*** Improper input parameters ***')
8010 FORMAT (3X, 'The value of H(1) or H(2) is negative or zero')
8020 FORMAT (3X, 'The value of M is negative or zero')
8030 FORMAT (3X, 'The value of ICON is less than 1 or larger than 4')
8040 FORMAT (3X, 'The value of NR is negative or zero')
8050 FORMAT (3X, 'The function is not analytic in the region specified',  &
                 ' by LV and H') 
8060 FORMAT (3X, 'The value of EPS_STOP is negative')
8070 FORMAT (/3X,'WARNING -- Improper input value of NUMABS or NUMREL',   &
             /3X,'I will use the default values',                         &
                 ' NUMABS = 0.07 and NUMREL = 0.0')
8080 FORMAT (/3X,'WARNING -- Improper input value of NEWTONZ or NEWTONF', &
             /3X,'I will use the following values: ',                     &
             /3X,'  NEWTONZ = ',G22.15,                                   &
             /3X,'  NEWTONF = ',G22.15)

     END SUBROUTINE CHECK_INPUT

!-----------------------------------------------------------------------
     SUBROUTINE ERROR_EXIT
!-----------------------------------------------------------------------
!**PURPOSE
!  Check the possibility of error exit.
!  
!-----------------------------------------------------------------------

     IF ( IFAIL .EQ. 0 ) THEN
       STOP
     ELSE
       FAIL = .TRUE.
       INFORMATION = INFO
       RETURN
     END IF 

     END SUBROUTINE ERROR_EXIT


END MODULE Zout 
