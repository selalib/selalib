!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Split_Module
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
!  This modules is mainly composed by two subroutines. Subroutine INBOX
!  calculates the number of zeros that are conatined inside a given box.
!  Subroutine SPLITBOX splits a given region into two subregions and
!  calcualtes the number of zeros within each subregion.
!
!**END PROLOGUE Split_Module
!-----------------------------------------------------------------------

MODULE Split_Module

     USE Precision_Module
     USE Integration_Input_Module
     USE Function_Input_Module
     USE Quad_Module
     USE Error_Module

     IMPLICIT NONE

     LOGICAL :: FAIL1, FAIL2, FAIL3, FAIL4

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
!     PRIVATE
     PUBLIC :: INBOX, SPLITBOX, TRAP1, TRAP2, TRAP3, TRAP4


CONTAINS


     SUBROUTINE INBOX(POINT,STEP,INT,ERR,NUMBER,FLAG)
!-----------------------------------------------------------------------
!**PURPOSE
!  Compute the total number NUMBER of zeros that lie inside the 
!  rectangular region specified by POINT and STEP. If necessary,
!  the edges of the region are shifted. 
!
!-----------------------------------------------------------------------
!  Parameters
!
     LOGICAL, INTENT(OUT) :: FLAG
     INTEGER, INTENT(OUT) :: NUMBER 
     REAL(KIND=DP), DIMENSION(2), INTENT(INOUT)  :: POINT, STEP
     REAL(KIND=DP), DIMENSION(6), INTENT(OUT)    :: INT, ERR

!-----------------------------------------------------------------------
!  Local variables needed by DQAG
!
     INTEGER                        :: KEY, IER, NEVAL, LAST, EDGE
     INTEGER, PARAMETER             :: LIMIT = 1000 
     INTEGER, PARAMETER             :: LENW = 4*LIMIT
     INTEGER, DIMENSION(LIMIT)      :: IWORK
     REAL(KIND=DP)                  :: RESULT, ABSERR
     REAL(KIND=DP), DIMENSION(LENW) :: WORK 

!-----------------------------------------------------------------------
!  Other local variables
!
     REAL(KIND=DP) :: TOLX, TOLY
     REAL(KIND=DP) :: DXOLD, DYOLD, DXLIM, DYLIM, DXINIT, DYINIT
     REAL(KIND=DP), DIMENSION(4) :: INTTMP, ERRTMP 

!-----------------------------------------------------------------------
     INTRINSIC NINT
     EXTERNAL DQAG

!
!  Specify which quadrature formula is to be used.
!
     KEY = 1

!  Specify the amount by which an edge will be shifted in case DQAG
!  reports that the numerical integration cannot be performed in a 
!  satisfactory way.
!
     TOLX = STEP(1)/100.0_DP
     TOLY = STEP(2)/100.0_DP
!
!  The numerical integration is along the edges of the box whose
!  geometry is described by the global variables X0, Y0, DX and DY.
!
     X0 = POINT(1)
     Y0 = POINT(2)
     DX = STEP(1)
     DY = STEP(2)
!
!  We store the initial stepsizes of the box and specify the
!  maximal amount by which the edges are allowed to be shifted.
!
     DXINIT = DX
     DYINIT = DY

     DXLIM = DX/16.0_DP
     DYLIM = DY/16.0_DP

     FLAG = .FALSE.

     EDGE = 1

     EXTDO: DO

       SELECT CASE (EDGE)
!
!  We start with the lower edge.
!
       CASE (1)

         INTDO1: DO
           CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,       &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
           IF ( IER .EQ. 0 ) EXIT INTDO1
           Y0 = Y0 - TOLY
           DY = DY + TOLY
           FLAG = .TRUE.
           POINT(2) = Y0
           STEP(2) = DY
           IF ( .NOT. VALREG(POINT,STEP) ) THEN
              INFO = 2
             RETURN
           END IF
           IF ( (DY - DYINIT) .GT. DYLIM ) THEN
             INFO = 2
             RETURN
           END IF
         END DO INTDO1

         EDGE = 2

         INTTMP(1) = RESULT
         ERRTMP(1) = ABSERR

         DXOLD = DX
!
!  Then the right edge.
!
       CASE (2)

         INTDO2: DO
           CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
           IF ( IER .EQ. 0 ) EXIT INTDO2
           DX = DX + TOLX
           FLAG = .TRUE.
           STEP(1) = DX
           IF ( .NOT. VALREG(POINT,STEP) ) THEN
             INFO = 2
             RETURN
           END IF
           IF ( (DX - DXINIT) .GT. DXLIM ) THEN
             INFO = 2
             RETURN
           END IF
         END DO INTDO2

         INTTMP(2) = RESULT
         ERRTMP(2) = ABSERR
!
!  If the right edge has been moved, then we have to calculate an 
!  additional integral along a short segment of the lower edge. 
!  If this fails, then the lower edge has to be moved.
!
         IF ( DXOLD .NE. DX ) THEN
           X0 = X0 + DXOLD 
           DX = DX - DXOLD

           CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,ABSERR,  &
                     NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

           X0 = POINT(1)
           DX = STEP(1)

           IF ( IER .GT. 0 ) THEN   
             Y0 = Y0 - TOLY
             DY = DY + TOLY
             POINT(2) = Y0
             STEP(2) = DY
             IF ( .NOT. VALREG(POINT,STEP) ) THEN
               INFO = 2
               RETURN
             END IF
             IF ( (DY - DYINIT) .GT. DYLIM ) THEN
               INFO = 2
               RETURN
             END IF
             EDGE = 1
             ELSE
               INTTMP(1) = INTTMP(1) + RESULT
               ERRTMP(1) = ERRTMP(1) + ABSERR
             END IF
         END IF

         IF ( EDGE .EQ. 2 ) THEN
           EDGE = 3
           DYOLD = DY
         END IF
!
! Then the upper edge.
!
       CASE (3)

         INTDO3: DO
           CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
           IF ( IER .EQ. 0 ) EXIT INTDO3
           DY = DY + TOLY
           FLAG = .TRUE.
           STEP(2) = DY
           IF ( .NOT. VALREG(POINT,STEP) ) THEN
             INFO = 2
             RETURN
           END IF
           IF ( (DY - DYINIT) .GT. DYLIM ) THEN
             INFO = 2
             RETURN
           END IF
         END DO INTDO3

         INTTMP(3) = RESULT
         ERRTMP(3) = ABSERR
!
!  If the upper edge has been moved, then we have to calculate an 
!  additional integral along a short segment of the right edge. 
!  If this fails, then the right edge has to be moved.
!
         IF ( DYOLD .NE. DY ) THEN 
           Y0 = Y0 + DYOLD
           DY = DY - DYOLD

           CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
      
           Y0 = POINT(2)
           DY = STEP(2)

           IF ( IER .GT. 0 ) THEN
             DXOLD = DX
             DX = DX + TOLX
             STEP(1) = DX
             IF ( .NOT. VALREG(POINT,STEP) ) THEN
               INFO = 2
               RETURN
             END IF
             IF ( (DX - DXINIT) .GT. DXLIM ) THEN
               INFO = 2
               RETURN
             END IF
             EDGE = 2
           ELSE
             INTTMP(2) = INTTMP(2) + RESULT
             ERRTMP(2) = ERRTMP(2) + ABSERR
           END IF
         END IF

         IF ( EDGE .EQ. 3 ) THEN
           EDGE = 4
           DXOLD = DX
         END IF
!
!  Then the left edge.
!
       CASE (4)

         INTDO4: DO
           CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
           IF ( IER .EQ. 0 ) EXIT INTDO4
           X0 = X0 - TOLX
           DX = DX + TOLX
           FLAG = .TRUE.
           POINT(1) = X0
           STEP(1) = DX
           IF ( .NOT. VALREG(POINT,STEP) ) THEN
             INFO = 2
             RETURN
           END IF
           IF ( (DX - DXINIT) .GT. DXLIM ) THEN
             INFO = 2
             RETURN
           END IF
         END DO INTDO4

         EDGE = 5

         INTTMP(4) = RESULT
         ERRTMP(4) = ABSERR
!
!  If the left edge has been moved, then we have to calculate an 
!  additional integral along a short segment of the upper edge
!  and also along a short segment of the lower edge. We start with
!  the upper edge. If this fails, then the upper edge has to be moved. 
!  Else, we proceed with the lower edge. If this fails, then the lower 
!  edge has to be moved.
!
         IF ( DXOLD .NE. DX ) THEN
           DX = DX - DXOLD

           CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

           DX = STEP(1)

           IF ( IER .GT. 0 ) THEN
             DYOLD = DY
             DY = DY + TOLY
             STEP(2) = DY
             IF ( .NOT. VALREG(POINT,STEP) ) THEN
               INFO = 2
               RETURN
             END IF
             IF ( (DY - DYINIT) .GT. DYLIM ) THEN
               INFO = 2
               RETURN
             END IF
             EDGE = 3
           ELSE
             INTTMP(3) = INTTMP(3) + RESULT
             ERRTMP(3) = ERRTMP(3) + ABSERR 

             DX = DX - DXOLD

             CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                       ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

             DX = STEP(1)

             IF ( IER .GT. 0 ) THEN
               Y0 = Y0 - TOLY
               DY = DY + TOLY
               POINT(2) = Y0
               STEP(2) = DY
               IF ( .NOT. VALREG(POINT,STEP) ) THEN
                 INFO = 2
                 RETURN
               END IF
               IF ( (DY - DYINIT) .GT. DYLIM ) THEN
                 INFO = 2
                 RETURN
               END IF
               EDGE = 1
             ELSE
               INTTMP(1) = INTTMP(1) + RESULT
               ERRTMP(1) = ERRTMP(1) + ABSERR 
             END IF
           END IF
         END IF
!
!  Exit case: all integral have been calculated.
!
       CASE (5)

         EXIT EXTDO

       END SELECT

     END DO EXTDO
!
!  We calculate the number of zeros inside the box.
!
     NUMBER = NINT(INTTMP(1) + INTTMP(2) + INTTMP(3) + INTTMP(4))
!
!  We check if the (estimated) error is sufficiently small.
!
     IF ( ERRTMP(1)+ERRTMP(2)+ERRTMP(3)+ERRTMP(4) .GT. INTACC ) THEN
       INFO = 2
       RETURN
     END IF
!
!  We calculate the integrals along the halves of the longest edges.
!     
     IF ( NUMBER .GE. 1 ) THEN
       IF ( STEP(2) .GT. STEP(1) ) THEN
         INT(1) = INTTMP(1)
         ERR(1) = ERRTMP(1)
         INT(4) = INTTMP(3)
         ERR(4) = ERRTMP(3)

         DY = DY/2.0_DP
         CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT(2) = RESULT
         ERR(2) = ABSERR
         INT(3) = INTTMP(2) - RESULT
         ERR(3) = ERRTMP(2) + ABSERR
         CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT(6) = RESULT
         ERR(6) = ABSERR
         INT(5) = INTTMP(4) - RESULT
         ERR(5) = ERRTMP(4) + ABSERR
       ELSE
         INT(3) = INTTMP(2)
         ERR(3) = ERRTMP(2)
         INT(6) = INTTMP(4)
         ERR(6) = ERRTMP(4)

         DX = DX/2.0_DP
         CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT(1) = RESULT
         ERR(1) = ABSERR
         INT(2) = INTTMP(1) - RESULT
         ERR(2) = ERRTMP(1) + ABSERR
         CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT(5) = RESULT
         ERR(5) = ABSERR
         INT(4) = INTTMP(3) - RESULT
         ERR(4) = ERRTMP(3) + ABSERR
       END IF
     END IF
     END SUBROUTINE INBOX
!-----------------------------------------------------------------------


     SUBROUTINE SPLITBOX(POINT,STEP,INT,ERR,NUMBER,       &
                         POINT1,STEP1,INT1,ERR1,NUMBER1,  &
                         POINT2,STEP2,INT2,ERR2,NUMBER2)
!-----------------------------------------------------------------------
!**PURPOSE
!
!  Given a rectangular region specified by the variables POINT
!  and STEP, and given that the total number of zeros that lie
!  inside this region is equal to NUMBER, split this region into
!  two subregions, POINT1/STEP1 and POINT2/STEP2, respectively,
!  and compute the total number NUMBER1 and NUMBER2 of zeros 
!  that lie inside these subregions.
!
!-----------------------------------------------------------------------
!  Parameters
!
     INTEGER, INTENT(IN)  :: NUMBER 
     INTEGER, INTENT(OUT) :: NUMBER1, NUMBER2
     REAL(KIND=DP), DIMENSION(2), INTENT(IN)   :: POINT, STEP
     REAL(KIND=DP), DIMENSION(6), INTENT(IN)   :: INT, ERR
     REAL(KIND=DP), DIMENSION(2), INTENT(OUT)  :: POINT1, STEP1 
     REAL(KIND=DP), DIMENSION(2), INTENT(OUT)  :: POINT2, STEP2
     REAL(KIND=DP), DIMENSION(6), INTENT(OUT)  :: INT1, ERR1
     REAL(KIND=DP), DIMENSION(6), INTENT(OUT)  :: INT2, ERR2

!-----------------------------------------------------------------------
!  Local variables needed by DQAG
!
     INTEGER                        :: KEY, IER, NEVAL, LAST
     INTEGER, PARAMETER             :: LIMIT = 50000
     INTEGER, PARAMETER             :: LENW = 4*LIMIT
     INTEGER, DIMENSION(LIMIT)      :: IWORK
     REAL(KIND=DP)                  :: RESULT, ABSERR
     REAL(KIND=DP), DIMENSION(LENW) :: WORK

!-----------------------------------------------------------------------
!  Other local variables
!
     REAL(KIND=DP) :: TOLX, TOLY
     REAL(KIND=DP) :: DXOLD, DYOLD, DXLIM, DYLIM, DXINIT, DYINIT
     REAL(KIND=DP), DIMENSION(4) :: INTTM1, ERRTM1
     REAL(KIND=DP), DIMENSION(4) :: INTTM2, ERRTM2

!-----------------------------------------------------------------------
     INTRINSIC NINT
     EXTERNAL DQAG

!
!  Specify which quadrature formula is to be used.
!
     KEY = 1
!
!  N.B.: The numerical integration is along the edges of the box whose 
!        geometry is described by the variables X0, Y0, DX and DY.
!
!  Specify the amount by which the inner edge will be shifted 
!  (to the right or up) in case DQAG reports that the numerical
!  integration cannot be performed in a satisfactory way.
!
     TOLX = STEP(1)/100.0_DP
     TOLY = STEP(2)/100.0_DP

     IF ( STEP(2) .GT. STEP(1) ) THEN
       POINT1(1) = POINT(1)
       POINT1(2) = POINT(2)
       STEP1(1) = STEP(1)
       STEP1(2) = STEP(2)/2.0_DP

       X0 = POINT1(1)
       Y0 = POINT1(2)
       DX = STEP1(1)
       DY = STEP1(2)
!
!  We store the initial stepsize along the Y-direction and specify the 
!  maximal amount by which the inner edge is allowed to be shifted up.
!
       DYINIT = DY
       DYLIM = DY/16.0_DP

       INTTM1(1) = INT(1)
       INTTM1(2) = INT(2)
       INTTM1(4) = INT(6)
       ERRTM1(1) = ERR(1)
       ERRTM1(2) = ERR(2)
       ERRTM1(4) = ERR(6)

       INTTM2(2) = INT(3)
       INTTM2(3) = INT(4)
       INTTM2(4) = INT(5)
       ERRTM2(2) = ERR(3)
       ERRTM2(3) = ERR(4)
       ERRTM2(4) = ERR(5)

       DYOLD = DY

       DO
         CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,       &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         IF ( IER .EQ. 0 ) EXIT
         DY = DY + TOLY
         STEP1(2) = DY
         IF ( (DY - DYINIT) .GT. DYLIM ) THEN
           INFO = 3
           RETURN
         END IF
       ENDDO

       INTTM1(3) = RESULT
       ERRTM1(3) = ABSERR

       INTTM2(1) = -RESULT
       ERRTM2(1) = ABSERR

       IF ( DYOLD .NE. DY ) THEN
         Y0 = Y0 + DYOLD
         DY = DY - DYOLD

         CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         INTTM1(2) = INTTM1(2) + RESULT
         ERRTM1(2) = ERRTM1(2) + ABSERR

         INTTM2(2) = INTTM2(2) - RESULT
         ERRTM2(2) = ERRTM2(2) + ABSERR

         CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         INTTM1(4) = INTTM1(4) + RESULT
         ERRTM1(4) = ERRTM1(4) + ABSERR

         INTTM2(4) = INTTM2(4) - RESULT
         ERRTM2(4) = ERRTM2(4) + ABSERR
        
         Y0 = POINT1(2)
         DY = STEP1(2)
       END IF

       POINT2(1) = POINT(1)
       POINT2(2) = POINT(2) + STEP1(2)
       STEP2(1) = STEP(1) 
       STEP2(2) = STEP(2) - STEP1(2) 
!
!  We calculate the number of zeros inside the upper and lower box.
!
       NUMBER1 = NINT(INTTM1(1) + INTTM1(2) + INTTM1(3) + INTTM1(4))
       NUMBER2 = NINT(INTTM2(1) + INTTM2(2) + INTTM2(3) + INTTM2(4))
!
!  We check if the (estimated) error is sufficiently small.
!
       IF ( ERRTM1(1)+ERRTM1(2)+ERRTM1(3)+ERRTM1(4) .GT. INTACC .OR.   &
            ERRTM2(1)+ERRTM2(2)+ERRTM2(3)+ERRTM2(4) .GT. INTACC ) THEN
         INFO = 3
         RETURN
       END IF
!
!  We check if the boxes add up.
!
       IF ( NUMBER .NE. NUMBER1 + NUMBER2 ) THEN
         INFO = 3
         RETURN
       END IF
     ELSE
       POINT1(1) = POINT(1)
       POINT1(2) = POINT(2)
       STEP1(1) = STEP(1)/2.0_DP
       STEP1(2) = STEP(2)

       X0 = POINT1(1)
       Y0 = POINT1(2)
       DX = STEP1(1)
       DY = STEP1(2)
!
!  We store the initial stepsize along the X-direction and specify
!  the maximal amount by which the inner edge is allowed to be
!  shifted to the right.
!
       DXINIT = DX
       DXLIM = DX/16.0_DP

       INTTM1(1) = INT(1)
       INTTM1(3) = INT(5)
       INTTM1(4) = INT(6)
       ERRTM1(1) = ERR(1)
       ERRTM1(3) = ERR(5)
       ERRTM1(4) = ERR(6)

       INTTM2(1) = INT(2)
       INTTM2(2) = INT(3)
       INTTM2(3) = INT(4)
       ERRTM2(1) = ERR(2)
       ERRTM2(2) = ERR(3)
       ERRTM2(3) = ERR(4)

       DXOLD = DX

       DO
         CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         IF ( IER .EQ. 0 ) EXIT
         DX = DX + TOLX
         STEP1(1) = DX
         IF ( (DX - DXINIT) .GT. DXLIM ) THEN
           INFO = 3
           RETURN
         END IF
       ENDDO

       INTTM1(2) = RESULT
       ERRTM1(2) = ABSERR

       INTTM2(4) = -RESULT
       ERRTM2(4) = ABSERR

       IF ( DXOLD .NE. DX ) THEN
         X0 = X0 + DXOLD
         DX = DX - DXOLD

         CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         INTTM1(1) = INTTM1(1) + RESULT
         ERRTM1(1) = ERRTM1(1) + ABSERR

         INTTM2(1) = INTTM2(1) - RESULT
         ERRTM2(1) = ERRTM2(1) + ABSERR

         CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         INTTM1(3) = INTTM1(3) + RESULT
         ERRTM1(3) = ERRTM1(3) + ABSERR

         INTTM2(3) = INTTM2(3) - RESULT
         ERRTM2(3) = ERRTM2(3) + ABSERR
        
         X0 = POINT1(1)
         DX = STEP1(1)
       END IF

       POINT2(1) = POINT(1) + STEP1(1)  
       POINT2(2) = POINT(2) 
       STEP2(1) = STEP(1) - STEP1(1) 
       STEP2(2) = STEP(2) 
!
!  We calculate the number of zeros inside the left and right box.
!
       NUMBER1 = NINT(INTTM1(1) + INTTM1(2) + INTTM1(3) + INTTM1(4))
       NUMBER2 = NINT(INTTM2(1) + INTTM2(2) + INTTM2(3) + INTTM2(4))
!
!  We check if the (estimated) error is sufficiently small.
!
       IF ( ERRTM1(1)+ERRTM1(2)+ERRTM1(3)+ERRTM1(4) .GT. INTACC .OR.   &
            ERRTM2(1)+ERRTM2(2)+ERRTM2(3)+ERRTM2(4) .GT. INTACC ) THEN
         INFO = 3
         RETURN
       END IF
!
!  We check if the boxes add up.
!
       IF ( NUMBER .NE. NUMBER1 + NUMBER2 ) THEN
         INFO = 3
         RETURN
       END IF
     END IF
!
!  We calculate the integrals along the halves of the longest edges.
!
     IF ( NUMBER1 .GE. 1 ) THEN
       DX = STEP1(1)
       DY = STEP1(2)
       IF ( STEP1(2) .GT. STEP1(1) ) THEN
         INT1(1) = INTTM1(1)
         ERR1(1) = ERRTM1(1)
         INT1(4) = INTTM1(3)
         ERR1(4) = ERRTM1(3)

         DY = DY/2.0_DP
         CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT1(2) = RESULT
         ERR1(2) = ABSERR
         INT1(3) = INTTM1(2) - RESULT
         ERR1(3) = ERRTM1(2) + ABSERR
         CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT1(6) = RESULT
         ERR1(6) = ABSERR
         INT1(5) = INTTM1(4) - RESULT
         ERR1(5) = ERRTM1(4) + ABSERR
       ELSE
         INT1(3) = INTTM1(2)
         ERR1(3) = ERRTM1(2)
         INT1(6) = INTTM1(4)
         ERR1(6) = ERRTM1(4)

         DX = DX/2.0_DP
         CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

         INT1(1) = RESULT
         ERR1(1) = ABSERR
         INT1(2) = INTTM1(1) - RESULT
         ERR1(2) = ERRTM1(1) + ABSERR
         CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT1(5) = RESULT
         ERR1(5) = ABSERR
         INT1(4) = INTTM1(3) - RESULT
         ERR1(4) = ERRTM1(3) + ABSERR
       END IF
     END IF

     IF ( NUMBER2 .GE. 1 ) THEN
       X0 = POINT2(1)
       Y0 = POINT2(2)
       DX = STEP2(1)
       DY = STEP2(2)
       IF ( STEP2(2) .GT. STEP2(1) ) THEN
         INT2(1) = INTTM2(1)
         ERR2(1) = ERRTM2(1)
         INT2(4) = INTTM2(3)
         ERR2(4) = ERRTM2(3)

         DY = DY/2.0_DP
         CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT2(2) = RESULT
         ERR2(2) = ABSERR
         INT2(3) = INTTM2(2) - RESULT
         ERR2(3) = ERRTM2(2) + ABSERR
         CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT2(6) = RESULT
         ERR2(6) = ABSERR
         INT2(5) = INTTM2(4) - RESULT
         ERR2(5) = ERRTM2(4) + ABSERR
       ELSE
         INT2(3) = INTTM2(2)
         ERR2(3) = ERRTM2(2)
         INT2(6) = INTTM2(4)
         ERR2(6) = ERRTM2(4)

         DX = DX/2.0_DP
         CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT2(1) = RESULT
         ERR2(1) = ABSERR
         INT2(2) = INTTM2(1) - RESULT
         ERR2(2) = ERRTM2(1) + ABSERR
         CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                   ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         INT2(5) = RESULT
         ERR2(5) = ABSERR
         INT2(4) = INTTM2(3) - RESULT
         ERR2(4) = ERRTM2(3) + ABSERR
       END IF
     END IF

     END SUBROUTINE SPLITBOX

!----------------------------------------------------------------------

     SUBROUTINE TRAP1()

       INTEGER                        :: KEY, IER, NEVAL, LAST, EDGE
       INTEGER, PARAMETER             :: LIMIT = 1000 
       INTEGER, PARAMETER             :: LENW = 4*LIMIT
       INTEGER, DIMENSION(LIMIT)      :: IWORK
       REAL(KIND=DP)                  :: RESULT, ABSERR
       REAL(KIND=DP), DIMENSION(LENW) :: WORK 



       KEY = 1
       FAIL1 = .FALSE.

       CALL DQAG(F1,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

       if (IER .NE. 0) then 
          FAIL1 = .TRUE.
       END if
     end subroutine TRAP1


!---------------------------------------------------------------------
     SUBROUTINE TRAP2()

       INTEGER                        :: KEY, IER, NEVAL, LAST, EDGE
       INTEGER, PARAMETER             :: LIMIT = 1000 
       INTEGER, PARAMETER             :: LENW = 4*LIMIT
       INTEGER, DIMENSION(LIMIT)      :: IWORK
       REAL(KIND=DP)                  :: RESULT, ABSERR
       REAL(KIND=DP), DIMENSION(LENW) :: WORK 



       KEY = 1
       FAIL2 = .FALSE.

       CALL DQAG(F2,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

       if (IER .NE. 0) then 
          FAIL2 = .TRUE.
       END if
     end subroutine TRAP2


!-------------------------------------------------------------------

     SUBROUTINE TRAP3()

       INTEGER                        :: KEY, IER, NEVAL, LAST, EDGE
       INTEGER, PARAMETER             :: LIMIT = 1000 
       INTEGER, PARAMETER             :: LENW = 4*LIMIT
       INTEGER, DIMENSION(LIMIT)      :: IWORK
       REAL(KIND=DP)                  :: RESULT, ABSERR
       REAL(KIND=DP), DIMENSION(LENW) :: WORK 



       KEY = 1
       FAIL3 = .FALSE.

       CALL DQAG(F3,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

       if (IER .NE. 0) then 
          FAIL3 = .TRUE.
       END if
     end subroutine TRAP3

!----------------------------------------------------------------------



     SUBROUTINE TRAP4()

       INTEGER                        :: KEY, IER, NEVAL, LAST, EDGE
       INTEGER, PARAMETER             :: LIMIT = 1000 
       INTEGER, PARAMETER             :: LENW = 4*LIMIT
       INTEGER, DIMENSION(LIMIT)      :: IWORK
       REAL(KIND=DP)                  :: RESULT, ABSERR
       REAL(KIND=DP), DIMENSION(LENW) :: WORK 



       KEY = 1
       FAIL4 = .FALSE.

       CALL DQAG(F4,ZERO,ONE,NUMABS,NUMREL,KEY,RESULT,      &
                     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

       if (IER .NE. 0) then 
          FAIL4 = .TRUE.
       END if
     end subroutine TRAP4

!------------------------------------------------------------------------------
     subroutine plot(Xmin,Xmax, Y, pas, tabu,tabv, tabux, tabuy, tabpsix)
       REAL (KIND=DP) :: Xmin, Xmax, Y
       REAL(KIND=DP) :: h
       REAL(KIND=DP), dimension(101) :: tabu, tabv, tabpsix, tabux, tabuy
       integer :: j, pas

       h = (Xmax - Xmin)/pas
       do j = 1, pas+1
          tabu(j)=U(Xmin + (j-1)*h , Y)
          tabv(j)=V(Xmin + (j-1)*h, Y)
          tabux(j)=UX(Xmin + (j-1)*h , Y)
          tabuy(j)=UY(Xmin + (j-1)*h, Y)
          tabpsix =PSIX(Xmin + (j-1)*h, Y)
       end do 
       
     end subroutine plot











!-----------------------------------------------------------------------
     FUNCTION U(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: U
     COMPLEX (KIND=DP)          :: F, DF
     INTRINSIC CMPLX, REAL

     CALL FDF(CMPLX(X,Y,DP),F,DF)
     U = REAL(F,DP)

     END FUNCTION U
!-----------------------------------------------------------------------
     FUNCTION V(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: V
     COMPLEX (KIND=DP)          :: F, DF
     INTRINSIC CMPLX, AIMAG

     CALL FDF(CMPLX(X,Y,DP),F,DF)
     V = AIMAG(F)

     END FUNCTION V
!-----------------------------------------------------------------------
     FUNCTION UX(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: UX 
     COMPLEX (KIND=DP)          :: F, DF
     INTRINSIC CMPLX, REAL

     CALL FDF(CMPLX(X,Y,DP),F,DF)
     UX = REAL(DF,DP)

     END FUNCTION UX
!-----------------------------------------------------------------------
     FUNCTION UY(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: UY 
     COMPLEX (KIND=DP)          :: F, DF
     INTRINSIC CMPLX, AIMAG

     CALL FDF(CMPLX(X,Y,DP),F,DF)
     UY = -AIMAG(DF)

     END FUNCTION UY
!-----------------------------------------------------------------------
     FUNCTION PSIX(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: PSIX, UU, VV, UUX, VVX

     UU = U(X,Y)
     VV = V(X,Y)
     UUX = UU/(UU**2 + VV**2)
     UUX =  UUX*(-UY(X,Y))
     VVX = VV/(UU**2 + VV**2)
     VVX = VVX*UX(X,Y)
     PSIX = UUX - VVX
!     PSIX = (UU*(-UY(X,Y)) - VV*UX(X,Y))/(UU**2 + VV**2)

     END FUNCTION PSIX
!-----------------------------------------------------------------------
     FUNCTION PSIY(X,Y)
     REAL (KIND=DP), INTENT(IN) :: X, Y
     REAL (KIND=DP)             :: PSIY, UU, VV, UUX, VVX

     UU = U(X,Y)
     VV = V(X,Y)
     UUX = UU/(UU**2 + VV**2)
     UUX =  UUX*(UX(X,Y))
     VVX = VV/(UU**2 + VV**2)
     VVX = VVX*UY(X,Y)
     PSIY = UUX-VVX
!     PSIY = (UU*UX(X,Y) - VV*UY(X,Y))/(UU**2 + VV**2)

     END FUNCTION PSIY
!-----------------------------------------------------------------------
     FUNCTION F1(T)
     REAL (KIND=DP), INTENT(IN) :: T
     REAL (KIND=DP)             :: F1, TWOPI
     INTRINSIC ATAN

     TWOPI = 8.0_DP*ATAN(1.0_DP)
     F1 = DX*PSIX(X0+T*DX,Y0)/TWOPI 

     END FUNCTION F1
!-----------------------------------------------------------------------
     FUNCTION F2(T)
     REAL (KIND=DP), INTENT(IN) :: T
     REAL (KIND=DP)             :: F2, TWOPI
     INTRINSIC ATAN      

     TWOPI = 8.0_DP*ATAN(1.0_DP)
     F2 = DY*PSIY(X0+DX,Y0+T*DY)/TWOPI

     END FUNCTION F2
!-----------------------------------------------------------------------
     FUNCTION F3(T)
     REAL (KIND=DP), INTENT(IN) :: T
     REAL (KIND=DP)             :: F3, TWOPI
     INTRINSIC ATAN      

     TWOPI = 8.0_DP*ATAN(1.0_DP)
     F3 = -DX*PSIX(X0+T*DX,Y0+DY)/TWOPI

     END FUNCTION F3
!-----------------------------------------------------------------------
     FUNCTION F4(T)
     REAL (KIND=DP), INTENT(IN) :: T
     REAL (KIND=DP)             :: F4, TWOPI
     INTRINSIC ATAN      

     TWOPI = 8.0_DP*ATAN(1.0_DP)
     F4 = -DY*PSIY(X0,Y0+T*DY)/TWOPI

     END FUNCTION F4
!-----------------------------------------------------------------------


END MODULE Split_Module
