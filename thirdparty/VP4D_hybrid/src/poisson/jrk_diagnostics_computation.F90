MODULE JRK_DIAGNOSTICS_COMPUTATION
  USE JOREK_GLOB_DEF
  USE INDICES_DEF
  USE TYPEDEF
  USE JOREK_PARAM
  USE JRK_MODEL_DEF
  USE MATRIX_DEF
  USE FIELD
  IMPLICIT NONE

CONTAINS

 ! .........................................
SUBROUTINE Assembly_Diagnostics(ptr_matrix, ao_BBox2D, ao_GBox2D,nstep)
IMPLICIT NONE
CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ijg, kg  , nstep
REAL(KIND=RK) ::wVol
REAL(kind = rk), DIMENSION(3,1)     :: my_uh
real(kind = rk)                     :: uh
real(kind = rk), DIMENSION(1,1:2)                          :: lpr_x
  IF(nstep .ge. 0) THEN
     
       ijg   = ao_BBox2D % ijg
       wVol  = DABS(ao_BBox2D % wVol(ijg))

       g_diagnostics(1,nstep+1) = g_diagnostics(1,nstep+1) + &
            & ao_GBox2D%VarN_0(1, ijg) * &
            & wVol

       g_diagnostics(2,nstep+1) = g_diagnostics(2,nstep+1) + &
            & ao_GBox2D%VarN_0(1, ijg) * &
            & ao_GBox2D%VarN_0(1, ijg) * &
            & wVol
      
       g_diagnostics(3,nstep+1) = g_diagnostics(3,nstep+1) + &
            & ao_GBox2D%VarN_x1(1, ijg) * &
            & ao_GBox2D%VarN_x1(1, ijg) * &
            & wVol  + &
            & ao_GBox2D%VarN_x2(1, ijg) * &
            & ao_GBox2D%VarN_x2(1, ijg) * &
            & wVol 
     
       g_diagnostics(4,nstep+1) = g_diagnostics(4,nstep+1) + &
            & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
            & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
            & wVol 
       
       g_diagnostics(5,nstep+1) = g_diagnostics(5,nstep+1) + &
            & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
            & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
            & wVol + &
            & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
            & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
            & wVol 
    END IF

  RETURN 
  
END SUBROUTINE Assembly_Diagnostics
! .........................................
SUBROUTINE Assembly_Diagnostics_Field_Evaluate(ptr_matrix, ao_BBox2D, ao_GBox2D,nstep)
IMPLICIT NONE
CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ijg, kg  , nstep
REAL(KIND=RK) ::wVol
REAL(kind = rk), DIMENSION(3,1)     :: my_field_eval
real(kind = rk)                     :: error_in_eval
real(kind = rk), DIMENSION(1,1:2)                          :: lpr_x
 IF(nstep .ge. 0) THEN
     
        ijg   = ao_BBox2D % ijg
        wVol  = DABS(ao_BBox2D % wVol(ijg))


    lpr_x(1,1) = ao_BBox2D% ptr_quad % opr_points(1,ijg)
    lpr_x(1,2) = ao_BBox2D% ptr_quad % opr_points(2,ijg)
    CALL FIELD_EVALUATE(ptr_field,&
         ao_BBox2D%oi_elmt_id,&
         lpr_x,&
         my_field_eval)

    g_diagnostics(1,nstep+1) =g_diagnostics(1,nstep+1) &
         +  (ao_GBox2D%VarN_x1(1, ijg)-my_field_eval(2,1))**2 &
         +  (ao_GBox2D%VarN_x2(1, ijg)-my_field_eval(3,1))**2
!    print *, ao_GBox2D%VarN_x1(1, ijg), ao_GBox2D%VarN_x2(1, ijg)
!    print *, my_field_eval(2,1) , my_field_eval(3,1)

   
   END IF
  RETURN 
  
END SUBROUTINE Assembly_Diagnostics_Field_Evaluate
  
SUBROUTINE Plot_Diagnostics(nstep)
  IMPLICIT NONE
  INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
  REAL(KIND=RK) :: lr_dt, lr_time
  INTEGER, PARAMETER :: file_stream_norm=333 
  INTEGER :: i_iter
  
  CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

    IF(nstep .ge. 0) THEN
  
     g_diagnostics(2,nstep+1) = SQRT(g_diagnostics(2,nstep+1))
     g_diagnostics(3,nstep+1) = SQRT(g_diagnostics(3,nstep+1))

     mr_e_l2_norm = g_diagnostics(3,nstep+1) 
     
     PRINT *, "======      Masse ======"
     PRINT *,'Masse :',g_diagnostics(1,nstep+1)
     PRINT *, "======      NORMS ======"
     PRINT *,'Norm L2 :',g_diagnostics(2,nstep+1)
     PRINT *,'Norm H1 :',g_diagnostics(3,nstep+1)
     PRINT *, "====== DIFF-NORMS ======"
     PRINT *,'Error L2 :',g_diagnostics(4,nstep+1)/g_diagnostics(2,nstep+1)
     PRINT *,'Error H1 :',g_diagnostics(5,nstep+1)/g_diagnostics(3,nstep+1)
     PRINT *, "========================"

!     PRINT *, "========================"
    open(unit=file_stream_norm, file='diagnostics.dat', status='unknown')
       DO i_iter = 1, nstep
          lr_time = i_iter * lr_dt
          write(file_stream_norm,*) lr_time, g_diagnostics(3,i_iter+1)
       END DO
    close(file_stream_norm)
!     PRINT *, "========================"

      END IF
  RETURN 
END SUBROUTINE Plot_Diagnostics
! SUBROUTINE Plot_Diagnostics_Eval_Field(nstep)
!   IMPLICIT NONE
!   INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
!   REAL(KIND=RK) :: lr_dt, lr_time
!   INTEGER, PARAMETER :: file_stream_norm=333 
!   INTEGER :: i_iter
  
!   CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

!     IF(nstep .ge. 0) THEN
  
!      g_diagnostics(1,nstep+1) = SQRT(g_diagnostics(1,nstep+1))
     
!      PRINT *,'Error  :',g_diagnostics(1,nstep+1)
     
!   END IF
!   RETURN 
! END SUBROUTINE Plot_Diagnostics_Eval_Field
! ! .........................................


! SUBROUTINE Diagnostics_write_grad_in_file(nstep)
!   IMPLICIT NONE
!   INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
!   REAL(KIND=RK) :: lr_dt, lr_time
!   INTEGER, PARAMETER :: file_stream_norm=333 
!   INTEGER :: i_iter
  
!   CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

!     IF(nstep .ge. 0) THEN
  
!      g_diagnostics(2,nstep+1) = SQRT(g_diagnostics(2,nstep+1))
!      g_diagnostics(3,nstep+1) = SQRT(g_diagnostics(3,nstep+1))

!      mr_e_l2_norm = g_diagnostics(3,nstep+1) 
! !     PRINT *, "========================"
!     open(unit=file_stream_norm, file='diagnostics.dat', status='unknown')
!        DO i_iter = 0, nstep
!           lr_time = i_iter * lr_dt
!           write(file_stream_norm,*) lr_time, (g_diagnostics(3,i_iter+1) )
!        END DO
!     close(file_stream_norm)

! !     PRINT *, "========================"   

!       END IF
!   RETURN 
! END SUBROUTINE Diagnostics_write_grad_in_file
! ! .........................................
SUBROUTINE write_diagnostics_errornorm_in_file(nstep)
  IMPLICIT NONE
  INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
  REAL(KIND=RK) :: lr_dt, lr_time
  INTEGER, PARAMETER :: file_stream_norm=333 
  INTEGER :: i_iter
  open(unit=file_stream_norm, file='diagnostics.dat', status='unknown')
  DO i_iter = 0, nstep
     lr_time = i_iter * lr_dt
     write(file_stream_norm,*) lr_time, (g_diagnostics(4,i_iter+1) )
  END DO
  close(file_stream_norm)
END SUBROUTINE write_diagnostics_errornorm_in_file
END MODULE JRK_DIAGNOSTICS_COMPUTATION
