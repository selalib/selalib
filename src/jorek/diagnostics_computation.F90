MODULE DIAGNOSTICS_COMPUTATION
  USE JOREK_GLOB_DEF
  USE INDICES_DEF
  USE TYPEDEF
  USE JOREK_PARAM
  USE MODEL_DEF
  USE TESTCASE
  IMPLICIT NONE

CONTAINS

 ! .........................................
SUBROUTINE Assembly_Diagnostics(ao_BBox2D, ao_GBox2D,nstep)
IMPLICIT NONE
TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ijg, kg  , nstep
REAL(KIND=RK) ::wVol

  IF(nstep .ge. 0) THEN

       ijg   = ao_BBox2D % ijg
       wVol  = ao_BBox2D % wVol(ijg)
 
       fem_model % opr_diagnostics(1,nstep+1) = fem_model % opr_diagnostics(1,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       fem_model % opr_diagnostics(2,nstep+1) = fem_model % opr_diagnostics(2,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       fem_model % opr_diagnostics(3,nstep+1) = fem_model % opr_diagnostics(3,nstep+1) + &
                  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & wVol  + &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & wVol 

       !... Diff norms
       fem_model % opr_diagnostics(4,nstep+1) = fem_model % opr_diagnostics(4,nstep+1) + &
                  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & wVol 

       fem_model % opr_diagnostics(5,nstep+1) = fem_model % opr_diagnostics(5,nstep+1) + &
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

SUBROUTINE Plot_Diagnostics(nstep)
  IMPLICIT NONE
  INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
  REAL(KIND=RK) :: lr_dt, lr_time
  
  CALL JOREK_Param_GETInt(INT_COORDINATES_POLOIDAL_ID,mi_COORDINATES_POLOIDAL,ierr)
  CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

    IF(nstep .ge. 0) THEN
  
     fem_model % opr_diagnostics(2,nstep+1) = SQRT(fem_model % opr_diagnostics(2,nstep+1))
     fem_model % opr_diagnostics(3,nstep+1) = SQRT(fem_model % opr_diagnostics(3,nstep+1))

     fem_model % opr_diagnostics(4,nstep+1) = SQRT(fem_model % opr_diagnostics(4,nstep+1))
     fem_model % opr_diagnostics(5,nstep+1) = SQRT(fem_model % opr_diagnostics(5,nstep+1))
     
     PRINT *, "======      Masse ======"
     PRINT *,'Masse :',fem_model % opr_diagnostics(1,nstep+1)
     PRINT *, "======      NORMS ======"
     PRINT *,'Norm L2 :',fem_model % opr_diagnostics(2,nstep+1)
     PRINT *,'Norm H1 :',fem_model % opr_diagnostics(3,nstep+1)
     PRINT *, "====== DIFF-NORMS ======"
     PRINT *,'Error L2 :',fem_model % opr_diagnostics(4,nstep+1)/fem_model % opr_diagnostics(2,nstep+1)
     PRINT *,'Error H1 :',fem_model % opr_diagnostics(5,nstep+1)/fem_model % opr_diagnostics(3,nstep+1)
     PRINT *, "========================"
     
   
        IF(nstep .eq. 0) THEN
           write(li_file_stream_norm,*) '# ,time, Var_id, masse,  norm_L2, semi_norm_H1, diff_norm_L2, diff_semi_norm_H1'
        END IF
        lr_time=nstep*lr_dt
        write(li_file_stream_norm,*) lr_time, 1, fem_model % opr_diagnostics(1,nstep+1), fem_model % opr_diagnostics(2,nstep+1), &
             fem_model % opr_diagnostics(3,nstep+1), fem_model % opr_diagnostics(4,nstep+1),fem_model % opr_diagnostics(5,nstep+1)
     END IF

  
  RETURN 
END SUBROUTINE Plot_Diagnostics
! .........................................

END MODULE DIAGNOSTICS_COMPUTATION
