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
SUBROUTINE Assembly_Diagnostics(ai_ijg, ai_kg,ar_weight,nstep)
IMPLICIT NONE
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ai_ijg, ai_kg  , nstep
REAL(KIND=8) ::ar_weight

IF(nstep .ge. 0) THEN
 
           Diagnostics(1,nstep+1) = Diagnostics(1,nstep+1) + &
		  & VarN_0(1, ai_ijg, 1) * &
		  & ar_weight

           Diagnostics(2,nstep+1) = Diagnostics(2,nstep+1) + &
		  & VarN_0(1, ai_ijg, 1) * &
		  & VarN_0(1, ai_ijg, 1) * &
		  & ar_weight

           Diagnostics(3,nstep+1) = Diagnostics(3,nstep+1) + &
                  & VarN_R(1, ai_ijg, 1) * &
		  & VarN_R(1, ai_ijg, 1) * &
		  & ar_weight  + &
		  & VarN_Z(1, ai_ijg, 1) * &
		  & VarN_Z(1, ai_ijg, 1) * &
		  & ar_weight 

	   !... Diff norms
           Diagnostics(4,nstep+1) = Diagnostics(4,nstep+1) + &
                & ( VarN_0(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 1) ) * &
		  & ( VarN_0(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 1) ) * &
		  & ar_weight 

          Diagnostics(5,nstep+1) = Diagnostics(5,nstep+1) + &
		  & ( VarN_R(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 2) ) * &
		  & ( VarN_R(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 2) ) * &
		  & ar_weight + &
		  & ( VarN_Z(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 3) ) * &
		  & ( VarN_Z(1, ai_ijg, 1) - Sol_analytical_2D(ai_ijg, 1, 3) ) * &
		  & ar_weight 

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
  
     Diagnostics(2,nstep+1) = SQRT(Diagnostics(2,nstep+1))
     Diagnostics(3,nstep+1) = SQRT(Diagnostics(3,nstep+1))

     Diagnostics(4,nstep+1) = SQRT(Diagnostics(4,nstep+1))
     Diagnostics(5,nstep+1) = SQRT(Diagnostics(5,nstep+1))
     
     PRINT *, "======      Masse ======"
     PRINT *,'Masse :',Diagnostics(1,nstep+1)
     PRINT *, "======      NORMS ======"
     PRINT *,'Norm L2 :',Diagnostics(2,nstep+1)
     PRINT *,'Norm H1 :',Diagnostics(3,nstep+1)
     PRINT *, "====== DIFF-NORMS ======"
     PRINT *,'Error L2 :',Diagnostics(4,nstep+1)/Diagnostics(2,nstep+1)
     PRINT *,'Error H1 :',Diagnostics(5,nstep+1)/Diagnostics(3,nstep+1)
     PRINT *, "========================"
     
   
        IF(nstep .eq. 0) THEN
           write(li_file_stream_norm,*) '# ,time, Var_id, masse,  norm_L2, semi_norm_H1, diff_norm_L2, diff_semi_norm_H1'
        END IF
        lr_time=nstep*lr_dt
        write(li_file_stream_norm,*) lr_time, 1, Diagnostics(1,nstep+1), Diagnostics(2,nstep+1), &
             Diagnostics(3,nstep+1), Diagnostics(4,nstep+1),Diagnostics(5,nstep+1)
     END IF

  
  RETURN 
END SUBROUTINE Plot_Diagnostics
! .........................................

END MODULE DIAGNOSTICS_COMPUTATION
