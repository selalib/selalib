MODULE ELEMENT_CONTRIBUTION 
   USE JOREK_GLOB_DEF
   USE INDICES_DEF
   USE TYPEDEF
   USE MODEL_DEF
   USE TESTCASE
   IMPLICIT NONE

CONTAINS
    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi()
    IMPLICIT NONE
       REAL(KIND=RK) :: f_rhs
       REAL(KIND=RK) :: contribution
  
       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS()

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(i_Vu_Rho)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi
    ! ............................................................... 
  
    ! ............................................................... 
    SUBROUTINE Matrix_for_Vi_Vj()
    IMPLICIT NONE
       REAL(KIND=RK) :: contribution
  
       ! ... ADD STIFFNESS CONTRIBUTION
       contribution = ( Vi_R*Vj_R + Vi_Z*Vj_Z ) * wVol

       Matrix_Contribution(i_Vu_Rho, i_Vu_Rho) =  contribution 
       ! ...

    END SUBROUTINE Matrix_for_Vi_Vj
    ! ............................................................... 

END MODULE ELEMENT_CONTRIBUTION 

