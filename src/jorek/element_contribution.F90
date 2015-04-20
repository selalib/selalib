MODULE ELEMENT_CONTRIBUTION 
   USE JOREK_GLOB_DEF
   USE INDICES_DEF
   USE TYPEDEF
   USE MODEL_DEF
   USE TESTCASE
   USE MATRIX_DEF
   IMPLICIT NONE

CONTAINS
    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi(ptr_matrix, ao_BBox2Di, ao_GBox2D)
    IMPLICIT NONE
       CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
       ! LOCAL
       REAL(KIND=RK) :: f_rhs
       REAL(KIND=RK) :: contribution
       INTEGER       :: ijg
       REAL(KIND=RK) :: wVol
       REAL(KIND=RK) :: Vi_0
       REAL(KIND=RK) :: Vi_R
       REAL(KIND=RK) :: Vi_Z
       REAL(KIND=RK) :: Vi_RR
       REAL(KIND=RK) :: Vi_RZ
       REAL(KIND=RK) :: Vi_ZZ
       REAL(KIND=RK), DIMENSION(:), POINTER :: Rhs_Contribution

       Rhs_Contribution => ao_GBox2D % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)
  
       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(i_Vu_Rho)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi
    ! ............................................................... 
  
    ! ............................................................... 
    SUBROUTINE Matrix_for_Vi_Vj(ptr_matrix, ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
    IMPLICIT NONE
       CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Dj
       TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
       ! LOCAL
       REAL(KIND=RK) :: contribution
       INTEGER       :: ijg
       REAL(KIND=RK) :: wVol
       REAL(KIND=RK) :: Vi_0
       REAL(KIND=RK) :: Vi_R
       REAL(KIND=RK) :: Vi_Z
       REAL(KIND=RK) :: Vi_RR
       REAL(KIND=RK) :: Vi_RZ
       REAL(KIND=RK) :: Vi_ZZ
       REAL(KIND=RK) :: Vj_0
       REAL(KIND=RK) :: Vj_R
       REAL(KIND=RK) :: Vj_Z
       REAL(KIND=RK) :: Vj_RR
       REAL(KIND=RK) :: Vj_RZ
       REAL(KIND=RK) :: Vj_ZZ
       REAL(KIND=RK), DIMENSION(:,:), POINTER :: Matrix_Contribution

       Matrix_Contribution => ao_GBox2D % Matrix_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)    ; Vj_0  = ao_BBox2Dj % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)   ; Vj_R  = ao_BBox2Dj % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)   ; Vj_Z  = ao_BBox2Dj % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg) ; Vj_RR = ao_BBox2Dj % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg) ; Vj_RZ = ao_BBox2Dj % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg) ; Vj_ZZ = ao_BBox2Dj % B_x2x2(ijg)

       ! ... ADD STIFFNESS CONTRIBUTION
       contribution = ( Vi_R*Vj_R + Vi_Z*Vj_Z ) * wVol

       Matrix_Contribution(i_Vu_Rho, i_Vu_Rho) =  contribution 
       ! ...

    END SUBROUTINE Matrix_for_Vi_Vj
    ! ............................................................... 

END MODULE ELEMENT_CONTRIBUTION 

