MODULE JRK_ELEMENT_CONTRIBUTION 
   USE BLACKBOX_DEF
   USE GREENBOX_DEF
   USE JOREK_GLOB_DEF
   USE MATRIX_DEF
   USE JRK_TESTCASE
   IMPLICIT NONE

CONTAINS
    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_1(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs = ANALYTICAL_RHS_1(ao_BBox2Di)

       contribution         = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_1
    ! ............................................................... 



    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_100(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs = sin(Vi_R)*sin(Vi_Z)

       contribution         = Vi_0 * wVol * f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

     END SUBROUTINE RHS_for_Vi_100
    ! ............................................................... 



    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_2(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS_2(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_2
    ! ............................................................... 

    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_3(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS_3(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_3
    ! ............................................................... 

    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_4(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS_4(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_4
    ! ............................................................... 

    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_5(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS_5(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_5
    ! ............................................................... 

    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_6(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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

       Rhs_Contribution => ptr_matrix % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)

       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS_6(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(1)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi_6
    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi_ONE(ptr_matrix, ao_BBox2Di, ao_GBox2D)
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
       REAL(KIND=RK), DIMENSION(:), POINTER :: Rhs_Contribution
       Rhs_Contribution => ptr_matrix % Rhs_Contribution       
       f_rhs       = ANALYTICAL_RHS_ONE(ao_BBox2Di)
       

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)

       ! ... ADD L2 contribution
       

       contribution                = Vi_0 *wVol
       Rhs_Contribution(1)  =  contribution 
       ! ...

     END SUBROUTINE RHS_for_Vi_ONE
    ! ............................................................... 


    ! ............................................................... 
    SUBROUTINE Stiffness_Matrix_for_Vi_Vj(ptr_matrix, ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
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

       Matrix_Contribution => ptr_matrix % Matrix_Contribution

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

       Matrix_Contribution(1, 1) =  contribution 
       ! ...

    END SUBROUTINE Stiffness_Matrix_for_Vi_Vj
    ! ............................................................... 

    ! ............................................................... 
    SUBROUTINE Mass_Matrix_for_Vi_Vj(ptr_matrix, ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
    IMPLICIT NONE
       CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Dj
       TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
       REAL(KIND=RK)               :: contribution_Mass  
       INTEGER       :: ijg
       REAL(KIND=RK) :: wVol
       REAL(KIND=RK) :: Vi_0
       REAL(KIND=RK) :: Vj_0
       REAL(KIND=RK), DIMENSION(:,:), POINTER :: Matrix_Contribution

       Matrix_Contribution => ptr_matrix % Matrix_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vj_0  = ao_BBox2Dj % B_0(ijg)
       
       contribution_Mass = Vi_0 * Vj_0 * wVol
       Matrix_Contribution(1, 1) =  contribution_Mass
  
    END SUBROUTINE  Mass_Matrix_for_Vi_Vj
    ! ............................................................... 

END MODULE JRK_ELEMENT_CONTRIBUTION
  

