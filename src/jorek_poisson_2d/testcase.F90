MODULE TESTCASE 
  USE JOREK_GLOB_DEF
  USE INDICES_DEF
  USE TYPEDEF
  USE MODEL_DEF
  IMPLICIT NONE

CONTAINS

  !! Test case 1 :
  !!             - mesh: square, collela or square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Dirichet or periodic boundary condition
  !!             - f(x,y)= 8pi**2 * SIN(2pi*R)*SIN(2pi*Z)

   !! Test case 2 : validate this case

  !! Test case 3 : validate this case
  
  !! Test case 4 :
  !!             - mesh: square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Periodic boundary condition
  !!             - f(x,y)= 8pi**2 * cos(2pi*R)*cos(2pi*Z)

  
  FUNCTION ANALYTICAL_RHS()
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS 
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     INTEGER       :: li_testcase
     INTEGER       :: ierr
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    

     lr_R   = Xp_0(1,ijg)
     lr_Z   = Xp_0(2,ijg)
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     CALL JOREK_Param_GETInt(INT_TESTCASE_ID, li_testcase, ierr)
     IF(li_testcase .eq. 1) THEN
        ANALYTICAL_RHS = (lr_k1**2 + lr_k2**2) * SIN(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
     ENDIF    
     IF(li_testcase .eq. 2) THEN
        ANALYTICAL_RHS = 4.0 * ( lr_R**2 + lr_Z**2 ) * SIN ( 1.0 - lr_R**2 - lr_Z**2 ) &
		& + 4.0 * COS ( 1.0 - lr_R**2 - lr_Z**2 ) 
     ENDIF
     IF(li_testcase .eq. 3) THEN
        ANALYTICAL_RHS =8*lr_Z**2/lr_a**2 - 4 + 2*(2*lr_R - 2*lr_R0)**2/lr_a**2 &
		&+ 4*(lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2 &
		& + 4*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2 
     ENDIF

     ! periodic but not dirichet homogeneous
      IF(li_testcase .eq. 4) THEN
        ANALYTICAL_RHS = (lr_k1**2 + lr_k2**2) * cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
     ENDIF 
     
     
  
  END FUNCTION ANALYTICAL_RHS
  ! ..................................................
  
  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL(apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     REAL(KIND=8), DIMENSION(n_dimension), INTENT(IN)  :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1), INTENT(OUT) :: apr_v
     REAL(KIND=8), DIMENSION(10), INTENT(IN)  :: apr_info
     INTEGER, DIMENSION(10), INTENT(IN) :: api_info
     INTEGER                            :: n_variable,n_dimension
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     INTEGER       :: li_testcase
     INTEGER       :: ierr
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter

     lr_R = apr_x(1)
     lr_Z = apr_x(2)

     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)

     CALL JOREK_Param_GETInt(INT_TESTCASE_ID, li_testcase, ierr)
     IF(li_testcase .eq. 1) THEN
        ! ... u
        apr_v(1, 1) = SIN(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
        ! ... u_R
        apr_v(1, 2) = lr_k1*COS(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
        ! ... u_Z
        apr_v(1, 3) = lr_k2*SIN(lr_k1*lr_R)*COS(lr_k2*lr_Z)
        ! ...
     ENDIF    
     IF(li_testcase .eq. 2) THEN
        ! ... u
        apr_v(1, 1) = SIN ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ... u_R
        apr_v(1, 2) = - 2.0 * lr_R * COS ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ... u_Z
        apr_v(1, 3) = - 2.0 * lr_Z * COS ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ...
     ENDIF
     IF(li_testcase .eq. 3) THEN
        ! ... u
        apr_v(1, 1) = (1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2)*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)
        ! ... u_R
        apr_v(1, 2) = (1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2)*(2*lr_R - 2*lr_R0) &
		& - (2*lr_R - 2*lr_R0)*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2
        ! ... u_Z
        apr_v(1, 3) = 2*lr_Z*(1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2) &
		& - 2*lr_Z*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2
        ! ...
     ENDIF

     ! periodic but not dirichet homogeneous
     IF(li_testcase .eq. 4) THEN
        ! ... u
        apr_v(1, 1) = cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
        ! ... u_R
        apr_v(1, 2) = -lr_k1*sin(lr_k1*lr_R)*cos(lr_k2*lr_Z)
        ! ... u_Z
        apr_v(1, 3) = -lr_k2*cos(lr_k1*lr_R)*sin(lr_k2*lr_Z)
        ! ...
     ENDIF    
     

  END SUBROUTINE ANALYTICAL_MODEL
  ! ..................................................
END MODULE TESTCASE 

