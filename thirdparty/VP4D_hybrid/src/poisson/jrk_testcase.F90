MODULE JRK_TESTCASE
  USE JOREK_GLOB_DEF
  USE BLACKBOX_DEF
  USE JRK_MODEL_DEF

  IMPLICIT NONE

CONTAINS

  !! Test case 1 :
  !!             - mesh: square, collela or square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Dirichet or periodic boundary condition
  !!             - f(x,y)= 8pi**2 * SIN(2pi*R)*SIN(2pi*Z)

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_1(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_1
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     ijg = ao_BBox2D%ijg

     lr_R   = ao_BBox2D%Xp_0(1,ijg)
     lr_Z   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     ANALYTICAL_RHS_1 = (lr_k1**2 + lr_k2**2) / (4.0d0 * pi)**2 * &
          SIN(lr_k1*lr_R/(4.0d0*pi))*SIN(lr_k2*lr_Z/(4.0d0*pi))
  
  END FUNCTION ANALYTICAL_RHS_1
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_1(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = SIN(lr_k1*lr_R/(4.0d0*pi))*SIN(lr_k2*lr_Z/(4.0d0*pi))
     ! ... u_R
     apr_v(1, 2) = lr_k1/(4.0d0*pi)*COS(lr_k1*lr_R/(4.0d0*pi))*SIN(lr_k2*lr_Z/(4.0d0*pi))
     ! ... u_Z     
     apr_v(1, 3) = lr_k2/(4.0d0*pi)*SIN(lr_k1*lr_R/(4.0d0*pi))*COS(lr_k2*lr_Z/(4.0d0*pi))
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_1
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_100(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = SIN(lr_R)*SIN(lr_Z)
     ! ... u_R
     apr_v(1, 2) = COS(lr_R)*SIN(lr_Z)
     ! ... u_Z                       
     apr_v(1, 3) = SIN(lr_R)*COS(lr_Z)
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_100
  ! ..................................................

  !! Test case 2 : validate this case

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_2(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_2
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     ijg = ao_BBox2D%ijg

     lr_R   = ao_BBox2D%Xp_0(1,ijg)
     lr_Z   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     ANALYTICAL_RHS_2 = 4.0 * ( lr_R**2 + lr_Z**2 ) * SIN ( 1.0 - lr_R**2 - lr_Z**2 ) &
		& + 4.0 * COS ( 1.0 - lr_R**2 - lr_Z**2 ) 
  
  END FUNCTION ANALYTICAL_RHS_2
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_2(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter

     lr_R = apr_x(1)
     lr_Z = apr_x(2)

     lr_k1 = 2.0 * PI * FLOAT(1)
     lr_k2 = 2.0 * PI * FLOAT(1)

     ! ... u
     apr_v(1, 1) = SIN ( 1.0 - lr_R**2 - lr_Z**2 )
     ! ... u_R
     apr_v(1, 2) = - 2.0 * lr_R * COS ( 1.0 - lr_R**2 - lr_Z**2 )
     ! ... u_Z
     apr_v(1, 3) = - 2.0 * lr_Z * COS ( 1.0 - lr_R**2 - lr_Z**2 )
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_2

  !! Test case 3 :
  !!             - mesh: square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Periodic boundary condition
  !!             - f(x,y)= 

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_3(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_3
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     ijg = ao_BBox2D%ijg

     lr_R   = ao_BBox2D%Xp_0(1,ijg)
     lr_Z   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     ANALYTICAL_RHS_3 = 0.01 * cos(2.0*PI *lr_R) 

  END FUNCTION ANALYTICAL_RHS_3
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_3(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = 0.01 / (4.0*PI**2) * cos(2.0*PI *lr_R) 
     ! ... u_R
     apr_v(1, 2) =- 0.01 / (2.0*PI) * sin(2.0*PI *lr_R) 
     ! ... u_Z
     apr_v(1, 3) = 0.0
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_3
  ! ..................................................
  
  !! Test case 4 :
  !!             - mesh: square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Periodic boundary condition
  !!             - f(x,y)= 8pi**2 * cos(2pi*R)*cos(2pi*Z)

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_4(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_4
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     ijg = ao_BBox2D%ijg

     lr_R   = ao_BBox2D%Xp_0(1,ijg)
     lr_Z   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     ANALYTICAL_RHS_4 =  (lr_k1**2 + lr_k2**2) * cos(lr_k1*lr_R)*cos(lr_k2*lr_Z) 
     
!     ANALYTICAL_RHS_4 = 0.01 * cos(2.0*PI *lr_R) !+1
  END FUNCTION ANALYTICAL_RHS_4
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_4(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
!     apr_v(1, 1) = (lr_k1**2 + lr_k2**2) * cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
!     apr_v(1, 1) = 1.0 + 0.01 * cos(0.5 *lr_R)
     ! ... u_R
     apr_v(1, 2) = -lr_k1*sin(lr_k1*lr_R)*cos(lr_k2*lr_Z)
     ! ... u_Z
     apr_v(1, 3) = -lr_k2*cos(lr_k1*lr_R)*sin(lr_k2*lr_Z)
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_4
  ! ..................................................

  
  !! Test case 5 :
  !!             - mesh: square-xper
  !!             - R0=Z0=1 (lenght square)
  !!             - Cartesian coordindate
  !!             - Periodic in x direction
  !!             - Homogeneous Dirichlet in y direction 
  !!             - f(x,y)= 

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_5(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_5
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     ! LOCAL
     REAL(KIND=RK) :: lr_x
     REAL(KIND=RK) :: lr_y
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_Rmin
     REAL(KIND=RK) :: lr_Rmax
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     lr_Rmin = 0.1
     lr_Rmax = 1.0

     ijg = ao_BBox2D%ijg

     lr_x   = ao_BBox2D%Xp_0(1,ijg)
     lr_y   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = PI * FLOAT(mi_mode_m1)
     lr_k2 = PI * FLOAT(mi_mode_n1)

     ANALYTICAL_RHS_5 =  (lr_k1**2 + lr_k2**2) * cos(lr_k1*lr_x)*sin(lr_k2*lr_y) 
  
  END FUNCTION ANALYTICAL_RHS_5
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_5(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_x
     REAL(KIND=RK) :: lr_y
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2

     lr_x = apr_x(1)
     lr_y = apr_x(2)

     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)

     ! ... u
     apr_v(1, 1) = cos(lr_k1*lr_x)*sin(lr_k2*lr_y)
     ! ... u_R
     apr_v(1, 2) = -lr_k1*sin(lr_k1*lr_x)*sin(lr_k2*lr_y)
     ! ... u_Z
     apr_v(1, 3) = lr_k2*cos(lr_k1*lr_x)*cos(lr_k2*lr_y)
     ! ...


  END SUBROUTINE ANALYTICAL_MODEL_5
  ! ..................................................


  
  !! Test case 6 :
  !!             - mesh: annular 
  !!             - R0=Z0=1 (lenght square)
  !!             - Cartesian coordindate
  !!             - Homogeneous Dirichlet BC 
  !!             - f(x,y)= 

  ! ..................................................
  FUNCTION ANALYTICAL_RHS_6(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_6
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     ! LOCAL
     REAL(KIND=RK) :: lr_x
     REAL(KIND=RK) :: lr_y
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_c
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_Rmin
     REAL(KIND=RK) :: lr_Rmax
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     lr_Rmin = 0.1
     lr_Rmax = 1.0

     ijg = ao_BBox2D%ijg

     lr_x   = ao_BBox2D%Xp_0(1,ijg)
     lr_y   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = PI * FLOAT(mi_mode_m1)
     lr_c = lr_k1  / (lr_Rmax**2 - lr_Rmin**2) 
    
     ANALYTICAL_RHS_6 = 4.0 * lr_c**2 * ( lr_x**2 + lr_y**2 ) * SIN( lr_c * (lr_x**2 + lr_y**2 - lr_Rmin**2) ) &
                    & - 4.0 * lr_c * COS( lr_c * (lr_x**2 + lr_y**2 - lr_Rmin**2) )
    
  
  END FUNCTION ANALYTICAL_RHS_6
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_6(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_x
     REAL(KIND=RK) :: lr_y
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_c
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_Rmin
     REAL(KIND=RK) :: lr_Rmax
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_Rmin = 0.1
     lr_Rmax = 1.0

     lr_x = apr_x(1)
     lr_y = apr_x(2)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = PI * FLOAT(mi_mode_m1)
     lr_c = lr_k1  / (lr_Rmax**2 - lr_Rmin**2) 

     ! ... u
     apr_v(1, 1) = SIN ( lr_c * (lr_x**2 + lr_y**2 - lr_Rmin**2) ) 
     ! ... u_R
     apr_v(1, 2) = 2.0 * lr_c * lr_x * COS ( lr_c * (lr_x**2 + lr_y**2 - lr_Rmin**2) ) 
     ! ... u_Z
     apr_v(1, 3) = 2.0 * lr_c * lr_y * COS ( lr_c * (lr_x**2 + lr_y**2 - lr_Rmin**2) ) 
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_6
  ! ..................................................
  ! ..................................................
  FUNCTION ANALYTICAL_RHS_ONE(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS_ONE
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     !LOCAL

         
     ANALYTICAL_RHS_ONE = 1.0    
  
   END FUNCTION ANALYTICAL_RHS_ONE
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_MASS_1(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
    ! testcase : accumulation without mesh rescaling
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = - 0.01 * cos(2*PI*lr_R)
     ! ... u_R
     apr_v(1, 2) =   0.01*2*PI*sin(2*PI*lr_R)
     ! ... u_Z
     apr_v(1, 3) = 0.0
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_MASS_1
  ! ..................................................
  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_MASS_2(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
    ! testcase : accumulation with rescaling 
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
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

     ! ... u
     apr_v(1, 1) = - 0.01 * cos(lr_R)
     ! ... u_R
     apr_v(1, 2) =   0.01 * sin(lr_R)
     ! ... u_Z
     apr_v(1, 3) = 0.0
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_MASS_2
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_MASS_3(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     CLASS(DEF_FIELD_2D), POINTER                       :: self
     INTEGER                                            :: n_variable
     INTEGER                                            :: n_dimension
     REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
     REAL(KIND=8), DIMENSION(10)                        :: apr_info
     INTEGER, DIMENSION(10)                             :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter

     lr_R = apr_x(1)
     lr_Z = apr_x(2)

     lr_k1 = 2.0 * PI !* FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI !* FLOAT(mi_mode_n1)

     ! ... u
     apr_v(1, 1) =- SIN(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
     ! ... u_R
     apr_v(1, 2) =- lr_k1*COS(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
     ! ... u_Z
     apr_v(1, 3) =- lr_k2*SIN(lr_k1*lr_R)*COS(lr_k2*lr_Z)
     ! ...

  END SUBROUTINE ANALYTICAL_MODEL_MASS_3
  ! ..................................................

END MODULE JRK_TESTCASE 

