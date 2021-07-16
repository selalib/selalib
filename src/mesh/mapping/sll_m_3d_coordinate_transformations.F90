!> @ingroup mesh
!> @brief
!> definition of analytical coordinate transformations and their jacobi matrix, inverse jacobi matrix and jacobian
!> the transformation should be defined in a way, that the jacobian is strictly positive as singular points have to be excluded 
!> @author
!> Benedikt Perse
module sll_m_3d_coordinate_transformations
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
       sll_p_pi, &
       sll_p_twopi

  implicit none

  public :: &
       sll_f_colella_x1,&
       sll_f_colella_x2,&
       sll_f_colella_x3,&
       sll_f_colella_jac11,&
       sll_f_colella_jac12,&
       sll_f_colella_jac13,&
       sll_f_colella_jac21,&
       sll_f_colella_jac22,&
       sll_f_colella_jac23,&
       sll_f_colella_jac31,&
       sll_f_colella_jac32,&
       sll_f_colella_jac33,&
       sll_f_colella_jacobian,&
       sll_f_colbound_x1,&
       sll_f_colbound_x2,&
       sll_f_colbound_x3,&
       sll_f_colbound_jac11,&
       sll_f_colbound_jac12,&
       sll_f_colbound_jac13,&
       sll_f_colbound_jac21,&
       sll_f_colbound_jac22,&
       sll_f_colbound_jac23,&
       sll_f_colbound_jac31,&
       sll_f_colbound_jac32,&
       sll_f_colbound_jac33,&
       sll_f_colbound_jacobian,&
       sll_f_parallelogram_x1,&
       sll_f_parallelogram_x2,&
       sll_f_parallelogram_x3,&
       sll_f_parallelogram_jac11,&
       sll_f_parallelogram_jac12,&
       sll_f_parallelogram_jac13,&
       sll_f_parallelogram_jac21,&
       sll_f_parallelogram_jac22,&
       sll_f_parallelogram_jac23,&
       sll_f_parallelogram_jac31,&
       sll_f_parallelogram_jac32,&
       sll_f_parallelogram_jac33,&
       sll_f_parallelogram_jacobian,&
       sll_f_orthogonal_x1,&
       sll_f_orthogonal_x2,&
       sll_f_orthogonal_x3,&
       sll_f_orthogonal_jac11,&
       sll_f_orthogonal_jac12,&
       sll_f_orthogonal_jac13,&
       sll_f_orthogonal_jac21,&
       sll_f_orthogonal_jac22,&
       sll_f_orthogonal_jac23,&
       sll_f_orthogonal_jac31,&
       sll_f_orthogonal_jac32,&
       sll_f_orthogonal_jac33,&
       sll_f_orthogonal_jacobian,&
       sll_f_polynomial_x1,&
       sll_f_polynomial_x2,&
       sll_f_polynomial_x3,&
       sll_f_polynomial_jac11,&
       sll_f_polynomial_jac12,&
       sll_f_polynomial_jac13,&
       sll_f_polynomial_jac21,&
       sll_f_polynomial_jac22,&
       sll_f_polynomial_jac23,&
       sll_f_polynomial_jac31,&
       sll_f_polynomial_jac32,&
       sll_f_polynomial_jac33,&
       sll_f_polynomial_jacobian,&
       sll_f_rotation_x1,&
       sll_f_rotation_x2,&
       sll_f_rotation_x3,&
       sll_f_rotation_jac11,&
       sll_f_rotation_jac12,&
       sll_f_rotation_jac13,&
       sll_f_rotation_jac21,&
       sll_f_rotation_jac22,&
       sll_f_rotation_jac23,&
       sll_f_rotation_jac31,&
       sll_f_rotation_jac32,&
       sll_f_rotation_jac33,&
       sll_f_rotation_jacobian,&
       sll_f_cylindrical_x1,&
       sll_f_cylindrical_x2,&
       sll_f_cylindrical_x3,&
       sll_f_cylindrical_xi1,&
       sll_f_cylindrical_xi2,&
       sll_f_cylindrical_xi3,&
       sll_f_cylindrical_jac11,&
       sll_f_cylindrical_jac12,&
       sll_f_cylindrical_jac13,&
       sll_f_cylindrical_jac21,&
       sll_f_cylindrical_jac22,&
       sll_f_cylindrical_jac23,&
       sll_f_cylindrical_jac31,&
       sll_f_cylindrical_jac32,&
       sll_f_cylindrical_jac33,&
       sll_f_cylindrical_jacobian,&
       sll_f_cylindrical_sqrt_x1,&
       sll_f_cylindrical_sqrt_x2,&
       sll_f_cylindrical_sqrt_x3,&
       sll_f_cylindrical_sqrt_xi1,&
       sll_f_cylindrical_sqrt_xi2,&
       sll_f_cylindrical_sqrt_xi3,&
       sll_f_cylindrical_sqrt_jac11,&
       sll_f_cylindrical_sqrt_jac12,&
       sll_f_cylindrical_sqrt_jac13,&
       sll_f_cylindrical_sqrt_jac21,&
       sll_f_cylindrical_sqrt_jac22,&
       sll_f_cylindrical_sqrt_jac23,&
       sll_f_cylindrical_sqrt_jac31,&
       sll_f_cylindrical_sqrt_jac32,&
       sll_f_cylindrical_sqrt_jac33,&
       sll_f_cylindrical_sqrt_jacobian,&
       sll_f_elliptical_x1,&
       sll_f_elliptical_x2,&
       sll_f_elliptical_x3,&
       sll_f_elliptical_jac11,&
       sll_f_elliptical_jac12,&
       sll_f_elliptical_jac13,&
       sll_f_elliptical_jac21,&
       sll_f_elliptical_jac22,&
       sll_f_elliptical_jac23,&
       sll_f_elliptical_jac31,&
       sll_f_elliptical_jac32,&
       sll_f_elliptical_jac33,&
       sll_f_elliptical_jacobian,&
       sll_f_Dshaped_singular_x1,&
       sll_f_Dshaped_singular_x2,&
       sll_f_Dshaped_singular_x3,&
       sll_f_Dshaped_singular_jac11,&
       sll_f_Dshaped_singular_jac12,&
       sll_f_Dshaped_singular_jac13,&
       sll_f_Dshaped_singular_jac21,&
       sll_f_Dshaped_singular_jac22,&
       sll_f_Dshaped_singular_jac23,&
       sll_f_Dshaped_singular_jac31,&
       sll_f_Dshaped_singular_jac32,&
       sll_f_Dshaped_singular_jac33,&
       sll_f_Dshaped_singular_jacobian,&
       sll_f_Dshaped_singular_pseudoinv11,&
       sll_f_Dshaped_singular_pseudoinv12,&
       sll_f_Dshaped_singular_pseudoinv13,&
       sll_f_Dshaped_singular_pseudoinv21,&
       sll_f_Dshaped_singular_pseudoinv22,&
       sll_f_Dshaped_singular_pseudoinv23,&
       sll_f_Dshaped_singular_pseudoinv31,&
       sll_f_Dshaped_singular_pseudoinv32,&
       sll_f_Dshaped_singular_pseudoinv33,&
       sll_f_rotated_singular_x1,&
       sll_f_rotated_singular_x2,&
       sll_f_rotated_singular_x3,&
       sll_f_rotated_singular_jac11,&
       sll_f_rotated_singular_jac12,&
       sll_f_rotated_singular_jac13,&
       sll_f_rotated_singular_jac21,&
       sll_f_rotated_singular_jac22,&
       sll_f_rotated_singular_jac23,&
       sll_f_rotated_singular_jac31,&
       sll_f_rotated_singular_jac32,&
       sll_f_rotated_singular_jac33,&
       sll_f_rotated_singular_jacobian,&
       sll_f_rotated_singular_pseudoinv11,&
       sll_f_rotated_singular_pseudoinv12,&
       sll_f_rotated_singular_pseudoinv13,&
       sll_f_rotated_singular_pseudoinv21,&
       sll_f_rotated_singular_pseudoinv22,&
       sll_f_rotated_singular_pseudoinv23,&
       sll_f_rotated_singular_pseudoinv31,&
       sll_f_rotated_singular_pseudoinv32,&
       sll_f_rotated_singular_pseudoinv33,&
       sll_f_toroidal_cylinder_x1,&
       sll_f_toroidal_cylinder_x2,&
       sll_f_toroidal_cylinder_x3,&
       sll_f_toroidal_cylinder_jac11,&
       sll_f_toroidal_cylinder_jac12,&
       sll_f_toroidal_cylinder_jac13,&
       sll_f_toroidal_cylinder_jac21,&
       sll_f_toroidal_cylinder_jac22,&
       sll_f_toroidal_cylinder_jac23,&
       sll_f_toroidal_cylinder_jac31,&
       sll_f_toroidal_cylinder_jac32,&
       sll_f_toroidal_cylinder_jac33,&
       sll_f_toroidal_cylinder_jacobian,&
       sll_f_spherical_x1,&
       sll_f_spherical_x2,&
       sll_f_spherical_x3,&
       sll_f_spherical_jac11,&
       sll_f_spherical_jac12,&
       sll_f_spherical_jac13,&
       sll_f_spherical_jac21,&
       sll_f_spherical_jac22,&
       sll_f_spherical_jac23,&
       sll_f_spherical_jac31,&
       sll_f_spherical_jac32,&
       sll_f_spherical_jac33,&
       sll_f_spherical_jacobian,&
       sll_f_identity_x1,&
       sll_f_identity_x2,&
       sll_f_identity_x3,&
       sll_f_identity_xi1,&
       sll_f_identity_xi2,&
       sll_f_identity_xi3,&
       sll_f_identity_jac11,&
       sll_f_identity_jac12,&
       sll_f_identity_jac13,&
       sll_f_identity_jac21,&
       sll_f_identity_jac22,&
       sll_f_identity_jac23,&
       sll_f_identity_jac31,&
       sll_f_identity_jac32,&
       sll_f_identity_jac33,&
       sll_f_affine_x1,&
       sll_f_affine_x2,&
       sll_f_affine_x3,&
       sll_f_affine_jac11,&
       sll_f_affine_jac12,&
       sll_f_affine_jac13,&
       sll_f_affine_jac21,&
       sll_f_affine_jac22,&
       sll_f_affine_jac23,&
       sll_f_affine_jac31,&
       sll_f_affine_jac32,&
       sll_f_affine_jac33,&
       sll_f_scaling_x1,&
       sll_f_scaling_x2,&
       sll_f_scaling_x3,&
       sll_f_scaling_xi1,&
       sll_f_scaling_xi2,&
       sll_f_scaling_xi3,&
       sll_f_scaling_jac11,&
       sll_f_scaling_jac12,&
       sll_f_scaling_jac13,&
       sll_f_scaling_jac21,&
       sll_f_scaling_jac22,&
       sll_f_scaling_jac23,&
       sll_f_scaling_jac31,&
       sll_f_scaling_jac32,&
       sll_f_scaling_jac33,&
       sll_f_scaling_jacobian,&
       sll_f_periodic_x1,&
       sll_f_periodic_x2,&
       sll_f_periodic_x3,&
       sll_f_periodic_jac11,&
       sll_f_periodic_jac12,&
       sll_f_periodic_jac13,&
       sll_f_periodic_jac21,&
       sll_f_periodic_jac22,&
       sll_f_periodic_jac23,&
       sll_f_periodic_jac31,&
       sll_f_periodic_jac32,&
       sll_f_periodic_jac33,&
       sll_f_periodic_jacobian
      
       

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  sll_real64, parameter :: a = 0.5_f64*sll_p_pi!sll_p_twopi !
  sll_real64, parameter :: b = sll_p_twopi !sll_p_pi !
contains

    ! ***************************************************************************
  !
  ! <b> "Colella transformation with Boundary"; </b>
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !  X1 = L1 * (\xi(1) + \alpha1 *  \xi(2) ) 
  !  X2 = L2 * \xi(2)
  !  X3 = L3 * \xi(3)  
  ! 
  ! The parameters are:
  !    - L1     = params(1)
  !    - L2     = params(2)
  !    - L3     = params(3)
  !    - alpha1 = params(4)
  !
  ! ***************************************************************************

  
  !> direct mapping
  function sll_f_parallelogram_x1( xi, params )
    sll_real64 :: sll_f_parallelogram_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1

    SLL_ASSERT(size(params) >= 4)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_parallelogram_x1 = L1*(xi(1)+alpha1*xi(2) ) 
  end function sll_f_parallelogram_x1

  !> direct mapping
  function sll_f_parallelogram_x2( xi, params )
    sll_real64 :: sll_f_parallelogram_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2

    SLL_ASSERT(size(params) >= 4)
    L2     = params(2)
    alpha2 = params(4)
    sll_f_parallelogram_x2 = L2*(xi(2))
  end function sll_f_parallelogram_x2

  !> direct mapping
  function sll_f_parallelogram_x3( xi, params )
    sll_real64 :: sll_f_parallelogram_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    
    SLL_ASSERT(size(params) >= 4)
    L3     = params(3)
    sll_f_parallelogram_x3 = L3*xi(3)
  end function sll_f_parallelogram_x3

  !> jacobian matrix
  function sll_f_parallelogram_jac11 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 4)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_parallelogram_jac11 = L1
  end function sll_f_parallelogram_jac11

  !> jacobian matrix
  function sll_f_parallelogram_jac12 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha1, L1
    
    SLL_ASSERT(size(params) >= 4)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_parallelogram_jac12 = L1*alpha1
  end function sll_f_parallelogram_jac12

  !> jacobian matrix
  function sll_f_parallelogram_jac13 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    
    sll_f_parallelogram_jac13 =0._f64
  end function sll_f_parallelogram_jac13

  !> jacobian matrix
  function sll_f_parallelogram_jac21 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha2, L2
    
    SLL_ASSERT(size(params) >= 4)
    L2     = params(2)
    alpha2 = params(4)
    sll_f_parallelogram_jac21 = 0._f64
  end function sll_f_parallelogram_jac21

  !> jacobian matrix
  function sll_f_parallelogram_jac22 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2
    
    SLL_ASSERT(size(params) >= 4)
    L2     = params(2)
    alpha2 = params(4)
    sll_f_parallelogram_jac22 =L2
  end function sll_f_parallelogram_jac22

  !> jacobian matrix
  function sll_f_parallelogram_jac23 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_parallelogram_jac23 = 0._f64
  end function sll_f_parallelogram_jac23

  !> jacobian matrix
  function sll_f_parallelogram_jac31 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_parallelogram_jac31 = 0._f64
  end function sll_f_parallelogram_jac31

  !> jacobian matrix
  function sll_f_parallelogram_jac32 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_parallelogram_jac32 = 0._f64
  end function sll_f_parallelogram_jac32

  !> jacobian matrix
  function sll_f_parallelogram_jac33 ( xi, params )
    sll_real64  :: sll_f_parallelogram_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 4)
    L3 = params(3)
    sll_f_parallelogram_jac33 = L3
  end function sll_f_parallelogram_jac33

    !> jacobian 
  function sll_f_parallelogram_jacobian ( xi, params )
    sll_real64  :: sll_f_parallelogram_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2, L3
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 4)
    L1 =     params(1)
    L2 =     params(2)
    L3 =     params(3)
    alpha1 = params(4)
    sll_f_parallelogram_jacobian = L3*L1*L2
  end function sll_f_parallelogram_jacobian

  
  ! ***************************************************************************
  !
  ! <b> "Colella transformation with Boundary"; </b>
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !  X1 = L1 * (\xi(1) + \alpha1 * sin(a * \xi(1)) * sin(b * \xi(2)) ) 
  !  X2 = L2 * (\xi(2) + \alpha2 * sin(a * \xi(1)) * sin(b * \xi(2)) ) 
  !  X3 = L3 * (\xi(3) + \alpha3 * sin(2*\pi * \xi(3)) ) 
  ! 
  ! The parameters are:
  !    - L1     = params(1)
  !    - L2     = params(2)
  !    - L3     = params(3)
  !    - alpha1 = params(4)
  !    - alpha2 = params(5)
  !    - alpha3 = params(6)
  !    - a = \pi/2
  !    - b = 2*\pi
  !
  ! ***************************************************************************

  
  !> direct mapping
  function sll_f_colbound_x1( xi, params )
    sll_real64 :: sll_f_colbound_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1

    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colbound_x1 = L1*(xi(1)+alpha1*sin(a*xi(1))*sin(b*xi(2)))
  end function sll_f_colbound_x1

  !> direct mapping
  function sll_f_colbound_x2( xi, params )
    sll_real64 :: sll_f_colbound_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2

    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colbound_x2 = L2*(xi(2)+alpha2*sin(a*xi(1))*sin(b*xi(2)))
  end function sll_f_colbound_x2

  !> direct mapping
  function sll_f_colbound_x3( xi, params )
    sll_real64 :: sll_f_colbound_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L3     = params(3)
    alpha3 = params(6)
    sll_f_colbound_x3 = L3*(xi(3)+alpha3*sin(sll_p_twopi*xi(3)))
  end function sll_f_colbound_x3

  !> jacobian matrix
  function sll_f_colbound_jac11 ( xi, params )
    sll_real64  :: sll_f_colbound_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colbound_jac11 = L1*(1._f64+alpha1*a*cos(a*xi(1))*sin(b*xi(2)))
  end function sll_f_colbound_jac11

  !> jacobian matrix
  function sll_f_colbound_jac12 ( xi, params )
    sll_real64  :: sll_f_colbound_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha1, L1
    
    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colbound_jac12 = L1*alpha1*sin(a*xi(1))*cos(b*xi(2))*b
  end function sll_f_colbound_jac12

  !> jacobian matrix
  function sll_f_colbound_jac13 ( xi, params )
    sll_real64  :: sll_f_colbound_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    
    sll_f_colbound_jac13 =0._f64
  end function sll_f_colbound_jac13

  !> jacobian matrix
  function sll_f_colbound_jac21 ( xi, params )
    sll_real64  :: sll_f_colbound_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha2, L2
    
    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colbound_jac21 = L2*alpha2*a*cos(a*xi(1))*sin(b*xi(2))
  end function sll_f_colbound_jac21

  !> jacobian matrix
  function sll_f_colbound_jac22 ( xi, params )
    sll_real64  :: sll_f_colbound_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colbound_jac22 =L2*(1._f64+alpha2*sin(a*xi(1))*cos(b*xi(2))*b)
  end function sll_f_colbound_jac22

  !> jacobian matrix
  function sll_f_colbound_jac23 ( xi, params )
    sll_real64  :: sll_f_colbound_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colbound_jac23 = 0._f64
  end function sll_f_colbound_jac23

  !> jacobian matrix
  function sll_f_colbound_jac31 ( xi, params )
    sll_real64  :: sll_f_colbound_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_colbound_jac31 = 0._f64
  end function sll_f_colbound_jac31

  !> jacobian matrix
  function sll_f_colbound_jac32 ( xi, params )
    sll_real64  :: sll_f_colbound_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colbound_jac32 = 0._f64
  end function sll_f_colbound_jac32

  !> jacobian matrix
  function sll_f_colbound_jac33 ( xi, params )
    sll_real64  :: sll_f_colbound_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L3 = params(3)
    alpha3 = params(6)
    sll_f_colbound_jac33 = L3*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_colbound_jac33

    !> jacobian 
  function sll_f_colbound_jacobian ( xi, params )
    sll_real64  :: sll_f_colbound_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2, L3
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    L3 =     params(3)
    alpha1 = params(4)
    alpha2 = params(5)
    alpha3 = params(6)
    sll_f_colbound_jacobian = L3*L1*L2*(1._f64+b*alpha2*sin(a*xi(1))*cos(b*xi(2))+a*alpha1*cos(a*xi(1))*sin(b*xi(2)))*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_colbound_jacobian

  

  ! ***************************************************************************
  !
  ! <b> "Colella transformation"; </b>
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !  X1 = L1 * (\xi(1) + \alpha1 * sin(2*\pi * \xi(1)) * sin(2*\pi * \xi(2)) ) 
  !  X2 = L2 * (\xi(2) + \alpha2 * sin(2*\pi * \xi(1)) * sin(2*\pi * \xi(2)) ) 
  !  X3 = L3 * (\xi(3) + \alpha3 * sin(2*\pi * \xi(3)) ) 
  ! 
  ! The parameters are:
  !    - L1     = params(1)
  !    - L2     = params(2)
  !    - L3     = params(3)
  !    - alpha1 = params(4)
  !    - alpha2 = params(5)
  !    - alpha3 = params(6)
  !
  ! ***************************************************************************

  
  !> direct mapping
  function sll_f_colella_x1( xi, params )
    sll_real64 :: sll_f_colella_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1

    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colella_x1 = L1*(xi(1)+alpha1*sin(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2)))
  end function sll_f_colella_x1

  !> direct mapping
  function sll_f_colella_x2( xi, params )
    sll_real64 :: sll_f_colella_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2

    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colella_x2 = L2*(xi(2)+alpha2*sin(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2)))
  end function sll_f_colella_x2

  !> direct mapping
  function sll_f_colella_x3( xi, params )
    sll_real64 :: sll_f_colella_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L3     = params(3)
    alpha3 = params(6)
    sll_f_colella_x3 = L3*(xi(3)+alpha3*sin(sll_p_twopi*xi(3)))
  end function sll_f_colella_x3

  !> jacobian matrix
  function sll_f_colella_jac11 ( xi, params )
    sll_real64  :: sll_f_colella_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colella_jac11 = L1*(1._f64+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))*sll_p_twopi)
  end function sll_f_colella_jac11

  !> jacobian matrix
  function sll_f_colella_jac12 ( xi, params )
    sll_real64  :: sll_f_colella_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha1, L1
    
    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    sll_f_colella_jac12 = L1*alpha1*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_colella_jac12

  !> jacobian matrix
  function sll_f_colella_jac13 ( xi, params )
    sll_real64  :: sll_f_colella_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    
    sll_f_colella_jac13 =0._f64
  end function sll_f_colella_jac13

  !> jacobian matrix
  function sll_f_colella_jac21 ( xi, params )
    sll_real64  :: sll_f_colella_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: alpha2, L2
    
    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colella_jac21 = L2*alpha2*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_colella_jac21

  !> jacobian matrix
  function sll_f_colella_jac22 ( xi, params )
    sll_real64  :: sll_f_colella_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    sll_f_colella_jac22 =L2*(1._f64+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))*sll_p_twopi)
  end function sll_f_colella_jac22

  !> jacobian matrix
  function sll_f_colella_jac23 ( xi, params )
    sll_real64  :: sll_f_colella_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jac23 = 0._f64
  end function sll_f_colella_jac23

  !> jacobian matrix
  function sll_f_colella_jac31 ( xi, params )
    sll_real64  :: sll_f_colella_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_colella_jac31 = 0._f64
  end function sll_f_colella_jac31

  !> jacobian matrix
  function sll_f_colella_jac32 ( xi, params )
    sll_real64  :: sll_f_colella_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jac32 = 0._f64
  end function sll_f_colella_jac32

  !> jacobian matrix
  function sll_f_colella_jac33 ( xi, params )
    sll_real64  :: sll_f_colella_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L3 = params(3)
    alpha3 = params(6)
    sll_f_colella_jac33 = L3*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_colella_jac33

    !> jacobian 
  function sll_f_colella_jacobian ( xi, params )
    sll_real64  :: sll_f_colella_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2, L3
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    L3 =     params(3)
    alpha1 = params(4)
    alpha2 = params(5)
    alpha3 = params(6)
    sll_f_colella_jacobian = L3*L1*L2*(1._f64+sll_p_twopi*alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))+sll_p_twopi*alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2)))*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_colella_jacobian

!> jacobian matrix inverse
  function sll_f_colella_jacinv11 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L1, L2
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    alpha1 = params(4)
    alpha2 = params(5)
    sll_f_colella_jacinv11 = (1._f64/sll_p_twopi+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2)))/(L1*(1._f64/sll_p_twopi+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))))
  end function sll_f_colella_jacinv11

  !> jacobian matrix inverse
  function sll_f_colella_jacinv12 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    alpha1 = params(4)
    alpha2 = params(5)
    sll_f_colella_jacinv12 = - alpha1*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))/(L2*(1._f64/sll_p_twopi+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))))
  end function sll_f_colella_jacinv12

  !> jacobian matrix inverse
  function sll_f_colella_jacinv13 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jacinv13 = 0._f64
  end function sll_f_colella_jacinv13

  !> jacobian matrix inverse
  function sll_f_colella_jacinv21 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    alpha1 = params(4)
    alpha2 = params(5)
    sll_f_colella_jacinv21 = -alpha2*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))/(L1*(1._f64/sll_p_twopi+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))))
  end function sll_f_colella_jacinv21

  !> jacobian matrix inverse
  function sll_f_colella_jacinv22 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L1, L2
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    alpha1 = params(4)
    alpha2 = params(5)
    sll_f_colella_jacinv22 = (1._f64/sll_p_twopi+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2)))/(L2*(1._f64/sll_p_twopi+alpha2*sin(sll_p_twopi*xi(1))*cos(sll_p_twopi*xi(2))+alpha1*cos(sll_p_twopi*xi(1))*sin(sll_p_twopi*xi(2))))
  end function sll_f_colella_jacinv22

  !> jacobian matrix inverse
  function sll_f_colella_jacinv23 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jacinv23 = 0._f64
  end function sll_f_colella_jacinv23

  !> jacobian matrix inverse
  function sll_f_colella_jacinv31 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jacinv31 =0._f64
  end function sll_f_colella_jacinv31

  !> jacobian matrix inverse
  function sll_f_colella_jacinv32 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_colella_jacinv32 = 0._f64
  end function sll_f_colella_jacinv32

  !> jacobian matrix inverse
  function sll_f_colella_jacinv33 ( xi, params )
    sll_real64  :: sll_f_colella_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    !local variables
    sll_real64 :: L3
    SLL_ASSERT(size(params) >= 6)
    
    L3 = params(3)
    sll_f_colella_jacinv33 = 1._f64/L3
  end function sll_f_colella_jacinv33


  ! ***************************************************************************
  !
  ! <b> "Orthogonal transformation"; </b>
  !
  !  X1 = L1 * (\xi(1) + \alpha1 * sin(2*\pi * \xi(1)) ) 
  !  X2 = L2 * (\xi(2) + \alpha2 * sin(2*\pi * \xi(2)) ) 
  !  X3 = L3 * (\xi(3) + \alpha3 * sin(2*\pi * \xi(3)) ) 
  ! 
  ! The parameters are:
  !    - L1     = params(1)
  !    - L2     = params(2)
  !    - L3     = params(3)
  !    - alpha1 = params(4)
  !    - alpha2 = params(5)
  !    - alpha3 = params(6)
  !
  ! ***************************************************************************

  
  !> direct mapping
  function sll_f_orthogonal_x1( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1

    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    res = L1*(xi(1)+alpha1*sin(sll_p_twopi*xi(1)))
  end function sll_f_orthogonal_x1

  !> direct mapping
  function sll_f_orthogonal_x2( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2

    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    res = L2*(xi(2)+alpha2*sin(sll_p_twopi*xi(2)))
  end function sll_f_orthogonal_x2

  !> direct mapping
  function sll_f_orthogonal_x3( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
        
    SLL_ASSERT(size(params) >= 6)
    L3     = params(3)
    alpha3 = params(6)
    res= L3*(xi(3)+alpha3*sin(sll_p_twopi*xi(3)))
  end function sll_f_orthogonal_x3

  !> jacobian matrix
  function sll_f_orthogonal_jac11 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 6)
    L1     = params(1)
    alpha1 = params(4)
    res = L1*(1._f64+alpha1*cos(sll_p_twopi*xi(1))*sll_p_twopi)
  end function sll_f_orthogonal_jac11

  !> jacobian matrix
  function sll_f_orthogonal_jac12 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0.0_f64
  end function sll_f_orthogonal_jac12

  !> jacobian matrix
  function sll_f_orthogonal_jac13 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res =0._f64
  end function sll_f_orthogonal_jac13

  !> jacobian matrix
  function sll_f_orthogonal_jac21 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jac21

  !> jacobian matrix
  function sll_f_orthogonal_jac22 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L2     = params(2)
    alpha2 = params(5)
    res = L2*(1._f64+alpha2*cos(sll_p_twopi*xi(2))*sll_p_twopi)
  end function sll_f_orthogonal_jac22

  !> jacobian matrix
  function sll_f_orthogonal_jac23 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jac23

  !> jacobian matrix
  function sll_f_orthogonal_jac31 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jac31

  !> jacobian matrix
  function sll_f_orthogonal_jac32 ( xi, params )result( res )
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jac32

  !> jacobian matrix
  function sll_f_orthogonal_jac33 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L3 = params(3)
    alpha3 = params(6)
    res = L3*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_orthogonal_jac33

  !> jacobian 
  function sll_f_orthogonal_jacobian ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1, L2, L3
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    L2 =     params(2)
    L3 =     params(3)
    alpha1 = params(4)
    alpha2 = params(5)
    alpha3 = params(6)
    res = L3*L2*L1*(1._f64+alpha1*cos(sll_p_twopi*xi(1))*sll_p_twopi)*(1._f64+alpha2*cos(sll_p_twopi*xi(2))*sll_p_twopi)*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi)
  end function sll_f_orthogonal_jacobian
  
  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv11 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1
    
    SLL_ASSERT(size(params) >= 6)
    L1 =     params(1)
    alpha1 = params(4)
    
    res = 1._f64/(L1*(1._f64+alpha1*cos(sll_p_twopi*xi(1))*sll_p_twopi))
  end function sll_f_orthogonal_jacinv11

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv12 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0.0_f64
  end function sll_f_orthogonal_jacinv12

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv13 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output

    res =0._f64
  end function sll_f_orthogonal_jacinv13

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv21 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jacinv21

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv22 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha2
    
    SLL_ASSERT(size(params) >= 6)
    L2 =     params(2)
    alpha2 = params(5)
    res = 1._f64/(L2*(1._f64+alpha2*cos(sll_p_twopi*xi(2))*sll_p_twopi))
  end function sll_f_orthogonal_jacinv22

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv23 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jacinv23

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv31 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jacinv31

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv32 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_orthogonal_jacinv32

  !> jacobian matrix inverse
  function sll_f_orthogonal_jacinv33 ( xi, params ) result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha3
    SLL_ASSERT(size(params) >= 6)
    
    L3 = params(3)
    alpha3 = params(6)
    res = 1._f64/(L3*(1._f64+alpha3*cos(sll_p_twopi*xi(3))*sll_p_twopi))
  end function sll_f_orthogonal_jacinv33


  ! ***************************************************************************
  !
  ! <b> "Polynomial transformation"; </b>
  !
  !  X1 = L1 * (\xi(1) + \alpha11*\xi(1)**2) ) + \alpha12 *\xi(1)**3)
  !  X2 = L2 * (\xi(2) + \alpha21*\xi(2)**2) ) + \alpha22 *\xi(2)**3)
  !  X3 = L3 * (\xi(3) + \alpha31*\xi(3)**2) ) + \alpha32 *\xi(3)**3)
  ! 
  ! The parameters are:
  !    - L1     = params(1)
  !    - L2     = params(2)
  !    - L3     = params(3)
  !    - alpha11 = params(4)
  !    - alpha21 = params(5)
  !    - alpha31 = params(6)
  !    - alpha12 = params(7)
  !    - alpha22 = params(8)
  !    - alpha32 = params(9)
  !
  ! ***************************************************************************

  
  !> direct mapping
  function sll_f_polynomial_x1( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1, alpha2

    SLL_ASSERT(size(params) >= 9)
    L1     = params(1)
    alpha1 = params(4)
    alpha2 = params(7)
    res = L1 * (xi(1) + alpha1*xi(1)**2 + alpha2 *xi(1)**3)
  end function sll_f_polynomial_x1

  !> direct mapping
  function sll_f_polynomial_x2( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha1, alpha2

    SLL_ASSERT(size(params) >= 9)
    L2     = params(2)
    alpha1 = params(5)
    alpha2 = params(8)
  
    res = L2 * (xi(2) + alpha1*xi(2)**2 + alpha2 *xi(2)**3)
  end function sll_f_polynomial_x2

  !> direct mapping
  function sll_f_polynomial_x3( xi, params )result(res)
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha1, alpha2
        
    SLL_ASSERT(size(params) >= 9)
    L3     = params(3)
    alpha1 = params(6)
    alpha2 = params(9)
 
    res = L3 * (xi(3) + alpha1*xi(3)**2 + alpha2*xi(3)**3)
  end function sll_f_polynomial_x3
  
  !> jacobian matrix
  function sll_f_polynomial_jac11 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 9)
    L1     = params(1)
    alpha1 = params(4)
    alpha2 = params(7)

    res = L1*(1._f64+alpha1*2._f64*xi(1)+alpha2*3._f64*xi(1)**2)
  end function sll_f_polynomial_jac11

  !> jacobian matrix
  function sll_f_polynomial_jac12 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0.0_f64
  end function sll_f_polynomial_jac12

  !> jacobian matrix
  function sll_f_polynomial_jac13 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res =0._f64
  end function sll_f_polynomial_jac13

  !> jacobian matrix
  function sll_f_polynomial_jac21 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jac21

  !> jacobian matrix
  function sll_f_polynomial_jac22 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 9)
    L2     = params(2)
    alpha1 = params(5)
    alpha2 = params(8)

    res = L2*(1._f64+alpha1*2._f64*xi(2)+alpha2*3._f64*xi(2)**2)
  end function sll_f_polynomial_jac22

  !> jacobian matrix
  function sll_f_polynomial_jac23 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jac23

  !> jacobian matrix
  function sll_f_polynomial_jac31 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jac31

  !> jacobian matrix
  function sll_f_polynomial_jac32 ( xi, params )result( res )
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jac32

  !> jacobian matrix
  function sll_f_polynomial_jac33 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha1, alpha2
    
    SLL_ASSERT(size(params) >= 9)
    L3 = params(3)
    alpha1 = params(6)
    alpha2 = params(9)
 
    res = L3*(1._f64+alpha1*2._f64*xi(3)+alpha2*3._f64*xi(3)**2)
  end function sll_f_polynomial_jac33

  !> jacobian 
  function sll_f_polynomial_jacobian ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L1, L2, L3
    sll_real64 :: alpha11, alpha21, alpha31, alpha12, alpha22, alpha32
    
    SLL_ASSERT(size(params) >= 9)
    L1 =     params(1)
    L2 =     params(2)
    L3 =     params(3)
    alpha11 = params(4)
    alpha21 = params(5)
    alpha31 = params(6)
    alpha12 = params(7)
    alpha22 = params(8)
    alpha32 = params(9)
  

    res = L3*L2*L1*(1._f64+alpha11*2._f64*xi(1)+alpha12*3._f64*xi(1)**2)*(1._f64+alpha21*2._f64*xi(2)+alpha22*3._f64*xi(2)**2)*(1._f64+alpha31*2._f64*xi(3)+alpha32*3._f64*xi(3)**2)
  end function sll_f_polynomial_jacobian
  
  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv11 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
     sll_real64 :: L1
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 12)
    L1     = params(1)
    alpha1 = params(4)
    alpha2 = params(7)
    alpha3 = params(10)
    res = 1._f64/(L1*(alpha1+alpha2*xi(1)+alpha3*0.5_f64*xi(1)**2))
  end function sll_f_polynomial_jacinv11

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv12 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0.0_f64
  end function sll_f_polynomial_jacinv12

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv13 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output

    res =0._f64
  end function sll_f_polynomial_jacinv13

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv21 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jacinv21

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv22 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L2
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 12)
    L2     = params(2)
    alpha1 = params(5)
    alpha2 = params(8)
    alpha3 = params(11)
    res = 1._f64/(L2*(alpha1+alpha2*xi(2)+alpha3*0.5_f64*xi(2)**2))
  end function sll_f_polynomial_jacinv22

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv23 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jacinv23

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv31 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jacinv31

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv32 ( xi, params )result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    
    res = 0._f64
  end function sll_f_polynomial_jacinv32

  !> jacobian matrix inverse
  function sll_f_polynomial_jacinv33 ( xi, params ) result(res)
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res !Output
    !local variables
    sll_real64 :: L3
    sll_real64 :: alpha1, alpha2, alpha3
    
    SLL_ASSERT(size(params) >= 12)
    L3 = params(3)
    alpha1 = params(6)
    alpha2 = params(9)
    alpha3 = params(12)
    res = 1._f64/(L3*(alpha1+alpha2*xi(3)+alpha3*0.5_f64*xi(3)**2))
  end function sll_f_polynomial_jacinv33

  

  ! ***************************************************************************
  !
  ! Cylindrical coordinate transformation:
  !
  ! X1 = (Rmin + (Rmax-Rmin)*xi(1))*cos(2*pi*xi(2))
  ! X2 = (Rmin + (Rmax-Rmin)*xi(1))*sin(2*pi*xi(2))
  ! X3 =  L*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (Rmin, Rmax, L). Typically:
  !
  ! Rmin = 0.0
  ! Rmax = 2 *pi 
  ! L    = 2 *pi 
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_cylindrical_x1( xi, params )
    sll_real64 :: sll_f_cylindrical_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_x1 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_x1

  !> direct mapping
  function sll_f_cylindrical_x2( xi, params )
    sll_real64 :: sll_f_cylindrical_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_x2 = (r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_x2

  !> direct mapping
  function sll_f_cylindrical_x3( xi, params )
    sll_real64 :: sll_f_cylindrical_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L=params(3)
    sll_f_cylindrical_x3 = L*xi(3)
  end function sll_f_cylindrical_x3

  !> direct mapping
  function sll_f_cylindrical_xi1( x, params )
    sll_real64 :: sll_f_cylindrical_xi1
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_xi1 = (sqrt(x(1)**2+x(2)**2)-r1)/(r2-r1)
  end function sll_f_cylindrical_xi1

  !> direct mapping
  function sll_f_cylindrical_xi2( x, params )
    sll_real64 :: sll_f_cylindrical_xi2
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_xi2 = modulo(atan2(x(2), x(1))/sll_p_twopi, 1._f64) !atan2 gives a value between [-pi,pi], which we transform to [0,1] 
  end function sll_f_cylindrical_xi2

  !> direct mapping
  function sll_f_cylindrical_xi3( x, params )
    sll_real64 :: sll_f_cylindrical_xi3
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L = params(3)
    sll_f_cylindrical_xi3 = x(3)/L
  end function sll_f_cylindrical_xi3

  !> jacobian matrix
  function sll_f_cylindrical_jac11 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac11 = (r2-r1)*cos(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_jac11

  !> jacobian matrix
  function sll_f_cylindrical_jac12 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac12 = -(r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_cylindrical_jac12

  !> jacobian matrix
  function sll_f_cylindrical_jac13 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac13 =0._f64
  end function sll_f_cylindrical_jac13

  !> jacobian matrix
  function sll_f_cylindrical_jac21 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac21 = (r2-r1)*sin(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_jac21

  !> jacobian matrix
  function sll_f_cylindrical_jac22 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac22 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_cylindrical_jac22

  !> jacobian matrix
  function sll_f_cylindrical_jac23 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac23 = 0._f64
  end function sll_f_cylindrical_jac23

  !> jacobian matrix
  function sll_f_cylindrical_jac31 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac31 =0._f64
  end function sll_f_cylindrical_jac31

  !> jacobian matrix
  function sll_f_cylindrical_jac32 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jac32 = 0._f64
  end function sll_f_cylindrical_jac32

  !> jacobian matrix
  function sll_f_cylindrical_jac33 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L  = params(3)
    sll_f_cylindrical_jac33 = L
  end function sll_f_cylindrical_jac33

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv11 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv11 = cos(sll_p_twopi*xi(2))/(r2-r1)
  end function sll_f_cylindrical_jacinv11

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv12 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv12 = sin(sll_p_twopi*xi(2))/(r2-r1)
  end function sll_f_cylindrical_jacinv12

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv13 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv13 =0._f64
  end function sll_f_cylindrical_jacinv13

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv21 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv21 = -sin(sll_p_twopi*xi(2))/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))
  end function sll_f_cylindrical_jacinv21

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv22 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv22 = cos(sll_p_twopi*xi(2))/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))
  end function sll_f_cylindrical_jacinv22

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv23 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv23 = 0._f64
  end function sll_f_cylindrical_jacinv23

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv31 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv31 =0._f64
  end function sll_f_cylindrical_jacinv31

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv32 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_jacinv32 = 0._f64
  end function sll_f_cylindrical_jacinv32

  !> jacobian matrix inverse
  function sll_f_cylindrical_jacinv33 ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L  = params(3)
    sll_f_cylindrical_jacinv33 = 1._f64/L
  end function sll_f_cylindrical_jacinv33

  !> jacobian 
  function sll_f_cylindrical_jacobian ( xi, params )
    sll_real64  :: sll_f_cylindrical_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L  = params(3)
    sll_f_cylindrical_jacobian = sll_p_twopi*(r2-r1)*(r1 + (r2-r1)*xi(1))*L
  end function sll_f_cylindrical_jacobian

  ! ***************************************************************************
  !
  ! Cylindrical coordinate transformation:
  !
  ! X1 = (Rmin + (Rmax-Rmin)*sqrt(xi(1)))*cos(2*pi*xi(2))
  ! X2 = (Rmin + (Rmax-Rmin)*sqrt(xi(1)))*sin(2*pi*xi(2))
  ! X3 =  L*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (Rmin, Rmax, L). Typically:
  !
  ! Rmin = 0.0
  ! Rmax = 2 *pi 
  ! L    = 2 *pi 
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_cylindrical_sqrt_x1( xi, params )
    sll_real64 :: sll_f_cylindrical_sqrt_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_x1 = (r1 + (r2-r1)*sqrt(xi(1)))*cos(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_sqrt_x1

  !> direct mapping
  function sll_f_cylindrical_sqrt_x2( xi, params )
    sll_real64 :: sll_f_cylindrical_sqrt_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_x2 = (r1 + (r2-r1)*sqrt(xi(1)))*sin(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_sqrt_x2

  !> direct mapping
  function sll_f_cylindrical_sqrt_x3( xi, params )
    sll_real64 :: sll_f_cylindrical_sqrt_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L=params(3)
    sll_f_cylindrical_sqrt_x3 = L*xi(3)
  end function sll_f_cylindrical_sqrt_x3

  !> direct mapping
  function sll_f_cylindrical_sqrt_xi1( x, params )
    sll_real64 :: sll_f_cylindrical_sqrt_xi1
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_xi1 = ((sqrt(x(1)**2+x(2)**2)-r1)/(r2-r1))**2
  end function sll_f_cylindrical_sqrt_xi1

  !> direct mapping
  function sll_f_cylindrical_sqrt_xi2( x, params )
    sll_real64 :: sll_f_cylindrical_sqrt_xi2
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_xi2 = modulo(atan2(x(2), x(1))/sll_p_twopi, 1._f64)!atan2 gives a value between [-pi,pi], which we transform to [0,1] 
  end function sll_f_cylindrical_sqrt_xi2

  !> direct mapping
  function sll_f_cylindrical_sqrt_xi3( x, params )
    sll_real64 :: sll_f_cylindrical_sqrt_xi3
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L = params(3)
    sll_f_cylindrical_sqrt_xi3 = x(3)/L
  end function sll_f_cylindrical_sqrt_xi3

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac11 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac11 = (r2-r1)/( 2._f64*sqrt(xi(1)) )*cos(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_sqrt_jac11

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac12 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac12 = -(r1 + (r2-r1)*sqrt(xi(1)))*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_cylindrical_sqrt_jac12

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac13 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac13 =0._f64
  end function sll_f_cylindrical_sqrt_jac13

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac21 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac21 = (r2-r1)/( 2._f64*sqrt(xi(1)) )*sin(sll_p_twopi*xi(2))
  end function sll_f_cylindrical_sqrt_jac21

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac22 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac22 = (r1 + (r2-r1)*sqrt(xi(1)))*cos(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_cylindrical_sqrt_jac22

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac23 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac23 = 0._f64
  end function sll_f_cylindrical_sqrt_jac23

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac31 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac31 =0._f64
  end function sll_f_cylindrical_sqrt_jac31

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac32 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    sll_f_cylindrical_sqrt_jac32 = 0._f64
  end function sll_f_cylindrical_sqrt_jac32

  !> jacobian matrix
  function sll_f_cylindrical_sqrt_jac33 ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L  = params(3)
    sll_f_cylindrical_sqrt_jac33 = L
  end function sll_f_cylindrical_sqrt_jac33

  !> jacobian 
  function sll_f_cylindrical_sqrt_jacobian ( xi, params )
    sll_real64  :: sll_f_cylindrical_sqrt_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    L  = params(3)
    sll_f_cylindrical_sqrt_jacobian = sll_p_pi*(r2-r1)*(r1/sqrt(xi(1)) + (r2-r1))*L
  end function sll_f_cylindrical_sqrt_jacobian


  

  ! ***************************************************************************
  !
  ! Elliptical coordinate transformation:
  !
  ! X1 = R*cosh(xi(1)+rmin)*cos(2*pi*xi(2))
  ! X2 = R*sinh(xi(1)+rmin)*sin(2*pi*xi(2))
  ! X3 = L*xi(3)
  ! Where xi are defined in the interval [0,1]. The 'params' array
  ! contains the information (rmin, R, L). Typically:
  !
  ! rmin = 0.0
  ! R    = 2 *pi 
  ! L    = 2 *pi 
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_elliptical_x1( xi, params )
    sll_real64 :: sll_f_elliptical_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
   
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_x1 = R*cosh(xi(1)+rmin)*cos(sll_p_twopi*xi(2))
  end function sll_f_elliptical_x1

  !> direct mapping
  function sll_f_elliptical_x2( xi, params )
    sll_real64 :: sll_f_elliptical_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
   
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_x2 = R*sinh(xi(1)+rmin)*sin(sll_p_twopi*xi(2))
  end function sll_f_elliptical_x2

  !> direct mapping
  function sll_f_elliptical_x3( xi, params )
    sll_real64 :: sll_f_elliptical_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    L=params(3)
    sll_f_elliptical_x3 = L*xi(3)
  end function sll_f_elliptical_x3


  !> jacobian matrix
  function sll_f_elliptical_jac11 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
   
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_jac11 = R*sinh(xi(1)+rmin)*cos(sll_p_twopi*xi(2))
  end function sll_f_elliptical_jac11

  !> jacobian matrix
  function sll_f_elliptical_jac12 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
  
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_jac12 = - R*cosh(xi(1)+rmin)*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_elliptical_jac12

  !> jacobian matrix
  function sll_f_elliptical_jac13 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
 
    sll_f_elliptical_jac13 = 0._f64
  end function sll_f_elliptical_jac13

  !> jacobian matrix
  function sll_f_elliptical_jac21 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
   
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_jac21 = R*cosh(xi(1)+rmin)*sin(sll_p_twopi*xi(2))
  end function sll_f_elliptical_jac21

  !> jacobian matrix
  function sll_f_elliptical_jac22 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
   
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    sll_f_elliptical_jac22 = R*sinh(xi(1)+rmin)*cos(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_elliptical_jac22

  !> jacobian matrix
  function sll_f_elliptical_jac23 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
 
    sll_f_elliptical_jac23 = 0._f64
  end function sll_f_elliptical_jac23

  !> jacobian matrix
  function sll_f_elliptical_jac31 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_elliptical_jac31 = 0._f64
  end function sll_f_elliptical_jac31

  !> jacobian matrix
  function sll_f_elliptical_jac32 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
 
    sll_f_elliptical_jac32 = 0._f64
  end function sll_f_elliptical_jac32

  !> jacobian matrix
  function sll_f_elliptical_jac33 ( xi, params )
    sll_real64  :: sll_f_elliptical_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
  
    L  = params(3)
    sll_f_elliptical_jac33 = L
  end function sll_f_elliptical_jac33

  

  !> jacobian 
  function sll_f_elliptical_jacobian ( xi, params )
    sll_real64  :: sll_f_elliptical_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: rmin, R
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 3)
    rmin = params(1)
    R = params(2)
    L  = params(3)
    sll_f_elliptical_jacobian = L*R**2*sll_p_twopi*(sinh(xi(1)+rmin )**2+sin(sll_p_twopi*xi(2))**2 )
  end function sll_f_elliptical_jacobian


    ! ***************************************************************************
  !
  ! Dshaped_Singular coordinate transformation:
  !
  ! X1 = 1/eps*(1 - sqrt(1+eps*(eps+2*xi(1))*cos(2*pi*xi(2)) )
  ! X2 = y_0 + (e*(1/sqrt(1-eps^2/4))*xi(1)*sin(2*pi*xi(2)))/(1+eps*X1(xi))
  ! X3 =  L*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information ( eps, e, Lz, y_0 ). Typically:
  !
  ! eps = 0.3
  ! e = 1.4
  ! L = 2 *pi
  ! y_0 = 0
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_Dshaped_singular_x1( xi, params )
    sll_real64 :: sll_f_Dshaped_singular_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e

    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    sll_f_Dshaped_singular_x1 = (1-sqrt(1+eps*(eps+2._f64*xi(1)*cos( sll_p_twopi*xi(2) ) ) ) )/eps
  end function sll_f_Dshaped_singular_x1

  !> direct mapping
  function sll_f_Dshaped_singular_x2( xi, params )
    sll_real64 :: sll_f_Dshaped_singular_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e
    sll_real64 :: y_0

    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    y_0 = params(4)
    sll_f_Dshaped_singular_x2 = y_0+ e*(1._f64/sqrt(1._f64-eps**2/4._f64) )*xi(1)*sin(sll_p_twopi*xi(2)) / (1+eps*sll_f_Dshaped_singular_x1( xi, params ) )
  end function sll_f_Dshaped_singular_x2

  !> direct mapping
  function sll_f_Dshaped_singular_x3( xi, params )
    sll_real64 :: sll_f_Dshaped_singular_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 4)
    L=params(3)
    sll_f_Dshaped_singular_x3 = L*xi(3)
  end function sll_f_Dshaped_singular_x3

 
  !> jacobian matrix
  function sll_f_Dshaped_singular_jac11 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps

    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    sll_f_Dshaped_singular_jac11 = -cos(sll_p_twopi*xi(2))/(1-eps*sll_f_Dshaped_singular_x1( xi, params ) )
  end function sll_f_Dshaped_singular_jac11

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac12 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    sll_f_Dshaped_singular_jac12 = sll_p_twopi*sin(sll_p_twopi*xi(2))/(1-eps*sll_f_Dshaped_singular_x1( xi, params ) )
  end function sll_f_Dshaped_singular_jac12

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac13 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
   
    sll_f_Dshaped_singular_jac13 =0._f64
  end function sll_f_Dshaped_singular_jac13

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac21 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e
    sll_real64 :: x1
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    x1 = sll_f_Dshaped_singular_x1( xi, params )
    sll_f_Dshaped_singular_jac21 = (e*(1/sqrt(1-eps**2/4._f64)) )/(1+eps*x1  ) * ( sin(sll_p_twopi*xi(2)) + (eps*xi(1)*sin(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(2)) ) / (1-eps**2 * x1**2 ) )
  end function sll_f_Dshaped_singular_jac21

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac22 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e
    sll_real64 :: x1
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    x1 = sll_f_Dshaped_singular_x1( xi, params )
    sll_f_Dshaped_singular_jac22 = (sll_p_twopi*e*(1/sqrt(1-eps**2/4._f64)) )/(1+eps*x1  ) * ( xi(1)*cos(sll_p_twopi*xi(2)) - (eps*xi(1)**2 * sin(sll_p_twopi*xi(2))**2 ) / (1-eps**2 * x1**2 ) )
  end function sll_f_Dshaped_singular_jac22

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac23 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_jac23 = 0._f64
  end function sll_f_Dshaped_singular_jac23

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac31 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_jac31 =0._f64
  end function sll_f_Dshaped_singular_jac31

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac32 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_jac32 = 0._f64
  end function sll_f_Dshaped_singular_jac32

  !> jacobian matrix
  function sll_f_Dshaped_singular_jac33 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 4)
    L  = params(3)
    sll_f_Dshaped_singular_jac33 = L
  end function sll_f_Dshaped_singular_jac33

 

  !> jacobian 
  function sll_f_Dshaped_singular_jacobian ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    L  = params(3)
    sll_f_Dshaped_singular_jacobian = L*sll_p_twopi*xi(1)/(eps*sll_f_Dshaped_singular_x1( xi, params ) -1._f64)
  end function sll_f_Dshaped_singular_jacobian

   !> jacobian matrix
  function sll_f_Dshaped_singular_pseudoinv11 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps

    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    sll_f_Dshaped_singular_pseudoinv11 = eps*sll_f_Dshaped_singular_x1( xi, params ) -1._f64
  end function sll_f_Dshaped_singular_pseudoinv11

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv12 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv12 = 0._f64
  end function sll_f_Dshaped_singular_pseudoinv12

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv13 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv13 =0._f64
  end function sll_f_Dshaped_singular_pseudoinv13

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv21 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e
    sll_real64 :: x1
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    x1 = sll_f_Dshaped_singular_x1( xi, params )
    sll_f_Dshaped_singular_pseudoinv21 = eps*xi(1) * sin(sll_p_twopi*xi(2)) /(1._f64 + eps*sll_f_Dshaped_singular_x1( xi, params ) )
  end function sll_f_Dshaped_singular_pseudoinv21

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv22 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: eps
    sll_real64 :: e
    sll_real64 :: x1
    SLL_ASSERT(size(params) >= 4)
    eps = params(1)
    e = params(2)
    x1 = sll_f_Dshaped_singular_x1( xi, params )
    sll_f_Dshaped_singular_pseudoinv22 = (1._f64 + eps*sll_f_Dshaped_singular_x1( xi, params ) )/ (e*(1/sqrt(1-eps**2/4._f64)) )
  end function sll_f_Dshaped_singular_pseudoinv22

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv23 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv23 = 0._f64
  end function sll_f_Dshaped_singular_pseudoinv23

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv31 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv31 =0._f64
  end function sll_f_Dshaped_singular_pseudoinv31

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv32 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv32 = 0._f64
  end function sll_f_Dshaped_singular_pseudoinv32

  !> pseudoinvobian matrix
  function sll_f_Dshaped_singular_pseudoinv33 ( xi, params )
    sll_real64  :: sll_f_Dshaped_singular_pseudoinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_Dshaped_singular_pseudoinv33 = 1._f64
  end function sll_f_Dshaped_singular_pseudoinv33


    ! ***************************************************************************
  !
  ! Rotated rotated singular coordinate transformation:
  !
  ! X1 = x_0 + (1-kappa)*xi(1))*cos(2*pi*xi(2)) - delta xi(1)^2
  ! X2 = y_0 + (1+kappa)*xi(1))*sin(2*pi*xi(2))
  ! X3 =  Lz*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (kappa, delta, Lz,  x_0, y_0). Typically:
  !
  ! kappa = 0.3
  ! delta = 0.2
  ! Lz    = 2 *pi
  ! x_0 = 0
  ! y_0 = 0
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_rotated_singular_x1( xi, params )
    sll_real64 :: sll_f_rotated_singular_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta
    sll_real64 :: x_0

    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    x_0 = params(4)
    sll_f_rotated_singular_x1 = x_0 + (1-kappa)*xi(1)*cos(sll_p_twopi*xi(2))-delta*xi(1)**2
  end function sll_f_rotated_singular_x1

  !> direct mapping
  function sll_f_rotated_singular_x2( xi, params )
    sll_real64 :: sll_f_rotated_singular_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa
    sll_real64 :: y_0

    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    y_0 = params(5)
    sll_f_rotated_singular_x2 = y_0+(1+kappa)*xi(1)*sin(sll_p_twopi*xi(2))
  end function sll_f_rotated_singular_x2

  !> direct mapping
  function sll_f_rotated_singular_x3( xi, params )
    sll_real64 :: sll_f_rotated_singular_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 5)
    L=params(3)
    sll_f_rotated_singular_x3 = L*xi(3)
  end function sll_f_rotated_singular_x3

  !> jacobian matrix
  function sll_f_rotated_singular_jac11 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta
 
    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    sll_f_rotated_singular_jac11 = (1-kappa)*xi(1)*cos(sll_p_twopi*xi(2))-2._f64*delta*xi(1)
  end function sll_f_rotated_singular_jac11

  !> jacobian matrix
  function sll_f_rotated_singular_jac12 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta

    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    sll_f_rotated_singular_jac12 =  -(1-kappa)*xi(1)*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_rotated_singular_jac12

  !> jacobian matrix
  function sll_f_rotated_singular_jac13 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
   
    sll_f_rotated_singular_jac13 =0._f64
  end function sll_f_rotated_singular_jac13

  !> jacobian matrix
  function sll_f_rotated_singular_jac21 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa
  
    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    sll_f_rotated_singular_jac21 = (1+kappa)*sin(sll_p_twopi*xi(2))
  end function sll_f_rotated_singular_jac21

  !> jacobian matrix
  function sll_f_rotated_singular_jac22 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa
  
    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    sll_f_rotated_singular_jac22 = (1+kappa)*xi(1)*cos(sll_p_twopi*xi(2)) *sll_p_twopi
  end function sll_f_rotated_singular_jac22

  !> jacobian matrix
  function sll_f_rotated_singular_jac23 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
   
    sll_f_rotated_singular_jac23 = 0._f64
  end function sll_f_rotated_singular_jac23

  !> jacobian matrix
  function sll_f_rotated_singular_jac31 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
  
    sll_f_rotated_singular_jac31 =0._f64
  end function sll_f_rotated_singular_jac31

  !> jacobian matrix
  function sll_f_rotated_singular_jac32 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_rotated_singular_jac32 = 0._f64
  end function sll_f_rotated_singular_jac32

  !> jacobian matrix
  function sll_f_rotated_singular_jac33 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: L
    SLL_ASSERT(size(params) >= 5)
    
    L  = params(3)
    sll_f_rotated_singular_jac33 = L
  end function sll_f_rotated_singular_jac33

 
  !> jacobian 
  function sll_f_rotated_singular_jacobian ( xi, params )
    sll_real64  :: sll_f_rotated_singular_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta
    sll_real64 :: L

    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    L  = params(3)
    sll_f_rotated_singular_jacobian = sll_p_twopi*xi(1)*(1+kappa) *((1-kappa)-2*delta*xi(1)*cos(sll_p_twopi*xi(2)))*L
  end function sll_f_rotated_singular_jacobian

  !> jacobian matrix
  function sll_f_rotated_singular_pseudoinv11 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta
 
    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    sll_f_rotated_singular_pseudoinv11 = 1._f64/(1-kappa-2._f64*delta*xi(1)*cos(sll_p_twopi*xi(2)) )
  end function sll_f_rotated_singular_pseudoinv11

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv12 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa, delta

    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    delta = params(2)
    sll_f_rotated_singular_pseudoinv12 =  2._f64*delta*xi(1)*sin(sll_p_twopi*xi(2)) /((1-kappa-2._f64*delta*xi(1)*cos(sll_p_twopi*xi(2)) )*(1+kappa) )
  end function sll_f_rotated_singular_pseudoinv12

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv13 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
   
    sll_f_rotated_singular_pseudoinv13 =0._f64
  end function sll_f_rotated_singular_pseudoinv13

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv21 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_rotated_singular_pseudoinv21 = 0._f64
  end function sll_f_rotated_singular_pseudoinv21

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv22 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: kappa
  
    SLL_ASSERT(size(params) >= 5)
    kappa = params(1)
    sll_f_rotated_singular_pseudoinv22 = 1._f64/(1+kappa)
  end function sll_f_rotated_singular_pseudoinv22

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv23 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
   
    sll_f_rotated_singular_pseudoinv23 = 0._f64
  end function sll_f_rotated_singular_pseudoinv23

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv31 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
  
    sll_f_rotated_singular_pseudoinv31 =0._f64
  end function sll_f_rotated_singular_pseudoinv31

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv32 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_rotated_singular_pseudoinv32 = 0._f64
  end function sll_f_rotated_singular_pseudoinv32

  !> pseudoinvobian matrix
  function sll_f_rotated_singular_pseudoinv33 ( xi, params )
    sll_real64  :: sll_f_rotated_singular_pseudoinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    sll_f_rotated_singular_pseudoinv33 = 1._f64
  end function sll_f_rotated_singular_pseudoinv33


  ! ***************************************************************************
  !
  ! Rotation coordinate transformation:
  !
  ! X1 = L1*(xi(1)*cos(thxi)-xi(2)*sin(thxi))
  ! X2 = L2*(xi(1)*sin(thxi)+xi(2)*cos(thxi))
  ! X3 = L3*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (L1,L2,L3,thxi).
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_rotation_x1( xi, params )
    sll_real64 :: sll_f_rotation_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_x1 = params(1)*modulo(xi(1)*cos(params(4))-xi(2)*sin(params(4)),1._f64)
  end function sll_f_rotation_x1

  !> direct mapping
  function sll_f_rotation_x2( xi, params )
    sll_real64 :: sll_f_rotation_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_x2 = params(2)*modulo(xi(1)*sin(params(4))+xi(2)*cos(params(4)),1._f64)
  end function sll_f_rotation_x2

  !> direct mapping
  function sll_f_rotation_x3( xi, params )
    sll_real64 :: sll_f_rotation_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_x3 = params(3)*modulo(xi(3),1._f64)
  end function sll_f_rotation_x3

  !> jacobian matrix
  function sll_f_rotation_jac11 ( xi, params )
    sll_real64  :: sll_f_rotation_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_rotation_jac11 =params(1)*cos(params(4))
  end function sll_f_rotation_jac11

  !> jacobian matrix
  function sll_f_rotation_jac12 ( xi, params )
    sll_real64  :: sll_f_rotation_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_rotation_jac12 = params(1)*sin(params(4))
  end function sll_f_rotation_jac12

  !> jacobian matrix
  function sll_f_rotation_jac13 ( xi, params )
    sll_real64  :: sll_f_rotation_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac13 =0._f64
  end function sll_f_rotation_jac13

  !> jacobian matrix
  function sll_f_rotation_jac21 ( xi, params )
    sll_real64  :: sll_f_rotation_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac21 = params(2)*sin(params(4))
  end function sll_f_rotation_jac21

  !> jacobian matrix
  function sll_f_rotation_jac22 ( xi, params )
    sll_real64  :: sll_f_rotation_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac22 = params(2)*cos(params(4))
  end function sll_f_rotation_jac22

  !> jacobian matrix
  function sll_f_rotation_jac23 ( xi, params )
    sll_real64  :: sll_f_rotation_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac23 = 0._f64
  end function sll_f_rotation_jac23

  !> jacobian matrix
  function sll_f_rotation_jac31 ( xi, params )
    sll_real64  :: sll_f_rotation_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac31 =0._f64
  end function sll_f_rotation_jac31

  !> jacobian matrix
  function sll_f_rotation_jac32 ( xi, params )
    sll_real64  :: sll_f_rotation_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac32 = 0._f64
  end function sll_f_rotation_jac32

  !> jacobian matrix
  function sll_f_rotation_jac33 ( xi, params )
    sll_real64  :: sll_f_rotation_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jac33 = params(3)
  end function sll_f_rotation_jac33

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv11 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv11 = cos(params(4))/params(1)
  end function sll_f_rotation_jacinv11

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv12 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv12 = sin(params(4))/params(1)
  end function sll_f_rotation_jacinv12

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv13 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv13 =0._f64
  end function sll_f_rotation_jacinv13

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv21 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv21 = -sin(params(4))/params(2)
  end function sll_f_rotation_jacinv21

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv22 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv22 = cos(params(4))/params(2)
  end function sll_f_rotation_jacinv22

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv23 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv23 = 0._f64
  end function sll_f_rotation_jacinv23

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv31 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv31 =0._f64
  end function sll_f_rotation_jacinv31

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv32 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv32 = 0._f64
  end function sll_f_rotation_jacinv32

  !> jacobian matrix inverse
  function sll_f_rotation_jacinv33 ( xi, params )
    sll_real64  :: sll_f_rotation_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacinv33 = 1._f64/params(3)
  end function sll_f_rotation_jacinv33

  !> jacobian 
  function sll_f_rotation_jacobian ( xi, params )
    sll_real64  :: sll_f_rotation_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params

    SLL_ASSERT(size(params) >= 4)
    sll_f_rotation_jacobian = params(1)*params(2)*params(3)
  end function sll_f_rotation_jacobian

  
  ! ***************************************************************************
  !
  ! Toroidal_Cylinder coordinate transformation:
  !
  ! X1 = (R0 + (Rmin + (Rmax-Rmin)*xi(1))*cos(2*pi*xi(3))) * cos(xi(2))
  ! X2 = (R0 + (Rmin + (Rmax-Rmin)*xi(1))*cos(2*pi*xi(3))) * sin(xi(2))
  ! X3 =       (Rmin + (Rmax-Rmin)*xi(1))*sin(2*pi*xi(3))
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (Rmin, Rmax, R0). Typically:
  ! Rmin <= Rmax < R0
  ! R0 = 2.0  
  ! R1 = 0.0
  ! R2 = 1.0
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_toroidal_cylinder_x1( xi, params )
    sll_real64 :: sll_f_toroidal_cylinder_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_x1 = (r0+(r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(3))) * cos(sll_p_twopi*xi(2))
  end function sll_f_toroidal_cylinder_x1

  !> direct mapping
  function sll_f_toroidal_cylinder_x2( xi, params )
    sll_real64 :: sll_f_toroidal_cylinder_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_x2 = (r0+(r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(3)))*sin(sll_p_twopi*xi(2)) 
  end function sll_f_toroidal_cylinder_x2

  !> direct mapping
  function sll_f_toroidal_cylinder_x3( xi, params )
    sll_real64 :: sll_f_toroidal_cylinder_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_x3 = (r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(3))
  end function sll_f_toroidal_cylinder_x3

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac11 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac11 = (r2-r1)*cos(sll_p_twopi*xi(3))*cos(sll_p_twopi*xi(2))
  end function sll_f_toroidal_cylinder_jac11

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac12 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac12 = -(r0+(r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(3)))*sin(sll_p_twopi*xi(2)) *sll_p_twopi
  end function sll_f_toroidal_cylinder_jac12

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac13 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac13 = -(r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(3))*cos(sll_p_twopi*xi(2))*sll_p_twopi 
  end function sll_f_toroidal_cylinder_jac13

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac21 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac21 = (r2-r1)*cos(sll_p_twopi*xi(3))*sin(sll_p_twopi*xi(2))
  end function sll_f_toroidal_cylinder_jac21

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac22 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac22 = (r0+(r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(3)))*cos(sll_p_twopi*xi(2))*sll_p_twopi 
  end function sll_f_toroidal_cylinder_jac22

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac23 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac23 = -(r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(3))*sin(sll_p_twopi*xi(2))*sll_p_twopi 
  end function sll_f_toroidal_cylinder_jac23

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac31 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac31 = (r2-r1)*sin(sll_p_twopi*xi(3))
  end function sll_f_toroidal_cylinder_jac31

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac32 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac32 = 0._f64
  end function sll_f_toroidal_cylinder_jac32

  !> jacobian matrix
  function sll_f_toroidal_cylinder_jac33 ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jac33 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(3))*sll_p_twopi 
  end function sll_f_toroidal_cylinder_jac33

    !> jacobian 
  function sll_f_toroidal_cylinder_jacobian ( xi, params )
    sll_real64  :: sll_f_toroidal_cylinder_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    sll_real64 :: r0

    SLL_ASSERT(size(params) >= 3)
    r1 = params(1)
    r2 = params(2)
    r0 = params(3)
    sll_f_toroidal_cylinder_jacobian = r0*(r1 + (r2-r1)*xi(1))*sll_p_twopi+(r1 + (r2-r1)*xi(1))**2*cos(sll_p_twopi*xi(3))*sll_p_twopi**2 
    !r*R0+xi(1)^2 cos(xi(3))
  end function sll_f_toroidal_cylinder_jacobian
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv11 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv11
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv11 = cos(sll_p_twopi*xi(2))/(r2-r1)
!!$  end function sll_f_toroidal_cylinder_jacinv11
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv12 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv12
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv12 = sin(sll_p_twopi*xi(2))/(r2-r1)
!!$  end function sll_f_toroidal_cylinder_jacinv12
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv13 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv13
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv13 =0._f64
!!$  end function sll_f_toroidal_cylinder_jacinv13
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv21 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv21
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv21 = -sin(sll_p_twopi*xi(2))/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))
!!$  end function sll_f_toroidal_cylinder_jacinv21
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv22 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv22
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv22 = cos(sll_p_twopi*xi(2))/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))
!!$  end function sll_f_toroidal_cylinder_jacinv22
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv23 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv23
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv23 = 0._f64
!!$  end function sll_f_toroidal_cylinder_jacinv23
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv31 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv31
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv31 =0._f64
!!$  end function sll_f_toroidal_cylinder_jacinv31
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv32 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv32
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv32 = 0._f64
!!$  end function sll_f_toroidal_cylinder_jacinv32
!!$
!!$  !> jacobian matrix inverse
!!$  function sll_f_toroidal_cylinder_jacinv33 ( xi, params )
!!$    sll_real64  :: sll_f_toroidal_cylinder_jacinv33
!!$    sll_real64, intent(in)   :: xi(3)
!!$    sll_real64, dimension(:), intent(in) :: params
!!$    sll_real64 :: r1
!!$    sll_real64 :: r2
!!$    SLL_ASSERT(size(params) >= 2)
!!$    r1 = params(1)
!!$    r2 = params(2)
!!$    sll_f_toroidal_cylinder_jacinv33 = 1._f64
!!$  end function sll_f_toroidal_cylinder_jacinv33





  ! ***************************************************************************
  !
  ! Spherical coordinate transformation:
  !
  ! X1 = (Rmin + (Rmax-Rmin)*xi(1))*sin(2*pi*xi(2))*cos(2*pi*xi(3))
  ! X2 = (Rmin + (Rmax-Rmin)*xi(1))*sin(2*pi*xi(2))*sin(2*pi*xi(3))
  ! X3 = (Rmin + (Rmax-Rmin)*xi(1))*cos(2*pi*xi(2))
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (Rmin, Rmax). Typically:
  !
  ! Rmin = 0.0
  ! Rmax = 1.0
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_spherical_x1( xi, params )
    sll_real64 :: sll_f_spherical_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_x1 = (r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))
  end function sll_f_spherical_x1

  !> direct mapping
  function sll_f_spherical_x2( xi, params )
    sll_real64 :: sll_f_spherical_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_x2= (r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))
  end function sll_f_spherical_x2

  !> direct mapping
  function sll_f_spherical_x3( xi, params )
    sll_real64 :: sll_f_spherical_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_x3 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(2))
  end function sll_f_spherical_x3

  !> jacobian matrix
  function sll_f_spherical_jac11 ( xi, params )
    sll_real64  :: sll_f_spherical_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac11 = (r2-r1)*sin(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))
  end function sll_f_spherical_jac11

  !> jacobian matrix
  function sll_f_spherical_jac12 ( xi, params )
    sll_real64  :: sll_f_spherical_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac12 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))*sll_p_twopi
  end function sll_f_spherical_jac12

  !> jacobian matrix
  function sll_f_spherical_jac13 ( xi, params )
    sll_real64  :: sll_f_spherical_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac13 = -(r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))*sll_p_twopi
  end function sll_f_spherical_jac13

  !> jacobian matrix
  function sll_f_spherical_jac21 ( xi, params )
    sll_real64  :: sll_f_spherical_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac21 = (r2-r1)*sin(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))
  end function sll_f_spherical_jac21

  !> jacobian matrix
  function sll_f_spherical_jac22 ( xi, params )
    sll_real64  :: sll_f_spherical_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac22 = (r1 + (r2-r1)*xi(1))*cos(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))*sll_p_twopi
  end function sll_f_spherical_jac22

  !> jacobian matrix
  function sll_f_spherical_jac23 ( xi, params )
    sll_real64  :: sll_f_spherical_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac23 = (r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))*sll_p_twopi
  end function sll_f_spherical_jac23

  !> jacobian matrix
  function sll_f_spherical_jac31 ( xi, params )
    sll_real64  :: sll_f_spherical_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac31 = (r2-r1)*cos(sll_p_twopi*xi(2))
  end function sll_f_spherical_jac31

  !> jacobian matrix
  function sll_f_spherical_jac32 ( xi, params )
    sll_real64  :: sll_f_spherical_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac32 = -(r1 + (r2-r1)*xi(1))*sin(sll_p_twopi*xi(2))*sll_p_twopi
  end function sll_f_spherical_jac32

  !> jacobian matrix
  function sll_f_spherical_jac33 ( xi, params )
    sll_real64  :: sll_f_spherical_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jac33 = 0._f64
  end function sll_f_spherical_jac33

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv11 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv11 = sin(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))/(r2-r1) 
  end function sll_f_spherical_jacinv11

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv12 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv12 = sin(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))/(r2-r1) 
  end function sll_f_spherical_jacinv12

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv13 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv13 = cos(sll_p_twopi*xi(2))/(r2-r1) 
  end function sll_f_spherical_jacinv13

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv21 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv21 = 1._f64/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))*cos(sll_p_twopi*xi(2))*cos(sll_p_twopi*xi(3))
  end function sll_f_spherical_jacinv21

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv22 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv22 = 1._f64/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))*cos(sll_p_twopi*xi(2))*sin(sll_p_twopi*xi(3))
  end function sll_f_spherical_jacinv22

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv23 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv23 = -1._f64/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))*sin(sll_p_twopi*xi(2))
  end function sll_f_spherical_jacinv23

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv31 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv31 = -1._f64/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))*sin(sll_p_twopi*xi(3))/sin(sll_p_twopi*xi(2))
  end function sll_f_spherical_jacinv31

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv32 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv32 = 1._f64/(sll_p_twopi*(r1 + (r2-r1)*xi(1)))*cos(sll_p_twopi*xi(3))/sin(sll_p_twopi*xi(2))
  end function sll_f_spherical_jacinv32

  !> jacobian matrix inverse
  function sll_f_spherical_jacinv33 ( xi, params )
    sll_real64  :: sll_f_spherical_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacinv33 = 0._f64
  end function sll_f_spherical_jacinv33

  !> jacobian
  function sll_f_spherical_jacobian ( xi, params )
    sll_real64  :: sll_f_spherical_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: r1
    sll_real64 :: r2
    SLL_ASSERT(size(params) >= 2)
    r1 = params(1)
    r2 = params(2)
    sll_f_spherical_jacobian =(r1 + (r2-r1)*xi(1))**2._f64 *sin(sll_p_twopi*xi(2))*sll_p_twopi**2._f64*(r2-r1)
  end function sll_f_spherical_jacobian



  ! ***************************************************************************
  !
  !Identity mapping:
  !
  ! X1 = xi(1)
  ! X2 = xi(2)
  ! X3 = xi(3)
  ! Where xi(3) are defined in the interval [0,1]. 
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_identity_x1( xi, params )
    sll_real64 :: sll_f_identity_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_x1 = xi(1)
  end function sll_f_identity_x1

  !> direct mapping
  function sll_f_identity_x2( xi, params )
    sll_real64 :: sll_f_identity_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_x2= xi(2)
  end function sll_f_identity_x2

  !> direct mapping
  function sll_f_identity_x3( xi, params )
    sll_real64 :: sll_f_identity_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_x3 = xi(3)
  end function sll_f_identity_x3

    !> inverse mapping
  function sll_f_identity_xi1( x, params )
    sll_real64 :: sll_f_identity_xi1
    sll_real64, intent(in) ::  x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_xi1 = x(1)
  end function sll_f_identity_xi1

  !> inverse mapping
  function sll_f_identity_xi2( x, params )
    sll_real64 :: sll_f_identity_xi2
    sll_real64, intent(in) ::  x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_xi2= x(2)
  end function sll_f_identity_xi2

  !> inverse mapping
  function sll_f_identity_xi3( x, params )
    sll_real64 :: sll_f_identity_xi3
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_xi3 = x(3)
  end function sll_f_identity_xi3

  !> jacobian matrix
  function sll_f_identity_jac11 ( xi, params )
    sll_real64  :: sll_f_identity_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac11 = 1._f64
  end function sll_f_identity_jac11

  !> jacobian matrix
  function sll_f_identity_jac12 ( xi, params )
    sll_real64  :: sll_f_identity_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac12 = 0._f64
  end function sll_f_identity_jac12

  !> jacobian matrix
  function sll_f_identity_jac13 ( xi, params )
    sll_real64  :: sll_f_identity_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac13 = 0._f64
  end function sll_f_identity_jac13

  !> jacobian matrix
  function sll_f_identity_jac21 ( xi, params )
    sll_real64  :: sll_f_identity_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac21 = 0._f64
  end function sll_f_identity_jac21

  !> jacobian matrix
  function sll_f_identity_jac22 ( xi, params )
    sll_real64  :: sll_f_identity_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac22 = 1._f64
  end function sll_f_identity_jac22

  !> jacobian matrix
  function sll_f_identity_jac23 ( xi, params )
    sll_real64  :: sll_f_identity_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac23 = 0._f64
  end function sll_f_identity_jac23

  !> jacobian matrix
  function sll_f_identity_jac31 ( xi, params )
    sll_real64  :: sll_f_identity_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac31 = 0._f64
  end function sll_f_identity_jac31

  !> jacobian matrix
  function sll_f_identity_jac32 ( xi, params )
    sll_real64  :: sll_f_identity_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac32 = 0._f64
  end function sll_f_identity_jac32

  !> jacobian matrix
  function sll_f_identity_jac33 ( xi, params )
    sll_real64  :: sll_f_identity_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_identity_jac33 = 1._f64
  end function sll_f_identity_jac33


  ! ***************************************************************************
  !
  ! Affine  coordinate transformation:
  !
  ! X1 = a1*xi(1)+a2*xi(2)+a3*xi(3)+b1
  ! X2 = a4*xi(1)+a5*xi(2)+a6*xi(3)+b2
  ! X3 = a7*xi(1)+a8*xi(2)+a9*xi(3)+b3
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (a1,a2,a3,a4,a5,a,a7,a8,a9,b1,b2,b3).
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_affine_x1( xi, params )
    sll_real64 :: sll_f_affine_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_x1 = params(1)*xi(1)+params(2)*xi(2)+params(3)*xi(3)+params(10)
  end function sll_f_affine_x1

  !> direct mapping
  function sll_f_affine_x2( xi, params )
    sll_real64 :: sll_f_affine_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_x2=  params(4)*xi(1)+params(5)*xi(2)+params(6)*xi(3)+params(11)
  end function sll_f_affine_x2

  !> direct mapping
  function sll_f_affine_x3( xi, params )
    sll_real64 :: sll_f_affine_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_x3 =  params(7)*xi(1)+params(8)*xi(2)+params(9)*xi(3)+params(12)
  end function sll_f_affine_x3

  !> jacobian matrix
  function sll_f_affine_jac11 ( xi, params )
    sll_real64  :: sll_f_affine_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac11 = params(1)
  end function sll_f_affine_jac11

  !> jacobian matrix
  function sll_f_affine_jac12 ( xi, params )
    sll_real64  :: sll_f_affine_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac12 = params(2)
  end function sll_f_affine_jac12

  !> jacobian matrix
  function sll_f_affine_jac13 ( xi, params )
    sll_real64  :: sll_f_affine_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac13 = params(3)
  end function sll_f_affine_jac13

  !> jacobian matrix
  function sll_f_affine_jac21 ( xi, params )
    sll_real64  :: sll_f_affine_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac21 = params(4)
  end function sll_f_affine_jac21

  !> jacobian matrix
  function sll_f_affine_jac22 ( xi, params )
    sll_real64  :: sll_f_affine_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac22 = params(5)
  end function sll_f_affine_jac22

  !> jacobian matrix
  function sll_f_affine_jac23 ( xi, params )
    sll_real64  :: sll_f_affine_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac23 = params(6)
  end function sll_f_affine_jac23

  !> jacobian matrix
  function sll_f_affine_jac31 ( xi, params )
    sll_real64  :: sll_f_affine_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac31 = params(7)
  end function sll_f_affine_jac31

  !> jacobian matrix
  function sll_f_affine_jac32 ( xi, params )
    sll_real64  :: sll_f_affine_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    sll_f_affine_jac32 = params(8)
  end function sll_f_affine_jac32

  !> jacobian matrix
  function sll_f_affine_jac33 ( xi, params )
    sll_real64  :: sll_f_affine_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 12)
    sll_f_affine_jac33 = params(9)
  end function sll_f_affine_jac33


  ! ***************************************************************************
  !
  ! Scaling  coordinate transformation:
  !
  ! X1 = a1*xi(1)+b1
  ! X2 = a2*xi(2)+b2
  ! X3 = a3*xi(3)+b3
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (a1,a2,a3,b1,b2,b3).
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_scaling_x1( xi, params )
    sll_real64 :: sll_f_scaling_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_x1 = params(1)*xi(1)+params(4)
  end function sll_f_scaling_x1

  !> direct mapping
  function sll_f_scaling_x2( xi, params )
    sll_real64 :: sll_f_scaling_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_x2=  params(2)*xi(2)+params(5)
  end function sll_f_scaling_x2

  !> direct mapping
  function sll_f_scaling_x3( xi, params )
    sll_real64 :: sll_f_scaling_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_x3 = params(3)*xi(3)+params(6)
  end function sll_f_scaling_x3

  !> direct mapping
  function sll_f_scaling_xi1( x, params )
    sll_real64 :: sll_f_scaling_xi1
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_xi1 = (x(1)-params(4))/params(1)
  end function sll_f_scaling_xi1

    !> direct mapping
  function sll_f_scaling_xi2( x, params )
    sll_real64 :: sll_f_scaling_xi2
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_xi2 = (x(2)-params(5))/params(2)
  end function sll_f_scaling_xi2

  !> direct mapping
  function sll_f_scaling_xi3( x, params )
    sll_real64 :: sll_f_scaling_xi3
    sll_real64, intent(in) :: x(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_xi3 = (x(3)-params(6))/params(3)
  end function sll_f_scaling_xi3
  
  !> jacobian matrix
  function sll_f_scaling_jac11 ( xi, params )
    sll_real64  :: sll_f_scaling_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac11 = params(1)
  end function sll_f_scaling_jac11

  !> jacobian matrix
  function sll_f_scaling_jac12 ( xi, params )
    sll_real64  :: sll_f_scaling_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac12 = 0._f64
  end function sll_f_scaling_jac12

  !> jacobian matrix
  function sll_f_scaling_jac13 ( xi, params )
    sll_real64  :: sll_f_scaling_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac13 = 0._f64
  end function sll_f_scaling_jac13

  !> jacobian matrix
  function sll_f_scaling_jac21 ( xi, params )
    sll_real64  :: sll_f_scaling_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac21 = 0._f64
  end function sll_f_scaling_jac21

  !> jacobian matrix
  function sll_f_scaling_jac22 ( xi, params )
    sll_real64  :: sll_f_scaling_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac22 = params(2)
  end function sll_f_scaling_jac22

  !> jacobian matrix
  function sll_f_scaling_jac23 ( xi, params )
    sll_real64  :: sll_f_scaling_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac23 = 0._f64
  end function sll_f_scaling_jac23

  !> jacobian matrix
  function sll_f_scaling_jac31 ( xi, params )
    sll_real64  :: sll_f_scaling_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac31 = 0._f64
  end function sll_f_scaling_jac31

  !> jacobian matrix
  function sll_f_scaling_jac32 ( xi, params )
    sll_real64  :: sll_f_scaling_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac32 = 0._f64
  end function sll_f_scaling_jac32

  !> jacobian matrix
  function sll_f_scaling_jac33 ( xi, params )
    sll_real64  :: sll_f_scaling_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jac33 = params(3)
  end function sll_f_scaling_jac33

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv11 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv11 = 1._f64/params(1)
  end function sll_f_scaling_jacinv11

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv12 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv12 = 0._f64
  end function sll_f_scaling_jacinv12

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv13 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv13 = 0._f64
  end function sll_f_scaling_jacinv13

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv21 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv21 = 0._f64
  end function sll_f_scaling_jacinv21

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv22 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv22 = 1._f64/params(2)
  end function sll_f_scaling_jacinv22

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv23 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv23 = 0._f64
  end function sll_f_scaling_jacinv23

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv31 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv31 = 0._f64
  end function sll_f_scaling_jacinv31

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv32 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv32 = 0._f64
  end function sll_f_scaling_jacinv32

  !> jacobian matrix inverse
  function sll_f_scaling_jacinv33 ( xi, params )
    sll_real64  :: sll_f_scaling_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacinv33 = 1._f64/params(3)
  end function sll_f_scaling_jacinv33

  !> jacobian 
  function sll_f_scaling_jacobian ( xi, params )
    sll_real64  :: sll_f_scaling_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 6)
    sll_f_scaling_jacobian = params(1)*params(2)*params(3)
  end function sll_f_scaling_jacobian


  ! ***************************************************************************
  !
  ! Periodic  coordinate transformation:
  !
  ! X1 = L1*xi(1)
  ! X2 = L2*xi(2)+alpha*sin(2*\pi*xi(1))
  ! X3 = L3*xi(3)
  ! Where xi(3) are defined in the interval [0,1]. The 'params' array
  ! contains the information (L1,L2,L3,alpha).
  !
  ! ***************************************************************************


  !> direct mapping
  function sll_f_periodic_x1( xi, params )
    sll_real64 :: sll_f_periodic_x1
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_x1 = modulo(params(1)*xi(1), params(1))
  end function sll_f_periodic_x1

  !> direct mapping
  function sll_f_periodic_x2( xi, params )
    sll_real64 :: sll_f_periodic_x2
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_x2=  modulo(params(2)*xi(2)+params(4)*sin(sll_p_twopi*xi(1)),params(2))
  end function sll_f_periodic_x2

  !> direct mapping
  function sll_f_periodic_x3( xi, params )
    sll_real64 :: sll_f_periodic_x3
    sll_real64, intent(in) :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_x3 = modulo(params(3)*xi(3), params(3))
  end function sll_f_periodic_x3

  !> jacobian matrix
  function sll_f_periodic_jac11 ( xi, params )
    sll_real64  :: sll_f_periodic_jac11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac11 = params(1)
  end function sll_f_periodic_jac11

  !> jacobian matrix
  function sll_f_periodic_jac12 ( xi, params )
    sll_real64  :: sll_f_periodic_jac12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac12 = 0._f64
  end function sll_f_periodic_jac12

  !> jacobian matrix
  function sll_f_periodic_jac13 ( xi, params )
    sll_real64  :: sll_f_periodic_jac13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac13 = 0._f64
  end function sll_f_periodic_jac13

  !> jacobian matrix
  function sll_f_periodic_jac21 ( xi, params )
    sll_real64  :: sll_f_periodic_jac21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac21 = params(4)*sll_p_twopi*cos(sll_p_twopi*xi(1))
  end function sll_f_periodic_jac21

  !> jacobian matrix
  function sll_f_periodic_jac22 ( xi, params )
    sll_real64  :: sll_f_periodic_jac22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac22 = params(2)
  end function sll_f_periodic_jac22

  !> jacobian matrix
  function sll_f_periodic_jac23 ( xi, params )
    sll_real64  :: sll_f_periodic_jac23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac23 = 0._f64
  end function sll_f_periodic_jac23

  !> jacobian matrix
  function sll_f_periodic_jac31 ( xi, params )
    sll_real64  :: sll_f_periodic_jac31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac31 = 0._f64
  end function sll_f_periodic_jac31

  !> jacobian matrix
  function sll_f_periodic_jac32 ( xi, params )
    sll_real64  :: sll_f_periodic_jac32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac32 = 0._f64
  end function sll_f_periodic_jac32

  !> jacobian matrix
  function sll_f_periodic_jac33 ( xi, params )
    sll_real64  :: sll_f_periodic_jac33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jac33 = params(3)
  end function sll_f_periodic_jac33

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv11 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv11
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv11 = 1._f64/params(1)
  end function sll_f_periodic_jacinv11

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv12 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv12
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv12 = 0._f64
  end function sll_f_periodic_jacinv12

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv13 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv13
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv13 = 0._f64
  end function sll_f_periodic_jacinv13

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv21 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv21
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv21 = - params(4)*sll_p_twopi*cos(sll_p_twopi*xi(1))/params(1)
  end function sll_f_periodic_jacinv21

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv22 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv22
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv22 = 1._f64/params(2)
  end function sll_f_periodic_jacinv22

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv23 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv23
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv23 = 0._f64
  end function sll_f_periodic_jacinv23

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv31 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv31
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv31 = 0._f64
  end function sll_f_periodic_jacinv31

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv32 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv32
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv32 = 0._f64
  end function sll_f_periodic_jacinv32

  !> jacobian matrix inverse
  function sll_f_periodic_jacinv33 ( xi, params )
    sll_real64  :: sll_f_periodic_jacinv33
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacinv33 = 1._f64/params(3)
  end function sll_f_periodic_jacinv33

  !> jacobian 
  function sll_f_periodic_jacobian ( xi, params )
    sll_real64  :: sll_f_periodic_jacobian
    sll_real64, intent(in)   :: xi(3)
    sll_real64, dimension(:), intent(in) :: params
    SLL_ASSERT(size(params) >= 4)
    sll_f_periodic_jacobian = params(1)*params(2)*params(3)
  end function sll_f_periodic_jacobian

  




end module sll_m_3d_coordinate_transformations
