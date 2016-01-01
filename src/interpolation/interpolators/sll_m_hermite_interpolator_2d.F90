!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup interpolators
!>@brief
!>Hermite interpolation in 2d
!!@details
!!derivatives are given with finite stencil formulae of order p
!!which can be arbitrary in each direction
!!- If p is odd, the reconstruction has discontinuous derivatives
!!- If p is even, the reconstruction has continuous derivatives
!!- p=6 should be quite similar to cubic splines
!!
!!do not hesitate to take large odd p, like the favourite p=17
!!@todo
!! for the moment only in implementation for the case DIRICHLET x PERIODIC
module sll_m_hermite_interpolator_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_hermite_interpolation_2d, only: &
    sll_s_compute_interpolants_hermite_2d, &
    sll_f_interpolate_value_hermite_2d, &
    sll_f_new_hermite_interpolation_2d, &
    sll_t_hermite_interpolation_2d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  implicit none

  public :: &
    sll_f_new_hermite_interpolator_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> The hermite-based interpolator is only a wrapper around the capabilities
!! of the hermite interpolation. 
!> @details
!! All interpolators share a common interface with
!! respect to their use, as described by the interpolator_2d_base class.
!! Where the diverse interpolators diverge is in the way to initialize them.
!! We basically copy the analog for cubic splines
  type, extends(sll_c_interpolator_2d) :: sll_hermite_interpolator_2d
    !> PLEASE ADD DOCUMENTATION
    type(sll_t_hermite_interpolation_2d), pointer :: hermite
    !> PLEASE ADD DOCUMENTATION
    sll_int32 :: npts1
    !> PLEASE ADD DOCUMENTATION
    sll_int32 :: npts2
  contains
    !> PLEASE ADD DOCUMENTATION
    procedure, pass(interpolator) :: initialize=>initialize_hermite_interpolator_2d
    !> PLEASE ADD DOCUMENTATION
    procedure :: compute_interpolants => wrap_compute_interpolants_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure :: interpolate_from_interpolant_value => wrap_interpolate_value_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure :: interpolate_from_interpolant_derivative_eta1 => wrap_interpolate_deriv1_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure :: interpolate_from_interpolant_derivative_eta2 => wrap_interpolate_deriv2_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: interpolate_array => wrap_interpolate_array_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: interpolate_array_disp => wrap_interpolate2d_disp_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: set_coefficients => wrap_set_coefficients_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: get_coefficients => wrap_get_coefficients_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: coefficients_are_set => wrap_coefficients_are_set_hermite_2d
    !> PLEASE ADD DOCUMENTATION
    procedure, pass :: delete => delete_sll_hermite_interpolator_2d
  end type sll_hermite_interpolator_2d


contains

    !> PLEASE ADD DOCUMENTATION
  function sll_f_new_hermite_interpolator_2d( &
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    degree1, &
    degree2, &
    eta1_hermite_continuity, &
    eta2_hermite_continuity, &
    eta1_bc_type, &
    eta2_bc_type, &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes ) &   
    result(interpolator)

    type(sll_hermite_interpolator_2d), pointer :: interpolator
    sll_int32, intent(in) :: npts1
    sll_int32, intent(in) :: npts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: degree1
    sll_int32, intent(in) :: degree2    
    sll_int32, intent(in) :: eta1_hermite_continuity
    sll_int32, intent(in) :: eta2_hermite_continuity
    sll_int32, intent(in) :: eta1_bc_type
    sll_int32, intent(in) :: eta2_bc_type
    sll_real64, intent(in), optional :: const_eta1_min_slope
    sll_real64, intent(in), optional :: const_eta1_max_slope
    sll_real64, intent(in), optional :: const_eta2_min_slope
    sll_real64, intent(in), optional :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interpolator,ierr)
    
    interpolator%npts1 = npts1
    interpolator%npts2 = npts2
    
    call interpolator%initialize( &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      degree1, &
      degree2, &
      eta1_hermite_continuity, &
      eta2_hermite_continuity, &
      eta1_bc_type, &
      eta2_bc_type, &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes )    

     
  end function  sll_f_new_hermite_interpolator_2d
  
  subroutine initialize_hermite_interpolator_2d( &
    interpolator, &    
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    degree1, &
    degree2, &
    eta1_hermite_continuity, &
    eta2_hermite_continuity, &
    eta1_bc_type, &
    eta2_bc_type, &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes )    

    class(sll_hermite_interpolator_2d), intent(inout) :: interpolator
    sll_int32, intent(in) :: npts1
    sll_int32, intent(in) :: npts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: degree1
    sll_int32, intent(in) :: degree2    
    sll_int32, intent(in) :: eta1_hermite_continuity
    sll_int32, intent(in) :: eta2_hermite_continuity
    sll_int32, intent(in) :: eta1_bc_type
    sll_int32, intent(in) :: eta2_bc_type
    sll_real64, intent(in), optional :: const_eta1_min_slope
    sll_real64, intent(in), optional :: const_eta1_max_slope
    sll_real64, intent(in), optional :: const_eta2_min_slope
    sll_real64, intent(in), optional :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes
       
    interpolator%hermite  => sll_f_new_hermite_interpolation_2d( &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      degree1, &
      degree2, &
      eta1_hermite_continuity, &
      eta2_hermite_continuity, &
      eta1_bc_type, &
      eta2_bc_type, &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes )    
      
  end subroutine initialize_hermite_interpolator_2d
  
  
  subroutine wrap_compute_interpolants_hermite_2d( &
    interpolator, &
    data_array, &
    eta1_coords, &
    size_eta1_coords, &
    eta2_coords, &
    size_eta2_coords )
    class(sll_hermite_interpolator_2d), intent(inout) :: interpolator
    sll_real64, dimension(:,:), intent(in) :: data_array
    sll_real64, dimension(:), intent(in),optional   :: eta1_coords
    sll_real64, dimension(:), intent(in),optional   :: eta2_coords
    sll_int32, intent(in), optional                 :: size_eta1_coords
    sll_int32, intent(in),optional                  :: size_eta2_coords

    if(present(eta1_coords))then
      !print *,'#Warning eta1_coords not used'
    endif
    if(present(eta2_coords))then
      !print *,'#Warning eta2_coords not used'
    endif
    if(present(size_eta1_coords))then
      !print *,'#Warning size_eta1_coords not used'
    endif
    if(present(size_eta2_coords))then
      !print *,'#Warning size_eta2_coords not used'
    endif    
    call sll_s_compute_interpolants_hermite_2d( interpolator%hermite, data_array )
  end subroutine wrap_compute_interpolants_hermite_2d
  
  function wrap_interpolate_value_hermite_2d( interpolator, eta1, eta2 ) result(val)
    class(sll_hermite_interpolator_2d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = sll_f_interpolate_value_hermite_2d( eta1, eta2, interpolator%hermite )
      
  end function wrap_interpolate_value_hermite_2d
  
  function wrap_interpolate_deriv1_hermite_2d( interpolator, eta1, eta2 ) result(val)
    class(sll_hermite_interpolator_2d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = 0._f64+eta1+eta2
    print *,'#wrap_interpolate_deriv1_hermite_2d'
    print *,'#not implemented for the moment'
    stop
    SLL_ASSERT(interpolator%npts1>0)
    !interpolate_x1_derivative_2D(eta1,eta2,interpolator%spline)
  end function wrap_interpolate_deriv1_hermite_2d

  function wrap_interpolate_deriv2_hermite_2d( interpolator, eta1, eta2 ) result(val)
    class(sll_hermite_interpolator_2d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = 0._f64+eta1+eta2
    print *,'#wrap_interpolate_deriv1_hermite_2d'
    print *,'#not implemented for the moment'
    stop
    SLL_ASSERT(interpolator%npts1>0)
    !interpolate_x1_derivative_2D(eta1,eta2,interpolator%spline)
  end function wrap_interpolate_deriv2_hermite_2d

  subroutine wrap_interpolate_array_hermite_2d( &
    this, &
    num_points1, &
    num_points2, &
    data_in, &
    eta1, &
    eta2, &
    data_out)
    class(sll_hermite_interpolator_2d),  intent(in)  :: this
    sll_int32,  intent(in)                           :: num_points1
    sll_int32,  intent(in)                           :: num_points2
    sll_real64, dimension(:,:), intent(in)           :: eta1
    sll_real64, dimension(:,:), intent(in)           :: eta2
    sll_real64, dimension(:,:), intent(in)           :: data_in
    sll_real64,                 intent(out)          :: data_out(num_points1,num_points2) 
    sll_int32 :: i
    sll_int32 :: j
    call sll_s_compute_interpolants_hermite_2d( this%hermite, data_in )
    do j = 1, num_points2
      do i = 1, num_points1
        data_out(i,j) = this%interpolate_from_interpolant_value(eta1(i,j),eta2(i,j))
      end do
    end do
    !print *,'#wrap_interpolate_array_hermite_2d'
    !print *,'#not implemented for the moment'
    !stop
  end subroutine wrap_interpolate_array_hermite_2d

  subroutine wrap_interpolate2d_disp_hermite_2d( &
       this, &
       num_points1, &
       num_points2, &
       data_in, &
       alpha1, &
       alpha2, &
       data_out)
    class(sll_hermite_interpolator_2d), intent(in) :: this
    sll_int32,  intent(in)                         :: num_points1
    sll_int32,  intent(in)                         :: num_points2
    sll_real64, dimension(:,:), intent(in)         :: alpha1
    sll_real64, dimension(:,:), intent(in)         :: alpha2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_real64,                 intent(out)        :: data_out(num_points1,num_points2)
    print *,'#wrap_interpolate2d_disp_hermite_2d'
    print *,'#not implemented for the moment'
    data_out = 0.0_f64 * data_in
    stop
    SLL_ASSERT(size(alpha1,1) == num_points1)
    SLL_ASSERT(size(alpha2,1) == num_points1)
    SLL_ASSERT(this%npts1     == num_points1)
  end subroutine wrap_interpolate2d_disp_hermite_2d

  
  subroutine wrap_set_coefficients_hermite_2d( &
    interpolator, &
    coeffs_1d, &
    coeffs_2d, &
    coeff2d_size1, &
    coeff2d_size2, &
    knots1, &
    size_knots1, &
    knots2, &
    size_knots2)
    class(sll_hermite_interpolator_2d),  intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in), optional :: coeffs_1d
    sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
    sll_int32, intent(in), optional :: coeff2d_size1
    sll_int32, intent(in), optional :: coeff2d_size2
    sll_real64, dimension(:), intent(in), optional   :: knots1
    sll_real64, dimension(:), intent(in), optional   :: knots2
    sll_int32, intent(in), optional :: size_knots1
    sll_int32, intent(in), optional :: size_knots2
    print *, 'wrap_set_coefficients_hermite_2d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%npts1
    if(present(coeffs_1d))then
      print *,'coeffs_1d present but not used'
    endif     
    if(present(coeffs_2d))then
      print *,'coeffs_2d present but not used'
    endif     
    if(present(coeff2d_size1))then
      print *,'coeff2d_size1 present but not used'
    endif     
    if(present(coeff2d_size2))then
      print *,'coeff2d_size2 present but not used'
    endif     
    if(present(knots1))then
      print *,'knots1 present but not used'
    endif     
    if(present(knots2))then
      print *,'knots2 present but not used'
    endif     
    if(present(size_knots1))then
      print *,'size_knots1 present but not used'
    endif     
    if(present(size_knots2))then
      print *,'size_knots2 present but not used'
    endif     
    stop
  end subroutine wrap_set_coefficients_hermite_2d

  function wrap_get_coefficients_hermite_2d(interpolator) result(res)
    class(sll_hermite_interpolator_2d), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer            :: res     
    
    print *, 'wrap_get_coefficients_hermite_2d: ERROR: This function has not been ', &
         'implemented yet.'
    res => null()
    print *,interpolator%npts1    
    stop      
  end function wrap_get_coefficients_hermite_2d

  function wrap_coefficients_are_set_hermite_2d( interpolator ) result(res)
    class(sll_hermite_interpolator_2d), intent(in) :: interpolator
    logical :: res
    res = .false.
    print *, 'wrap_coefficients_are_set_hermite_2d(): '
    print *, 'this function has not been implemented yet.'
    print *,'#',interpolator%npts1
    !stop
  end function wrap_coefficients_are_set_hermite_2d

  subroutine delete_sll_hermite_interpolator_2d( interpolator )
    class(sll_hermite_interpolator_2d), intent(inout) :: interpolator    
    print *,'#warning delete_sll_hermite_interpolator_2d'
    print *,'#not implemented for the moment'     
    SLL_ASSERT(interpolator%npts1>0)
  end subroutine delete_sll_hermite_interpolator_2d  



end module sll_m_hermite_interpolator_2d
