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

module sll_cubic_spline_interpolator_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

#ifndef STDF95
  use sll_module_interpolators_2d_base
#endif
  use sll_cubic_splines
  implicit none
  
  ! The spline-based interpolator is only a wrapper around the capabilities
  ! of the cubic splines. All interpolators share a common interface with
  ! respect to their use, as described by the interpolator_2d_base class.
  !
  ! Where the diverse interpolators diverge is in the way to initialize them.
#ifdef STDF95
  type                                :: cubic_spline_2d_interpolator
#else
  type, extends(sll_interpolator_2d_base) :: cubic_spline_2d_interpolator
#endif
     sll_int32                           :: npts1
     sll_int32                           :: npts2
     type(sll_cubic_spline_2D), pointer        :: spline
     sll_int32                           :: bc_type1
     sll_int32                           :: bc_type2
     sll_real64, dimension(:,:), pointer :: interpolation_points 
#ifdef STDF95
#else
   contains
     procedure, pass(interpolator) :: initialize=>initialize_cs2d_interpolator
     procedure :: compute_interpolants => compute_interpolants_cs2d
     procedure :: interpolate_value => interpolate_value_cs2d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_cs2d
     procedure :: interpolate_derivative_eta2 => interpolate_deriv2_cs2d
     procedure, pass :: interpolate_array => spline_interpolate2d
     procedure, pass :: interpolate_array_disp => spline_interpolate2d_disp
     procedure, pass :: set_coefficients => set_coefficients_cs2d
     procedure, pass :: get_coefficients => get_coefficients_cs2d
     procedure, pass :: coefficients_are_set => coefficients_are_set_cs2d
     procedure, pass :: delete => delete_cubic_spline_2d_interpolator
    ! procedure, pass :: compute_spline_coefficients => compute_spl_coeff_cs2d
#endif
  end type cubic_spline_2d_interpolator
  
  interface delete
     module procedure delete_cubic_spline_2d_interpolator
  end interface delete

contains

  subroutine delete_cubic_spline_2d_interpolator( interpolator )
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    call delete(interpolator%spline)
  end subroutine delete_cubic_spline_2d_interpolator
  
  function new_cubic_spline_2d_interpolator( &
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_bc_type,   &
    eta2_bc_type,   &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes ) &
    result(interpolator)

    type(cubic_spline_2d_interpolator), pointer :: interpolator
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    sll_real64, intent(in)                        :: eta1_min
    sll_real64, intent(in)                        :: eta1_max
    sll_real64, intent(in)                        :: eta2_min
    sll_real64, intent(in)                        :: eta2_max
    sll_int32, intent(in), optional               :: eta1_bc_type
    sll_int32, intent(in), optional               :: eta2_bc_type
    sll_real64, intent(in), optional              :: const_eta1_min_slope
    sll_real64, intent(in), optional              :: const_eta1_max_slope
    sll_real64, intent(in), optional              :: const_eta2_min_slope
    sll_real64, intent(in), optional              :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interpolator,ierr)
    
    call interpolator%initialize( &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      eta1_bc_type,   &
      eta2_bc_type,   &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes )

     
  end function  new_cubic_spline_2d_interpolator


  ! We allow to use the enumerators of the splines module in this interpolator
  ! because:
  ! a. This is just a wrapper and is intimately related to the underlying
  !    cubic splines module.
  ! b. There is no uniform interface for the initialization anyway.
  ! The underlying implementation with the splines module could be hidden but
  ! I can't see a compelling reason why.
#ifdef STDF95
  subroutine cubic_spline_2d_initialize( &
#else
  subroutine initialize_cs2d_interpolator( &
#endif
    interpolator, &
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_bc_type,   &
    eta2_bc_type,   &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes )

#ifdef STDF95
    type(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#else
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#endif
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    sll_real64, intent(in)                        :: eta1_min
    sll_real64, intent(in)                        :: eta1_max
    sll_real64, intent(in)                        :: eta2_min
    sll_real64, intent(in)                        :: eta2_max
    sll_int32, intent(in), optional               :: eta1_bc_type
    sll_int32, intent(in), optional               :: eta2_bc_type
    sll_real64, intent(in), optional              :: const_eta1_min_slope
    sll_real64, intent(in), optional              :: const_eta1_max_slope
    sll_real64, intent(in), optional              :: const_eta2_min_slope
    sll_real64, intent(in), optional              :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes

    interpolator%npts1 = npts1
    interpolator%npts2 = npts2
    interpolator%bc_type1 = eta1_bc_type
    interpolator%bc_type2 = eta2_bc_type
    interpolator%spline => new_spline_2D( &
         npts1, &
         npts2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         eta1_bc_type, &
         eta2_bc_type, &
         const_slope_x1_min=const_eta1_min_slope, &
         const_slope_x1_max=const_eta1_max_slope, &
         const_slope_x2_min=const_eta2_min_slope, &
         const_slope_x2_max=const_eta2_max_slope, &
         x1_min_slopes=eta1_min_slopes, &
         x1_max_slopes=eta1_max_slopes, &
         x2_min_slopes=eta2_min_slopes, &
         x2_max_slopes=eta2_max_slopes )
  end subroutine

#ifdef STDF95
  subroutine cubic_spline_2d_compute_interpolants( &
#else
  subroutine compute_interpolants_cs2d( &
#endif
       interpolator, &
       data_array, &
       eta1_coords, &
       size_eta1_coords, &
       eta2_coords, &
       size_eta2_coords )
#ifdef STDF95
    type(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#else
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:,:), intent(in) :: data_array
    sll_real64, dimension(:), intent(in),optional   :: eta1_coords
    sll_real64, dimension(:), intent(in),optional   :: eta2_coords
    sll_int32, intent(in), optional                 :: size_eta1_coords
    sll_int32, intent(in),optional                  :: size_eta2_coords
    call compute_spline_2D( data_array, interpolator%spline )
  end subroutine

#ifdef STDF95
  function cubic_spline_2d_interpolate_value( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_value_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_value_2D( eta1, eta2, interpolator%spline )
  end function

#ifdef STDF95
  function cubic_spline_2d_interpolate_derivative_eta1( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_deriv1_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_x1_derivative_2D(eta1,eta2,interpolator%spline)
  end function

#ifdef STDF95
  function cubic_spline_2d_interpolate_derivative_eta2( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_deriv2_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2

    val = interpolate_x2_derivative_2D(eta1,eta2,interpolator%spline)

  end function

#ifdef STDF95
  function cubic_spline_2d_interpolate_array(this, num_points1, num_points2, data_in, &
                                eta1, eta2) &
       result(data_out)
    type(cubic_spline_2d_interpolator),  intent(in)       :: this
#else
  function spline_interpolate2d(this, num_points1, num_points2, data_in, &
                                eta1, eta2) &
       result(data_out)

    class(cubic_spline_2d_interpolator),  intent(in) :: this
#endif
    sll_int32,  intent(in)                           :: num_points1
    sll_int32,  intent(in)                           :: num_points2
    sll_real64, dimension(:,:), intent(in)           :: eta1
    sll_real64, dimension(:,:), intent(in)           :: eta2
    sll_real64, dimension(:,:), intent(in)           :: data_in
    sll_real64, dimension(num_points1,num_points2)   :: data_out
    ! local variables
    sll_int32 :: i,j
    ! compute the interpolating spline coefficients
    call compute_spline_2D( data_in, this%spline )
    do j = 1, num_points2
    do i = 1, num_points1
#ifdef STDF95
        data_out(i,j) = cubic_spline_2d_interpolate_value(this,eta1(i,j),eta2(i,j))     
#else
        data_out(i,j) = this%interpolate_value(eta1(i,j),eta2(i,j))
#endif
    end do
    end do

  end function !spline_interpolate2d

#ifdef STDF95

  function cubic_spline_interpolate2d_disp(this,        &
                                           num_points1, &
                                           num_points2, &
                                           data_in,     &
                                           alpha1,      &
                                           alpha2) result(data_out)

    type(cubic_spline_2d_interpolator), intent(in) :: this

#else

  function spline_interpolate2d_disp(this,        &
                                     num_points1, &
                                     num_points2, &
                                     data_in,     &
                                     alpha1,      &
                                     alpha2) result(data_out)

    class(cubic_spline_2d_interpolator),  intent(in) :: this

#endif

    sll_int32,  intent(in)                         :: num_points1
    sll_int32,  intent(in)                         :: num_points2
    sll_real64, dimension(:,:), intent(in)         :: alpha1
    sll_real64, dimension(:,:), intent(in)         :: alpha2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_real64, dimension(num_points1,num_points2) :: data_out
    sll_real64                                     :: eta1
    sll_real64                                     :: eta1_min
    sll_real64                                     :: eta1_max
    sll_real64                                     :: delta_eta1
    sll_real64                                     :: eta2
    sll_real64                                     :: eta2_min
    sll_real64                                     :: eta2_max
    sll_real64                                     :: delta_eta2
    sll_int32                                      :: i
    sll_int32                                      :: j

    eta1_min   = this%spline%x1_min 
    eta1_max   = this%spline%x1_max 
    eta2_min   = this%spline%x2_min 
    eta2_max   = this%spline%x2_max 
    delta_eta1 = this%spline%x1_delta  
    delta_eta2 = this%spline%x2_delta  
    
    call compute_spline_2D( data_in, this%spline )

    if(this%bc_type1 == SLL_PERIODIC .and. &
       this%bc_type2 == SLL_PERIODIC ) then
       
       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1
             eta2 = eta2_min + (j-1)*delta_eta2
             eta1 = eta1_min + &
                  modulo(eta1-eta1_min-alpha1(i,j),eta1_max-eta1_min)
             eta2 = eta2_min + &
                  modulo(eta2-eta2_min-alpha2(i,j),eta2_max-eta2_min)
#ifdef STDF95
          data_out(i,j) = cubic_spline_2d_interpolate_value(this,eta1,eta2)     
#else
             data_out(i,j) = this%interpolate_value(eta1,eta2)
#endif
          end do
       end do

    else if(this%bc_type1 == SLL_HERMITE .and. &
            this%bc_type2 == SLL_HERMITE ) then
       
       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1
             eta2 = eta2_min + (j-1)*delta_eta2
             eta1 = min(eta1,eta1_max)
             eta2 = min(eta2,eta2_max)
             eta1 = max(eta1,eta1_min)
             eta2 = max(eta2,eta2_min)
#ifdef STDF95
             data_out(i,j) = cubic_spline_2d_interpolate_value(this,eta1,eta2)
#else
             data_out(i,j) = this%interpolate_value(eta1,eta2)
#endif
          end do
       end do
       
    else

       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1 - alpha1(i,j)
             eta2 = eta2_min + (j-1)*delta_eta2 - alpha2(i,j)
             SLL_ASSERT(eta1_min <= eta1 .and. eta1 <= eta1_max)
             SLL_ASSERT(eta2_min <= eta2 .and. eta2 <= eta2_max)
#ifdef STDF95
             data_out(i,j) = cubic_spline_2d_interpolate_value(this,eta1,eta2)
#else
             data_out(i,j) = this%interpolate_value(eta1,eta2)
#endif
          end do
       end do
    end if
  end function 

#ifdef STDF95
  subroutine cubic_spline_2d_set_coefficients(&
       interpolator,&
       coeffs_1d,&
       coeffs_2d,&
       coeff2d_size1,&
       coeff2d_size2,&
       knots1,&
       size_knots1,&
       knots2,&
       size_knots2)
    type (cubic_spline_2d_interpolator),  intent(inout) :: interpolator
#else
  subroutine set_coefficients_cs2d( &
       interpolator,&
       coeffs_1d,&
       coeffs_2d,&
       coeff2d_size1,&
       coeff2d_size2,&
       knots1,&
       size_knots1,&
       knots2,&
       size_knots2)
    class(cubic_spline_2d_interpolator),  intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in), optional :: coeffs_1d
    sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
    ! size coeffs 2D 
    sll_int32, intent(in), optional :: coeff2d_size1
    sll_int32, intent(in), optional :: coeff2d_size2
    sll_real64, dimension(:), intent(in), optional   :: knots1
    sll_real64, dimension(:), intent(in), optional   :: knots2
    sll_int32, intent(in), optional :: size_knots1
    sll_int32, intent(in), optional :: size_knots2
    print *, 'set_coefficients_cs2d(): ERROR: This function has not been ', &
         'implemented yet.'
    stop
  end subroutine !set_coefficients_cs2d
  
!!$  subroutine compute_spl_coeff_cs2d(interpolator, &
!!$       data_array, &
!!$       eta1_coords, &
!!$       size_eta1_coords, &
!!$       eta2_coords, &
!!$       size_eta2_coords )
!!$    class(cubic_spline_2d_interpolator), intent(inout)  :: interpolator
!!$    sll_real64, dimension(:,:), intent(in)     :: data_array
!!$    sll_real64, dimension(:), intent(in),optional       :: eta1_coords
!!$    sll_real64, dimension(:), intent(in),optional       :: eta2_coords
!!$    sll_int32, intent(in), optional                     :: size_eta1_coords
!!$    sll_int32, intent(in),optional                      :: size_eta2_coords
!!$
!!$    print *, 'compute_coefficients_cs2d(): ERROR: This function has not been',&
!!$         'implemented yet.'
!!$    stop
!!$  end subroutine compute_spl_coeff_cs2d
  

#ifdef STDF95
  function cubic_spline_2d_get_coefficients(interpolator)
    type (cubic_spline_2d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer            :: cubic_spline_2d_get_coefficients     
    
    print *, 'cubic_spline_2d_get_coefficients(): ERROR: This function has not been ', &
         'implemented yet.' 
  end function cubic_spline_2d_get_coefficients
#else  
  function get_coefficients_cs2d(interpolator)
    class(cubic_spline_2d_interpolator), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer            :: get_coefficients_cs2d     
    
    print *, 'get_coefficients_cs2d(): ERROR: This function has not been ', &
         'implemented yet.' 
  end function get_coefficients_cs2d
#endif

  function coefficients_are_set_cs2d( interpolator ) result(res)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
    logical :: res
    print *, 'coefficients_are_set_cs2d(): this function has not been implemented yet.'
  end function coefficients_are_set_cs2d

end module sll_cubic_spline_interpolator_2d
