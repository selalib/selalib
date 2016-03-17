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
!> @brief
!> Class for the bspline inmplementation of sll_c_interpolator_2d
!> @details
!> Implements the sll_c_interpolator_2d interface
module sll_m_bspline_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_bspline_interpolation, only: &
    sll_o_compute_bspline_2d, &
    sll_f_interpolate_value_2d, &
    sll_t_bspline_interpolation_2d, &
    sll_s_bspline_interpolation_2d_init

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  implicit none

  public :: &
    sll_f_new_bspline_interpolator_2d, &
    sll_t_bspline_interpolator_2d, &
    sll_o_delete

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> @brief
!> The spline-based interpolator is only a wrapper around the capabilities
!> of the bplines. 
!> @details
!> All interpolators share a common interface with
!> respect to their use, as described by the interpolator_2d_base class.
!> Where the diverse interpolators diverge is in the way to initialize them.
type, extends(sll_c_interpolator_2d) :: sll_t_bspline_interpolator_2d

  sll_int32                      :: npts1    !< Number of points along x direction
  sll_int32                      :: npts2    !< Number of points along y direction
  type(sll_t_bspline_2d), pointer  :: spline   !< The spline object to store coefficients
  sll_int32                      :: spline_degree1 !< Boundady condition for x
  sll_int32                      :: spline_degree2 !< Boundary condition for y
  sll_int32                      :: bc_type1       !< Boundady condition for x
  sll_int32                      :: bc_type2       !< Boundary condition for y
  sll_real64, pointer            :: interpolation_points(:,:) !< positions

contains

  !> Allocate data, set dimensions and boundary conditions
  procedure, pass(interpolator) :: initialize=>initialize_bs2d_interpolator
  !> Compute bspline coefficients
  procedure :: compute_interpolants => compute_interpolants_bs2d
  !> Interpolate single value from last interpolants computed
  procedure :: interpolate_from_interpolant_value => interpolate_value_bs2d
  !> Interpolate first derivative from last inteprolants computed
  procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_bs2d
  !> Interpolate first derivative from last inteprolants computed
  procedure :: interpolate_from_interpolant_derivative_eta2 => interpolate_deriv2_bs2d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: interpolate_array => spline_interpolate2d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: interpolate_array_disp => spline_interpolate2d_disp
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: set_coefficients => set_coefficients_bs2d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: get_coefficients => get_coefficients_bs2d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: coefficients_are_set => coefficients_are_set_bs2d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: delete => delete_bspline_interpolator_2d
  ! procedure, pass :: compute_spline_coefficients => compute_spl_coeff_bs2d

end type sll_t_bspline_interpolator_2d

!> Pointer to this interpolator derived type
type :: sll_bspline_interpolator_2d_ptr
  type(sll_t_bspline_interpolator_2d), pointer :: interp
end type sll_bspline_interpolator_2d_ptr
  
!> Deallocate the interpolator object
interface sll_o_delete
  module procedure delete_bspline_interpolator_2d
end interface sll_o_delete


contains

  subroutine delete_bspline_interpolator_2d( interpolator )
    class(sll_t_bspline_interpolator_2d), intent(inout) :: interpolator
!PN    call sll_o_delete(interpolator%spline)
    
    SLL_ASSERT(interpolator%npts1>0)
  end subroutine delete_bspline_interpolator_2d
  
  !> Create a pointer to a 2d interpolator using bsplines.
  function sll_f_new_bspline_interpolator_2d( &
    npts1,                              &
    npts2,                              &
    eta1_min,                           &
    eta1_max,                           &
    eta2_min,                           &
    eta2_max,                           &
    spline_degree1,                     &
    spline_degree2,                     &
    eta1_bc_type,                       &
    eta2_bc_type,                       &
    const_eta1_min_slope,               &
    const_eta1_max_slope,               &
    const_eta2_min_slope,               &
    const_eta2_max_slope,               &
    eta1_min_slopes,                    &
    eta1_max_slopes,                    &
    eta2_min_slopes,                    &
    eta2_max_slopes )                   &
    result(interpolator)

    type(sll_t_bspline_interpolator_2d),    pointer  :: interpolator
    sll_int32,                intent(in)           :: npts1
    sll_int32,                intent(in)           :: npts2
    sll_real64,               intent(in)           :: eta1_min
    sll_real64,               intent(in)           :: eta1_max
    sll_real64,               intent(in)           :: eta2_min
    sll_real64,               intent(in)           :: eta2_max
    sll_int32,                intent(in)           :: spline_degree1
    sll_int32,                intent(in)           :: spline_degree2
    sll_int32,                intent(in)           :: eta1_bc_type
    sll_int32,                intent(in)           :: eta2_bc_type
    sll_real64,               intent(in), optional :: const_eta1_min_slope
    sll_real64,               intent(in), optional :: const_eta1_max_slope
    sll_real64,               intent(in), optional :: const_eta2_min_slope
    sll_real64,               intent(in), optional :: const_eta2_max_slope
    sll_real64, dimension(:), intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:), intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:), intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:), intent(in), optional :: eta2_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interpolator,ierr)
    
    call interpolator%initialize( &
      npts1,                      &
      npts2,                      &
      eta1_min,                   &
      eta1_max,                   &
      eta2_min,                   &
      eta2_max,                   &
      spline_degree1,             &
      spline_degree2,             &
      eta1_bc_type,               &
      eta2_bc_type,               &
      const_eta1_min_slope,       &
      const_eta1_max_slope,       &
      const_eta2_min_slope,       &
      const_eta2_max_slope,       &
      eta1_min_slopes,            &
      eta1_max_slopes,            &
      eta2_min_slopes,            &
      eta2_max_slopes )
     
  end function sll_f_new_bspline_interpolator_2d

  subroutine initialize_bs2d_interpolator( &
    interpolator,                          &
    npts1,                                 &
    npts2,                                 &
    eta1_min,                              &
    eta1_max,                              &
    eta2_min,                              &
    eta2_max,                              &
    spline_degree1,                        &
    spline_degree2,                        &
    eta1_bc_type,                          &
    eta2_bc_type,                          &
    const_eta1_min_slope,                  &
    const_eta1_max_slope,                  &
    const_eta2_min_slope,                  &
    const_eta2_max_slope,                  &
    eta1_min_slopes,                       &
    eta1_max_slopes,                       &
    eta2_min_slopes,                       &
    eta2_max_slopes )

    class(sll_t_bspline_interpolator_2d), intent(inout) :: interpolator
    sll_int32,  intent(in)                            :: npts1
    sll_int32,  intent(in)                            :: npts2
    sll_real64, intent(in)                            :: eta1_min
    sll_real64, intent(in)                            :: eta1_max
    sll_real64, intent(in)                            :: eta2_min
    sll_real64, intent(in)                            :: eta2_max
    sll_int32,  intent(in)                            :: spline_degree1
    sll_int32,  intent(in)                            :: spline_degree2
    sll_int32,  intent(in),              optional     :: eta1_bc_type
    sll_int32,  intent(in),              optional     :: eta2_bc_type
    sll_real64, intent(in),              optional     :: const_eta1_min_slope
    sll_real64, intent(in),              optional     :: const_eta1_max_slope
    sll_real64, intent(in),              optional     :: const_eta2_min_slope
    sll_real64, intent(in),              optional     :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional     :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional     :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional     :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional     :: eta2_max_slopes

    interpolator%npts1    = npts1
    interpolator%npts2    = npts2
    interpolator%bc_type1 = eta1_bc_type
    interpolator%bc_type2 = eta2_bc_type

    interpolator%spline => sll_f_new_bspline_2d( &
         npts1,                            &
         spline_degree1,                   &
         eta1_min,                         &
         eta1_max,                         &
         eta1_bc_type,                     &
         npts2,                            &
         spline_degree2,                   &
         eta2_min,                         &
         eta2_max,                         &
         eta2_bc_type,                     &
         const_eta1_min_slope,             &
         const_eta1_max_slope,             &
         const_eta2_min_slope,             &
         const_eta2_max_slope              &
    )

    return
    SLL_ASSERT(present(eta1_min_slopes))
    SLL_ASSERT(present(eta1_max_slopes))
    SLL_ASSERT(present(eta2_min_slopes))
    SLL_ASSERT(present(eta2_max_slopes))

  end subroutine initialize_bs2d_interpolator

  subroutine compute_interpolants_bs2d( &
    interpolator,                       &
    data_array,                         &
    eta1_coords,                        &
    size_eta1_coords,                   &
    eta2_coords,                        &
    size_eta2_coords )

    class(sll_t_bspline_interpolator_2d), intent(inout)        :: interpolator
    sll_real64, dimension(:,:),         intent(in)           :: data_array
    sll_real64, dimension(:),           intent(in), optional :: eta1_coords
    sll_real64, dimension(:),           intent(in), optional :: eta2_coords
    sll_int32,                          intent(in), optional :: size_eta1_coords
    sll_int32,                          intent(in), optional :: size_eta2_coords

    if (present(eta1_coords))      print*, '#Warning eta1_coords not used'
    if (present(eta2_coords))      print*, '#Warning eta2_coords not used'
    if (present(size_eta1_coords)) print*, '#Warning size_eta1_coords not used'
    if (present(size_eta2_coords)) print*, '#Warning size_eta2_coords not used'

    call sll_o_compute_bspline_2d( interpolator%spline, data_array )

  end subroutine

  function interpolate_value_bs2d( interpolator, eta1, eta2 ) result(val)

    class(sll_t_bspline_interpolator_2d), intent(in) :: interpolator
    sll_real64,                         intent(in) :: eta1
    sll_real64,                         intent(in) :: eta2

    sll_real64 :: val

    val = sll_f_interpolate_value_2d( interpolator%spline, eta1, eta2, 0, 0 )

  end function

  function interpolate_deriv1_bs2d( interpolator, eta1, eta2 ) result(val)

    class(sll_t_bspline_interpolator_2d), intent(in) :: interpolator
    sll_real64,                         intent(in) :: eta1
    sll_real64,                         intent(in) :: eta2

    sll_real64 :: val

    val = sll_f_interpolate_value_2d( interpolator%spline, eta1, eta2, 1, 0 )

  end function

  function interpolate_deriv2_bs2d( interpolator, eta1, eta2 ) result(val)

    class(sll_t_bspline_interpolator_2d), intent(in) :: interpolator
    sll_real64,                         intent(in) :: eta1
    sll_real64,                         intent(in) :: eta2

    sll_real64 :: val

    val = sll_f_interpolate_value_2d( interpolator%spline, eta1, eta2, 0, 1 )

  end function

  subroutine spline_interpolate2d(this,              &
                                num_points1,       &
                                num_points2,       &
                                data_in,           &
                                eta1,              &
                                eta2,              &
                                data_out) 

    class(sll_t_bspline_interpolator_2d),  intent(in) :: this
    sll_int32,  intent(in)                          :: num_points1
    sll_int32,  intent(in)                          :: num_points2
    sll_real64, dimension(:,:),          intent(in) :: eta1
    sll_real64, dimension(:,:),          intent(in) :: eta2
    sll_real64, dimension(:,:),          intent(in) :: data_in
    sll_real64,                          intent(out):: data_out(num_points1,num_points2)

    sll_int32 :: i
    sll_int32 :: j

    call sll_o_compute_bspline_2d( this%spline, data_in )

    do j = 1, num_points2
      do i = 1, num_points1
        data_out(i,j) = this%interpolate_from_interpolant_value(eta1(i,j),eta2(i,j))
      end do
    end do

  end subroutine spline_interpolate2d

  subroutine spline_interpolate2d_disp(this,        &
                                     num_points1, &
                                     num_points2, &
                                     data_in,     &
                                     alpha1,      &
                                     alpha2,      &
                                     data_out)

    class(sll_t_bspline_interpolator_2d),  intent(in) :: this

    sll_int32,  intent(in)                         :: num_points1
    sll_int32,  intent(in)                         :: num_points2
    sll_real64, dimension(:,:), intent(in)         :: alpha1
    sll_real64, dimension(:,:), intent(in)         :: alpha2
    sll_real64, dimension(:,:), intent(in)         :: data_in
    sll_real64,                 intent(out)        :: data_out(num_points1,num_points2)
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

!PN    eta1_min   = get_x1_min( this%spline ) !this%spline%x1_min 
!PN    eta1_max   = get_x1_max( this%spline ) !this%spline%x1_max 
!PN    eta2_min   = get_x1_min( this%spline ) !this%spline%x2_min 
!PN    eta2_max   = get_x2_max( this%spline ) !this%spline%x2_max 
!PN    delta_eta1 = get_x1_delta( this%spline ) !this%spline%x1_delta  
!PN    delta_eta2 = get_x2_delta( this%spline ) !this%spline%x2_delta  
    
    call sll_o_compute_bspline_2d( this%spline, data_in )

    if(this%bc_type1 == sll_p_periodic .and. &
       this%bc_type2 == sll_p_periodic ) then
       
       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1
             eta2 = eta2_min + (j-1)*delta_eta2
             eta1 = eta1_min + &
                  modulo(eta1-eta1_min+alpha1(i,j),eta1_max-eta1_min)
             eta2 = eta2_min + &
                  modulo(eta2-eta2_min+alpha2(i,j),eta2_max-eta2_min)
             data_out(i,j) = this%interpolate_from_interpolant_value(eta1,eta2)
          end do
       end do

    else if(this%bc_type1 == sll_p_hermite .and. &
            this%bc_type2 == sll_p_hermite ) then
       
       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1 + alpha1(i,j)
             eta2 = eta2_min + (j-1)*delta_eta2 + alpha2(i,j)
             eta1 = min(eta1,eta1_max)
             eta2 = min(eta2,eta2_max)
             eta1 = max(eta1,eta1_min)
             eta2 = max(eta2,eta2_min)
             data_out(i,j) = this%interpolate_from_interpolant_value(eta1,eta2)
          end do
       end do
      
       
    else

       do j = 1, num_points2
          do i = 1, num_points1
             eta1 = eta1_min + (i-1)*delta_eta1 + alpha1(i,j)
             eta2 = eta2_min + (j-1)*delta_eta2 + alpha2(i,j)
             SLL_ASSERT(eta1_min <= eta1 .and. eta1 <= eta1_max)
             SLL_ASSERT(eta2_min <= eta2 .and. eta2 <= eta2_max)
             data_out(i,j) = this%interpolate_from_interpolant_value(eta1,eta2)
          end do
       end do
    end if
  end subroutine spline_interpolate2d_disp

  subroutine set_coefficients_bs2d(  &
       interpolator,                 &
       coeffs_1d,                    &
       coeffs_2d,                    &
       coeff2d_size1,                &
       coeff2d_size2,                &
       knots1,                       &
       size_knots1,                  &
       knots2,                       &
       size_knots2)

    class(sll_t_bspline_interpolator_2d), intent(inout)        :: interpolator
    sll_real64, dimension(:),           intent(in), optional :: coeffs_1d
    sll_real64, dimension(:,:),         intent(in), optional :: coeffs_2d
    sll_int32,                          intent(in), optional :: coeff2d_size1
    sll_int32,                          intent(in), optional :: coeff2d_size2
    sll_real64, dimension(:),           intent(in), optional :: knots1
    sll_real64, dimension(:),           intent(in), optional :: knots2
    sll_int32,                          intent(in), optional :: size_knots1
    sll_int32,                          intent(in), optional :: size_knots2

    print*, 'set_coefficients_bs2d(): This function has not been implemented yet.'
    print*, interpolator%npts1

    if(present(coeffs_1d))     print*,'coeffs_1d present but not used'
    if(present(coeffs_2d))     print*,'coeffs_2d present but not used'
    if(present(coeff2d_size1)) print*,'coeff2d_size1 present but not used'
    if(present(coeff2d_size2)) print*,'coeff2d_size2 present but not used'
    if(present(knots1))        print*,'knots1 present but not used'
    if(present(knots2))        print*,'knots2 present but not used'
    if(present(size_knots1))   print*,'size_knots1 present but not used'
    if(present(size_knots2))   print*,'size_knots2 present but not used'

    stop

  end subroutine !set_coefficients_bs2d
  
!!$  subroutine compute_spl_coeff_bs2d(interpolator, &
!!$       data_array, &
!!$       eta1_coords, &
!!$       size_eta1_coords, &
!!$       eta2_coords, &
!!$       size_eta2_coords )
!!$    class(sll_t_bspline_interpolator_2d), intent(inout)  :: interpolator
!!$    sll_real64, dimension(:,:), intent(in)     :: data_array
!!$    sll_real64, dimension(:), intent(in),optional       :: eta1_coords
!!$    sll_real64, dimension(:), intent(in),optional       :: eta2_coords
!!$    sll_int32, intent(in), optional                     :: size_eta1_coords
!!$    sll_int32, intent(in),optional                      :: size_eta2_coords
!!$
!!$    print *, 'compute_coefficients_bs2d(): ERROR: This function has not been',&
!!$         'implemented yet.'
!!$    stop
!!$  end subroutine compute_spl_coeff_bs2d
  

  function get_coefficients_bs2d(interpolator)
    class(sll_t_bspline_interpolator_2d), intent(in)    :: interpolator
    sll_real64, dimension(:,:), pointer            :: get_coefficients_bs2d     
    
    print *, 'get_coefficients_bs2d(): ERROR: This function has not been ', &
         'implemented yet.'
    get_coefficients_bs2d => null()
    print *,interpolator%npts1    
    stop      
  end function get_coefficients_bs2d

  function coefficients_are_set_bs2d( interpolator ) result(res)
    class(sll_t_bspline_interpolator_2d), intent(in) :: interpolator
    logical :: res
    res = .false.
    print *, 'coefficients_are_set_bs2d(): this function has not been implemented yet.'
    print *,'#',interpolator%npts1
    !stop
  end function coefficients_are_set_bs2d

end module sll_m_bspline_interpolator_2d
