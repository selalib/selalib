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
!> Class for the nufft inmplementation of sll_c_interpolator_2d
!> @details
!> Implements the sll_c_interpolator_2d interface
module sll_m_nufft_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_nufft_interpolation, only:        & 
 sll_t_nufft_2d,                            &
 sll_s_nufft_2d_compute_fft,                &
 sll_s_nufft_2d_init,                       &
 sll_s_nufft_2d_free,                       &
 sll_f_nufft_2d_interpolate_value_from_fft, &
 sll_s_nufft_2d_interpolate_array_values

use sll_m_interpolators_2d_base, only: sll_c_interpolator_2d

implicit none

public :: &
  sll_s_nufft_interpolator_2d_init, &
  sll_s_nufft_interpolator_2d_free, &
  sll_t_nufft_interpolator_2d

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> @brief
!> The nufft-based interpolator is only a wrapper around the capabilities
!> of the nufft package. 
!> @details
!> All interpolators share a common interface with
!> respect to their use, as described by the interpolator_2d_base class.
!> Where the diverse interpolators diverge is in the way to initialize them.
type, extends(sll_c_interpolator_2d) :: sll_t_nufft_interpolator_2d

  sll_int32                      :: num_cells1    !< Number of cells along x
  sll_int32                      :: num_cells2    !< Number of cells along y 
  type(sll_t_nufft_2d)           :: nufft         !< The nufft object
  sll_real64                     :: eta1_min
  sll_real64                     :: eta1_max
  sll_real64                     :: eta2_min
  sll_real64                     :: eta2_max

contains

  !> Allocate data, set dimensions and boundary conditions
  procedure, pass(interpolator) :: initialize => sll_s_nufft_interpolator_2d_init
  !> Compute fft 
  procedure :: compute_interpolants => compute_interpolants_nufft2d
  !> Interpolate single value from last interpolants computed
  procedure :: interpolate_from_interpolant_value => interpolate_value_nufft2d
  !> Interpolate first derivative from last inteprolants computed
  procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_nufft2d
  !> Interpolate first derivative from last inteprolants computed
  procedure :: interpolate_from_interpolant_derivative_eta2 => interpolate_deriv2_nufft2d
  !> Interpolate a 2d array
  procedure, pass :: interpolate_array => nufft_interpolate2d
  !> Interpolate a 2d array after a displacement of original points
  procedure, pass :: interpolate_array_disp => nufft_interpolate2d_disp
  !> Set coefficients (non relevant)
  procedure, pass :: set_coefficients => set_coefficients_nufft2d
  !> Get coefficients (non relevant)
  procedure, pass :: get_coefficients => get_coefficients_nufft2d
  !> Return true if the coefficients are set (non relevant)
  procedure, pass :: coefficients_are_set => coefficients_are_set_nufft2d
  !> Free the memory
  procedure, pass :: delete => sll_s_nufft_interpolator_2d_free

end type sll_t_nufft_interpolator_2d

!> Pointer to this interpolator derived type
type :: sll_nufft_interpolator_2d_ptr
  type(sll_t_nufft_interpolator_2d), pointer :: interp
end type sll_nufft_interpolator_2d_ptr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_nufft_interpolator_2d_free( interpolator )
class(sll_t_nufft_interpolator_2d), intent(inout) :: interpolator

call sll_s_nufft_2d_free(interpolator%nufft)
    
end subroutine sll_s_nufft_interpolator_2d_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
subroutine sll_s_nufft_interpolator_2d_init( &
  interpolator,                              &
  npts1,                                     &
  npts2,                                     &
  eta1_min,                                  &
  eta1_max,                                  &
  eta2_min,                                  &
  eta2_max                                   )

class(sll_t_nufft_interpolator_2d), intent(inout) :: interpolator
sll_int32,  intent(in)                            :: npts1
sll_int32,  intent(in)                            :: npts2
sll_real64, intent(in)                            :: eta1_min
sll_real64, intent(in)                            :: eta1_max
sll_real64, intent(in)                            :: eta2_min
sll_real64, intent(in)                            :: eta2_max

interpolator%num_cells1 = npts1-1
interpolator%num_cells2 = npts2-1
interpolator%eta1_min   = eta1_min
interpolator%eta1_max   = eta1_max
interpolator%eta2_min   = eta2_min
interpolator%eta2_max   = eta2_max

call sll_s_nufft_2d_init(              &
     interpolator%nufft,               &
     npts1-1,                          &
     eta1_min,                         &
     eta1_max,                         &
     npts2-1,                          &
     eta2_min,                         &
     eta2_max)

end subroutine sll_s_nufft_interpolator_2d_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_interpolants_nufft2d( &
     interpolator,                       &
     data_array,                         &
     eta1_coords,                        &
     size_eta1_coords,                   &
     eta2_coords,                        &
     size_eta2_coords )

class(sll_t_nufft_interpolator_2d), intent(inout)        :: interpolator
sll_real64, dimension(:,:),         intent(in)           :: data_array
sll_real64, dimension(:),           intent(in), optional :: eta1_coords
sll_real64, dimension(:),           intent(in), optional :: eta2_coords
sll_int32,                          intent(in), optional :: size_eta1_coords
sll_int32,                          intent(in), optional :: size_eta2_coords
sll_int32 :: nc1, nc2

nc1 = interpolator%num_cells1
nc2 = interpolator%num_cells2

call sll_s_nufft_2d_compute_fft( interpolator%nufft, &
                                 data_array(1:nc1,1:nc2) )

end subroutine compute_interpolants_nufft2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function interpolate_value_nufft2d( interpolator, eta1, eta2 ) result(val)

class(sll_t_nufft_interpolator_2d), intent(in) :: interpolator
sll_real64,                         intent(in) :: eta1
sll_real64,                         intent(in) :: eta2

sll_real64 :: val 

val = sll_f_nufft_2d_interpolate_value_from_fft( interpolator%nufft, eta1, eta2)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function interpolate_deriv1_nufft2d( interpolator, eta1, eta2 ) result(val)

class(sll_t_nufft_interpolator_2d), intent(in) :: interpolator
sll_real64,                         intent(in) :: eta1
sll_real64,                         intent(in) :: eta2

sll_real64 :: val 

val = 0.0_f64
SLL_ERROR( 'interpolate_deriv1_nufft2d', 'not implemented' )

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function interpolate_deriv2_nufft2d( interpolator, eta1, eta2 ) result(val)

class(sll_t_nufft_interpolator_2d), intent(in) :: interpolator
sll_real64,                         intent(in) :: eta1
sll_real64,                         intent(in) :: eta2

sll_real64 :: val

val = 0.0_f64
SLL_ERROR( 'interpolate_deriv2_nufft2d', 'not implemented' )

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nufft_interpolate2d(this,              &
                               num_points1,       &
                               num_points2,       &
                               data_in,           &
                               eta1,              &
                               eta2,              &
                               data_out) 

class(sll_t_nufft_interpolator_2d),  intent(in) :: this
sll_int32,  intent(in)                          :: num_points1
sll_int32,  intent(in)                          :: num_points2
sll_real64, dimension(:,:),          intent(in) :: eta1
sll_real64, dimension(:,:),          intent(in) :: eta2
sll_real64, dimension(:,:),          intent(in) :: data_in
sll_real64,                          intent(out):: data_out(num_points1,num_points2)

sll_int32 :: nc1, nc2

nc1 = this%num_cells1
nc2 = this%num_cells2

call sll_s_nufft_2d_interpolate_array_values( this%nufft,           &
                                              data_in(1:nc1,1:nc2), &
                                              eta1(1:nc1,1:nc2),    &
                                              eta2(1:nc1,1:nc2),    &
                                              data_out(1:nc1,1:nc2) )

data_out(nc1+1,:) = data_out(1,:)
data_out(:,nc2+1) = data_out(:,1)

end subroutine nufft_interpolate2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nufft_interpolate2d_disp(this,        &
                                    num_points1, &
                                    num_points2, &
                                    data_in,     &
                                    alpha1,      &
                                    alpha2,      &
                                    data_out)

class(sll_t_nufft_interpolator_2d),  intent(in) :: this

sll_int32,  intent(in)                         :: num_points1
sll_int32,  intent(in)                         :: num_points2
sll_real64, dimension(:,:), intent(in)         :: alpha1
sll_real64, dimension(:,:), intent(in)         :: alpha2
sll_real64, dimension(:,:), intent(in)         :: data_in
sll_real64,                 intent(out)        :: data_out(num_points1,num_points2)

sll_real64, dimension(:,:), allocatable        :: eta1
sll_real64, dimension(:,:), allocatable        :: eta2

sll_real64                                     :: eta1_min
sll_real64                                     :: eta1_max
sll_real64                                     :: delta_eta1
sll_real64                                     :: eta2_min
sll_real64                                     :: eta2_max
sll_real64                                     :: delta_eta2
sll_int32                                      :: i
sll_int32                                      :: j
sll_int32                                      :: error
sll_int32                                      :: nc1
sll_int32                                      :: nc2

nc1 = this%num_cells1
nc2 = this%num_cells2

eta1_min   = this%eta1_min
eta1_max   = this%eta1_max
eta2_min   = this%eta2_min
eta2_max   = this%eta2_max
delta_eta1 = (this%eta1_max-this%eta1_min)/real(this%num_cells1,f64)
delta_eta2 = (this%eta2_max-this%eta2_min)/real(this%num_cells2,f64)

SLL_ALLOCATE(eta1(1:nc1,1:nc2), error)
SLL_ALLOCATE(eta2(1:nc1,1:nc2), error)
do j = 1, nc2
  do i = 1, nc1
    eta1(i,j) = (i-1)*delta_eta1
    eta2(i,j) = (j-1)*delta_eta2
    eta1(i,j) = eta1_min + modulo(eta1(i,j)+alpha1(i,j),eta1_max-eta1_min)
    eta2(i,j) = eta2_min + modulo(eta2(i,j)+alpha2(i,j),eta2_max-eta2_min)
  end do
end do

call sll_s_nufft_2d_interpolate_array_values( this%nufft, &
  data_in(1:nc1,1:nc2), eta1, eta2, data_out(1:nc1,1:nc2) )

data_out(nc1+1,:) = data_out(1,:)
data_out(:,nc2+1) = data_out(:,1)

end subroutine nufft_interpolate2d_disp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_coefficients_nufft2d(  &
        interpolator,                 &
        coeffs_1d,                    &
        coeffs_2d,                    &
        coeff2d_size1,                &
        coeff2d_size2,                &
        knots1,                       &
        size_knots1,                  &
        knots2,                       &
        size_knots2)

class(sll_t_nufft_interpolator_2d), intent(inout)          :: interpolator
sll_real64, dimension(:),           intent(in),   optional :: coeffs_1d
sll_real64, dimension(:,:),         intent(in),   optional :: coeffs_2d
sll_int32,                          intent(in),   optional :: coeff2d_size1
sll_int32,                          intent(in),   optional :: coeff2d_size2
sll_real64, dimension(:),           intent(in),   optional :: knots1
sll_real64, dimension(:),           intent(in),   optional :: knots2
sll_int32,                          intent(in),   optional :: size_knots1
sll_int32,                          intent(in),   optional :: size_knots2

print*, interpolator%num_cells1, interpolator%num_cells2

if(present(coeffs_1d))     print*,'coeffs_1d present but not used'
if(present(coeffs_2d))     print*,'coeffs_2d present but not used'
if(present(coeff2d_size1)) print*,'coeff2d_size1 present but not used'
if(present(coeff2d_size2)) print*,'coeff2d_size2 present but not used'
if(present(knots1))        print*,'knots1 present but not used'
if(present(knots2))        print*,'knots2 present but not used'
if(present(size_knots1))   print*,'size_knots1 present but not used'
if(present(size_knots2))   print*,'size_knots2 present but not used'

SLL_ERROR('set_coefficients','Not implemented for nufft2d')

end subroutine !set_coefficients_nufft2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_coefficients_nufft2d(interpolator) result(res)
class(sll_t_nufft_interpolator_2d), intent(in)    :: interpolator
sll_real64, pointer :: res(:,:)

res = 0.0_f64
   
SLL_ERROR('get_coefficients','Not implemented for nufft2d')

end function get_coefficients_nufft2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function coefficients_are_set_nufft2d( interpolator ) result(res) 
class(sll_t_nufft_interpolator_2d), intent(in) :: interpolator

logical :: res

res = .false.

SLL_ERROR('coefficients_are_set','Not implemented for nufft2d')

end function coefficients_are_set_nufft2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module sll_m_nufft_interpolator_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
