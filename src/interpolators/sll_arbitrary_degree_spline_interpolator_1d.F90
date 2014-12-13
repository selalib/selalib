!**************************************************************
!  Copyright INRIA
!  Authors :
!     Aurore
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

!> Class interpolator and methods for arbitrary degree spline 1D interpolator
module sll_module_arbitrary_degree_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

use sll_module_interpolators_1d_base

implicit none
private

type deboor_type
   sll_int32 :: ilo = 1
   sll_int32 :: j   = 1
   sll_real64, dimension ( 20 ) :: deltal
   sll_real64, dimension ( 20 ) :: deltar
end type deboor_type

!> Class for arbitrary degree spline 1d interpolator
type, public, extends(sll_interpolator_1d_base) :: sll_arbitrary_degree_spline_interpolator_1d
  sll_int32                         :: num_pts
  sll_real64                        :: eta_min
  sll_real64                        :: eta_max
  sll_int32                         :: bc_left
  sll_int32                         :: bc_right
  sll_int32                         :: spline_degree
  sll_real64, dimension(:), pointer :: eta
  sll_real64, dimension(:), pointer :: t
  sll_real64, dimension(:), pointer :: work
  sll_int32                         :: size_t
  sll_int64                         :: bc_selector 
  sll_real64, dimension(:), pointer :: coeff_splines
  sll_int32                         :: size_coeffs
  sll_real64                        :: slope_left
  sll_real64                        :: slope_right
  sll_real64                        :: value_left
  sll_real64                        :: value_right
  logical                           :: compute_slope_left = .TRUE.
  logical                           :: compute_slope_right= .TRUE.
  logical                           :: compute_value_left = .TRUE.
  logical                           :: compute_value_right= .TRUE.
  type(deboor_type)                 :: deboor
contains
  procedure :: initialize=>initialize_ad1d_interpolator
  procedure :: set_coefficients => set_coefficients_ad1d
  procedure :: compute_interpolants => compute_interpolants_ad1d
  procedure :: interpolate_value => interpolate_value_ad1d
  procedure :: interpolate_array_values => interpolate_values_ad1d
  procedure :: interpolate_pointer_values => interpolate_pointer_values_ad1d
  procedure :: interpolate_derivative_eta1 => interpolate_derivative_ad1d
  procedure :: interpolate_array_derivatives => interpolate_derivatives_ad1d
  procedure :: interpolate_pointer_derivatives =>interpolate_pointer_derivatives_ad1d
  procedure :: interpolate_array => interpolate_array_ad1d
  procedure :: interpolate_array_disp => interpolate_1d_array_disp_ad1d
  procedure :: get_coefficients => get_coefficients_ad1d
  procedure :: reconstruct_array
  procedure :: delete => delete_arbitrary_degree_1d_interpolator
end type sll_arbitrary_degree_spline_interpolator_1d

!> Deallocate
interface sll_delete
   module procedure delete_arbitrary_degree_1d_interpolator
end interface sll_delete

public sll_delete 
public new_arbitrary_degree_1d_interpolator
public set_values_at_boundary1d
public initialize_ad1d_interpolator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief delete interpolator arbitrary degree splines.
!> @details
!>
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d

subroutine delete_arbitrary_degree_1d_interpolator( interpolator )
  class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator
  sll_int32 :: ierr
  SLL_DEALLOCATE(interpolator%t,ierr)
  SLL_DEALLOCATE(interpolator%coeff_splines,ierr)
end subroutine delete_arbitrary_degree_1d_interpolator


!> @brief Initialization of a pointer interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d such as a pointer
!>
!> The parameters are
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @return the type interpolator arbitrary degree splines 1d

function new_arbitrary_degree_1d_interpolator( num_pts,       &
                                               eta_min,       &
                                               eta_max,       &
                                               bc_left,       &
                                               bc_right,      &
                                               spline_degree) &
result(interpolator)

class(sll_arbitrary_degree_spline_interpolator_1d),pointer :: interpolator
sll_int32,  intent(in) :: num_pts
sll_real64, intent(in) :: eta_min
sll_real64, intent(in) :: eta_max
sll_int32,  intent(in) :: bc_left
sll_int32,  intent(in) :: bc_right
sll_int32,  intent(in) :: spline_degree

sll_int32              :: ierr

SLL_ALLOCATE(interpolator,ierr)

call initialize_ad1d_interpolator( interpolator, &
                                   num_pts,      &
                                   eta_min,      &
                                   eta_max,      &
                                   bc_left,      &
                                   bc_right,     &
                                   spline_degree)

end function new_arbitrary_degree_1d_interpolator

!> @brief Initialization of interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d
!>
!> The parameters are
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d

subroutine initialize_ad1d_interpolator( interpolator, &
                                         num_pts,      &
                                         eta_min,      &
                                         eta_max,      &
                                         bc_left,      &
                                         bc_right,     &
                                         spline_degree)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator
sll_int32, intent(in) :: num_pts
sll_real64, intent(in) :: eta_min
sll_real64, intent(in) :: eta_max
sll_int32, intent(in) :: bc_left
sll_int32, intent(in) :: bc_right
sll_int32, intent(in) :: spline_degree
sll_int32 :: ierr
sll_int32 :: tmp
sll_int64 :: bc_selector
sll_int32 :: i, k
sll_real64 :: delta_eta

! do some argument checking...
if(((bc_left == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC))) then
   print *, 'initialize_arbitrary_degree_1d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in first direction.'
end if

bc_selector = 0

if( bc_left  == SLL_DIRICHLET ) bc_selector = bc_selector + 1
if( bc_left  == SLL_NEUMANN   ) bc_selector = bc_selector + 2
if( bc_left  == SLL_HERMITE   ) bc_selector = bc_selector + 4
if( bc_right == SLL_DIRICHLET ) bc_selector = bc_selector + 8
if( bc_right == SLL_NEUMANN   ) bc_selector = bc_selector + 16
if( bc_right == SLL_HERMITE   ) bc_selector = bc_selector + 32

interpolator%spline_degree = spline_degree
interpolator%eta_min       = eta_min
interpolator%eta_max       = eta_max
interpolator%bc_left       = bc_left
interpolator%bc_right      = bc_right
interpolator%bc_selector   = bc_selector
interpolator%num_pts       = num_pts

delta_eta = (eta_max - eta_min)/(num_pts-1)
SLL_ALLOCATE(interpolator%eta(num_pts),ierr)
do i = 1, num_pts
   interpolator%eta(i) = eta_min + delta_eta*(i-1)
end do

k = spline_degree+1

select case (bc_selector)

case (0) ! 1. periodic

  SLL_ALLOCATE(interpolator%coeff_splines(num_pts),ierr)
  SLL_ALLOCATE(interpolator%t(num_pts+k),ierr)
  interpolator%t(1:k) = eta_min
  interpolator%t(num_pts+1:num_pts+k) = eta_max
  
  if ( mod(k,2) == 0 ) then
    do i = k+1,num_pts
      interpolator%t(i) = interpolator%eta(i-k/2) 
    end do
  else
    do i = k+1, num_pts
      interpolator%t(i) = 0.5*(interpolator%eta(i-(k-1)/2)+    &
                               interpolator%eta(i-1-(k-1)/2))
    end do
  end if
  SLL_ALLOCATE(interpolator%work(num_pts*(2*k-1)),ierr)

case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%value_left  = 0.0_f64
  interpolator%value_right = 0.0_f64

  SLL_ALLOCATE(interpolator%coeff_splines(num_pts),ierr)
  SLL_ALLOCATE(interpolator%t(num_pts+k),ierr)
  interpolator%t(1:k) = eta_min
  interpolator%t(num_pts+1:num_pts+k) = eta_max
  
  if ( mod(k,2) == 0 ) then
    do i = k+1,num_pts
      interpolator%t(i) = interpolator%eta(i-k/2) 
    end do
  else
    do i = k+1, num_pts
      interpolator%t(i) = 0.5*(interpolator%eta(i-(k-1)/2)+    &
                               interpolator%eta(i-1-(k-1)/2))
    end do
  end if
  SLL_ALLOCATE(interpolator%work(num_pts*(2*k-1)),ierr)

case(10) ! Neumann - Dirichlet
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_left  = 0.0_f64
   interpolator%value_right = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(12) ! Hermite- Dirichlet
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_left = 0.0_f64
   interpolator%value_right = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(17) ! Dirichlet-Neumann
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%value_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(18) ! Neumann - Neumann
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%slope_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(20) ! Hermite - Neumann
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%slope_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(33) ! Dirichlet - Hermite
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%value_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(34) ! Neumann- Hermite
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%slope_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case(36) ! Hermitte - Hermite
   tmp = num_pts*num_pts
   SLL_ALLOCATE( interpolator%coeff_splines(tmp),ierr)
   interpolator%slope_right = 0.0_f64
   interpolator%slope_left = 0.0_f64
   SLL_ALLOCATE( interpolator%t(num_pts*num_pts),ierr)
   interpolator%t(:) = 0.0_f64

case default
   print *, 'initialize_ad1d_interpolator: BC combination not implemented.'
end select

interpolator%coeff_splines(:) = 0.0_f64
end subroutine initialize_ad1d_interpolator


!> Initialization of the boundary for interpolator arbitrary degree splines 1d.
!> The parameters are
!> @param[in]  value_left  contains the value in the left
!> @param[in]  value_right contains the value in the right
!> @param[in]  slope_left  contains the value in the left for derivative
!> @param[in]  slope_right contains the value in the right for derivative
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_values_at_boundary1d( interpolator, &
                                     value_left,   &
                                     value_right,  &
                                     slope_left,   &
                                     slope_right)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout) :: interpolator

sll_real64, intent(in), optional :: value_left
sll_real64, intent(in), optional :: value_right
sll_real64, intent(in), optional :: slope_left
sll_real64, intent(in), optional :: slope_right

sll_int32 :: bc_left
sll_int32 :: bc_right
sll_int64 :: bc_selector

bc_left = interpolator%bc_left
bc_right= interpolator%bc_right
bc_selector = interpolator%bc_selector

if (present(value_left)) then
  interpolator%value_left = value_left
  interpolator%compute_value_left = .FALSE.
end if

if (present(value_right)) then
  interpolator%value_right = value_right
  interpolator%compute_value_right = .FALSE.
end if

if (present(slope_left)) then
  interpolator%slope_left = slope_left
  interpolator%compute_slope_left = .FALSE.
end if

if (present(slope_right)) then
  interpolator%slope_right = slope_right
  interpolator%compute_slope_right = .FALSE.
end if

end subroutine set_values_at_boundary1d

!> @brief computing the coefficients spline with a given
!>  data_array 1D cooresponding at the values of a function
!> @details 
!>  on eta_coords of size size_eta_coords
!>  if the eta_coords and eta_coords is not given
!>  we consider that the values of the function is on the points in the mesh_1d
!>
!> The parameters are
!> @param[in]  data_array the 1d arrays corresponding at the values of a function
!> @param[in]  eta_coords the 1d arrays
!> @param[in]  size_eta_coords the size of eta_coords
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine compute_interpolants_ad1d( interpolator,    &
                                      data_array,      &
                                      eta_coords,      &
                                      size_eta_coords)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout)  :: interpolator

sll_real64, dimension(:), intent(in)          :: data_array
sll_real64, dimension(:),pointer              :: data_array_derivative
sll_real64, dimension(:), intent(in),optional :: eta_coords
sll_int32, intent(in),optional                :: size_eta_coords

sll_real64, dimension(:),pointer              :: point_locate_eta
sll_int32, dimension(:),pointer               :: point_locate_eta_derivative
sll_real64 :: delta_eta
sll_int32  :: sz,sz_deriv
sll_real64 :: period
sll_int32  :: order
sll_int32  :: ierr
sll_int32  :: i
logical    :: user_coords

if(present(eta_coords) .and. (.not. present(size_eta_coords))) then
   print *, 'compute_interpolants_ad1d(), ERROR: if eta_coords is ', &
        'passed, its size must be specified as well through ', &
        'size_eta_coords.'
   SLL_ERROR(' ')
end if

if( present(eta_coords) ) then
  user_coords = .true.
  sz = size(data_array)
  SLL_ALLOCATE(point_locate_eta(sz),ierr)
  point_locate_eta = eta_coords
else
  user_coords = .false.
  sz = interpolator%num_pts
  SLL_ALLOCATE(point_locate_eta(sz),ierr)
  point_locate_eta = interpolator%eta
end if

if (interpolator%compute_slope_left .eqv. .true.) then
  interpolator%slope_left = forward_fd_5pt( data_array,point_locate_eta)
end if
if (interpolator%compute_slope_right .eqv. .true.) then
  interpolator%slope_right = backward_fd_5pt( data_array,point_locate_eta,sz)
end if

if (interpolator%compute_value_left .eqv. .true.) then
  interpolator%value_left = data_array(1)
end if
if (interpolator%compute_value_right .eqv. .true.) then
  interpolator%value_right = data_array(sz)
end if

SLL_ASSERT(sz .le. interpolator%num_pts* interpolator%num_pts)
SLL_ASSERT(size(data_array) .le. sz)
SLL_ASSERT(size(point_locate_eta)  .ge. sz)

order  = interpolator%spline_degree + 1
period = interpolator%eta_max - interpolator%eta_min

select case (interpolator%bc_selector)
case (0) ! periodic

  interpolator%size_coeffs = sz !+ 1
  interpolator%size_t = order + sz !+ 1
  interpolator%deboor%ilo = 1

  call splint ( interpolator%deboor,        &
                interpolator%eta,           &
                data_array,                 &
                interpolator%t,             &
                interpolator%num_pts,       &
                order,                      &
                interpolator%work,          &
                interpolator%coeff_splines, &
                ierr)
    
case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%size_coeffs = sz
  interpolator%size_t = order + sz
    
  interpolator%coeff_splines  = data_array

  call splint ( interpolator%deboor,        &
                interpolator%eta,           &
                data_array,                 &
                interpolator%t,             &
                interpolator%num_pts,       &
                order,                      &
                interpolator%work,          &
                interpolator%coeff_splines, &
                ierr)
    
  interpolator%coeff_splines(1)  = interpolator%value_left
  interpolator%coeff_splines(sz) = interpolator%value_right

case(10) ! Neumann - Dirichlet
       ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       ! test dirichlet non homogene
       interpolator%coeff_splines(sz+sz_deriv) = interpolator%value_right

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(12) ! Hermite - Dirichlet

       ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       ! test dirichlet non homogene
       interpolator%coeff_splines(sz+sz_deriv) = interpolator%value_right
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(17) ! Dirichlet - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       ! test dirichlet non homogene
       interpolator%coeff_splines(1) = interpolator%value_left
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
    case(18) ! Neumann - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64
       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(20) ! Hermite - Neumann

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = 0.0_f64

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(33) ! Dirichlet - Hermite

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       ! test dirichlet non homogene
       interpolator%coeff_splines(1) = interpolator%value_left
       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)
    case(34) ! Neumann - Hermite

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------

       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

       SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = 0.0_f64
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    case(36) ! Hermite - Hermite

        ! -----------------------------------
       !!! It is only for cubic spline !!!!
       ! -----------------------------------
       sz_deriv = 2
       interpolator%size_coeffs = sz + sz_deriv
       interpolator%size_t = order + sz + sz_deriv

        SLL_ALLOCATE(point_locate_eta_derivative(sz_deriv),ierr)
       SLL_ALLOCATE(data_array_derivative(sz_deriv),ierr)

       point_locate_eta_derivative(1) = 1
       data_array_derivative(1)       = interpolator%slope_left
       point_locate_eta_derivative(2) = sz
       data_array_derivative(2)       = interpolator%slope_right

       call spli1d_der(sz,sz_deriv,order,&
            point_locate_eta,&
            data_array,&
            point_locate_eta_derivative,&
            data_array_derivative,&
            interpolator%coeff_splines(1:sz+sz_deriv),&
            interpolator%t(1:sz+order+sz_deriv))

       SLL_DEALLOCATE(point_locate_eta_derivative,ierr)
       SLL_DEALLOCATE(data_array_derivative,ierr)

    end select
    SLL_DEALLOCATE(point_locate_eta,ierr)

end subroutine compute_interpolants_ad1d

!> @brief Interpolation on the points eta using
!> the arbitrary degree splines interpolator 1d
!> @details computing the values with the interpolator arbitrary degree splines 1d
!>  on the points eta of arbitrary degree splines 1d
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
!> @param[in] eta1 the point
!> @return val the values on the points eta
function interpolate_value_ad1d( interpolator, eta1) result(val)

class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)  :: interpolator

sll_real64, intent(in)          :: eta1
sll_real64                      :: val
sll_int32                       :: size_coeffs
sll_real64                      :: res
sll_real64,dimension(:),pointer :: knot_tmp
sll_real64,dimension(:),pointer :: coef_tmp

size_coeffs = interpolator%size_coeffs

res = eta1

select case (interpolator%bc_selector)
case (0) ! periodic

  if( res < interpolator%eta_min ) then
     res = res+interpolator%eta_max-interpolator%eta_min
  else if( res >  interpolator%eta_max ) then
     res = res+interpolator%eta_min-interpolator%eta_max
  end if

end select

SLL_ASSERT( res >= interpolator%eta_min )
SLL_ASSERT( res <= interpolator%eta_max )

val = bvalue( interpolator%deboor,                       &
              interpolator%t(1:interpolator%size_t),     &
              interpolator%coeff_splines(1:size_coeffs), &
              size_coeffs,                               &
              interpolator%spline_degree+1,              &
              res,                                       &
              0)

end function interpolate_value_ad1d


  !> @brief First derivative interpolation on the point eta
  !> @details computing the values of the first derivative
  !> with the interpolator arbitrary degree splines 1d
  !> on the points eta of arbitrary degree splines 1d
  !>
  !> The parameters are
  !> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
  !> @param[in] eta1 the point
  !> @return val the values on the point eta of the first derivative

  function interpolate_derivative_ad1d( &
    interpolator, &
    eta1 ) result(val)

    class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)  :: interpolator
    sll_real64, intent(in)         :: eta1
    sll_real64                     :: val
    sll_int32 :: size_coeffs
    sll_real64 :: res
    sll_real64, dimension(:),pointer :: knot_tmp
    sll_real64, dimension(:),pointer :: coef_tmp

    SLL_ASSERT( eta1 .ge. interpolator%eta_min )
    SLL_ASSERT( eta1 .le. interpolator%eta_max )

    size_coeffs = interpolator%size_coeffs

    res = eta1

    select case (interpolator%bc_selector)
    case (0) ! periodic
!
       if( res < interpolator%eta_min ) then
          res = res+interpolator%eta_max-interpolator%eta_min
       else if( res >  interpolator%eta_max ) then
          res = res+interpolator%eta_min-interpolator%eta_max
       end if
    case (9) ! 2. dirichlet-left, dirichlet-right
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(10) ! Neumann - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(12) ! Hermite - Dirichlet
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(17) ! Dirichlet - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(18) ! Neumann - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(20) ! Hermite - Neumann
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(33) ! Dirichlet - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(34) ! Neumann - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if
    case(36) ! Hermite - Hermite
       SLL_ASSERT( res >= interpolator%eta_min )
       SLL_ASSERT( res <= interpolator%eta_max )
       if ( res > interpolator%eta_max) then
          print*, 'problem  x > eta_max'
          stop
       end if
       if ( res < interpolator%eta_min) then
          print*, 'problem  x < eta_min'
          stop
       end if

    end select

    knot_tmp =>  interpolator%t(1:interpolator%size_t)
    coef_tmp => interpolator%coeff_splines(1:size_coeffs)
  
    val = bvalue( &
         interpolator%deboor, &
         knot_tmp,&
         coef_tmp, &
         size_coeffs, &
         interpolator%spline_degree+1, &
         res, &
         1)
  end function interpolate_derivative_ad1d

  function interpolate_array_ad1d( &
       this, num_points, data, coordinates)&
    result(data_out)
    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in)       :: this
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out

    print *, 'interpolate_array_ad1d: not implemented'
    data_out = -1000000._f64*data*coordinates*this%spline_degree
  end function interpolate_array_ad1d

  function interpolate_1d_array_disp_ad1d( &
       this,        &
       num_points, &
       data,     &
       alpha) result(res)

    class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)    :: this
    sll_int32, intent(in)                          :: num_points
    sll_real64, dimension(:), intent(in)         :: data
    sll_real64, intent(in)         :: alpha
    sll_real64, dimension(num_points) :: res

    print *, 'interpolate_1d_array_disp_ad1d: not implemented.'
    res = -1000000._f64*alpha*data*this%spline_degree
  end function interpolate_1d_array_disp_ad1d


  function get_coefficients_ad1d(interpolator)
    class(sll_arbitrary_degree_spline_interpolator_1d), intent(in)    :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_ad1d

    get_coefficients_ad1d => interpolator%coeff_splines
  end function get_coefficients_ad1d

  subroutine interpolate_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )

    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    sll_int32 :: idx

    SLL_ASSERT(num_pts==size(vals_to_interpolate))
    do idx=1,num_pts
          output_array(idx)=interpolate_value_ad1d( &
                                interpolator, &
                                vals_to_interpolate(idx))
    enddo

!    print*, 'interpolate_values_ad1d NOT iMPLEMENTED YET'
!    output_array = -1000000._f64*num_pts&
!         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_values_ad1d

  subroutine interpolate_pointer_values_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )

    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    sll_int32 :: idx

    SLL_ASSERT(num_pts==size(vals_to_interpolate))
    do idx=1,num_pts
          output(idx)=interpolate_value_ad1d( &
                                interpolator, &
                                vals_to_interpolate(idx))
    enddo
!    print*, 'interpolate_pointer_values_ad1d NOT iMPLEMENTED YET'
!    output = -1000000._f64*num_pts&
!         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_pointer_values_ad1d

  subroutine interpolate_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )

    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in) :: vals_to_interpolate
    sll_real64, dimension(:), intent(out) :: output_array
    sll_int32 :: idx

    SLL_ASSERT(num_pts==size(vals_to_interpolate))
    do idx=1,num_pts
          output_array(idx)=interpolate_derivative_ad1d( &
                                interpolator, &
                                vals_to_interpolate(idx))
    enddo
!    print*, 'interpolate_derivatives_ad1d NOT iMPLEMENTED YET'
!    output_array = -1000000._f64*num_pts &
!         *vals_to_interpolate*interpolator%spline_degree
  end subroutine interpolate_derivatives_ad1d

  subroutine interpolate_pointer_derivatives_ad1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )
    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    sll_int32 :: idx

    SLL_ASSERT(num_pts==size(vals_to_interpolate))
    do idx=1,num_pts
          output(idx)=interpolate_derivative_ad1d( &
                                interpolator, &
                                vals_to_interpolate(idx))
    enddo
!    output = -1000000.0_f64*num_pts&
!         *vals_to_interpolate*interpolator%spline_degree
!    print*, 'interpolate_pointer_derivatives_ad1d NOT iMPLEMENTED YET'
  end subroutine interpolate_pointer_derivatives_ad1d

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure

    class(sll_arbitrary_degree_spline_interpolator_1d),  intent(in) :: this
    sll_int32, intent(in)                :: num_points! size of output array
    sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated
    sll_real64, dimension(num_points)    :: res
    res(:) = -1000000.0_f64*data*this%spline_degree
    print*, 'reconstruct_array 1d not implemented yet'
  end function reconstruct_array

  ! The following two functions are wrong, the stencil to compute the
  ! derivatives is valid only for uniform spacing between the points, the
  ! specific coefficients chosen here using the eta coordinates is not to
  ! be used. ABSOLUTE FIXME!!
  function forward_fd_5pt( data,eta) result(res)
    sll_real64, dimension(:), intent(in) :: data
    sll_real64, dimension(:), intent(in) :: eta
    sll_real64 :: res

    res = (-(25.0_f64/12.0_f64)*data(1)*(eta(2) - eta(1)) &
                      + 4.0_f64*data(2)*(eta(3) - eta(2)) &
                      - 3.0_f64*data(3)*(eta(4) - eta(3)) &
            + (4.0_f64/3.0_f64)*data(4)*(eta(5) - eta(4)) &
                      -0.25_f64*data(5)*(eta(6) - eta(5)))
  end function forward_fd_5pt

  function backward_fd_5pt( data,eta,li)result(res)
    sll_real64, dimension(:), intent(in) :: data
    sll_real64, dimension(:), intent(in) :: eta
    sll_int32, intent(in)  :: li  ! last index of the array
    sll_real64 :: res


    res = (0.25_f64*data(li-4)*(eta(li-5) - eta(li-4)) -&
         (4.0_f64/3.0_f64)*  data(li-3)*(eta(li-4) - eta(li-3)) + &
          3.0_f64*           data(li-2)*(eta(li-3) - eta(li-2)) - &
          4.0_f64*           data(li-1)*(eta(li-2) - eta(li-1)) + &
         (25.0_f64/12.0_f64)*data(li)*  (eta(li-1) - eta(li)) )
  end function backward_fd_5pt


subroutine interv( deboor, xt, lxt, x, left, mflag )
    
type(deboor_type)                       :: deboor
sll_real64, dimension(:), intent(in)    :: xt!(lxt)
sll_int32,                intent(in)    :: lxt
sll_real64,               intent(in)    :: x
sll_int32,                intent(out)   :: left
sll_int32,                intent(out)   :: mflag
sll_int32                               :: ihi

sll_int32 :: istep
sll_int32 :: middle

ihi = deboor%ilo + 1
if ( lxt <= ihi ) then
   if ( xt(lxt) <= x ) go to 110
   if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
   end if
   deboor%ilo = lxt - 1
   ihi = lxt
end if
if ( xt(ihi) <= x ) go to 20
if ( xt(deboor%ilo) <= x ) then
   mflag = 0
   left = deboor%ilo
   return
end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
istep = 1
    
10  continue
    
ihi = deboor%ilo
deboor%ilo = ihi - istep

if ( 1 < deboor%ilo ) then
   if ( xt(deboor%ilo) <= x ) go to 50
   istep = istep * 2
   go to 10
end if

deboor%ilo = 1

if ( x < xt(1) ) then
   mflag = -1
   left = 1
   return
end if

go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20  continue
    
istep = 1
    
30  continue
    
deboor%ilo = ihi
ihi = deboor%ilo + istep

if ( ihi < lxt ) then
   if ( x < xt(ihi) ) go to 50
   istep = istep * 2
   go to 30
end if

if ( xt(lxt) <= x ) go to 110
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
ihi = lxt

50  continue
    
do
   
   middle = ( deboor%ilo + ihi ) / 2
   if ( middle == deboor%ilo ) then
      mflag = 0
      left = deboor%ilo
      return
   end if
   !
   !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
   !
   if ( xt(middle) <= x ) then
      deboor%ilo = middle
   else
      ihi = middle
   end if
   
end do
!
!  Set output and return.
!

110 continue
    
mflag = 1

if ( x == xt(lxt) ) mflag = 0

do left = lxt, 1, -1
   if ( xt(left) < xt(lxt) ) return
end do

end subroutine interv


subroutine bsplvb ( deboor, t, jhigh, index, x, left, biatx )

    implicit none
    
    
    type(deboor_type) :: deboor
    sll_int32:: jhigh
    
    sll_real64,dimension(jhigh):: biatx !(jhigh)
    sll_int32:: i
    sll_int32:: index
    sll_int32:: left
    sll_real64:: saved
    sll_real64,dimension(left+jhigh):: t!() left+jhigh
    sll_real64:: term
    sll_real64:: x
    
    if ( index == 1 ) then
       deboor%j = 1
       biatx(1) = 1.0_8
       if ( jhigh <= deboor%j ) then
          return
       end if
    end if
    
    if ( t(left+1) <= t(left) ) then
       print*,'x=',x
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'BSPLVB - Fatal error!'
       write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
       write ( *, '(a,i8)' ) '  But LEFT = ', left
       write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
       write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
       stop
    end if
    
    do
       
       deboor%deltar(deboor%j) = t(left+deboor%j) - x
       deboor%deltal(deboor%j) = x - t(left+1-deboor%j)
       
       saved = 0.0_f64
       do i = 1, deboor%j
          term = biatx(i) / ( deboor%deltar(i) + deboor%deltal(deboor%j+1-i) )
          biatx(i) = saved + deboor%deltar(i) * term
          saved = deboor%deltal(deboor%j+1-i) * term
       end do
    
       biatx(deboor%j+1) = saved
       deboor%j = deboor%j + 1
       
       if ( jhigh <= deboor%j ) then
          
          exit
       end if
    
    end do
    
    return
  end subroutine bsplvb
  


  subroutine bsplvd ( deboor, t, k, x, left, a, dbiatx, nderiv )

    implicit none
    
    sll_int32 :: k
    sll_int32 :: left
    sll_int32 :: nderiv
    
    sll_real64, dimension(k,k):: a!(k,k)
    sll_real64,dimension(k,nderiv), intent(out) :: dbiatx!(k,nderiv)
    sll_real64:: factor
    sll_real64:: fkp1mm
    sll_int32 :: i
    sll_int32 :: ideriv
    sll_int32 :: il
    sll_int32 :: j
    sll_int32 :: jlow
    sll_int32 :: jp1mid
    sll_int32 :: ldummy
    sll_int32 :: m
    sll_int32 :: mhigh
    !  sll_real64 sum1  ! this one is not used...
    sll_real64,dimension(left+k):: t ! (left+k)
    sll_real64:: x
    type(deboor_type) :: deboor
    
    
    mhigh = max ( min ( nderiv, k ), 1 )
    !
    !  MHIGH is usually equal to NDERIV.
    !
    call bsplvb ( deboor, t, k+1-mhigh, 1, x, left, dbiatx )
    
    if ( mhigh == 1 ) then
       return
    end if
    !
    !  The first column of DBIATX always contains the B-spline values
    !  for the current order.  These are stored in column K+1-current
    !  order before BSPLVB is called to put values for the next
    !  higher order on top of it.
    !
    ideriv = mhigh
    do m = 2, mhigh
       jp1mid = 1
       do j = ideriv, k
          dbiatx(j,ideriv) = dbiatx(jp1mid,1)
          jp1mid = jp1mid + 1
          
       end do
       ideriv = ideriv - 1
       
       call bsplvb ( deboor, t, k+1-ideriv, 2, x, left, dbiatx )
       
    end do
    !
    !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
    !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
    !
    !  In particular, the first column of DBIATX is already in final form.
    !
    !  To obtain corresponding derivatives of B-splines in subsequent columns,
    !  generate their B-representation by differencing, then evaluate at X.
    !
    jlow = 1
    do i = 1, k
       a(jlow:k,i) = 0.0D+00
       jlow = i
       a(i,i) = 1.0D+00
    end do
    !
    !  At this point, A(.,J) contains the B-coefficients for the J-th of the
    !  K B-splines of interest here.
    !
    do m = 2, mhigh
       
       fkp1mm = real ( k + 1 - m, kind = 8 )
       il = left
       i = k
       !
       !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
       !  B-splines from those for preceding derivative by differencing
       !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
       !  I < J is used.
       !
       do ldummy = 1, k+1-m
          
          factor = fkp1mm / ( t(il+k+1-m) - t(il) )
          !
          !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
          !  in FACTOR nonzero.
          !
          a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
          
          il = il - 1
          i = i - 1
       
       end do
       !
       !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
       !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
       !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
       !
       !  Storage of this value over the value of a B-spline
       !  of order M there is safe since the remaining B-spline derivatives
       !  of the same order do not use this value due to the fact
       !  that  A(J,I) = 0  for J < I.
       !
       do i = 1, k
          
          jlow = max ( i, m )
          
          dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )
          
       end do
    
    end do
    return
  end subroutine bsplvd

function bvalue( deboor, t, bcoef, n, k, x, jderiv )
    
implicit none

type(deboor_type)        :: deboor
sll_real64               :: bvalue

sll_real64, dimension(:) :: t     !(n+k)
sll_real64, dimension(:) :: bcoef !(n)
sll_int32                :: n
sll_int32                :: k
sll_real64               :: x
sll_int32                :: jderiv

sll_real64, dimension(k) :: aj
sll_real64, dimension(k) :: dl
sll_real64, dimension(k) :: dr

sll_int32 :: i
sll_int32 :: ilo
sll_int32 :: j
sll_int32 :: jc
sll_int32 :: jcmax
sll_int32 :: jcmin
sll_int32 :: jj
sll_int32 :: mflag
sll_int32 :: ierr

bvalue = 0.0_8

aj(:)=0.0_8
dl(:)=0.0_8
dr(:)=0.0_8

if ( k <= jderiv ) return

!  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
!
!  If no such I can be found, X lies outside the support of the
!  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
!  right continuous.

call interv ( deboor, t, n+k, x, i, mflag )

if ( mflag /= 0 ) return
!
!  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
!
if ( k <= 1 ) then
   bvalue = bcoef(i)
   return
end if
!
!  Store the K B-spline coefficients relevant for the knot interval
!  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
!  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
!  from input to zero.
!
!  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
!
jcmin = 1

if ( k <= i ) then
   
   do j = 1, k-1
      dl(j) = x - t(i+1-j)
   end do
   
else
   
   jcmin = 1 - ( i - k )
   do j = 1, i
      dl(j) = x - t(i+1-j)
   end do
   do j = i, k-1
      aj(k-j) = 0.0_8
      dl(j) = dl(i)
   end do
   
end if

jcmax = k

if ( n < i ) then
   
   jcmax = k + n - i
   do j = 1, k + n - i
      dr(j) = t(i+j) - x
   end do
   
   do j = k+n-i, k-1
      aj(j+1) = 0.0_8
      dr(j) = dr(k+n-i)
   end do
   
else
   
   do j = 1, k-1
      dr(j) = t(i+j) - x
   end do
   
end if

do jc = jcmin, jcmax
   aj(jc) = bcoef(i-k+jc)
end do
!
!  Difference the coefficients JDERIV times.
!
do j = 1, jderiv
   
   ilo = k - j
   do jj = 1, k - j
      aj(jj) = ((aj(jj+1)-aj(jj))/(dl(ilo)+dr(jj)))*(k-j)
      ilo = ilo - 1
   end do
   
end do
!
!  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
!  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
!
do j = jderiv+1, k-1
   ilo = k-j
   do jj = 1, k-j
      aj(jj) = (aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
      ilo = ilo - 1
   end do
end do

bvalue = aj(1)

return

end function bvalue
  




   subroutine spli1d_der(&
       ai_nx,&
       ai_nx_der,&
       ai_kx,&
       apr_taux,&
       apr_g,&
       apr_taux_der,&
       apr_g_der,&
       apr_Bcoef,&
       apr_tx)
    implicit none
    ! INPUT
    sll_int32  :: ai_nx, ai_kx,ai_nx_der
    sll_real64, dimension ( ai_nx) :: apr_taux
    sll_real64, dimension ( ai_nx) :: apr_g
    sll_int32, dimension ( ai_nx_der) :: apr_taux_der
    sll_real64, dimension ( ai_nx_der) :: apr_g_der
    ! OUTPUT
    sll_real64, dimension ( ai_nx + ai_nx_der) :: apr_Bcoef
    sll_real64, dimension ( ai_nx + ai_nx_der + ai_kx) :: apr_tx
    ! LOCAL VARIABLES		
    sll_real64, dimension ( (ai_nx+ai_nx_der) *( 2*ai_kx-1) ) :: lpr_work
    !sll_real64, dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
    sll_int32 :: iflag
    
    apr_tx = 0.0_f64
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx+ ai_nx_der + 1: ai_nx + ai_nx_der + ai_kx ) = apr_taux ( ai_nx )
  
    
    if (ai_nx + ai_nx_der + ai_kx == ai_nx + 2*(ai_kx-1)) then
       apr_tx (ai_kx+1: ai_nx+ ai_nx_der) = apr_taux(2:ai_nx-1)
       
    else
       print*, 'problem with construction of knots' 
    end if


    call splint_der ( &
         apr_taux,&
         apr_g,&
         apr_taux_der,&
         apr_g_der,&
         apr_tx,&
         ai_nx,ai_nx_der,ai_kx,&
         lpr_work, apr_Bcoef, iflag )
    
    
  end subroutine spli1d_der
  
  


  subroutine splint ( deboor, tau, gtau, t, n, k, q, bcoef, iflag )
    
    implicit none
    
    sll_int32 ::  n
    
    type(deboor_type) :: deboor
    sll_real64,dimension(:):: bcoef ! (n)
    sll_real64,dimension(:)::  gtau ! (n)
    sll_int32 :: i
    sll_int32 :: iflag
    sll_int32 :: ilp1mx
    sll_int32 :: j
    sll_int32 :: jj
    sll_int32 :: k
    sll_int32 :: kpkm2
    sll_int32 :: left
    sll_real64, dimension((2*k-1)*n) :: q!((2*k-1)*n)
    sll_real64,dimension(n+k) ::  t!(n+k)
    sll_real64,dimension(n) ::  tau!!(n)
    sll_real64:: taui
  
    kpkm2 = 2 * ( k - 1 )
    left = k
    q(1:(2*k-1)*n) = 0.0_f64
    !
    !  Loop over I to construct the N interpolation equations.
    !
    do i = 1, n
       
       taui = tau(i)
       ilp1mx = min ( i + k, n + 1 )
       !
       !  Find LEFT in the closed interval (I,I+K-1) such that
       !
       !    T(LEFT) <= TAU(I) < T(LEFT+1)
       !
       !  The matrix is singular if this is not possible.
       !
       left = max ( left, i )
       
       if ( taui < t(left) ) then
          iflag = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINT - Fatal Error!'
          write ( *, '(a)' ) '  The linear system is not invertible!'
          return
       end if
   
       do while ( t(left+1) <= taui )
          
          left = left + 1
          
          if ( left < ilp1mx ) then
             cycle
          end if
          
          left = left - 1
          
          if ( t(left+1) < taui ) then
             iflag = 2
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'SPLINT - Fatal Error!'
             write ( *, '(a)' ) '  The linear system is not invertible!'
             return
          end if
          
          exit
          
       end do
       !
       !  The I-th equation enforces interpolation at TAUI, hence for all J,
       !
       !    A(I,J) = B(J,K,T)(TAUI).
       !
       !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
       !
       !These K numbers are returned, in BCOEF 
       ! (used for temporary storage here),
       !  by the following.
       !
       call bsplvb ( deboor, t, k, 1, taui, left, bcoef )
       !
       !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
       !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
       !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
       !  as a two-dimensional array, with  2*K-1 rows.  See comments in
       !  BANFAC.
       !
       !  In the present program, we treat Q as an equivalent
       !  one-dimensional array, because of fortran restrictions on
       !  dimension statements.
       !
       !  We therefore want  BCOEF(J) to go into the entry of Q with index:
       !
       !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
       !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
       !
       jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
       
       do j = 1, k
          jj = jj + kpkm2
          q(jj) = bcoef(j)
       end do
       
    end do
    !
    !  Obtain factorization of A, stored again in Q.
    !
    call banfac ( q, k+k-1, n, k-1, k-1, iflag )
    
    if ( iflag == 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'SPLINT - Fatal Error!'
       write ( *, '(a)' ) '  The linear system is not invertible!'
       return
    end if
    !
    !  Solve 
    !
    !    A * BCOEF = GTAU
    !
    !  by back substitution.
    !
    bcoef(1:n) = gtau(1:n)
    
    call banslv ( q, k+k-1, n, k-1, k-1, bcoef )
    
    return
  end subroutine splint


 subroutine splint_der( tau,gtau,tau_der,gtau_der,t,n,m,k, q, bcoef_spline, iflag )
    
    implicit none
    
    type(deboor_type) :: deboor
    sll_int32 ::  n
    sll_int32 ::  m
    sll_real64,dimension(n+m):: bcoef ! (n)
    sll_real64,dimension(n+m):: bcoef_spline 
    sll_real64,dimension(n)::  gtau ! (n)
    sll_real64,dimension(m)::  gtau_der ! (n)
    sll_int32 :: i
    sll_int32 :: iflag,mflag
    sll_int32 :: j,l
    sll_int32 :: jj
    sll_int32 :: k
    sll_int32 :: kpkm2
    sll_int32 :: left
    sll_real64, dimension((2*k-1)*(n+m)) :: q!((2*k-1)*n)
    sll_real64,dimension(n+k+m) ::  t!(n+k)
    sll_real64,dimension(n) ::  tau!!(n)
    sll_int32,dimension(m) ::  tau_der!!(n)
    sll_real64:: taui,taui_der
    sll_real64, dimension(k,k):: a
    sll_real64,dimension(k,2) :: bcoef_der
  
    kpkm2 = 2 * ( k - 1 )
    left = k
    q(1:(2*k-1)*(n+m)) = 0.0_f64
    a(1:k,1:k) = 0.0_f64
    bcoef_der(1:k,1:2) = 0.0_f64

    ! we must suppose that m is <= than n 
    if (m > n) then
       print*, 'problem m must be < = at n'
       print*, 'value m =', m, 'value n =', n
       stop
    end if
    l = 1 ! index for the derivative
    !
    !  Loop over I to construct the N interpolation equations.
    !
    do i = 1, n-1
       
       taui = tau(i)
       
       !
       !  Find LEFT in the closed interval (I,I+K-1) such that
       !
       !    T(LEFT) <= TAU(I) < T(LEFT+1)
       !
       !  The matrix is singular if this is not possible.
       !  With help of the Schoenberg-Whitney theorem 
       !  we can prove that if the diagonal of the 
       !  matrix B_j(x_i) is null, we have a non-inversible matrix.  

       call interv( deboor, t, n+m+k, taui, left, mflag )

       !

       !
       !  The I-th equation enforces interpolation at TAUI, hence for all J,
       !
       !    A(I,J) = B(J,K,T)(TAUI).
       !
       !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
       !
       !These K numbers are returned, in BCOEF 
       ! (used for temporary storage here),
       !  by the following.
       !
       
       call bsplvb ( deboor, t, k, 1, taui, left, bcoef )
       
          
       !
       !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
       !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
       !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
       !  as a two-dimensional array, with  2*K-1 rows.  See comments in
       !  BANFAC.
       !
       !  In the present program, we treat Q as an equivalent
       !  one-dimensional array, because of fortran restrictions on
       !  dimension statements.
       !
       !  We therefore want  BCOEF(J) to go into the entry of Q with index:
       !
       !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
       !    =  begin_ligne +  (begin_col -1) * number_coef_different_0
       !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
       !
       jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
       
       do j = 1, k
          jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
          q(jj) = bcoef(j)
       end do

       bcoef_spline(i+ l-1) = gtau(i)
       if ( tau_der(l) == i ) then   
          taui_der = taui
          
          call bsplvd( deboor, t, k, taui_der, left, a, bcoef_der, 2)

          l = l + 1
          jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
       
          do j = 1, k
             jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
             q(jj) = bcoef_der(j,2)
          end do
       bcoef_spline(i+ l-1) = gtau_der(l-1)
       end if
       

    end do


    taui = tau(n)
    call interv( deboor, t, n+m+k, taui, left, mflag )
    if ( tau_der(l)== n ) then   
          taui_der = taui
          
          call bsplvd( deboor, t, k, taui_der, left, a, bcoef_der, 2)

          
          jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
       
          do j = 1, k
             jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
             q(jj) = bcoef_der(j,2)
          end do
          bcoef_spline(n+ l-1) = gtau_der(l)
          l = l + 1
          
       end if
       

     call bsplvb ( deboor, t, k, 1, taui, left, bcoef )
     jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
       
     do j = 1, k
        jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
        q(jj) = bcoef(j)
     end do
     bcoef_spline(n+l-1) = gtau(n)

    !
    !  Obtain factorization of A, stored again in Q.
    !
    call banfac ( q, k+k-1, n+m, k-1, k-1, iflag )
    
    if ( iflag == 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'SPLINT - Fatal Error!'
       write ( *, '(a)' ) '  The linear system is not invertible!'
       return
    end if
    !
    !  Solve 
    !
    !    A * BCOEF = GTAU
    !
    !  by back substitution.
    !
    


    call banslv ( q, k+k-1, n+m, k-1, k-1, bcoef_spline )
    
    return
  end subroutine splint_der

  !> @brief initializing the coefficients of splines.
  !> @details  initializing the coefficients of splines
  !>  fot the arbitrary degree splines interpolator 1d
  !> The parameters are
  !> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
  !> @param[in] coeffs the 1d arrays corresponding of the splines coefficients
  !> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d

  subroutine set_coefficients_ad1d( &
   interpolator, &
   coeffs)

   class(sll_arbitrary_degree_spline_interpolator_1d), intent(inout)  :: interpolator
   sll_real64, dimension(:), intent(in), optional :: coeffs
   sll_int32 :: sp_deg
   sll_int32 :: num_cells
   sll_int32 :: tmp
   sll_int32 :: i
   sll_real64 :: eta_min, eta_max
   sll_real64 :: delta
   sll_int32  ::  nb_spline_eta
   sll_real64 :: eta


   sp_deg    = interpolator%spline_degree
   num_cells = interpolator%num_pts - 1
   eta_min   = interpolator%eta_min
   eta_max   = interpolator%eta_max
   delta     = (eta_max - eta_min)/num_cells

   tmp = (sp_deg + 1)/2

   ! The interpretation and further filling of the spline coefficients array
   ! depends on the boundary conditions.
   select case (interpolator%bc_selector)
   case(0) ! periodic

      interpolator%size_coeffs =  num_cells + sp_deg
      interpolator%size_t = 2*sp_deg + num_cells +1
      if ( size(coeffs) .ne.  num_cells + 1 ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',num_cells + 1
         print*, 'and not =', size(coeffs)
         stop
      endif
      ! allocation and definition of knots
      do i = -sp_deg, num_cells + sp_deg
         interpolator%t(i+sp_deg+1)=eta_min+i*delta
      end do
      do i = 1,num_cells
         interpolator%coeff_splines(i) = coeffs( i )
      end do
      do i = 1, sp_deg
         interpolator%coeff_splines(num_cells + i ) = coeffs(i )
      end do
      do i= 1,sp_deg
         interpolator%coeff_splines(num_cells +  i ) = &
              interpolator%coeff_splines(sp_deg-(i-1))
      end do

   case (9) ! 2. dirichlet-left, dirichlet-right
      interpolator%size_coeffs=  num_cells + sp_deg
      interpolator%size_t = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 2
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      ! allocation and definition of knots

      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo


      do i = 1,nb_spline_eta

         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do

      interpolator%coeff_splines(1)               = interpolator%value_left
      interpolator%coeff_splines(nb_spline_eta+2) = interpolator%value_right

   case(10) ! Neumann - Dirichlet

      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 1

      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif

      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do
      interpolator%coeff_splines(nb_spline_eta+2) =  interpolator%value_right
   case(12) ! Hermitte- Dirichlet
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 1
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do
      interpolator%coeff_splines(nb_spline_eta+2) = interpolator%value_right
   case(17) ! Dirichlet-Neumann

      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 1
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif

      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do
      interpolator%coeff_splines(1) = interpolator%value_left
   case(18) ! Neumann - Neumann
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i ) =  coeffs(i)
      end do
   case(20) ! Hermite - Neumann
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i) =  coeffs(i)
      end do
   case(33) ! Dirichlet - Hermite
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg - 1
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i + 1 ) =  coeffs(i)
      end do
      interpolator%coeff_splines(1) = interpolator%value_left
   case(34) ! Neumann- Hermite
      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i) =  coeffs(i)
      end do
   case(36)! Hermite - Hermite

      interpolator%size_coeffs= num_cells + sp_deg
      interpolator%size_t     = 2*sp_deg + num_cells + 1
      nb_spline_eta = num_cells + sp_deg
      if ( size(coeffs) .ne.  nb_spline_eta ) then
         print*, 'problem in set_coeff_1d_arb_deg_spline '
         print*, 'size coeffs must be equal to ',nb_spline_eta
         print*, 'and not =', size(coeffs)
         stop
      endif
      do i = 1, sp_deg + 1
         interpolator%t(i) = eta_min
      enddo
      eta = eta_min
      do i = sp_deg + 2, num_cells + 1 + sp_deg
         eta = eta + delta
         interpolator%t(i) = eta
      enddo
      do i = num_cells + sp_deg + 2, num_cells + 1 + 2*sp_deg
         interpolator%t(i) = eta
      enddo

      do i = 1,nb_spline_eta
         interpolator%coeff_splines(i ) =  coeffs(i)
      end do
   case default
      print *, 'arbitrary_degree_spline_1d() error: set_spline_coefficients ',&
           'not recognized.'
      stop
   end select
 end subroutine set_coefficients_ad1d

  
end module sll_module_arbitrary_degree_spline_interpolator_1d
