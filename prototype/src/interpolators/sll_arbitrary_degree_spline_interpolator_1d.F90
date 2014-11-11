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

#define sll_interpolator class(sll_arbitrary_degree_spline_interpolator_1d)

!> @ingroup interpolators
!> @brief
!> Class interpolator and methods for arbitrary degree spline 1D interpolator
!> @details
!> Arbitrary degree splines are implemented in this class that derived from
!> sll_interpolator_1d_base.
module sll_module_arbitrary_degree_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

use sll_module_interpolators_1d_base
use sll_module_deboor_splines_1d

implicit none
private

!> Class for arbitrary degree spline 1d interpolator
type, public, extends(sll_interpolator_1d_base) :: &
                 sll_arbitrary_degree_spline_interpolator_1d
  private
  !> PLEASE ADD DOCUMENTATION
  sll_real64, public, dimension(:), pointer :: bcoef
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: num_pts
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: eta_min
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: eta_max
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: bc_left
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: bc_right
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: spline_degree
  !> PLEASE ADD DOCUMENTATION
  sll_real64, dimension(:), pointer :: t
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: size_t
  !> PLEASE ADD DOCUMENTATION
  sll_int64  :: bc_selector 
  !> PLEASE ADD DOCUMENTATION
  sll_int32  :: size_coeffs
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: slope_left
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: slope_right
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: value_left
  !> PLEASE ADD DOCUMENTATION
  sll_real64 :: value_right
  !> PLEASE ADD DOCUMENTATION
  logical    :: compute_slope_left = .TRUE.
  !> PLEASE ADD DOCUMENTATION
  logical    :: compute_slope_right= .TRUE.
  !> PLEASE ADD DOCUMENTATION
  logical    :: compute_value_left = .TRUE.
  !> PLEASE ADD DOCUMENTATION
  logical    :: compute_value_right= .TRUE.

contains

  !> PLEASE ADD DOCUMENTATION
  procedure, pass(interpolator) :: initialize => initialize_ad1d_interpolator
  !> PLEASE ADD DOCUMENTATION
  procedure, pass :: set_coefficients => set_coefficients_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: compute_interpolants => compute_interpolants_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_value => interpolate_value_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array_values => interpolate_values_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_pointer_values => interpolate_pointer_values_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_derivative_eta1 => interpolate_derivative_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array_derivatives => interpolate_derivatives_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_pointer_derivatives =>interpolate_pointer_derivatives_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass:: interpolate_array => interpolate_array_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass:: interpolate_array_disp => interpolate_1d_array_disp_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass:: get_coefficients => get_coefficients_ad1d
  !> PLEASE ADD DOCUMENTATION
  procedure, pass:: reconstruct_array

end type sll_arbitrary_degree_spline_interpolator_1d

!> Deallocate
interface sll_delete
  module procedure delete_arbitrary_degree_1d_interpolator
end interface sll_delete

public sll_delete 
public new_arbitrary_degree_1d_interpolator
public set_values_at_boundary1d
public initialize_ad1d_interpolator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief delete interpolator arbitrary degree splines.
!> @details
!> The parameters are
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine delete_arbitrary_degree_1d_interpolator( interpolator )
  sll_interpolator, intent(inout) :: interpolator
  sll_int32 :: ierr
  DEALLOCATE(interpolator%t)
  DEALLOCATE(interpolator%bcoef)
end subroutine delete_arbitrary_degree_1d_interpolator

!> @brief Initialization of a pointer interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d such as a pointer.
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @return the type interpolator arbitrary degree splines 1d
function new_arbitrary_degree_1d_interpolator(     &
       num_pts,                                    &
       eta_min,                                    &
       eta_max,                                    &
       bc_left,                                    &
       bc_right,                                   &
       spline_degree) result(interpolator)

  sll_interpolator, pointer :: interpolator

  sll_int32,  intent(in) :: num_pts
  sll_real64, intent(in) :: eta_min
  sll_real64, intent(in) :: eta_max
  sll_int32,  intent(in) :: bc_left
  sll_int32,  intent(in) :: bc_right
  sll_int32,  intent(in) :: spline_degree
  sll_int32              :: ierr

  SLL_ALLOCATE(interpolator,ierr)

  call initialize_ad1d_interpolator(interpolator, &
                                    num_pts,      &
                                    eta_min,      &
                                    eta_max,      &
                                    bc_left,      &
                                    bc_right,     &
                                    spline_degree)

end function new_arbitrary_degree_1d_interpolator

!> @brief Initialization of interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d.
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

  sll_interpolator, intent(inout) :: interpolator

  sll_int32,  intent(in) :: num_pts
  sll_real64, intent(in) :: eta_min
  sll_real64, intent(in) :: eta_max
  sll_int32,  intent(in) :: bc_left
  sll_int32,  intent(in) :: bc_right
  sll_int32,  intent(in) :: spline_degree
  sll_int32 :: ierr
  sll_int64 :: bc_selector

  if(bc_left  == SLL_PERIODIC .or. bc_right == SLL_PERIODIC) then
    SLL_ASSERT(bc_right == SLL_PERIODIC .and. bc_left  == SLL_PERIODIC)
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

  SLL_CLEAR_ALLOCATE(interpolator%bcoef(1:num_pts*num_pts),ierr)
  SLL_CLEAR_ALLOCATE(interpolator%t(1:num_pts*num_pts),ierr)
  interpolator%value_left  = 0.0_f64
  interpolator%value_right = 0.0_f64

end subroutine initialize_ad1d_interpolator

!> @brief initializing the coefficients of splines.
!> @details  initializing the coefficients of splines
!> for the arbitrary degree splines interpolator 1d.
!> @param interpolator the type sll_arbitrary_degree_spline_interpolator_1d
!> @param[in] coeffs the 1d arrays corresponding of the splines coefficients
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_coefficients_ad1d(interpolator, coeffs)

  sll_interpolator, intent(inout)  :: interpolator
  sll_real64, intent(in), optional :: coeffs(:)
  sll_int32                        :: sp_deg
  sll_int32                        :: num_cells
  sll_int32                        :: i
  sll_real64                       :: eta_min
  sll_real64                       :: eta_max
  sll_real64                       :: delta
  sll_int32                        :: nb_spline_eta
  sll_real64                       :: eta

  sp_deg    = interpolator%spline_degree
  num_cells = interpolator%num_pts - 1
  eta_min   = interpolator%eta_min
  eta_max   = interpolator%eta_max
  delta     = (eta_max - eta_min)/num_cells

  interpolator%size_coeffs = num_cells + sp_deg
  interpolator%size_t      = 2*sp_deg + num_cells +1

  ! The interpretation and further filling of the spline coefficients array
  ! depends on the boundary conditions.
  select case (interpolator%bc_selector)
  case(0) ! periodic

    do i = -sp_deg, num_cells + sp_deg
      interpolator%t(i+sp_deg+1) = eta_min + i*delta
    end do

    if (present(coeffs)) then

      SLL_ASSERT ( size(coeffs) ==  num_cells + sp_deg )
      interpolator%bcoef(1:num_cells+sp_deg) = coeffs
      do i = 1, sp_deg
        interpolator%bcoef(num_cells+i) = coeffs(i)
      end do

    else

      do i= 1,sp_deg
        interpolator%bcoef(num_cells+i ) = &
                interpolator%bcoef(sp_deg-(i-1))
      end do

    end if

  case (9) ! 2. dirichlet-left, dirichlet-right

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

    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg - 2
      SLL_ASSERT(size(coeffs) ==  nb_spline_eta )
      do i = 1,nb_spline_eta
        interpolator%bcoef(i + 1 ) =  coeffs(i)
      end do
    endif

    interpolator%bcoef(1)               = interpolator%value_left
    interpolator%bcoef(nb_spline_eta+2) = interpolator%value_right

  case(10) ! Neumann - Dirichlet

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
 
    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg - 1
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i + 1 ) =  coeffs(i)
      end do
    end if

    interpolator%bcoef(nb_spline_eta+2) =  interpolator%value_right
 
  case(12) ! Hermitte- Dirichlet
 
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
 
    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg - 1
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i + 1 ) =  coeffs(i)
      end do
    end if
    interpolator%bcoef(nb_spline_eta+2) = interpolator%value_right
 
  case(17) ! Dirichlet-Neumann

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

    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg - 1
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i + 1 ) =  coeffs(i)
      end do
    end if
    interpolator%bcoef(1) = interpolator%value_left

  case(18) ! Neumann - Neumann

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

    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i ) =  coeffs(i)
      end do
    end if

  case(20) ! Hermite - Neumann

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

    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i) =  coeffs(i)
      end do
    end if

  case(33) ! Dirichlet - Hermite

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

    if (present(coeffs)) then
      nb_spline_eta = num_cells + sp_deg - 1
      SLL_ASSERT(size(coeffs) == nb_spline_eta)
      do i = 1,nb_spline_eta
        interpolator%bcoef(i + 1 ) =  coeffs(i)
      end do
    end if
    interpolator%bcoef(1) = interpolator%value_left

   case(34) ! Neumann- Hermite

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

     if (present(coeffs)) then
       nb_spline_eta = num_cells + sp_deg
       SLL_ASSERT(size(coeffs) == nb_spline_eta)
       do i = 1,nb_spline_eta
         interpolator%bcoef(i) =  coeffs(i)
       end do
     end if

   case(36)! Hermite - Hermite

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

     if (present(coeffs)) then
       nb_spline_eta = num_cells + sp_deg
       SLL_ASSERT(size(coeffs) == nb_spline_eta)
       do i = 1,nb_spline_eta
         interpolator%bcoef(i ) =  coeffs(i)
       end do
     end if

   case default

     SLL_ERROR('set_spline_coefficients, BC not recognized.')
    
  end select

end subroutine set_coefficients_ad1d

!> Initialization of the boundary for interpolator arbitrary degree splines 1d.
!> The parameters are
!> @param[in] value_left  contains the value in the left
!> @param[in] value_right contains the value in the right
!> @param[in] slope_left  contains the value in the left for derivative
!> @param[in] slope_right contains the value in the right for derivative
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_values_at_boundary1d( interpolator, &
                                     value_left,   &
                                     value_right,  &
                                     slope_left,   &
                                     slope_right)

  sll_interpolator, intent(inout)  :: interpolator
  sll_real64, intent(in), optional :: value_left
  sll_real64, intent(in), optional :: value_right
  sll_real64, intent(in), optional :: slope_left
  sll_real64, intent(in), optional :: slope_right

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

  sll_interpolator, intent(inout)                :: interpolator
  sll_real64, dimension(:), intent(in)           :: data_array
  sll_real64, dimension(:), intent(in), optional :: eta_coords
  sll_int32,  intent(in),   optional             :: size_eta_coords

  sll_real64, dimension(:), pointer              :: point_locate_eta

  sll_int32  :: sz
  sll_real64 :: delta_eta
  sll_real64 :: period
  sll_int32  :: order
  sll_int32  :: ierr
  sll_int32  :: i

  ! ----------------------------!
  ! It is only for cubic spline !
  ! ----------------------------!
  sll_int32, parameter            :: sz_deriv = 2
  sll_real64, dimension(sz_deriv) :: data_array_derivative
  sll_int32,  dimension(sz_deriv) :: point_locate_eta_derivative

  if( present(eta_coords) ) then

    SLL_ASSERT(present(size_eta_coords))

    sz = size(data_array)
    SLL_ALLOCATE(point_locate_eta(sz),ierr)
    point_locate_eta = eta_coords
    SLL_ASSERT(sz .le. interpolator%num_pts*interpolator%num_pts)

  else 

    sz = interpolator%num_pts
    SLL_ASSERT(size(data_array) .le. sz)

    delta_eta = (interpolator%eta_max - interpolator%eta_min)&
          /(interpolator%num_pts -1)
    SLL_ALLOCATE(point_locate_eta(sz),ierr)

    do i = 1,sz
      point_locate_eta(i) = interpolator%eta_min + delta_eta*(i-1)
    end do

  end if

  if (interpolator%compute_slope_left) then
    interpolator%slope_left = forward_fd_5pt( data_array,point_locate_eta)
  end if
  if (interpolator%compute_slope_right) then
    interpolator%slope_right = backward_fd_5pt( data_array,point_locate_eta,sz)
  end if
  if (interpolator%compute_value_left) then
    interpolator%value_left = data_array(1)
  end if
  if (interpolator%compute_value_right) then
    interpolator%value_right = data_array(sz)
  end if

  order  = interpolator%spline_degree + 1
  period = interpolator%eta_max - interpolator%eta_min


  select case (interpolator%bc_selector)

  case (0) ! periodic

    interpolator%size_coeffs = sz 
    interpolator%size_t = order + sz 
    call spli1d_per( period,                           &
                     sz,                               &
                     order,                            &
                     point_locate_eta,                 &
                     data_array,                       &
                     interpolator%bcoef(1:sz), &
                     interpolator%t(1:order + sz ))

  case (9) ! 2. dirichlet-left, dirichlet-right

    interpolator%size_coeffs = sz
    interpolator%size_t      = order + sz

    call spli1d_dir( sz,                               &
                     order,                            &
                     point_locate_eta,                 &
                     data_array,                       &
                     interpolator%bcoef(1:sz), &
                     interpolator%t(1:sz+order) )

    interpolator%bcoef(1)  = interpolator%value_left
    interpolator%bcoef(sz) = interpolator%value_right

  case(10) ! Neumann - Dirichlet

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = 0.0_f64
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = interpolator%slope_right

    call spli1d_der(sz,                             &
         sz_deriv,                                  &
         order,                                     &
         point_locate_eta,                          &
         data_array,                                &
         point_locate_eta_derivative,               &
         data_array_derivative,                     &
         interpolator%bcoef(1:sz+sz_deriv), &
         interpolator%t(1:sz+order+sz_deriv))

    ! test dirichlet non homogene
    interpolator%bcoef(sz+sz_deriv) = interpolator%value_right

 case(17) ! Dirichlet - Neumann

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = interpolator%slope_left
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = 0.0_f64

    call spli1d_der(sz,sz_deriv,order,              &
         point_locate_eta,                          &
         data_array,                                &
         point_locate_eta_derivative,               &
         data_array_derivative,                     &
         interpolator%bcoef(1:sz+sz_deriv), &
         interpolator%t(1:sz+order+sz_deriv))

    ! test dirichlet non homogene
    interpolator%bcoef(1) = interpolator%value_left

 case(12) ! Hermite - Dirichlet

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = interpolator%slope_left
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = interpolator%slope_right

    call spli1d_der(sz,sz_deriv,order,              &
         point_locate_eta,                          &
         data_array,                                &
         point_locate_eta_derivative,               &
         data_array_derivative,                     &
         interpolator%bcoef(1:sz+sz_deriv), &
         interpolator%t(1:sz+order+sz_deriv))

    ! test dirichlet non homogene
    interpolator%bcoef(sz+sz_deriv) = interpolator%value_right

 case(18) ! Neumann - Neumann

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = 0.0_f64
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = 0.0_f64
    call spli1d_der(sz,sz_deriv,order,              &
         point_locate_eta,                          &
         data_array,                                &
         point_locate_eta_derivative,               &
         data_array_derivative,                     &
         interpolator%bcoef(1:sz+sz_deriv), & 
         interpolator%t(1:sz+order+sz_deriv))

 case(20) ! Hermite - Neumann

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = interpolator%slope_left
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = 0.0_f64

    call spli1d_der(sz,sz_deriv,order,              &
         point_locate_eta,                          &
         data_array,                                &
         point_locate_eta_derivative,               &
         data_array_derivative,                     &
         interpolator%bcoef(1:sz+sz_deriv), &
         interpolator%t(1:sz+order+sz_deriv))

 case(33) ! Dirichlet - Hermite

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = interpolator%slope_left
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = interpolator%slope_right

    call spli1d_der(sz,sz_deriv,order,                &
         point_locate_eta,                            &
         data_array,                                  &
         point_locate_eta_derivative,                 &
         data_array_derivative,                       &
         interpolator%bcoef(1:sz+sz_deriv),   &
         interpolator%t(1:sz+order+sz_deriv))

    interpolator%bcoef(1) = interpolator%value_left

 case(34) ! Neumann - Hermite

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = 0.0_f64
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = interpolator%slope_right

    call spli1d_der(sz,sz_deriv,order,               &
         point_locate_eta,                           &
         data_array,                                 &
         point_locate_eta_derivative,                &
         data_array_derivative,                      &
         interpolator%bcoef(1:sz+sz_deriv),  &
         interpolator%t(1:sz+order+sz_deriv))

 case(36) ! Hermite - Hermite

    interpolator%size_coeffs = sz + sz_deriv
    interpolator%size_t = order + sz + sz_deriv

    point_locate_eta_derivative(1) = 1
    data_array_derivative(1)       = interpolator%slope_left
    point_locate_eta_derivative(2) = sz
    data_array_derivative(2)       = interpolator%slope_right

    call spli1d_der(sz,sz_deriv,order,               &
         point_locate_eta,                           &
         data_array,                                 &
         point_locate_eta_derivative,                &
         data_array_derivative,                      &
         interpolator%bcoef(1:sz+sz_deriv),  &
         interpolator%t(1:sz+order+sz_deriv))

  end select

  DEALLOCATE(point_locate_eta)

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

  sll_interpolator, intent(in) :: interpolator
  sll_real64, intent(in)       :: eta1
  sll_real64                   :: val
  sll_real64                   :: res

  res = eta1

  if (interpolator%bc_selector == 0) then ! periodic

    if( res < interpolator%eta_min ) then
      res = res+interpolator%eta_max-interpolator%eta_min
    else if( res >  interpolator%eta_max ) then
      res = res+interpolator%eta_min-interpolator%eta_max
    end if

  else

      SLL_ASSERT( res >= interpolator%eta_min )
      SLL_ASSERT( res <= interpolator%eta_max )

  end if

  val = bvalue( interpolator%t(1:interpolator%size_t),                  &
                interpolator%bcoef(1:interpolator%size_coeffs), &
                interpolator%size_coeffs,                               &
                interpolator%spline_degree+1,                           &
                res,                                                    &
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

function interpolate_derivative_ad1d( interpolator, eta1 ) result(val)

  sll_interpolator, intent(in) :: interpolator
  sll_real64, intent(in)       :: eta1
  sll_real64                   :: val
  sll_int32                    :: size_coeffs
  sll_real64                   :: res

  SLL_ASSERT( eta1 .ge. interpolator%eta_min )
  SLL_ASSERT( eta1 .le. interpolator%eta_max )

  size_coeffs = interpolator%size_coeffs

  res = eta1

  if(interpolator%bc_selector == 0) then

    if( res < interpolator%eta_min ) then
      res = res+interpolator%eta_max-interpolator%eta_min
    else if( res >  interpolator%eta_max ) then
      res = res+interpolator%eta_min-interpolator%eta_max
    end if

  else

    SLL_ASSERT( res >= interpolator%eta_min )
    SLL_ASSERT( res <= interpolator%eta_max )

  end if

  val = bvalue( interpolator%t(1:interpolator%size_t),     &
                interpolator%bcoef(1:size_coeffs), &
                size_coeffs,                               &
                interpolator%spline_degree+1,              &
                res,                                       &
                1 )
    
end function interpolate_derivative_ad1d

function interpolate_array_ad1d( this,        &
                                 num_points,  &
                                 data,        &
                                 coordinates) result(data_out)

  sll_interpolator,         intent(in) :: this
  sll_int32,                intent(in) :: num_points
  sll_real64, dimension(:), intent(in) :: coordinates
  sll_real64, dimension(:), intent(in) :: data
  sll_real64, dimension(num_points)    :: data_out

  SLL_WARNING('interpolate_array_ad1d: not implemented')
  data_out = -1000000._f64*data*coordinates*this%spline_degree

end function interpolate_array_ad1d

function interpolate_1d_array_disp_ad1d( this,        &
                                         num_points,  &
                                         data,        &
                                         alpha) result(res)

  sll_interpolator,         intent(in) :: this
  sll_int32,                intent(in) :: num_points
  sll_real64, dimension(:), intent(in) :: data
  sll_real64,               intent(in) :: alpha
  sll_real64, dimension(num_points)    :: res

  SLL_WARNING('interpolate_1d_array_disp_ad1d: not implemented')
  res = -1000000._f64*alpha*data*this%spline_degree

end function interpolate_1d_array_disp_ad1d

function get_coefficients_ad1d(interpolator)

   sll_interpolator, intent(in)      :: interpolator
   sll_real64, dimension(:), pointer :: get_coefficients_ad1d

   get_coefficients_ad1d = interpolator%bcoef

end function get_coefficients_ad1d

subroutine interpolate_values_ad1d( interpolator,        &
                                    num_pts,             &
                                    vals_to_interpolate, &
                                    output_array )

  sll_interpolator,         intent(in)  :: interpolator
  sll_int32,                intent(in)  :: num_pts
  sll_real64, dimension(:), intent(in)  :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array
  sll_int32 :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))

  do idx=1,num_pts
    output_array(idx)=interpolate_value_ad1d(            &
                                interpolator,            &
                                vals_to_interpolate(idx))
  enddo

end subroutine interpolate_values_ad1d

subroutine interpolate_pointer_values_ad1d( interpolator,        &
                                            num_pts,             &
                                            vals_to_interpolate, &
                                            output )

  sll_interpolator,      intent(in) :: interpolator
  sll_int32,             intent(in) :: num_pts
  sll_real64, dimension(:), pointer :: vals_to_interpolate
  sll_real64, dimension(:), pointer :: output
  sll_int32 :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))
  do idx=1,num_pts
    output(idx)=interpolate_value_ad1d( interpolator, vals_to_interpolate(idx))
  enddo

end subroutine interpolate_pointer_values_ad1d

subroutine interpolate_derivatives_ad1d( interpolator,        &
                                         num_pts,             &
                                         vals_to_interpolate, &
                                         output_array )

  sll_interpolator,         intent(in)  :: interpolator
  sll_int32,                intent(in)  :: num_pts
  sll_real64, dimension(:), intent(in)  :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array
  sll_int32                             :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))
  do idx=1,num_pts
    output_array(idx)=interpolate_derivative_ad1d( &
                        interpolator, vals_to_interpolate(idx))
  enddo

end subroutine interpolate_derivatives_ad1d

subroutine interpolate_pointer_derivatives_ad1d( interpolator,        &
                                                 num_pts,             &
                                                 vals_to_interpolate, &
                                                 output )

  sll_interpolator,  intent(in)     :: interpolator
  sll_int32,         intent(in)     :: num_pts
  sll_real64, dimension(:), pointer :: vals_to_interpolate
  sll_real64, dimension(:), pointer :: output
  sll_int32                         :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))
  do idx=1,num_pts
    output(idx)=interpolate_derivative_ad1d( interpolator, &
                                             vals_to_interpolate(idx))
  enddo

end subroutine interpolate_pointer_derivatives_ad1d

function reconstruct_array(this, num_points, data) result(res)

  sll_interpolator,  intent(in)        :: this
  sll_int32, intent(in)                :: num_points
  sll_real64, dimension(:), intent(in) :: data   
  sll_real64, dimension(num_points)    :: res

  res(:) = -1000000.0_f64*data*this%spline_degree
  SLL_WARNING('reconstruct_array 1d not implemented yet')

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

end module sll_module_arbitrary_degree_spline_interpolator_1d
