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

module sll_m_particle_lbf_utilities

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d

  implicit none

  public :: &
    sll_t_int_list_element, &
    sll_t_int_list_element_ptr, &
    sll_f_add_element_in_int_list, &
    sll_s_apply_periodic_bc_x, &
    sll_s_apply_periodic_bc_y, &
    sll_f_eval_hat_function, &
    sll_s_convert_4d_index_to_1d, &
    sll_s_convert_1d_index_to_4d, &
    sll_s_update_closest_particle_arrays, &
    sll_s_get_1d_cell_containing_point


  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! linked lists of integers
  type sll_t_int_list_element
    sll_int32 :: value
    type(sll_t_int_list_element), pointer :: next
  end type sll_t_int_list_element

  type sll_t_int_list_element_ptr ! we need this type for arrays of lists
    type(sll_t_int_list_element), pointer :: pointed_element
  end type sll_t_int_list_element_ptr

contains

  function sll_f_add_element_in_int_list(head, new_element) result( new_list )
    type( sll_t_int_list_element ), pointer :: head, new_element
    type( sll_t_int_list_element ), pointer :: new_list

    new_element%next => head
    new_list => new_element
  end function

  ! name was sll_f_marker_index_from_initial_position_on_cartesian_grid
  subroutine sll_s_convert_4d_index_to_1d(  &
    k, &
    j_x, j_y, j_vx, j_vy,  &
    n_parts_x, n_parts_y, n_parts_vx, n_parts_vy   &
  )
    sll_int32, intent(out):: k
    sll_int32, intent(in) :: j_x
    sll_int32, intent(in) :: j_y
    sll_int32, intent(in) :: j_vx
    sll_int32, intent(in) :: j_vy
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy

    if(         j_x <= 0  .or. j_x > n_parts_x      &
        .or.    j_y <= 0  .or. j_y > n_parts_y      &
        .or.    j_vx <= 0 .or. j_vx > n_parts_vx    &
        .or.    j_vy <= 0 .or. j_vy > n_parts_vy  )then
        k = 0
    else
        k = 1+ (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y
    end if

    SLL_ASSERT( k <= n_parts_x * n_parts_y * n_parts_vx * n_parts_vy )

  end subroutine

  ! name was sll_s_get_initial_position_on_cartesian_grid_from_marker_index
  subroutine sll_s_convert_1d_index_to_4d(  &
    j_x, j_y, j_vx, j_vy,  &
    k, &
    n_parts_x, n_parts_y, n_parts_vx, n_parts_vy  &
  )
    sll_int32, intent(out) :: j_x
    sll_int32, intent(out) :: j_y
    sll_int32, intent(out) :: j_vx
    sll_int32, intent(out) :: j_vy
    sll_int32, intent(in) :: k
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy
    sll_int32              :: k_aux

    k_aux = k-1
    ! here, k_aux = (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y
    j_vy = mod(k_aux, n_parts_vy) + 1

    k_aux = (k_aux - (j_vy-1)) / n_parts_vy
    ! here, k_aux = (j_vx-1) + (j_y-1) * n_parts_vx + (j_x-1) * n_parts_vx * n_parts_y
    j_vx = mod(k_aux, n_parts_vx) + 1

    k_aux = (k_aux - (j_vx-1)) / n_parts_vx
    ! here, k_aux = (j_y-1) + (j_x-1) * n_parts_y
    j_y = mod(k_aux, n_parts_y) + 1

    k_aux = (k_aux - (j_y-1)) / n_parts_y
    ! here, k_aux = (j_x-1)
    j_x = k_aux + 1

    SLL_ASSERT(j_x <= n_parts_x)

  end subroutine

  !> compute the index of the (1d) cell containing a point
  !> in a grid with minimum coordinate eta_min and first cell indexed by 1. Hence cell index j is such that
  !>    eta_min + (j-1)*h <= eta < eta_min + j*h
  !> bounds check should be performed outside this routine
  subroutine sll_s_get_1d_cell_containing_point( j, x, x_min, h )

    sll_int32,  intent(out) :: j
    sll_real64, intent(in)  :: x
    sll_real64, intent(in)  :: x_min
    sll_real64, intent(in)  :: h

    j = int( (x - x_min) / h ) + 1

  end subroutine

  !> sll_f_pic_shape(degree, x, y, vx, vy, inv_hx, inv_hy, inv_hvx, inv_hvy)
  !! computes the value of the 4d B-spline particle shape (no particle transformation here)
  !!
  !! note:
  !!  - inv_hx, inv_hy, inv_hvx, inv_hvy are the inverse
  !!    of the inter-particle distances (hx, hy, hvx, hvy) on the initial (or remapping) grid
  !!
  !!  - x, y, vx, vy are phase-space coordinates relative to the particle center
  !!    and they are not scaled: for instance if degree = 1 (resp. if degree = 3),
  !!    and |x| >= hx (resp. |x| >= 2*hx), then sll_f_pic_shape returns 0
  !!    (and similarly for y, vx, vy)

  !function sll_f_pic_shape( &
  !  degree,             &
  !  x, y, vx, vy,       &
  !  inv_hx, inv_hy, inv_hvx, inv_hvy   &
  !  ) &
  !  result(shape)
  !
  !  sll_int32 :: degree
  !  sll_real64 :: x, y, vx, vy
  !  sll_real64 :: inv_hx, inv_hy, inv_hvx, inv_hvy
  !  sll_real64 :: shape
  !!    sll_int32 :: nc_eta2
  !!    sll_int32 :: ierr
  !
  !  ! uses [[sll_b_spline]]
  !  shape =   inv_hx  * sll_b_spline(degree,inv_hx  * x ) &
  !          * inv_hy  * sll_b_spline(degree,inv_hy  * y ) &
  !          * inv_hvx * sll_b_spline(degree,inv_hvx * vx) &
  !          * inv_hvy * sll_b_spline(degree,inv_hvy * vy)
  !end function sll_f_pic_shape


  !> added by MCP
  !! 4d-reference B-spline shape function (independent of the grid resolution), function support is [-part_cp, part_cp]^4
  !! with part_cp = (degree +1)/2
  !function sll_ref_pic_shape( &
  !  degree,     &
  !  x,y,vx,vy &
  !  ) &
  !  result(ref_shape)
  !
  !  sll_int32 :: degree
  !  sll_real64 :: x, y, vx, vy
  !  sll_real64 :: ref_shape
  !!    sll_int32 :: nc_eta2
  !!    sll_int32 :: ierr
  !
  !  ref_shape = sll_b_spline(degree,x) * sll_b_spline(degree,y) * sll_b_spline(degree,vx) * sll_b_spline(degree, vy)
  !end function sll_ref_pic_shape


  ! added by MCP
  ! <<sll_b_spline>>
  ! univariate centered (and reference, ie independent of the grid resolution) B-splines. Support is ( -(degree+1)/2, (degree+1)/2 )
  !function sll_b_spline( &
  !  degree, &
  !  x       &
  !  ) &
  !  result(res)
  !  sll_int32 :: degree
  !  sll_real64 :: x
  !  sll_real64 :: x_aux
  !  sll_real64 :: res
  !
  !  res = 0.0_f64
  !  if ( degree == 1 ) then
  !      ! pw affine spline
  !     if ( x <= -1 .or. x >= 1 ) then
  !          return
  !      end if
  !      if ( x < 0 ) then
  !          res = x+1.0_f64
  !          return
  !      end if
  !      res = 1.0_f64-x
  !      return
  !  end if
  !  if ( degree == 3 ) then
  !      x_aux = x+2
  !      if ( x_aux <= 0 .or. x_aux >= 4 ) then
  !          return
  !      end if
  !      if ( x_aux < 1 ) then
  !          res = 0.16666666666666666*(x_aux**3)
  !          return
  !      end if
  !      if ( x_aux < 2 ) then
  !          res = 0.6666666666666666 - 2.*x_aux + 2.*(x_aux**2) - 0.5*(x_aux**3)
  !          return
  !      end if
  !      if ( x_aux < 3 ) then
  !          res = -7.333333333333333 + 10.*x_aux - 4.*(x_aux**2) + 0.5*(x_aux**3)
  !          return
  !      end if
  !      res = 10.66666666666666 - 8.*x_aux + 2.*(x_aux**2) - 0.16666666666666666*(x_aux**3)
  !      return
  !  end if
  !  if ( degree == 5 ) then
  !      x_aux = x+3
  !      if ( x_aux <= 0 .or. x_aux >= 6 ) then
  !          res = 0.0_f64
  !          return
  !      end if
  !      if ( x_aux < 1 ) then
  !          res = 0.00833333333333333*(x_aux**5)
  !          return
  !      end if
  !      if ( x_aux < 2 ) then
  !          res = -0.0416666666667 *(x_aux**5) + 0.25 *(x_aux**4) - 0.5 *(x_aux**3) + 0.5 *(x_aux**2) - 0.25 *(x_aux) + 0.05
  !          return
  !      end if
  !      if ( x_aux < 3 ) then
  !          res = 0.0833333333333 *(x_aux**5) - 1.0 *(x_aux**4) + 4.5 *(x_aux**3) - 9.5 *(x_aux**2) + 9.75 *(x_aux) - 3.95
  !          return
  !      end if
  !      if ( x_aux < 4 ) then
  !          res = -0.0833333333333 *(x_aux**5) + 1.5 *(x_aux**4) - 10.5 *(x_aux**3) + 35.5 *(x_aux**2) - 57.75 *(x_aux) + 36.55
  !          return
  !      end if
  !      if ( x_aux < 5 ) then
  !          res = 0.0416666666667 *(x_aux**5) - 1.0 *(x_aux**4) + 9.5 *(x_aux**3) - 44.5 *(x_aux**2) + 102.25 *(x_aux) - 91.45
  !          return
  !      end if
  !      res = -0.008333333333333333 * ((x_aux-6.)**5)
  !      return
  !  end if
  !  print *, 'sll_b_spline(): ERROR, invalid value of argument degree ', degree
  !  STOP
  !  return
  !end function sll_b_spline

  subroutine sll_s_apply_periodic_bc_x( mesh, x)

    type(sll_t_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
  end subroutine sll_s_apply_periodic_bc_x


  subroutine sll_s_apply_periodic_bc_y( mesh, y)

    type(sll_t_cartesian_mesh_2d), pointer, intent(in) :: mesh
    sll_real64, intent(inout) :: y

    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine sll_s_apply_periodic_bc_y

  ! update the arrays closest_particle and closest_particle_distance with the index of the given particle
  ! if closer to what had been stored up to now.

  ! x_aux : x_particle - flow_grid_x_min   and  similarly for y, vx, vy
  subroutine sll_s_update_closest_particle_arrays(k_part,    &
                                            x_aux, y_aux, vx_aux, vy_aux,   &
                                            i, j, l, m,                     &
                                            h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy,   &
                                            closest_particle,               &
                                            closest_particle_distance)

      sll_int32, intent(in) :: k_part
      sll_int32, intent(in) :: i, j, l, m
      sll_real64, intent(in) :: x_aux, y_aux, vx_aux, vy_aux
      sll_real64, intent(in) :: h_virtual_cell_x, h_virtual_cell_y, h_virtual_cell_vx, h_virtual_cell_vy

      sll_int32, dimension(:,:,:,:), intent(inout)  :: closest_particle
      sll_real64, dimension(:,:,:,:), intent(inout) :: closest_particle_distance

      sll_real64 :: square_dist_to_cell_center

      square_dist_to_cell_center = (x_aux  - (i + 0.5) * h_virtual_cell_x )**2.    &
                                 + (y_aux  - (j + 0.5) * h_virtual_cell_y )**2.    &
                                 + (vx_aux - (l + 0.5) * h_virtual_cell_vx )**2.    &
                                 + (vy_aux - (m + 0.5) * h_virtual_cell_vy )**2.

      ! if new particle is closer to center, keep the new one
      if(closest_particle(i,j,l,m) == 0 .or. square_dist_to_cell_center < closest_particle_distance(i,j,l,m)) then
         closest_particle(i,j,l,m) = k_part
         closest_particle_distance(i,j,l,m) = square_dist_to_cell_center
      end if

  end subroutine sll_s_update_closest_particle_arrays

  function sll_f_eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, shift, x, y, vx, vy)
    sll_real64 :: x0,y0,vx0,vy0             ! centers of the hat
    sll_real64 :: r_x,r_y,r_vx,r_vy         ! radii of the hat
    sll_real64 :: basis_height, shift
    sll_real64 :: x, y, vx, vy
    sll_real64 :: sll_f_eval_hat_function

    sll_real64 :: inv_r_x
    sll_real64 :: inv_r_y
    sll_real64 :: inv_r_vx
    sll_real64 :: inv_r_vy

    inv_r_x  = 1./r_x
    inv_r_y  = 1./r_y
    inv_r_vx = 1./r_vx
    inv_r_vy = 1./r_vy

    sll_f_eval_hat_function = basis_height + shift * max(0._f64, 1. - inv_r_x*abs(x-x0) )             &
                                                 * max(0._f64, 1. - inv_r_y*abs(y-y0) )             &
                                                 * max(0._f64, 1. - inv_r_vx*abs(vx-vx0) )          &
                                                 * max(0._f64, 1. - inv_r_vy*abs(vy-vy0) )
  end function sll_f_eval_hat_function

end module  sll_m_particle_lbf_utilities
