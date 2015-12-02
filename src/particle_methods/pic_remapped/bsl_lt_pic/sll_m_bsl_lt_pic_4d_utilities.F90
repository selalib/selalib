!**************************************************************
!  Copyright INRIA
!  Authors : MCP,ALH
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

module sll_m_bsl_lt_pic_4d_utilities

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_m_accumulators.h"
!#include "particle_representation.h"

  use sll_m_constants, only: sll_pi
  use sll_m_cartesian_meshes
  use sll_m_working_precision

  use sll_m_timer
  use sll_m_utilities, only: sll_new_file_id, int2string
  use sll_m_gnuplot

implicit none 

contains

! <<get_particle_index_from_initial_position_on_cartesian_grid>>


  !! transforms a standard particle position (x,y) in (i_cell_x, i_cell_y, dx, dy)
  !! -> here the indices i_cell_x and i_cell_y do not need to be within [1, m2d%num_cells1] or [1, m2d%num_cells2]
  !!    so that in non-periodic domains we can track outside particles (markers)
  !!    (in periodic domains I don't think this is useful)
  subroutine global_to_cell_offset_extended( x, y, &
                      m2d,      &
                      i_cell_x, &
                      i_cell_y, &
                      offset_x, offset_y )

    sll_real64, intent(in)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(out) :: i_cell_x       !! not necessarily in [1, m2d%num_cells1], see comments above
    sll_int32,  intent(out) :: i_cell_y       !! not necessarily in [1, m2d%num_cells2], see comments above
    sll_real32, intent(out) :: offset_x, offset_y
    sll_real64              :: temp

    temp = (x - m2d%eta1_min) / m2d%delta_eta1
    i_cell_x  = 1 + int(floor(temp))
    offset_x = real(temp - real(i_cell_x - 1,f64), f32)

    temp = (y - m2d%eta2_min) / m2d%delta_eta2
    i_cell_y  = 1 + int(floor(temp))
    offset_y = real(temp - real(i_cell_y - 1,f64), f32)
    SLL_ASSERT(offset_x >= 0)
    SLL_ASSERT(offset_x <= 1 )
    SLL_ASSERT(offset_y >= 0)
    SLL_ASSERT(offset_y <= 1 )

    !! note: the (integer) index of the Poisson cell (within space computational domain) is then obtained with get_poisson_cell_index

  end subroutine global_to_cell_offset_extended

  !> <<cell_offset_to_global_extended>> performs the inverse transformation of global_to_cell_offset_extended above
  !!
  subroutine cell_offset_to_global_extended ( offset_x, offset_y, &
                                   i_cell_x, i_cell_y, m2d, &
                                   x, y )
  ! transforms sll_type of a particle (i_cell, dx, dy) into the standard
  ! particle position (x,y)
    sll_real64, intent(out)  :: x, y
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: i_cell_x, i_cell_y
    sll_real32, intent(in) :: offset_x, offset_y

    x = m2d%eta1_min + m2d%delta_eta1*( real(offset_x + i_cell_x-1, f64) )
    y = m2d%eta2_min + m2d%delta_eta2*( real(offset_y + i_cell_y-1, f64) )

  end subroutine cell_offset_to_global_extended


  subroutine get_poisson_cell_index( m2d, i_cell_x, i_cell_y, i_cell )
    type(sll_cartesian_mesh_2d), intent(in) :: m2d
    sll_int32,  intent(in) :: i_cell_x       !! not necessarily in [1, m2d%num_cells1]
    sll_int32,  intent(in) :: i_cell_y      !! not necessarily in [1, m2d%num_cells2]
    sll_int32,  intent(out) :: i_cell        !! in [1, m2d%num_cells1 * m2d%num_cells2]

    i_cell = 1 + modulo(i_cell_x-1,  m2d%num_cells1) + modulo(i_cell_y-1,  m2d%num_cells2) * m2d%num_cells1

    SLL_ASSERT( i_cell >= 1)
    SLL_ASSERT( i_cell <= m2d%num_cells1 * m2d%num_cells2 )

  end subroutine get_poisson_cell_index


  subroutine get_particle_index_from_initial_position_on_cartesian_grid (j_x, j_y, j_vx, j_vy,                            &
                                                                       n_parts_x, n_parts_y, n_parts_vx, n_parts_vy,    &
                                                                       k                                                &
                                                                       )
    sll_int32, intent(in) :: j_x
    sll_int32, intent(in) :: j_y
    sll_int32, intent(in) :: j_vx
    sll_int32, intent(in) :: j_vy
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy
    sll_int32, intent(out) :: k

    if(         j_x <= 0  .or. j_x > n_parts_x      &
        .or.    j_y <= 0  .or. j_y > n_parts_y      &
        .or.    j_vx <= 0 .or. j_vx > n_parts_vx    &
        .or.    j_vy <= 0 .or. j_vy > n_parts_vy  )then
        k = 0
    else
        k = 1+ (j_vy-1) + (j_vx-1) * n_parts_vy + (j_y-1) * n_parts_vy * n_parts_vx + (j_x-1) * n_parts_vy * n_parts_vx * n_parts_y
    end if

    SLL_ASSERT( k <= n_parts_x * n_parts_y * n_parts_vx * n_parts_vy )

  end subroutine get_particle_index_from_initial_position_on_cartesian_grid


  ! <<get_initial_position_on_cartesian_grid_from_particle_index>>

  subroutine get_initial_position_on_cartesian_grid_from_particle_index (k,                                               &
                                                                       n_parts_x, n_parts_y, n_parts_vx, n_parts_vy,    &
                                                                       j_x, j_y, j_vx, j_vy                             &
                                                                       )
    sll_int32, intent(in) :: k
    sll_int32, intent(in) :: n_parts_x
    sll_int32, intent(in) :: n_parts_y
    sll_int32, intent(in) :: n_parts_vx
    sll_int32, intent(in) :: n_parts_vy
    sll_int32, intent(out) :: j_x
    sll_int32, intent(out) :: j_y
    sll_int32, intent(out) :: j_vx
    sll_int32, intent(out) :: j_vy
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

    return
    SLL_ASSERT(n_parts_x>0)

  end subroutine


  ! <<sll_pic_shape>>
  !> sll_pic_shape(degree, x, y, vx, vy, inv_hx, inv_hy, inv_hvx, inv_hvy)
  !! computes the value of the 4d B-spline particle shape (no particle transformation here)
  !!
  !! note:
  !!  - inv_hx, inv_hy, inv_hvx, inv_hvy are the inverse
  !!    of the inter-particle distances (hx, hy, hvx, hvy) on the initial (or remapping) grid
  !!
  !!  - x, y, vx, vy are phase-space coordinates relative to the particle center
  !!    and they are not scaled: for instance if degree = 1 (resp. if degree = 3),
  !!    and |x| >= hx (resp. |x| >= 2*hx), then sll_pic_shape returns 0
  !!    (and similarly for y, vx, vy)

  function sll_pic_shape( &
    degree,             &
    x, y, vx, vy,       &
    inv_hx, inv_hy, inv_hvx, inv_hvy   &
    ) &
    result(shape)

    sll_int32 :: degree
    sll_real64 :: x, y, vx, vy
    sll_real64 :: inv_hx, inv_hy, inv_hvx, inv_hvy
    sll_real64 :: shape
!    sll_int32 :: nc_eta2
!    sll_int32 :: ierr

    ! uses [[sll_b_spline]]
    shape =   inv_hx  * sll_b_spline(degree,inv_hx  * x ) &
            * inv_hy  * sll_b_spline(degree,inv_hy  * y ) &
            * inv_hvx * sll_b_spline(degree,inv_hvx * vx) &
            * inv_hvy * sll_b_spline(degree,inv_hvy * vy)
  end function sll_pic_shape


  !> added by MCP
  !! 4d-reference B-spline shape function (independent of the grid resolution), function support is [-part_cp, part_cp]^4
  !! with part_cp = (degree +1)/2
  function sll_ref_pic_shape( &
    degree,     &
    x,y,vx,vy &
    ) &
    result(ref_shape)

    sll_int32 :: degree
    sll_real64 :: x, y, vx, vy
    sll_real64 :: ref_shape
!    sll_int32 :: nc_eta2
!    sll_int32 :: ierr

    ref_shape = sll_b_spline(degree,x) * sll_b_spline(degree,y) * sll_b_spline(degree,vx) * sll_b_spline(degree, vy)
  end function sll_ref_pic_shape


  ! added by MCP
  ! <<sll_b_spline>>
  ! univariate centered (and reference, ie independent of the grid resolution) B-splines. Support is ( -(degree+1)/2, (degree+1)/2 )
  function sll_b_spline( &
    degree, &
    x       &
    ) &
    result(res)
    sll_int32 :: degree
    sll_real64 :: x
    sll_real64 :: x_aux
    sll_real64 :: res

    res = 0.0_f64
    if ( degree == 1 ) then
        ! pw affine spline
       if ( x <= -1 .or. x >= 1 ) then
            return
        end if
        if ( x < 0 ) then
            res = x+1.0_f64
            return
        end if
        res = 1.0_f64-x
        return
    end if
    if ( degree == 3 ) then
        x_aux = x+2
        if ( x_aux <= 0 .or. x_aux >= 4 ) then
            return
        end if
        if ( x_aux < 1 ) then
            res = 0.16666666666666666*(x_aux**3)
            return
        end if
        if ( x_aux < 2 ) then
            res = 0.6666666666666666 - 2.*x_aux + 2.*(x_aux**2) - 0.5*(x_aux**3)
            return
        end if
        if ( x_aux < 3 ) then
            res = -7.333333333333333 + 10.*x_aux - 4.*(x_aux**2) + 0.5*(x_aux**3)
            return
        end if
        res = 10.66666666666666 - 8.*x_aux + 2.*(x_aux**2) - 0.16666666666666666*(x_aux**3)
        return
    end if
    if ( degree == 5 ) then
        x_aux = x+3
        if ( x_aux <= 0 .or. x_aux >= 6 ) then
            res = 0.0_f64
            return
        end if
        if ( x_aux < 1 ) then
            res = 0.00833333333333333*(x_aux**5)
            return
        end if
        if ( x_aux < 2 ) then
            res = -0.0416666666667 *(x_aux**5) + 0.25 *(x_aux**4) - 0.5 *(x_aux**3) + 0.5 *(x_aux**2) - 0.25 *(x_aux) + 0.05
            return
        end if
        if ( x_aux < 3 ) then
            res = 0.0833333333333 *(x_aux**5) - 1.0 *(x_aux**4) + 4.5 *(x_aux**3) - 9.5 *(x_aux**2) + 9.75 *(x_aux) - 3.95
            return
        end if
        if ( x_aux < 4 ) then
            res = -0.0833333333333 *(x_aux**5) + 1.5 *(x_aux**4) - 10.5 *(x_aux**3) + 35.5 *(x_aux**2) - 57.75 *(x_aux) + 36.55
            return
        end if
        if ( x_aux < 5 ) then
            res = 0.0416666666667 *(x_aux**5) - 1.0 *(x_aux**4) + 9.5 *(x_aux**3) - 44.5 *(x_aux**2) + 102.25 *(x_aux) - 91.45
            return
        end if
        res = -0.008333333333333333 * ((x_aux-6.)**5)
        return
    end if
    print *, 'sll_b_spline(): ERROR, invalid value of argument degree ', degree
    STOP
    return
  end function sll_b_spline

  subroutine apply_periodic_bc_x( mesh, x)

!    use sll_m_cartesian_meshes
    ! [[file:../working_precision/sll_m_working_precision.h]]
!    use sll_m_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
  end subroutine apply_periodic_bc_x


  subroutine apply_periodic_bc_y( mesh, y)

!    use sll_m_cartesian_meshes
    ! [[file:../working_precision/sll_m_working_precision.h]]
!    use sll_m_working_precision

    type(sll_cartesian_mesh_2d), pointer, intent(in) :: mesh
    sll_real64, intent(inout) :: y

    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc_y

  ! update the arrays closest_particle and closest_particle_distance with the index of the given particle
  ! if closer to what had been stored up to now.

  ! x_aux : x_particle - x_min_virtual_mesh   and  similarly for y, vx, vy
  subroutine update_closest_particle_arrays(k_part,                         &
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

  end subroutine update_closest_particle_arrays

  function eval_landau_fx(alpha, kx, x)
    sll_real64 :: alpha, kx, x
    sll_real64 :: eval_landau_fx
    eval_landau_fx = 1._f64 + alpha * cos(kx * x)
  end function eval_landau_fx

  function eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, basis_height, hat_shift, x, y, vx, vy)
    sll_real64 :: x0,y0,vx0,vy0             ! centers of the hat
    sll_real64 :: r_x,r_y,r_vx,r_vy         ! radii of the hat
    sll_real64 :: basis_height, hat_shift
    sll_real64 :: x, y, vx, vy
    sll_real64 :: eval_hat_function

    sll_real64 :: inv_r_x
    sll_real64 :: inv_r_y
    sll_real64 :: inv_r_vx
    sll_real64 :: inv_r_vy

    inv_r_x  = 1./r_x
    inv_r_y  = 1./r_y
    inv_r_vx = 1./r_vx
    inv_r_vy = 1./r_vy

    eval_hat_function = basis_height + hat_shift * max(0._f64, 1. - inv_r_x*abs(x-x0) )             &
                                                 * max(0._f64, 1. - inv_r_y*abs(y-y0) )             &
                                                 * max(0._f64, 1. - inv_r_vx*abs(vx-vx0) )          &
                                                 * max(0._f64, 1. - inv_r_vy*abs(vy-vy0) )
  end function eval_hat_function

end module  sll_m_bsl_lt_pic_4d_utilities
