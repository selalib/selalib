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

program unit_test_poisson_2d_polar_parallel_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gnuplot_parallel, only: &
    sll_o_gnuplot_2d_parallel

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

  use sll_m_poisson_2d_polar_parallel_solver, only: &
    sll_f_new_poisson_2d_polar_parallel_solver

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  class(sll_c_poisson_2d_base), pointer :: poisson 
  sll_int32 :: psize
  sll_int32 :: prank
  type(sll_t_layout_2d), pointer :: layout_r ! sequential in r direction
  type(sll_t_layout_2d), pointer :: layout_a ! sequential in theta direction
  sll_int32 :: nr_loc
  sll_int32 :: na_loc
  sll_real64, dimension(:,:), allocatable :: rhs
  sll_real64, dimension(:,:), allocatable :: phi
  sll_real64, dimension(:,:), allocatable :: phi_cos
  sll_real64, dimension(:,:), allocatable :: phi_sin
  sll_real64, dimension(:),   allocatable :: r
  sll_real64, dimension(:),   allocatable :: a
  sll_real64, dimension(:,:), allocatable :: x
  sll_real64, dimension(:,:), allocatable :: y


  sll_int32, dimension(2)   :: global
  sll_int32                 :: gi, gj
  sll_int32  :: i, j
  sll_int32  :: nc_r, nr
  sll_int32  :: nc_a, na
  sll_real64 :: r_min, r_max, delta_r
  sll_real64 :: a_min, a_max, delta_a
  sll_int32  :: error

  sll_int32, parameter  :: n = 4

  print*,'#Testing the Poisson solver in 2D, polar coordinate'

  r_min   = 1.0_f64
  r_max   = 2.0_f64

  a_min   = 0.0_f64
  a_max   = 2.0_f64 * sll_p_pi

  nc_r    = 32 !256
  nc_a    = 64 !1024
  nr      = nc_r+1
  na      = nc_a+1
  delta_r = (r_max-r_min)/real(nr-1,f64)
  delta_a = 2.0_f64*sll_p_pi/real(na-1,f64)

  
  
  
  
  

  !Boot parallel environment
  call sll_s_boot_collective()


  psize  = sll_f_get_collective_size(sll_v_world_collective)
  prank  = sll_f_get_collective_rank(sll_v_world_collective)

  layout_r => sll_f_new_layout_2d( sll_v_world_collective )
  layout_a => sll_f_new_layout_2d( sll_v_world_collective )

  call sll_o_initialize_layout_with_distributed_array( &
    nr, &
    na, &
    1, &
    psize, &
    layout_r)

  call sll_o_initialize_layout_with_distributed_array( &
    nr, &
    na, &
    psize, &
    1, &
    layout_a)


  call sll_o_compute_local_sizes(layout_a, nr_loc, na_loc )


  SLL_CLEAR_ALLOCATE(rhs(1:nr_loc,1:na_loc),error)
  SLL_CLEAR_ALLOCATE(phi(1:nr_loc,1:na_loc),error)
  SLL_CLEAR_ALLOCATE(phi_cos(1:nr_loc,1:na_loc),error)
  SLL_CLEAR_ALLOCATE(phi_sin(1:nr_loc,1:na_loc),error)
  SLL_CLEAR_ALLOCATE(r(1:nr_loc),error)
  SLL_CLEAR_ALLOCATE(a(1:na_loc),error)
  SLL_CLEAR_ALLOCATE(x(1:nr_loc,1:na_loc),error)
  SLL_CLEAR_ALLOCATE(y(1:nr_loc,1:na_loc),error)

  do i = 1, nr_loc
    global = sll_o_local_to_global( layout_a, (/i, 1/))
    gi = global(1)
    r(i)=r_min+(gi-1)*delta_r
  end do

  do j = 1, na_loc
    global = sll_o_local_to_global( layout_a, (/1, j/))
    gj = global(2)
    a(j)=(gj-1)*delta_a
  end do

  do j=1,na_loc
    do i=1,nr_loc
      phi_cos(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*a(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*sin(n*a(j))*r(i)
      x(i,j)   = r(i)*cos(a(j))
      y(i,j)   = r(i)*sin(a(j))
    end do
  end do







  
  poisson =>sll_f_new_poisson_2d_polar_parallel_solver( &
    layout_r, &
    layout_a, &
    r_min, &
    r_max, &
    nc_r, &
    nc_a, &
   (/sll_p_dirichlet, sll_p_dirichlet/))

  do i =1,nr_loc
    do j=1,na_loc
      rhs(i,j) = - f_sin(r(i), a(j))
      !rhs(i,j) = - f_cos(r(i), a(j))
    end do
  end do


  
  call poisson%compute_phi_from_rho( phi, rhs )

  print *,maxval(phi),minval(phi)
  

  call sll_o_gnuplot_2d_parallel(x, y, phi_sin, 'phi_sin',  1, error)
  call sll_o_gnuplot_2d_parallel(x, y, phi,     'solution', 1, error)

  call error_max(phi_sin,phi,1e-4_f64)
  call sll_s_halt_collective()

  contains

  subroutine error_max(phi, phi_exact, tolmax)

    sll_real64, intent(in), dimension(:,:) :: phi
    sll_real64, intent(in), dimension(:,:) :: phi_exact
    sll_real64 :: errmax, tolmax 

    errmax = maxval(abs(phi_exact-phi))
    write(*,201) errmax
    if ( errmax > tolmax ) then
      print*,'FAILED'
    else
      print*,'PASSED'
    end if

    201 format(' maximum error  =  ',e10.3)

  end subroutine error_max


  sll_real64 function f_cos( r, theta )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta

   f_cos = -(r-r_max)*(r-r_min)*n*n*cos(n*theta)/r &
           + ((r-r_max)*(r-r_min)*cos(n*theta)  &
           + (r-r_max)*r*cos(n*theta) + (r-r_min)*r*cos(n*theta) &
           + 2.0_f64*((r-r_max)*cos(n*theta) + (r-r_min)*cos(n*theta) &
           + r*cos(n*theta))*r)/r


  end function f_cos

  sll_real64 function f_sin( r, theta)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta
   
   f_sin = -(r-r_max)*(r-r_min)*n*n*sin(n*theta)/r &
         + ((r-r_max)*(r-r_min)*sin(n*theta) &
         + (r-r_max)*r*sin(n*theta) + (r-r_min)*r*sin(n*theta) &
         + 2.0_f64*((r-r_max)*sin(n*theta) + (r-r_min)*sin(n*theta)  &
         + r*sin(n*theta))*r)/r

  end function f_sin


end program
