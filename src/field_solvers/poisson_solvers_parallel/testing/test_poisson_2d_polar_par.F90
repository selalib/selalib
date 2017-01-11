program test_poisson_polar_parallel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

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

  use sll_m_poisson_2d_polar_par, only: &
    sll_t_poisson_2d_polar_par, &
    sll_s_poisson_2d_polar_par_init, &
    sll_s_poisson_2d_polar_par_solve

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_view_lims

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_poisson_2d_polar_par) :: poisson
  sll_real64, allocatable :: rhs(:,:)
  sll_real64, allocatable :: phi(:,:)
  sll_real64, allocatable :: phi_cos(:,:)
  sll_real64, allocatable :: phi_sin(:,:)
  sll_real64, allocatable :: r(:)
  sll_real64, allocatable :: a(:)
  sll_real64, allocatable :: x(:,:)
  sll_real64, allocatable :: y(:,:)

  type(sll_t_layout_2d), pointer :: layout_r ! sequential in r direction
  type(sll_t_layout_2d), pointer :: layout_a ! sequential in theta direction

  sll_int32 :: global(2)
  sll_int32 :: gi, gj
  sll_int32 :: i, j
  sll_int32 :: nc_r, nr
  sll_int32 :: nc_a, na
  sll_real64 :: rmin, rmax, delta_r
  sll_real64 :: amin, amax, delta_a
  sll_int32 :: error
  sll_int32 :: psize
  sll_int32 :: prank
  sll_int32 :: nr_loc
  sll_int32 :: na_loc

  sll_int32, parameter  :: n = 4

  write(*,'(a)') 'Testing the 2D Poisson solver in polar coordinates'

  rmin   = 1.0_f64
  rmax   = 2.0_f64

  amin   = 0.0_f64
  amax   = 2.0_f64 * sll_p_pi

  nc_r    = 32 !256
  nc_a    = 64 !1024
  nr      = nc_r+1
  na      = nc_a+1
  delta_r = (rmax-rmin)/real( nr-1, f64 )
  delta_a = 2.0_f64*sll_p_pi/real( na-1, f64 )

  !Boot parallel environment
  call sll_s_boot_collective()

  psize = sll_f_get_collective_size(sll_v_world_collective)
  prank = sll_f_get_collective_rank(sll_v_world_collective)

  layout_r => sll_f_new_layout_2d( sll_v_world_collective )
  layout_a => sll_f_new_layout_2d( sll_v_world_collective )

  call sll_o_initialize_layout_with_distributed_array( nr,na,1,psize,layout_r )

  call sll_o_initialize_layout_with_distributed_array( nr,na,psize,1,layout_a )

  flush( output_unit )

  if (prank == 0) then
    call sll_o_view_lims(layout_a)
    call sll_o_view_lims(layout_a)
  end if

  flush( output_unit )

  call sll_o_compute_local_sizes( layout_a, nr_loc, na_loc )

  allocate( rhs(1:nr_loc,1:na_loc) )
  allocate( phi(1:nr_loc,1:na_loc) )
  allocate( phi_cos(1:nr_loc,1:na_loc) )
  allocate( phi_sin(1:nr_loc,1:na_loc) )
  allocate( r(1:nr_loc) )
  allocate( a(1:na_loc) )
  allocate( x(1:nr_loc,1:na_loc) )
  allocate( y(1:nr_loc,1:na_loc) )

  do i = 1, nr_loc
    global = sll_o_local_to_global( layout_a, (/i, 1/))
    gi = global(1)
    r(i) = rmin+(gi-1)*delta_r
  end do

  do j = 1, na_loc
    global = sll_o_local_to_global( layout_a, (/1, j/))
    gj = global(2)
    a(j) = (gj-1)*delta_a
  end do

  do j=1,na_loc
     do i=1,nr_loc
       phi_cos(i,j) = (r(i)-rmin)*(r(i)-rmax)*cos(n*a(j))*r(i)
       phi_sin(i,j) = (r(i)-rmin)*(r(i)-rmax)*sin(n*a(j))*r(i)
       x(i,j) = r(i)*cos(a(j))
       y(i,j) = r(i)*sin(a(j))
     end do
  end do

  call sll_s_poisson_2d_polar_par_init( solver = poisson, &
       layout_r = layout_r, &
       layout_a = layout_a, &
       rmin = rmin, &
       rmax = rmax, &
       nr = nc_r, &
       ntheta = nc_a, &
       bc_rmin = sll_p_dirichlet, &
       bc_rmax = sll_p_dirichlet )

  do i= 1, nr_loc
     do j = 1, na_loc
       rhs(i,j) = -rhs_sin( r(i), a(j) )
     end do
  end do

  call sll_s_poisson_2d_polar_par_solve(poisson, rhs, phi)

  call sll_o_gnuplot_2d_parallel( x, y, phi_sin, 'phi_sin',  1, error )
  call sll_o_gnuplot_2d_parallel( x, y, phi, 'solution', 1, error )

  call error_max( phi, phi_sin, 1e-4_f64 )
  call sll_s_halt_collective()

contains

  subroutine error_max(phi, phi_exact, tolmax)

    sll_real64, intent(in), dimension(:,:) :: phi
    sll_real64, intent(in), dimension(:,:) :: phi_exact
    sll_real64 :: errmax, tolmax

    errmax = maxval(abs(phi-phi_exact))
    write(*,'(a,e10.3)') 'maximum error = ', errmax
    if ( errmax > tolmax ) then
      write(*,'(a)') 'FAILED'
    else
      write(*,'(a)') 'PASSED'
    end if

  end subroutine error_max


  sll_real64 function rhs_cos( r, theta )

    ! Assume r>=1 and r<=2
    ! phi = r*(r-rmin)*(r-rmax)*cos(n*theta)
    ! return: d2phi_dr + d1phi_dr/r + d2phi_dtheta/r^2

    sll_real64 :: r, theta
    sll_real64 :: d1phi_dr, d2phi_dr, d2phi_dtheta

    d1phi_dr = cos(n*theta)*(3.0_f64*r**2-2.0_f64*(rmin+rmax)*r+rmin*rmax)
    d2phi_dr = cos(n*theta)*(6.0_f64*r-2.0_f64*(rmin+rmax))
    d2phi_dtheta = -r*(r-rmin)*(r-rmax)*n**2*cos(n*theta)

    rhs_cos = d2phi_dr + d1phi_dr/r + d2phi_dtheta/r**2

!    rhs_cos = - (r-rmax)*(r-rmin)*n**2*cos(n*theta)/r &
!              + ((r-rmax)*(r-rmin)*cos(n*theta)  &
!              + (r-rmax)*r*cos(n*theta) + (r-rmin)*r*cos(n*theta) &
!              + 2.0_f64*((r-rmax)*cos(n*theta) + (r-rmin)*cos(n*theta) &
!              + r*cos(n*theta))*r)/r

  end function rhs_cos

  sll_real64 function rhs_sin( r, theta)

    ! Assume r>=1 and r<=2
    ! phi = r*(r-rmin)*(r-rmax)*sin(n*theta)
    ! return: d2phi_dr + d1phi_dr/r + d2phi_dtheta/r^2

    sll_real64 :: r, theta
    sll_real64 :: d1phi_dr, d2phi_dr, d2phi_dtheta

    d1phi_dr = sin(n*theta)*(3.0_f64*r**2-2.0_f64*(rmin+rmax)*r+rmin*rmax)
    d2phi_dr = sin(n*theta)*(6.0_f64*r-2.0_f64*(rmin+rmax))
    d2phi_dtheta = -r*(r-rmin)*(r-rmax)*n**2*sin(n*theta)

    rhs_sin = d2phi_dr + d1phi_dr/r + d2phi_dtheta/r**2

  end function rhs_sin

end program test_poisson_polar_parallel
