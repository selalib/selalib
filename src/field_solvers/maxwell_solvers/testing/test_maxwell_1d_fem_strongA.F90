!>$L_x$ domain dimensions and M is an integer.
!>$
!>$
!>$
!>B_z(x,y,t) =   \cos(\frac{2 M \pi}{L_x} x)  \cos(\frac{2 M \pi}{L_x} t)
!>$
!>$
!>E_y(x,y,t) = \sin(\frac{2 M \pi}{L_x} x)  \sin(\frac{2 M \pi}{L_x} t)
!>$
!
!  Contact : Katharina Kormann
!
program test_maxwell_1d_fem_strongA
  !------------------------------------------------------------------------
  !  test 1D Maxwell spline finite element solver on a periodic grid with strong Ampere
  !------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_eval_uniform_periodic_spline_curve

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_maxwell_1d_base, only: &
    sll_s_plot_two_fields_1d

  use sll_m_maxwell_1d_fem, only: &
    sll_t_maxwell_1d_fem

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_bsplines), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta1

  sll_int32  :: nc_eta1
  sll_int32  :: error
  sll_real64 :: tol

  type(sll_t_maxwell_1d_fem)  :: maxwell_1d

  sll_real64, dimension(:), allocatable :: ex
  sll_real64, dimension(:), allocatable :: ey
  sll_real64, dimension(:), allocatable :: bz
  sll_real64, dimension(:), allocatable :: bz_exact
  sll_real64, dimension(:), allocatable :: ex_exact
  sll_real64, dimension(:), allocatable :: ey_exact
  sll_real64, dimension(:), allocatable :: rho
  sll_real64, dimension(:), allocatable :: sval

  sll_int32                               :: i!, j
  sll_real64                              :: time
  sll_int32                               :: istep, nstep
  sll_real64                              :: err_ex
  sll_real64                              :: err_ex2
  sll_real64                              :: err_ey
  sll_real64                              :: err_bz
  sll_real64                              :: err_l2norm
  sll_real64                              :: dt
  !sll_real64                              :: cfl = 0.5_f64
  sll_real64                              :: Lx
  sll_real64                              :: xi
  sll_real64, dimension(2)                :: domain
  sll_int32                               :: deg
  sll_int32,  parameter                   :: mode = 2
  sll_real64                              :: l2norm

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta1 = 128
  Lx = eta1_max-eta1_min
  delta_eta1 = Lx/real(nc_eta1,f64)
  domain = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3

  ! Initialise maxwell FEM object
  call maxwell_1d%init(domain, nc_eta1, deg, 0.1_f64, strong_ampere=.true.)

  ! Allocate arrays
  SLL_CLEAR_ALLOCATE(ex(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(bz(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(bz_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ex_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(ey_exact(1:nc_eta1),error)
  SLL_CLEAR_ALLOCATE(rho(1:nc_eta1),error)

  SLL_CLEAR_ALLOCATE(sval(1:nc_eta1),error)


  ! Test Poisson
  !-------------
  ! Set exact solution
  do i = 1, nc_eta1
     xi = eta1_min + real(i-1,f64)*delta_eta1
     ex_exact(i) =   sin_k(xi)/(2.0_f64*mode*sll_p_pi/Lx)
     !rho(i) = cos_k( xi )
     rho(i) = (sin_k( xi ) - sin_k( xi - delta_eta1))/(2.0_f64*mode*sll_p_pi/Lx)/delta_eta1
  end do
  
  !call maxwell_1d%compute_rhs_from_function( cos_k, deg, rho)

  call maxwell_1d%compute_e_from_rho( rho, ex ) 

  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_eval_uniform_periodic_spline_curve(deg, ex, sval)
  err_ex = maxval(abs(sval-ex_exact))
  print*, 'error Poisson',  err_ex

  
  !call sll_s_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)
 
  ! Test Ampere
  !-------------
  ! Set time step
  dt = .5_f64 * delta_eta1
  ! Set exact solution
  do i = 1, nc_eta1
     xi = eta1_min + real(i-1,f64)*delta_eta1
     ex_exact(i) =   -cos_k(xi)*dt
  end do

  !call maxwell_1d%compute_rhs_from_function(cos_k, deg, rho)
  ex = 0.0_f64
  call maxwell_1d%compute_E_from_j(-ex_exact, 1, ex)!dt*rho, 1, ex )

  ! Evaluate spline curve at grid points and compute error
  ! Ex is a 1-form, i.e. one spline degree lower
  call sll_s_eval_uniform_periodic_spline_curve(deg, ex, sval)
  err_ex2 = maxval(abs(sval-ex_exact))
  write(12,*) sval
  write(13,*) ex_exact
  write(14,*) ex
  print*, 'error Ampere',  err_ex2
  !call sll_plot_two_fields_1d('ex',nc_eta1,sval,ex_exact,0,0.0_f64)

  l2norm =  maxwell_1d%l2norm_squared(ex,deg)
  err_l2norm = abs(l2norm - dt**2*sll_p_pi)
  print*, 'error l2 norm', err_l2norm


  ! Test Maxwell on By and Ez 
  !--------------------------
  ! Set time stepping parameters
  time  = 0.0_f64
  dt    = .5_f64 * delta_eta1
  nstep = 10

  ! Compute initial fields 
  ex = 0.0_f64 ! 0-form -> splines of degree deg
  ey = 0.0_f64
  call maxwell_1d%L2projection(cos_k, deg, bz) ! 1-form -> splines of degree deg

  ! Time loop. Second order Strang splitting
  do istep = 1, nstep 
     call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     call maxwell_1d%compute_e_from_b(         dt, bz, ey)
     call maxwell_1d%compute_b_from_e( 0.5_f64*dt, ey, bz)
     
     time = time + dt

     do i = 1, nc_eta1
        xi = eta1_min + real(i-1,f64)*delta_eta1
        ey_exact(i) =   sin(mode*2*sll_p_pi*xi/Lx) * sin(mode*2*sll_p_pi*time/Lx)
        bz_exact(i) =   cos(mode*2*sll_p_pi*xi/Lx) * cos(mode*2*sll_p_pi*time/Lx)
     end do

     call sll_s_eval_uniform_periodic_spline_curve(deg-1, ey, sval)
     err_ey = maxval(abs(sval-ey_exact))
     call sll_s_eval_uniform_periodic_spline_curve(deg, bz, sval)
     err_bz = maxval(abs(sval-bz_exact))

     write(*,"(10x,' istep = ',I6)",advance="no") istep
     write(*,"(' time = ',g12.3,' sec')",advance="no") time
     write(*,"(' erreur L2 = ',2g15.5)") err_ey, err_bz

     !call sll_s_plot_two_fields_1d('bz',nc_eta1,sval,bz_exact,istep,time)

  end do ! next time step

  tol = 1.0d-3
  if ((err_bz < tol) .and. (err_ey < tol) .and. (err_ex < 1.0d-14) .and. (err_ex2 < 1.0d-16)) then! .and. (err_l2norm < 1.0d-13)) then
     print*,'PASSED'
  endif

  call maxwell_1d%free()

  DEALLOCATE(ex)
  DEALLOCATE(ey)
  DEALLOCATE(bz)
  DEALLOCATE(bz_exact)
  DEALLOCATE(ex_exact)
  DEALLOCATE(ey_exact)
  DEALLOCATE(rho)

contains

  function cos_k(x)

    sll_real64             :: cos_k
    sll_real64, intent(in) :: x

    cos_k = cos(mode*2*sll_p_pi*x/Lx) 

  end function cos_k

  function sin_k(x)

    sll_real64             :: sin_k
    sll_real64, intent(in) :: x

    sin_k = sin(mode*2*sll_p_pi*x/Lx) 

  end function sin_k

end program test_maxwell_1d_fem_strongA
