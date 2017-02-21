!***************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d.F90
!
!> @brief 
!> Selalib poisson solvers (3D) unit test
!> Last modification: March 22, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!                                  
!***************************************************************************

program test_poisson_3d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_3d_periodic, only: &
    sll_t_poisson_3d_periodic, &
    sll_s_poisson_3d_periodic_init, &
    sll_s_poisson_3d_periodic_solve, &
    sll_s_poisson_3d_periodic_free

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32  :: nx, ny, nz
  sll_real64 :: Lx, Ly, Lz
  sll_int32  :: error
  sll_real64 :: dx, dy, dz
  sll_real64, dimension(:,:,:), allocatable :: x, y, z
  sll_real64, dimension(:,:,:), allocatable :: rho, phi_an, phi
  sll_int32                                 :: i, j, k
  type (sll_t_poisson_3d_periodic)          :: poisson
  sll_real64                                :: average_err
  sll_real64                                :: time_0, time_1, time_2
  sll_int32                                 :: i_test

  call cpu_time(time_0)
  nx = 128
  ny = 64
  nz = 32

  Lx = 2.0_f64*sll_p_pi
  Ly = 2.0_f64*sll_p_pi
  Lz = 2.0_f64*sll_p_pi

  dx = Lx/real(nx,f64)
  dy = Ly/real(ny,f64)
  dz = Lz/real(nz,f64)

  print*, 'Initialize Poisson 3D solver '
  SLL_ALLOCATE(x(nx,ny,nz),error)
  SLL_ALLOCATE(y(nx,ny,nz),error)
  SLL_ALLOCATE(z(nx,ny,nz),error)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           x(i,j,k) = real(i-1,f64)*dx
           y(i,j,k) = real(j-1,f64)*dy
           z(i,j,k) = real(k-1,f64)*dz
        end do
     end do
  end do

  !call write_mesh(x,y,z,nx,ny,nz,"mesh3d")

  SLL_ALLOCATE(rho(nx,ny,nz),error)
  SLL_ALLOCATE(phi(nx,ny,nz),error)
  SLL_ALLOCATE(phi_an(nx,ny,nz),error)

  call sll_s_poisson_3d_periodic_init(poisson, nx, ny, nz, Lx, Ly, Lz)

  print*, ' '
  call cpu_time(time_1)
  print*,' CPU time = ', time_1-time_0

  do i_test = 1, 2

     print*, '----------------------------------------'
     print*, 'poisson_3d in sequential test ', i_test
     print*, '----------------------------------------'

     if (i_test == 1) then
        phi_an = cos(x)*sin(y)*cos(z)
        rho = 3._f64 * phi_an
     else if (i_test == 2) then
        phi_an = (4.0_f64/(sll_p_pi*sqrt(sll_p_pi)*Lx*Ly*Lz)) *  exp(-.5_f64 & 
                 *(x-Lx/2)**2) * exp(-.5_f64*(y-Ly/2.0_f64)**2) * sin(z)
        rho    = phi_an * (3.0_f64 - ((x-Lx/2.0_f64)**2 + &
                 (y-Ly/2.0_f64)**2))
     end if

     call cpu_time(time_1)
     call sll_s_poisson_3d_periodic_solve(poisson, rho, phi)
     call cpu_time(time_2)

     average_err = sum( abs(phi_an-phi) ) / real(nx*ny*nz,f64)
     print*, 'Average error:', average_err
     print*, 'dx*dy*dz =', dx*dy*dz
     print*, 'CPU time = ', time_2-time_1

     if ( average_err > dx*dx*dy) then
        print*, 'i_test = ', i_test
        stop 'sll_m_poisson_3d_periodic test not passed'
     else
        print*, 'sll_m_poisson_3d_periodic test: PASSED'
     end if

  end do

  call sll_s_poisson_3d_periodic_free(poisson)

  call cpu_time(time_2)
  print*, 'Total CPU time : ', time_2-time_0


end program test_poisson_3d_periodic
