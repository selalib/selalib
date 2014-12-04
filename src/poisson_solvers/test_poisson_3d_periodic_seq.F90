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

program test_poisson_3d_periodic_seq
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_poisson_solvers.h"
  use sll_constants
  use sll_poisson_3d_periodic_seq

  implicit none

  sll_int32  :: nx, ny, nz
  sll_real64 :: Lx, Ly, Lz
  sll_int32  :: error
  sll_real64 :: dx, dy, dz
  sll_real64, dimension(:,:,:), allocatable :: x, y, z
  sll_real64, dimension(:,:,:), allocatable :: rho, phi_an, phi
  sll_int32                                 :: i, j, k
  type (poisson_3d_periodic_plan_seq), pointer  :: plan
  sll_real64                                :: average_err
  sll_real64                                :: time_0, time_1, time_2
  sll_int32                                 :: i_test

  call cpu_time(time_0)
  nx = 128
  ny = 64
  nz = 32

  Lx = 2.0_f64*sll_pi
  Ly = 2.0_f64*sll_pi
  Lz = 2.0_f64*sll_pi

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  print*, 'Initialize Poisson 3D solver plan'
  SLL_ALLOCATE(x(nx,ny,nz),error)
  SLL_ALLOCATE(y(nx,ny,nz),error)
  SLL_ALLOCATE(z(nx,ny,nz),error)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           x(i,j,k) = (i-1)*dx
           y(i,j,k) = (j-1)*dy
           z(i,j,k) = (k-1)*dz
        end do
     end do
  end do

  !call write_mesh(x,y,z,nx,ny,nz,"mesh3d")

  SLL_ALLOCATE(rho(nx,ny,nz),error)
  SLL_ALLOCATE(phi(nx,ny,nz),error)
  SLL_ALLOCATE(phi_an(nx,ny,nz),error)

  plan => new_poisson_3d_periodic_plan_seq(nx, ny, nz, Lx, Ly, Lz)

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
        phi_an = (4/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) *  exp(-.5 & 
                 *(x-Lx/2)**2) * exp(-.5*(y-Ly/2)**2) * sin(z)
        rho    = phi_an * (3.0_f64 - ((x-Lx/2.0_f64)**2 + &
                 (y-Ly/2.0_f64)**2))
     end if

     call cpu_time(time_1)
     call solve_poisson_3d_periodic_seq(plan, rho, phi)
     call cpu_time(time_2)

     average_err = sum( abs(phi_an-phi) ) / (nx*ny*nz)
     print*, 'Average error:', average_err
     print*, 'dx*dy*dz =', dx*dy*dz
     print*, 'CPU time = ', time_2-time_1

     if ( average_err > dx*dx*dy) then
        print*, ' '
        print*, 'sll_poisson_3d_periodic_seq test not passed'
     end if


     !call write_vec1d(rho,nx,ny,nz,"rho"//char(i_test+48),"mesh3d",0)
     !call write_vec1d(phi,nx,ny,nz,"phi"//char(i_test+48),"mesh3d",0)

  end do

  call delete_poisson_3d_periodic_plan_seq(plan)


  print*, 'sll_poisson_3d_periodic_seq test: PASSED'

  call cpu_time(time_2)
  print*, 'Total CPU time : ', time_2-time_0


end program test_poisson_3d_periodic_seq
