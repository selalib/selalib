!***************************************************************************
!
! Selalib 2012     
! File : test_poisson_3d.F90
!
!> @brief 
!> Selalib poisson solvers (3D) unit test
!> Last modification: March 06, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!> Pierre NAVARO (navaro@math.unistra.fr)
!                                  
!***************************************************************************

program test_poisson_3d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_poisson_3d_periodic_util
  use sll_poisson_3d_periodic_seq

  implicit none

  sll_int32  :: nx, ny, nz
  sll_real64 :: Lx, Ly, Lz
  sll_int32  :: ierr
  sll_real64 :: dx, dy, dz
  sll_real64                                :: x, y, z
  sll_real64, dimension(:,:,:), allocatable :: rho, phi_an, phi
  sll_int32                                 :: i, j, k
  type (poisson_3d_periodic_plan), pointer  :: plan
  sll_real64                                :: average_err
  sll_real64                                :: time_0, time_1, time_2
  sll_int32                                 :: i_test


  call cpu_time(time_0)
  nx = 128
  ny = 128
  nz = 128

  Lx = 2*sll_pi
  Ly = 2*sll_pi
  Lz = 2*sll_pi

  dx = Lx/nx
  dy = Ly/ny
  dz = Lz/nz

  print*, 'Initialize Poisson 3D solver plan'
  SLL_ALLOCATE(rho(nx,ny,nz),ierr)
  SLL_ALLOCATE(phi(nx,ny,nz),ierr)
  SLL_ALLOCATE(phi_an(nx,ny,nz),ierr)

  plan => new_poisson_3d_periodic_plan(cmplx(rho, 0_f64, kind=f64), Lx, Ly, Lz)

  print*, ' '
  call cpu_time(time_1)
  print*,' CPU time = ', time_1-time_0


  do i_test = 1, 2

     print*, 'poisson_3d in sequential test ', i_test
     print*, '----------------------------------------'
     do k=1,nz
        z = (k-1)*dz
        do j=1,ny
           y = (j-1)*dy
           do i=1,nx
              x = (i-1)*dx
              if (i_test == 1) then
                 phi_an(i,j,k) = cos(x)*sin(y)*cos(z)
                 rho(i,j,k) = 3._f64 * phi_an(i,j,k)
              else if (i_test == 2) then
                 phi_an(i,j,k) = (4/(sll_pi*sqrt(sll_pi)*Lx*Ly*Lz)) * &
                               exp(-.5*(x-Lx/2)**2) * exp(-.5*(y-Ly/2)**2) * sin(z)
                 rho(i,j,k) = phi_an(i,j,k) * &
                              (3.0_f64 - ((x-Lx/2.0_f64)**2 + (y-Ly/2.0_f64)**2))
              end if
           end do
        end do
     end do
   
     call cpu_time(time_1)
     call solve_poisson_3d_periodic_seq(plan, rho, phi)
     call cpu_time(time_2)
   
     average_err = sum( abs(phi_an-phi) ) / (nx*ny*nz)
     print*, 'Average error:', average_err
     print*, 'dx*dy*dz =', dx*dy*dz
     print*, 'CPU time = ', time_2-time_1
   
     if ( average_err > dx*dx*dy) then
        print*, ' '
        print*, 'Test stoppped by sll_poisson_3d_periodic_seq test'
        stop
     end if
   
  end do
   
  call delete_poisson_3d_periodic_plan(plan)
   

  print*, 'sll_poisson_3d_periodic_seq test: PASS'
   
  call cpu_time(time_2)
  print*, 'Total CPU time : ', time_2-time_0
     
   
  end program test_poisson_3d
