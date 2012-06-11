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

program test_poisson_3d_seq

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
#include "sll_remap.h"

  use numeric_constants
  use sll_poisson_3d_periodic_seq
  use sll_diagnostics

  implicit none

  sll_int32  :: nx, ny, nz
  sll_real64 :: Lx, Ly, Lz
  sll_int32  :: error
  sll_real64 :: dx, dy, dz
  sll_real64, dimension(:,:,:), allocatable    :: x, y, z
  sll_real64                                   :: xx, yy, zz
  sll_real64, dimension(:,:,:), allocatable    :: rho, phi_an, phi
  sll_int32                                    :: i, j, k
  type (poisson_3d_periodic_plan_seq), pointer :: plan
  sll_real64                                   :: average_err
  sll_real64                                   :: time_0, time_1, time_2
  sll_int32                                    :: i_test

  call cpu_time(time_0)
  nx = 128
  ny = 64
  nz = 32

  Lx = 1.0_f64
  Ly = 1.0_f64
  Lz = 1.0_f64

  dx = Lx/(nx-1)
  dy = Ly/(ny-1)
  dz = Lz/(nz-1)

  print*, 'Initialize Poisson 3D solver plan'
  SLL_ALLOCATE(x(nx,ny,nz),error)
  SLL_ALLOCATE(y(nx,ny,nz),error)
  SLL_ALLOCATE(z(nx,ny,nz),error)

  call write_mesh(x,y,z,nx,ny,nz,"mesh3d")

  SLL_ALLOCATE(rho(nx,ny,nz),error)
  SLL_ALLOCATE(phi(nx,ny,nz),error)
  SLL_ALLOCATE(phi_an(nx,ny,nz),error)

  plan => new_poisson_3d_periodic_plan_seq(nx, ny, nz, Lx, Ly, Lz)

  print*, ' '
  call cpu_time(time_1)
  print*,' CPU time = ', time_1-time_0

  do k=1,nz
     zz = (k-1)*dz
     do j=1,ny
        yy = (j-1)*dy
        do i=1,nx
           xx = (i-1)*dx
           rho(i,j,k) = 10000 * exp(-( (xx-Lx/2)**2 + (yy-Ly/2)**2 + &
                                                       (zz-Lz/2)**2 ))
        enddo
     enddo
   enddo

   call cpu_time(time_1)
   call solve_poisson_3d_periodic_seq(plan, rho, phi)
   call cpu_time(time_2)

   print*, 'CPU time = ', time_2-time_1

   call write_vec1d(rho,nx,ny,nz,"rho"//char(i_test+48),"mesh3d",0)
   call write_vec1d(phi,nx,ny,nz,"phi"//char(i_test+48),"mesh3d",0)

  call delete_poisson_3d_periodic_plan_seq(plan)

  call cpu_time(time_2)
  print*, 'Total CPU time : ', time_2-time_0


end program test_poisson_3d_seq
