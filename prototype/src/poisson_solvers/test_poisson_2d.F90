program test_poisson_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_poisson_solvers.h"
#include "sll_constants.h"

   use sll_poisson_2D_periodic

   implicit none

   sll_real64  :: eta1_max, eta1_min, eta2_max, eta2_min
   sll_int32   :: nc_eta1, nc_eta2
   sll_int32   :: error


   sll_real64, dimension(:,:),allocatable      :: ex
   sll_real64, dimension(:,:),allocatable      :: ey
   sll_real64, dimension(:,:),allocatable      :: ex_exact
   sll_real64, dimension(:,:),allocatable      :: ey_exact
   sll_real64, dimension(:,:),allocatable      :: rho
   sll_real64, dimension(:,:),allocatable      :: rhs
   sll_real64, dimension(:,:),allocatable      :: phi
   sll_real64, dimension(:,:),allocatable      :: phi_exact
   type(poisson_2d_periodic)                   :: poisson

   sll_real64                         :: x1, x2
   sll_int32                          :: mode
   sll_int32                          :: i, j

   eta1_min = .0_f64; eta1_max = 2.0_f64*sll_pi
   eta2_min = .0_f64; eta2_max = 2.0_f64*sll_pi

   nc_eta1 = 127; nc_eta2 = 127

   SLL_CLEAR_ALLOCATE(ex(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(ey(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(ex_exact(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(ey_exact(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(rhs(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(rho(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(phi(nc_eta1+1,nc_eta2+1),error)
   SLL_CLEAR_ALLOCATE(phi_exact(nc_eta1+1,nc_eta2+1),error)

   write(*,*) " eta1_min, eta1_max, nc_eta1 ", eta1_min, eta1_max, nc_eta1
   write(*,*) " eta2_min, eta2_max, nc_eta2 ", eta2_min, eta2_max, nc_eta2

   call initialize( poisson, eta1_min, eta1_max, nc_eta1, &
                    eta2_min, eta2_max, nc_eta2, error) 

   open(14, file="test_poisson_2d_rho.dat")
   mode = 2
   do i = 1, nc_eta1+1
      do j = 1, nc_eta2+1
         x1 = (i-1)*(eta1_max-eta1_min)/nc_eta1
         x2 = (j-1)*(eta2_max-eta2_min)/nc_eta2
         phi_exact(i,j) = mode * sin(mode*x1) * cos(mode*x2)
         ex_exact(i,j)  =  1_f64*mode**2*cos(mode*x1)*cos(mode*x2)
         ey_exact(i,j)  = -1_f64*mode**2*sin(mode*x1)*sin(mode*x2)
         rho(i,j) = -2_f64 * mode**3 * sin(mode*x1)*cos(mode*x2)
         write(14,*) x1, x2, rho(i,j)
      end do
   end do

   rhs = rho
   call solve( poisson, phi, rhs)
   write(*,*) " Po Error = " , maxval(abs(phi_exact+phi))
   rhs = rho
   call solve( poisson, phi, rhs)
   write(*,*) " Po Error = " , maxval(abs(phi_exact+phi))
   rhs = rho
   call solve( poisson, ex, ey, rhs)
   write(*,*) " Ex Error = " , maxval(abs(ex_exact-ex))
   write(*,*) " Ey Error = " , maxval(abs(ey_exact-ey))
   rhs = rho
   call solve( poisson, ex, ey, rhs)
   write(*,*) " Ex Error = " , maxval(abs(ex_exact-ex))
   write(*,*) " Ey Error = " , maxval(abs(ey_exact-ey))

end program test_poisson_2d
