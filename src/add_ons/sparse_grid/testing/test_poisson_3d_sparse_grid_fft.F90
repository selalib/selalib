program test_poisson_3d_sparsegrid_fft
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

   !-------------------------------------------------------------------
   !  test 3D Poisson solver based on FFT on sparse grid
   !-------------------------------------------------------------------

   use sll_m_poisson_3d_sparse_grid_fft
   use sll_m_sparse_grid_3d
   use sll_m_constants, only : sll_pi

   implicit none

   sll_real64                    :: eta_max(3), eta_min(3)
   sll_int32                    :: error

   type(sparse_grid_interpolator_3d) :: interp

   sll_real64, allocatable      :: ex(:)
   sll_real64, allocatable      :: ey(:)
   sll_real64, allocatable      :: ez(:)
   sll_real64, allocatable      :: ex_exact(:)
   sll_real64, allocatable      :: ey_exact(:)
   sll_real64, allocatable      :: ez_exact(:)
   sll_real64, allocatable      :: rho(:), phi_guess(:)
   sll_real64, allocatable      :: phi(:)
   sll_real64, allocatable      :: phi_exact(:)
   sll_real64, allocatable      :: rho_hsp(:), phi_hsp(:)
   type(sll_fft3d_derivative)   :: poisson

   sll_real64                   :: x(3)
   sll_int32                    :: i
   sll_int32                    :: order
   sll_int32                    :: levels(3)
   sll_real64                   :: err(3)
   sll_real64                   :: tol(3)


   order = 2;
   levels(1) = 8; levels(2) = 6; levels(3) = 7;

   eta_min(1) = .0_f64; eta_max(1) = 4.0_f64*sll_pi
   eta_min(2) = .0_f64; eta_max(2) = 4.0_f64*sll_pi
   eta_min(3) = .0_f64; eta_max(3) = 4.0_f64*sll_pi

   call interp%initialize(levels,order, order+1,0,eta_min,eta_max,0,0);

   SLL_ALLOCATE(ex(interp%size_basis),error)
   SLL_ALLOCATE(ey(interp%size_basis),error)
   SLL_ALLOCATE(ez(interp%size_basis),error)
   SLL_ALLOCATE(ex_exact(interp%size_basis),error)
   SLL_ALLOCATE(ey_exact(interp%size_basis),error)
   SLL_ALLOCATE(ez_exact(interp%size_basis),error)
   SLL_ALLOCATE(rho(interp%size_basis),error)
   SLL_ALLOCATE(phi(interp%size_basis),error)
   SLL_ALLOCATE(phi_exact(interp%size_basis),error)

   SLL_ALLOCATE(rho_hsp(interp%size_basis),error)
   SLL_ALLOCATE(phi_hsp(interp%size_basis),error) 
   SLL_ALLOCATE(phi_guess(interp%size_basis),error)


   do i=1,interp%size_basis
      x = interp%hierarchy(i)%coordinate
      phi_exact(i) = sin(x(2))/(2.0_f64+cos(x(1)))
      phi_guess(i) = phi_exact(i)
      ex_exact(i)  =  -sin(x(2))*sin(x(1))/(2.0_f64+cos(x(1)))**2
      ey_exact(i)  =  -cos(x(2))/(2.0_f64+cos(x(1)))
      ez_exact(i)  =  0.0_f64
      rho(i) = - sin(x(2))*(2.0_f64*cos(x(1))+1.0_f64+(sin(x(1)))**2)/(2.0_f64+cos(x(1)))**3+&
           sin(x(2))/(2.0_f64+cos(x(1)));
   end do


   call poisson%initialize( interp ) 
   call poisson%solve(interp,rho,ex,ey,ez)

   err(1) = maxval(abs(ex_exact-ex))
   err(2) = maxval(abs(ey_exact-ey))
   err(3) = maxval(abs(ez_exact-ez))

   write(*,*) " Ex Error = " , err(1)
   write(*,*) " Ey Error = " , err(2)
   write(*,*) " Ez Error = " , err(3)

   tol(1) = 1D-3
   tol(2) = 1D-4
   tol(3) = 1D-14

   if ((err(1) < tol(1)) .AND. (err(2) < tol(2)) .AND. (err(3) < tol(3))) then
      print*, 'PASSED'
   else
      print*, 'FAILED'
   end if

  
 end program test_poisson_3d_sparsegrid_fft
