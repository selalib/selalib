program test_poisson_2d_sparse_grid_fft
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

   !-------------------------------------------------------------------
   !  test 2D Poisson solver based on FFT on sparse grid
   !-------------------------------------------------------------------

   use sll_m_poisson_2d_sparse_grid_fft
   use sll_m_sparse_grid_2d
   use sll_m_constants, only : sll_pi

   implicit none

   sll_real64, dimension(2)  :: eta_max, eta_min
   sll_int32   :: nc_eta1, nc_eta2
   sll_int32   :: error

   type(sparse_grid_interpolator_2d), target :: interp

   sll_real64, dimension(:),allocatable      :: ex
   sll_real64, dimension(:),allocatable      :: ey
   sll_real64, dimension(:),allocatable      :: ex_exact
   sll_real64, dimension(:),allocatable      :: ey_exact
   sll_real64, dimension(:),allocatable      :: rho, phi_guess
   sll_real64, dimension(:),allocatable      :: phi
   sll_real64, dimension(:),allocatable      :: phi_exact
   sll_real64, dimension(:),allocatable        :: rho_hsp, phi_hsp
   type(sll_fft_derivative) :: poisson

   sll_real64                         :: x1, x2
   sll_int32                          :: i
   sll_int32                          :: order
   sll_int32, dimension(2)            :: levels
   sll_real64                         :: tol(2), err(2)

   order = 2;
   levels(1) = 9; levels(2) = 8;

   eta_min(1) = .0_f64; eta_max(1) = 4.0_f64*sll_pi
   eta_min(2) = .0_f64; eta_max(2) = 4.0_f64*sll_pi

   nc_eta1 = 20; nc_eta2 = 20

   call interp%initialize(levels,order, order+1,0,eta_min,eta_max,0,0);

   SLL_ALLOCATE(ex(interp%size_basis),error)
   SLL_ALLOCATE(ey(interp%size_basis),error)
   SLL_ALLOCATE(ex_exact(interp%size_basis),error)
   SLL_ALLOCATE(ey_exact(interp%size_basis),error)
   SLL_ALLOCATE(rho(interp%size_basis),error)
   SLL_ALLOCATE(phi(interp%size_basis),error)
   SLL_ALLOCATE(phi_exact(interp%size_basis),error)

   SLL_ALLOCATE(rho_hsp(interp%size_basis),error)
   SLL_ALLOCATE(phi_hsp(interp%size_basis),error) 
   SLL_ALLOCATE(phi_guess(interp%size_basis),error)


   do i=1,interp%size_basis
      x1 = interp%hierarchy(i)%coordinate(1)
      x2 = interp%hierarchy(i)%coordinate(2)
      phi_exact(i) = sin(x2)/(2.0_f64+cos(x1))
      phi_guess(i) = phi_exact(i)
      ex_exact(i)  =  -sin(x2)*sin(x1)/(2.0_f64+cos(x1))**2
      ey_exact(i)  =  -cos(x2)/(2.0_f64+cos(x1))
      rho(i) = - sin(x2)*(2.0_f64*cos(x1)+1.0_f64+(sin(x1))**2)/(2.0_f64+cos(x1))**3+&
           sin(x2)/(2.0_f64+cos(x1));
   end do

   call poisson%initialize( interp ) 
   call poisson%solve(interp,rho,ex,ey)

   err(1) =  maxval(abs(ex_exact-ex))
   err(2) =  maxval(abs(ey_exact-ey))

   write(*,*) " Ex Error = " , err(1)
   write(*,*) " Ey Error = " , err(2)

   tol(1) = 1.0D-8
   tol(2) = 1.0D-9

   if ((err(1) < tol(1)) .AND. (err(2) < tol(2))) then
      print*, 'PASSED'   
   else
      print*, 'FAILED'
   end if

 end program test_poisson_2d_sparse_grid_fft
