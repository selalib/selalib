! input is the data with size NP
! output is the out: the derivatives of data with size NP
! default that the size of output and input data are the same = NP

! to compile
! gfortran -O3 -ffree-line-length-none sll_memory.F90 sll_assert.F90 sll_working_precision.F90 WENO_recon.F90 unit_test.F90 -o unit_test

!
! the option -O3 is for optimization
! the option -ffree-line-length-none: 
!    -ffree-line-length-n
!    Set column after which characters are ignored in typical free-form lines in the source file. For free-form, the default value is 132. n may be none, meaning that the entire line is meaningful. -ffree-line-length-0 means the same thing as -ffree-line-length-none. 
!    -fmax-identifier-length=n 


program WENO_recon_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use WENO_recon
  implicit none
#define NP 40

  sll_int32 :: err    ! indicator for allocating data array
  sll_int32 :: i
    sll_int32 :: i_weno ! indicator for weno or linear reconstruction
  sll_int32 :: order
  type(WENO_recon_1D), pointer :: sp1
  sll_real64, allocatable, dimension(:) :: data       ! data at cell centers with size NP
  sll_real64, allocatable, dimension(:) :: out        ! data at coordinates_i with size NP
  sll_real64, allocatable, dimension(:) :: data_deri  ! exact function value at interpolated locations
  sll_real64 :: dx
  sll_real64 :: accumulator1
  sll_real64 :: sll_pi

  
  ! constants
  accumulator1 = 0.0_f64
  sll_pi = 4.0_f64*atan(1.0_f64)
  order = 5
  i_weno = 1
  
  print *, 'WENO recon module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NP), err)
  SLL_ALLOCATE(out(NP), err)
  SLL_ALLOCATE(data_deri(NP), err)
    
  print *, 'initialize data and coordinates array'
  dx = sll_pi*2.0_f64/NP
 
  do i=1,NP      
     data(i)  = sin((i-0.5_f64)*dx)+cos((i-0.5_f64)*dx) 
     data_deri(i) =cos((i-0.5_f64)*dx)-sin((i-0.5_f64)*dx)
  enddo
  print *, 'proceed to allocate the data for WENO reconstruction...'
  sp1 =>  new_WENO_recon_1D(NP, 0.0_f64, sll_pi*2.0_f64, order, i_weno)  
  ! NP is the number of data points, the second and the third argument is the min and max of coordinates
  ! set up the basic information for data: np, xmin, xmax, delta (mesh size), rdelta (reciprocal of delta)
  call FD_WENO_recon_1D(sp1, NP, data, out)

  print *, 'NP=', NP
  print *, 'Contents of the sp1:'
  print *, 'sp1%xmin=', sp1%xmin
  print *, 'sp1%xmax=',sp1%xmax
  print *, 'dx=',sp1%delta
  print *, '1/dx=',sp1%rdelta
  print *, 'order=', order
  print *, 'cumulative errors: '
  print *, 'periodic case, NP points: '
  print *, 'derivatives of function values'
  do i=1, NP
     accumulator1 = accumulator1 + abs(data_deri(i) - out(i))
     if(i.eq.1)then
        print *, abs(data_deri(i) - out(i)), data_deri(i) , out(i)
     endif
  end do


  print *, '----------------------------------------------------'
  print *, 'RESULTS: '
  print *, 'Periodic case: '
  print *, 'average error at the nodes (single values) = '
  print *, accumulator1/real(NP,f64)

  call delete_WENO_recon_1D(sp1)
 
end program WENO_recon_tester
