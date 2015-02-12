! input is the data and its corresponding coordinates_d with size NP
! output is the out at the coordinates_o with size NP
! default that the size of output and input data are the same = NP

! to compile
! gfortran -O3 -ffree-line-length-none sll_memory.F90 sll_assert.F90 sll_working_precision.F90 WENO_interp.F90 unit_test.F90 -o unit_test

!
! the option -O3 is for optimization
! the option -ffree-line-length-none: 
!    -ffree-line-length-n
!    Set column after which characters are ignored in typical free-form lines in the source file. For free-form, the default value is 132. n may be none, meaning that the entire line is meaningful. -ffree-line-length-0 means the same thing as -ffree-line-length-none. 
!    -fmax-identifier-length=n 


program SLWENO_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use WENO_interp
  implicit none
#define NP 60

  sll_int32 :: err    ! indicator for allocating data array
  sll_int32 :: i
  sll_int32 :: i_weno ! indicator for weno(1) or linear interpolation

  type(WENO_interp_1d), pointer :: sp1
  sll_real64, allocatable, dimension(:) :: data       ! data at coordinates_d with size NP
  sll_real64, allocatable, dimension(:) :: coordinates_d  ! coordinates for data with size NP
  sll_real64, allocatable, dimension(:) :: out        ! data at coordinates_i with size NP
  sll_real64, allocatable, dimension(:) :: coordinates_o  ! coordinates for interpolated locations with size NP
  sll_real64, allocatable, dimension(:) :: data_interp  ! exact function value at interpolated locations
  sll_real64 :: dx
  sll_real64 :: accumulator1
  sll_real64 :: sll_pi
  sll_int32  :: order

  
  ! constants
  accumulator1 = 0.0_f64
  sll_pi = 4.0_f64*atan(1.0_f64)
  order = 5
  i_weno = 1
  
  print *, 'SL WENO module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NP), err)
  SLL_ALLOCATE(out(NP), err)
  SLL_ALLOCATE(coordinates_d(NP), err)
  SLL_ALLOCATE(coordinates_o(NP), err)
  SLL_ALLOCATE(data_interp(NP), err)
    
  print *, 'initialize data and coordinates array'
  dx = sll_pi*2.0_f64/NP
  
  do i=1,NP
     coordinates_d(i) = (i-0.5_f64)*dx
     coordinates_o(i) = coordinates_d(i) - dx/2.0_f64
     data(i)        = (sin(coordinates_d(i)) + cos(coordinates_d(i)))
     data_interp(i) = (sin(coordinates_o(i)) + cos(coordinates_o(i)))
  enddo
  print *, 'proceed to allocate the data for WENO interpolation...'
  sp1 =>  new(NP, 0.0_f64, sll_pi*2.0_f64, order, i_weno)  
  
  ! NP is the number of data points, the second and the third argument is the min and max of coordinates
  ! set up the basic information for data: np, xmin, xmax, delta (mesh size), rdelta (reciprocal of delta)
  ! order is the order of interpolation, and i_weno is the index for linear or weno interopolation
  ! when i_weno = 1, it is the weno; otherwise, it is the linear.
  
  call interpolate_WENO_1D( sp1, NP, data, coordinates_o, out)
  
  print *, 'NP=', NP
  print *, 'Contents of the sp1:'
  print *, 'sp1%xmin=', sp1%xmin
  print *, 'sp1%xmax=',sp1%xmax
  print *, 'dx=',sp1%delta
  print *, '1/dx=',sp1%rdelta
  print *, 'left b. =', sp1%xmin - sp1%delta/2.0_f64
  print *, 'right b= ', sp1%xmax + sp1%delta/2.0_f64 
  print *, 'cumulative errors: '
  print *, 'periodic case, NP points: '
  print *, 'interpolating individual values from 1 to NP:'
  do i=1, NP
     accumulator1 = accumulator1 + abs(data_interp(i) - out(i))
     if(i.eq.1)then
        print *, abs(data_interp(i) - out(i)), data_interp(i) , out(i)
     endif
  end do


  print *, '----------------------------------------------------'
  print *, 'RESULTS: '
  print *, 'Periodic case: '
  print *, 'average error at the nodes (single values) = '
  print *, accumulator1/real(NP,f64)

  call delete(sp1)
 
end program SLWENO_tester
