! input is the data and its corresponding coordinates_d with size NP
! output is the out at the coordinates_o with size NP
! default that the size of output and input data are the same = NP
! Tests WENO interpolation and reconstruction
! This code is made for cell centered values

program WENO_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use WENO_interp
  use WENO_recon
  use sll_constants
  implicit none
#define NP 30

  sll_int32 :: err    ! indicator for allocating data array
  sll_int32 :: i
  type(WENO_interp_1d) :: WENO_int
  type(WENO_recon_1d)  :: WENO_rec
  ! data at coordinates_d with size NP
  sll_real64, allocatable, dimension(:) :: data 
  ! coordinates for data with size NP
  sll_real64, allocatable, dimension(:) :: coordinates_d
  ! data at coordinates_i with size NP 
  sll_real64, allocatable, dimension(:) :: out        
  ! coordinates for interpolated locations with size NP
  sll_real64, allocatable, dimension(:) :: coordinates_o 
  ! exact function value at interpolated locations
  sll_real64, allocatable, dimension(:) :: data_interp 
  ! exact derivative value at initial locations
  sll_real64, allocatable, dimension(:) :: data_deriv  

  sll_real64 :: dx
  sll_real64 :: accumulator1

  ! constants
  accumulator1 = 0.0_f64

  print *, 'WENO module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NP), err)
  SLL_ALLOCATE(out(NP), err)
  SLL_ALLOCATE(coordinates_d(NP), err)
  SLL_ALLOCATE(coordinates_o(NP), err)
  SLL_ALLOCATE(data_interp(NP), err)
  SLL_ALLOCATE(data_deriv(NP), err)

  print *, 'initialize data and coordinates array'
  dx = sll_pi*2.0_f64/NP

  do i=1,NP
     coordinates_d(i) = (i-0.5_f64)*dx
     coordinates_o(i) = modulo(coordinates_d(i) - dx/3.0_f64, 2*sll_pi)
     data(i)        = 2.0_f64*(sin(coordinates_d(i)) + 2.5_f64 &
          + cos(coordinates_d(i)))
     data_interp(i) = 2.0_f64*(sin(coordinates_o(i)) + 2.5_f64 &
          + cos(coordinates_o(i)))
     data_deriv(i) = 2.0_f64*(cos(coordinates_d(i))  &
          - sin(coordinates_d(i)))
  enddo
  print *, 'Test WENO interpolation...'
  weno_int =  new_WENO_1D(NP, coordinates_d(1), coordinates_d(NP))  
  ! NP is the number of data points, the second and the third argument 
  ! is the min and max of coordinates
  ! set up the basic information for data: np, xmin, xmax, delta (mesh size),
  ! rdelta (reciprocal of delta)
  call interpolate_WENO_1D( weno_int, NP, data, coordinates_o, 4, out )

  print *, 'Contents of the weno:'
  print *, 'weno%xmin=', weno_int%xmin
  print *, 'weno%xmax=',weno_int%xmax
  print *, 'dx=',weno_int%delta
  print *, '1/dx=',weno_int%rdelta
  print *, 'left b. =', weno_int%xmin !- weno%delta/2.0_f64
  print *, 'right b= ', weno_int%xmax !+ weno%delta/2.0_f64 
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
  print *, 'WENO interpolation. Periodic case: '
  print *, 'average error at the nodes (single values) = '
  print *, accumulator1/real(NP,f64)

  call delete_WENO_1D(weno_int)

  print *, 'Test WENO reconstruction...'
  weno_rec =  new_WENO_recon_1D(NP, 0.0_f64, sll_pi*2.0_f64)   
  ! NP is the number of data points, the second and the third argument 
  ! is the min and max of coordinates
  ! set up the basic information for data: np, xmin, xmax, delta (mesh size), 
  ! rdelta (reciprocal of delta)
   call FD_WENO_recon_1D( weno_rec, NP, data,  4, out )

  print *, 'Contents of the weno:'
  print *, 'weno%xmin=', weno_rec%xmin
  print *, 'weno%xmax=',weno_rec%xmax
  print *, 'dx=',weno_rec%delta
  print *, '1/dx=',weno_rec%rdelta
  print *, 'left b. =', weno_rec%xmin !- weno%delta/2.0_f64
  print *, 'right b= ', weno_rec%xmax !+ weno%delta/2.0_f64 
  print *, 'cumulative errors: '
  print *, 'periodic case, NP points: '
  print *, 'reconstructing individual values from 1 to NP:'
  accumulator1 = 0.0_f64
  do i=1, NP
     accumulator1 = accumulator1 + abs(data_deriv(i) - out(i))
     if(i.eq.1)then
        print *, abs(data_deriv(i) - out(i)), data_deriv(i) , out(i)
     endif
  end do


  print *, '----------------------------------------------------'
  print *, 'RESULTS: '
  print *, 'WENO reconstruction. Periodic case: '
  print *, 'average error (single values) = '
  print *, accumulator1/real(NP,f64)

end program WENO_tester
