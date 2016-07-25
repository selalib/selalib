! Test the module sll_m_intial_distribution
! author: Katharina Kormann, IPP

program test_initial_distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_initial_distribution

  use sll_m_io_utilities, only : &
       sll_s_concatenate_filename_and_path
  
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_distribution_params), allocatable :: params
  sll_int32                                     :: file_id
  character(len=256)                            :: filename
  sll_real64                                    :: xi(2), vi(2), val, val_ref


  call sll_s_concatenate_filename_and_path( "initial_distribution.nml", __FILE__,&
       filename)
  open(newunit=file_id, file=trim(filename))
  call sll_s_initial_distribution_new( "cossum_twogaussian", [2,2], file_id, params )
  close(file_id)

  ! Evaluate params
  xi = [0.1_f64, 2.5_f64]
  vi = [-1.2_f64, 0.5_f64]
  val = params%eval( xi, vi )
  val_ref = 1.0711384771936658E-003_f64
  if ( abs(val-val_ref) > 1E-13_f64 ) then
     print*, 'Error in %eval', val-val_ref
     print*, 'FAILED.'
     stop
  end if
  val = params%evalx( xi )
  val_ref = 1.0026749882862458_f64
  if ( abs(val-val_ref) > 1E-13_f64 ) then
     print*, 'Error in %evalx', val-val_ref
     print*, 'FAILED.'
     stop
  end if
    val = params%evalv( vi )
  val_ref = 1.0682808384643529E-003_f64
  if ( abs(val-val_ref) > 1E-13_f64 ) then
     print*, 'Error in %evalv', val-val_ref
     print*, 'FAILED.'
     stop
  end if

  ! Test passed if we git here
  print*, 'PASSED.'

end program test_initial_distribution
