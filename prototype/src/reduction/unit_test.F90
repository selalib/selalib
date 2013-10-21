program test_reduction
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_reduction_module
implicit none

  type(sll_logical_mesh_1d), pointer    :: m_x4
  sll_real64, dimension(:,:,:,:), allocatable   :: data_4d
  sll_real64, dimension(:,:,:), allocatable   :: data_3d
  type(sll_logical_mesh_4d), pointer :: m4d
  sll_int32 :: ierr
  sll_real64 :: err
  sll_int32 :: N(4) = (/32,32,32,64/)
  
  m_x4 => new_logical_mesh_1d(N(4))
  m4d => new_logical_mesh_4d(N(1),N(2),N(3),N(4))
  SLL_ALLOCATE(data_4d(N(1)+1,N(2)+1,N(3)+1,N(4)+1),ierr)
  SLL_ALLOCATE(data_3d(N(1)+1,N(2)+1,N(3)+1),ierr)
  
  
  
  data_4d = 1._f64
  
  call compute_reduction_4d_to_3d(&
    m_x4, &
    data_4d, &
    data_3d, &
    m4d=m4d)

  err=maxval(abs(data_3d-1._f64))
  
  
  print *,'#err=',err

  if(err==0._f64)then
    print *,'#PASSED'
  endif

end program test_reduction