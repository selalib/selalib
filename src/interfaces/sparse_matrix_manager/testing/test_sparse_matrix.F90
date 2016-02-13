program unit_test_sparse_matrix

#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_sparse_matrix
#ifdef PASTIX
use sll_m_collective
#endif /* PASTIX */

implicit none

type(sll_t_csr_matrix),     pointer     :: mat
sll_int32                               :: num_rows
sll_int32                               :: num_elements
sll_int32,  dimension(:,:), allocatable :: local_to_global_row
sll_real64, dimension(:),   allocatable :: b
sll_real64, dimension(:),   allocatable :: x
sll_int32                               :: num_local_dof_row
sll_int32                               :: ierr
sll_int32                               :: i
sll_real64                              :: val

#ifdef PASTIX

sll_int32                               :: prank
sll_int32                               :: psize

call sll_s_boot_collective()
prank = sll_f_get_collective_rank( sll_v_world_collective )
psize = sll_f_get_collective_size( sll_v_world_collective )
write(*,*) " #", prank, " of ", psize

#endif /* PASTIX */

num_rows = 2000
num_elements = num_rows
num_local_dof_row = 1

SLL_ALLOCATE(local_to_global_row(num_local_dof_row,num_elements),ierr)
SLL_ALLOCATE(b(num_rows),ierr)
SLL_ALLOCATE(x(num_rows),ierr)

local_to_global_row(1:num_local_dof_row,1:num_elements) = 0
do i=1,num_elements
  local_to_global_row(1,i)=i
enddo

b = 1._f64

print *,'#begin sll_f_new_csr_matrix'

mat => sll_f_new_csr_matrix(  num_rows,            &
                              num_rows,            &
                              num_elements,        &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              local_to_global_row, &
                              num_local_dof_row)

print *,'#end sll_f_new_csr_matrix'

!do identity matrix
print *,'#begin add_to_csr_matrix'
val = 1._f64
do i=1,num_rows
  call sll_s_add_to_csr_matrix(mat, val, i, i)
enddo
print *,'#end add_to_csr_matrix'

x = 0.0_f64
print *,'#begin sll_s_solve_csr_matrix'
call sll_s_solve_csr_matrix(mat, b, x)
print *,'#end sll_s_solve_csr_matrix'
  
if(maxval(abs(x-1.0_f64))<1.e-13)then
  print *,'#PASSED'
else
  print *,'#system is not good resolved',maxval(abs(x))
endif

call sll_s_delete_csr_matrix(mat)

#ifdef PASTIX
call sll_s_halt_collective()
#endif /* PASTIX */

print *,'#begin sll_f_new_csr_matrix with umfpack'

#ifdef UMFPACK
mat  => sll_f_new_csr_matrix( num_rows,            &
                              num_rows,            &
                              num_elements,        &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              local_to_global_row, &
                              num_local_dof_row,   &
                              sll_p_umfpack )

print *,'#end sll_f_new_csr_matrix with umfpack'

!do identity matrix
print *,'#begin add_to_csr_matrix with umfpack'
val = 1._f64
do i=1,num_rows
  call sll_s_add_to_csr_matrix(mat, val, i, i)
enddo
print *,'#end add_to_csr_matrix with umfpack'

call sll_s_factorize_csr_matrix(mat)
b = 1._f64
x = 0.0_f64
print *,'#begin sll_s_solve_csr_matrix with umfpack'
call sll_s_solve_csr_matrix(mat, b, x)
print *,'#end sll_s_solve_csr_matrix with umfpack'
  
if(maxval(abs(x-1.0_f64))<1.e-13)then
  print *,'#PASSED'
else
  print *,'#system is not good resolved with umfpack',maxval(abs(x))
endif

call sll_s_delete_csr_matrix(mat)
#endif /* UMFPACK */

#ifdef PASTIX
call sll_s_halt_collective()
#endif /* PASTIX */

contains

!> @brief
!> Test function to initialize a CSR matrix
!> @details
!> Fill a matrix in CSR format corresponding to a constant coefficient
!> five-point stencil on a square grid. This function comes from AGMG
!> test program.

subroutine uni_laplace_2d(m,f,a,ja,ia)
sll_int32  :: m
sll_real64 :: f(:)
sll_real64 :: a(:)
sll_int32  :: ia(:)
sll_int32  :: ja(:)

sll_int32  :: k,l,i,j
sll_real64, parameter :: zero =  0.0_f64
sll_real64, parameter :: cx   = -1.0_f64
sll_real64, parameter :: cy   = -1.0_f64
sll_real64, parameter :: cd   =  4.0_f64

k=0
l=0
ia(1)=1
do i=1,m
  do j=1,m
    k=k+1
    l=l+1
    a(l)=cd
    ja(l)=k
    f(k)=zero
    if(j < m) then
      l=l+1
      a(l)=cx
      ja(l)=k+1
    else
      f(k)=f(k)-cx
    end if
    if(i < m) then
      l=l+1
      a(l)=cy
      ja(l)=k+m
    else
      f(k)=f(k)-cy
    end if
    if(j > 1) then
      l=l+1
      a(l)=cx
      ja(l)=k-1
    else
      f(k)=f(k)-cx
    end if
    if(i >  1) then
      l=l+1
      a(l)=cy
      ja(l)=k-m
    else
      f(k)=f(k)-cy
    end if
    ia(k+1)=l+1
  end do
end do

end subroutine uni_laplace_2D

end program
