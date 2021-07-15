program test_linear_operator_maxwell_eb_schur
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants
  
  use sll_m_linear_operator_maxwell_eb_schur

  use sll_m_maxwell_3d_fem_fft
  
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type(sll_t_linear_operator_maxwell_eb_schur) :: linear_op_eb_schur
  type(sll_t_maxwell_3d_fem_fft), pointer :: maxwell_solver      !< Maxwell solver

  sll_int32 :: n_dofs(3), degree
  sll_real64 :: domain(3,2)

  sll_real64, allocatable :: rhs(:), lhs(:)
  sll_int32 :: j, r_size
  sll_int32, allocatable :: seed(:)
  
  degree = 3
  n_dofs = [16,8,8]
  domain(1,:) = [0.0_f64, sll_p_twopi]
  domain(2,:) = [0.0_f64, sll_p_twopi]
  domain(3,:) = [0.0_f64, sll_p_twopi]
  
  allocate( maxwell_solver )
  call maxwell_solver%init( domain, n_dofs, [degree,degree,degree] )
  
  call linear_op_eb_schur%create( maxwell_solver )
  linear_op_eb_schur%factor = 1.0_f64

  allocate( rhs( product(n_dofs)*3 ) )
  allocate( lhs( product(n_dofs)*3 ) )

  lhs = 0.0_f64

  call random_seed(size=r_size)
  allocate(seed(r_size))
  seed = 0
  call random_seed(put=seed)
  do j=1,3*product(n_dofs)
     call random_number(lhs(j) )
  end do

!!$  ind=0
!!$  do k=1,n_dofs(3)
!!$     do j=1,n_dofs(2)
!!$        do i=1,n_dofs(1)
!!$           ind = ind+1
!!$           lhs(ind) = cos(sll_p_twopi*(j-1)/real(n_dofs(2),f64))
!!$        end do
!!$     end do
!!$  end do

  write(34,*) lhs
  
  call linear_op_eb_schur%dot( lhs, rhs )

  write(22,*) rhs
  print*, rhs(1)
  
end program test_linear_operator_maxwell_eb_schur
