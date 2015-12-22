program test_mudpack_cartesian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_2d

  use sll_m_mudpack, only: &
    sll_s_delete_mudpack_cartesian, &
    sll_s_initialize_mudpack_cartesian, &
    sll_t_mudpack_solver, &
    sll_s_solve_mudpack_cartesian

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

integer :: nc_eta1
integer :: nc_eta2
type(sll_t_mudpack_solver) :: periodic
type(sll_t_mudpack_solver) :: dirichlet
real(8), allocatable :: sol(:,:)
real(8), allocatable :: phi(:,:)
real(8), allocatable :: rhs(:,:)
real(8), allocatable :: eta1(:,:)
real(8), allocatable :: eta2(:,:)

real(8) :: eta1_min, eta1_max, eta2_min, eta2_max
real(8) :: delta_eta1, delta_eta2

integer :: i, j, error

nc_eta1 = 64
nc_eta2 = 64

SLL_CLEAR_ALLOCATE(eta1(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(eta2(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(sol( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(phi( 1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(rhs( 1:nc_eta1+1,1:nc_eta2+1),error)

!set end points of solution rectangle in (x,y) space
eta1_min = 0.0_8
eta1_max = 4.0_8
eta2_min = 0.0_8
eta2_max = 4.0_8

delta_eta1 = (eta1_max-eta1_min)/dble(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/dble(nc_eta2)
do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      eta1(i,j) = eta1_min+dble(i-1)*delta_eta1
      eta2(i,j) = eta2_min+dble(j-1)*delta_eta2
   end do
end do

!Poisson periodic

call sll_s_initialize_mudpack_cartesian(periodic,  &
                eta1_min, eta1_max, nc_eta1, &
                eta2_min, eta2_max, nc_eta2, &
                sll_p_periodic, sll_p_periodic,  &
                sll_p_periodic, sll_p_periodic)


sol  = sin(2*sll_p_pi*eta1)*sin(2*sll_p_pi*eta2)
rhs  = -8*sll_p_pi**2 * sol + 1.

call sll_s_solve_mudpack_cartesian(periodic, phi, rhs)

call sll_o_gnuplot_2d(eta1_min, eta1_max, nc_eta1+1, &
                    eta2_min, eta2_max, nc_eta2+1, &
                    phi, "sinsin", 1, error)

!compute and print maximum norm of error
write(*,201) maxval(abs(phi-sol))

call sll_s_delete_mudpack_cartesian(periodic)

print*,"PASSED"

eta1_min = -5.0_8
eta1_max =  5.0_8
eta2_min = -5.0_8
eta2_max =  5.0_8

delta_eta1 = (eta1_max-eta1_min)/dble(nc_eta1)
delta_eta2 = (eta2_max-eta2_min)/dble(nc_eta2)
do i=1,nc_eta1+1
   do j=1,nc_eta2+1
      eta1(i,j) = eta1_min+dble(i-1)*delta_eta1
      eta2(i,j) = eta2_min+dble(j-1)*delta_eta2
   end do
end do
     
call sll_s_initialize_mudpack_cartesian(dirichlet,  &
                eta1_min, eta1_max, nc_eta1,  &
                eta2_min, eta2_max, nc_eta2,  &
                sll_p_dirichlet, sll_p_dirichlet, &
                sll_p_dirichlet, sll_p_dirichlet)


sol = exp(-(eta1*eta1+eta2*eta2))

call sll_o_gnuplot_2d(eta1_min, eta1_max, nc_eta1+1, &
                    eta2_min, eta2_max, nc_eta2+1, &
                    sol, "sol_dirichlet", 1, error)

do j=2,nc_eta2
   do i=2,nc_eta1
      rhs(i,j) = (sol(i-1,j)-2.*sol(i,j)+sol(i+1,j))/(delta_eta1*delta_eta1) &
               + (sol(i,j-1)-2.*sol(i,j)+sol(i,j+1))/(delta_eta2*delta_eta2)
   end do
end do

!rhs = 4 * sol * (eta1*eta1 + eta2*eta2 - 1)

call sll_o_gnuplot_2d(eta1_min, eta1_max, nc_eta1+1, &
                    eta2_min, eta2_max, nc_eta2+1, &
                    rhs, "rhs_dirichlet", 1, error)

!rhs = 4.0_8
phi(:,1) = sol(:,1)
phi(:,nc_eta2+1) = sol(:,nc_eta2+1)
phi(1,:) = sol(1,:)
phi(nc_eta2+1,:) = sol(nc_eta2+1,:)

call sll_s_solve_mudpack_cartesian(dirichlet, phi, rhs)

call sll_o_gnuplot_2d(eta1_min, eta1_max, nc_eta1+1, &
                    eta2_min, eta2_max, nc_eta2+1, &
                    phi, "phi_dirichlet", 1, error)

!compute and print maximum norm of error
write(*,201) maxval(abs(phi-sol))

call sll_s_delete_mudpack_cartesian(dirichlet)

print*,"PASSED"

201 format(' maximum error  =  ',e10.3)
     
end program test_mudpack_cartesian
