
!***************************************************************************
!
! Selalib 2012     
! Module: sll_qns2d_angular_spect_method_seq.F90
!
!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver with angular spectral method
!> Start date: April 10, 2012
!> Last modification: May 04, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!***************************************************************************

module sll_qn_solver_2d
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
  use sll_constants
  use sll_fft
  use sll_tridiagonal
  use sll_boundary_condition_descriptors

  implicit none

  type qn_solver_2d
     sll_int32                   :: BC ! Boundary_conditions
     sll_int32                   :: NP_r!Number of points in r-direction
     sll_int32                   :: NP_theta!Number of points in theta-direction
     sll_real64                  :: rmin
     sll_real64                  :: rmax
     type(sll_fft_plan), pointer :: fft_plan
     type(sll_fft_plan), pointer :: inv_fft_plan
  end type qn_solver_2d


interface new
module procedure new_qn_solver_2d
end interface new
interface solve
module procedure solve_qn_solver_2d
end interface solve
interface delete
module procedure delete_qn_solver_2d
end interface delete

contains


  function new_qn_solver_2d(BC, rmin, rmax, NP_r, NP_theta) &
                                                                  result (plan)

    sll_int32                                     :: BC ! Boundary_conditions
    sll_real64                                    :: rmin
    sll_real64                                    :: rmax
    sll_comp64, dimension(:),   allocatable       :: x
    sll_int32                                     :: NP_r, NP_theta, ierr
    type(qn_solver_2d), pointer :: plan

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE( x(NP_theta), ierr )

    plan%BC     = BC
    plan%NP_r     = NP_r
    plan%NP_theta = NP_theta
    plan%rmin   = rmin
    plan%rmax   = rmax

    ! For FFTs in theta-direction
    !plan%fft_plan => fft_new_plan_c2c_1d( NP_theta, x, x, FFT_FORWARD )
    plan%fft_plan => fft_new_plan( NP_theta, x, x, FFT_FORWARD )

    ! For inverse FFTs in theta-direction
    !plan%inv_fft_plan => fft_new_plan_c2c_1d( NP_theta, x, x, FFT_INVERSE )
    plan%inv_fft_plan => fft_new_plan( NP_theta, x, x, FFT_INVERSE )

    SLL_DEALLOCATE_ARRAY( x, ierr )

  end function new_qn_solver_2d


  subroutine solve_qn_solver_2d(plan, rho, c, Te, f, g, Zi, phi)

    type(qn_solver_2d),pointer   :: plan
    sll_real64                                     :: dr, dtheta, Zi
    sll_int32                                      :: NP_r, NP_theta, i, j, k
    sll_real64, dimension(:,:)                     :: rho, phi
    sll_comp64, dimension(plan%NP_r,plan%NP_theta) :: hat_rho, hat_phi
    sll_real64, dimension(:)                       :: c, Te, f, g
    sll_comp64, dimension(plan%NP_theta)           :: hat_f, hat_g
    sll_int32, dimension(plan%NP_r)                :: ipiv
    sll_real64, dimension(3*plan%NP_r)             :: a_resh ! 3*n
    sll_real64, dimension(7*plan%NP_r)             :: cts  ! 7*n allocation

    NP_r = plan%NP_r
    NP_theta = plan%NP_theta
    if (plan%BC==SLL_NEUMANN) then
       dr = (plan%rmax-plan%rmin)/(NP_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(NP_r+1)
    endif
    dtheta = 2*sll_pi/NP_theta

    hat_rho = cmplx(rho, 0_f64, kind=f64)
    hat_f = cmplx(f, 0_f64, kind=f64)
    hat_g = cmplx(g, 0_f64, kind=f64)

    !call fft_apply_plan_c2c_1d( plan%fft_plan, hat_f, hat_f )
    !call fft_apply_plan_c2c_1d( plan%fft_plan, hat_g, hat_g )
    !
    !call fft_apply_plan_c2c_1d( plan%fft_plan, hat_rho(1,:), hat_rho(1,:) )
    !call fft_apply_plan_c2c_1d( plan%fft_plan, hat_rho(NP_r,:), hat_rho(NP_r,:) )
    call fft_apply_plan( plan%fft_plan, hat_f, hat_f )
    call fft_apply_plan( plan%fft_plan, hat_g, hat_g )
    
    call fft_apply_plan( plan%fft_plan, hat_rho(1,:), hat_rho(1,:) )
    call fft_apply_plan( plan%fft_plan, hat_rho(NP_r,:), hat_rho(NP_r,:) )

    if (plan%BC==SLL_NEUMANN) then
       hat_rho(1,:)  = hat_rho(1,:) + (c(1)-2/dr)*hat_f 
       hat_rho(NP_r,:) = hat_rho(NP_r,:) + (c(NP_r)+2/dr)*hat_g
    else ! 'dirichlet'
       hat_rho(1,:)  = hat_rho(1,:) + (1/dr**2 - c(1)/(2*dr))*hat_f
       hat_rho(NP_r,:) = hat_rho(NP_r,:) + (1/dr**2 + c(NP_r)/(2*dr))*hat_g
    endif

    do i=2,NP_r-1
       !call fft_apply_plan_c2c_1d( plan%fft_plan, hat_rho(i,:), hat_rho(i,:) ) 
       call fft_apply_plan( plan%fft_plan, hat_rho(i,:), hat_rho(i,:) ) 
    enddo

    do j=1,NP_theta
       if (j<=NP_theta/2) then
          k = j-1
       else
          k = NP_theta-(j-1)
       endif
       if (plan%BC==SLL_NEUMANN) then
          call neumann_matrix(plan%np_r, plan%rmin, plan%rmax, &
                              plan%np_theta, k, c, Te, Zi, a_resh)
       else ! 'dirichlet'
          call dirichlet_matrix(plan%np_r, plan%rmin, plan%rmax, &
                                plan%np_theta, k, c, Te, Zi, a_resh)
       endif
       ! b is given by taking the FFT in the theta-direction of rho_{r,theta}      
       ! Solving the linear system: Ax = b  
       call setup_cyclic_tridiag( a_resh, NP_r, cts, ipiv )
       call solve_cyclic_tridiag( cts, ipiv, hat_rho(:,j), NP_r, hat_phi(:,j))         
    enddo

    ! Solution phi of the Quasi-neutral equation is given by taking the inverse
    ! FFT in the k-direction of Tild_phi (storaged in phi)  
    do i=1,NP_r
       !call fft_apply_plan_c2c_1d( plan%inv_fft_plan, hat_phi(i,:), hat_phi(i,:) ) 
       call fft_apply_plan( plan%inv_fft_plan, hat_phi(i,:), hat_phi(i,:) ) 
    enddo

    phi = real(hat_phi, f64)/NP_theta

  end subroutine solve_qn_solver_2d


  subroutine delete_qn_solver_2d(plan)

    type (qn_solver_2d), pointer :: plan
    sll_int32                                      :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call fft_delete_plan(plan%fft_plan)
    call fft_delete_plan(plan%inv_fft_plan)

    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qn_solver_2d


  subroutine dirichlet_matrix(np_r, rmin, rmax, &
                                             np_theta, k, c, Te, Zi, a_resh)

    sll_real64, dimension(:)                      :: c, Te
    sll_real64                                    :: dr, dtheta, Zi
    sll_real64                                    :: r, rmin, rmax
    sll_int32                                     :: i, k, NP_r, np_theta
    sll_real64, dimension(:)                      :: a_resh

    dr = (rmax-rmin)/(NP_r+1)
    dtheta = 2*sll_pi / NP_theta        

    a_resh = 0.d0

    do i=1,NP_r
       r = rmin + i*dr
       if (i>1) then
          a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       endif
       a_resh(3*(i-1)+2) = 2/dr**2 + 1/(Zi*Te(i)) + (k/r)**2  
       if (i<NP_r) then
          a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
       endif
    enddo

  end subroutine dirichlet_matrix


  subroutine neumann_matrix(np_r, rmin, rmax, np_theta, &
                            k, c, Te, Zi, a_resh)

    sll_real64, dimension(:)                      :: c, Te
    sll_real64                                    :: dr, dtheta, Zi
    sll_real64                                    :: rmin, rmax, r
    sll_int32                                     :: i, k, NP_r, np_theta
    sll_real64, dimension(:)                      :: a_resh

    dr = (rmax-rmin)/(NP_r-1)
    dtheta = 2*sll_pi / NP_theta

    a_resh = 0.d0

    a_resh(2) = 2/dr**2 + 1/(Zi*Te(1)) + (k/rmin)**2
    a_resh(3) = -2/dr**2

    do i=2,NP_r-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       a_resh(3*(i-1)+2) = 2/dr**2 + 1/(Zi*Te(i)) + (k/r)**2
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(NP_r-1)+1) = -2/dr**2
    a_resh(3*(NP_r-1)+2) = 2/dr**2 + 1/(Zi*Te(NP_r)) + (k/rmax)**2

  end subroutine neumann_matrix


end module sll_qn_solver_2d

