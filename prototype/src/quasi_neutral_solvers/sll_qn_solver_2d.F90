!> @brief 
!> Selalib 2D (r, theta) quasi-neutral solver 
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
module sll_qn_solver_2d

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_fft
use sll_tridiagonal
use sll_collective
use sll_boundary_condition_descriptors

implicit none

!> QNS solver object 2d
type qn_solver_2d
   sll_int32                   :: BC       !< boundary_conditions
   sll_int32                   :: np_r     !< points in r-direction
   sll_int32                   :: np_theta !< points in theta-direction
   sll_real64                  :: rmin     !< r min
   sll_real64                  :: rmax     !< r max

   type(sll_fft_plan), pointer :: fw_fft   !< forward fft plan
   type(sll_fft_plan), pointer :: bw_fft   !< backward fft plan

   sll_comp64, dimension(:,:), pointer :: hat_rho
   sll_comp64, dimension(:,:), pointer :: hat_phi
   sll_comp64, dimension(:),   pointer :: hat_f
   sll_comp64, dimension(:),   pointer :: hat_g
   sll_int32,  dimension(:),   pointer :: ipiv
   sll_real64, dimension(:),   pointer :: a_resh 
   sll_real64, dimension(:),   pointer :: cts  

end type qn_solver_2d

interface new
module procedure  new_qn_solver_2d
end interface new

interface solve
module procedure  solve_qn_solver_2d
end interface solve

interface delete
module procedure  delete_qn_solver_2d
end interface delete

sll_int32, private :: i, j, k

contains

  !> Initialialize a new QNS solver in 2D space
  function new_qn_solver_2d(bc, rmin, rmax, np_r,  np_theta) result (plan)

    type(qn_solver_2d), pointer :: plan !< qns solver

    sll_int32  :: BC       !< Boundary_conditions
    sll_real64 :: rmin     !< r min
    sll_real64 :: rmax     !< r max
    sll_int32  :: np_r     !< points in r
    sll_int32  :: np_theta !< points in theta
    sll_int32  :: ierr     !< error code

    SLL_ALLOCATE(plan, ierr)
    SLL_CLEAR_ALLOCATE(plan%hat_rho(1:np_r,1:np_theta),ierr) 
    SLL_CLEAR_ALLOCATE(plan%hat_phi(1:np_r,1:np_theta),ierr) 

    SLL_CLEAR_ALLOCATE(plan%hat_f(1:np_theta),ierr)
    SLL_CLEAR_ALLOCATE(plan%hat_g(1:np_theta),ierr)
    SLL_ALLOCATE(plan%ipiv(1:np_r),ierr); plan%ipiv = 0
    SLL_CLEAR_ALLOCATE(plan%a_resh(1:3*np_r),ierr) 
    SLL_CLEAR_ALLOCATE(plan%cts(1:7*np_r),ierr)

    plan%BC       = BC
    plan%np_r     = np_r
    plan%np_theta = np_theta
    plan%rmin     = rmin
    plan%rmax     = rmax

    ! For FFTs in theta-direction
    plan%fw_fft => fft_new_plan( np_theta, plan%hat_f, plan%hat_g, FFT_FORWARD )

    plan%bw_fft => fft_new_plan( np_theta, plan%hat_g, plan%hat_f, FFT_INVERSE )

  end function new_qn_solver_2d


  subroutine solve_qn_solver_2d(plan, rho, c, Te, f, g, Zi, phi)

    type(qn_solver_2d), pointer :: plan
    sll_real64                                     :: dr
    sll_real64                                     :: dtheta
    sll_real64                                     :: Zi 
    sll_int32                                      :: np_r
    sll_int32                                      :: np_theta
    sll_real64, dimension(:,:)                     :: rho, phi
    sll_real64, dimension(:)                       :: c, Te, f, g

    np_r = plan%np_r
    np_theta = plan%np_theta
    if (plan%BC==SLL_NEUMANN) then
       dr = (plan%rmax-plan%rmin)/(np_r-1)
    else ! 'dirichlet'
       dr = (plan%rmax-plan%rmin)/(np_r+1)
    endif
    dtheta = 2*sll_pi/np_theta

    plan%hat_rho = cmplx(rho, 0_f64, kind=f64)
    plan%hat_f   = cmplx(f,   0_f64, kind=f64)
    plan%hat_g   = cmplx(g,   0_f64, kind=f64)

    call fft_apply_plan( plan%fw_fft, plan%hat_f, plan%hat_f )
    call fft_apply_plan( plan%fw_fft, plan%hat_g, plan%hat_g )
    
    call fft_apply_plan( plan%fw_fft, plan%hat_rho(1,:), plan%hat_rho(1,:) )
    call fft_apply_plan( plan%fw_fft, plan%hat_rho(np_r,:), plan%hat_rho(np_r,:) )

    if (plan%bc==SLL_NEUMANN) then
       plan%hat_rho(1,:)  = plan%hat_rho(1,:) + (c(1)-2/dr)*plan%hat_f 
       plan%hat_rho(np_r,:) = plan%hat_rho(np_r,:) + (c(np_r)+2/dr)*plan%hat_g
    else ! dirichlet
       plan%hat_rho(1,:)  = plan%hat_rho(1,:) + (1/dr**2 - c(1)/(2*dr))*plan%hat_f
       plan%hat_rho(np_r,:) = plan%hat_rho(np_r,:) + (1/dr**2 + c(np_r)/(2*dr))*plan%hat_g
    endif

    do i=2,np_r-1
       call fft_apply_plan( plan%fw_fft, plan%hat_rho(i,:), plan%hat_rho(i,:) ) 
    enddo

    do j=1,np_theta

       if (j<=np_theta/2) then
           k = j-1
       else
           k = np_theta-(j-1)
       endif
       
       if (plan%BC==SLL_NEUMANN) then
          call neumann_matrix_resh(np_r, plan%rmin, plan%rmax, np_theta, k, c, Te, Zi, plan%a_resh)
       else 
          call dirichlet_matrix_resh(np_r, plan%rmin, plan%rmax, np_theta, k, c, Te, Zi,  plan%a_resh)
       endif
       call setup_cyclic_tridiag( plan%a_resh, np_r, plan%cts, plan%ipiv )
       call solve_cyclic_tridiag( plan%cts, plan%ipiv, plan%hat_rho(:,j), &
                                  np_r, plan%hat_phi(:,j))         
    enddo

    ! Solution phi of the Quasi-neutral equation is given by taking the inverse
    ! FFT in the k-direction of Tild_phi (storaged in phi)  
    do i=1,np_r
       call fft_apply_plan( plan%bw_fft, plan%hat_phi(i,:), plan%hat_phi(i,:) ) 
    enddo

    phi = real(plan%hat_phi, f64)/np_theta

  end subroutine solve_qn_solver_2d


  subroutine delete_qn_solver_2d(plan)

     type (qn_solver_2d), pointer :: plan
     sll_int32                                       :: ierr

     ! Fixme: some error checking, whether the poisson pointer is associated
     ! for instance
     SLL_ASSERT( associated(plan) )

     call fft_delete_plan(plan%fw_fft)
     call fft_delete_plan(plan%bw_fft)

     SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_qn_solver_2d


  subroutine dirichlet_matrix_resh(np_r, rmin, rmax, np_theta, j, c, Te, Zi, a_resh)

    sll_real64, dimension(:)                       :: c, Te
    sll_real64                                     :: dr, dtheta, Zi
    sll_real64                                     :: r, rmin, rmax
    sll_int32                                      :: i, j, np_r, np_theta
    sll_real64, dimension(:)                       :: a_resh

    dr = (rmax-rmin)/(np_r+1)
    dtheta = 2*sll_pi / np_theta

    a_resh = 0.d0

    do i=1,np_r
       r = rmin + i*dr
       if (i>1) then
          a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       endif
      ! a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2*(1-cos(j*dtheta)) &
      !                                               + 1/(Zi*Te(i))
       a_resh(3*(i-1)+2) = 2/dr**2 + 1/(Zi*Te(i)) + (j/r)**2
                                           
       if (i<np_r) then
          a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
       endif
    enddo

  end subroutine dirichlet_matrix_resh


  subroutine neumann_matrix_resh(np_r, rmin, rmax, np_theta, j, c, Te, Zi, a_resh)

    sll_real64, dimension(:)                       :: c, Te
    sll_real64                                     :: dr, dtheta, Zi
    sll_real64                                     :: rmin, rmax, r
    sll_int32                                      :: i, j, np_r, np_theta
    sll_real64, dimension(:)                       :: a_resh

    dr = (rmax-rmin)/(np_r-1)
    dtheta = 2*sll_pi / np_theta

    a_resh = 0.d0

   ! a_resh(2) = 2/dr**2 + 2/(rmin*dtheta)**2*(1-cos(j*dtheta)) &
   !                                          + 1/(Zi*Te(1))

    a_resh(2) = 2/dr**2 + 1/(Zi*Te(1)) + (k/rmin)**2
    a_resh(3) = -2/dr**2

    do i=2,np_r-1
       r = rmin + (i-1)*dr
       a_resh(3*(i-1)+1) = c(i)/(2*dr) - 1/dr**2
       a_resh(3*(i-1)+2) = 2/dr**2 + 1/(Zi*Te(i)) + (j/r)**2
      ! a_resh(3*(i-1)+2) = 2/dr**2 + 2/(r*dtheta)**2* &
     !                  (1-cos(j*dtheta)) + 1/(Zi*Te(i))
       a_resh(3*(i-1)+3) = -( 1/dr**2 + c(i)/(2*dr) )
    enddo

    a_resh(3*(np_r-1)+1) = -2/dr**2
    !a_resh(3*(np_r-1)+2) = 2/dr**2 + 2/(rmax*dtheta)**2* &
    !                   (1-cos(j*dtheta)) + 1/(Zi*Te(np_r))
    a_resh(3*(np_r-1)+2) = 2/dr**2 + 1/(Zi*Te(np_r)) + (j/rmax)**2

  end subroutine neumann_matrix_resh


end module sll_qn_solver_2d
