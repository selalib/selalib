subroutine mgdrelax(sxm,exm,sym,eym,szm,ezm,phi,cof,iters,  &
     &              comm3dp,neighbor,bd,phibc,planetype,IOUT)

use mpi
implicit none
#include "mgd3.h"
integer :: sxm,exm,sym,eym,szm,ezm,iters,IOUT
integer :: comm3dp,neighbor(26),bd(26),planetype(3)
real(8) :: phi(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
real(8) :: cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8),phibc(6)
!------------------------------------------------------------------------
! Gauss-Seidel point relaxation with Red & Black ordering. Works for
! periodic, Neumann, and Dirichlet boundary conditions.
! 
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : mgdbdry, gxch1pla, MPI_WAITALL (non-blocking version)
!------------------------------------------------------------------------
integer rb,rbs,it,ipass,i,j,k
integer ireq,req(52)

integer status(MPI_STATUS_SIZE,52),ierr

!
! do iters sweeps in the subdomain; impose the wall derivative
! BC after each half-sweep
!
do it=1,iters
  rb=mod(sxm,2)
  do ipass=1,2
    do k=szm,ezm
      rbs=rb
      do j=sym,eym
        do i=sxm+rb,exm,2
          phi(i,j,k)=(cof(i,j,k,8)-(cof(i,j,k,1)*phi(i-1,j,k)   &
     &                                   +cof(i,j,k,2)*phi(i+1,j,k)   &
     &                                   +cof(i,j,k,3)*phi(i,j-1,k)   &
     &                                   +cof(i,j,k,4)*phi(i,j+1,k)   &
     &                                   +cof(i,j,k,5)*phi(i,j,k-1)   &
     &                                   +cof(i,j,k,6)*phi(i,j,k+1))) &
     &                     /cof(i,j,k,7)
        end do
        rb=1-rb
      end do
      rb=1-rbs
    end do
    rb=1-mod(sxm,2)
# if WMGD
!
! new version: impose Neumann and Dirichlet boundary conditions
!
    call mgdbdry(sxm,exm,sym,eym,szm,ezm,phi,bd,phibc,IOUT)
# endif
  end do
end do
!
! Exchange the plane boundary values only once at the end. Since the 
! number of relaxation sweeps at each level is characteristically 
! small (1 or 2 are common values), this does not damage the 
! convergence rate too badly. Overall, I have found a significant 
! reduction in execution time. This also imposes the periodic BCs.
!

ireq=0

call gxch1pla(sxm,exm,sym,eym,szm,ezm,phi,comm3dp,neighbor, &
     &              bd,planetype,req,ireq,IOUT)

call MPI_WAITALL(ireq,req,status,ierr)

# if WMGD
!
! new version: impose Neumann and Dirichlet boundary conditions
!
call mgdbdry(sxm,exm,sym,eym,szm,ezm,phi,bd,phibc,IOUT)
# endif
!
return
end
