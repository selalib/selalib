subroutine mgdrelax(sxm,exm,sym,eym,phi,cof,iters,comm2d,myid, &
                    neighbor,bd,phibc,itype,jtype)

#include "sll_working_precision.h"
use gxch1_2d

sll_int32 :: sxm,exm,sym,eym,iters
sll_int32 :: comm2d,myid,neighbor(8),bd(8),itype,jtype
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),phibc(4)
!------------------------------------------------------------------------
! Gauss-Seidel point relaxation with Red & Black ordering. Works for
! periodic, Neumann, and Dirichlet boundary conditions.
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdkcyc
! Calls     : mgdbdry, gxch1lin
!------------------------------------------------------------------------
sll_int32 :: rb,it,ipass,i,j
!
! do iters sweeps in the subdomain
!
do it=1,iters
  rb=mod(sxm,2)
  do ipass=1,2
    do j=sym,eym
      do i=sxm+rb,exm,2
        phi(i,j)=(cof(i,j,6)-(cof(i,j,1)*phi(i-1,j)               &
                             +cof(i,j,2)*phi(i+1,j)               &
                             +cof(i,j,3)*phi(i,j-1)               &
                             +cof(i,j,4)*phi(i,j+1)))/cof(i,j,5)
      end do
      rb=1-rb
    end do
    rb=1-mod(sxm,2)
# if WMGD
!
! new version: impose Neumann and Dirichlet boundary conditions
!
    call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
# endif
  end do
end do
!
! Exchange boundary data only once at the end. Since the number
! of relaxation sweeps at each level is characteristically small
! (1 or 2 are common values), this does not damage the convergence
! rate too badly. Overall, I have found a significant reduction
! in execution time. This also imposes the periodic BCs.
!
call gxch1lin(phi,comm2d,sxm,exm,sym,eym,neighbor,bd,itype,jtype)
# if WMGD
!
! new version: impose Neumann and Dirichlet boundary conditions
!
call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
# endif 

end subroutine
