!> Enforce the Neumann and Dirichlet boundary conditions
!>
!> Code      : mgd2, 2-D parallel multigrid solver
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!> Called in : mgdrelax, mgdsolver
!> Calls     : --
subroutine mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
#include "sll_working_precision.h"

sll_int32  :: sxm,exm,sym,eym,bd(8)
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1),phibc(4)
sll_int32 :: i,j

if (bd(1).eq.1) then
   do j=sym-1,eym+1
      phi(exm+1,j)=phi(exm,j)
   end do
else if (bd(1).eq.2) then
   do j=sym-1,eym+1
      phi(exm+1,j)=2.0d0*phibc(1)-phi(exm,j)
   end do
end if
if (bd(5).eq.1) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=phi(sxm,j)
   end do
else if (bd(5).eq.2) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=2.0d0*phibc(3)-phi(sxm,j)
   end do
end if
if (bd(3).eq.1) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=phi(i,sym)
   end do
else if (bd(3).eq.2) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=2.0d0*phibc(2)-phi(i,sym)
   end do
end if
if (bd(7).eq.1) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=phi(i,eym)
   end do
else if (bd(7).eq.2) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=2.0d0*phibc(4)-phi(i,eym)
   end do
end if

end subroutine
