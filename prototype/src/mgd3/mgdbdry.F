      subroutine mgdbdry(sxm,exm,sym,eym,szm,ezm,phi,bd,phibc,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,szm,ezm,bd(26),IOUT
      REALN phi(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1),phibc(6)
c------------------------------------------------------------------------
c Enforce the Neumann and Dirichlet boundary conditions
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdrelax, mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j,k
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      if (bd(1).eq.1) then
        do k=szm-1,ezm+1
          do i=sxm-1,exm+1
            phi(i,eym+1,k)=phi(i,eym,k)
          end do
        end do
      else if (bd(1).eq.2) then
        do k=szm-1,ezm+1
          do i=sxm-1,exm+1
            phi(i,eym+1,k)=2.0d0*phibc(1)-phi(i,eym,k)
          end do
        end do
      end if
      if (bd(3).eq.1) then
        do j=sym-1,eym+1
          do i=sxm-1,exm+1
            phi(i,j,szm-1)=phi(i,j,szm)
          end do
        end do
      else if (bd(3).eq.2) then
        do j=sym-1,eym+1
          do i=sxm-1,exm+1
            phi(i,j,szm-1)=2.0d0*phibc(2)-phi(i,j,szm)
          end do
        end do
      end if
      if (bd(5).eq.1) then
        do k=szm-1,ezm+1
          do i=sxm-1,exm+1
            phi(i,sym-1,k)=phi(i,sym,k)
          end do
        end do
      else if (bd(5).eq.2) then
        do k=szm-1,ezm+1
          do i=sxm-1,exm+1
            phi(i,sym-1,k)=2.0d0*phibc(3)-phi(i,sym,k)
          end do
        end do
      end if
      if (bd(9).eq.1) then
        do k=szm-1,ezm+1
          do j=sym-1,eym+1
            phi(sxm-1,j,k)=phi(sxm,j,k)
          end do
        end do
      else if (bd(9).eq.2) then
        do k=szm-1,ezm+1
          do j=sym-1,eym+1
            phi(sxm-1,j,k)=2.0d0*phibc(5)-phi(sxm,j,k)
          end do
        end do
      end if
      if (bd(18).eq.1) then
        do k=szm-1,ezm+1
          do j=sym-1,eym+1
            phi(exm+1,j,k)=phi(exm,j,k)
          end do
        end do
      else if (bd(18).eq.2) then
        do k=szm-1,ezm+1
          do j=sym-1,eym+1
            phi(exm+1,j,k)=2.0d0*phibc(6)-phi(exm,j,k)
          end do
        end do
      end if
c
# if cdebug
      timing(93)=timing(93)+MPI_WTIME()-tinitial
# endif
      return
      end
