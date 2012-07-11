      subroutine mgderr(relmax,sxm,exm,sym,eym,szm,ezm,phio,phin,
     1                  comm3d,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,szm,ezm,comm3d,IOUT
      REALN relmax
      REALN phio(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
      REALN phin(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
c------------------------------------------------------------------------
c Calculate the error between the new and old iterates of phi and 
c save the new iterate into the phio array.
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : MPI_ALLREDUCE
c------------------------------------------------------------------------
      REALN locval(2),totval(2)
      integer i,j,k,ierr
# if cdebug
      double precision tinitial,tmpi
      tinitial=MPI_WTIME()
# endif
c
c calculate local error
c
      locval(1)=0.0d0
      locval(2)=0.0d0
      do k=szm,ezm
        do j=sym,eym
          do i=sxm,exm
            locval(1)=max(locval(1),abs(phin(i,j,k)))
            locval(2)=max(locval(2),abs(phin(i,j,k)-phio(i,j,k)))
          end do
        end do
      end do
c
c global reduce across all processes
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_ALLREDUCE(locval,totval,2,MPI_DOUBLE_PRECISION,MPI_MAX,
     1                   comm3d,ierr)
# if cdebug
      nallreduce=nallreduce+1
      tallreduce=tallreduce+MPI_WTIME()-tmpi
# endif
      if (totval(1).gt.0.0) then
        relmax=totval(2)/totval(1)
      else
        relmax=0.0d0
      end if
c
c save new values into ouput array
c
      do k=szm-1,ezm+1
        do j=sym-1,eym+1
          do i=sxm-1,exm+1
            phio(i,j,k)=phin(i,j,k)
          end do
        end do
      end do
c
# if cdebug
      timing(94)=timing(94)+MPI_WTIME()-tinitial
# endif
      return
      end

