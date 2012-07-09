      subroutine mgdrelax(sxm,exm,sym,eym,szm,ezm,phi,cof,iters,
     1                    comm3dp,neighbor,bd,phibc,planetype,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,szm,ezm,iters,IOUT
      integer comm3dp,neighbor(26),bd(26),planetype(3)
      REALN phi(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
      REALN cof(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1,8),phibc(6)
c------------------------------------------------------------------------
c Gauss-Seidel point relaxation with Red & Black ordering. Works for
c periodic, Neumann, and Dirichlet boundary conditions.
c 
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : mgdbdry, gxch1pla, MPI_WAITALL (non-blocking version)
c------------------------------------------------------------------------
      integer rb,rbs,it,ipass,i,j,k
      integer ireq,req(52)
# if NBLOCKGR
      integer status(MPI_STATUS_SIZE,52),ierr
# endif
# if cdebug
      double precision tinitial
# if NBLOCKGR
      double precision tmpi
# endif
      tinitial=MPI_WTIME()
# endif
c
c do iters sweeps in the subdomain; impose the wall derivative
c BC after each half-sweep
c
      do it=1,iters
        rb=mod(sxm,2)
        do ipass=1,2
          do k=szm,ezm
            rbs=rb
            do j=sym,eym
              do i=sxm+rb,exm,2
                phi(i,j,k)=(cof(i,j,k,8)-(cof(i,j,k,1)*phi(i-1,j,k)
     1                                   +cof(i,j,k,2)*phi(i+1,j,k)
     2                                   +cof(i,j,k,3)*phi(i,j-1,k)
     3                                   +cof(i,j,k,4)*phi(i,j+1,k)
     4                                   +cof(i,j,k,5)*phi(i,j,k-1)
     5                                   +cof(i,j,k,6)*phi(i,j,k+1)))
     6                     /cof(i,j,k,7)
              end do
              rb=1-rb
            end do
            rb=1-rbs
          end do
          rb=1-mod(sxm,2)
# if WMGD
c
c new version: impose Neumann and Dirichlet boundary conditions
c
          call mgdbdry(sxm,exm,sym,eym,szm,ezm,phi,bd,phibc,IOUT)
# endif
        end do
      end do
c
c Exchange the plane boundary values only once at the end. Since the 
c number of relaxation sweeps at each level is characteristically 
c small (1 or 2 are common values), this does not damage the 
c convergence rate too badly. Overall, I have found a significant 
c reduction in execution time. This also imposes the periodic BCs.
c
# if NBLOCKGR
      ireq=0
# endif
      call gxch1pla(sxm,exm,sym,eym,szm,ezm,phi,comm3dp,neighbor,
     1              bd,planetype,req,ireq,IOUT)
# if NBLOCKGR
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nwaitall=nwaitall+1
      twaitall=twaitall+MPI_WTIME()-tmpi
# endif
# endif
# if WMGD
c
c new version: impose Neumann and Dirichlet boundary conditions
c
      call mgdbdry(sxm,exm,sym,eym,szm,ezm,phi,bd,phibc,IOUT)
# endif
c
# if cdebug
      timing(90)=timing(90)+MPI_WTIME()-tinitial
# endif
      return
      end
