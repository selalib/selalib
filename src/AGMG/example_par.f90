! COPYRIGHT (c) 2012 Universite Libre de Bruxelles (ULB)
!
! This file is part of AGMG software package
! Release 3.2.1-aca built on "Mar 20 2014" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE,
! PLEASE CONTACT ynotay@ulb.ac.be
!
! In particular, if you have a free academic license:
!
! (1) You must be a member of an educational, academic or research institution.
!     The license agreement automatically terminates once you no longer fulfill
!     this requirement.
!
! (2) You are obliged to cite AGMG in any publication or report as:
!     "Yvan Notay, AGMG software and documentation;
!      see http://homepages.ulb.ac.be/~ynotay/AGMG".
!
! (3) You may not make available to others the software in any form, either
!     as source or as a precompiled object.
!
! (4) You may not use AGMG for the benefit of any third party or for any
!     commercial purposes. Note that this excludes the use within the
!     framework of a contract with an industrial partner.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! See the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG> for
! release information and possible upgrade.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DICLAIMER:
!    AGMG is provided on an "AS IS" basis, without any explicit or implied
!    WARRANTY; see the file "LICENSE" for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If you use AGMG for research, please observe that your work benefits
!   our past research efforts that allowed the development of AGMG.
!   Hence, even if you do not see it directly, the results obtained thanks
!   to the use of AGMG depend on the results in publications [1-3] below,
!   where the main algorithms used in AGMG are presented and justified.
!   It is then a normal duty to cite these publications (besides citing
!   AGMG itself) in any scientific work depending on the usage of AGMG,
!   as you would do with any former research result you are using.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        program example_par
!
!  Solves the discrete Laplacian on the unit square by simple call to agmg.
!  The right-hand-side is such that the exact solution is the vector of all 1.
!  Uses a strip partitioning of the domain, with internal boundaries parallel
!       to the x direction.
!
        implicit none
        include 'mpif.h'
        real (kind(0d0)),allocatable :: a(:),f(:),x(:)
        integer,allocatable :: ja(:),ia(:),listrank(:)
        integer :: n,iter,iprint,nhinv,NPROC,IRANK,mx,my,ifirstlistrank,ierr
        real (kind(0d0)) :: tol
        character*10 filename
!
!       set inverse of the mesh size (feel free to change)
        nhinv=1000
!
!       maximal number of iterations
        iter=50
!
!       tolerance on relative residual norm
        tol=1.e-6
!
!  Initialize MPI
!
       call MPI_INIT(ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,ierr)
!
!       unit number for output messages (alternative: iprint=10+IRANK)
        iprint=10
        filename(1:8)='res.out_'
        write (filename(9:10),'(i2.2)') IRANK        ! processor dependent
        open(iprint,file=filename,form='formatted')
!
!    calculate local grid size
!
       mx=nhinv-1
       my=(nhinv-1)/NPROC
       if (IRANK < mod(nhinv-1,NPROC)) my=my+1
!
!       generate the matrix in required format (CSR)
!
!         first allocate the vectors with correct size
            N=mx*my
            allocate (a(5*N),ja(5*N),ia(N+1),f(N),x(n),listrank(2*mx))
!           external nodes connected with local ones on top and bottom
!           internal boundaries will receive numbers [N+1,...,N+2*mx]
            ifirstlistrank=N+1
!         next call subroutine to set entries
!           before, initialize listrank to zero so that entries
!           that do not correspond to a nonlocal variable present
!           in ja are anyway properly defined
         listrank(1:2*mx)=0
         call uni2dstrip(mx,my,f,a,ja,ia,IRANK,NPROC,listrank,ifirstlistrank)
!
!       call agmg
!         argument 5 (ijob)  is 0 because we want a complete solve
!         argument 7 (nrest) is 1 because we want to use flexible CG
!                            (the matrix is symmetric positive definite)
!
         call dagmgpar(N,a,ja,ia,f,x,0,iprint,1,iter,tol,         &
                       MPI_COMM_WORLD,listrank,ifirstlistrank)
!
!      uncomment the following to write solution on disk for checking
!
!        filename(1:8)='sol.out_'
!        write (filename(9:10),'(i2.2)') IRANK         ! processor dependent
!        open(11,file=filename,form='formatted')
!        write(11,'(e22.15)') x(1:n)
!        close(11)
      END program example_par
!----------------------------------------------------------------------
    subroutine uni2dstrip(mx,my,f,a,ja,ia,IRANK,NPROC,listrank,ifirstlistrank)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a rectangular grid
! Bottom boundary is an internal boundary if IRANK > 0, and
!    top boundary is an internal boundary if IRANK < NPROC-1
!
      implicit none
      real (kind(0d0)) :: f(*),a(*)
      integer :: mx,my,ia(*),ja(*),ifirstlistrank,listrank(ifirstlistrank:*)
      integer :: IRANK,NPROC,k,l,i,j
      real (kind(0d0)), parameter :: zero=0.0d0,cx=-1.0d0,cy=-1.0d0, cd=4.0d0
!
      k=0
      l=0
      ia(1)=1
      do i=1,my
        do j=1,mx
          k=k+1
          l=l+1
          a(l)=cd
          ja(l)=k
          f(k)=zero
          if(j < mx) then
             l=l+1
             a(l)=cx
             ja(l)=k+1
            else
             f(k)=f(k)-cx
          end if
          if(i < my) then
             l=l+1
             a(l)=cy
             ja(l)=k+mx
            else if (IRANK == NPROC-1) then
             f(k)=f(k)-cy             !real boundary
            else
             l=l+1                  !internal boundary (top)
             a(l)=cy                ! these external nodes are given the
             ja(l)=k+mx             !  numbers [mx*my+1,...,mx*(my+1)]
             listrank(k+mx)=IRANK+1 !Thus listrank(mx*my+1:mx*(my+1))=IRANK+1
          end if
          if(j > 1) then
             l=l+1
             a(l)=cx
             ja(l)=k-1
            else
             f(k)=f(k)-cx
          end if
          if(i >  1) then
             l=l+1
             a(l)=cy
             ja(l)=k-mx
            else if (IRANK == 0) then
             f(k)=f(k)-cy             !real boundary
            else
             l=l+1                    !internal boundary (bottom)
             a(l)=cy                  ! these external nodes are given the
             ja(l)=k+mx*(my+1)        ! numbers [mx*(my+1)+1,...,mx*(my+2)]
             listrank(k+mx*(my+1))=IRANK-1
                              !Thus listrank(mx*(my+1)+1:mx*(my+2))=IRANK-1
          end if
          ia(k+1)=l+1
        end do
      end do
!
      return
    end subroutine uni2dstrip
