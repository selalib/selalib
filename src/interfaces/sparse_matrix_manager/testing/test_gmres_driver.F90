!*************************************************************************
!**                 DRIVER EXAMPLE FOR THE GMRes CODE
!*************************************************************************
program validation

integer, parameter :: lda = 500, ldstrt = 40
integer, parameter :: lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1
integer, parameter :: matvec=1, precondLeft=2, precondRight=3, dotProd=4
real(8), parameter :: ZERO = 0.0_8, ONE = 1.0_8

integer :: i, j, n, m
integer :: revcom, colx, coly, colz, nbscal
integer :: irc(5), icntl(8), info(3)
integer :: nout

real(8) ::  a(lda,lda), work(lwork)
real(8) ::  cntl(5), rinfo(2)

!***************************************************************
!** Generate the test matrix a and set the right-hand side
!** in positions (n+1) to 2n of the array work.
!** The right-hand side is chosen such that the exact solution
!** is the vector of all ones.
!***************************************************************

write(*,*) '***********************************************'
write(*,*) 'This code is an example of use of GMRES'
write(*,*) 'Results are written in output files'
write(*,*)  'fort.20','  ','sol_dTest'
write(*,*) '***********************************************'
write(*,*)
write(*,*) 'Matrix size < ', lda
n = 100
if (n.gt.lda) then
  write(*,*) 'You are asking for a too large matrix'
  goto 100
endif

a = ZERO
do i = 1,n
  a(i,i) = 4.0_8
enddo
do i = 1,n-1
  a(i,i+1) = -2.0_8
  a(i+1,i) = -1.0_8
enddo

!*********************************
!** Choose the restart parameter
!*********************************

write(*,*) 'Restart  <', ldstrt
m = 1

!*******************************************************
!** Initialize the control parameters to default value
!*******************************************************

call init_dgmres(icntl,cntl)

!*************************
!*c Tune some parameters
!*************************

! Save the convergence history on standard output
icntl(3) = 6
! Maximum number of iterations
icntl(7) = 1000 
! preconditioner location
icntl(4) = 0
! orthogonalization scheme
icntl(5)=0
!initial guess
icntl(6) = 0
!residual calculation strategy at restart
icntl(8) = 1
!Initialise the right hand side
print*, 'n=', n, one, size(work)
work(1:n) = ONE

call dgemv('N',n,n,ONE,A,lda,work(1),1,ZERO,work(n+1),1)
do j = 1,n
  work(j) = ONE/2.0_8
enddo
!*****************************************
!** Reverse communication implementation
!*****************************************

10 continue

call drive_dgmres(n,n,m,lwork,work, &
&         irc,icntl,cntl,info,rinfo)

revcom = irc(1)
colx   = irc(2)
coly   = irc(3)
colz   = irc(4)
nbscal = irc(5)

if (revcom.eq.matvec) then
  ! perform the matrix vector product work(colz) <-- A * work(colx)
  call dgemv('N',n,n,ONE,a,lda,work(colx),1,ZERO,work(colz),1)
  goto 10

else if (revcom.eq.precondLeft) then
  ! perform the left preconditioning work(colz) <-- M^{-1} * work(colx)
  call dcopy(n,work(colx),1,work(colz),1)
  call dtrsm('L','L','N','N',n,1,ONE,A,lda,work(colz),n)
  goto 10

else if (revcom.eq.precondRight) then
  !perform the right preconditioning
  call dcopy(n,work(colx),1,work(colz),1)
  call dtrsm('L','U','N','N',n,1,ONE,A,lda,work(colz),n)
  goto 10

else if (revcom.eq.dotProd) then
  ! perform the scalar product work(colz) <-- work(colx) work(coly)

  call dgemv('C',n,nbscal,ONE,work(colx),n,work(coly),1,ZERO,work(colz),1)
  goto 10
endif

!*******************************
!* dump the solution on a file
!*******************************

nout = 6
open(nout,FILE='sol_dTest',STATUS='unknown')
if (icntl(5).eq.0) then
  write(nout,*) 'Orthogonalisation : MGS'
elseif (icntl(5).eq.1) then
  write(nout,*) 'Orthogonalisation : IMGS'
elseif (icntl(5).eq.2) then
  write(nout,*) 'Orthogonalisation : CGS'
elseif (icntl(5).eq.3) then
  write(nout,*) 'Orthogonalisation : ICGS'
endif
write(nout,*) 'Restart : ', m
write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
write(nout,*) 'Optimal workspace = ', info(3)
write(nout,*) ' ************************************************ '
write(nout,*)
close(nout)

100 continue

stop
end
