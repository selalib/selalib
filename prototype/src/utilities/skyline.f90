! fortran structure for solving linear systems
! the matrix can be symmetric or not.
! It is stored in the skyline format
! into three vectors "vdiag", "vlow", "vsup" for representing
! respectively its diagonal, lower part and upper part.
! The shape of the non-zero part of the matrix is symmetric.
! The profile "prof" of the matrix is the number of nonzero terms 
! above the diagonal for each column (or the number of nonzero terms
! to the left of the diagonal for each line).
! A fast access is provided to the matrix with the array "mkld" pointing
! in "vsup" to the beginning of each column 
! (or in "vlow" to the beginning of each line).
! Subroutines are provided for accessing (i,j) terms of the matrix,
! computing matrix vector product,
! computing LU decomposition and solving linear system
module skyline
  implicit none

  type skyline_matrix
     integer :: neq=0 ! size of the linear system
     logical :: sym=.false. ! true if symmetric system    
     integer :: memsize ! size in memory of vlow and vsup
     real*8,dimension(:), allocatable :: vdiag,vlow,vsup
     integer,dimension(:),allocatable :: prof,mkld
  end type skyline_matrix


contains

  ! create the matrix
  subroutine create(sky,n)
    implicit none
    integer :: n
    type(skyline_matrix) :: sky

    sky%neq=n
    allocate(sky%prof(n))
    sky%prof=0
    allocate(sky%mkld(n+1))
    sky%mkld=0
    allocate(sky%vdiag(n))
    sky%vdiag=0

  end subroutine create

  ! allocate the matrix from the profile
  subroutine init(sky)
    implicit none
    integer :: i
    type(skyline_matrix) :: sky

    ! array of column start
    sky%mkld(1)=1
    do i=1,sky%neq
       sky%mkld(i+1)=sky%mkld(i)+sky%prof(i)
    enddo

    sky%memsize=sky%mkld(sky%neq+1)-1

    allocate(sky%vsup(sky%memsize))
    sky%vsup=0
    allocate(sky%vlow(sky%memsize))
    sky%vlow=0
  end subroutine init

  ! add a value to  a term of the matrix
  subroutine add(sky,val,i,j)

    implicit none
    integer :: i,j,ll
    real*8 :: val
    type(skyline_matrix) :: sky

    if (i-j.gt.sky%prof(i).or.j-i.gt.sky%prof(j)) then
       write(*,*) '(',i,',',j,') out of matrix profile'
       write(*,*) sky%prof(i)
       write(*,*) sky%prof(j)
       stop
    end if
    if (i.eq.j) then
       sky%vdiag(i)=sky%vdiag(i)+val
    else if (j.gt.i) then
       ll=sky%mkld(j+1)-j+i
       sky%vsup(ll)=sky%vsup(ll)+val
    else 
       ll=sky%mkld(i+1)-i+j
       sky%vlow(ll)=sky%vlow(ll)+val
    end if

  end subroutine add

  ! get the (i,j) term of the matrix
  subroutine get(sky,val,i,j)

    implicit none
    integer :: i,j,ll
    real*8 :: val
    type(skyline_matrix) :: sky

    if (abs(i-j).gt.sky%prof(i).or.abs(i-j).gt.sky%prof(j)) then
       val=0
    else if (i.eq.j) then
       val=sky%vdiag(i)
    else if (j.gt.i) then
       ll=sky%mkld(j+1)-j+i
       val=sky%vsup(ll)
    else 
       ll=sky%mkld(i+1)-i+j
       val=sky%vlow(ll)
    end if

  end subroutine get

  ! scalar produc
  function scal(x,y,n)
    implicit none
    integer :: n,i
    real*8 :: x(n),y(n),scal

    scal=0
    do i=1,n
       scal=scal+x(i)*y(i)
    end do
  end function scal

  ! compute (in place) the LU decomposition
  subroutine computeLU(sky)
    implicit none
    type(skyline_matrix) :: sky

    integer :: nsym=1,mp=6,ifac=1,isol=0,ier
    real*8 :: energ,vu,vfg

    call sol(sky%vsup,sky%vdiag,sky%vlow,   &
         vfg,sky%mkld,vu,sky%neq,mp,ifac,isol, &
         nsym,energ,ier,sky%memsize)  
  end subroutine computeLU

  ! solve a linear system
  subroutine solve(sky,rhs,solution)
    implicit none
    type(skyline_matrix) :: sky

    integer :: nsym=1,mp=6,ifac=0,isol=1,ier
    real*8 :: energ,solution(sky%neq),rhs(sky%neq)

    call sol(sky%vsup,sky%vdiag,sky%vlow,   &
         rhs,sky%mkld,solution,sky%neq,mp,ifac,isol, &
         nsym,energ,ier,sky%memsize)  
  end subroutine solve

  ! compute matrix vector product (real)
  subroutine mul(sky,x,Ax)
    implicit none
    type(skyline_matrix) :: sky

    integer :: nsym=1
    real*8 :: x(sky%neq),Ax(sky%neq)

    call mulku(sky%vsup,sky%vdiag,sky%vlow,   &
         sky%mkld,x,sky%neq, &
         nsym,Ax,sky%memsize)  

  end subroutine mul

  ! compute matrix vector product (complex)
  subroutine cmul()
  end subroutine cmul

end module skyline


program test_sky

  use skyline
  implicit none

  type(skyline_matrix)  :: sky
  integer,parameter :: n=10
  integer :: i,j
  real*8,dimension(n,n) :: A=0
  real*8,dimension(n) :: b1,x1,x2,b2


  A(1,1) = 5
  A(1,2) = -2
  A(1,3) = -1
  A(2,1) = -1
  A(2,2) = 5
  A(2,3) = -2
  A(2,4) = -1
  A(3,1) = -1
  A(3,2) = -1
  A(3,3) = 5
  A(3,4) = -2
  A(3,5) = -1
  A(3,8) = 1
  A(4,2) = -1
  A(4,3) = -1
  A(4,4) = 5
  A(4,5) = -2
  A(4,6) = -1
  A(5,3) = -1
  A(5,4) = -1
  A(5,5) = 5
  A(5,6) = -2
  A(5,7) = -1
  A(6,4) = -1
  A(6,5) = -1
  A(6,6) = 5
  A(6,7) = -2
  A(6,8) = -1
  A(7,5) = -1
  A(7,6) = -1
  A(7,7) = 5
  A(7,8) = -2
  A(7,9) = -1
  A(8,6) = -1
  A(8,7) = -1
  A(8,8) = 5
  A(8,9) = -2
  A(8,10) = -1
  A(9,7) = -1
  A(9,8) = -1
  A(9,9) = 5
  A(9,10) = -2
  A(10,8) = -1
  A(10,9) = -1
  A(10,10) = 5

  x1(1) = 4
  x1(2) = 2
  x1(3) = -5
  x1(4) = -9
  x1(5) = 1
  x1(6) = -7
  x1(7) = -5
  x1(8) = -6
  x1(9) = 6
  x1(10) = -3


  write(*,*) 'A=',A
  write(*,*) 'x1=',x1
  b1=matmul(A,x1)
  write(*,*) 'b1=',b1

  call create(sky,n)

  sky%prof(1) = 0
  sky%prof(2) = 1
  sky%prof(3) = 2
  sky%prof(4) = 2
  sky%prof(5) = 2
  sky%prof(6) = 2
  sky%prof(7) = 2
  sky%prof(8) = 5
  sky%prof(9) = 2
  sky%prof(10) = 2

!!$  sky%prof(1)=0
!!$  sky%prof(2)=1
!!$  sky%prof(3)=2

  call init(sky)

  do i=1,n
     do j=1,n
        if (A(i,j).ne.0) then
           call add(sky,A(i,j),i,j)
        end if
     end do
  end do

  call mul(sky,x1,b2)
  write(*,*) 'b2=',b2

  call computeLU(sky)

  call solve(sky,b1,x2)
  write(*,*) 'x2=',x2


end program test_sky

