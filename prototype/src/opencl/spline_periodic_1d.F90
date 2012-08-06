module spline_periodic_1d

  implicit none

  type :: spline_t
     integer                          :: n
     real, dimension (:), allocatable :: p2
     real, dimension (:), allocatable :: xi
     real, dimension (:), allocatable :: fi
  end type

contains

subroutine test()
  
  integer, parameter :: n=20, m=100
  integer :: i, j, k
  real, dimension (n+1) :: xi, fi
  real, dimension (m) :: xj, fj
  real :: h
  type(spline_t) :: spline
!
! read in data points xi and and data fi
!
  open (unit=7, file="xy.data", status="old", action="read")
  do i = 1, n+1
   read (7,*) xi(i), fi(i)
   write(*,*) xi(i), fi(i)
  end do

  call new_spline( spline, xi, n)

  call compute_spline(spline, fi, n)

  h = (xi(n+1)-xi(1))/m
  xj(1) = xi(1)
  do i = 1, m-1
    xj(i+1) = xj(i) + h
  end do

  call interpolate_array_values(spline, xj, fj )

  do i = 1, m
     write(10,*) xj(i), fj(i)
  end do

end subroutine test

subroutine new_spline( spline, xi, n )
integer, intent(in) :: n
type(spline_t) :: spline
real(4), dimension(n+1) :: xi
spline%n = n
allocate(spline%xi(n+1))
allocate(spline%fi(n+1))
allocate(spline%p2(n+1))

spline%xi = xi

end subroutine new_spline

subroutine compute_spline (spline, fi, n)
!
! function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  implicit none
  type(spline_t) :: spline
  integer :: i
  integer, intent (in) :: n
  real, intent (in), dimension (n+1):: fi
  real, dimension (n):: g, h
  real, dimension (n-1):: d, b, c
!
! assign the intervals and function differences
!
  spline%fi = fi

  do i = 1, n
    h(i) = spline%xi(i+1) - spline%xi(i)
    g(i) = spline%fi(i+1) - spline%fi(i)
  end do
!
! evaluate the coefficient matrix elements
  do i = 1, n-1
    d(i) = 2*(h(i+1)+h(i))
    b(i) = 6*(g(i+1)/h(i+1)-g(i)/h(i))
    c(i) = h(i+1)
  end do
!
! obtain the second-order derivatives
!
  call tridiagonal_linear_eq (n-1, d, c, c, b, g)

  spline%p2(1) = 0
  spline%p2(n+1) = 0
  do i = 2, n 
    spline%p2(i) = g(i-1)
  end do

end subroutine compute_spline
!

subroutine interpolate_array_values( spline, xj, fj )

implicit none
type(spline_t) :: spline
integer  :: m
real(4), dimension(:) :: xj
real(4), dimension(:) :: fj
real :: x, dx, alpha, beta, gamma, eta
integer :: i, k
!
! find the approximation of the function
!
  m = size(xj)

  x = xj(1)
  do i = 1, m-1
    x = xj(i)
!
! find the interval that x resides
    k = 1
    dx = x-spline%xi(1)
    do while (dx .ge. 0)
      k = k + 1
      dx = x-spline%xi(k)
    end do
    k = k - 1
!
! find the value of function f(x)
    dx = spline%xi(k+1) - spline%xi(k)
    alpha = spline%p2(k+1)/(6*dx)
    beta = -spline%p2(k)/(6*dx)
    gamma = spline%fi(k+1)/dx - dx*spline%p2(k+1)/6
    eta = dx*spline%p2(k)/6 - spline%fi(k)/dx


    fj(i) =  alpha*(x-spline%xi(k))*(x-spline%xi(k))*(x-spline%xi(k)) &
       +beta*(x-spline%xi(k+1))*(x-spline%xi(k+1))*(x-spline%xi(k+1)) &
       +gamma*(x-spline%xi(k))+eta*(x-spline%xi(k+1))

  end do

  fj(m) = fj(1)

end subroutine interpolate_array_values

!
subroutine tridiagonal_linear_eq (l, d, e, c, b, z)
!
! functione to solve the tridiagonal linear equation set.
!
  implicit none
  integer, intent (in) :: l
  integer :: i
  real, intent (in), dimension (l):: d, e, c, b
  real, intent (out), dimension (l):: z
  real, dimension (l):: y, w
  real, dimension (l-1):: v, t
!
! evaluate the elements in the lu decomposition
!
  w(1) = d(1)
  v(1) = c(1)
  t(1) = e(1)/w(1)
  do i = 2, l - 1
    w(i) = d(i)-v(i-1)*t(i-1)
    v(i) = c(i)
    t(i) = e(i)/w(i)
  end do
  w(l) = d(l)-v(l-1)*t(l-1)
!
! forward substitution to obtain y
!
  y(1) = b(1)/w(1)
  do i = 2, l
    y(i) = (b(i)-v(i-1)*y(i-1))/w(i)
  end do
!
! backward substitution to obtain z
  z(l) = y(l)
  do i = l-1, 1, -1
    z(i) = y(i) - t(i)*z(i+1)
  end do
end subroutine tridiagonal_linear_eq

end module spline_periodic_1d
