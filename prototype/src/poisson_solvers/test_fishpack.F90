program test_fishpack

use sll_fishpack

implicit none
integer :: l, m, n, mp1, np1
integer :: i, j, k

real(8), allocatable, dimension(:,:)   :: f_2d
real(8), allocatable, dimension(:,:)   :: f_polar
real(8), allocatable, dimension(:,:,:) :: f_3d
real(8) :: pi, piby2, pisq, err, z
real(8) :: xs, xf, ys, yf, zs, zf

pi    = 4.0*atan(1.0)
pisq  = pi*pi
piby2 = 0.5*pi

m = 40; n = 80
mp1 = m + 1; np1 = n + 1
xs =  0.; xf =  2.; ys = -1.; yf =  3.
allocate(f_2d(mp1,np1)); f_2d = 0.0
call test_cartesian_2d(f_2d,xs,xf,m,2,ys,yf,n,0)


m = 50; n = 48
mp1 = m + 1; np1 = n + 1
xs = 0.; xf = 1.; ys = 0.; yf = 0.5*pi
allocate(f_polar(mp1,np1)); f_2d = 0.0
call test_polar_2d(f_polar,xs,xf,m,5,ys,yf,n,3)

xs = 0.; xf = 1.
ys = 0.; yf = 2.*pi
zs = 0.; zf = pi/2.
l = 10
m = 40
n = 15

allocate(f_3d(l+1,m+1,n+1)); f_3d = 0.0
call test_cartesian_3d(f_3d,xs,xf,l,1,ys,yf,m,0,zs,zf,n,2)

deallocate(f_2d)
deallocate(f_polar)
deallocate(f_3d)

contains

subroutine test_cartesian_2d(field,                               &
                             eta1_min, eta1_max, nc_eta1, bc_eta1, &
                             eta2_min, eta2_max, nc_eta2, bc_eta2  )
implicit none
integer, intent(in)                     :: nc_eta1, nc_eta2
integer, intent(in)                     :: bc_eta1, bc_eta2 
real(8), dimension(nc_eta1+1,nc_eta2+1) :: field
real(8), dimension(nc_eta1+1)           :: eta1
real(8), dimension(nc_eta2+1)           :: eta2
real(8), intent(in)                     :: eta1_min, eta1_max
real(8), intent(in)                     :: eta2_min, eta2_max

type(fishpack_2d) :: poisson

call poisson%initialize(CARTESIAN_2D, &
         eta1_min, eta1_max, nc_eta1, bc_eta1,&
         eta2_min, eta2_max, nc_eta2, bc_eta2  )

poisson%bc_eta1 = bc_eta1
poisson%bc_eta2 = bc_eta2
poisson%elmbda = -4.

!     generate and store grid points for the purpose of computing
!     boundary data and the right side of the helmholtz equation.

do i = 1, nc_eta1+1
   eta1(i) = eta1_min+(eta1_max-eta1_min)*float(i-1)/nc_eta1
end do

do j = 1, nc_eta2+1
   eta2(j) = eta2_min+(eta2_max-eta2_min)*float(j-1)/nc_eta2
end do

!     generate boundary data.
!     bda, bdc, and bdd are dummy variables.

poisson%bdb(:) = 4.*cos((eta2(:)+1.)*piby2)

!     generate right side of equation.

do i = 2, nc_eta1+1
   do j = 1, nc_eta2+1
      field(i,j) = (2. - (4. + pisq/4.)*eta1(i)**2)*cos((eta2(j)+1.)*piby2)
   end do
end do

call poisson%solve(field)

!     compute discretization error.  the exact solution is
!                u(x,y) = x**2*cos((y+1)*piby2)

err = 0.
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      z = abs(field(i,j)-eta1(i)**2*cos((eta2(j)+1.)*piby2))
      err = max(z,err)
   end do
end do

!     print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer

write (*, *) '    fishpack 2d cartesian test run *** '
write (*, *) '    previous 64 bit floating point arithmetic result '
write (*, *) '    ierror = 0,  discretization error = 5.36508-4'
write (*, *) '    previous 32 bit floating point arithmetic result '
write (*, *) '    ierror = 0,  discretization error = 4.9305e-4'
write (*, *) '    the output from your computer is: '
write (*, *) '    ierror =', poisson%error, ' discretization error = ', err

end subroutine test_cartesian_2d

subroutine test_polar_2d(field, &
                         eta1_min, eta1_max, nc_eta1, bc_eta1,&
                         eta2_min, eta2_max, nc_eta2, bc_eta2)
implicit none
integer, intent(in)                     :: nc_eta1, nc_eta2
integer, intent(in)                     :: bc_eta1, bc_eta2
real(8), dimension(nc_eta1+1,nc_eta2+1) :: field
real(8), dimension(nc_eta1+1)           :: eta1
real(8), dimension(nc_eta2+1)           :: eta2
real(8), intent(in)                     :: eta1_min, eta1_max
real(8), intent(in)                     :: eta2_min, eta2_max

type(fishpack_2d) :: poisson


!--------------------------------------------------------------------------
!									  !
!     program to illustrate the use of subroutine hwsplr to solve	  !
!     the equation							  !
!									  !
!     (1/r)(d/dr)(r*(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta) = 16*r**2	  !
!									  !
!     on the quarter-disk 0 .lt. r .lt. 1, 0 .lt. theta .lt. pi/2 with	  !
!     with the boundary conditions					  !
!									  !
!     u(1,theta) = 1 - cos(4*theta), 0 .le. theta .le. 1		  !
!									  !
!     and								  !
!									  !
!     (du/dtheta)(r,0) = (du/dtheta)(r,pi/2) = 0,  0 .le. r .le. 1.	  !
!									  !
!     (note that the solution u is unspecified at r = 0.)		  !
!          the r-interval will be divided into 50 panels and the	  !
!     theta-interval will be divided into 48 panels.			  !
!									  !
!     from dimension statement we get value of idimf.			  !
!									  !
!-------------------------------------------------------------------------!

call poisson%initialize(POLAR_2D, &
         eta1_min, eta1_max, nc_eta1, bc_eta1,&
         eta2_min, eta2_max, nc_eta2, bc_eta2   )

poisson%elmbda = 0.

do i = 1, nc_eta1+1
   eta1(i) = eta1_min+(eta1_max-eta1_min)*float(i-1)/nc_eta1
end do

do j = 1, nc_eta2+1
   eta2(j) = eta2_min+(eta2_max-eta2_min)*float(j-1)/nc_eta2
end do

poisson%bdc(:) = 0.
poisson%bdd(:) = 0.

do j = 1, nc_eta2+1
   field(nc_eta1+1,j) = 1. - cos(4.*eta2(J))
end do

do i = 1, nc_eta1
   field(i,:nc_eta2+1) = 16.*eta1(i)**2
end do

call poisson%solve(field)

err = 0.
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      z = abs(field(i,j)-eta1(i)**4*(1.-cos(4.*eta2(j))))
      err = max(z,err)
   end do
end do

write (*, *) '    HWSPLR TEST RUN *** '
write (*, *) '    Previous 64 bit floating point arithmetic result '
write (*, *) '    IERROR = 0,  Discretization Error = 6.19134E-4'
write (*, *) '    Previous 32 bit floating point arithmetic result '
write (*, *) '    IERROR = 0,  Discretization Error = 6.20723E-4'
write (*, *) '    The output from your computer is: '
write (*, *) '    IERROR =', poisson%error, ' Discretization Error = ', ERR

end subroutine test_polar_2d



subroutine test_cartesian_3d(field, &
                             eta1_min, eta1_max, nc_eta1, bc_eta1,&
                             eta2_min, eta2_max, nc_eta2, bc_eta2,&
                             eta3_min, eta3_max, nc_eta3, bc_eta3)

implicit none

integer, intent(in)                     :: nc_eta1, nc_eta2, nc_eta3
integer, intent(in)                     :: bc_eta1, bc_eta2, bc_eta3
real(8), intent(in)                     :: eta1_min, eta1_max
real(8), intent(in)                     :: eta2_min, eta2_max
real(8), intent(in)                     :: eta3_min, eta3_max

real(8), dimension(nc_eta1+1,nc_eta2+1,nc_eta3+1) :: field

type(fishpack_3d) :: poisson

real(8) , allocatable, dimension(:) :: eta1
real(8) , allocatable, dimension(:) :: eta2
real(8) , allocatable, dimension(:) :: eta3
real(8) :: t, delta_eta1, delta_eta2, delta_eta3

call poisson%initialize(CARTESIAN_3D,     &
         eta1_min, eta1_max, nc_eta1, bc_eta1,&
         eta2_min, eta2_max, nc_eta2, bc_eta2,&
         eta3_min, eta3_max, nc_eta3, bc_eta3)
 
poisson%elmbda = -3.

delta_eta1 = (eta1_max - eta1_min)/float(nc_eta1)
allocate(eta1(nc_eta1+1))
do i = 1, nc_eta1+1
   eta1(i) = eta1_min + float(i - 1)*delta_eta1
end do
delta_eta2 = (eta2_max - eta2_min)/float(nc_eta2)
allocate(eta2(nc_eta2+1))
do j = 1, nc_eta2+1
   eta2(j) = eta2_min + float(j - 1)*delta_eta2
end do
delta_eta3 = (eta3_max - eta3_min)/float(nc_eta3)
allocate(eta3(nc_eta3+1))
do k = 1, nc_eta3+1
   eta3(k) = eta3_min + float(k - 1)*delta_eta3
end do

do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      poisson%bdf(i,j) = -eta1(i)**4*sin(eta2(j))
   end do
end do

do j = 1, nc_eta2+1
   do k = 1, nc_eta3+1
      field(1,j,k) = 0.
      field(nc_eta1+1,j,k) = sin(eta2(j))*cos(eta3(k))
   end do
end do
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      field(i,j,1) = eta1(i)**4*sin(eta2(j))
   end do
end do

do i = 2, nc_eta1
   do j = 1, nc_eta2+1
      do k = 2, nc_eta3+1
         field(i,j,k) = 4.*eta1(i)**2*(3. - eta1(i)**2)*sin(eta2(j))*cos(eta3(k))
      end do
   end do
end do

call poisson%solve(field)


!     compute discretization error.  the exact solution to the
!     problem is
!
!        u(x,y,z) = x**4*sin(y)*cos(z)
!
err = 0.
do i = 1, nc_eta1+1
   do j = 1, nc_eta2+1
      do k = 1, nc_eta3+1
         t = abs(field(i,j,k)-eta1(i)**4*sin(eta2(j))*cos(eta3(k)))
         err = max(t,err)
      end do
   end do
end do
!     print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
write (*, *) '    hw3crt test run *** '
write (*, *) '    previous 64 bit floating point arithmetic result '
write (*, *) '    ierror = 0,  discretization error = 9.6480e-3'
write (*, *) '    previous 32 bit floating point arithmetic result '
write (*, *) '    ierror = 0,  discretization error = 9.6480e-3'
write (*, *) '    the output from your computer is: '
write (*, *) '    ierror =', poisson%error, ' discretization error = ', err

end subroutine test_cartesian_3d
 
end program test_fishpack
