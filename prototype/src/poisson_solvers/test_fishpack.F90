program test_fishpack

use fishpack

implicit none
!-----------------------------------------------
!   l o c a l   v a r i a b l e s
!-----------------------------------------------
integer :: idimf, m, n, mp1, np1
integer :: i, j, ierror

real(8), allocatable, dimension(:,:) :: f
real(8) :: a,b,c,d
real(8) :: err, z, w
real(8) :: pi, piby2, pisq

pi    = 4.0*atan(1.0)
pisq  = pi*pi
piby2 = 0.5*pi

m = 40; n = 80
mp1 = m + 1; np1 = n + 1
a =  0.; b =  2.; c = -1.; d =  3.
allocate(f(mp1,np1)); f = 0.0
call test_cartesian_2d(f,a,b,m,c,d,n)
deallocate(f)


m = 50; n = 48
mp1 = m + 1; np1 = n + 1
a = 0.; b = 1.; c = 0.; d = 0.5*pi
allocate(f(mp1,np1)); f = 0.0
call test_polar_2d(f,a,b,m,c,d,n)
deallocate(f)

contains

subroutine test_cartesian_2d(field, eta1_min, eta1_max, nc_eta1, &
                            eta2_min, eta2_max, nc_eta2)
implicit none
integer, intent(in)                     :: nc_eta1, nc_eta2
real(8), dimension(nc_eta1+1,nc_eta2+1) :: field
real(8), dimension(nc_eta1+1)           :: eta1
real(8), dimension(nc_eta2+1)           :: eta2
real(8), intent(in)                     :: eta1_min, eta1_max
real(8), intent(in)                     :: eta2_min, eta2_max

type(poisson_fishpack) :: poisson

call new_poisson_2d_fishpack(poisson, CARTESIAN_2D, &
                             eta1_min, eta1_max, nc_eta1, &
                             eta2_min, eta2_max, nc_eta2  )

poisson%mbdcnd = 2
poisson%nbdcnd = 0
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

call solve_poisson_fishpack(poisson, field)

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

subroutine test_polar_2d(field, eta1_min, eta1_max, nc_eta1, &
                        eta2_min, eta2_max, nc_eta2)
implicit none
integer, intent(in)                     :: nc_eta1, nc_eta2
real(8), dimension(nc_eta1+1,nc_eta2+1) :: field
real(8), dimension(nc_eta1+1)           :: bdc, bdd
real(8)                                 :: bda, bdb 
real(8), dimension(nc_eta1+1)           :: eta1
real(8), dimension(nc_eta2+1)           :: eta2
real(8), intent(in)                     :: eta1_min, eta1_max
real(8), intent(in)                     :: eta2_min, eta2_max
integer                                 :: mbdcnd, nbdcnd
real(8)                                 :: elmbda, pertrb

type(poisson_fishpack) :: poisson


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

call new_poisson_2d_fishpack(poisson, POLAR_2D, &
                             eta1_min, eta1_max, nc_eta1, &
                             eta2_min, eta2_max, nc_eta2  )


poisson%mbdcnd = 5
poisson%nbdcnd = 3
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

call solve_poisson_fishpack (poisson, field)

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
 
end program test_fishpack
