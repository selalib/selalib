! new definition of dpx,dpy
! cubic spline basis
module m_particules
#include "sll_working_precision.h"
#include "sll_memory.h"
use m_zone
use m_quietstart

implicit none

type particle
  sll_real64, pointer :: dpx(:)
  sll_real64, pointer :: dpy(:)
  sll_int32 , pointer :: idx(:)
  sll_int32 , pointer :: idy(:)
  sll_real64, pointer :: vpx(:)
  sll_real64, pointer :: vpy(:)
  sll_real64, pointer :: epx(:)
  sll_real64, pointer :: epy(:)
  sll_real64, pointer :: p(:)
end type particle

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plasma( ele )

type (particle) :: ele
sll_real64 :: speed, theta, vth, n
sll_real64 :: xi, yi, zi
sll_real64 :: temm, z(3)
sll_int32  :: k, error

vth    = 1.0_f64
n      = 1.0_f64/real(nbpart,f64) ! int(f(t=0,x,v)) = dimx*dimy

SLL_ALLOCATE(ele%dpx(nbpart),error)
SLL_ALLOCATE(ele%dpy(nbpart),error)
SLL_ALLOCATE(ele%idx(nbpart),error)
SLL_ALLOCATE(ele%idy(nbpart),error)
SLL_ALLOCATE(ele%vpx(nbpart),error)
SLL_ALLOCATE(ele%vpy(nbpart),error)
SLL_ALLOCATE(ele%epx(nbpart),error)
SLL_ALLOCATE(ele%epy(nbpart),error)
SLL_ALLOCATE(ele%p(nbpart),error)

do k=0,nbpart-1

  speed = vth * sqrt(-2.0_f64 * log( (k+0.5)*n ))

  theta = trinary_reversing( k ) * 2.0_f64 *  sll_p_pi

  ele%vpx(k+1) = speed * cos(theta)
  ele%vpy(k+1) = speed * sin(theta)

  ele%p(k+1)   = dimx * dimy * n

enddo

!--add by Xiaofei : rejection sampling--
k = 1
do while (k<=nbpart)
  call random_number(z)
  xi   = dimx*z(1)
  yi   = dimy*z(2)
  zi   = (2.0_f64+alpha)*z(3)
  temm = 1.0_f64+sin(yi)+alpha*cos(kx*xi)
  if (temm >= zi) then
    ele%idx(k) = floor(xi/dimx*nx)
    ele%idy(k) = floor(yi/dimy*ny)
    ele%dpx(k) = real(xi/dx - ele%idx(k), f64)
    ele%dpy(k) = real(yi/dy - ele%idy(k), f64)
    k=k+1
  endif
enddo

end subroutine plasma

subroutine interpol_eb_cic( f, ele )

type(particle)    :: ele
type(mesh_fields) :: f
sll_real32 :: a1, a2, a3, a4
sll_int32  :: k
sll_int32  :: i, j
sll_real64 :: dpx
sll_real64 :: dpy

! i,j+1_________i+1,j+1
!  |     |        |
!  | a2     a1    |
!  |____ P _______|
!  |              |
!  | a3  |  a4    |
!  |     |        |
! i,j____|______i+1,j

do k=1,nbpart

   i   = ele%idx(k)
   j   = ele%idy(k)
   dpx = ele%dpx(k)
   dpy = ele%dpy(k)

   a1 = (1.0-dpx) * (1.0-dpy)
   a2 = (    dpx) * (1.0-dpy)
   a3 = (    dpx) * (    dpy)
   a4 = (1.0-dpx) * (    dpy)

   ele%epx(k) = a1 * f%ex(i  ,j  ) + a2 * f%ex(i+1,j  ) &
            & + a3 * f%ex(i+1,j+1) + a4 * f%ex(i  ,j+1)
   ele%epy(k) = a1 * f%ey(i  ,j  ) + a2 * f%ey(i+1,j  ) &
            & + a3 * f%ey(i+1,j+1) + a4 * f%ey(i  ,j+1)

end do

end subroutine interpol_eb_cic



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_rho_cic( ele, f )

type(particle)    :: ele
type(mesh_fields) :: f
sll_real64        :: a1, a2, a3, a4, weight
sll_real64        :: rho_total
sll_int32         :: k
sll_int32         :: i, j
sll_real64        :: dpx
sll_real64        :: dpy

f%r0 = 0.d0

!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do k=1,nbpart

  i             = ele%idx(k)
  j             = ele%idy(k)
  dpx           = ele%dpx(k)
  dpy           = ele%dpy(k)
  weight        = ele%p(k)
  a1            = (1.0_f64-dpx) * (1.0_f64-dpy) * weight
  a2            = (dpx)         * (1.0_f64-dpy) * weight
  a3            = (dpx)         * (dpy)         * weight
  a4            = (1.0_f64-dpx) * (dpy)         * weight
  f%r0(i,j)     = f%r0(i,j)     + a1
  f%r0(i+1,j)   = f%r0(i+1,j)   + a2
  f%r0(i+1,j+1) = f%r0(i+1,j+1) + a3
  f%r0(i,j+1)   = f%r0(i,j+1)   + a4

end do

f%r0 = f%r0 / (dx*dy)

rho_total = sum(f%r0(0:nx-1,0:ny-1))*dx*dy
f%r0 = f%r0 - rho_total/dimx/dimy

f%r0(0:nx-1,ny)  = f%r0(0:nx-1,0)
f%r0(nx,0:ny-1)  = f%r0(0,0:ny-1)
f%r0(nx,ny)      = f%r0(0,0)

end subroutine calcul_rho_cic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> M4 function from Monhagan (SPH method)
!> Cubic spline
!>             |  1 - 1.5 q^2 + 0.75 q^3  for 0 <= q <= 1
!> M4(x) = 2/3 |  1/4 (2−q)^3,            for 1 <= q <= 2
!>             |  0                       for q > 2.

function f_m4( x )
sll_real32, intent(in) :: x
sll_real32             :: f_m4

if ( x < 1. ) then
   f_m4 = 1 - 1.5*x*x+0.75*x*x*x
else if ( x >= 1. .and. x < 2. ) then
   f_m4 = 0.25*(2.-x)**3
else
   f_m4 = 0.0
end if

f_m4 = f_m4 * 2.0 / 3.0

return
end function f_m4

subroutine interpol_eb_m4( f, ele )

type(particle)       :: ele
type(mesh_fields)    :: f
sll_int32            :: k
sll_int32            :: i, j
sll_int32            :: im1, im2, ip1, ip2
sll_int32            :: jm1, jm2, jp1, jp2
sll_real32           :: dpx
sll_real32           :: dpy
sll_real32           :: cm2x
sll_real32           :: cp2x
sll_real32           :: cm1x
sll_real32           :: cp1x
sll_real32           :: cx
sll_real32           :: cy
sll_real32           :: cm2y
sll_real32           :: cp2y
sll_real32           :: cm1y
sll_real32           :: cp1y

do k=1,nbpart

   i   = ele%idx(k)
   j   = ele%idy(k)

   im2 = modulo(i-2,nx)
   im1 = modulo(i-1,nx)
   ip1 = modulo(i+1,nx)
   ip2 = modulo(i+2,nx)
   jm2 = modulo(j-2,ny)
   jm1 = modulo(j-1,ny)
   jp1 = modulo(j+1,ny)
   jp2 = modulo(j+2,ny)

   dpx = ele%dpx(k)
   dpy = ele%dpy(k)

   cm2x = f_m4(2.+dpx)
   cp2x = f_m4(2.-dpx)
   cm1x = f_m4(1.+dpx)
   cp1x = f_m4(1.-dpx)
   cx   = f_m4(dpx)
   cy   = f_m4(dpy)
   cp1y = f_m4(1.-dpy)
   cm1y = f_m4(1.+dpy)
   cp2y = f_m4(2.-dpy)
   cm2y = f_m4(2.+dpy)

   ele%epx(k) =                                 &
   &             + cm2x * cm2y * f%ex(im2,jm2)   &
   &             + cm2x * cm1y * f%ex(im2,jm1)   &
   &             + cm2x * cy   * f%ex(im2,j  )   &
   &             + cm2x * cp1y * f%ex(im2,jp1)   &
   &             + cm2x * cp2y * f%ex(im2,jp2)   &
   &             + cm1x * cm2y * f%ex(im1,jm2)   &
   &             + cm1x * cm1y * f%ex(im1,jm1)   &
   &             + cm1x * cy   * f%ex(im1,j  )   &
   &             + cm1x * cp1y * f%ex(im1,j+1)   &
   &             + cm1x * cp2y * f%ex(im1,jp2)   &
   &             + cx   * cm2y * f%ex(i  ,jm2)   &
   &             + cx   * cm1y * f%ex(i  ,jm1)   &
   &             + cx   * cy   * f%ex(i  ,j  )   &
   &             + cx   * cp1y * f%ex(i  ,jp1)   &
   &             + cx   * cp2y * f%ex(i  ,jp2)   &
   &             + cp1x * cm2y * f%ex(ip1,jm2)   &
   &             + cp1x * cm1y * f%ex(ip1,jm1)   &
   &             + cp1x * cy   * f%ex(ip1,j  )   &
   &             + cp1x * cp1y * f%ex(ip1,jp1)   &
   &             + cp1x * cp2y * f%ex(ip1,jp2)   &
   &             + cp2x * cm2y * f%ex(ip2,jm2)   &
   &             + cp2x * cm1y * f%ex(ip2,jm1)   &
   &             + cp2x * cy   * f%ex(ip2,j  )   &
   &             + cp2x * cp1y * f%ex(ip2,jp1)   &
   &             + cp2x * cp2y * f%ex(ip2,jp2)

   ele%epy(k) =                                 &
   &             + cm2x * cm2y * f%ey(im2,jm2)   &
   &             + cm2x * cm1y * f%ey(im2,jm1)   &
   &             + cm2x * cy   * f%ey(im2,j  )   &
   &             + cm2x * cp1y * f%ey(im2,jp1)   &
   &             + cm2x * cp2y * f%ey(im2,jp2)   &
   &             + cm1x * cm2y * f%ey(im1,jm2)   &
   &             + cm1x * cm1y * f%ey(im1,jm1)   &
   &             + cm1x * cy   * f%ey(im1,j  )   &
   &             + cm1x * cp1y * f%ey(im1,jp1)   &
   &             + cm1x * cp2y * f%ey(im1,jp2)   &
   &             + cx   * cm2y * f%ey(i  ,jm2)   &
   &             + cx   * cm1y * f%ey(i  ,jm1)   &
   &             + cx   * cy   * f%ey(i  ,j  )   &
   &             + cx   * cp1y * f%ey(i  ,jp1)   &
   &             + cx   * cp2y * f%ey(i  ,jp2)   &
   &             + cp1x * cm2y * f%ey(ip1,jm2)   &
   &             + cp1x * cm1y * f%ey(ip1,jm1)   &
   &             + cp1x * cy   * f%ey(ip1,j  )   &
   &             + cp1x * cp1y * f%ey(ip1,jp1)   &
   &             + cp1x * cp2y * f%ey(ip1,jp2)   &
   &             + cp2x * cm2y * f%ey(ip2,jm2)   &
   &             + cp2x * cm1y * f%ey(ip2,jm1)   &
   &             + cp2x * cy   * f%ey(ip2,j  )   &
   &             + cp2x * cp1y * f%ey(ip2,jp1)   &
   &             + cp2x * cp2y * f%ey(ip2,jp2)

end do

end subroutine interpol_eb_m4

subroutine calcul_rho_m4( ele, f )

type(particle)    :: ele
type(mesh_fields) :: f

sll_int32  :: i, j, k
sll_int32  :: im1, im2, ip1, ip2
sll_int32  :: jm1, jm2, jp1, jp2
sll_real32 :: dpx, dpy
sll_real32 :: cx, cm1x, cm2x, cp1x, cp2x
sll_real32 :: cy, cm1y, cm2y, cp1y, cp2y
sll_real64 :: weight, rho_total

f%r0 = 0.0_f64

do k = 1, nbpart

  i      = ele%idx(k)
  j      = ele%idy(k)
  dpx    = ele%dpx(k)
  dpy    = ele%dpy(k)
  weight = ele%p(k)

  im2 = modulo(i-2,nx)
  im1 = modulo(i-1,nx)
  ip1 = modulo(i+1,nx)
  ip2 = modulo(i+2,nx)
  jm2 = modulo(j-2,ny)
  jm1 = modulo(j-1,ny)
  jp1 = modulo(j+1,ny)
  jp2 = modulo(j+2,ny)

  cm2x = f_m4(2.+dpx)
  cp2x = f_m4(2.-dpx)
  cm1x = f_m4(1.+dpx)
  cp1x = f_m4(1.-dpx)
  cx   = f_m4(dpx)
  cy   = f_m4(dpy)
  cm2y = f_m4(2.+dpy)
  cp2y = f_m4(2.-dpy)
  cm1y = f_m4(1.+dpy)
  cp1y = f_m4(1.-dpy)

  f%r0(im2,jm2) = f%r0(im2,jm2) + cm2x * cm2y * weight
  f%r0(im2,jm1) = f%r0(im2,jm1) + cm2x * cm1y * weight
  f%r0(im2,j  ) = f%r0(im2,j  ) + cm2x * cy   * weight
  f%r0(im2,jp1) = f%r0(im2,jp1) + cm2x * cp1y * weight
  f%r0(im2,jp2) = f%r0(im2,jp2) + cm2x * cp2y * weight

  f%r0(im1,jm2) = f%r0(im1,jm2) + cm1x * cm2y * weight
  f%r0(im1,jm1) = f%r0(im1,jm1) + cm1x * cm1y * weight
  f%r0(im1,j  ) = f%r0(im1,j  ) + cm1x * cy   * weight
  f%r0(im1,jp1) = f%r0(im1,jp1) + cm1x * cp1y * weight
  f%r0(im1,jp2) = f%r0(im1,jp2) + cm1x * cp2y * weight

  f%r0(i  ,jm2) = f%r0(i  ,jm2) + cx   * cm2y * weight
  f%r0(i  ,jm1) = f%r0(i  ,jm1) + cx   * cm1y * weight
  f%r0(i  ,j  ) = f%r0(i  ,j  ) + cx   * cy   * weight
  f%r0(i  ,jp1) = f%r0(i  ,jp1) + cx   * cp1y * weight
  f%r0(i  ,jp2) = f%r0(i  ,jp2) + cx   * cp2y * weight

  f%r0(ip1,jm2) = f%r0(ip1,jm2) + cp1x * cm2y * weight
  f%r0(ip1,jm1) = f%r0(ip1,jm1) + cp1x * cm1y * weight
  f%r0(ip1,j  ) = f%r0(ip1,j  ) + cp1x * cy   * weight
  f%r0(ip1,jp1) = f%r0(ip1,jp1) + cp1x * cp1y * weight
  f%r0(ip1,jp2) = f%r0(ip1,jp2) + cp1x * cp2y * weight

  f%r0(ip2,jm2) = f%r0(ip2,jm2) + cp2x * cm2y * weight
  f%r0(ip2,jm1) = f%r0(ip2,jm1) + cp2x * cm1y * weight
  f%r0(ip2,j  ) = f%r0(ip2,j  ) + cp2x * cy   * weight
  f%r0(ip2,jp1) = f%r0(ip2,jp1) + cp2x * cp1y * weight
  f%r0(ip2,jp2) = f%r0(ip2,jp2) + cp2x * cp2y * weight

end do

!periodic boundary conditions

f%r0 = f%r0 / (dx*dy)

rho_total = sum(f%r0(0:nx-1,0:ny-1))*dx*dy
f%r0 = f%r0 - rho_total/dimx/dimy

f%r0(0:nx-1,ny)  = f%r0(0:nx-1,0)
f%r0(nx,0:ny-1)  = f%r0(0,0:ny-1)
f%r0(nx,ny)      = f%r0(0,0)

end subroutine calcul_rho_m4

!Bspline weighting function
!        | 0.5*(1.5+x)^2        -1.5 <  x < -0.5
! W(x) = | 0.75* x^2            -0.5 <= x <  0.5
!        | 0.5*(1.5-x)^2         0.5 <= x <  1.5

#define FONCTION1( X ) (0.75-(X)*(X))

#define FONCTION2( X ) (0.5 * ( 1.5 - (X) )**2)

#define BSPLINES(X,Y) \
c_1x = FONCTION2(1+X); \
c1x  = FONCTION1(X)  ; \
c2x  = FONCTION2(1-X); \
c_1y = FONCTION2(1+Y); \
c1y  = FONCTION1(Y)  ; \
c2y  = FONCTION2(1-Y)


subroutine calcul_rho_m3( ele, f)
type(particle)    :: ele
type(mesh_fields) :: f

sll_int32  :: i, j, k
sll_real32 :: dpx, dpy
sll_real64 :: weight, rho_total
sll_real32 :: c1x, c_1x, c2x
sll_real32 :: c1y, c_1y, c2y


f%r0 = 0.0_f64

do k = 1, nbpart

  i      = ele%idx(k)
  j      = ele%idy(k)
  dpx    = ele%dpx(k)
  dpy    = ele%dpy(k)
  weight = ele%p(k)

  BSPLINES(dpx,dpy)

  f%r0( i  ,j  ) = f%r0( i  ,j  ) + c1x*c1y   * weight
  f%r0( i  ,j-1) = f%r0( i  ,j-1) + c1x*c_1y  * weight
  f%r0( i  ,j+1) = f%r0( i  ,j+1) + c1x*c2y   * weight
  f%r0( i-1,j  ) = f%r0( i-1,j  ) + c_1x*c1y  * weight
  f%r0( i-1,j-1) = f%r0( i-1,j-1) + c_1x*c_1y * weight
  f%r0( i-1,j+1) = f%r0( i-1,j+1) + c_1x*c2y  * weight
  f%r0( i+1,j  ) = f%r0( i+1,j  ) + c2x*c1y   * weight
  f%r0( i+1,j-1) = f%r0( i+1,j-1) + c2x*c_1y  * weight
  f%r0( i+1,j+1) = f%r0( i+1,j+1) + c2x*c2y   * weight

end do

!periodic boundary conditions
f%r0 = f%r0 / (dx*dy)
rho_total = sum(f%r0(0:nx-1,0:ny-1))*dx*dy
f%r0 = f%r0 - rho_total/dimx/dimy
f%r0(0:nx-1,ny)  = f%r0(0:nx-1,0)
f%r0(nx,0:ny-1)  = f%r0(0,0:ny-1)
f%r0(nx,ny)      = f%r0(0,0)

end subroutine calcul_rho_m3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> M6 function 
!> Quintic spline
!>               |  (3-q)^5-6(2-q)^5+15(1-q)^5  for 0 <= q < 1
!> M6(x) = 1/120 |  (3-q)^5-6(2-q)^5            for 1 <= q < 2
!>               |  (3-q)^5                     for 2 <= q < 3
!>               |  0                           for q >= 3

function f_m6( q )
sll_real64, intent(in) :: q
sll_real64             :: f_m6

f_m6 = 0.0d0
if ( q < 3.0d0 ) f_m6 = (3.0d0-q)**5
if ( q < 2.0d0 ) f_m6 = f_m6 - 6.0d0*(2.0d0-q)**5
if ( q < 1.0d0 ) f_m6 = f_m6 + 15.0d0*(1.0d0-q)**5

f_m6 = f_m6 * 0.008333333333333333_f64 !/ 120.0d0

return
end function f_m6

subroutine calcul_rho_m6( ele, tm )

type(particle)    :: ele
type(mesh_fields) :: tm

sll_int32  :: i, j, k
sll_int32  :: im1, im2, im3, ip1, ip2, ip3
sll_int32  :: jm1, jm2, jm3, jp1, jp2, jp3
sll_real64 :: dpx, dpy
sll_real64 :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
sll_real64 :: cy, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y
sll_real64 :: weight, rho_total

tm%r0 = 0.0_f64

do k = 1, nbpart

  i      = ele%idx(k)
  j      = ele%idy(k)
  dpx    = ele%dpx(k)
  dpy    = ele%dpy(k)
  weight = ele%p(k)

  im3 = modulo(i-3,nx)
  im2 = modulo(i-2,nx)
  im1 = modulo(i-1,nx)
  ip1 = modulo(i+1,nx)
  ip2 = modulo(i+2,nx)
  ip3 = modulo(i+3,nx)
  jm3 = modulo(j-3,ny)
  jm2 = modulo(j-2,ny)
  jm1 = modulo(j-1,ny)
  jp1 = modulo(j+1,ny)
  jp2 = modulo(j+2,ny)
  jp3 = modulo(j+3,ny)

  cm3x = f_m6(3.0d0+dpx)
  cp3x = f_m6(3.0d0-dpx)
  cm2x = f_m6(2.0d0+dpx)
  cp2x = f_m6(2.0d0-dpx)
  cm1x = f_m6(1.0d0+dpx)
  cp1x = f_m6(1.0d0-dpx)
  cx   = f_m6(dpx)
  cy   = f_m6(dpy)
  cp1y = f_m6(1.0d0-dpy)
  cm1y = f_m6(1.0d0+dpy)
  cp2y = f_m6(2.0d0-dpy)
  cm2y = f_m6(2.0d0+dpy)
  cp3y = f_m6(3.0d0-dpy)
  cm3y = f_m6(3.0d0+dpy)

  tm%r0(im3,jm3) = tm%r0(im3,jm3) + cm3x * cm3y * weight
  tm%r0(im3,jm2) = tm%r0(im3,jm2) + cm3x * cm2y * weight
  tm%r0(im3,jm1) = tm%r0(im3,jm1) + cm3x * cm1y * weight
  tm%r0(im3,j  ) = tm%r0(im3,j  ) + cm3x * cy   * weight
  tm%r0(im3,jp1) = tm%r0(im3,jp1) + cm3x * cp1y * weight
  tm%r0(im3,jp2) = tm%r0(im3,jp2) + cm3x * cp2y * weight
  tm%r0(im3,jp3) = tm%r0(im3,jp3) + cm3x * cp3y * weight

  tm%r0(im2,jm3) = tm%r0(im2,jm3) + cm2x * cm3y * weight
  tm%r0(im2,jm2) = tm%r0(im2,jm2) + cm2x * cm2y * weight
  tm%r0(im2,jm1) = tm%r0(im2,jm1) + cm2x * cm1y * weight
  tm%r0(im2,j  ) = tm%r0(im2,j  ) + cm2x * cy   * weight
  tm%r0(im2,jp1) = tm%r0(im2,jp1) + cm2x * cp1y * weight
  tm%r0(im2,jp2) = tm%r0(im2,jp2) + cm2x * cp2y * weight
  tm%r0(im2,jp3) = tm%r0(im2,jp3) + cm2x * cp3y * weight

  tm%r0(im1,jm3) = tm%r0(im1,jm3) + cm1x * cm3y * weight
  tm%r0(im1,jm2) = tm%r0(im1,jm2) + cm1x * cm2y * weight
  tm%r0(im1,jm1) = tm%r0(im1,jm1) + cm1x * cm1y * weight
  tm%r0(im1,j  ) = tm%r0(im1,j  ) + cm1x * cy   * weight
  tm%r0(im1,jp1) = tm%r0(im1,jp1) + cm1x * cp1y * weight
  tm%r0(im1,jp2) = tm%r0(im1,jp2) + cm1x * cp2y * weight
  tm%r0(im1,jp3) = tm%r0(im1,jp3) + cm1x * cp3y * weight

  tm%r0(i  ,jm3) = tm%r0(i  ,jm3) + cx   * cm3y * weight
  tm%r0(i  ,jm2) = tm%r0(i  ,jm2) + cx   * cm2y * weight
  tm%r0(i  ,jm1) = tm%r0(i  ,jm1) + cx   * cm1y * weight
  tm%r0(i  ,j  ) = tm%r0(i  ,j  ) + cx   * cy   * weight
  tm%r0(i  ,jp1) = tm%r0(i  ,jp1) + cx   * cp1y * weight
  tm%r0(i  ,jp2) = tm%r0(i  ,jp2) + cx   * cp2y * weight
  tm%r0(i  ,jp3) = tm%r0(i  ,jp3) + cx   * cp3y * weight

  tm%r0(ip1,jm3) = tm%r0(ip1,jm3) + cp1x * cm3y * weight
  tm%r0(ip1,jm2) = tm%r0(ip1,jm2) + cp1x * cm2y * weight
  tm%r0(ip1,jm1) = tm%r0(ip1,jm1) + cp1x * cm1y * weight
  tm%r0(ip1,j  ) = tm%r0(ip1,j  ) + cp1x * cy   * weight
  tm%r0(ip1,jp1) = tm%r0(ip1,jp1) + cp1x * cp1y * weight
  tm%r0(ip1,jp2) = tm%r0(ip1,jp2) + cp1x * cp2y * weight
  tm%r0(ip1,jp3) = tm%r0(ip1,jp3) + cp1x * cp3y * weight

  tm%r0(ip2,jm3) = tm%r0(ip2,jm3) + cp2x * cm3y * weight
  tm%r0(ip2,jm2) = tm%r0(ip2,jm2) + cp2x * cm2y * weight
  tm%r0(ip2,jm1) = tm%r0(ip2,jm1) + cp2x * cm1y * weight
  tm%r0(ip2,j  ) = tm%r0(ip2,j  ) + cp2x * cy   * weight
  tm%r0(ip2,jp1) = tm%r0(ip2,jp1) + cp2x * cp1y * weight
  tm%r0(ip2,jp2) = tm%r0(ip2,jp2) + cp2x * cp2y * weight
  tm%r0(ip2,jp3) = tm%r0(ip2,jp3) + cp2x * cp3y * weight

  tm%r0(ip3,jm3) = tm%r0(ip3,jm3) + cp3x * cm3y * weight
  tm%r0(ip3,jm2) = tm%r0(ip3,jm2) + cp3x * cm2y * weight
  tm%r0(ip3,jm1) = tm%r0(ip3,jm1) + cp3x * cm1y * weight
  tm%r0(ip3,j  ) = tm%r0(ip3,j  ) + cp3x * cy   * weight
  tm%r0(ip3,jp1) = tm%r0(ip3,jp1) + cp3x * cp1y * weight
  tm%r0(ip3,jp2) = tm%r0(ip3,jp2) + cp3x * cp2y * weight
  tm%r0(ip3,jp3) = tm%r0(ip3,jp3) + cp3x * cp3y * weight

end do


tm%r0 = tm%r0 / (dx*dy)

rho_total = sum(tm%r0(0:nx-1,0:ny-1))*dx*dy
tm%r0 = tm%r0 - rho_total/dimx/dimy

tm%r0(0:nx-1,ny)  = tm%r0(0:nx-1,0)
tm%r0(nx,0:ny-1)  = tm%r0(0,0:ny-1)
tm%r0(nx,ny)      = tm%r0(0,0)

end subroutine calcul_rho_m6

subroutine interpol_eb_m6( tm, ele )

type(particle)       :: ele
type(mesh_fields)    :: tm
sll_int32            :: k
sll_int32            :: i, j
sll_int32            :: im1, im2, im3, ip1, ip2, ip3
sll_int32            :: jm1, jm2, jm3, jp1, jp2, jp3
sll_real64           :: dpx
sll_real64           :: dpy
sll_real64           :: cm3x
sll_real64           :: cp3x
sll_real64           :: cm2x
sll_real64           :: cp2x
sll_real64           :: cm1x
sll_real64           :: cp1x
sll_real64           :: cx
sll_real64           :: cy
sll_real64           :: cm3y
sll_real64           :: cp3y
sll_real64           :: cm2y
sll_real64           :: cp2y
sll_real64           :: cm1y
sll_real64           :: cp1y
sll_real64           :: cm3x_cm3y
sll_real64           :: cm3x_cm2y
sll_real64           :: cm3x_cm1y
sll_real64           :: cm3x_cy  
sll_real64           :: cm3x_cp1y
sll_real64           :: cm3x_cp2y
sll_real64           :: cm3x_cp3y
sll_real64           :: cm2x_cm3y
sll_real64           :: cm2x_cm2y
sll_real64           :: cm2x_cm1y
sll_real64           :: cm2x_cy  
sll_real64           :: cm2x_cp1y
sll_real64           :: cm2x_cp2y
sll_real64           :: cm2x_cp3y
sll_real64           :: cm1x_cm3y
sll_real64           :: cm1x_cm2y
sll_real64           :: cm1x_cm1y
sll_real64           :: cm1x_cy  
sll_real64           :: cm1x_cp1y
sll_real64           :: cm1x_cp2y
sll_real64           :: cm1x_cp3y
sll_real64           :: cx_cm3y  
sll_real64           :: cx_cm2y  
sll_real64           :: cx_cm1y  
sll_real64           :: cx_cy    
sll_real64           :: cx_cp1y  
sll_real64           :: cx_cp2y  
sll_real64           :: cx_cp3y  
sll_real64           :: cp1x_cm3y
sll_real64           :: cp1x_cm2y
sll_real64           :: cp1x_cm1y
sll_real64           :: cp1x_cy  
sll_real64           :: cp1x_cp1y
sll_real64           :: cp1x_cp2y
sll_real64           :: cp1x_cp3y
sll_real64           :: cp2x_cm3y
sll_real64           :: cp2x_cm2y
sll_real64           :: cp2x_cm1y
sll_real64           :: cp2x_cy  
sll_real64           :: cp2x_cp1y
sll_real64           :: cp2x_cp2y
sll_real64           :: cp2x_cp3y
sll_real64           :: cp3x_cm3y
sll_real64           :: cp3x_cm2y
sll_real64           :: cp3x_cm1y
sll_real64           :: cp3x_cy  
sll_real64           :: cp3x_cp1y
sll_real64           :: cp3x_cp2y
sll_real64           :: cp3x_cp3y


do k=1,nbpart

   i   = ele%idx(k)
   j   = ele%idy(k)

   im3 = modulo(i-3,nx)
   im2 = modulo(i-2,nx)
   im1 = modulo(i-1,nx)
   ip1 = modulo(i+1,nx)
   ip2 = modulo(i+2,nx)
   ip3 = modulo(i+3,nx)
   jm3 = modulo(j-3,ny)
   jm2 = modulo(j-2,ny)
   jm1 = modulo(j-1,ny)
   jp1 = modulo(j+1,ny)
   jp2 = modulo(j+2,ny)
   jp3 = modulo(j+3,ny)

   dpx = ele%dpx(k)
   dpy = ele%dpy(k)

   cm3x = f_m6(3.0d0+dpx)
   cp3x = f_m6(3.0d0-dpx)
   cm2x = f_m6(2.0d0+dpx)
   cp2x = f_m6(2.0d0-dpx)
   cm1x = f_m6(1.0d0+dpx)
   cp1x = f_m6(1.0d0-dpx)
   cx   = f_m6(dpx)
   cy   = f_m6(dpy)
   cm3y = f_m6(3.0d0+dpy)
   cp3y = f_m6(3.0d0-dpy)
   cm2y = f_m6(2.0d0+dpy)
   cp2y = f_m6(2.0d0-dpy)
   cm1y = f_m6(1.0d0+dpy)
   cp1y = f_m6(1.0d0-dpy)

   cm3x_cm3y = cm3x * cm3y 
   cm3x_cm2y = cm3x * cm2y 
   cm3x_cm1y = cm3x * cm1y 
   cm3x_cy   = cm3x * cy   
   cm3x_cp1y = cm3x * cp1y 
   cm3x_cp2y = cm3x * cp2y 
   cm3x_cp3y = cm3x * cp3y 
   cm2x_cm3y = cm2x * cm3y 
   cm2x_cm2y = cm2x * cm2y 
   cm2x_cm1y = cm2x * cm1y 
   cm2x_cy   = cm2x * cy   
   cm2x_cp1y = cm2x * cp1y 
   cm2x_cp2y = cm2x * cp2y 
   cm2x_cp3y = cm2x * cp3y 
   cm1x_cm3y = cm1x * cm3y 
   cm1x_cm2y = cm1x * cm2y 
   cm1x_cm1y = cm1x * cm1y 
   cm1x_cy   = cm1x * cy   
   cm1x_cp1y = cm1x * cp1y 
   cm1x_cp2y = cm1x * cp2y 
   cm1x_cp3y = cm1x * cp3y 
   cx_cm3y   = cx   * cm3y 
   cx_cm2y   = cx   * cm2y 
   cx_cm1y   = cx   * cm1y 
   cx_cy     = cx   * cy   
   cx_cp1y   = cx   * cp1y 
   cx_cp2y   = cx   * cp2y 
   cx_cp3y   = cx   * cp3y 
   cp1x_cm3y = cp1x * cm3y 
   cp1x_cm2y = cp1x * cm2y 
   cp1x_cm1y = cp1x * cm1y 
   cp1x_cy   = cp1x * cy   
   cp1x_cp1y = cp1x * cp1y 
   cp1x_cp2y = cp1x * cp2y 
   cp1x_cp3y = cp1x * cp3y 
   cp2x_cm3y = cp2x * cm3y 
   cp2x_cm2y = cp2x * cm2y 
   cp2x_cm1y = cp2x * cm1y 
   cp2x_cy   = cp2x * cy   
   cp2x_cp1y = cp2x * cp1y 
   cp2x_cp2y = cp2x * cp2y 
   cp2x_cp3y = cp2x * cp3y 
   cp3x_cm3y = cp3x * cm3y 
   cp3x_cm2y = cp3x * cm2y 
   cp3x_cm1y = cp3x * cm1y 
   cp3x_cy   = cp3x * cy   
   cp3x_cp1y = cp3x * cp1y 
   cp3x_cp2y = cp3x * cp2y 
   cp3x_cp3y = cp3x * cp3y


   ele%epx(k) =                               &
   &             + cm3x_cm3y * tm%ex(im3,jm3)   &
   &             + cm3x_cm2y * tm%ex(im3,jm2)   &
   &             + cm3x_cm1y * tm%ex(im3,jm1)   &
   &             + cm3x_cy   * tm%ex(im3,j  )   &
   &             + cm3x_cp1y * tm%ex(im3,jp1)   &
   &             + cm3x_cp2y * tm%ex(im3,jp2)   &
   &             + cm3x_cp3y * tm%ex(im3,jp3)   &
   &             + cm2x_cm3y * tm%ex(im2,jm3)   &
   &             + cm2x_cm2y * tm%ex(im2,jm2)   &
   &             + cm2x_cm1y * tm%ex(im2,jm1)   &
   &             + cm2x_cy   * tm%ex(im2,j  )   &
   &             + cm2x_cp1y * tm%ex(im2,jp1)   &
   &             + cm2x_cp2y * tm%ex(im2,jp2)   &
   &             + cm2x_cp3y * tm%ex(im2,jp3)   &
   &             + cm1x_cm3y * tm%ex(im1,jm3)   &
   &             + cm1x_cm2y * tm%ex(im1,jm2)   &
   &             + cm1x_cm1y * tm%ex(im1,jm1)   &
   &             + cm1x_cy   * tm%ex(im1,j  )   &
   &             + cm1x_cp1y * tm%ex(im1,jp1)   &
   &             + cm1x_cp2y * tm%ex(im1,jp2)   &
   &             + cm1x_cp3y * tm%ex(im1,jp3)   &
   &             + cx_cm3y   * tm%ex(i  ,jm3)   &
   &             + cx_cm2y   * tm%ex(i  ,jm2)   &
   &             + cx_cm1y   * tm%ex(i  ,jm1)   &
   &             + cx_cy     * tm%ex(i  ,j  )   &
   &             + cx_cp1y   * tm%ex(i  ,jp1)   &
   &             + cx_cp2y   * tm%ex(i  ,jp2)   &
   &             + cx_cp3y   * tm%ex(i  ,jp3)   &
   &             + cp1x_cm3y * tm%ex(ip1,jm3)   &
   &             + cp1x_cm2y * tm%ex(ip1,jm2)   &
   &             + cp1x_cm1y * tm%ex(ip1,jm1)   &
   &             + cp1x_cy   * tm%ex(ip1,j  )   &
   &             + cp1x_cp1y * tm%ex(ip1,jp1)   &
   &             + cp1x_cp2y * tm%ex(ip1,jp2)   &
   &             + cp1x_cp3y * tm%ex(ip1,jp3)   &
   &             + cp2x_cm3y * tm%ex(ip2,jm3)   &
   &             + cp2x_cm2y * tm%ex(ip2,jm2)   &
   &             + cp2x_cm1y * tm%ex(ip2,jm1)   &
   &             + cp2x_cy   * tm%ex(ip2,j  )   &
   &             + cp2x_cp1y * tm%ex(ip2,jp1)   &
   &             + cp2x_cp2y * tm%ex(ip2,jp2)   &
   &             + cp2x_cp3y * tm%ex(ip2,jp3)   &
   &             + cp3x_cm3y * tm%ex(ip3,jm3)   &
   &             + cp3x_cm2y * tm%ex(ip3,jm2)   &
   &             + cp3x_cm1y * tm%ex(ip3,jm1)   &
   &             + cp3x_cy   * tm%ex(ip3,j  )   &
   &             + cp3x_cp1y * tm%ex(ip3,jp1)   &
   &             + cp3x_cp2y * tm%ex(ip3,jp2)   &
   &             + cp3x_cp3y * tm%ex(ip3,jp3) 

   ele%epy(k) =                               &
   &             + cm3x_cm3y * tm%ey(im3,jm3)   &
   &             + cm3x_cm2y * tm%ey(im3,jm2)   &
   &             + cm3x_cm1y * tm%ey(im3,jm1)   &
   &             + cm3x_cy   * tm%ey(im3,j  )   &
   &             + cm3x_cp1y * tm%ey(im3,jp1)   &
   &             + cm3x_cp2y * tm%ey(im3,jp2)   &
   &             + cm3x_cp3y * tm%ey(im3,jp3)   &
   &             + cm2x_cm3y * tm%ey(im2,jm3)   &
   &             + cm2x_cm2y * tm%ey(im2,jm2)   &
   &             + cm2x_cm1y * tm%ey(im2,jm1)   &
   &             + cm2x_cy   * tm%ey(im2,j  )   &
   &             + cm2x_cp1y * tm%ey(im2,jp1)   &
   &             + cm2x_cp2y * tm%ey(im2,jp2)   &
   &             + cm2x_cp3y * tm%ey(im2,jp3)   &
   &             + cm1x_cm3y * tm%ey(im1,jm3)   &
   &             + cm1x_cm2y * tm%ey(im1,jm2)   &
   &             + cm1x_cm1y * tm%ey(im1,jm1)   &
   &             + cm1x_cy   * tm%ey(im1,j  )   &
   &             + cm1x_cp1y * tm%ey(im1,jp1)   &
   &             + cm1x_cp2y * tm%ey(im1,jp2)   &
   &             + cm1x_cp3y * tm%ey(im1,jp3)   &
   &             + cx_cm3y   * tm%ey(i  ,jm3)   &
   &             + cx_cm2y   * tm%ey(i  ,jm2)   &
   &             + cx_cm1y   * tm%ey(i  ,jm1)   &
   &             + cx_cy     * tm%ey(i  ,j  )   &
   &             + cx_cp1y   * tm%ey(i  ,jp1)   &
   &             + cx_cp2y   * tm%ey(i  ,jp2)   &
   &             + cx_cp3y   * tm%ey(i  ,jp3)   &
   &             + cp1x_cm3y * tm%ey(ip1,jm3)   &
   &             + cp1x_cm2y * tm%ey(ip1,jm2)   &
   &             + cp1x_cm1y * tm%ey(ip1,jm1)   &
   &             + cp1x_cy   * tm%ey(ip1,j  )   &
   &             + cp1x_cp1y * tm%ey(ip1,jp1)   &
   &             + cp1x_cp2y * tm%ey(ip1,jp2)   &
   &             + cp1x_cp3y * tm%ey(ip1,jp3)   &
   &             + cp2x_cm3y * tm%ey(ip2,jm3)   &
   &             + cp2x_cm2y * tm%ey(ip2,jm2)   &
   &             + cp2x_cm1y * tm%ey(ip2,jm1)   &
   &             + cp2x_cy   * tm%ey(ip2,j  )   &
   &             + cp2x_cp1y * tm%ey(ip2,jp1)   &
   &             + cp2x_cp2y * tm%ey(ip2,jp2)   &
   &             + cp2x_cp3y * tm%ey(ip2,jp3)   &
   &             + cp3x_cm3y * tm%ey(ip3,jm3)   &
   &             + cp3x_cm2y * tm%ey(ip3,jm2)   &
   &             + cp3x_cm1y * tm%ey(ip3,jm1)   &
   &             + cp3x_cy   * tm%ey(ip3,j  )   &
   &             + cp3x_cp1y * tm%ey(ip3,jp1)   &
   &             + cp3x_cp2y * tm%ey(ip3,jp2)   &
   &             + cp3x_cp3y * tm%ey(ip3,jp3) 

end do

end subroutine interpol_eb_m6

end module m_particules
