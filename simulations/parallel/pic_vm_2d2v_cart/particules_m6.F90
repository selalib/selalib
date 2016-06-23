module particules_m6
#include "sll_working_precision.h"
#include "sll_memory.h"

use particules

implicit none

private

public calcul_rho_m6
public interpol_eb_m6

contains

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

if ( q < 1.0d0 ) then
  f_m6 = (3.0d0-q)**5-6.0d0*(2.0d0-q)**5+15.0d0*(1.0d0-q)**5
else if ( q >= 1.0d0 .and. q < 2.0d0 ) then
  f_m6 = (3.0d0-q)**5-6.0d0*(2.0d0-q)**5
else if ( q >= 2.0d0 .and. q < 3.0d0 ) then
  f_m6 = (3.0d0-q)**5
else
  f_m6 = 0.0d0
end if

f_m6 = f_m6 / 120.0d0

return
end function f_m6

subroutine calcul_rho_m6( ele, tm )

type(particle)       :: ele
type(tm_mesh_fields) :: tm

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

tm%r0(0:nx-1,ny)  = tm%r0(0:nx-1,0)
tm%r0(nx,0:ny-1)  = tm%r0(0,0:ny-1)
tm%r0(nx,ny)      = tm%r0(0,0)

tm%r0 = tm%r0 / (dx*dy)

rho_total = sum(tm%r0(0:nx-1,0:ny-1))*dx*dy
!print*,'rho total',rho_total
tm%r0 = tm%r0 - rho_total/dimx/dimy

end subroutine calcul_rho_m6

subroutine interpol_eb_m6( tm1, ele )

type(particle)       :: ele
type(tm_mesh_fields) :: tm1
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
sll_real64           :: e(-2:2,-2:2)

sll_real64, allocatable :: ehx(:,:)
sll_real64, allocatable :: ehy(:,:)
!sll_real64, allocatable :: bhz(:,:)

e(-2,-2) =    1.0_f64/14400.0_f64
e(-1,-2) =   26.0_f64/14400.0_f64
e( 0,-2) =   66.0_f64/14400.0_f64
e(+1,-2) =   26.0_f64/14400.0_f64
e(+2,-2) =    1.0_f64/14400.0_f64

e(-2,-1) =   26.0_f64/14400.0_f64
e(-1,-1) =  676.0_f64/14400.0_f64
e( 0,-1) = 1716.0_f64/14400.0_f64
e(+1,-1) =  676.0_f64/14400.0_f64
e(+2,-1) =   26.0_f64/14400.0_f64

e(-2, 0) =   66.0_f64/14400.0_f64
e(-1, 0) = 1716.0_f64/14400.0_f64
e( 0, 0) = 4356.0_f64/14400.0_f64
e(+1, 0) = 1716.0_f64/14400.0_f64
e(+2, 0) =   66.0_f64/14400.0_f64

e(-2,+1) =   26.0_f64/14400.0_f64
e(-1,+1) =  676.0_f64/14400.0_f64
e( 0,+1) = 1716.0_f64/14400.0_f64
e(+1,+1) =  676.0_f64/14400.0_f64
e(+2,+1) =   26.0_f64/14400.0_f64

e(-2,+2) =    1.0_f64/14400.0_f64
e(-1,+2) =   26.0_f64/14400.0_f64
e( 0,+2) =   66.0_f64/14400.0_f64
e(+1,+2) =   26.0_f64/14400.0_f64
e(+2,+2) =    1.0_f64/14400.0_f64

allocate(ehx(0:nx-1,0:ny-1))
allocate(ehy(0:nx-1,0:ny-1))
!allocate(bhz(0:nx-1,0:ny-1))

do j = 0,ny-1

  jm2 = modulo(j-2,ny)
  jm1 = modulo(j-1,ny)
  jp1 = modulo(j+1,ny)
  jp2 = modulo(j+2,ny)

  do i = 0,nx-1

     im2 = modulo(i-2,nx)
     im1 = modulo(i-1,nx)
     ip1 = modulo(i+1,nx)
     ip2 = modulo(i+2,nx)

     ehx(i,j) = e(-2,-2)*tm1%ex(im2,jm2) & 
              + e(-2,-1)*tm1%ex(im2,jm1) & 
              + e(-2, 0)*tm1%ex(im2,j  ) &
              + e(-2,+1)*tm1%ex(im2,jp1) &
              + e(-2,+2)*tm1%ex(im2,jp2) &
              + e(-1,-2)*tm1%ex(im1,jm2) & 
              + e(-1,-1)*tm1%ex(im1,jm1) & 
              + e(-1, 0)*tm1%ex(im1,j  ) &
              + e(-1,+1)*tm1%ex(im1,jp1) &
              + e(-1,+2)*tm1%ex(im1,jp2) &
              + e( 0,-2)*tm1%ex(  i,jm2) &
              + e( 0,-1)*tm1%ex(  i,jm1) &
              + e( 0, 0)*tm1%ex(  i,j  ) &
              + e( 0,+1)*tm1%ex(  i,jp1) &
              + e( 0,+2)*tm1%ex(  i,jp2) &
              + e(+1,-2)*tm1%ex(ip1,jm2) &
              + e(+1,-1)*tm1%ex(ip1,jm1) &
              + e(+1, 0)*tm1%ex(ip1,j  ) &
              + e(+1,+1)*tm1%ex(ip1,jp1) &
              + e(+1,+2)*tm1%ex(ip1,jp2) &
              + e(+2,-2)*tm1%ex(ip2,jm2) &
              + e(+2,-1)*tm1%ex(ip2,jm1) &
              + e(+2, 0)*tm1%ex(ip2,j  ) &
              + e(+2,+1)*tm1%ex(ip2,jp1) &
              + e(+2,+2)*tm1%ex(ip2,jp2)

     ehy(i,j) = e(-2,-2)*tm1%ey(im2,jm2) & 
              + e(-2,-1)*tm1%ey(im2,jm1) & 
              + e(-2, 0)*tm1%ey(im2,j  ) &
              + e(-2,+1)*tm1%ey(im2,jp1) &
              + e(-2,+2)*tm1%ey(im2,jp2) &
              + e(-1,-2)*tm1%ey(im1,jm2) & 
              + e(-1,-1)*tm1%ey(im1,jm1) & 
              + e(-1, 0)*tm1%ey(im1,j  ) &
              + e(-1,+1)*tm1%ey(im1,jp1) &
              + e(-1,+2)*tm1%ey(im1,jp2) &
              + e( 0,-2)*tm1%ey(  i,jm2) &
              + e( 0,-1)*tm1%ey(  i,jm1) &
              + e( 0, 0)*tm1%ey(  i,j  ) &
              + e( 0,+1)*tm1%ey(  i,jp1) &
              + e( 0,+2)*tm1%ey(  i,jp2) &
              + e(+1,-2)*tm1%ey(ip1,jm2) &
              + e(+1,-1)*tm1%ey(ip1,jm1) &
              + e(+1, 0)*tm1%ey(ip1,j  ) &
              + e(+1,+1)*tm1%ey(ip1,jp1) &
              + e(+1,+2)*tm1%ey(ip1,jp2) &
              + e(+2,-2)*tm1%ey(ip2,jm2) &
              + e(+2,-1)*tm1%ey(ip2,jm1) &
              + e(+2, 0)*tm1%ey(ip2,j  ) &
              + e(+2,+1)*tm1%ey(ip2,jp1) &
              + e(+2,+2)*tm1%ey(ip2,jp2)

     !bhz(i,j) = e(-2,-2)*tm1%bz(im2,jm2) &
!              + e(-2,-1)*tm1%bz(im2,jm1) & 
!              + e(-2, 0)*tm1%bz(im2,j  ) &
!              + e(-2,+1)*tm1%bz(im2,jp1) &
!              + e(-2,+2)*tm1%bz(im2,jp2) &
!              + e(-1,-2)*tm1%bz(im1,jm2) & 
!              + e(-1,-1)*tm1%bz(im1,jm1) & 
!              + e(-1, 0)*tm1%bz(im1,j  ) &
!              + e(-1,+1)*tm1%bz(im1,jp1) &
!              + e(-1,+2)*tm1%bz(im1,jp2) &
!              + e( 0,-2)*tm1%bz(  i,jm2) &
!              + e( 0,-1)*tm1%bz(  i,jm1) &
!              + e( 0, 0)*tm1%bz(  i,j  ) &
!              + e( 0,+1)*tm1%bz(  i,jp1) &
!              + e( 0,+2)*tm1%bz(  i,jp2) &
!              + e(+1,-2)*tm1%bz(ip1,jm2) &
!              + e(+1,-1)*tm1%bz(ip1,jm1) &
!              + e(+1, 0)*tm1%bz(ip1,j  ) &
!              + e(+1,+1)*tm1%bz(ip1,jp1) &
!              + e(+1,+2)*tm1%bz(ip1,jp2) &
!              + e(+2,-2)*tm1%bz(ip2,jm2) &
!              + e(+2,-1)*tm1%bz(ip2,jm1) &
!              + e(+2, 0)*tm1%bz(ip2,j  ) &
!              + e(+2,+1)*tm1%bz(ip2,jp1) &
!              + e(+2,+2)*tm1%bz(ip2,jp2)
  end do
end do

!print*, sum(abs(ehx)), sum(abs(tm1%ex(0:nx-1,0:ny-1)))
!print*, sum(abs(ehy)), sum(abs(tm1%ey(0:nx-1,0:ny-1)))
!print*, sum(abs(bhz)), sum(abs(tm1%bz(0:nx-1,0:ny-1)))

do k=1,nbpart

   i   = ele%idx(k)
   j   = ele%idy(k)
    i=modulo(i,nx)
    j=modulo(j,ny)

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

   ele%epx(k) =                                 &
   &             + cm3x * cm3y * ehx(im3,jm3)   &
   &             + cm3x * cm2y * ehx(im3,jm2)   &
   &             + cm3x * cm1y * ehx(im3,jm1)   &
   &             + cm3x * cy   * ehx(im3,j  )   &
   &             + cm3x * cp1y * ehx(im3,jp1)   &
   &             + cm3x * cp2y * ehx(im3,jp2)   &
   &             + cm3x * cp3y * ehx(im3,jp3)   &
   &             + cm2x * cm3y * ehx(im2,jm3)   &
   &             + cm2x * cm2y * ehx(im2,jm2)   &
   &             + cm2x * cm1y * ehx(im2,jm1)   &
   &             + cm2x * cy   * ehx(im2,j  )   &
   &             + cm2x * cp1y * ehx(im2,jp1)   &
   &             + cm2x * cp2y * ehx(im2,jp2)   &
   &             + cm2x * cp3y * ehx(im2,jp3)   &
   &             + cm1x * cm3y * ehx(im1,jm3)   &
   &             + cm1x * cm2y * ehx(im1,jm2)   &
   &             + cm1x * cm1y * ehx(im1,jm1)   &
   &             + cm1x * cy   * ehx(im1,j  )   &
   &             + cm1x * cp1y * ehx(im1,jp1)   &
   &             + cm1x * cp2y * ehx(im1,jp2)   &
   &             + cm1x * cp3y * ehx(im1,jp3)   &
   &             + cx   * cm3y * ehx(i  ,jm3)   &
   &             + cx   * cm2y * ehx(i  ,jm2)   &
   &             + cx   * cm1y * ehx(i  ,jm1)   &
   &             + cx   * cy   * ehx(i  ,j  )   &
   &             + cx   * cp1y * ehx(i  ,jp1)   &
   &             + cx   * cp2y * ehx(i  ,jp2)   &
   &             + cx   * cp3y * ehx(i  ,jp3)   &
   &             + cp1x * cm3y * ehx(ip1,jm3)   &
   &             + cp1x * cm2y * ehx(ip1,jm2)   &
   &             + cp1x * cm1y * ehx(ip1,jm1)   &
   &             + cp1x * cy   * ehx(ip1,j  )   &
   &             + cp1x * cp1y * ehx(ip1,jp1)   &
   &             + cp1x * cp2y * ehx(ip1,jp2)   &
   &             + cp1x * cp3y * ehx(ip1,jp3)   &
   &             + cp2x * cm3y * ehx(ip2,jm3)   &
   &             + cp2x * cm2y * ehx(ip2,jm2)   &
   &             + cp2x * cm1y * ehx(ip2,jm1)   &
   &             + cp2x * cy   * ehx(ip2,j  )   &
   &             + cp2x * cp1y * ehx(ip2,jp1)   &
   &             + cp2x * cp2y * ehx(ip2,jp2)   &
   &             + cp2x * cp3y * ehx(ip2,jp3)   &
   &             + cp3x * cm3y * ehx(ip3,jm3)   &
   &             + cp3x * cm2y * ehx(ip3,jm2)   &
   &             + cp3x * cm1y * ehx(ip3,jm1)   &
   &             + cp3x * cy   * ehx(ip3,j  )   &
   &             + cp3x * cp1y * ehx(ip3,jp1)   &
   &             + cp3x * cp2y * ehx(ip3,jp2)   &
   &             + cp3x * cp3y * ehx(ip3,jp3) 

   ele%epy(k) =                                 &
   &             + cm3x * cm3y * ehy(im3,jm3)   &
   &             + cm3x * cm2y * ehy(im3,jm2)   &
   &             + cm3x * cm1y * ehy(im3,jm1)   &
   &             + cm3x * cy   * ehy(im3,j  )   &
   &             + cm3x * cp1y * ehy(im3,jp1)   &
   &             + cm3x * cp2y * ehy(im3,jp2)   &
   &             + cm3x * cp3y * ehy(im3,jp3)   &
   &             + cm2x * cm3y * ehy(im2,jm3)   &
   &             + cm2x * cm2y * ehy(im2,jm2)   &
   &             + cm2x * cm1y * ehy(im2,jm1)   &
   &             + cm2x * cy   * ehy(im2,j  )   &
   &             + cm2x * cp1y * ehy(im2,jp1)   &
   &             + cm2x * cp2y * ehy(im2,jp2)   &
   &             + cm2x * cp3y * ehy(im2,jp3)   &
   &             + cm1x * cm3y * ehy(im1,jm3)   &
   &             + cm1x * cm2y * ehy(im1,jm2)   &
   &             + cm1x * cm1y * ehy(im1,jm1)   &
   &             + cm1x * cy   * ehy(im1,j  )   &
   &             + cm1x * cp1y * ehy(im1,jp1)   &
   &             + cm1x * cp2y * ehy(im1,jp2)   &
   &             + cm1x * cp3y * ehy(im1,jp3)   &
   &             + cx   * cm3y * ehy(i  ,jm3)   &
   &             + cx   * cm2y * ehy(i  ,jm2)   &
   &             + cx   * cm1y * ehy(i  ,jm1)   &
   &             + cx   * cy   * ehy(i  ,j  )   &
   &             + cx   * cp1y * ehy(i  ,jp1)   &
   &             + cx   * cp2y * ehy(i  ,jp2)   &
   &             + cx   * cp3y * ehy(i  ,jp3)   &
   &             + cp1x * cm3y * ehy(ip1,jm3)   &
   &             + cp1x * cm2y * ehy(ip1,jm2)   &
   &             + cp1x * cm1y * ehy(ip1,jm1)   &
   &             + cp1x * cy   * ehy(ip1,j  )   &
   &             + cp1x * cp1y * ehy(ip1,jp1)   &
   &             + cp1x * cp2y * ehy(ip1,jp2)   &
   &             + cp1x * cp3y * ehy(ip1,jp3)   &
   &             + cp2x * cm3y * ehy(ip2,jm3)   &
   &             + cp2x * cm2y * ehy(ip2,jm2)   &
   &             + cp2x * cm1y * ehy(ip2,jm1)   &
   &             + cp2x * cy   * ehy(ip2,j  )   &
   &             + cp2x * cp1y * ehy(ip2,jp1)   &
   &             + cp2x * cp2y * ehy(ip2,jp2)   &
   &             + cp2x * cp3y * ehy(ip2,jp3)   &
   &             + cp3x * cm3y * ehy(ip3,jm3)   &
   &             + cp3x * cm2y * ehy(ip3,jm2)   &
   &             + cp3x * cm1y * ehy(ip3,jm1)   &
   &             + cp3x * cy   * ehy(ip3,j  )   &
   &             + cp3x * cp1y * ehy(ip3,jp1)   &
   &             + cp3x * cp2y * ehy(ip3,jp2)   &
   &             + cp3x * cp3y * ehy(ip3,jp3) 

!   ele%bpz(k) =                                 &
!   &             + cm3x * cm3y * bhz(im3,jm3)   &
!   &             + cm3x * cm2y * bhz(im3,jm2)   &
!   &             + cm3x * cm1y * bhz(im3,jm1)   &
!   &             + cm3x * cy   * bhz(im3,j  )   &
!   &             + cm3x * cp1y * bhz(im3,jp1)   &
!   &             + cm3x * cp2y * bhz(im3,jp2)   &
!   &             + cm3x * cp3y * bhz(im3,jp3)   &
!   &             + cm2x * cm3y * bhz(im2,jm3)   &
!   &             + cm2x * cm2y * bhz(im2,jm2)   &
!   &             + cm2x * cm1y * bhz(im2,jm1)   &
!   &             + cm2x * cy   * bhz(im2,j  )   &
!   &             + cm2x * cp1y * bhz(im2,jp1)   &
!   &             + cm2x * cp2y * bhz(im2,jp2)   &
!   &             + cm2x * cp3y * bhz(im2,jp3)   &
!   &             + cm1x * cm3y * bhz(im1,jm3)   &
!   &             + cm1x * cm2y * bhz(im1,jm2)   &
!   &             + cm1x * cm1y * bhz(im1,jm1)   &
!   &             + cm1x * cy   * bhz(im1,j  )   &
!   &             + cm1x * cp1y * bhz(im1,jp1)   &
!   &             + cm1x * cp2y * bhz(im1,jp2)   &
!   &             + cm1x * cp3y * bhz(im1,jp3)   &
!   &             + cx   * cm3y * bhz(i  ,jm3)   &
!   &             + cx   * cm2y * bhz(i  ,jm2)   &
!   &             + cx   * cm1y * bhz(i  ,jm1)   &
!   &             + cx   * cy   * bhz(i  ,j  )   &
!   &             + cx   * cp1y * bhz(i  ,jp1)   &
!   &             + cx   * cp2y * bhz(i  ,jp2)   &
!   &             + cx   * cp3y * bhz(i  ,jp3)   &
!   &             + cp1x * cm3y * bhz(ip1,jm3)   &
!   &             + cp1x * cm2y * bhz(ip1,jm2)   &
!   &             + cp1x * cm1y * bhz(ip1,jm1)   &
!   &             + cp1x * cy   * bhz(ip1,j  )   &
!   &             + cp1x * cp1y * bhz(ip1,jp1)   &
!   &             + cp1x * cp2y * bhz(ip1,jp2)   &
!   &             + cp1x * cp3y * bhz(ip1,jp3)   &
!   &             + cp2x * cm3y * bhz(ip2,jm3)   &
!   &             + cp2x * cm2y * bhz(ip2,jm2)   &
!   &             + cp2x * cm1y * bhz(ip2,jm1)   &
!   &             + cp2x * cy   * bhz(ip2,j  )   &
!   &             + cp2x * cp1y * bhz(ip2,jp1)   &
!   &             + cp2x * cp2y * bhz(ip2,jp2)   &
!   &             + cp2x * cp3y * bhz(ip2,jp3)   &
!   &             + cp3x * cm3y * bhz(ip3,jm3)   &
!   &             + cp3x * cm2y * bhz(ip3,jm2)   &
!   &             + cp3x * cm1y * bhz(ip3,jm1)   &
!   &             + cp3x * cy   * bhz(ip3,j  )   &
!   &             + cp3x * cp1y * bhz(ip3,jp1)   &
!   &             + cp3x * cp2y * bhz(ip3,jp2)   &
!   &             + cp3x * cp3y * bhz(ip3,jp3) 

end do

deallocate(ehx)
deallocate(ehy)
!deallocate(bhz)

end subroutine interpol_eb_m6

subroutine calcul_energy( ele, tm, energy)

type(particle)       :: ele
type(tm_mesh_fields) :: tm

sll_int32  :: i, j, k
sll_int32  :: im1, im2, im3, ip1, ip2, ip3
sll_int32  :: jm1, jm2, jm3, jp1, jp2, jp3
sll_real64 :: dpx, dpy
sll_real64 :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
sll_real64 :: cy, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y
sll_real64 :: weight, rho_total,velocity
real(8), intent(out)    :: energy
tm%r0 = 0.0_f64

do k = 1, nbpart

velocity=(ele%vpx(k))**2+(ele%vpy(k))**2!dcos(t)*ele%vpx(k)+dsin(t)*ele%vpy(k)!

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

tm%r0(im3,jm3) = tm%r0(im3,jm3) + cm3x * cm3y * weight*velocity
tm%r0(im3,jm2) = tm%r0(im3,jm2) + cm3x * cm2y * weight*velocity
tm%r0(im3,jm1) = tm%r0(im3,jm1) + cm3x * cm1y * weight*velocity
tm%r0(im3,j  ) = tm%r0(im3,j  ) + cm3x * cy   * weight*velocity
tm%r0(im3,jp1) = tm%r0(im3,jp1) + cm3x * cp1y * weight*velocity
tm%r0(im3,jp2) = tm%r0(im3,jp2) + cm3x * cp2y * weight*velocity
tm%r0(im3,jp3) = tm%r0(im3,jp3) + cm3x * cp3y * weight*velocity

tm%r0(im2,jm3) = tm%r0(im2,jm3) + cm2x * cm3y * weight*velocity
tm%r0(im2,jm2) = tm%r0(im2,jm2) + cm2x * cm2y * weight*velocity
tm%r0(im2,jm1) = tm%r0(im2,jm1) + cm2x * cm1y * weight*velocity
tm%r0(im2,j  ) = tm%r0(im2,j  ) + cm2x * cy   * weight*velocity
tm%r0(im2,jp1) = tm%r0(im2,jp1) + cm2x * cp1y * weight*velocity
tm%r0(im2,jp2) = tm%r0(im2,jp2) + cm2x * cp2y * weight*velocity
tm%r0(im2,jp3) = tm%r0(im2,jp3) + cm2x * cp3y * weight*velocity

tm%r0(im1,jm3) = tm%r0(im1,jm3) + cm1x * cm3y * weight*velocity
tm%r0(im1,jm2) = tm%r0(im1,jm2) + cm1x * cm2y * weight*velocity
tm%r0(im1,jm1) = tm%r0(im1,jm1) + cm1x * cm1y * weight*velocity
tm%r0(im1,j  ) = tm%r0(im1,j  ) + cm1x * cy   * weight*velocity
tm%r0(im1,jp1) = tm%r0(im1,jp1) + cm1x * cp1y * weight*velocity
tm%r0(im1,jp2) = tm%r0(im1,jp2) + cm1x * cp2y * weight*velocity
tm%r0(im1,jp3) = tm%r0(im1,jp3) + cm1x * cp3y * weight*velocity

tm%r0(i  ,jm3) = tm%r0(i  ,jm3) + cx   * cm3y * weight*velocity
tm%r0(i  ,jm2) = tm%r0(i  ,jm2) + cx   * cm2y * weight*velocity
tm%r0(i  ,jm1) = tm%r0(i  ,jm1) + cx   * cm1y * weight*velocity
tm%r0(i  ,j  ) = tm%r0(i  ,j  ) + cx   * cy   * weight*velocity
tm%r0(i  ,jp1) = tm%r0(i  ,jp1) + cx   * cp1y * weight*velocity
tm%r0(i  ,jp2) = tm%r0(i  ,jp2) + cx   * cp2y * weight*velocity
tm%r0(i  ,jp3) = tm%r0(i  ,jp3) + cx   * cp3y * weight*velocity

tm%r0(ip1,jm3) = tm%r0(ip1,jm3) + cp1x * cm3y * weight*velocity
tm%r0(ip1,jm2) = tm%r0(ip1,jm2) + cp1x * cm2y * weight*velocity
tm%r0(ip1,jm1) = tm%r0(ip1,jm1) + cp1x * cm1y * weight*velocity
tm%r0(ip1,j  ) = tm%r0(ip1,j  ) + cp1x * cy   * weight*velocity
tm%r0(ip1,jp1) = tm%r0(ip1,jp1) + cp1x * cp1y * weight*velocity
tm%r0(ip1,jp2) = tm%r0(ip1,jp2) + cp1x * cp2y * weight*velocity
tm%r0(ip1,jp3) = tm%r0(ip1,jp3) + cp1x * cp3y * weight*velocity

tm%r0(ip2,jm3) = tm%r0(ip2,jm3) + cp2x * cm3y * weight*velocity
tm%r0(ip2,jm2) = tm%r0(ip2,jm2) + cp2x * cm2y * weight*velocity
tm%r0(ip2,jm1) = tm%r0(ip2,jm1) + cp2x * cm1y * weight*velocity
tm%r0(ip2,j  ) = tm%r0(ip2,j  ) + cp2x * cy   * weight*velocity
tm%r0(ip2,jp1) = tm%r0(ip2,jp1) + cp2x * cp1y * weight*velocity
tm%r0(ip2,jp2) = tm%r0(ip2,jp2) + cp2x * cp2y * weight*velocity
tm%r0(ip2,jp3) = tm%r0(ip2,jp3) + cp2x * cp3y * weight*velocity

tm%r0(ip3,jm3) = tm%r0(ip3,jm3) + cp3x * cm3y * weight*velocity
tm%r0(ip3,jm2) = tm%r0(ip3,jm2) + cp3x * cm2y * weight*velocity
tm%r0(ip3,jm1) = tm%r0(ip3,jm1) + cp3x * cm1y * weight*velocity
tm%r0(ip3,j  ) = tm%r0(ip3,j  ) + cp3x * cy   * weight*velocity
tm%r0(ip3,jp1) = tm%r0(ip3,jp1) + cp3x * cp1y * weight*velocity
tm%r0(ip3,jp2) = tm%r0(ip3,jp2) + cp3x * cp2y * weight*velocity
tm%r0(ip3,jp3) = tm%r0(ip3,jp3) + cp3x * cp3y * weight*velocity

end do

tm%r0(0:nx-1,ny)  = tm%r0(0:nx-1,0)
tm%r0(nx,0:ny-1)  = tm%r0(0,0:ny-1)
tm%r0(nx,ny)      = tm%r0(0,0)

!tm%ex= tm%r0 / (dx*dy)

energy = sum(tm%r0(0:nx,0:ny))/2.0d0
do i=0,nx
do j=0,ny
energy=energy+(tm%ex(i,j)**2+tm%ey(i,j)**2)*dx*dy/2.0d0
enddo
enddo
!tm%r0 = 0.0_f64
!
!do k = 1, nbpart
!
!velocity=dcos(t)*ele%vpy(k)-dsin(t)*ele%vpx(k)
!
!i      = ele%idx(k)
!j      = ele%idy(k)
!dpx    = ele%dpx(k)
!dpy    = ele%dpy(k)
!weight = ele%p(k)
!
!im3 = modulo(i-3,nx)
!im2 = modulo(i-2,nx)
!im1 = modulo(i-1,nx)
!ip1 = modulo(i+1,nx)
!ip2 = modulo(i+2,nx)
!ip3 = modulo(i+3,nx)
!jm3 = modulo(j-3,ny)
!jm2 = modulo(j-2,ny)
!jm1 = modulo(j-1,ny)
!jp1 = modulo(j+1,ny)
!jp2 = modulo(j+2,ny)
!jp3 = modulo(j+3,ny)
!
!cm3x = f_m6(3.0d0+dpx)
!cp3x = f_m6(3.0d0-dpx)
!cm2x = f_m6(2.0d0+dpx)
!cp2x = f_m6(2.0d0-dpx)
!cm1x = f_m6(1.0d0+dpx)
!cp1x = f_m6(1.0d0-dpx)
!cx   = f_m6(dpx)
!cy   = f_m6(dpy)
!cp1y = f_m6(1.0d0-dpy)
!cm1y = f_m6(1.0d0+dpy)
!cp2y = f_m6(2.0d0-dpy)
!cm2y = f_m6(2.0d0+dpy)
!cp3y = f_m6(3.0d0-dpy)
!cm3y = f_m6(3.0d0+dpy)
!
!tm%r0(im3,jm3) = tm%r0(im3,jm3) + cm3x * cm3y * weight*velocity
!tm%r0(im3,jm2) = tm%r0(im3,jm2) + cm3x * cm2y * weight*velocity
!tm%r0(im3,jm1) = tm%r0(im3,jm1) + cm3x * cm1y * weight*velocity
!tm%r0(im3,j  ) = tm%r0(im3,j  ) + cm3x * cy   * weight*velocity
!tm%r0(im3,jp1) = tm%r0(im3,jp1) + cm3x * cp1y * weight*velocity
!tm%r0(im3,jp2) = tm%r0(im3,jp2) + cm3x * cp2y * weight*velocity
!tm%r0(im3,jp3) = tm%r0(im3,jp3) + cm3x * cp3y * weight*velocity
!
!tm%r0(im2,jm3) = tm%r0(im2,jm3) + cm2x * cm3y * weight*velocity
!tm%r0(im2,jm2) = tm%r0(im2,jm2) + cm2x * cm2y * weight*velocity
!tm%r0(im2,jm1) = tm%r0(im2,jm1) + cm2x * cm1y * weight*velocity
!tm%r0(im2,j  ) = tm%r0(im2,j  ) + cm2x * cy   * weight*velocity
!tm%r0(im2,jp1) = tm%r0(im2,jp1) + cm2x * cp1y * weight*velocity
!tm%r0(im2,jp2) = tm%r0(im2,jp2) + cm2x * cp2y * weight*velocity
!tm%r0(im2,jp3) = tm%r0(im2,jp3) + cm2x * cp3y * weight*velocity
!
!tm%r0(im1,jm3) = tm%r0(im1,jm3) + cm1x * cm3y * weight*velocity
!tm%r0(im1,jm2) = tm%r0(im1,jm2) + cm1x * cm2y * weight*velocity
!tm%r0(im1,jm1) = tm%r0(im1,jm1) + cm1x * cm1y * weight*velocity
!tm%r0(im1,j  ) = tm%r0(im1,j  ) + cm1x * cy   * weight*velocity
!tm%r0(im1,jp1) = tm%r0(im1,jp1) + cm1x * cp1y * weight*velocity
!tm%r0(im1,jp2) = tm%r0(im1,jp2) + cm1x * cp2y * weight*velocity
!tm%r0(im1,jp3) = tm%r0(im1,jp3) + cm1x * cp3y * weight*velocity
!
!tm%r0(i  ,jm3) = tm%r0(i  ,jm3) + cx   * cm3y * weight*velocity
!tm%r0(i  ,jm2) = tm%r0(i  ,jm2) + cx   * cm2y * weight*velocity
!tm%r0(i  ,jm1) = tm%r0(i  ,jm1) + cx   * cm1y * weight*velocity
!tm%r0(i  ,j  ) = tm%r0(i  ,j  ) + cx   * cy   * weight*velocity
!tm%r0(i  ,jp1) = tm%r0(i  ,jp1) + cx   * cp1y * weight*velocity
!tm%r0(i  ,jp2) = tm%r0(i  ,jp2) + cx   * cp2y * weight*velocity
!tm%r0(i  ,jp3) = tm%r0(i  ,jp3) + cx   * cp3y * weight*velocity
!
!tm%r0(ip1,jm3) = tm%r0(ip1,jm3) + cp1x * cm3y * weight*velocity
!tm%r0(ip1,jm2) = tm%r0(ip1,jm2) + cp1x * cm2y * weight*velocity
!tm%r0(ip1,jm1) = tm%r0(ip1,jm1) + cp1x * cm1y * weight*velocity
!tm%r0(ip1,j  ) = tm%r0(ip1,j  ) + cp1x * cy   * weight*velocity
!tm%r0(ip1,jp1) = tm%r0(ip1,jp1) + cp1x * cp1y * weight*velocity
!tm%r0(ip1,jp2) = tm%r0(ip1,jp2) + cp1x * cp2y * weight*velocity
!tm%r0(ip1,jp3) = tm%r0(ip1,jp3) + cp1x * cp3y * weight*velocity
!
!tm%r0(ip2,jm3) = tm%r0(ip2,jm3) + cp2x * cm3y * weight*velocity
!tm%r0(ip2,jm2) = tm%r0(ip2,jm2) + cp2x * cm2y * weight*velocity
!tm%r0(ip2,jm1) = tm%r0(ip2,jm1) + cp2x * cm1y * weight*velocity
!tm%r0(ip2,j  ) = tm%r0(ip2,j  ) + cp2x * cy   * weight*velocity
!tm%r0(ip2,jp1) = tm%r0(ip2,jp1) + cp2x * cp1y * weight*velocity
!tm%r0(ip2,jp2) = tm%r0(ip2,jp2) + cp2x * cp2y * weight*velocity
!tm%r0(ip2,jp3) = tm%r0(ip2,jp3) + cp2x * cp3y * weight*velocity
!
!tm%r0(ip3,jm3) = tm%r0(ip3,jm3) + cp3x * cm3y * weight*velocity
!tm%r0(ip3,jm2) = tm%r0(ip3,jm2) + cp3x * cm2y * weight*velocity
!tm%r0(ip3,jm1) = tm%r0(ip3,jm1) + cp3x * cm1y * weight*velocity
!tm%r0(ip3,j  ) = tm%r0(ip3,j  ) + cp3x * cy   * weight*velocity
!tm%r0(ip3,jp1) = tm%r0(ip3,jp1) + cp3x * cp1y * weight*velocity
!tm%r0(ip3,jp2) = tm%r0(ip3,jp2) + cp3x * cp2y * weight*velocity
!tm%r0(ip3,jp3) = tm%r0(ip3,jp3) + cp3x * cp3y * weight*velocity
!
!end do
!
!tm%r0(0:nx-1,ny)  = tm%r0(0:nx-1,0)
!tm%r0(nx,0:ny-1)  = tm%r0(0,0:ny-1)
!tm%r0(nx,ny)      = tm%r0(0,0)
!
!tm%ey= tm%r0 / (dx*dy)

end subroutine calcul_energy


subroutine calcul_momentum( ele, momentum)

type(particle)       :: ele

sll_int32  :: i, j, k,nv
sll_int32  :: im1, im2, im3, ip1, ip2, ip3
sll_int32  :: jm1, jm2, jm3, jp1, jp2, jp3
sll_real64 :: dpx, dpy, hv
sll_real64 :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
sll_real64 :: cy, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y
sll_real64 :: weight
real(8), intent(inout)   :: momentum(0:64,0:64)

nv=64
hv=8.0d0/64.0d0
momentum=0.0d0

do k = 1, nbpart
    do while ( ele%vpx(k) > 4.0d0 )
    ele%vpx(k) = ele%vpx(k) - 4.0d0
    enddo
    do while ( ele%vpx(k) < -4.0d0 )
    ele%vpx(k) = ele%vpx(k) + 4.0d0
    enddo
    do while ( ele%vpy(k) > 4.0d0 )
    ele%vpy(k) = ele%vpy(k) - 4.0d0
    enddo
    do while ( ele%vpy(k) < -4.0d0 )
    ele%vpy(k) = ele%vpy(k) + 4.0d0
    enddo
    ele%idx(k) = floor((ele%vpx(k)+4.0d0)/8.0d0*nv)
    ele%vpx(k) = real((ele%vpx(k)+4.0d0)/hv- ele%idx(k), f64)
    ele%idy(k) = floor((ele%vpy(k)+4.0d0)/8.0d0*nv)
    ele%vpy(k) = real((ele%vpy(k)+4.0d0)/hv- ele%idy(k), f64)
enddo

do k = 1, nbpart

i      = ele%idx(k)
j      = ele%idy(k)
dpx    = ele%vpx(k)
dpy    = ele%vpy(k)
weight = ele%p(k)

im3 = modulo(i-3,nv)
im2 = modulo(i-2,nv)
im1 = modulo(i-1,nv)
ip1 = modulo(i+1,nv)
ip2 = modulo(i+2,nv)
ip3 = modulo(i+3,nv)
jm3 = modulo(j-3,nv)
jm2 = modulo(j-2,nv)
jm1 = modulo(j-1,nv)
jp1 = modulo(j+1,nv)
jp2 = modulo(j+2,nv)
jp3 = modulo(j+3,nv)

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

momentum(im3,jm3) = momentum(im3,jm3) + cm3x * cm3y * weight
momentum(im3,jm2) = momentum(im3,jm2) + cm3x * cm2y * weight
momentum(im3,jm1) = momentum(im3,jm1) + cm3x * cm1y * weight
momentum(im3,j  ) = momentum(im3,j  ) + cm3x * cy   * weight
momentum(im3,jp1) = momentum(im3,jp1) + cm3x * cp1y * weight
momentum(im3,jp2) = momentum(im3,jp2) + cm3x * cp2y * weight
momentum(im3,jp3) = momentum(im3,jp3) + cm3x * cp3y * weight

momentum(im2,jm3) = momentum(im2,jm3) + cm2x * cm3y * weight
momentum(im2,jm2) = momentum(im2,jm2) + cm2x * cm2y * weight
momentum(im2,jm1) = momentum(im2,jm1) + cm2x * cm1y * weight
momentum(im2,j  ) = momentum(im2,j  ) + cm2x * cy   * weight
momentum(im2,jp1) = momentum(im2,jp1) + cm2x * cp1y * weight
momentum(im2,jp2) = momentum(im2,jp2) + cm2x * cp2y * weight
momentum(im2,jp3) = momentum(im2,jp3) + cm2x * cp3y * weight

momentum(im1,jm3) = momentum(im1,jm3) + cm1x * cm3y * weight
momentum(im1,jm2) = momentum(im1,jm2) + cm1x * cm2y * weight
momentum(im1,jm1) = momentum(im1,jm1) + cm1x * cm1y * weight
momentum(im1,j  ) = momentum(im1,j  ) + cm1x * cy   * weight
momentum(im1,jp1) = momentum(im1,jp1) + cm1x * cp1y * weight
momentum(im1,jp2) = momentum(im1,jp2) + cm1x * cp2y * weight
momentum(im1,jp3) = momentum(im1,jp3) + cm1x * cp3y * weight

momentum(i  ,jm3) = momentum(i  ,jm3) + cx   * cm3y * weight
momentum(i  ,jm2) = momentum(i  ,jm2) + cx   * cm2y * weight
momentum(i  ,jm1) = momentum(i  ,jm1) + cx   * cm1y * weight
momentum(i  ,j  ) = momentum(i  ,j  ) + cx   * cy   * weight
momentum(i  ,jp1) = momentum(i  ,jp1) + cx   * cp1y * weight
momentum(i  ,jp2) = momentum(i  ,jp2) + cx   * cp2y * weight
momentum(i  ,jp3) = momentum(i  ,jp3) + cx   * cp3y * weight

momentum(ip1,jm3) = momentum(ip1,jm3) + cp1x * cm3y * weight
momentum(ip1,jm2) = momentum(ip1,jm2) + cp1x * cm2y * weight
momentum(ip1,jm1) = momentum(ip1,jm1) + cp1x * cm1y * weight
momentum(ip1,j  ) = momentum(ip1,j  ) + cp1x * cy   * weight
momentum(ip1,jp1) = momentum(ip1,jp1) + cp1x * cp1y * weight
momentum(ip1,jp2) = momentum(ip1,jp2) + cp1x * cp2y * weight
momentum(ip1,jp3) = momentum(ip1,jp3) + cp1x * cp3y * weight

momentum(ip2,jm3) = momentum(ip2,jm3) + cp2x * cm3y * weight
momentum(ip2,jm2) = momentum(ip2,jm2) + cp2x * cm2y * weight
momentum(ip2,jm1) = momentum(ip2,jm1) + cp2x * cm1y * weight
momentum(ip2,j  ) = momentum(ip2,j  ) + cp2x * cy   * weight
momentum(ip2,jp1) = momentum(ip2,jp1) + cp2x * cp1y * weight
momentum(ip2,jp2) = momentum(ip2,jp2) + cp2x * cp2y * weight
momentum(ip2,jp3) = momentum(ip2,jp3) + cp2x * cp3y * weight

momentum(ip3,jm3) = momentum(ip3,jm3) + cp3x * cm3y * weight
momentum(ip3,jm2) = momentum(ip3,jm2) + cp3x * cm2y * weight
momentum(ip3,jm1) = momentum(ip3,jm1) + cp3x * cm1y * weight
momentum(ip3,j  ) = momentum(ip3,j  ) + cp3x * cy   * weight
momentum(ip3,jp1) = momentum(ip3,jp1) + cp3x * cp1y * weight
momentum(ip3,jp2) = momentum(ip3,jp2) + cp3x * cp2y * weight
momentum(ip3,jp3) = momentum(ip3,jp3) + cp3x * cp3y * weight

end do

momentum(0:nv-1,nv)  = momentum(0:nv-1,0)
momentum(nv,0:nv-1)  = momentum(0,0:nv-1)
momentum(nv,nv)      = momentum(0,0)

momentum = momentum / (hv**2)

end subroutine calcul_momentum

end module particules_m6
