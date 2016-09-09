module particules_m6
#include "sll_working_precision.h"
#include "sll_memory.h"

use particules

implicit none

private

public calcul_rho_m6
public interpol_eb_m6

sll_real64, allocatable :: ehx(:,:)
sll_real64, allocatable :: ehy(:,:)

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

f_m6 = 0.0d0
if ( q < 3.0d0 ) f_m6 = (3.0d0-q)**5
if ( q < 2.0d0 ) f_m6 = f_m6 - 6.0d0*(2.0d0-q)**5
if ( q < 1.0d0 ) f_m6 = f_m6 + 15.0d0*(1.0d0-q)**5

f_m6 = f_m6 * 0.008333333333333333_f64 !/ 120.0d0

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

if (.not. allocated(ehx)) allocate(ehx(0:nx-1,0:ny-1))
if (.not. allocated(ehy)) allocate(ehy(0:nx-1,0:ny-1))

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

  end do
end do

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
   &             + cm3x_cm3y * ehx(im3,jm3)   &
   &             + cm3x_cm2y * ehx(im3,jm2)   &
   &             + cm3x_cm1y * ehx(im3,jm1)   &
   &             + cm3x_cy   * ehx(im3,j  )   &
   &             + cm3x_cp1y * ehx(im3,jp1)   &
   &             + cm3x_cp2y * ehx(im3,jp2)   &
   &             + cm3x_cp3y * ehx(im3,jp3)   &
   &             + cm2x_cm3y * ehx(im2,jm3)   &
   &             + cm2x_cm2y * ehx(im2,jm2)   &
   &             + cm2x_cm1y * ehx(im2,jm1)   &
   &             + cm2x_cy   * ehx(im2,j  )   &
   &             + cm2x_cp1y * ehx(im2,jp1)   &
   &             + cm2x_cp2y * ehx(im2,jp2)   &
   &             + cm2x_cp3y * ehx(im2,jp3)   &
   &             + cm1x_cm3y * ehx(im1,jm3)   &
   &             + cm1x_cm2y * ehx(im1,jm2)   &
   &             + cm1x_cm1y * ehx(im1,jm1)   &
   &             + cm1x_cy   * ehx(im1,j  )   &
   &             + cm1x_cp1y * ehx(im1,jp1)   &
   &             + cm1x_cp2y * ehx(im1,jp2)   &
   &             + cm1x_cp3y * ehx(im1,jp3)   &
   &             + cx_cm3y   * ehx(i  ,jm3)   &
   &             + cx_cm2y   * ehx(i  ,jm2)   &
   &             + cx_cm1y   * ehx(i  ,jm1)   &
   &             + cx_cy     * ehx(i  ,j  )   &
   &             + cx_cp1y   * ehx(i  ,jp1)   &
   &             + cx_cp2y   * ehx(i  ,jp2)   &
   &             + cx_cp3y   * ehx(i  ,jp3)   &
   &             + cp1x_cm3y * ehx(ip1,jm3)   &
   &             + cp1x_cm2y * ehx(ip1,jm2)   &
   &             + cp1x_cm1y * ehx(ip1,jm1)   &
   &             + cp1x_cy   * ehx(ip1,j  )   &
   &             + cp1x_cp1y * ehx(ip1,jp1)   &
   &             + cp1x_cp2y * ehx(ip1,jp2)   &
   &             + cp1x_cp3y * ehx(ip1,jp3)   &
   &             + cp2x_cm3y * ehx(ip2,jm3)   &
   &             + cp2x_cm2y * ehx(ip2,jm2)   &
   &             + cp2x_cm1y * ehx(ip2,jm1)   &
   &             + cp2x_cy   * ehx(ip2,j  )   &
   &             + cp2x_cp1y * ehx(ip2,jp1)   &
   &             + cp2x_cp2y * ehx(ip2,jp2)   &
   &             + cp2x_cp3y * ehx(ip2,jp3)   &
   &             + cp3x_cm3y * ehx(ip3,jm3)   &
   &             + cp3x_cm2y * ehx(ip3,jm2)   &
   &             + cp3x_cm1y * ehx(ip3,jm1)   &
   &             + cp3x_cy   * ehx(ip3,j  )   &
   &             + cp3x_cp1y * ehx(ip3,jp1)   &
   &             + cp3x_cp2y * ehx(ip3,jp2)   &
   &             + cp3x_cp3y * ehx(ip3,jp3) 

   ele%epy(k) =                               &
   &             + cm3x_cm3y * ehy(im3,jm3)   &
   &             + cm3x_cm2y * ehy(im3,jm2)   &
   &             + cm3x_cm1y * ehy(im3,jm1)   &
   &             + cm3x_cy   * ehy(im3,j  )   &
   &             + cm3x_cp1y * ehy(im3,jp1)   &
   &             + cm3x_cp2y * ehy(im3,jp2)   &
   &             + cm3x_cp3y * ehy(im3,jp3)   &
   &             + cm2x_cm3y * ehy(im2,jm3)   &
   &             + cm2x_cm2y * ehy(im2,jm2)   &
   &             + cm2x_cm1y * ehy(im2,jm1)   &
   &             + cm2x_cy   * ehy(im2,j  )   &
   &             + cm2x_cp1y * ehy(im2,jp1)   &
   &             + cm2x_cp2y * ehy(im2,jp2)   &
   &             + cm2x_cp3y * ehy(im2,jp3)   &
   &             + cm1x_cm3y * ehy(im1,jm3)   &
   &             + cm1x_cm2y * ehy(im1,jm2)   &
   &             + cm1x_cm1y * ehy(im1,jm1)   &
   &             + cm1x_cy   * ehy(im1,j  )   &
   &             + cm1x_cp1y * ehy(im1,jp1)   &
   &             + cm1x_cp2y * ehy(im1,jp2)   &
   &             + cm1x_cp3y * ehy(im1,jp3)   &
   &             + cx_cm3y   * ehy(i  ,jm3)   &
   &             + cx_cm2y   * ehy(i  ,jm2)   &
   &             + cx_cm1y   * ehy(i  ,jm1)   &
   &             + cx_cy     * ehy(i  ,j  )   &
   &             + cx_cp1y   * ehy(i  ,jp1)   &
   &             + cx_cp2y   * ehy(i  ,jp2)   &
   &             + cx_cp3y   * ehy(i  ,jp3)   &
   &             + cp1x_cm3y * ehy(ip1,jm3)   &
   &             + cp1x_cm2y * ehy(ip1,jm2)   &
   &             + cp1x_cm1y * ehy(ip1,jm1)   &
   &             + cp1x_cy   * ehy(ip1,j  )   &
   &             + cp1x_cp1y * ehy(ip1,jp1)   &
   &             + cp1x_cp2y * ehy(ip1,jp2)   &
   &             + cp1x_cp3y * ehy(ip1,jp3)   &
   &             + cp2x_cm3y * ehy(ip2,jm3)   &
   &             + cp2x_cm2y * ehy(ip2,jm2)   &
   &             + cp2x_cm1y * ehy(ip2,jm1)   &
   &             + cp2x_cy   * ehy(ip2,j  )   &
   &             + cp2x_cp1y * ehy(ip2,jp1)   &
   &             + cp2x_cp2y * ehy(ip2,jp2)   &
   &             + cp2x_cp3y * ehy(ip2,jp3)   &
   &             + cp3x_cm3y * ehy(ip3,jm3)   &
   &             + cp3x_cm2y * ehy(ip3,jm2)   &
   &             + cp3x_cm1y * ehy(ip3,jm1)   &
   &             + cp3x_cy   * ehy(ip3,j  )   &
   &             + cp3x_cp1y * ehy(ip3,jp1)   &
   &             + cp3x_cp2y * ehy(ip3,jp2)   &
   &             + cp3x_cp3y * ehy(ip3,jp3) 

end do

end subroutine interpol_eb_m6

end module particules_m6
