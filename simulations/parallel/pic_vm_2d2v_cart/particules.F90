module particules
#include "sll_working_precision.h"
#include "sll_memory.h"
use zone
use quietstart

implicit none

type particle
  sll_real32, pointer :: dpx(:)
  sll_real32, pointer :: dpy(:)
  sll_int32 , pointer :: idx(:)
  sll_int32 , pointer :: idy(:)
  sll_real64, pointer :: vpx(:)
  sll_real64, pointer :: vpy(:)
  sll_real64, pointer :: epx(:)
  sll_real64, pointer :: epy(:)
  sll_real64, pointer :: bpz(:)
  sll_real64, pointer :: p(:)
end type particle

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_eb( tm1, ele )

type  (particle) :: ele
type(tm_mesh_fields) :: tm1
sll_real64 :: a1, a2, a3, a4
sll_int32  :: k 
sll_int32  :: i, j
sll_real32 :: dpx
sll_real32 :: dpy
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|


do k=1,nbpart

   i   = ele%idx(k)
   j   = ele%idy(k)
   dpx = ele%dpx(k)
   dpy = ele%dpy(k)

   a1 = (1.0-dpx) * (1.0-dpy) 
   a2 = (    dpx) * (1.0-dpy) 
   a3 = (    dpx) * (    dpy) 
   a4 = (1.0-dpx) * (    dpy) 

   ele%epx(k) = a1 * tm1%ex(i  ,j  ) + a2 * tm1%ex(i+1,j  ) &
            & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i  ,j+1) 
   ele%epy(k) = a1 * tm1%ey(i  ,j  ) + a2 * tm1%ey(i+1,j  ) &
            & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i  ,j+1) 
   ele%bpz(k) = a1 * tm1%bz(i  ,j  ) + a2 * tm1%bz(i+1,j  ) &
            & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i  ,j+1) 

end do

end subroutine interpol_eb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_vitesse( ele )

type (particle) :: ele
sll_real64 :: dum, u2
sll_real64 :: tantheta, sintheta
sll_real64 :: gamma
sll_int32  :: k 

do k = 1, nbpart

   !*** Changement de variable u = gamma*vit

   if( relativ ) then

      u2  = ele%vpx(k)*ele%vpx(k) &
          + ele%vpy(k)*ele%vpy(k)
      if ( u2 >= csq ) then 
         print*,'Erreur : u2 >= c2 dans le calcul de la vitesse'
         print*,'k = ',k,' vx = ',ele%vpx(k),' vy = ',ele%vpy(k)
         stop
      else
         gamma = 1./sqrt( 1. - u2/csq )
      endif

      ele%vpx(k) = gamma*ele%vpx(k)
      ele%vpy(k) = gamma*ele%vpy(k)

   else

      gamma=1.0_f64

   end if


   !*** Separation des effets electriques et magnetiques

   !*** On ajoute la moitie de l'effet champ electrique E

   dum = 0.5 * dt * q_sur_m
   ele%vpx(k) = ele%vpx(k) + dum*(ele%epx(k)+exext)
   ele%vpy(k) = ele%vpy(k) + dum*(ele%epy(k)+eyext)

   !*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (ele%bpz(k)+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vpx(k) = ele%vpx(k) + ele%vpy(k)*tantheta
   ele%vpy(k) = ele%vpy(k) - ele%vpx(k)*sintheta
   ele%vpx(k) = ele%vpx(k) + ele%vpy(k)*tantheta

   !*** Autre moitie de l'effet du champ electrique E

   ele%vpx(k) = ele%vpx(k) + dum*(ele%epx(k)+exext)
   ele%vpy(k) = ele%vpy(k) + dum*(ele%epy(k)+eyext)

   !*** On repasse a la vitesse (changement de variable inverse)

   if( relativ ) then

      u2 =   ele%vpx(k)*ele%vpx(k) &
           + ele%vpy(k)*ele%vpy(k)

      gamma = sqrt( 1. + u2/csq )
 
      ele%vpx(k) = ele%vpx(k) / gamma
      ele%vpy(k) = ele%vpy(k) / gamma

   end if

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
sll_real64 :: coef
sll_int32  :: k 
sll_real64 :: ppx
sll_real64 :: ppy

do k=1,nbpart     

  ppx = (ele%idx(k)+ele%dpx(k))*dx + ele%vpx(k)*dt*coef
  ppy = (ele%idy(k)+ele%dpy(k))*dy + ele%vpy(k)*dt*coef

  if (bcname == 'period') then

    if( ppx >= dimx ) then
      ppx = ppx - dimx
    else if( ppx < 0.d0 ) then
      ppx = ppx + dimx
    end if

    if( ppy >= dimy ) then
      ppy = ppy - dimy
    else if( ppy < 0.d0 ) then
      ppy = ppy + dimy
    end if

  end if

  ele%idx(k) = floor(ppx/dimx*nx)
  ele%idy(k) = floor(ppy/dimy*ny)
  ele%dpx(k) = real(ppx/dx - ele%idx(k), f32)
  ele%dpy(k) = real(ppy/dy - ele%idy(k), f32) 

end do

if (bcname /= 'period') then
  k = 1
  do while ( k <= nbpart )  
    if (      ele%idx(k) < 0 .or. ele%idx(k) >= nx      &	
         .or. ele%idy(k) < 0 .or. ele%idy(k) >= ny ) then
        
        ele%dpx(k) = ele%dpx(nbpart)
        ele%dpy(k) = ele%dpy(nbpart)
        ele%vpx(k) = ele%vpx(nbpart)
        ele%vpy(k) = ele%vpy(nbpart)
        ele%p(k)   = ele%p(nbpart)
        ele%idx(k) = ele%idx(nbpart)
        ele%idy(k) = ele%idy(nbpart)
        nbpart = nbpart - 1
    if (nbpart == 0) stop 'no more particle'
    else
      k = k + 1
    end if
  end do

end if

end subroutine avancee_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_rho( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
sll_real64 :: a1, a2, a3, a4, weight
sll_real64 :: rho_total
sll_int32  :: k 
sll_int32  :: i, j 
sll_real32 :: dpx
sll_real32 :: dpy

tm%r0 = 0.d0    
                
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do k=1,nbpart

  i      = ele%idx(k)
  j      = ele%idy(k)
  dpx    = ele%dpx(k)
  dpy    = ele%dpy(k)
  weight = ele%p(k) 
  a1  = (1.0-dpx) * (1.0-dpy) * weight
  a2  = (dpx)     * (1.0-dpy) * weight
  a3  = (dpx)     * (dpy)     * weight
  a4  = (1.0-dpx) * (dpy)     * weight
  tm%r0(i,j)     = tm%r0(i,j)     + a1 
  tm%r0(i+1,j)   = tm%r0(i+1,j)   + a2 
  tm%r0(i+1,j+1) = tm%r0(i+1,j+1) + a3
  tm%r0(i,j+1)   = tm%r0(i,j+1)   + a4

end do

tm%r0 = tm%r0 / (dx*dy)

if (bcname == 'period') then
  do i=0,nx
    tm%r0(i,0)  = tm%r0(i,0) + tm%r0(i,ny)
    tm%r0(i,ny) = tm%r0(i,0)
  end do
  do j=0,ny
    tm%r0(0,j)  = tm%r0(0,j) + tm%r0(nx,j)
    tm%r0(nx,j) = tm%r0(0,j)
  end do
end if

rho_total = sum(tm%r0(1:nx,1:ny))*dx*dy
print*,'rho total',rho_total
! Neutralisation du milieu
tm%r0 = tm%r0 - rho_total/dimx/dimy

end subroutine calcul_rho

subroutine plasma( ele )

type (particle) :: ele
sll_real64 :: speed, theta, vth, n
sll_real64 :: a, b, eps, R,temm,xi,yi,zi
sll_real64 :: ppx, ppy
sll_int32  :: k, error

eps = 1.d-12

vth =  1.0_f64
nbpart =150000 !100*(nx)*(ny)
n = 1.d0/nbpart

SLL_ALLOCATE(ele%dpx(nbpart),error)
SLL_ALLOCATE(ele%dpy(nbpart),error)
SLL_ALLOCATE(ele%idx(nbpart),error)
SLL_ALLOCATE(ele%idy(nbpart),error)
SLL_ALLOCATE(ele%vpx(nbpart),error)
SLL_ALLOCATE(ele%vpy(nbpart),error)
SLL_ALLOCATE(ele%epx(nbpart),error)
SLL_ALLOCATE(ele%epy(nbpart),error)
SLL_ALLOCATE(ele%bpz(nbpart),error)
SLL_ALLOCATE(ele%p(nbpart),error)

do k=0,nbpart-1

  speed = vth * sqrt(-2.0_f64 * log( (k+0.5)*n ))

  theta = trinary_reversing( k ) * 2.0_f64 * pi

!  a = 0.0_f64; b = dimx ! 2*pi/kx 
!  R = bit_reversing( k )
!  call dichotomie_x(a,b,R,eps) 
!  ppx = a
!  ppy = dimy * penta_reversing( k ) 

!  ele%idx(k+1) = floor(ppx/dimx*nx)
!  ele%idy(k+1) = floor(ppy/dimy*ny)
!  ele%dpx(k+1) = real(ppx/dx - ele%idx(k+1), f32)
!  ele%dpy(k+1) = real(ppy/dx - ele%idy(k+1), f32)
  ele%vpx(k+1) = speed * cos(theta)
  ele%vpy(k+1) = speed * sin(theta)

  ele%p(k+1) = poids * n

enddo

!--add by me: rejection sampling--
k=1
do while (k<=nbpart)
    call random_number(xi)
    xi=dimx*xi
    call random_number(yi)
    yi=dimy*yi
    call random_number(zi)
    zi=2.0d0+alpha
    temm=1.0d0+dsin(yi)+alpha*dcos(kx*xi)
    if (temm>=zi) then
        ele%idx(k) = floor(xi/dimx*nx)
        ele%idy(k) = floor(yi/dimy*ny)
        ele%dpx(k) = real(xi/dx - ele%idx(k), f32)
        ele%dpy(k) = real(yi/dx - ele%idy(k), f32)
        k=k+1
    endif
enddo

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> M4 function from Monhagan (SPH method)
!> Cubic spline
!>       |  1 - 1.5 q^2 + 0.75 q^3  for 0 <= q <= 1
!> M4(x) |  1/4 (2−q)^3,            for 1 <= q <= 2
!>       |  0                       for q > 2.

function f_m4( x )
sll_real32, intent(in) :: x
sll_real32, parameter  :: pi = 3.1415926535897932384626433832795
sll_real32             :: f_m4

if ( x < 1. ) then
   f_m4 = 1 - 1.5*x*x+0.75*x*x*x
else if ( x >= 1. .and. x <= 2. ) then
   f_m4 = 0.25*(2.-x)**3 
else
   f_m4 = 0.0
end if

f_m4 = f_m4 * 2.0 / 3.0

return
end function f_m4

subroutine interpol_eb_m4( tm1, ele )

type(particle)       :: ele
type(tm_mesh_fields) :: tm1
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
   cm2y = f_m4(2.+dpy)
   cp2y = f_m4(2.-dpy)
   cm1y = f_m4(1.+dpy)
   cp1y = f_m4(1.-dpy)

   ele%epx(k) =                                    &
   &             + cm2x * cm2y * tm1%ex(im2,jm2)   &
   &             + cm2x * cm1y * tm1%ex(im2,jm1)   &
   &             + cm2x * cy   * tm1%ex(im2,j  )   &
   &             + cm2x * cp1y * tm1%ex(im2,jp1)   &
   &             + cm2x * cp2y * tm1%ex(im2,jp2)   &
   &             + cm1x * cm2y * tm1%ex(im1,jm2)   &
   &             + cm1x * cm1y * tm1%ex(im1,jm1)   &
   &             + cm1x * cy   * tm1%ex(im1,j  )   &
   &             + cm1x * cp1y * tm1%ex(im1,j+1)   &
   &             + cm1x * cp2y * tm1%ex(im1,jp2)   &
   &             + cx   * cm2y * tm1%ex(i  ,jm2)   &
   &             + cx   * cm1y * tm1%ex(i  ,jm1)   &
   &             + cx   * cy   * tm1%ex(i  ,j  )   &
   &             + cx   * cp1y * tm1%ex(i  ,jp1)   &
   &             + cx   * cp2y * tm1%ex(i  ,jp2)   &
   &             + cp1x * cm2y * tm1%ex(ip1,jm2)   &
   &             + cp1x * cm1y * tm1%ex(ip1,jm1)   &
   &             + cp1x * cy   * tm1%ex(ip1,j  )   &
   &             + cp1x * cp1y * tm1%ex(ip1,jp1)   &
   &             + cp1x * cp2y * tm1%ex(ip1,jp2)   &
   &             + cp2x * cm2y * tm1%ex(ip2,jm2)   &
   &             + cp2x * cm1y * tm1%ex(ip2,jm1)   &
   &             + cp2x * cy   * tm1%ex(ip2,j  )   &
   &             + cp2x * cp1y * tm1%ex(ip2,jp1)   &
   &             + cp2x * cp2y * tm1%ex(ip2,jp2)   

   ele%epy(k) =                                    &
   &             + cm2x * cm2y * tm1%ey(im2,jm2)   &
   &             + cm2x * cm1y * tm1%ey(im2,jm1)   &
   &             + cm2x * cy   * tm1%ey(im2,j  )   &
   &             + cm2x * cp1y * tm1%ey(im2,jp1)   &
   &             + cm2x * cp2y * tm1%ey(im2,jp2)   &
   &             + cm1x * cm2y * tm1%ey(im1,jm2)   &
   &             + cm1x * cm1y * tm1%ey(im1,jm1)   &
   &             + cm1x * cy   * tm1%ey(im1,j  )   &
   &             + cm1x * cp1y * tm1%ey(im1,jp1)   &
   &             + cm1x * cp2y * tm1%ey(im1,jp2)   &
   &             + cx   * cm2y * tm1%ey(i  ,jm2)   &
   &             + cx   * cm1y * tm1%ey(i  ,jm1)   &
   &             + cx   * cy   * tm1%ey(i  ,j  )   &
   &             + cx   * cp1y * tm1%ey(i  ,jp1)   &
   &             + cx   * cp2y * tm1%ey(i  ,jp2)   &
   &             + cp1x * cm2y * tm1%ey(ip1,jm2)   &
   &             + cp1x * cm1y * tm1%ey(ip1,jm1)   &
   &             + cp1x * cy   * tm1%ey(ip1,j  )   &
   &             + cp1x * cp1y * tm1%ey(ip1,jp1)   &
   &             + cp1x * cp2y * tm1%ey(ip1,jp2)   &
   &             + cp2x * cm2y * tm1%ey(ip2,jm2)   &
   &             + cp2x * cm1y * tm1%ey(ip2,jm1)   &
   &             + cp2x * cy   * tm1%ey(ip2,j  )   &
   &             + cp2x * cp1y * tm1%ey(ip2,jp1)   &
   &             + cp2x * cp2y * tm1%ey(ip2,jp2)   

   ele%bpz(k) =                                    &
   &             + cm2x * cm2y * tm1%bz(im2,jm2)   &
   &             + cm2x * cm1y * tm1%bz(im2,jm1)   &
   &             + cm2x * cy   * tm1%bz(im2,j  )   &
   &             + cm2x * cp1y * tm1%bz(im2,jp1)   &
   &             + cm2x * cp2y * tm1%bz(im2,jp2)   &
   &             + cm1x * cm2y * tm1%bz(im1,jm2)   &
   &             + cm1x * cm1y * tm1%bz(im1,jm1)   &
   &             + cm1x * cy   * tm1%bz(im1,j  )   &
   &             + cm1x * cp1y * tm1%bz(im1,jp1)   &
   &             + cm1x * cp2y * tm1%bz(im1,jp2)   &
   &             + cx   * cm2y * tm1%bz(i  ,jm2)   &
   &             + cx   * cm1y * tm1%bz(i  ,jm1)   &
   &             + cx   * cy   * tm1%bz(i  ,j  )   &
   &             + cx   * cp1y * tm1%bz(i  ,jp1)   &
   &             + cx   * cp2y * tm1%bz(i  ,jp2)   &
   &             + cp1x * cm2y * tm1%bz(ip1,jm2)   &
   &             + cp1x * cm1y * tm1%bz(ip1,jm1)   &
   &             + cp1x * cy   * tm1%bz(ip1,j  )   &
   &             + cp1x * cp1y * tm1%bz(ip1,jp1)   &
   &             + cp1x * cp2y * tm1%bz(ip1,jp2)   &
   &             + cp2x * cm2y * tm1%bz(ip2,jm2)   &
   &             + cp2x * cm1y * tm1%bz(ip2,jm1)   &
   &             + cp2x * cy   * tm1%bz(ip2,j  )   &
   &             + cp2x * cp1y * tm1%bz(ip2,jp1)   &
   &             + cp2x * cp2y * tm1%bz(ip2,jp2)

end do 

end subroutine interpol_eb_m4


subroutine calcul_rho_m4( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm

sll_int32  :: i, j, k
sll_int32  :: im1, im2, ip1, ip2
sll_int32  :: jm1, jm2, jp1, jp2
sll_real32 :: dpx, dpy
sll_real32 :: cx, cm1x, cm2x, cp1x, cp2x
sll_real32 :: cy, cm1y, cm2y, cp1y, cp2y
sll_real64 :: weight, rho_total

tm%r0 = 0.0_f64

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

  tm%r0(im2,jm2) = tm%r0(im2,jm2) + cm2x * cm2y * weight
  tm%r0(im2,jm1) = tm%r0(im2,jm1) + cm2x * cm1y * weight
  tm%r0(im2,j  ) = tm%r0(im2,j  ) + cm2x * cy   * weight
  tm%r0(im2,jp1) = tm%r0(im2,jp1) + cm2x * cp1y * weight
  tm%r0(im2,jp2) = tm%r0(im2,jp2) + cm2x * cp2y * weight

  tm%r0(im1,jm2) = tm%r0(im1,jm2) + cm1x * cm2y * weight
  tm%r0(im1,jm1) = tm%r0(im1,jm1) + cm1x * cm1y * weight
  tm%r0(im1,j  ) = tm%r0(im1,j  ) + cm1x * cy   * weight
  tm%r0(im1,jp1) = tm%r0(im1,jp1) + cm1x * cp1y * weight
  tm%r0(im1,jp2) = tm%r0(im1,jp2) + cm1x * cp2y * weight

  tm%r0(i  ,jm2) = tm%r0(i  ,jm2) + cx   * cm2y * weight
  tm%r0(i  ,jm1) = tm%r0(i  ,jm1) + cx   * cm1y * weight
  tm%r0(i  ,j  ) = tm%r0(i  ,j  ) + cx   * cy   * weight
  tm%r0(i  ,jp1) = tm%r0(i  ,jp1) + cx   * cp1y * weight
  tm%r0(i  ,jp2) = tm%r0(i  ,jp2) + cx   * cp2y * weight

  tm%r0(ip1,jm2) = tm%r0(ip1,jm2) + cp1x * cm2y * weight
  tm%r0(ip1,jm1) = tm%r0(ip1,jm1) + cp1x * cm1y * weight
  tm%r0(ip1,j  ) = tm%r0(ip1,j  ) + cp1x * cy   * weight
  tm%r0(ip1,jp1) = tm%r0(ip1,jp1) + cp1x * cp1y * weight
  tm%r0(ip1,jp2) = tm%r0(ip1,jp2) + cp1x * cp2y * weight

  tm%r0(ip2,jm2) = tm%r0(ip2,jm2) + cp2x * cm2y * weight
  tm%r0(ip2,jm1) = tm%r0(ip2,jm1) + cp2x * cm1y * weight
  tm%r0(ip2,j  ) = tm%r0(ip2,j  ) + cp2x * cy   * weight
  tm%r0(ip2,jp1) = tm%r0(ip2,jp1) + cp2x * cp1y * weight
  tm%r0(ip2,jp2) = tm%r0(ip2,jp2) + cp2x * cp2y * weight

end do

if (bcname == 'period') then
  do i=0,nx
    tm%r0(i,0)  = tm%r0(i,0) + tm%r0(i,ny)
    tm%r0(i,ny) = tm%r0(i,0)
  end do
  do j=0,ny
    tm%r0(0,j)  = tm%r0(0,j) + tm%r0(nx,j)
    tm%r0(nx,j) = tm%r0(0,j)
  end do
end if

tm%r0 = tm%r0 / (dx*dy) 

rho_total = sum(tm%r0(1:nx,1:ny))*dx*dy
print*,'rho total',rho_total

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


subroutine calcul_rho_m3( ele, tm)
type(particle) :: ele
type(tm_mesh_fields) :: tm

sll_int32  :: i, j, k
sll_real32 :: dpx, dpy
sll_real64 :: weight
sll_real32 :: c1x, c_1x, c2x
sll_real32 :: c1y, c_1y, c2y


tm%r0 = 0.0_f64

do k = 1, nbpart

  i      = ele%idx(k)
  j      = ele%idy(k)
  dpx    = ele%dpx(k)
  dpy    = ele%dpy(k)
  weight = ele%p(k) 

  BSPLINES(dpx,dpy)

  tm%r0( i  ,j  ) = tm%r0( i  ,j  ) + c1x*c1y   * weight
  tm%r0( i  ,j-1) = tm%r0( i  ,j-1) + c1x*c_1y  * weight
  tm%r0( i  ,j+1) = tm%r0( i  ,j+1) + c1x*c2y   * weight
  tm%r0( i-1,j  ) = tm%r0( i-1,j  ) + c_1x*c1y  * weight
  tm%r0( i-1,j-1) = tm%r0( i-1,j-1) + c_1x*c_1y * weight
  tm%r0( i-1,j+1) = tm%r0( i-1,j+1) + c_1x*c2y  * weight
  tm%r0( i+1,j  ) = tm%r0( i+1,j  ) + c2x*c1y   * weight
  tm%r0( i+1,j-1) = tm%r0( i+1,j-1) + c2x*c_1y  * weight
  tm%r0( i+1,j+1) = tm%r0( i+1,j+1) + c2x*c2y   * weight

end do

tm%r0 = tm%r0 / (dx*dy) 

end subroutine calcul_rho_m3

end module particules
