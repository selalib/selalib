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
sll_real64 :: xp, yp, dum
sll_int32  :: k 
sll_int32  :: i, j
sll_real64 :: dpx
sll_real64 :: dpy
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

dum = 1./(dx*dy)

do k=1,nbpart

   i = ele%idx(k)
   j = ele%idy(k)
   dpx = ele%dpx(k)
   dpy = ele%dpy(k)

   a1 = (dx-dpx) * (dy-dpy) * dum
   a2 = (   dpx) * (dy-dpy) * dum
   a3 = (   dpx) * (   dpy) * dum
   a4 = (dx-dpx) * (   dpy) * dum

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

  ppx = ele%idx(k)*dx + ele%dpx(k) + ele%vpx(k)*dt*coef
  ppy = ele%idy(k)*dy + ele%dpy(k) + ele%vpy(k)*dt*coef

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
  ele%dpx(k) = ppx - ele%idx(k)*dx
  ele%dpy(k) = ppy - ele%idy(k)*dy

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
sll_real64 :: a1, a2, a3, a4, dum, xp, yp
sll_real64 :: rho_total
sll_int32  :: k 
sll_int32  :: i, j 
sll_real64 :: dpx
sll_real64 :: dpy

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

  i   = ele%idx(k)
  j   = ele%idy(k)
  dpx = ele%dpx(k)
  dpy = ele%dpy(k)
  dum = ele%p(k) 
  a1  = (dx-dpx) * (dy-dpy) * dum
  a2  = (dpx)    * (dy-dpy) * dum
  a3  = (dpx)    * (dpy)    * dum
  a4  = (dx-dpx) * (dpy)    * dum
  tm%r0(i,j)     = tm%r0(i,j)     + a1 
  tm%r0(i+1,j)   = tm%r0(i+1,j)   + a2 
  tm%r0(i+1,j+1) = tm%r0(i+1,j+1) + a3
  tm%r0(i,j+1)   = tm%r0(i,j+1)   + a4

end do

tm%r0 = tm%r0 / (dx*dy) / (dx*dy)

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
sll_real64 :: a, b, eps, R
sll_real64 :: ppx, ppy
sll_int32  :: k, error

eps = 1.d-12

vth =  1.0_f64
nbpart = 100*(nx)*(ny)
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

  a = 0.0_f64; b = dimx ! 2*pi/kx 
  R = bit_reversing( k )
  call dichotomie_x(a,b,R,eps) 
  ppx = a
  ppy = dimy * penta_reversing( k ) 

  ele%idx(k+1) = floor(ppx/dimx*nx)
  ele%idy(k+1) = floor(ppy/dimy*ny)
  
  ele%vpx(k+1) = speed * cos(theta)
  ele%vpy(k+1) = speed * sin(theta)

  ele%p(k+1) = poids * n

enddo

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
