module particules
#include "sll_working_precision.h"
#include "sll_memory.h"
use zone
use quietstart

implicit none


type particle
  sll_real64, pointer :: pos(:,:)
  sll_int32 , pointer :: idx(:)
  sll_int32 , pointer :: idy(:)
  sll_real64, pointer :: vit(:,:)
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
   xp = ele%pos(k,1)
   yp = ele%pos(k,2)

   a1 = ((i+1)*dx-xp) * ((j+1)*dy-yp) * dum
   a2 = (xp-(i)*dx) * ((j+1)*dy-yp) * dum
   a3 = (xp-(i)*dx) * (yp-(j)*dy) * dum
   a4 = ((i+1)*dx-xp) * (yp-(j)*dy) * dum

   ele%epx(k) = a1 * tm1%ex(i,j) + a2 * tm1%ex(i+1,j) &
        & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i,j+1) 
   ele%epy(k) = a1 * tm1%ey(i,j) + a2 * tm1%ey(i+1,j) &
        & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i,j+1) 
   ele%bpz(k) =  a1 * tm1%bz(i,j) + a2 * tm1%bz(i+1,j) &
        & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i,j+1) 
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

      u2  = ele%vit(k,1)*ele%vit(k,1) &
          + ele%vit(k,2)*ele%vit(k,2)
      if ( u2 >= csq ) then 
         print*,'Erreur : u2 >= c2 dans le calcul de la vitesse'
         print*,'k = ',k,' vx = ',ele%vit(k,1),' vy = ',ele%vit(k,2)
         stop
      else
         gamma = 1./sqrt( 1. - u2/csq )
      endif

      ele%vit(k,1) = gamma*ele%vit(k,1)
      ele%vit(k,2) = gamma*ele%vit(k,2)

   else

      gamma=1.0_f64

   end if


   !*** Separation des effets electriques et magnetiques

   !*** On ajoute la moitie de l'effet champ electrique E

   dum = 0.5 * dt * q_sur_m
   ele%vit(k,1) = ele%vit(k,1) + dum*(ele%epx(k)+exext)
   ele%vit(k,2) = ele%vit(k,2) + dum*(ele%epy(k)+eyext)

   !*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (ele%bpz(k)+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vit(k,1) = ele%vit(k,1) + ele%vit(k,2)*tantheta
   ele%vit(k,2) = ele%vit(k,2) - ele%vit(k,1)*sintheta
   ele%vit(k,1) = ele%vit(k,1) + ele%vit(k,2)*tantheta

   !*** Autre moitie de l'effet du champ electrique E

   ele%vit(k,1) = ele%vit(k,1) + dum*(ele%epx(k)+exext)
   ele%vit(k,2) = ele%vit(k,2) + dum*(ele%epy(k)+eyext)

   !*** On repasse a la vitesse (changement de variable inverse)

   if( relativ ) then

      u2 =   ele%vit(k,1)*ele%vit(k,1) &
           + ele%vit(k,2)*ele%vit(k,2)

      gamma = sqrt( 1. + u2/csq )
 
      ele%vit(k,1) = ele%vit(k,1) / gamma
      ele%vit(k,2) = ele%vit(k,2) / gamma

   end if

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
sll_real64 :: coef
sll_int32  :: k 

do k=1,nbpart     
   ele%pos(k,1) = ele%pos(k,1) + ele%vit(k,1)*dt*coef
   ele%pos(k,2) = ele%pos(k,2) + ele%vit(k,2)*dt*coef
enddo

!*** Mise a jour des "ids"

do k=1,nbpart
   ele%idx(k) = floor(ele%pos(k,1)/dimx*nx)
   ele%idy(k) = floor(ele%pos(k,2)/dimy*ny)
end do

!*** Traitement de la sortie des particules

if (bcname == 'period') then

  do k=1,nbpart
    if( ele%pos(k,1) >= dimx ) then
      ele%pos(k,1) = ele%pos(k,1) - dimx
      ele%idx(k) = floor(ele%pos(k,1)/dimx*nx)
    end if
    if( ele%pos(k,2) >= dimy ) then
      ele%pos(k,2) = ele%pos(k,2) - dimy
      ele%idy(k) = floor(ele%pos(k,2)/dimy*ny)
    end if
    if( ele%pos(k,1) < 0.d0 ) then
      ele%pos(k,1) = ele%pos(k,1) + dimx
      ele%idx(k) = floor(ele%pos(k,1)/dimx*nx)
    end if
    if( ele%pos(k,2) < 0.d0 ) then
      ele%pos(k,2) = ele%pos(k,2) + dimy
      ele%idy(k) = floor(ele%pos(k,2)/dimy*ny)
    end if
  end do   

else
   k = 1
   do while ( k <= nbpart )  
      if (      ele%pos(k,1) < 0.0        &
           .or. ele%pos(k,1) >= dimx      &	
           .or. ele%pos(k,2) < 0.         &
           .or. ele%pos(k,2) >= dimy ) then
         !*** Recuperation du trou laisse par la particule sortie
         ele%pos(k,1) = ele%pos(nbpart,1)
         ele%pos(k,2) = ele%pos(nbpart,2)
         ele%vit(k,1) = ele%vit(nbpart,1)
         ele%vit(k,2) = ele%vit(nbpart,2)
         ele%p(k)     = ele%p(nbpart)
         ele%idx(k) = floor(ele%pos(k,1)/dimx*nx)
         ele%idy(k) = floor(ele%pos(k,2)/dimy*ny)
         nbpart = nbpart - 1
         if (nbpart == 0) stop 'plus de particule'
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
  xp  = ele%pos(k,1)
  yp  = ele%pos(k,2)
  dum = ele%p(k) 
  a1  = ((i+1)*dx-xp) * ((j+1)*dy-yp) * dum
  a2  = (xp-(i)*dx)   * ((j+1)*dy-yp) * dum
  a3  = (xp-(i)*dx)   * (yp-(j)*dy)   * dum
  a4  = ((i+1)*dx-xp) * (yp-(j)*dy)   * dum
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
sll_int32  :: i, j, k, error

eps = 1.d-12

vth =  1.0_f64
nbpart = 100*(nx)*(ny)
n = 1.d0/nbpart

SLL_ALLOCATE(ele%pos(nbpart,2),error)
SLL_ALLOCATE(ele%idx(nbpart),error)
SLL_ALLOCATE(ele%idy(nbpart),error)
SLL_ALLOCATE(ele%vit(nbpart,2),error)
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
  ele%pos(k+1,1) = a
  ele%pos(k+1,2) = dimy * penta_reversing( k ) 

  i = 0
  do while (ele%pos(k+1,1) >= i*dx) 
     i=i+1
  enddo
  ele%idx(k+1) = i-1
  
  j = 0
  do while (ele%pos(k+1,2) >= j*dy) 
     j=j+1 
  enddo
  ele%idy(k+1) = j-1
  
  ele%vit(k+1,1) = speed * cos(theta)  !
  ele%vit(k+1,2) = speed * sin(theta)  !

  ele%p(k+1) = poids * n

enddo

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
