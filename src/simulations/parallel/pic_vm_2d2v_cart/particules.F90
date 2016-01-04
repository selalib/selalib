module particules
#include "sll_working_precision.h"
#include "sll_memory.h"
use zone
use marsaglia
use quietstart

implicit none

sll_int32, private :: ipart 
sll_int32, private :: i, j, k
sll_int32, private :: ivarx, ivary
sll_real64, private :: varx, vary

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine creapa( ele, time )

type  (particle) :: ele
sll_real64 :: time

select case (nomcas)

case ("viry__")

   nbpart = 1
   allocate(ele%pos(nbpart,2))
   allocate(ele%case(nbpart,2))
   allocate(ele%vit(nbpart,2))
   allocate(ele%epx(nbpart))
   allocate(ele%epy(nbpart))
   allocate(ele%bpz(nbpart))
   allocate(ele%p(nbpart))
   
   ele%pos(1,1) = 0.1*dimx 
   ele%pos(1,2) = 0.1*dimy
   
   ele%vit(1,1) = 0.0d0         ! 3d0*c
   ele%vit(1,2) = 0.d0

   i = 0
   do while (ele%pos(1,1) >= x(i)) 
      i=i+1
   enddo
   ele%case(1,1) = i-1
   
   j = 0
   do while (ele%pos(1,2) >= y(j)) 
      j=j+1 
   enddo
   ele%case(1,2) = j-1

   ele%p(1) = poids


case("faisce")

   call faisceau_aleatoire( ele, time )

case("gaussv")

   call gaussienne( ele, time )

case("plasma")

   call plasma( ele, time )

end select

end subroutine creapa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_eb( tm1, ele )

type  (particle) :: ele
type(tm_mesh_fields) :: tm1
real(kind = prec) :: a1, a2, a3, a4
real(kind = prec) :: xp, yp, dum
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)

   dum = 1./(hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum

   ele%epx(ipart) = a1 * tm1%ex(i,j) + a2 * tm1%ex(i+1,j) &
        & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i,j+1) 
   ele%epy(ipart) = a1 * tm1%ey(i,j) + a2 * tm1%ey(i+1,j) &
        & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i,j+1) 
   ele%bpz(ipart) =  a1 * tm1%bz(i,j) + a2 * tm1%bz(i+1,j) &
        & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i,j+1) 
end do

end subroutine interpol_eb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_vitesse( ele )

type (particle) :: ele
sll_real64 :: dum, u2
sll_real64 :: tantheta, sintheta
sll_real64 :: gamma

do ipart = 1, nbpart

   !*** Changement de variable u = gamma*vit

   if( relativ ) then

      u2  = ele%vit(ipart,1)*ele%vit(ipart,1) &
          + ele%vit(ipart,2)*ele%vit(ipart,2)
      if ( u2 >= csq ) then 
         print*,'Erreur : u2 >= c2 dans le calcul de la vitesse'
         print*,'ipart = ',ipart,' vx = ',ele%vit(ipart,1),' vy = ',ele%vit(ipart,2)
         stop
      else
         gamma = 1./sqrt( 1. - u2/csq )
      endif

      ele%vit(ipart,1) = gamma*ele%vit(ipart,1)
      ele%vit(ipart,2) = gamma*ele%vit(ipart,2)

   else

      gamma=1.

   end if


   !*** Separation des effets electriques et magnetiques

   !*** On ajoute la moitie de l'effet champ electrique E

   dum = 0.5 * dt * q_sur_m
   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*(ele%epx(ipart)+exext)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*(ele%epy(ipart)+eyext)

   !*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (ele%bpz(ipart)+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta
   ele%vit(ipart,2) = ele%vit(ipart,2) - ele%vit(ipart,1)*sintheta
   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta

   !*** Autre moitie de l'effet du champ electrique E

   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*(ele%epx(ipart)+exext)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*(ele%epy(ipart)+eyext)

   !*** On repasse a la vitesse (changement de variable inverse)

   if( relativ ) then

      u2 =   ele%vit(ipart,1)*ele%vit(ipart,1) &
           + ele%vit(ipart,2)*ele%vit(ipart,2)

      gamma = sqrt( 1. + u2/csq )
 
      ele%vit(ipart,1) = ele%vit(ipart,1) / gamma
      ele%vit(ipart,2) = ele%vit(ipart,2) / gamma

   end if

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
sll_real64 :: coef

do ipart=1,nbpart     
   ele%pos(ipart,1) = ele%pos(ipart,1) + ele%vit(ipart,1)*dt*coef
   ele%pos(ipart,2) = ele%pos(ipart,2) + ele%vit(ipart,2)*dt*coef
enddo

!*** Mise a jour des "cases"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(ipart,1) >= x(i) .and. ele%pos(ipart,1)<dimx) 
      i=i+1
   enddo
   if ( ele%pos(ipart,1) >= dimx ) i=nx+1
   ele%case(ipart,1) = i-1
   j = 0
   do while (ele%pos(ipart,2) >= y(j) .and. ele%pos(ipart,2)<dimy) 
      j=j+1 
   enddo
   if ( ele%pos(ipart,2) >= dimy ) j=ny+1
   ele%case(ipart,2) = j-1
end do

end subroutine avancee_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sortie_part( ele )

type(particle) :: ele

!*** Traitement de la sortie des particules

if (bcname == 'period') then

   do ipart=1,nbpart
      if( ele%pos(ipart,1) >= dimx ) ele%pos(ipart,1) = ele%pos(ipart,1) - dimx
      if( ele%pos(ipart,2) >= dimy ) ele%pos(ipart,2) = ele%pos(ipart,2) - dimy
      if( ele%pos(ipart,1) < 0.d0 )  ele%pos(ipart,1) = ele%pos(ipart,1) + dimx
      if( ele%pos(ipart,2) < 0.d0 )  ele%pos(ipart,2) = ele%pos(ipart,2) + dimy
   end do   

else
   ipart = 1
   do while ( ipart <= nbpart )  
      if (      ele%pos(ipart,1) < 0.0        &
           .or. ele%pos(ipart,1) >= dimx      &	
           .or. ele%pos(ipart,2) < 0.         &
           .or. ele%pos(ipart,2) >= dimy ) then
         !*** Recuperation du trou laisse par la particule sortie
         ele%pos(ipart,1) = ele%pos(nbpart,1)
         ele%pos(ipart,2) = ele%pos(nbpart,2)
         ele%vit(ipart,1) = ele%vit(nbpart,1)
         ele%vit(ipart,2) = ele%vit(nbpart,2)
         ele%p(ipart)     = ele%p(nbpart)
         nbpart = nbpart - 1
         if (nbpart == 0) stop 'plus de particule'
      else
         ipart = ipart + 1
      end if
   end do
end if


!*** Mise a jour des "cases"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(ipart,1) >= x(i) .and. ele%pos(ipart,1)<=dimx) 
      i=i+1
   enddo
   ele%case(ipart,1) = i-1
   j = 0
   do while (ele%pos(ipart,2) >= y(j) .and. ele%pos(ipart,2)<=dimy) 
      j=j+1 
   enddo
   ele%case(ipart,2) = j-1
end do

end subroutine sortie_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_rho( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
sll_real64 :: a1, a2, a3, a4, dum, xp, yp
sll_real64 :: rho_total

tm%r0 = tm%r1   
tm%r1 = 0.d0    
                
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   tm%r1(i,j)     = tm%r1(i,j)     + a1/(hhx(i)*hhy(j)) !charge unite = 1
   tm%r1(i+1,j)   = tm%r1(i+1,j)   + a2/(hhx(i+1)*hhy(j)) 
   tm%r1(i+1,j+1) = tm%r1(i+1,j+1) + a3/(hhx(i+1)*hhy(j+1))
   tm%r1(i,j+1)   = tm%r1(i,j+1)   + a4/(hhx(i)*hhy(j+1))
end do

if (bcname == 'period') then
   do i=0,nx
      tm%r1(i,0)  = tm%r1(i,0) + tm%r1(i,ny)
      tm%r1(i,ny) = tm%r1(i,0)
   end do
   do j=0,ny
      tm%r1(0,j)  = tm%r1(0,j) + tm%r1(nx,j)
      tm%r1(nx,j) = tm%r1(0,j)
   end do
end if

if (nomcas == 'gaussv' .or. nomcas == 'plasma' ) then
   rho_total = 0.d0
   do i=0,nx-1
      do j=0,ny-1
         rho_total = rho_total + tm%r1(i,j)*hhx(i)*hhy(j)
      enddo
   enddo
   print*,'rho total',rho_total 
   ! Neutralisation du milieu
   tm%r1 = tm%r1 - rho_total/dimx/dimy
endif

end subroutine calcul_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_j_cic( ele, tm )

type(particle) :: ele
type(tm_mesh_fields) :: tm
sll_real64 :: a1, a2, a3, a4, dum, xp, yp
sll_real64, dimension(0:nx,0:ny) :: jx, jy

jx = 0.d0
jy = 0.d0

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   dum = ele%vit(ipart,1) / (hx(i)*hy(j)) !charge unite = 1
   jx(i,j)     = jx(i,j)     + a1*dum  
   jx(i+1,j)   = jx(i+1,j)   + a2*dum 
   jx(i+1,j+1) = jx(i+1,j+1) + a3*dum 
   jx(i,j+1)   = jx(i,j+1)   + a4*dum 
   dum = ele%vit(ipart,2) / (hx(i)*hy(j)) 
   jy(i,j)     = jy(i,j)     + a1*dum  
   jy(i+1,j)   = jy(i+1,j)   + a2*dum 
   jy(i+1,j+1) = jy(i+1,j+1) + a3*dum 
   jy(i,j+1)   = jy(i,j+1)   + a4*dum 
end do

if (bcname == 'period') then
   do i=0,nx
      jx(i,0)  = jx(i,0) + jx(i,ny)
      jx(i,ny) = jx(i,0)
      jy(i,0)  = jy(i,0) + jy(i,ny)
      jy(i,ny) = jy(i,0)
   end do
   do j=0,ny
      jx(0,j)  = jx(0,j) + jx(nx,j)
      jx(nx,j) = jx(0,j)
      jy(0,j)  = jy(0,j) + jy(nx,j)
      jy(nx,j) = jy(0,j)
   end do
end if


do i=0,nx-1
do j=0,ny
   tm%jx(i,j) = 0.5 * (jx(i,j)+jx(i+1,j))
end do
end do

do i=0,nx
do j=0,ny-1
   tm%jy(i,j) = 0.5 * (jy(i,j)+jy(i,j+1))
end do
end do

end subroutine calcul_j_cic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faisceau_aleatoire( ele, time )

type (particle) :: ele
sll_real64 :: time, vitesse, intensite, y0, y1
sll_int32, parameter :: maxnbpart=100000
sll_int32 :: nbnewpart, error
sll_real64 :: aux

intensite = 1.d0
vitesse = 0.2
y0 = 0.4*dimy 
y1 = 0.6*dimy
nbnewpart = 30 !int( intensite/nx/vitesse )
!print*,'nb de nouvelles particules = ',nbnewpart

if( time == 0. ) then

   SLL_ALLOCATE(ele%pos(maxnbpart,2),error)
   SLL_ALLOCATE(ele%case(maxnbpart,2),error)
   SLL_ALLOCATE(ele%vit(maxnbpart,2),error)
   SLL_ALLOCATE(ele%epx(maxnbpart),error)
   SLL_ALLOCATE(ele%epy(maxnbpart),error)
   SLL_ALLOCATE(ele%bpz(maxnbpart),error)
   SLL_ALLOCATE(ele%p(maxnbpart),error)

end if

do i = 1,nbnewpart
   nbpart=nbpart+1
   aux= modulo(kiss(),int(1e6))*1.e-6
   ele%pos(nbpart,1) = 0.
   ele%pos(nbpart,2) = y0+aux*(y1-y0)
   ele%case(nbpart,1) = 0
   ele%vit(nbpart,1) = 0.0 !vitesse 
   ele%vit(nbpart,2) = 0.0
   ele%p(nbpart) = charge
   j = 0
   do while (ele%pos(nbpart,2) >= y(j) .and. ele%pos(nbpart,2)<dimy)
      j=j+1
   end do
   ele%case(nbpart,2) = j-1
end do

end subroutine faisceau_aleatoire

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faisceau_regulier( ele, time )

type (particle) :: ele
sll_real64 :: time, vitesse
sll_int32, parameter :: maxnbpart=100000
sll_int32 :: nbnewpart, error

vitesse = 0.5d0
nbnewpart = 0
print*,'nb de nouvelles particules = ',nbnewpart

if( time == 0. ) then

   SLL_ALLOCATE(ele%pos(maxnbpart,2),error)
   SLL_ALLOCATE(ele%case(maxnbpart,2),error)
   SLL_ALLOCATE(ele%vit(maxnbpart,2),error)
   SLL_ALLOCATE(ele%epx(maxnbpart),error)
   SLL_ALLOCATE(ele%epy(maxnbpart),error)
   SLL_ALLOCATE(ele%bpz(maxnbpart),error)
   SLL_ALLOCATE(ele%p(maxnbpart),error)

end if

do j = 16,24
   nbpart=nbpart+1
   nbnewpart=nbnewpart+1
   ele%pos(nbpart,1) = 0.
   ele%pos(nbpart,2) = y(j)
   ele%case(nbpart,1) = 0
   ele%case(nbpart,2) = j
   ele%vit(nbpart,1) = vitesse 
   ele%vit(nbpart,2) = 0.0
   ele%p(nbpart) = poids

end do

print*,'nb de nouvelles particules = ',nbnewpart

end subroutine faisceau_regulier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gaussienne( ele, time )

type (particle) :: ele
sll_real64 :: time, speed, theta, vth, n, aux
sll_int32 :: k, error

if( time == 0 ) then
 
   vth = 1.
   nbpart = 100*nx*ny !je joue sur le nb de part/maille
   n = 1.d0/nbpart

   SLL_ALLOCATE(ele%pos(nbpart,2),error)
   SLL_ALLOCATE(ele%case(nbpart,2),error)
   SLL_ALLOCATE(ele%vit(nbpart,2),error)
   SLL_ALLOCATE(ele%epx(nbpart),error)
   SLL_ALLOCATE(ele%epy(nbpart),error)
   SLL_ALLOCATE(ele%bpz(nbpart),error)
   SLL_ALLOCATE(ele%p(nbpart),error)


   do k=0,nbpart-1

      speed = vth * sqrt(-2 * log( (k+0.5)*n ))
      aux = trinary_reversing( k )
      theta = aux * 2 * pi

      aux =  bit_reversing( k )
      ele%pos(k+1,1) = dimx * aux
      aux = penta_reversing( k ) 
      ele%pos(k+1,2) = dimy * aux 

      i = 0
      do while (ele%pos(k+1,1) >= x(i)) 
         i=i+1
      enddo
      ele%case(k+1,1) = i-1

      j = 0
      do while (ele%pos(k+1,2) >= y(j)) 
         j=j+1 
      enddo
      ele%case(k+1,2) = j-1
      
      ele%vit(k+1,1) = speed * cos(theta)
      ele%vit(k+1,2) = speed * sin(theta) 

      ele%p(k+1) = poids * n

   enddo

end if

end subroutine gaussienne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plasma( ele, time )

type (particle) :: ele
sll_real64 :: time, speed, theta, vth, n
sll_real64 :: a, b, eps, R
sll_int32 :: k, error

if( time == 0 ) then

   eps = 1.d-12

   vth =  1.
   nbpart = 100*(nx)*(ny)
   n = 1.d0/nbpart

   SLL_ALLOCATE(ele%pos(nbpart,2),error)
   SLL_ALLOCATE(ele%case(nbpart,2),error)
   SLL_ALLOCATE(ele%vit(nbpart,2),error)
   SLL_ALLOCATE(ele%epx(nbpart),error)
   SLL_ALLOCATE(ele%epy(nbpart),error)
   SLL_ALLOCATE(ele%bpz(nbpart),error)
   SLL_ALLOCATE(ele%p(nbpart),error)

   do k=0,nbpart-1

      speed = vth * sqrt(-2 * log( (k+0.5)*n ))

      theta = trinary_reversing( k ) * 2 * pi

      a = 0; b = dimx ! 2*pi/kx 
      R = bit_reversing( k )
      call dichotomie_x(a,b,R,eps) 
      ele%pos(k+1,1) = a
      ele%pos(k+1,2) = dimy * penta_reversing( k ) 

      i = 0
      do while (ele%pos(k+1,1) >= x(i)) 
         i=i+1
      enddo
      ele%case(k+1,1) = i-1
      
      j = 0
      do while (ele%pos(k+1,2) >= y(j)) 
         j=j+1 
      enddo
      ele%case(k+1,2) = j-1
      
      ele%vit(k+1,1) = speed * cos(theta)  !
      ele%vit(k+1,2) = speed * sin(theta)  !

      ele%p(k+1) = poids * n

   enddo


end if

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
