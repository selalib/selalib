module Maxwell
#include "sll_working_precision.h"
use zone
implicit none

integer, private :: i, j

sll_real64, private :: dex_dx, dey_dy
sll_real64, private :: dex_dy, dey_dx
sll_real64, private :: dbz_dx, dbz_dy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( tm )

type( tm_mesh_fields ) :: tm

   !*** On utilise l'equation de Faraday sur un demi pas
   !*** de temps pour le calcul du champ magnetique  Bz 
   !*** a l'instant n puis n+1/2 apres deplacement des
   !*** particules

   do i=0,nx-1
   do j=0,ny-1
      dex_dy     = (tm%ex(i,j+1)-tm%ex(i,j)) / hy(j)
      dey_dx     = (tm%ey(i+1,j)-tm%ey(i,j)) / hx(i)
      tm%bz(i,j) = tm%bz(i,j) + 0.5 * dt * (dex_dy - dey_dx)
   end do
   end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere( tm )

type( tm_mesh_fields ) :: tm

   !*** Calcul du champ electrique E au temps n+1
   !*** sur les points internes du maillage
   !*** Ex aux points (i+1/2,j)
   !*** Ey aux points (i,j+1/2)

   do i=0,nx-1
   do j=1,ny-1
      dbz_dy = (tm%bz(i,j)-tm%bz(i,j-1)) / hhy(j)
      tm%ex(i,j) = tm%ex(i,j) + csq * dt * dbz_dy - dt * tm%jx(i,j)/e0
   end do
   end do

   do i=1,nx-1
   do j=0,ny-1
      dbz_dx = (tm%bz(i,j)-tm%bz(i-1,j)) / hhx(i)
      tm%ey(i,j) = tm%ey(i,j) - csq * dt * dbz_dx - dt * tm%jy(i,j)/e0
   end do
   end do

end subroutine ampere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine conditions_limites( tm, time )

type( tm_mesh_fields ) :: tm
sll_real64 :: a11,a12,a21,a22,b1,b2,dis
sll_real64 :: time, alpha, omega
integer :: mm=1 !parametre du cas entran

select case ( bcname )

case ("period")

do i=0,nx-1
   dbz_dy = (tm%bz(i,0)-tm%bz(i,ny-1)) / hhy(0)
   tm%ex(i,0)  = tm%ex(i,0) + csq * dt * dbz_dy - dt * tm%jx(i,0)/e0
   tm%ex(i,ny) = tm%ex(i,0)
end do

do j=0,ny-1
   dbz_dx = (tm%bz(0,j)-tm%bz(nx-1,j)) / hhx(0)
   tm%ey(0,j)  = tm%ey(0,j) - csq * dt * dbz_dx - dt * tm%jy(0,j)/e0
   tm%ey(nx,j) = tm%ey(0,j)
end do

case ("silver")

!Frontiere Ouest : Ey + c Bz = 0
do j=0,ny-1

   a11 = 1
   a12 = -csq*dt/hx(0)
   a21 = 1
   a22 = c
   b1  = tm%ey(0,j) - csq*dt/hx(0) * tm%bz(0,j) - dt * tm%jy(0,j)/e0
   b2  = - tm%ey(0,j) - c*tm%bz(0,j)

   dis = a11*a22-a21*a12 

   tm%ey(0,j) = (b1*a22-b2*a12) / dis

end do

!Frontiere Est : -Ey + c Bz = 0
do j=0,ny-1

   a11 = 1
   a12 = csq*dt/hx(nx-1)
   a21 = 1
   a22 = -c
   b1  = tm%ey(nx,j) + csq*dt/hx(nx-1) * tm%bz(nx-1,j) - dt * tm%jy(nx,j)/e0
   b2  = - tm%ey(nx,j) + c*tm%bz(nx-1,j)

   dis = a11*a22-a21*a12 

   tm%ey(nx,j) = (b1*a22-b2*a12) / dis

end do

!Frontiere Sud : -Ex + c Bz = 0
do i=0,nx-1

   a11 = 1
   a12 = csq*dt/hy(0)
   a21 = 1
   a22 = -c
   b1  = tm%ex(i,0) + csq*dt/hy(0) * tm%bz(i,0) - dt * tm%jx(i,0)/e0
   b2  = - tm%ex(i,0) + c*tm%bz(i,0)

   dis = a11*a22-a21*a12 

   tm%ex(i,0) = (b1*a22-b2*a12) / dis

end do

!Frontiere Nord : Ex + c Bz = 0
do i=0,nx-1

   a11 = 1
   a12 = -csq*dt/hy(ny-1)
   a21 = 1
   a22 = c
   b1  = tm%ex(i,ny) - csq*dt/hy(ny-1) * tm%bz(i,ny-1) - dt * tm%jx(i,ny)/e0
   b2  = - tm%ex(i,ny) - c*tm%bz(i,ny-1)

   dis = a11*a22-a21*a12 

   tm%ex(i,ny) = (b1*a22-b2*a12) / dis

end do

case("conduc")

   !Frontieres Nord et Sud : Ex = 0
   do i=0,nx-1
      tm%ex(i,0)  = 0.
      tm%ex(i,ny) = 0.
   end do
   
   !Frontieres Est et Ouest : Ey = 0
   do j=0,ny-1
      tm%ey(0,j)  = 0.
      tm%ey(nx,j) = 0.
   end do
   
case ("entran")
   
   !Frontieres Nord et Sud : Ex = 0
   do i=0,nx-1
      tm%ex(i,0)  = 0.
      tm%ex(i,ny) = 0.
   end do
   
   !Frontiere Ouest : Ey + Bz = donnees !!!onde entrant et non S-M

   alpha = 2*mm*pi/dimx
   omega = alpha

   do j=0,ny-1
      
      a11 = 1
      a12 = -csq*dt/hx(0)
      a21 = 1
      a22 = c
      b1  = tm%ey(0,j) - csq*dt/hx(0) * tm%bz(0,j) - dt * tm%jy(0,j)/e0
      b2  = - tm%ey(0,j) - c*tm%bz(0,j) + 4*cos(omega*(time+0.5*dt))
      
      dis = a11*a22-a21*a12 
      
      tm%ey(0,j) = (b1*a22-b2*a12) / dis
      
   end do


   !Frontiere Est : -Ey + Bz = 0
   do j=0,ny-1
      
      a11 = 1
      a12 = csq*dt/hx(nx-1)
      a21 = 1
      a22 = -c
      b1  = tm%ey(nx,j) + csq*dt/hx(nx-1) * tm%bz(nx-1,j) - dt * tm%jy(nx,j)/e0
      b2  = - tm%ey(nx,j) + c*tm%bz(nx-1,j)
      
      dis = a11*a22-a21*a12 
      
      tm%ey(nx,j) = (b1*a22-b2*a12) / dis
      
   end do

case ("Eincom")

   !Frontieres Nord et Sud : Ex = 0
   do i=0,nx-1
      tm%ex(i,0)  = 0.
      tm%ex(i,ny) = 0.
   end do

   
   !Frontiere Ouest : Ey + c Bz = 0
   
   alpha = 2*mm*pi/dimx
   
   do j=0,ny-1
      
      a11 = 1
      a12 = -csq*dt/hx(0)
      a21 = 1
      a22 = c
      b1  = tm%ey(0,j) - csq*dt/hx(0) * tm%bz(0,j) - dt * tm%jy(0,j)/e0
      b2  = - tm%ey(0,j) - c*tm%bz(0,j) 
      
      dis = a11*a22-a21*a12 
      
      tm%ey(0,j) = (b1*a22-b2*a12) / dis
      
   end do

   !Frontiere Est : -Ey + c Bz = donnee
   do j=0,ny-1
      
      a11 = 1
      a12 = csq*dt/hx(nx-1)
      a21 = 1
      a22 = -c
      b1  = tm%ey(nx,j) + csq*dt/hx(nx-1) * tm%bz(nx-1,j) - dt * tm%jy(nx,j)/e0
      b2  = - tm%ey(nx,j) + c*tm%bz(nx-1,j) +  4*cos(alpha*(time+0.5*dt+dimx))
      
      dis = a11*a22-a21*a12 
      
      tm%ey(nx,j) = (b1*a22-b2*a12) / dis
      
   end do

case ("Sincom")

   !Frontieres Est et Ouest : Ey = 0
   do j=0,ny-1
      tm%ey(0,j)  = 0.
      tm%ey(nx,j) = 0.
   end do
   
   !Frontiere Sud : -Ex + c Bz = donnee
   
   alpha = 2*mm*pi/dimy
   
   do i=0,nx-1
      
      a11 = 1
      a12 = csq*dt/hy(0)
      a21 = 1
      a22 = -c
      b1  = tm%ex(i,0) + csq*dt/hy(0) * tm%bz(i,0) - dt * tm%jx(i,0)/e0
      b2  = - tm%ex(i,0) + c*tm%bz(i,0) + 4*cos(alpha*(time+0.5*dt))
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,0) = (b1*a22-b2*a12) / dis
      
   end do
   
   !Frontiere Nord : Ex + c Bz = 0
   do i=0,nx-1
      
      a11 = 1
      a12 = -csq*dt/hy(ny-1)
      a21 = 1
      a22 = c
      b1  = tm%ex(i,ny) - csq*dt/hy(ny-1) * tm%bz(i,ny-1) - dt * tm%jx(i,ny)/e0
      b2  = - tm%ex(i,ny) - c*tm%bz(i,ny-1)
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,ny) = (b1*a22-b2*a12) / dis
      
   end do
   
case ("Nincom")

   !Frontieres Est et Ouest : Ey = 0
   do j=0,ny-1
      tm%ey(0,j)  = 0.
      tm%ey(nx,j) = 0.
   end do
   
   !Frontiere Sud : -Ex + c Bz = 0
   
   alpha = 2*mm*pi/dimy
   
   do i=0,nx-1
      
      a11 = 1
      a12 = csq*dt/hy(0)
      a21 = 1
      a22 = -c
      b1  = tm%ex(i,0) + csq*dt/hy(0) * tm%bz(i,0) - dt * tm%jx(i,0)/e0
      b2  = - tm%ex(i,0) + c*tm%bz(i,0) 
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,0) = (b1*a22-b2*a12) / dis
      
   end do
   
   !Frontiere Nord : Ex + c Bz = donnee
   do i=0,nx-1
      
      a11 = 1
      a12 = -csq*dt/hy(ny-1)
      a21 = 1
      a22 = c
      b1  = tm%ex(i,ny) - csq*dt/hy(ny-1) * tm%bz(i,ny-1) - dt * tm%jx(i,ny)/e0
      b2  = - tm%ex(i,ny) - c*tm%bz(i,ny-1) + 4*cos(alpha*(time+0.5*dt+dimy))
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,ny) = (b1*a22-b2*a12) / dis
      
   end do

case("faisce")

   !Frontiere Sud : -Ex + c Bz = 0 Silver Muller
   do i=0,nx-1
      
      a11 = 1
      a12 = csq*dt/hy(0)
      a21 = 1
      a22 = -c
      b1  = tm%ex(i,0) + csq*dt/hy(0) * tm%bz(i,0) - dt * tm%jx(i,0)/e0
      b2  = - tm%ex(i,0) + c*tm%bz(i,0)
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,0) = (b1*a22-b2*a12) / dis
      
   end do
   
   !Frontiere Nord : Ex + c Bz = 0 Silver Muller
   do i=0,nx-1
      
      a11 = 1
      a12 = -csq*dt/hy(ny-1)
      a21 = 1
      a22 = c
      b1  = tm%ex(i,ny) - csq*dt/hy(ny-1) * tm%bz(i,ny-1) - dt * tm%jx(i,ny)/e0
      b2  = - tm%ex(i,ny) - c*tm%bz(i,ny-1)
      
      dis = a11*a22-a21*a12 
      
      tm%ex(i,ny) = (b1*a22-b2*a12) / dis
      
   end do

   !** Est et Ouest : conducteur parfait **
   do j = 0, ny-1
      tm%ey(nx,j) = 0.d0
      tm%ey(0,j) = 0.d0
   end do

end select

end subroutine conditions_limites

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine decalage( tm, tm1 )

type(tm_mesh_fields) :: tm, tm1

!*** Calcul des composantes des champs 
!*** sur les noeuds du maillage de rho
!*** par interpolation lineaire

do i=1,nx-1
   do j=1,ny-1
      tm1%ex(i,j) = ( hx(i)*tm%ex(i-1,j) + hx(i-1)*tm%ex(i,j) ) &
           & / (hx(i)+hx(i-1))
      tm1%ey(i,j) = ( hy(j)*tm%ey(i,j-1) + hy(j-1)*tm%ey(i,j) ) &
           & / (hy(j)+hy(j-1))
      tm1%bz(i,j) = ( ( hx(i)*tm%bz(i-1,j-1) + hx(i-1)*tm%bz(i,j-1) ) &
           & * hy(j) + ( hx(i)*tm%bz(i-1,j) + hx(i-1)*tm%bz(i,j) ) &
           & * hy(j-1) ) / ( (hx(i)+hx(i-1)) * (hy(j)+hy(j-1)) )
   end do
end do

if (bcname == 'period') then

   do i=1,nx-1 !Sud et Nord
      tm1%ex(i,0) = ( hx(i)*tm%ex(i-1,0) + hx(i-1)*tm%ex(i,0) ) &
              & / (hx(i)+hx(i-1))
      tm1%ey(i,0) = ( hy(0)*tm%ey(i,ny-1) + hy(ny-1)*tm%ey(i,0) ) &
              & / (hy(0)+hy(ny-1))
      tm1%bz(i,0) = ( ( hx(i)*tm%bz(i-1,ny-1) + hx(i-1)*tm%bz(i,ny-1) ) &
              & * hy(0) + ( hx(i)*tm%bz(i-1,0) + hx(i-1)*tm%bz(i,0) ) &
              & * hy(ny-1) ) / ( (hx(i)+hx(i-1)) * (hy(0)+hy(ny-1)) )
      tm1%ex(i,ny) = tm1%ex(i,0) 
      tm1%ey(i,ny) = tm1%ey(i,0) 
      tm1%bz(i,ny) = tm1%bz(i,0) 
   end do

   do j=1,ny-1 !Ouest et Est
      tm1%ex(0,j) = ( hx(0)*tm%ex(nx-1,j) + hx(nx-1)*tm%ex(0,j) ) &
           & / (hx(0)+hx(nx-1))
      tm1%ey(0,j) = ( hy(j)*tm%ey(0,j-1) + hy(j-1)*tm%ey(0,j) ) &
           & / (hy(j)+hy(j-1))
      tm1%bz(0,j) = ( ( hx(0)*tm%bz(nx-1,j-1) + hx(nx-1)*tm%bz(0,j-1) ) &
           & * hy(j) + ( hx(0)*tm%bz(nx-1,j) + hx(nx-1)*tm%bz(0,j) ) &
           & * hy(j-1) ) / ( (hx(0)+hx(nx-1)) * (hy(j)+hy(j-1)) )
      tm1%ex(nx,j) = tm1%ex(0,j) 
      tm1%ey(nx,j) = tm1%ey(0,j) 
      tm1%bz(nx,j) = tm1%bz(0,j) 
   end do

   !Coins
   tm1%ex(0,0) = ( hx(0)*tm%ex(nx-1,0) + hx(nx-1)*tm%ex(0,0) ) &
           & / (hx(0)+hx(nx-1))
   tm1%ey(0,0) = ( hy(0)*tm%ey(0,ny-1) + hy(ny-1)*tm%ey(0,0) ) &
           & / (hy(0)+hy(ny-1))
   tm1%bz(0,0) = ( ( hx(0)*tm%bz(nx-1,ny-1) + hx(nx-1)*tm%bz(0,ny-1) ) &
           & * hy(0) + ( hx(0)*tm%bz(nx-1,0) + hx(nx-1)*tm%bz(0,0) ) &
           & * hy(ny-1) ) / ( (hx(0)+hx(nx-1)) * (hy(0)+hy(ny-1)) )

   tm1%ex(nx,0)  = tm1%ex(0,0) 
   tm1%ex(nx,ny) = tm1%ex(0,0) 
   tm1%ex(0,ny)  = tm1%ex(0,0) 

   tm1%ey(nx,0)  = tm1%ey(0,0) 
   tm1%ey(nx,ny) = tm1%ey(0,0) 
   tm1%ey(0,ny)  = tm1%ey(0,0) 

   tm1%bz(nx,0)  = tm1%bz(0,0) 
   tm1%bz(nx,ny) = tm1%bz(0,0) 
   tm1%bz(0,ny)  = tm1%bz(0,0) 

else

   do i=1,nx-1 
      ! Sud 
      tm1%ex(i,0) = ( hx(i)*tm%ex(i-1,0) + hx(i-1)*tm%ex(i,0) ) &
              & / (hx(i)+hx(i-1))
      tm1%ey(i,0) = ( -hy(0)*tm%ey(i,1) + (hy(1)+2*hy(0))*tm%ey(i,0) ) &
              & / (hy(0)+hy(1))
      tm1%bz(i,0) = ( ( -hy(0)*tm%bz(i-1,1) + (2*hy(0)+hy(1))*tm%bz(i-1,0) ) &
           & * hx(i-1) + ( -hy(0)*tm%bz(i,1) + (2*hy(0)+hy(1))*tm%bz(i,0) ) &
           & * hx(i) ) / ( (hx(i-1)+hx(i)) * (hy(0)+hy(1)) )
      ! Nord
      tm1%ex(i,ny) = ( hx(i)*tm%ex(i-1,ny) + hx(i-1)*tm%ex(i,ny) ) &
           & / (hx(i)+hx(i-1))
      tm1%ey(i,ny) = ( -hy(ny-1)*tm%ey(i,ny-2) + (hy(ny-2)+2*hy(ny-1)) &
           & *tm%ey(i,ny-1) ) / (hy(ny-1)+hy(ny-2))
      tm1%bz(i,ny) = ( ( -hy(ny-1)*tm%bz(i-1,ny-2) + (2*hy(ny-1)+hy(ny-2)) &
           & *tm%bz(i-1,ny-1) ) * hx(i-1) + ( -hy(ny-1)*tm%bz(i,ny-2) &
           & + (2*hy(ny-1)+hy(ny-2))*tm%bz(i,ny-1) )  * hx(i) ) &
           & / ( (hx(i-1)+hx(i)) * (hy(ny-1)+hy(ny-2)) )
   end do

   do j=1,ny-1 
      ! Ouest 
      tm1%ex(0,j) = ( -hx(0)*tm%ex(1,j) + (2*hx(0)+hx(1))*tm%ex(0,j) ) &
           & / (hx(0)+hx(1))
      tm1%ey(0,j) = ( hy(j)*tm%ey(0,j-1) + hy(j-1)*tm%ey(0,j) ) &
           & / (hy(j)+hy(j-1))
      tm1%bz(0,j) = ( ( -hx(0)*tm%bz(1,j-1) + (2*hx(0)+hx(1))*tm%bz(0,j-1) ) &
           & * hy(j) + ( -hx(0)*tm%bz(1,j) + (2*hx(0)+hx(1))*tm%bz(0,j) ) &
           & * hy(j-1) ) / ( (hx(0)+hx(1)) * (hy(j)+hy(j-1)) )
      ! Est
      tm1%ex(nx,j) = ( -hx(nx-1)*tm%ex(nx-2,j) + (2*hx(nx-1)+hx(nx-2)) &
           & *tm%ex(nx-1,j) ) / (hx(nx-1)+hx(nx-2))
      tm1%ey(nx,j) = ( hy(j)*tm%ey(nx,j-1) + hy(j-1)*tm%ey(nx,j) ) &
           & / (hy(j)+hy(j-1))
      tm1%bz(nx,j) =  ( ( -hx(nx-1)*tm%bz(nx-2,j-1) + (2*hx(nx-1)+hx(nx-2)) &
           & *tm%bz(nx-1,j-1) ) * hy(j) + ( -hx(nx-1)*tm%bz(nx-2,j) &
           & + (2*hx(nx-1)+hx(nx-2))*tm%bz(nx-1,j) ) * hy(j-1) ) &
           & / ( (hx(nx-1)+hx(nx-2)) * (hy(j)+hy(j-1)) )
   end do

   ! Coin SW
   tm1%ex(0,0) = ( -hx(0)*tm%ex(1,0) + (2*hx(0)+hx(1))*tm%ex(0,0) ) &
        & / (hx(0)+hx(1))
   tm1%ey(0,0) = ( -hy(0)*tm%ey(0,1) + (hy(1)+2*hy(0))*tm%ey(0,0) ) &
        & / (hy(0)+hy(1))
   tm1%bz(0,0) = ( ( -hx(0)*tm%bz(1,0) + (2*hx(0)+hx(1))*tm%bz(0,0) ) &
           & * (2*hy(0)+hy(1)) + ( -hx(0)*tm%bz(1,1) + (2*hx(0)+hx(1)) & 
           & *tm%bz(0,1) ) * (-hy(0)) ) / ( (hx(0)+hx(1)) * (hy(0)+hy(1)) )

   ! Coin NW
   tm1%ex(0,ny) = ( -hx(0)*tm%ex(1,ny) + (2*hx(0)+hx(1))*tm%ex(0,ny) ) &
           & / (hx(0)+hx(1))
   tm1%ey(0,ny) = ( -hy(ny-1)*tm%ey(0,ny-2) + (hy(ny-2)+2*hy(ny-1)) &
           & *tm%ey(0,ny-1) ) / (hy(ny-1)+hy(ny-2))
   tm1%bz(0,ny) = ( ( -hx(0)*tm%bz(1,ny-1) + (2*hx(0)+hx(1))*tm%bz(0,ny-1) ) &
           & * (2*hy(ny-1)+hy(ny-2)) + ( -hx(0)*tm%bz(1,ny-2) &
           & + (2*hx(0)+hx(1)) *tm%bz(0,ny-2) ) * (-hy(ny-1)) ) &
           & / ( (hx(0)+hx(1)) * (hy(ny-1)+hy(ny-2)) )

   ! Coin SE
   tm1%ex(nx,0) = ( -hx(nx-1)*tm%ex(nx-2,0) + (2*hx(nx-1)+hx(nx-2)) &
        & *tm%ex(nx-1,0) ) / (hx(nx-1)+hx(nx-2))
   tm1%ey(nx,0) = ( -hy(0)*tm%ey(nx,1) + (hy(1)+2*hy(0))*tm%ey(nx,0) ) &
        & / (hy(0)+hy(1))
   tm1%bz(nx,0) = ( ( -hx(nx-1)*tm%bz(nx-2,0) + (2*hx(nx-1)+hx(nx-2)) &
        & *tm%bz(nx-1,0) )  * (2*hy(0)+hy(1)) + ( -hx(nx-1)*tm%bz(nx-2,1) &
        & + (2*hx(nx-1)+hx(nx-2)) *tm%bz(nx-1,1) ) * (-hy(0)) ) &
        & / ( (hx(nx-1)+hx(nx-2)) * (hy(0)+hy(1)) )

   ! Coin NE
   tm1%ex(nx,ny) = ( -hx(nx-1)*tm%ex(nx-2,ny) + (2*hx(nx-1)+hx(nx-2)) &
        & *tm%ex(nx-1,ny) ) / (hx(nx-1)+hx(nx-2))
   tm1%ey(nx,ny) = ( -hy(ny-1)*tm%ey(nx,ny-2) + (hy(ny-2)+2*hy(ny-1)) &
        & *tm%ey(nx,ny-1) ) / (hy(ny-1)+hy(ny-2))
   tm1%bz(nx,0) = ( ( -hx(nx-1)*tm%bz(nx-2,ny-1) + (2*hx(nx-1)+hx(nx-2)) &
        & *tm%bz(nx-1,ny-1) )  * (2*hy(ny-1)+hy(ny-2)) &
        & + ( -hx(nx-1)*tm%bz(nx-2,ny-2) &
        & + (2*hx(nx-1)+hx(nx-2)) *tm%bz(nx-1,ny-2) ) * (-hy(ny-1)) ) &
        & / ( (hx(nx-1)+hx(nx-2)) * (hy(ny-1)+hy(ny-2)) )

end if

end subroutine decalage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Maxwell
