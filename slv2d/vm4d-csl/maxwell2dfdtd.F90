module maxwell2dfdtd_module 

#include "sll_working_precision.h"
use geometry_module
use fft_module

implicit none
sll_int32,  private :: i, j
sll_real64, private, parameter :: e0  = 1d0
sll_real64, private, parameter :: c = 1d0, csq = 1d0
sll_real64, private :: dex_dx, dey_dy, dex_dy, dey_dx
sll_real64, private :: dbz_dx, dbz_dy
public :: initialize, dealloc, solve_Ampere, solve_faraday
type, public:: maxwell2dfdtd
   type(geometry) :: geomx
   logical :: transpose
   sll_int32 :: jstartx,jendx,istartk,iendk
end type maxwell2dfdtd

interface initialize
   module procedure new_maxwell2dfdtd
end interface

interface dealloc
   module procedure dealloc_maxwell2dfdtd
end interface

interface solve_ampere
   module procedure solve_ampere2dfdtd
end interface

interface solve_faraday
   module procedure solve_faraday2dfdtd
end interface

contains

subroutine new_maxwell2dfdtd(this,geomx,iflag,jstartx,jendx)

type(maxwell2dfdtd),intent(out) :: this
type(geometry),intent(in)  :: geomx
sll_int32, intent(out) :: iflag
sll_int32, intent(in), optional ::  jstartx, jendx
sll_int32 :: ierr ! indicateur d'erreur
sll_int32 :: nxh1
! indicateur d'erreur
iflag = 0
this%transpose=.false.
! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
! le decoupage des tableaux exterieurs n'est pas gere par le module
! jstart et jend sont donnees en parametre d'entree dans le cas parallele
! on ne resout pas Maxwell en parallele pour l'instant
this%jstartx = 1
this%jendx = geomx%ny
! la taille totale de la zone en kx est nxh1
! ext et eyt. Ce sont des complexes. Pour cette raison leur taille est
! la moitie de celle des tableaux reels correspondants.
nxh1 = geomx%nx/2
this%istartk = 1  ! cas sequentiel
this%iendk = nxh1+1 ! cas sequentiel
!print*,'zones k ',this%istartk,this%iendk
!print*,'zones x ',this%jstartx,this%jendx

! initialisation de la geometrie
this%geomx=geomx

end subroutine new_maxwell2dfdtd

subroutine dealloc_maxwell2dfdtd(this)
type(maxwell2dfdtd),intent(out) :: this
end subroutine dealloc_maxwell2dfdtd

subroutine solve_faraday2dfdtd(this,ex,ey,bz,dt)
type(maxwell2dfdtd) :: this
sll_real64, dimension(:,:) :: ex,ey,bz
sll_real64, intent(in) :: dt	!pas de temps

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  bz(i,j) = Bz(x_{i+1/2},y_{j+1/2}) 
!*** a l'instant n puis n+1/2 
!*** on connait ex(i,j) = Ex(x_{i+1/2},y_j) et ey(i,j) = Ey(x_i,y_{j+1/2})

do i=1,this%geomx%nx-1
   do j=1,this%geomx%ny-1
      dex_dy  = (ex(i,j+1)-ex(i,j)) / this%geomx%dy
      dey_dx  = (ey(i+1,j)-ey(i,j)) / this%geomx%dx
      bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
   end do
end do

end subroutine solve_faraday2dfdtd

subroutine c_l_periodiques(this,ex,ey,bz,jx,jy,dt)
type(maxwell2dfdtd) :: this
sll_real64, intent(in) :: dt
sll_real64, dimension(:,:) :: ex,ey,bz,jx,jy
sll_int32 :: mx, my

mx = this%geomx%nx
my = this%geomx%ny
!Conditions limites periodiques

do i = 1, mx
   bz(i,my) = bz(i,1)
   dbz_dy = (bz(i,my)-bz(i,my-1)) / this%geomx%dy
!   ex(i,1) = ex(i,my) + dt*csq*dbz_dy - dt*jx(i,my)/e0
   ex(i,1) = ex(i,1) + dt*csq*dbz_dy - dt*jx(i,1)/e0
end do
     
do j = 1, my
   bz(mx,j) = bz(1,j)
   dbz_dx = (bz(mx,j)-bz(mx-1,j)) / this%geomx%dx
!   ey(1,j) = ey(mx,j) - dt*csq*dbz_dx - dt*jy(mx,j)/e0
   ey(1,j) = ey(1,j) - dt*csq*dbz_dx - dt*jy(1,j)/e0
end do

end subroutine c_l_periodiques

subroutine solve_ampere2dfdtd(this,ex,ey,bz,jx,jy,nrj,dt)

type(maxwell2dfdtd) :: this
sll_real64, dimension(:,:), intent(inout) :: ex,ey
sll_real64, dimension(:,:), intent(in)    :: bz,jx,jy
sll_real64, intent(in)  :: dt
sll_real64, intent(out) :: nrj
sll_real64, dimension(1:this%geomx%nx,1:this%geomx%ny) :: tmparray
sll_real64, parameter   :: csq = 1d0, e0 = 1d0
sll_real64 :: dx, dy

!*** Calcul du champ electrique E au temps n+1
!*** sur les points internes du maillage
!*** Ex aux points (i+1/2,j)
!*** Ey aux points (i,j+1/2)
dx = this%geomx%dx
dy = this%geomx%dy

do i=1,this%geomx%nx
   do j=2,this%geomx%ny
      dbz_dy  = (bz(i,j)-bz(i,j-1)) / dy
      ex(i,j) = ex(i,j) + dt*csq*dbz_dy - dt*jx(i,j)/e0
   end do
end do

do i=2,this%geomx%nx
   do j=1,this%geomx%ny
      dbz_dx  = (bz(i,j)-bz(i-1,j)) / dx
      ey(i,j) = ey(i,j) - dt*csq*dbz_dx - dt*jy(i,j)/e0
   end do
end do

!!$tmparray=ex
!!$do j=1,this%geomx%ny
!!$   do i=2,this%geomx%nx
!!$      tmparray(i,j)=0.5_8*(ex(i-1,j)+ex(i,j))
!!$   enddo
!!$   tmparray(1,j)=0.5_8*(ex(this%geomx%nx,j)+ex(1,j))
!!$enddo
!!$
!!$
!!$nrj=sum(tmparray*tmparray+ey*ey)*dx*dy
!!$nrj=0.5_wp*log(nrj)
!!$print *,'max ey',maxval(ey),maxval(bz)
!!$!nrj=sum(ex*ex+ey*ey)*dx*dy
!!$!nrj=0.5_wp*log(nrj)

end subroutine solve_ampere2dfdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_silver_muller(this, ex, ey, bz, jx, jy, dt)

type(maxwell2dfdtd) :: this
real(8) :: a11,a12,a21,a22,b1,b2,dis,dt
real(8), dimension(:,:) :: ex, ey, bz, jx, jy
sll_int32 :: mx, my

!Conditions de Silver-Muller
!------------------------------------
!Ey = -c Bz sur la frontiere ouest
!Ey =  c Bz sur la frontiere est
!Ex = -c Bz sur la frontiere nord
!Ex =  c Bz sur la frontiere sud

!On effectue le calcul de B sur les points fictifs du maillage
!simultanement avec la prise en compte des conditions limites sur
!E. Pour cela, on resout sur les points frontieres, l'equation de la
!condition limite en moyennant en temps pour E et en espace pour B puis
!l'equation d'Ampere

!On resout pour chaque frontiere
!---
! a11*x1 + a12*x2 = b1
! a21*x1 + a22*x2 = b2
!---
mx = this%geomx%nx
my = this%geomx%ny

!Frontiere Nord : Ex = -c Bz 
do i = 1, mx
      
   a11 = 1.;	a12 = + c
   a21 = 1./dt;	a22 = - csq / this%geomx%dy
   b1  = - ex(i,my) - c * bz(i,my-1)
   b2  =   ex(i,my)/dt - csq/this%geomx%dy*bz(i,my-1) - jx(i,my)/e0
      
   dis = a11*a22-a21*a12 
      
   !ex(i,my) = (b1*a22-b2*a12)/dis
   bz(i,my) = (a11*b2-a21*b1)/dis
      
end do
   
!Frontiere Sud : Ex =  c Bz
do i = 1, mx
      
   a11 = 1.;	a12 = - c
   a21 = 1./dt;	a22 = csq / this%geomx%dy
   b1  = - ex(i,1) + c * bz(i,1)
   b2  = ex(i,1)/dt + csq / this%geomx%dy * bz(i,1) - jx(i,1)/e0
         
   dis = a11*a22-a21*a12 
      
   ex(i,1) = (b1*a22-b2*a12)/dis
   !bz(i,0) = (a11*b2-a21*b1)/dis
      
end do
   
!Frontiere Est : Ey =  c Bz
do j = 1,my
      
   a11 = 1.;	a12 = - c
   a21 = 1./dt; a22 = + csq / this%geomx%dx
   b1  = - ey(mx,j) + c * bz(mx-1,j)
   b2  = ey(mx,j)/dt + csq/this%geomx%dx*bz(mx-1,j) - jy(mx,j)/e0
      
   dis = a11*a22-a21*a12 
      
   !ey(mx,j) = (b1*a22-b2*a12)/dis
   bz(mx,j) = (a11*b2-a21*b1)/dis
   
end do
   
!Frontiere Ouest : Ey = -c Bz
do j = 1, my
   
   a11 = 1.;	a12 = + c
   a21 = 1./dt;	a22 = - csq / this%geomx%dx
   b1  = - ey(1,j) - c * bz(1,j)
   b2  =   ey(1,j)/dt - csq/this%geomx%dx*bz(1,j) -jy(1,j)
   
   dis = a11*a22-a21*a12 

   ey(1,j) = (b1*a22-b2*a12)/dis
   !bz(0,j) = (a11*b2-a21*b1)/dis
      
end do

end subroutine cl_silver_muller

end module maxwell2dfdtd_module 
