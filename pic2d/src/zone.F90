module zone
#include "sll_working_precision.h"

sll_int32, parameter :: prec=8

type tm_mesh_fields
   sll_real64, dimension(:,:), pointer :: ex, ey
   sll_real64, dimension(:,:), pointer :: bz
   sll_real64, dimension(:,:), pointer :: r0, r1
   sll_real64, dimension(:,:), pointer :: jx, jy
end type tm_mesh_fields

type particle
   sll_real64   , pointer :: pos(:,:)
   sll_int32           , pointer :: case(:,:)
   sll_real64   , pointer :: vit(:,:)
   sll_real64   , pointer :: epx(:)
   sll_real64   , pointer :: epy(:)
   sll_real64   , pointer :: bpz(:)
   sll_real64   , pointer :: p(:)
end type particle

logical :: relativ 

sll_real64 :: pi 

character(len=6) :: nomcas
character(len=6) :: bcname 
character(len=6) :: jname


sll_int32 :: nx, ny
sll_int32 :: nstep, nstepmax
sll_int32 :: icrea, idiag
sll_int32 :: nbpart

sll_int32, private :: i, j

sll_real64 :: dt, alpha, kx, ky, c, csq, e0
sll_real64, private :: dx, dy, dx1, dy1, dx2, dy2
sll_real64, dimension(:), pointer :: x, y
sll_real64, dimension(:), allocatable :: hx, hy    ! les h_i+1/2
sll_real64, dimension(:), allocatable :: hhx, hhy  ! les h_i
sll_real64 :: dimx, dimy, dimx1, dimy1, dimx2, dimy2
sll_real64 :: cfl
sll_real64 :: tfinal

sll_real64 :: exext, eyext, bzext

sll_real64 :: charge, masse, poids
sll_real64 :: q_sur_m 


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readin( filename )

implicit none

character(len=*) :: filename

namelist/donnees/ dimx,  &      !dimensions du domaine
                  dimy,  & 
                    nx,  &      !nbre de pas
                    ny,  &
                  cfl,  &      !nbre de Courant
               tfinal,  &      !duree maxi
      nstepmax,  &	!nbre d'iterations maxi
                nomcas,  &      !nom du cas ce calcul
                 jname,  &      !calcul de j   
 icrea,  &	!frequence d'emission des particules
 idiag,  &	!frequence des diagnostics
 bcname, & 	!type de conditions limites
 exext,&	!champ electrique exterieur
 eyext,&	!champ electrique exterieur
 bzext,&	!champ magnetique exterieur
               charge, &	!charge d'une macroparticule
                masse,  &      !masse d'une macroparticule
                    c,  &      !vitesse de la lumiere
                   e0,  &      !permittivite du vide
                relativ        !calcul relativiste de la vitesse

!***Initialisation  des valeurs pas default

pi = 4. * atan(1.)

open(93,file=filename,status='old')
read(93,donnees) 
close(93)

csq = c*c
q_sur_m = charge / masse
poids = charge

if (nomcas == "plasma") then
   alpha = 0.1
   kx = 0.5
   ky = 0.
   dimx = 2*pi/kx
   poids = dimx * dimy ! car int(f0) = dimx*dimy
endif

!Creation du maillage

allocate(x(-1:nx+1))  !0:nx))
allocate(y(-1:ny+1))  !0:ny))
allocate(hx(-1:nx))
allocate(hy(-1:ny))
allocate(hhx(0:nx))
allocate(hhy(0:ny))

dx = dimx / nx
dy = dimy / ny

x(0) = 0.
y(0) = 0.

if (nomcas == "plasma") then
   do i=1,nx
      x(i) = (i*dx) *(i*dx+1)/(1+dimx)
   enddo
   do j=1,ny
      y(j) = (j*dy) *(j*dy+1)/(1+dimy)
   enddo
else
   do i=1,nx
      x(i) = i*dx 
   enddo
   do j=1,ny
      y(j) = j*dy
   enddo
end if


do i=0,nx-1
   hx(i) = x(i+1)-x(i)
end do
do j=0,ny-1
   hy(j) = y(j+1)-y(j)
end do
hx(nx) = hx(0)  ! CL periodiques
hx(-1) = hx(nx-1)
hy(ny) = hy(0)
hy(-1) = hy(ny-1)

x(-1)   = x(0) - hx(nx-1)  !points utiles pour le cas period
x(nx+1) = x(nx) + hx(0)
y(-1)   = y(0) - hy(ny-1)
y(ny+1) = y(ny) + hy(0)

hhx(0) =  0.5 * ( hx(0) + hx(nx-1) )  !0.5 * hx(0)
hhx(nx) =  0.5 * ( hx(0) + hx(nx-1) )   !0.5 * hx(nx-1)
do i=1,nx-1
   hhx(i) = 0.5 * ( hx(i) + hx(i-1) )
enddo
hhy(0) = 0.5 * ( hy(0) + hy(ny-1) )   !0.5 * hy(0)
hhy(ny) = 0.5 * ( hy(0) + hy(ny-1) )   !0.5 * hy(ny-1)
do j=1,ny-1
   hhy(j) = 0.5 * ( hy(j) + hy(j-1) )
enddo

dx = hx(0)  !on calcule le plus petit pas
do i=1,nx-1
   if (hx(i)<dx)  dx = hx(i)
end do

dy = hy(0)
do j=1,ny-1
   if (hy(j)<dy)  dy = hy(j)
end do

dt    = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

nstep = floor(tfinal/dt)

!write(*,*) " cfl = ", cfl
!write(*,*) " dx = ", dx, " dy = ", dy, " dt = ", dt
!if( nstep > nstepmax ) nstep = nstepmax
!write(*,*) " Nombre d'iteration nstep = ", nstep
!write(*,*)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
