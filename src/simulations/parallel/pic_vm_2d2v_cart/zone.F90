module zone
#include "sll_working_precision.h"
use sll_m_working_precision

type tm_mesh_fields
   sll_real64, dimension(:,:), pointer :: ex
   sll_real64, dimension(:,:), pointer :: ey
   sll_real64, dimension(:,:), pointer :: bz
   sll_real64, dimension(:,:), pointer :: r0
end type tm_mesh_fields

type particle
   sll_real64, pointer :: pos(:,:)
   sll_int32 , pointer :: case(:,:)
   sll_real64, pointer :: vit(:,:)
   sll_real64, pointer :: epx(:)
   sll_real64, pointer :: epy(:)
   sll_real64, pointer :: bpz(:)
   sll_real64, pointer :: p(:)
end type particle

logical :: relativ 

sll_real64 :: pi 

character(len=6) :: bcname = 'period' 

sll_int32 :: nx, ny
sll_int32 :: nstep, nstepmax
sll_int32 :: icrea, idiag
sll_int32 :: nbpart

sll_int32, private :: i, j

sll_real64 :: dt, alpha, kx, ky, c, csq, e0
sll_real64 :: dx, dy
sll_real64, dimension(:), pointer :: x, y

sll_real64 :: dimx, dimy
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

namelist/donnees/ dimx,  & !dimensions du domaine
                  dimy,  & 
                    nx,  & !nbre de pas
                    ny,  &
                   cfl,  & !nbre de Courant
                tfinal,  & !duree maxi
              nstepmax,  & !nbre d'iterations maxi
                 icrea,  & !frequence d'emission des particules
                 idiag,  & !frequence des diagnostics
                bcname,  & !type de conditions limites
                 exext,  & !champ electrique exterieur
                 eyext,  & !champ electrique exterieur
                 bzext,  & !champ magnetique exterieur
                charge,  & !charge d'une macroparticule
                 masse,  & !masse d'une macroparticule
                     c,  & !vitesse de la lumiere
                    e0,  & !permittivite du vide
               relativ     !calcul relativiste de la vitesse

!*** Initialisation des valeurs pas default

pi = 4.0_f64 * atan(1.)

write(*,*) " Input file name :"// filename
open(93,file=filename,status='old')
read(93,donnees) 
close(93)

csq     = c*c
q_sur_m = charge / masse
poids   = charge

alpha = 0.10_f64
kx    = 0.50_f64
ky    = 0.0_f64
dimx  = 2*pi/kx
poids = dimx * dimy ! car int(f0) = dimx*dimy

!Creation du maillage

allocate(x(-1:nx+1))  !0:nx))
allocate(y(-1:ny+1))  !0:ny))

dx = dimx / nx
dy = dimy / ny

x(0) = 0.0_f64
y(0) = 0.0_f64

do i=1,nx
  x(i) = i*dx 
enddo
do j=1,ny
  y(j) = j*dy
enddo

x(-1)   = x(0) - dx
x(nx+1) = x(nx) + dx
y(-1)   = y(0) - dy
y(ny+1) = y(ny) + dy

dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

nstep = floor(tfinal/dt)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
