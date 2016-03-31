module zone
#include "sll_working_precision.h"
use sll_m_working_precision

type tm_mesh_fields
   sll_real64, dimension(:,:), pointer :: ex
   sll_real64, dimension(:,:), pointer :: ey
   sll_real64, dimension(:,:), pointer :: bz
   sll_real64, dimension(:,:), pointer :: r0
end type tm_mesh_fields


logical :: relativ 

sll_real64 :: pi 

character(len=6) :: bcname = 'period' 

sll_int32  :: nx         ! Number of cells along x
sll_int32  :: ny         ! Number of cells along y
sll_real64 :: dimx       ! Domain length along x
sll_real64 :: dimy       ! Domain length along y
sll_real64 :: dx         ! dimx / nx
sll_real64 :: dy         ! dimy / ny

sll_int32  :: nstep      ! Time step number
sll_int32  :: nstepmax   ! Time step number (max)
sll_int32  :: idiag      ! Diagnostic interval
sll_int32  :: nbpart     ! Number of particles

sll_real64 :: dt         ! Time step
sll_real64 :: alpha      ! Perturbation amplitude
sll_real64 :: kx         ! Perturbation wave number along x
sll_real64 :: ky         ! Perturbation wave number along y
sll_real64 :: c          ! Speed of light
sll_real64 :: csq        ! c * c
sll_real64 :: e0         ! Electric conductivity

sll_real64 :: cfl        ! Courant-Friedrich-Levy coefficient
sll_real64 :: tfinal     ! time (max)
sll_real64 :: exext      ! External electric field (x)
sll_real64 :: eyext      ! External electric field (y)      
sll_real64 :: bzext      ! External magnetic field
sll_real64 :: charge     ! Particle charge
sll_real64 :: masse      ! Particle mass
sll_real64 :: poids      ! Size of particle
sll_real64 :: q_sur_m    ! charge / mass

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


nstepmax = 2000          ! nbre d'iterations maxi
dimx     = 1.0_f64
dimy     = 1.0_f64       ! dimensions du domaine 
nx       = 120           ! nombre de pts suivant x
ny       = 10            ! nombre de pts suivant y
cfl      = 0.9_f64       ! nombre de Courant-Friedrich-Levy
tfinal   = 10.0_f64      ! temps final
idiag    = 10            ! frequence des sorties graphiques
bcname   = "period"      ! type de conditions limites
exext    = 0.0_f64       ! champ electrique exterieur suivant x
eyext    = 0.0_f64       ! champ electrique exterieur suivant y
bzext    = 0.0_f64       ! champ magnetique exterieur
charge   = 1.0_f64       ! charge d'une macro particule
masse    = 1.0_f64       ! masse d'une macro particule
c        = 8.0_f64       ! vitesse de la lumiere
e0       = 1.0_f64       ! permittivite du vide
relativ  = .false.       ! relativistic pusher or not


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

dx = dimx / nx
dy = dimy / ny

dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

nstep = floor(tfinal/dt)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
