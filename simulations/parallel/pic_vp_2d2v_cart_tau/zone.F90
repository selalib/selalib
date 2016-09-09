module zone
#include "sll_working_precision.h"
use sll_m_working_precision

private

type, public :: tm_mesh_fields
  sll_real64, dimension(:,:), pointer :: ex
  sll_real64, dimension(:,:), pointer :: ey
  sll_real64, dimension(:,:), pointer :: bz
  sll_real64, dimension(:,:), pointer :: r0
end type tm_mesh_fields

sll_real64, public :: pi 

logical :: relativ 

sll_int32,  public :: nx         ! Number of cells along x
sll_int32,  public :: ny         ! Number of cells along y
sll_int32,  public :: nstep      ! Time step number
sll_int32,  public :: nbpart     ! Number of particles
sll_real64, public :: dimx       ! Domain length along x
sll_real64, public :: dimy       ! Domain length along y
sll_real64, public :: dx         ! dimx / nx
sll_real64, public :: dy         ! dimy / ny
sll_real64, public :: dt         ! Time step
sll_real64, public :: alpha      ! Perturbation amplitude
sll_real64, public :: kx         ! Perturbation wave number along x
sll_real64, public :: ky         ! Perturbation wave number along y
sll_real64, public :: poids      ! Size of particle
sll_int32,  public :: ntau
sll_real64, public :: ep

sll_int32  :: nstepmax   ! Time step number (max)
sll_int32  :: idiag      ! Diagnostic interval
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
sll_real64 :: q_sur_m    ! charge / mass

public readin

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
                 exext,  & !champ electrique exterieur
                 eyext,  & !champ electrique exterieur
                 bzext,  & !champ magnetique exterieur
                charge,  & !charge d'une macroparticule
                 masse,  & !masse d'une macroparticule
                     c,  & !vitesse de la lumiere
                    e0,  & !permittivite du vide
               relativ,  & !calcul relativiste de la vitesse
                  ntau,  & !nuber of tau iterations
                    ep,  & !epsilon
                nbpart,  & !number of particles
                    kx,  & !wave number along x
                    ky,  & !wave number along y
                 alpha     !amplitude of perturbation

!*** Initialisation des valeurs pas default


nstepmax = 2000          ! nbre d'iterations maxi
dimx     = 1.0_f64
dimy     = 1.0_f64       ! dimensions du domaine 
nx       = 120           ! nombre de pts suivant x
ny       = 10            ! nombre de pts suivant y
cfl      = 0.9_f64       ! nombre de Courant-Friedrich-Levy
tfinal   = 10.0_f64      ! temps final
idiag    = 10            ! frequence des sorties graphiques
exext    = 0.0_f64       ! champ electrique exterieur suivant x
eyext    = 0.0_f64       ! champ electrique exterieur suivant y
bzext    = 0.0_f64       ! champ magnetique exterieur
charge   = 1.0_f64       ! charge d'une macro particule
masse    = 1.0_f64       ! masse d'une macro particule
c        = 8.0_f64       ! vitesse de la lumiere
e0       = 1.0_f64       ! permittivite du vide
relativ  = .false.       ! relativistic pusher or not
ntau     = 32
nbpart   = 204800   
ep       = 0.1_f64/1._f64
alpha    = 0.05d0  !original it's 0.10_f64
kx       = 0.50_f64
ky       = 1.0d0   !original it's 0.0_f64


pi = 4.0_f64 * atan(1.0_f64)

write(*,*) " Input file name :"// filename
open(93,file=filename,status='old')
read(93,donnees) 
close(93)

csq     = c*c
q_sur_m = charge / masse
poids   = charge

dimx  = 2*pi/kx
dimy  = 2*pi/ky  ! original it's 1
poids = dimx * dimy ! car int(f0) = dimx*dimy

dx = dimx / nx
dy = dimy / ny

dt=1.0d-1/1.0d0

nstep = floor(tfinal/dt)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
