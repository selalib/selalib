!test example to visualize particles
module biot_savart
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
#include "sll_constants.h"
#include "sll_file_io.h"

implicit none

!> vortex circulation
real(8) :: gam0
!> angular velocity of vortex pair
real(8) :: gomeg
!> gaussian or constant distribution
logical :: gauss = .true.

contains

!> Compute particles velocities
subroutine vitesse (nbpart, xp, yp, op, up, vp, delta, time)
implicit none
integer,  intent(in)     :: nbpart      !< particles number
real(8), intent(in)     :: xp(nbpart)  !< particles x position
real(8), intent(in)     :: yp(nbpart)  !< particles y position
real(8), intent(in)     :: op(nbpart)  !< particles strength
real(8), intent(inout)  :: up(nbpart)  !< particles x velocity
real(8), intent(inout)  :: vp(nbpart)  !< particles y velocity
real(8), intent(in)     :: delta       !< smooth parameter
real(8) :: xo, yo, dx, dy, usum, vsum, dpi
integer  :: j, k
real(8) :: r2, r22, r2a1, r2a13, xm
real(8) :: a1, a12, a122
real(8) :: time

dpi   = 8.0 * atan( 1.0 )

a1    = delta
a12   = a1*a1
a122  = a12*a12

do k = 1, nbpart/2
   
   usum = 0.0; vsum = 0.0
   xo = xp(k); yo = yp(k)
  
   do j = 1 , nbpart/2
      if( j .ne. k ) then
         dx    = xp( j ) - xo
         dy    = yp( j ) - yo
         r2    = dx * dx + dy * dy
         r22   = r2 * r2
         r2a1  = r2 + a12
         r2a13 = r2a1 * r2a1 * r2a1
         xm    = (r22+3.0*a12*r2+4.0*a122) / r2a13
         usum  = usum + dy * op(j) * xm
         vsum  = vsum - dx * op(j) * xm
      end if
   end do

   up(k) = usum/dpi + gomeg * cos(gomeg*time)
   vp(k) = vsum/dpi - gomeg * sin(gomeg*time) 
   up(k+nbpart/2) = usum/dpi - gomeg * cos(gomeg*time)
   vp(k+nbpart/2) = vsum/dpi + gomeg * sin(gomeg*time)
   
end do

end subroutine vitesse

!> move particles
subroutine deplace (nbpart, xp, yp, up, vp, dt)
implicit none
integer, intent(in)  :: nbpart
real(8), intent(inout)  :: xp(nbpart) !< x particle position
real(8), intent(inout)  :: yp(nbpart) !< y particle position
real(8), intent(in)     :: up(nbpart) !< x particle velocity
real(8), intent(in)     :: vp(nbpart) !< y particle velocity
real(8), intent(in)     :: dt         !< time step
integer :: k

do k = 1, nbpart
   xp(k) = xp(k) + dt * up(k)
   yp(k) = yp(k) + dt * vp(k)
end do
   
end subroutine deplace

!> compute particles centers
subroutine centres(nbpart, xp, yp, op, time)
implicit none
integer, intent(in) :: nbpart
real(8), dimension(nbpart), intent(in) :: xp !< x position
real(8), dimension(nbpart), intent(in) :: yp !< y position
real(8), dimension(nbpart), intent(in) :: op !< particle weight
real(8) :: xc, yc, time

xc = sum(xp(1:nbpart/2)*op(1:nbpart/2))
yc = sum(yp(1:nbpart/2)*op(1:nbpart/2))
open(10, file="centres.dat",position="append")
if (time == 0.) rewind(10)
write(10,"(5f8.4)")time, xc/gam0, yc/gam0, &
      -sin(gomeg*time), cos(gomeg*time)
close(10)
end subroutine centres

!> compute real time
function getRealTimer()
implicit none
real(8) :: out, getRealTimer
sll_int64  :: count, count_rate
call system_clock(count, count_rate)
count = count - 1254348000*count_rate
out = count
out = out / count_rate
getRealTimer = out
end function getRealTimer

!> initialize particles positions
subroutine initialize( nstep, imov, xp, yp, op, delta, dt, nbpart ) 

namelist/donnees/nstep, dt, imov, amach, nray, r0, delta

integer, parameter :: nsec0 = 6
integer :: nstep, imov, idm, nbpart, nray, nsec, nr
real(8), dimension(:), pointer :: xp
real(8), dimension(:), pointer :: yp
real(8), dimension(:), pointer :: op

real(8), dimension(:), pointer :: rf
real(8), dimension(:), pointer :: zf
real(8), dimension(:), pointer :: gam

real(8) :: circ, al, ur, tau, aom, u0, r0
real(8) :: amach, delta, dt
integer  :: k
integer  :: error, file_id

nstep = 2000
dt    = 0.02
imov  = 1
amach = 0.1 !0.56
nray  = 20 !50
r0    = 0.5
delta = 0.01

!open(10, file = "input" )
!read(10,donnees)
!close(10)

u0 = amach 

gam0   = u0 * 2.0 * sll_pi / 0.7 * r0!gaussienne
!gam0   = 2. * sll_pi * r0 * u0	!constant
!gam0   = 2. * sll_pi / 10.0

aom    = gam0 / ( sll_pi * r0**2 )  ! Amplitude du vortex
tau    = 8.0 * sll_pi**2 / gam0     ! Periode de co-rotation
gomeg  = gam0/ (4.0*sll_pi)         ! Vitesse angulaire
ur     = gomeg                  ! Vitesse tangentielle du vortex
al     = 0.5 * tau              ! Longeur d'onde

call sll_new_file_id(file_id, error)
open(file_id, file="particles.out")
write(file_id,*) " iterations : ", nstep
write(file_id,*) " pas de temps : ", dt
write(file_id,*) " animation : ", imov, " steps "
write(file_id,*) " aom = ", aom
write(file_id,*) " r0 = ", r0
write(file_id,*) " Circulation = ", gam0
write(file_id,*) " Vitesse de rotation gomeg = ", gomeg
write(file_id,*) " ur = ", ur
write(file_id,*) " periode de corotation = ", tau
write(file_id,*) " --------------------------------------------- "

idm = 1
nsec = 0
do k = 1, nray
   nsec  = nsec + nsec0
   idm = idm + nsec
end do

write(*,*) "idm =", idm
SLL_ALLOCATE(rf(idm),error)
SLL_ALLOCATE(zf(idm),error)
SLL_ALLOCATE(gam(idm),error)
call distrib( rf, zf, gam, r0, idm, nray, nsec0, nr )
write(*,*) "nr =", nr

SLL_ASSERT(idm == nr)
nbpart = 2 * nr
SLL_ALLOCATE(xp(nbpart),error)
SLL_ALLOCATE(yp(nbpart),error)
SLL_ALLOCATE(op(nbpart),error)

do k = 1, nr
   xp( k    ) = rf( k ) 
   yp( k    ) = zf( k ) + 1.
   op( k    ) = gam( k )
   xp( k+nr ) = rf( k ) 
   yp( k+nr ) = zf( k ) - 1.
   op( k+nr ) = gam( k )
end do

!Calcul de la vitesse de propagation du systeme

circ = sum(op)

write(file_id,*) ' Nombre total de particules =',nbpart
write(file_id,*) ' Circulation totale  =',circ
close(file_id)

deallocate(rf,zf,gam)

end subroutine initialize

!---------------------------------------------------------------

!> set particles positions on a disc
subroutine distrib(rf, zf, cir, ray, idm, nray, nsec0, nr )

integer :: i, j, k
integer :: kd, nsec, nsec0, idm, nray, nr

real(8) :: rf( * ), zf( * ), cir( * ), ds( idm )
real(8) :: ssurf, q, sigma, teta, dss, r1, r2, s1, s2, eps
real(8) :: gamt, sgam, dteta, surf, ray, dray, r, dr

integer :: file_id, error


!     rf,zf : position de la particule
!     ds    : taille de l'element de surface
!     cir   : circulation de la particule
!     dr    : pas d'espace dans la direction radiale 
!     nray  : nb de pts ds la direction radiale.
!     dray  : rayon de la particule placee au centre
!     ray   : rayon de la section
!     gam0  : circulation totale 
!     surf  : surface de la section
!     nsec  : nombre de points dans la premiere couronne

dr      = ray / ( nray + 0.5 )
dray    = 0.5 * dr                !rayon de la section centrale
surf    = sll_pi * ray * ray 
dteta   = 2.0 * sll_pi / float( nsec0)

k       = 1
rf(  1) = 0.0
zf(  1) = 0.0
ds(  1) = sll_pi * dray * dray

if ( gauss ) then
   gamt = gam0 / ( 1. - exp( -1.0 ) )
   cir( 1) = gamt * ( 1.-exp(-(dray/ray)**2)) !gauss
else
   cir( 1) = gam0 * ds( 1 ) / surf   !uniforme
end if
sgam    = cir( 1 )

call sll_new_file_id(file_id, error)
open(file_id, file="particles_begin.out")
write(file_id,1000) rf(1), zf(1), cir(1), ds(1)

r1    = dray
s1    = sll_pi * r1**2
nsec  = 0

!cpn   *** parametre de l'ellipse ***
!c      eps = 0.01		! 0.0 --> disque
      eps = 0.0
!cpn   ******************************

do i = 1, nray

   nsec  = nsec + nsec0
   dteta = 2.0 * sll_pi / float(nsec)
   r     = float( i ) * dr 

   r2  = r + 0.5 * dr
   s2  = sll_pi * r2**2 
   dss = s2 - s1
   s1  = s2

   do j = 1, nsec

      k = k + 1
      if( k .gt. idm ) stop ' !! idm < Nr !! '

      teta    = float( j ) * dteta 
      sigma   = r * ( 1.0 + eps * cos( 2.0*teta ) )
      rf( k ) = sigma * cos( teta )
      zf( k ) = sigma * sin( teta )

      ds(  k ) = dss / float( nsec )

      if ( gauss ) then
         q       = 0.5 * (exp(-(r1/ray)**2)-exp(-(r2/ray)**2) )
         cir( k ) = gamt * dteta / sll_pi * q    ! gauss
      else
         cir( k ) = gam0 * ds( k ) / surf    ! uniforme
      end if

      sgam     = sgam + cir( k )

      write(file_id,1000) rf(k), zf(k), cir(k), ds(k)

    end do

    r1  = r2 

    kd = k - nsec + 1 
    write(file_id,1000) rf( kd ), zf( kd ), cir( kd ), ds(kd)
    write(file_id,1000) 

end do

close(file_id)

nr = k

ssurf = sum(ds(1:nr))

!write(*,*) 'Nb de pts sur la section :', nr,'(',idm,' max )'
!write(*,*) 'surface theorique - pratique :', (surf),' - ',(ssurf)
!if (gauss) then
!   write(*,*) 'circulation theorique - pratique :',(gam0),' ; ',(sgam)
!else
!   write(*,*) 'circulation theorique - pratique :',(gam0),' ; ',(sgam)
!end if

1000  format( F11.5, 2X, F11.5, 2X, F11.5, 1X, F11.5 )

end subroutine distrib

end module biot_savart
