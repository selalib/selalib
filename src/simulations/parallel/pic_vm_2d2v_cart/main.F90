program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use initialisation
use maxwell
use diagno

implicit none

type(tm_mesh_fields) :: f0
type(tm_mesh_fields) :: f1
type(particle)       :: p

sll_real64 :: time
sll_real64 :: xmin
sll_real64 :: xmax
sll_real64 :: ymin
sll_real64 :: ymax

sll_int32  :: istep
sll_int32  :: iplot
sll_int32  :: iargc
sll_int32  :: n
sll_int32  :: i
sll_int32  :: error

character(len=272) :: argv

n = iargc()
if (n == 0) stop 'Usage: ./bin/test_pic2d fichier-de-donnees.nml'
do i = 1, n
  call getarg( i, argv)
  write(*,'(i2, 1x, a)') i, argv
end do

call readin( trim(argv) )

SLL_ALLOCATE(f0%ex(0:nx-1,0:ny),  error)
SLL_ALLOCATE(f0%ey(0:nx,0:ny-1),  error)
SLL_ALLOCATE(f0%bz(0:nx-1,0:ny-1),error)
SLL_ALLOCATE(f0%jx(0:nx-1,0:ny),  error)
SLL_ALLOCATE(f0%jy(0:nx,0:ny-1),  error)
SLL_ALLOCATE(f0%r0(0:nx,0:ny),    error) !rho au temps n
SLL_ALLOCATE(f0%r1(0:nx,0:ny),    error) !rho au temps n+1
SLL_ALLOCATE(f1%ex(0:nx,0:ny),    error) !decales sur maillage de Maxwell
SLL_ALLOCATE(f1%ey(0:nx,0:ny),    error)
SLL_ALLOCATE(f1%bz(0:nx,0:ny),    error)

time  = 0.d0
iplot = 0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

call init( f0 )                 !initialisation des champs et densites

xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p, time ) !creation des particules

do istep = 1, nstep

  if (istep > 1) then
    call faraday( f0 )     !Calcul de B(n-1/2) --> B(n)			
  end if

  call decalage( f0, f1 )
  call interpol_eb( f1, p )

  call avancee_vitesse( p )

  call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
  call sortie_part( p )
  call calcul_j_cic( p, f0 )
  call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
  call sortie_part( p )
       
  call calcul_rho( p, f0 )

  call faraday( f0 )   !Calcul de B(n) --> B(n+1/2)
  call ampere( f0 )    !Calcul de E(n) --> E(n+1)
  call conditions_limites( f0, time )

  time = time + dt

  !call plot_particles_center( p )

  if ( istep==1 .or. mod(istep,idiag) == 0 .or. istep==nstep ) then
    iplot = iplot + 1
    ! call plot_particles_center( p, time)  
    ! call plot_particle_density( p, iplot)  
    ! call plot_particles_points3d(p, iplot)
    ! call plot_particles_xmdv( p, iplot, xmin, xmax, ymin, ymax)  
    ! call diag_coc( f0, p, time, iplot )
    ! call diag_champ_part( p, time, iplot )
    ! call plot_champ( f0, iplot, time )
    ! call plot_phases( p, iplot, time )
    ! call distribution_v( p, iplot, time )  
    ! call distribution_x( p, iplot, time )
  endif

  call modeE( f0, iplot, time )
  write(*,"('istep = ', i6, ' time = ')", advance='no') istep

end do

print*,'PASSED'

end program test_pic2d
