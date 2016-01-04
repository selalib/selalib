program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use initialisation
use poisson
use villasenor
use maxwell
use diagno

implicit none

type(tm_mesh_fields) :: f0, f1
type(particle) :: p

sll_real64 :: time, xmin, xmax, ymin, ymax
integer :: istep, iplot
integer :: iargc, n, i
character(len=72) :: argv
sll_int32 :: error

n = iargc()
if (n == 0) stop 'Usage: ./bin/test_pic2d fichier-de-donnees.nml'
do i = 1, n
   call getarg( i, argv); !write(*,'(i2, 1x, a)') i, argv
end do

call readin( trim(argv) )

SLL_ALLOCATE(f0%ex(0:nx-1,0:ny),error)
SLL_ALLOCATE(f0%ey(0:nx,0:ny-1),error)
SLL_ALLOCATE(f0%bz(0:nx-1,0:ny-1),error)
SLL_ALLOCATE(f0%jx(0:nx-1,0:ny),error)
SLL_ALLOCATE(f0%jy(0:nx,0:ny-1),error)
SLL_ALLOCATE(f0%r0(0:nx,0:ny),error)  !rho au temps n
SLL_ALLOCATE(f0%r1(0:nx,0:ny),error)  !rho au temps n+1
SLL_ALLOCATE(f1%ex(0:nx,0:ny),error) !decales sur maillage de Maxwell
SLL_ALLOCATE(f1%ey(0:nx,0:ny),error)
SLL_ALLOCATE(f1%bz(0:nx,0:ny),error)

time  = 0.d0
iplot = 0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

call init( f0 )                 !initialisation des champs et densites

xmin = 0.0; xmax = dimx
ymin = 0.0; ymax = dimy

do istep = 1, nstep

   if ((nomcas == "faisce") .or. istep == 1) then
      call creapa( p, time ) !creation des particules
   endif
   if (istep > 1) then
      call faraday( f0 )     !Calcul de B(n-1/2) --> B(n)			
   end if

   call decalage( f0, f1 )
   call interpol_eb( f1, p )

   call avancee_vitesse( p )

   if (jname == 'jcico1') then
      call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
      call sortie_part( p )
      call calcul_j_cic( p, f0 )
      call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
      call sortie_part( p )
   else if (jname == 'jcoco1') then
      call avancee_part( p, 1.d0 )
      call calcul_j_villa( p, f0 )
      call sortie_part( p )
   else
      call avancee_part( p, 1.d0 )
      call sortie_part( p )
   end if
        
   !call calcul_rho( p, f0 )

   call faraday( f0 )   !Calcul de B(n) --> B(n+1/2)
   call ampere( f0 )    !Calcul de E(n) --> E(n+1)
   call conditions_limites( f0, time )

   time = time + dt

   !call plot_particles_center( p )

   if ( istep==1 .or. mod(istep,idiag) == 0 .or. istep==nstep ) then
      iplot = iplot + 1
      !call plot_particles_center( p, time)  
      !call plot_particle_density( p, iplot)  
      call plot_particles_points3d(p, iplot)
      call plot_particles_xmdv( p, iplot, xmin, xmax, ymin, ymax)  
      !if (nomcas=='viry__') call plot_part( p, time, iplot )
     ! call diag_coc( f0, p, time, iplot )
     ! call diag_champ_part( p, time, iplot )
     !call plot_champ( f0, iplot, time )
     ! call plot_phases( p, iplot, time )
     ! call distribution_v( p, iplot, time )  
     ! call distribution_x( p, iplot, time )
      !if (nomcas == 'plasma') call modeE( f0, iplot, time )
   endif

end do

print*,'PASSED'

end program test_pic2d
