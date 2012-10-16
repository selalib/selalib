program VP2D
  !-------------------------------------------------------------------
  !  programme de simulation numerique d'un plasma electrostatique 2D
  !  modelise par les equations de Vlasov-Poisson
  !-------------------------------------------------------------------
#include "selalib.h"
use used_precision  
use geometry_module
use diagnostiques_module
use poisson2dpp_module
use vlasov2d_module
use vp2dinit

implicit none

type (geometry)    :: geomx ! geometrie dans l'espace physique
type (geometry)    :: geomv ! geometrie dans l'espace des vitesses
type (poisson2dpp) :: poiss2dpp !champ electrique
type (vlasov2d)    :: vlas2d ! vlasov
type (splinepx)    :: splx       ! vlasov1d
type (splinepy)    :: sply       ! vlasov1d

sll_real64, dimension(:,:,:,:), pointer :: f ! fonc de distribution
sll_real64, dimension(:,:),     pointer :: rho  ! densite de charge
sll_real64, dimension(:,:),     pointer :: ex,ey ! champ electrique

! donnees du probleme
sll_int32      :: nbiter   ! nombre d'iterations en temps
sll_real64     :: dt       ! pas de temps
sll_int32      :: fdiag, fthdiag    ! frequences des diagnostiques

sll_int32      :: iter ! variables de boucles       

sll_int32  :: jstartx, jendx, jstartv, jendv
sll_real64 :: nrj
sll_real64 :: tcpu1, tcpu2

! initialisation global
  tcpu1 = MPI_WTIME()
  call initialise_moduleMPI
  if (my_num.eq.0) then
     print*,'MPI Version of slv2d running on ',num_threads, ' processors'
  end if

call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  
if (my_num == MPI_MASTER) then
   ! write some run data
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i3,1x),6(g13.3,1x))") geomx%nx, geomx%ny, geomx%x0, &
                                     geomx%x0+(geomx%nx)*geomx%dx, &
                                     geomx%y0, geomx%y0+(geomx%ny)*geomx%dy, geomx%dx, geomx%dy   
   write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
   write(*,"(2(i3,1x),6(g13.3,1x))") geomv%nx, geomv%ny, geomv%x0, &
          geomv%x0+(geomv%nx-1)*geomv%dx, &
          geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy, geomv%dx, geomv%dy
   write(*,*) 'dt,nbiter,fdiag,fthdiag'
   write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
               f,rho,ex,ey,vlas2d,poiss2dpp,splx,sply)

iter = 0
call diagnostiques(f,rho,ex,ey,geomx,geomv,jstartx,jendx,jstartv,jendv,iter)
 
call advection_x(vlas2d,f,.5*dt)

do iter=1,nbiter

   call transposexv(vlas2d,f)

   call densite_charge(vlas2d,rho)

   call solve(poiss2dpp,ex,ey,rho,nrj)

   call advection_v(vlas2d,ex,ey,dt)

   call transposevx(vlas2d,f)

   if (mod(iter,fdiag) == 0) then 

       call advection_x(vlas2d,f,.5*dt)

       call diagnostiques(f,rho,ex,ey,geomx,geomv, &
              jstartx,jendx,jstartv,jendv,iter/fdiag)

       if (mod(iter,fthdiag).eq.0) then
          call thdiag(vlas2d,f,nrj,iter*dt)    
       end if

       call advection_x(vlas2d,f,.5*dt)

   else 
       call advection_x(vlas2d,f,dt)
   end if

   tcpu2 = MPI_WTIME()
   if (my_num == MPI_MASTER) &
      write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

end do

call termine_moduleMPI

print*,'PASSED'

end program VP2D
