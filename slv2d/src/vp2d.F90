!-------------------------------------------------------------------
!  programme de simulation numerique d'un plasma electrostatique 2D
!  modelise par les equations de Vlasov-Poisson
!-------------------------------------------------------------------
program VP2D

#include "selalib.h"
use used_precision  
use geometry_module
use diagnostiques_module
use poisson2dpp_seq
use vlasov2d_module
use vp2dinit

implicit none

type (geometry)    :: geomx 
type (geometry)    :: geomv 
type (poisson2dpp) :: poisson 
type (vlasov2d)    :: vlas2d 
type (splinepx)    :: splx 
type (splinepy)    :: sply

sll_real64, dimension(:,:,:,:), pointer :: f4d
sll_real64, dimension(:,:),     pointer :: rho
sll_real64, dimension(:,:),     pointer :: e_x
sll_real64, dimension(:,:),     pointer :: e_y 

sll_int32  :: nbiter  
sll_real64 :: dt     
sll_int32  :: fdiag, fthdiag  
sll_int32  :: iter 
sll_int32  :: sx, ex, sv, ev
sll_real64 :: nrj
sll_real64 :: tcpu1, tcpu2

! initialisation global
tcpu1 = MPI_WTIME()
call initialise_moduleMPI
if (my_num == MPI_MASTER) then
   print*,'MPI Version of slv2d running on ',num_threads, ' processors'
end if

call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  
if (my_num == MPI_MASTER) then
   ! write some run data
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i3,1x),6(g13.3,1x))") geomx%nx, geomx%ny, geomx%x0, &
                                     geomx%x0+(geomx%nx)*geomx%dx, &
                                     geomx%y0, geomx%y0+(geomx%ny)*geomx%dy, &
                                     geomx%dx, geomx%dy   
   write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
   write(*,"(2(i3,1x),6(g13.3,1x))") geomv%nx, geomv%ny, geomv%x0, &
                                     geomv%x0+(geomv%nx-1)*geomv%dx, &
                                     geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy, &
                                     geomv%dx, geomv%dy
   write(*,*) 'dt,nbiter,fdiag,fthdiag'
   write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,sv,ev,sx,ex, &
               f4d,rho,e_x,e_y,vlas2d,poisson,splx,sply)

iter = 0
call diagnostiques(f4d,rho,e_x,e_y,geomx,geomv,sx,ex,sv,ev,iter)
 
call advection_x(vlas2d,f4d,.5*dt)

do iter=1,nbiter

   call transposexv(vlas2d,f4d)

   call densite_charge(vlas2d,rho)

   call solve(poisson,e_x,e_y,rho,nrj)

   call advection_v(vlas2d,e_x,e_y,dt)

   call transposevx(vlas2d,f4d)

   if (mod(iter,fdiag) == 0) then 

       call advection_x(vlas2d,f4d,.5*dt)

       call diagnostiques(f4d,rho,e_x,e_y,geomx,geomv, &
              sx,ex,sv,ev,iter/fdiag)

       if (mod(iter,fthdiag).eq.0) then
          call thdiag(vlas2d,f4d,nrj,iter*dt)    
       end if

       call advection_x(vlas2d,f4d,.5*dt)

   else 

       call advection_x(vlas2d,f4d,dt)

   end if


end do

tcpu2 = MPI_WTIME()
if (my_num == MPI_MASTER) &
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

call termine_moduleMPI

print*,'PASSED'

end program VP2D
