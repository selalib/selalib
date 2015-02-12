program vp4d_transpose

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

use used_precision  
use geometry_module
use diagnostiques_module
use vlasov4d_plot
use poisson2dpp_seq
use vlasov2d_module
use splinepx_class
use splinepy_class

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
sll_int32  :: jstartx, jendx, jstartv, jendv
sll_real64 :: nrj
sll_real64 :: tcpu1, tcpu2

sll_int32 :: my_num, num_threads, comm

call sll_boot_collective()

my_num = sll_get_collective_rank(sll_world_collective)
num_threads = sll_get_collective_size(sll_world_collective)
comm   = sll_world_collective%comm

! initialisation global
tcpu1 = MPI_WTIME()
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
   write(*,"(g13.3,1x,3i6)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
               f4d,rho,e_x,e_y,vlas2d,poisson,splx,sply)

call plot_mesh4d(geomx,geomv,jstartx,jendx,jstartv,jendv)
 
call advection_x(vlas2d,f4d,.5*dt)

do iter=1,nbiter

   call transposexv(vlas2d,f4d)

   call densite_charge(vlas2d,rho)

   call solve(poisson,e_x,e_y,rho,nrj)

   if (mod(iter,fthdiag).eq.0) then
       call thdiag(vlas2d,f4d,nrj,iter*dt,jstartv)    
   end if

   call advection_v(vlas2d,e_x,e_y,dt)

   call transposevx(vlas2d,f4d)

   if (mod(iter,fdiag) == 0) then 

       call advection_x(vlas2d,f4d,.5*dt)

       call diagnostiques(f4d,rho,e_x,e_y,geomx,geomv, &
                          jstartx,jendx,jstartv,jendv,iter/fdiag)

       call plot_df(f4d,iter/fdiag,geomx,geomv,jstartx,jendx,jstartv,jendv,VXVY_VIEW)


       call advection_x(vlas2d,f4d,.5*dt)

   else 

       call advection_x(vlas2d,f4d,dt)

   end if

end do

tcpu2 = MPI_WTIME()
if (my_num == MPI_MASTER) &
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

call sll_halt_collective()

print*,'PASSED'



contains

 subroutine initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)

  type(geometry)  :: geomx       ! geometrie globale du probleme
  type(geometry)  :: geomv       ! geometrie globale du probleme
  sll_real64      :: dt          ! pas de temps
  sll_int32       :: nbiter      ! nombre d'iterations en temps
  sll_int32       :: fdiag       ! frequences des diagnostiques
  sll_int32       :: fthdiag     ! frequences des historiques en temps
  sll_int32       :: nx, ny      ! dimensions de l'espace physique
  sll_int32       :: nvx, nvy    ! dimensions de l'espace des vitesses
  sll_real64      :: x0, y0      ! coordonnees debut du maillage espace physique
  sll_real64      :: vx0, vy0    ! coordonnees debut du maillage espace vitesses
  sll_real64      :: x1, y1      ! coordonnees fin du maillage espace physique
  sll_real64      :: vx1, vy1    ! coordonnees fin du maillage espace vitesses
  sll_int32       :: iflag,ierr  ! indicateur d'erreur
   
  namelist /time/ dt, nbiter
  namelist /diag/ fdiag, fthdiag
  namelist /phys_space/ x0,x1,y0,y1,nx,ny
  namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
  sll_int32 :: my_num, num_threads, comm

  my_num = sll_get_collective_rank(sll_world_collective)
  num_threads = sll_get_collective_size(sll_world_collective)
  comm   = sll_world_collective%comm
  
  if (my_num == MPI_MASTER) then
   
    call fichinit()
    read(idata,NML=time)
    read(idata,NML=diag)
    read(idata,NML=phys_space)
    read(idata,NML=vel_space)

  end if

  call mpi_bcast(dt,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(nbiter,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fdiag,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fthdiag, 1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(x0,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(y0,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(x1,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(y1,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(nx,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(ny,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(vx0,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(vy0,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(vx1,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(vy1,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(nvx,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(nvy,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)

  call initialize(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")
  call initialize(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")

 end subroutine initglobal

 subroutine initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
                      f,rho,e_x,e_y,vlas2d,poisson,splx,sply)

  type(geometry) :: geomx
  type(geometry) :: geomv  
  sll_int32      :: jstartv
  sll_int32      :: jendv
  sll_int32      :: jstartx
  sll_int32      :: jendx

  sll_real64, dimension(:,:)    ,pointer :: rho
  sll_real64, dimension(:,:)    ,pointer :: e_x
  sll_real64, dimension(:,:)    ,pointer :: e_y
  sll_real64, dimension(:,:,:,:),pointer :: f

  type(vlasov2d)    :: vlas2d
  type(splinepx)    :: splx
  type(splinepy)    :: sply
  type(poisson2dpp) :: poisson
  
  sll_int32  :: ipiece_size_v
  sll_int32  :: ipiece_size_x

  sll_real64 :: xi,vx,vy,v2,x,y,eps,kx,ky
  sll_int32  :: i,j,iv,jv,iflag
  sll_int32  :: my_num, num_threads, comm

  my_num      = sll_get_collective_rank(sll_world_collective)
  num_threads = sll_get_collective_size(sll_world_collective)
  comm        = sll_world_collective%comm

  ! initialisation of size of parallel zones 
  ! the total size of the vx zone is nvx
  ! the total size of the y zone is ny
  ! ipiece_size = n/num_threads rounded up
  ipiece_size_v = (geomv%ny + num_threads-1) / num_threads
  ipiece_size_x = (geomx%ny + num_threads-1) / num_threads
  ! zone a traiter en fonction du numero de process
  jstartv = my_num * ipiece_size_v + 1
  jendv   = min(jstartv - 1 + ipiece_size_v, geomv%ny)
  jstartx = my_num * ipiece_size_x + 1
  jendx   = min(jstartx - 1 + ipiece_size_x, geomx%ny)
    
  SLL_ASSERT(jstartx<jendx)
  SLL_ASSERT(jstartv<jendv)
  print*,'init zone ',my_num,jstartx,jendx,ipiece_size_x, &
                             jstartv,jendv,ipiece_size_v

  SLL_ALLOCATE(f(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)

  SLL_ALLOCATE(rho(geomx%nx,geomx%ny),iflag)
  SLL_ALLOCATE(e_x(geomx%nx,geomx%ny),iflag)
  SLL_ALLOCATE(e_y(geomx%nx,geomx%ny),iflag)

  xi  = 0.90
  eps = 0.05
  kx  = 2*pi/((geomx%nx)*geomx%dx)
  ky  = 2*pi/((geomx%ny)*geomx%dy)
  do jv=jstartv,jendv
     vy = geomv%y0+(jv-1)*geomv%dy
     do iv=1,geomv%nx
        vx = geomv%x0+(iv-1)*geomv%dx
        v2 = vx*vx+vy*vy
        do j=1,geomx%ny
           y=geomx%y0+(j-1)*geomx%dy
           do i=1,geomx%nx
              x=geomx%x0+(i-1)*geomx%dx
              f(i,j,iv,jv)=(1+eps*cos(kx*x))*1/(2*pi)*exp(-.5*v2)
           end do
        end do
     end do
  end do

  call new(poisson, rho, geomx, iflag)
  call new(vlas2d, geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(splx,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(sply,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)

 end subroutine initlocal

end program vp4d_transpose
