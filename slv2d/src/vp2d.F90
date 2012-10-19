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

sll_int32  :: i, j, k, l
sll_int32  :: nbiter  
sll_real64 :: dt     
sll_int32  :: fdiag, fthdiag  
sll_int32  :: iter 
sll_int32  :: sx, ex, sy, ey, su, eu, sv, ev
sll_real64 :: nrj
sll_real64 :: tcpu1, tcpu2

sll_int32 :: my_num, num_threads, comm, error

sll_real64, dimension(:,:),pointer :: fxu
sll_real64, dimension(:,:),pointer :: fyv
sll_real64 :: sumloc

character(len=4)  :: prefix = "df4d"
integer(HID_T)    :: file_id
integer(HSSIZE_T) :: offset(2), offset_1d(1)
integer(HSIZE_T)  :: global_dims(2), global_dim(1)
character(len=4)  :: counter

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
   write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,sv,ev,sx,ex, &
               f4d,rho,e_x,e_y,vlas2d,poisson,splx,sply)

sy = 1
ey = geomx%ny
su = 1
eu = geomv%nx

SLL_ALLOCATE(fxu(sx:ex,su:eu),error)
SLL_ALLOCATE(fyv(sy:ey,sv:ev),error)

offset_1d = 0
call sll_hdf5_file_create("mesh4d.h5",file_id,error)
global_dim = geomx%nx
call sll_hdf5_write_array(file_id,global_dim,offset_1d,geomx%xgrid,"/x",error)
global_dim = geomx%ny
call sll_hdf5_write_array(file_id,global_dim,offset_1d,geomx%ygrid,"/y",error)
global_dim = geomv%nx
call sll_hdf5_write_array(file_id,global_dim,offset_1d,geomv%xgrid,"/u",error)
global_dim = geomv%ny
call sll_hdf5_write_array(file_id,global_dim,offset_1d,geomv%ygrid,"/v",error)
call sll_hdf5_file_close(file_id, error)

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

       offset = 0

       call int2string(iter/fdiag,counter)
       
       do k=su,eu
        do i=sx,ex
         sumloc= sum(f4d(i,sy:ey,k,sv:ev))
         call  mpi_reduce(sumloc,fxu(i,k),1,MPI_REAL8,MPI_SUM,0,comm,error)
        end do
       end do
       do l=sv,ev
        do j=sy,ey
         fyv(j,l)= sum(f4d(sx:ex,j,su:eu,l))
        end do
       end do

       call sll_hdf5_file_create(prefix//counter//".h5",file_id,error)
       global_dims = (/ex-sx+1,eu-su+1/)
       offset      = (/0, 0/)
       call sll_hdf5_write_array(file_id,global_dims,offset,fxu,"/fxvx",error)
       global_dims = (/ey-sy+1,ev-sv+1/)
       offset      = (/0, sv-1/)
       call sll_hdf5_write_array(file_id,global_dims,offset,fyv,"/fyvy",error)
       call sll_hdf5_file_close(file_id, error)

       call diagnostiques(f4d,rho,e_x,e_y,geomx,geomv,sx,ex,sv,ev,iter/fdiag)

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

call sll_halt_collective()

print*,'PASSED'

end program VP2D
