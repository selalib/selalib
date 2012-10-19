module vp2dinit

#include "selalib.h"

 use mpi
 use geometry_module
 use vlasov2d_module
 use splinepx_class
 use splinepy_class
 use poisson2dpp_seq
 use diagnostiques_module

 implicit none

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

  call new_geometry2(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")

  call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")

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

  call new(poisson, rho,   geomx, iflag)
  call new(vlas2d,  geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(splx,    geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(sply,    geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)

 end subroutine initlocal

end module vp2dinit
