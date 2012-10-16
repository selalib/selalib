module vm2dinit
#include "selalib.h"
use geometry_module
use vlasov2d_module
use splinepx_class
use splinepy_class
use diagnostiques_module
use poisson2dpp_module
use maxwell2dfdtd_module

implicit none

contains

subroutine initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
!-------------------------------------------------------
! sert a l'initialisation globale du programme VP2D
!-------------------------------------------------------
type(geometry)  :: geomx, geomv  ! geometrie globale du probleme
sll_real64     :: dt       ! pas de temps
sll_int32      :: nbiter   ! nombre d'iterations en temps
sll_int32      :: fdiag    ! frequences des diagnostiques
sll_int32      :: fthdiag    ! frequences des historiques en temps
sll_int32      :: nx, ny   ! dimensions de l'espace physique
sll_int32      :: nvx, nvy ! dimensions de l'espace des vitesses
sll_real64     :: x0, y0   ! coordonnees debut du maillage espace physique
sll_real64     :: vx0, vy0 ! coordonnees debut du maillage espace vitesses
sll_real64     :: x1, y1   ! coordonnees fin du maillage espace physique
sll_real64     :: vx1, vy1 ! coordonnees fin du maillage espace vitesses
sll_int32 :: iflag,ierr  ! indicateur d'erreur
! definition of namelists
namelist /time/ dt, nbiter
namelist /diag/ fdiag, fthdiag! freq. of diags and time hist diags in steps
namelist /phys_space/ x0,x1,y0,y1,nx,ny
namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
   
if (my_num == MPI_MASTER) then
   call fichinit
   read(idata,NML=time)
   read(idata,NML=diag)
   read(idata,NML=phys_space)
   read(idata,NML=vel_space)
end if

call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_bcast(dt,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(nbiter,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(fdiag,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(fthdiag,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(x0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(y0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(x1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(y1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(nx,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(ny,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(vx0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(vy0,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(vx1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(vy1,1,mpi_realtype,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(nvx,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)
call mpi_bcast(nvy,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,ierr)

call new(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")

call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")

end subroutine initglobal

subroutine initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
                     f,f1,rho,ex,ey,ex1,ey1,bz,bz1,jx,jy,vlas2d,maxw2dfdtd,  &
		     poiss2dpp,splx,sply)
!-------------------------------------------------------
! sert a l'initialisation en parallele du programme VP2D
!-------------------------------------------------------
type(geometry) :: geomx, geomv  ! geometrie globale du probleme
sll_int32 :: jstartv,jendv,jstartx,jendx ! definition de la bande du proc
sll_real64, dimension(:,:,:,:),pointer :: f,f1
sll_real64, dimension(:,:),pointer :: rho,ex,ey,ex1,ey1,bz,bz1,jx,jy
sll_int32 :: proc, mpierror
type(vlasov2d) :: vlas2d
type(splinepx) :: splx
type(splinepy) :: sply
type(maxwell2dfdtd) :: maxw2dfdtd
type(poisson2dpp) :: poiss2dpp
!variables locales
sll_int32 :: ipiece_size_v, ipiece_size_x
sll_real64 :: xi,vx,vy,v2,x,y,eps,kx,ky,nrj
sll_int32 :: i,j,iv,jv,iflag

! cas sequentiel
jstartv=1
jendv=geomv%nx
jstartx=1
jendx=geomx%ny

! initialisation of size of parallel zones 
! the total size of the vx zone is nvx
! the total size of the y zone is ny
! ipiece_size = n/num_threads rounded up
ipiece_size_v = (geomv%ny + num_threads - 1) / num_threads
ipiece_size_x = (geomx%ny + num_threads - 1) / num_threads
! zone a traiter en fonction du numero de process
jstartv = my_num * ipiece_size_v + 1
jendv = min(jstartv - 1 + ipiece_size_v, geomv%ny)
jstartx = my_num * ipiece_size_x + 1
jendx = min(jstartx - 1 + ipiece_size_x, geomx%ny)
    
SLL_ASSERT(jstartv < jendv) 
SLL_ASSERT(jstartx < jendx) 

print*,'init zone ',my_num,jstartx,jendx,ipiece_size_x, &
		           jstartv,jendv,ipiece_size_v

! allocation dynamique des tableaux
SLL_ALLOCATE(f(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)
SLL_ALLOCATE(f1(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)
!!$  SLL_ALLOCATE(rho(geomx%nx,jstartx:jendx),iflag)
!!$  if (iflag.ne.0) stop 'erreur dans l allocation de rho'
!!$  SLL_ALLOCATE(ex(geomx%nx,jstartx:jendx),iflag)
!!$  if (iflag.ne.0) stop 'erreur dans l allocation de ex'
!!$  SLL_ALLOCATE(ey(geomx%nx,jstartx:jendx),iflag)
!!$  if (iflag.ne.0) stop 'erreur dans l allocation de ey'
! Poisson n'est pas parallele pour l'instant
SLL_ALLOCATE(rho(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(ex(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(ey(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(bz(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(ex1(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(ey1(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(bz1(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(jx(geomx%nx,geomx%ny),iflag)
SLL_ALLOCATE(jy(geomx%nx,geomx%ny),iflag)

! initialisation parallele des tableaux globaux, 
! ce qui permet  de les distribuer sur les processeurs
! initialisation de la fonction de distribution 
xi=0.90
eps=0.05
kx=2*pi/((geomx%nx)*geomx%dx)
ky=2*pi/((geomx%ny)*geomx%dy)
do jv=jstartv,jendv
   vy = geomv%y0+(jv-1)*geomv%dy
   do iv=1,geomv%nx
      vx = geomv%x0+(iv-1)*geomv%dx
      v2 = vx*vx+vy*vy
      do j=1,geomx%ny
         y=geomx%y0+(j-1)*geomx%dy
         do i=1,geomx%nx
            x=geomx%x0+(i-1)*geomx%dx
!            f(i,j,iv,jv)=(1+eps*((cos(2*kx*x)+cos(3*kx*x))/1.2 &
!                 + cos(kx*x)))* &
!                 1/(2*pi)*((2-2*xi)/(3-2*xi))* &
!                 (1+.5*vx*vx/(1-xi))*exp(-.5*v2)
!             f(i,j,iv,jv)=(1+eps*cos(kx*x)*cos(ky*y))*1/(2*pi)*exp(-.5*v2)
            f(i,j,iv,jv)=(1._wp+eps*cos(kx*x))*(1._wp/(2._wp*pi))*exp(-0.5_wp*v2)
         end do
      end do
   end do
end do

! initialisation de ex,ey, bz, jx, jy

ex(:,:) = 0.0; ey(:,:)  = 0.0
ex1(:,:)= 0.0; ey1(:,:) = 0.0
jx(:,:) = 0.0; jy(:,:)  = 0.0 
bz(:,:) = 0.0; bz1(:,:) = 0.0; 
rho(:,:) = 0.0

!Initialisation du module poisson
call new(poiss2dpp,rho,geomx,iflag,jstartx,jendx)
!Initialisation du module vlasov
call new(vlas2d,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
!Intitialisation du champ electrique
call transposexv(vlas2d,f)		!Transposition
call densite_charge(vlas2d,rho)		!calcul de rho
call solve(poiss2dpp,ex,ey,rho,nrj)	!calcul de ex et ey
!call solve(poiss2dpp,ex,rho,nrj)	!calcul de ex (cas 1d)

open(15+my_num,file="r_init:"//char(48+my_num))
open(16+my_num,file="e_init:"//char(48+my_num))
do j = jstartx, jendx
   do i = 1, geomx%nx
      x = geomx%x0+(i-1)*geomx%dx
      y = geomx%y0+(j-1)*geomx%dy
      write(15+my_num,*) sngl(x),sngl(y),sngl(rho(i,j))
      write(16+my_num,*) sngl(x),sngl(y),sngl(ex(i,j)),sngl(ey(i,j))
   end do
   write(15+my_num,*); write(16+my_num,*)
end do 
close(15+my_num)
close(16+my_num)

! initialisation du calcul du champ magnetique
call new(maxw2dfdtd,geomx,iflag, jstartx, jendx)
call new(splx,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
call new(sply,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)

end subroutine initlocal

end module vm2dinit
