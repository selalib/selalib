program vm4d_transpose
!-------------------------------------------------------------------
!  programme de simulation numerique d'un plasma electromagnetique 2D
!  modelise par les equations de Vlasov-Maxwell
!-------------------------------------------------------------------

#define MPI_MASTER 0

#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use geometry_module
use maxwell2dfdtd_module
use poisson2dpp_seq
use diagnostiques_module
use vlasov4d_plot
use vlasov2d_module
use splinepx_class
use splinepy_class
use vlasov1d_module

implicit none

type(geometry)      :: geomx      ! geometrie dans l'espace physique
type(geometry)      :: geomv      ! geometrie dans l'espace des vitesses
type(maxwell2dfdtd) :: maxw2dfdtd ! champ electromagnetique
type(poisson2dpp)   :: poiss2dpp  ! potentiel pour la correction
type(vlasov2d)      :: vlas2d     ! vlasov
type(splinepx)      :: splx       ! vlasov1d
type(splinepy)      :: sply       ! vlasov1d

sll_real64, dimension(:,:,:,:), pointer :: f,f1     ! fonc de distribution
sll_real64, dimension(:,:),     pointer :: ex,ey ! champ electrique
sll_real64, dimension(:,:),     pointer :: ex1,ey1 ! champ electrique
sll_real64, dimension(:,:),     pointer :: jx,jy ! courant
sll_real64, dimension(:,:),     pointer :: bz,bz1    ! champ magnetique
sll_real64, dimension(:,:),     pointer :: rho   ! charge
sll_real64, dimension(:,:),     pointer :: Jx1,Jx2,Jy1,Jy2   ! courant partiel 
sll_real64, dimension(:,:),     pointer :: div   ! divergence E

! donnees du probleme
sll_int32  :: nbiter         ! nombre d'iterations en temps
sll_real64 :: dt             ! pas de temps
sll_int32  :: fdiag, fthdiag ! frequences des diagnostiques
sll_int32  :: iter,i,j       ! variables de boucles       

sll_int32  :: jstartx, jendx, jstartv, jendv
sll_real64 :: nrj
sll_int32  :: error, iplot
sll_int32  :: comm, my_num, num_threads

sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df
sll_real64 :: tcpu1, tcpu2

call sll_boot_collective()
num_threads  = sll_get_collective_size(sll_world_collective)
my_num       = sll_get_collective_rank(sll_world_collective)
comm         = sll_world_collective%comm

tcpu1 = MPI_WTIME()
if (my_num == MPI_MASTER) then
   print*,'MPI Version of slv2d running on ',num_threads, ' processors'
end if

call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  
if (my_num == MPI_MASTER) then
   ! write some run data
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i5,1x),6(g15.3,1x))")&
   geomx%nx, geomx%ny, geomx%x0,&
   geomx%x0+(geomx%nx-1)*geomx%dx,&
   geomx%y0, geomx%y0+(geomx%ny-1)*geomx%dy,&
   geomx%dx, geomx%dy   
   write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
   write(*,"(2(i6,1x),6(g15.3,1x))")&
   geomv%nx, geomv%ny, geomv%x0,&
   geomv%x0+(geomv%nx-1)*geomv%dx,&
   geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy,&
   geomv%dx, geomv%dy
   write(*,*) 'dt,nbiter,fdiag,fthdiag'
   write(*,"(g15.3,1x,3i6)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
               f,f1,rho,ex,ey,ex1,ey1,bz,bz1,jx,jy, &
               vlas2d,maxw2dfdtd,poiss2dpp,splx,sply)

! ecriture des resultats par le processeur 0 a l'instant initial
call mpi_barrier(comm,error)

iter = 0
call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv,&
jstartx,jendx,jstartv,jendv,iter)

call mpi_barrier(comm,error)

call transposevx(vlas2d,f)

SLL_ALLOCATE(div(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(Jx1(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(Jx2(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(Jy1(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(Jy2(geomx%nx,geomx%ny),error)

SLL_ALLOCATE(x1(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(x2(geomx%nx,geomx%ny),error)
SLL_ALLOCATE(df(geomx%nx,geomx%ny),error)

do j = 1, geomx%ny
   do i = 1, geomx%ny
      x1(i,j) = geomx%xgrid(i)
      x2(i,j) = geomx%ygrid(j)
      df(i,j) = sum(f(i,j,:,:))
   end do
end do 

do iter=1,nbiter

   ! advection d'un demi-pas de temps en espace
   !call advection_x(vlas2d,f,.5*dt)
   call advection1d_x(splx,f,0.5_f64*dt,Jx1)
   call advection1d_y(sply,f,0.5_f64*dt,Jy1)

   !Recopie de f^(n,2) dans f1
   f1=f;ex1=ex;ey1=ey;bz1=bz

   ! transposition de la fonction de distribution
   call transposexv(vlas2d,f)

   !###########################
   !Phase prediction a la VALIS
   !########################### 

   ! calcul de la densite de courant et de charge
   ! Les deux operations doivent etre faites dans cet ordre
   ! Chercher la raison ...
   call densite_charge(vlas2d,rho)
   call densite_courant(vlas2d,jx,jy)
     
   ! calcul du champ magnetique Bz(k+1/2)
   if (iter == 1) then
      call solve_faraday(maxw2dfdtd,ex,ey,bz,0.5_f64*dt)
   else
      call solve_faraday(maxw2dfdtd,ex,ey,bz,dt)
   end if

   call c_l_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   !call silver_muller(maxw2dfdtd,ex,ey,bz,jx,jy,dt)

   call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)

   ! advection d'un pas de temps en vitesse
   call advection_v(vlas2d,ex,ey,dt,bz)
!   call advection_v(vlas2d,ex,ey,dt)

   ! transposition de la fonction de distribution
   call transposevx(vlas2d,f)

   ! advection d'un demi-pas de temps en espace     
!      call advection_x(vlas2d,f,0.5_f64*dt)
   call advection1d_x(splx,f,0.5_f64*dt,Jx2)
   call advection1d_y(sply,f,0.5_f64*dt,Jy2)

   f=f1;ex=ex1;ey=ey1;

   !################
   !Phase correction
   !################
   call c_l_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   !call silver_muller(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
   call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)

   ! transposition de la fonction de distribution
   call transposexv(vlas2d,f)

   ! advection d'un pas de temps en vitesse
   call advection_v(vlas2d,ex,ey,dt,bz)

   ! transposition de la fonction de distribution
   call transposevx(vlas2d,f)

   ! advection d'un demi-pas de temps en espace     
!      call advection_x(vlas2d,f,0.5_f64*dt)
   call advection1d_x(splx,f,0.5_f64*dt,Jx1)
   call advection1d_y(sply,f,0.5_f64*dt,Jy1)

   !Diagnostic
   if (mod(iter,fdiag).eq.0) then 
      ! ecriture des resultats par le processeur 0
      call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv, &
                          jstartx,jendx,jstartv,jendv,iter/fdiag)
   endif

   if (mod(iter,fthdiag).eq.0) then 
      call thdiag(vlas2d,f,nrj,iter*dt,jstartv)
   endif

end do
tcpu2 = MPI_WTIME()
if (my_num == MPI_MASTER) &
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

print*,'PASSED'
call sll_halt_collective()

contains

subroutine plot_solution( )

   use sll_xdmf
   sll_int32 :: file_id
   character(len=4) :: cplot

   do j = 1, geomx%ny
      do i = 1, geomx%nx
         df(i,j) = sum(f(i,j,:,:))
      end do
   end do 

   iplot = iplot+1
   call int2string(iplot,cplot)

   call sll_xdmf_open("f"//cplot//".xmf","mesh",geomx%nx,geomx%ny,file_id,error)
   call sll_xdmf_write_array("mesh",x1,'x',error)
   call sll_xdmf_write_array("mesh",x2,'y',error)
   call sll_xdmf_write_array("fxy",df,"NodeVal",error,file_id,"Node")
   call sll_xdmf_close(file_id,error)

end subroutine plot_solution

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
sll_int32      :: iflag,ierr  ! indicateur d'erreur
sll_int32      :: comm
sll_int32 :: my_num, num_threads

! definition of namelists
namelist /time/ dt, nbiter
namelist /diag/ fdiag, fthdiag! freq. of diags and time hist diags in steps
namelist /phys_space/ x0,x1,y0,y1,nx,ny
namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy

num_threads  = sll_get_collective_size(sll_world_collective)
my_num = sll_get_collective_rank(sll_world_collective)
comm   = sll_world_collective%comm
   
if (my_num == MPI_MASTER) then
   call fichinit
   read(idata,NML=time)
   read(idata,NML=diag)
   read(idata,NML=phys_space)
   read(idata,NML=vel_space)
end if

call mpi_barrier(comm,ierr)
call mpi_bcast(dt,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(nbiter,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(fdiag,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(fthdiag,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(x0,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(y0,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(x1,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(y1,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(nx,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(ny,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(vx0,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(vy0,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(vx1,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(vy1,1,MPI_REAL8,MPI_MASTER,comm,ierr)
call mpi_bcast(nvx,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(nvy,1,MPI_INTEGER,MPI_MASTER,comm,ierr)

geomx = geometry(x0,  y0,  x1,  y1,  nx,  ny,  "perxy")
geomv = geometry(vx0, vy0, vx1, vy1, nvx, nvy, "natxy")

end subroutine initglobal

subroutine initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
                     f,f1,rho,ex,ey,ex1,ey1,bz,bz1,jx,jy,vlas2d,maxw2dfdtd,  &
                     poiss2dpp,splx,sply)

use splinepx_class
use splinepy_class

!-------------------------------------------------------
! sert a l'initialisation en parallele du programme VP2D
!-------------------------------------------------------
type(geometry) :: geomx, geomv  ! geometrie globale du probleme
sll_int32 :: jstartv,jendv,jstartx,jendx ! definition de la bande du proc
sll_real64, dimension(:,:,:,:),pointer :: f,f1
sll_real64, dimension(:,:),pointer :: rho,ex,ey,ex1,ey1,bz,bz1,jx,jy
type(vlasov2d) :: vlas2d
type(splinepx) :: splx
type(splinepy) :: sply
type(maxwell2dfdtd) :: maxw2dfdtd
type(poisson2dpp) :: poiss2dpp
!variables locales
sll_int32 :: ipiece_size_v, ipiece_size_x
sll_real64 :: xi,vx,vy,v2,x,y,eps,kx,ky,nrj
sll_int32 :: i,j,iv,jv,iflag, comm
sll_int32 :: my_num, num_threads

my_num = sll_get_collective_rank(sll_world_collective)
num_threads = sll_get_collective_size(sll_world_collective)

comm   = sll_world_collective%comm

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
SLL_CLEAR_ALLOCATE(rho(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ex(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ey(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(bz(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ex1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ey1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(bz1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(jx(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(jy(1:geomx%nx,1:geomx%ny),iflag)

! initialisation parallele des tableaux globaux, 
! ce qui permet  de les distribuer sur les processeurs
! initialisation de la fonction de distribution 
xi=0.90
eps=0.05
kx=2*sll_pi/((geomx%nx)*geomx%dx)
ky=2*sll_pi/((geomx%ny)*geomx%dy)
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
            f(i,j,iv,jv)=(1._wp+eps*cos(kx*x))*(1._wp/(2._wp*sll_pi))*exp(-0.5_wp*v2)
         end do
      end do
   end do
end do

!Initialisation du module poisson
call new(poiss2dpp,rho,geomx,iflag)
!Initialisation du module vlasov
call new(vlas2d,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
!Intitialisation du champ electrique
call transposexv(vlas2d,f)
call densite_charge(vlas2d,rho)
call solve(poiss2dpp,ex,ey,rho,nrj)
!call solve(poiss2dpp,ex,rho,nrj)

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

end program vm4d_transpose
