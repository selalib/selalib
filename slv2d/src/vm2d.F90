program VM2D
!-------------------------------------------------------------------
!  programme de simulation numerique d'un plasma electromagnetique 2D
!  modelise par les equations de Vlasov-Maxwell
!-------------------------------------------------------------------
#include "selalib.h"
use used_precision  
use geometry_module
use maxwell2dfdtd_module
use poisson2dpp_seq
use diagnostiques_module
use vlasov2d_module
use splinepx_class
use splinepy_class
use vlasov1d_module
use vm2dinit


implicit none

type (geometry)      :: geomx      ! geometrie dans l'espace physique
type (geometry)      :: geomv      ! geometrie dans l'espace des vitesses
type (maxwell2dfdtd) :: maxw2dfdtd ! champ electromagnetique
type (poisson2dpp)   :: poiss2dpp  ! potentiel pour la correction
type (vlasov2d)      :: vlas2d     ! vlasov
type (splinepx)      :: splx       ! vlasov1d
type (splinepy)      :: sply       ! vlasov1d


sll_real64, dimension(:,:,:,:), pointer :: f,f1     ! fonc de distribution
sll_real64, dimension(:,:),     pointer :: ex,ey ! champ electrique
sll_real64, dimension(:,:),     pointer :: ex1,ey1 ! champ electrique
sll_real64, dimension(:,:),     pointer :: jx,jy ! courant
sll_real64, dimension(:,:),     pointer :: bz,bz1    ! champ magnetique
sll_real64, dimension(:,:),     pointer :: rho   ! charge
sll_real64, dimension(:,:),     pointer :: Jx1,Jx2,Jy1,Jy2   ! courant partiel 
sll_real64, dimension(:,:),     pointer :: div   ! divergence E

! donnees du probleme
sll_int32      :: nbiter         ! nombre d'iterations en temps
sll_real64     :: dt             ! pas de temps
sll_int32      :: fdiag, fthdiag ! frequences des diagnostiques
sll_int32      :: iflag          ! indicateur d'erreur
sll_int32      :: iter,i,j       ! variables de boucles       

sll_int32  :: jstartx, jendx, jstartv, jendv
sll_real64 :: nrj
sll_int32       :: error

sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df
sll_real64 :: tcpu1, tcpu2

! initialisation global

call initialise_moduleMPI
tcpu1 = MPI_WTIME()
if (my_num == MPI_MASTER) then
   print*,'MPI Version of slv2d running on ',num_threads, ' processors'
end if

call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)
  
if (my_num == MPI_MASTER) then
   ! write some run data
   write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
   write(*,"(2(i3,1x),6(g13.3,1x))")&
   geomx%nx, geomx%ny, geomx%x0,&
   geomx%x0+(geomx%nx-1)*geomx%dx,&
   geomx%y0, geomx%y0+(geomx%ny-1)*geomx%dy,&
   geomx%dx, geomx%dy   
   write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
   write(*,"(2(i3,1x),6(g13.3,1x))")&
   geomv%nx, geomv%ny, geomv%x0,&
   geomv%x0+(geomv%nx-1)*geomv%dx,&
   geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy,&
   geomv%dx, geomv%dy
   write(*,*) 'dt,nbiter,fdiag,fthdiag'
   write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
endif

call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
               f,f1,rho,ex,ey,ex1,ey1,bz,bz1,jx,jy, &
               vlas2d,maxw2dfdtd,poiss2dpp,splx,sply)

! ecriture des resultats par le processeur 0 a l'instant initial
call mpi_barrier(MPI_COMM_WORLD,iflag)

iter = 0
call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv,&
jstartx,jendx,jstartv,jendv,iter)

call mpi_barrier(MPI_COMM_WORLD,iflag)

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

   call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
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
   call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)
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

   write(*,*) iter
   if (mod(iter,fthdiag).eq.0) then 
      call thdiag(vlas2d,f,nrj,iter*dt)
   endif

end do
tcpu2 = MPI_WTIME()
if (my_num == MPI_MASTER) &
   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

print*,'PASSED'
call termine_moduleMPI

contains

subroutine plot_solution( f )

   sll_real64, dimension(:,:,:,:), intent(in) :: f
   sll_int32 :: file_id
   sll_int32, save :: iplot = 0
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

end program VM2D
