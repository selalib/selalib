program vm4d_transpose
!-------------------------------------------------------------------
!  programme de simulation numerique d'un plasma electromagnetique 2D
!  modelise par les equations de Vlasov-Maxwell
!-------------------------------------------------------------------

#define MPI_MASTER 0

#include "selalib.h"

use geometry_module
use maxwell2dfdtd_module
use poisson2dpp_seq
use diagnostiques_module
use vlasov4d_plot
use vlasov2d_module
use splinepx_class
use splinepy_class
use vlasov1d_module
use init_functions
use sll_mudpack_cartesian

implicit none

type (geometry)      :: geomx      ! geometrie dans l'espace physique
type (geometry)      :: geomv      ! geometrie dans l'espace des vitesses
type (maxwell2dfdtd) :: maxw2dfdtd ! champ electromagnetique
type (poisson2dpp)   :: poiss2dpp  ! potentiel pour la correction
type (mudpack_2d)    :: poisson_mg
type (vlasov2d)      :: vlas2d     ! vlasov
type (splinepx)      :: splx       ! vlasov1d
type (splinepy)      :: sply       ! vlasov1d

sll_real64, dimension(:,:,:,:), pointer :: f,f1     ! fonc de distribution
sll_real64, dimension(:,:),     pointer :: ex,ey ! champ electrique
sll_real64, dimension(:,:),     pointer :: ex1,ey1 ! champ electrique
sll_real64, dimension(:,:),     pointer :: ex2,ey2 ! champ electrique
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
sll_real64 :: nrj,nrjex,nrjey,nrjbz
sll_real64 :: nrjtab(4)
sll_int32  :: error, iplot
sll_int32  :: comm, my_num, num_threads

sll_int32  :: va       = VA_VALIS
sll_int32  :: meth     = METH_CSL_PPM1
sll_int32  :: num_case = LANDAU_X_CASE
sll_int32  :: mud_case = SPECTRAL

sll_real64, allocatable, dimension(:,:) :: x1
sll_real64, allocatable, dimension(:,:) :: x2
sll_real64, allocatable, dimension(:,:) :: df
sll_real64 :: tcpu1, tcpu2

sll_real64, dimension(:,:),pointer :: phi, jxp

call sll_boot_collective()
num_threads  = sll_get_collective_size(sll_world_collective)
my_num       = sll_get_collective_rank(sll_world_collective)
comm         = sll_world_collective%comm

tcpu1 = MPI_WTIME()
if (my_num == MPI_MASTER) then
   print*,'MPI Version of slv2d running on ',num_threads, ' processors'
end if

call initglobal()

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

call initlocal( )


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

   if (va==VA_VALIS) then 
      div=rho;df=rho
      call transposexv(vlas2d,f)
      call densite_courant(vlas2d,jx,jy)
      jx1=jx;jy1=jy
      
      !compute Jx(x_{i+1/2},  y_j) and Jy(x_i, y_{j+1/2})
      call average_current(vlas2d,jx,jy,jx1,jy1)

      call transposevx(vlas2d,f)
      bz1=bz
      call solve_faraday(maxw2dfdtd,ex,ey,bz,dt)
      bz=0.5_8*(bz+bz1) !! bz = bz(n), bz1=bz(n-1/2)
      ex1=ex;ey1=ey    !! (ex1, ey1) = (ex, ey)(n) 
      jx1=0._8;jy1=0._8

      call c_l_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,0.5_8*dt)
      call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,0.5_8*dt)
      !(ex, ey)=(ex, ey)(n+1/2)
   endif

   call advection1d_x(splx,f,0.5_f64*dt,Jx1,meth)
   call advection1d_y(sply,f,0.5_f64*dt,Jy1,meth)

   ! transposition de la fonction de distribution
   call transposexv(vlas2d,f)

   if (va==VA_VLASOV_POISSON) then 
      call densite_charge(vlas2d,rho)

      if (mud_case==SPECTRAL) then 
         call solve(poiss2dpp,ex,ey,rho,nrj)
         call average(geomx,ex,ey)
      else
         call solve_poisson_mg(poisson_mg)
         !specific 1dx
         call poisson1d(vlas2d%geomx,rho,ex)
         ey=0._8;bz=0._f64;bz1=0._f64
      endif
         
   endif

   !compute (ex2, ey2) en i,j
   call average_int(vlas2d,ex,ey,ex2,ey2)

   call advection_v(vlas2d,ex2,ey2,dt,bz)

   call transposevx(vlas2d,f)

   call advection1d_y(sply,f,0.5_f64*dt,Jy2,meth)
   call advection1d_x(splx,f,0.5_f64*dt,Jx2,meth)

   Jx=0.5_8*(Jx1+Jx2)
   Jy=0.5_8*(Jy1+Jy2)

   if (va==VA_VALIS) then 
      call c_l_periodiques(maxw2dfdtd,ex1,ey1,bz,Jx,Jy,dt)
      call solve_ampere(maxw2dfdtd,ex1,ey1,bz,Jx,Jy,nrj,dt)
      ex=ex1;ey=ey1

      !verif charge conservation
      call transposexv(vlas2d,f)
      call densite_charge(vlas2d,rho)
      call transposevx(vlas2d,f)

      call verif_rho(vlas2d,div,df,rho,Jx,Jy)
      call verif_charge(vlas2d,div,ex,ey)

   else if (va==VA_VLASOV_POISSON) then 

      call transposexv(vlas2d,f)
      call densite_charge(vlas2d,rho)
      call transposevx(vlas2d,f)

      if (mud_case==SPECTRAL) then 
         call solve(poiss2dpp,ex,ey,rho,nrj)
      else
         call solve_poisson_mg(poisson_mg)
         !specific 1dx
         call poisson1d(geomx,rho,ex)
         ey=1.e-13
      endif

   endif
   
   !Diagnostic
   if (mod(iter,fdiag).eq.0) then 
      ! ecriture des resultats par le processeur 0
      call flush(6)
      call diagnostiquesm(f,rho,ex,ey,bz,jx,jy,geomx,geomv, &
                          jstartx,jendx,jstartv,jendv,iter/fdiag)
      call flush(6)
   endif

   if (mod(iter,fthdiag).eq.0) then 
      nrj=0.5_8*log(sum(ex*ex+ey*ey)*vlas2d%geomx%dx*vlas2d%geomx%dy)
      nrjex=sum(ex*ex)*vlas2d%geomx%dx*vlas2d%geomx%dy 
      nrjey=0.5_8*log(sum(ey*ey)*vlas2d%geomx%dx*vlas2d%geomx%dy)
      nrjbz=0.5_8*log(max(sum(bz*bz)*vlas2d%geomx%dx*vlas2d%geomx%dy,1.e-12))
      nrjtab(1)=nrjex;nrjtab(2)=nrjey;nrjtab(3)=nrjbz
      nrjtab(4)=0.5_8*sum(ex*ex+ey*ey+bz*bz)*vlas2d%geomx%dx*vlas2d%geomx%dy
      call thdiag_tab(vlas2d,f,nrjtab,iter*dt,jstartv)
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

subroutine initglobal()
!-------------------------------------------------------
! sert a l'initialisation globale du programme VP2D
!-------------------------------------------------------
sll_int32      :: nx, ny   ! dimensions de l'espace physique
sll_int32      :: nvx, nvy ! dimensions de l'espace des vitesses
sll_real64     :: x0, y0   ! coordonnees debut du maillage espace physique
sll_real64     :: vx0, vy0 ! coordonnees debut du maillage espace vitesses
sll_real64     :: x1, y1   ! coordonnees fin du maillage espace physique
sll_real64     :: vx1, vy1 ! coordonnees fin du maillage espace vitesses
sll_int32      :: iflag,ierr  ! indicateur d'erreur

! definition of namelists
namelist /time/ dt, nbiter
namelist /diag/ fdiag, fthdiag! freq. of diags and time hist diags in steps
namelist /phys_space/ x0,x1,y0,y1,nx,ny
namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
namelist /algo_charge/ va, meth
namelist /test_case/ num_case
namelist /solver_poisson/ mud_case

if (my_num == MPI_MASTER) then
   call fichinit
   read(idata,NML=time)
   read(idata,NML=diag)
   read(idata,NML=phys_space)
   read(idata,NML=vel_space)
   read(idata,NML=algo_charge)
   read(idata,NML=test_case)
   read(idata,NML=solver_poisson)
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
call mpi_bcast(va,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(meth,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(num_case,1,MPI_INTEGER,MPI_MASTER,comm,ierr)
call mpi_bcast(mud_case,1,MPI_INTEGER,MPI_MASTER,comm,ierr)

call initialize(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")

call initialize(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")


end subroutine initglobal

subroutine initlocal()

!variables locales
sll_int32 :: ipiece_size_v, ipiece_size_x
sll_real64 :: xi,vx,vy,v2,x,y,eps,kx,ky,nrj,vth,Tr,u
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
SLL_CLEAR_ALLOCATE(f(1:geomx%nx,1:geomx%ny,1:geomv%nx,jstartv:jendv),iflag)
SLL_CLEAR_ALLOCATE(f1(1:geomx%nx,1:geomx%ny,1:geomv%nx,jstartv:jendv),iflag)
SLL_CLEAR_ALLOCATE(rho(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ex(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ey(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(bz(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ex1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ey1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ex2(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(ey2(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(bz1(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(jx(1:geomx%nx,1:geomx%ny),iflag)
SLL_CLEAR_ALLOCATE(jy(1:geomx%nx,1:geomx%ny),iflag)

SLL_CLEAR_ALLOCATE(phi(1:geomx%nx+1,1:geomx%ny+1),iflag)
SLL_CLEAR_ALLOCATE(jxp(1:geomx%nx+1,1:geomx%ny+1),iflag)

! initialisation parallele des tableaux globaux, 
! ce qui permet  de les distribuer sur les processeurs
! initialisation de la fonction de distribution 
xi=0.90_8
eps=0.001_8
vth=0.02_8
Tr=12._8
u=0.1_8 !4.5_8
kx=2._8*sll_pi/(real(geomx%nx,8)*geomx%dx)
ky=2._8*sll_pi/(real(geomx%ny,8)*geomx%dy)
do jv=jstartv,jendv
   vy = geomv%y0+real(jv-1,8)*geomv%dy
   do iv=1,geomv%nx
      vx = geomv%x0+real(iv-1,8)*geomv%dx
      v2 = vx*vx+vy*vy
      do j=1,geomx%ny
         y=geomx%y0+real(j-1,8)*geomx%dy
         do i=1,geomx%nx
            x=geomx%x0+real(i-1,8)*geomx%dx

            select case(num_case)
                case(LANDAU_X_CASE)
                    f(i,j,iv,jv)= landau_1d(eps, kx, x, v2)
                case(LANDAU_Y_CASE)
                    f(i,j,iv,jv)= landau_1d(eps, ky, y, v2)
                case(LANDAU_COS_PROD_CASE)
                    f(i,j,iv,jv)= landau_cos_prod(eps, kx, ky, x, y, v2)
                case(LANDAU_COS_SUM_CASE)
                    f(i,j,iv,jv)= landau_cos_sum(eps, kx, ky, x, y, v2)
                case(TSI_CASE)
                    f(i,j,iv,jv)= tsi(eps, kx, x, vx, v2)
            end select
!BOT 1dx (to do)
!            f(i,j,iv,jv)=(1._wp+eps*cos(kx*x))*(0.9_8*exp(-0.5_8*vx*vx)+0.1_8*exp(-0.5_8*(vx-u)*(vx-u)))* &
!                 (1._wp/(2._wp*sll_pi))*exp(-0.5_wp*vy*vy)
!weibel (to do)
!            f(i,j,iv,jv)=(1._wp+eps*cos(kx*x)) &
!                 *(1._wp/(sll_pi*vth*sqrt(Tr)))*exp(-0.5_wp*(vx*vx+vy*vy/Tr)/vth)
         end do
      end do
   end do
end do

!Initialisation du module vlasov
call initialize(vlas2d,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)

!Intitialisation du champ electrique
call transposexv(vlas2d,f)
call densite_charge(vlas2d,rho)
call transposevx(vlas2d,f)
jx=rho

if (mud_case==SPECTRAL) then 
   call initialize(poiss2dpp,rho,geomx,iflag)
   call solve(poiss2dpp,ex,ey,rho,nrj)
   !compute ex(x_{i+1/2}, y_j), ey(x_i, y_{j+1/2})
   call average(geomx,ex,ey)
   print *,'MUD CASE',mud_case
else
   call initialize_mudpack_cartesian(poisson_mg, &
        geomx%x0, geomx%x1, geomx%nx, &
        geomx%y0, geomx%y1, geomx%ny, &
        SLL_PERIODIC, SLL_PERIODIC,   &
        SLL_PERIODIC, SLL_PERIODIC)
   !compute ex(x_{i+1/2}, y_j), ey(x_i, y_{j+1/2})
   call solve_poisson_mg(poisson_mg)
   print *,'MUD CASE',mud_case
   !specific 1dx
   call poisson1d(geomx,rho,ex)
   ey=0._8
   call verif_charge(vlas2d,rho,ex,ey)
endif

if (mud_case==SPECTRAL) then 
   call compute_rho_for_charge(vlas2d,rho,ex,ey)

   !compute f a partir de rho
   f=0._8
   do jv=jstartv,jendv
      vy = geomv%y0+real(jv-1,8)*geomv%dy
      do iv=1,geomv%nx
         vx = geomv%x0+real(iv-1,8)*geomv%dx
         v2 = vx*vx+vy*vy
         do j=1,geomx%ny
            y=geomx%y0+real(j-1,8)*geomx%dy
            do i=1,geomx%nx
               x=geomx%x0+real(i-1,8)*geomx%dx
               select case(num_case)
               case(TSI_CASE)
                  f(i,j,iv,jv)=rho(i,j)/(2._wp*sll_pi)*exp(-0.5_wp*v2)*vx*vx
               case default
                  f(i,j,iv,jv)=rho(i,j)/(2._wp*sll_pi)*exp(-0.5_wp*v2)
               end select
               !bot 1dx
               !            f(i,j,iv,jv)=rho(i,j)*(0.9_8*exp(-0.5_8*vx*vx)+0.1_8*exp(-(vx-u)*(vx-u)/2._8))* &
               !                 (1._wp/(2._wp*sll_pi))*exp(-0.5_wp*vy*vy)
               
               !TSI 1dx
               !            f(i,j,iv,jv)=rho(i,j)*(1._wp/(2._wp*sll_pi))*exp(-0.5_wp*v2)*vx*vx
            end do
         end do
      end do
   end do
   
endif

call transposexv(vlas2d,f)

call verif_charge(vlas2d,rho,ex,ey)

!weibel instability (to do)
do i=1,geomx%nx
   x=geomx%x0+real(i-0.5_8,8)*geomx%dx
   do j=1,geomx%ny
      bz(i,j)=0._8*eps*cos(kx*x)
   enddo
enddo


open(45+my_num,file="r_init:"//char(48+my_num))
open(46+my_num,file="e_init:"//char(48+my_num))
do j = jstartx, jendx
   do i = 1, geomx%nx
      x = geomx%x0+(i-1)*geomx%dx
      y = geomx%y0+(j-1)*geomx%dy
      write(45+my_num,*) sngl(x),sngl(y),sngl(rho(i,j))
      write(46+my_num,*) sngl(x),sngl(y),sngl(ex(i,j)),sngl(ey(i,j))
   end do
   write(45+my_num,*); write(16+my_num,*)
end do 
close(45+my_num)
close(46+my_num)

! initialisation du calcul du champ magnetique
call initialize(maxw2dfdtd,geomx,iflag, jstartx, jendx)
call initialize(splx,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)
call initialize(sply,geomx,geomv,iflag,jstartx,jendx,jstartv,jendv)

end subroutine initlocal

subroutine poisson1d(geomx,rho,ex)
  type(geometry) :: geomx
  sll_real64, dimension(:,:),pointer,intent(in)  :: rho
  sll_real64, dimension(:,:),pointer,intent(out) :: ex
  integer :: i,j
  sll_real64 :: mass

  do j=1,geomx%ny
     ex(1,j)=0._8
     mass=0._8
     do i=2,geomx%nx
        ex(i,j)=ex(i-1,j)+geomx%dx*(rho(i,j)-1._8)
        mass=mass+ex(i,j)
     enddo

     do i=1,geomx%nx
        ex(i,j)=ex(i,j)-mass/real(geomx%nx)
     enddo

  enddo

end subroutine poisson1d

subroutine average(geomx,ex,ey)
  type(geometry) :: geomx
  sll_real64, dimension(:,:),pointer :: ex,ey
  sll_real64, dimension(1:geomx%nx,1:geomx%ny) :: buf
  integer :: i,j

  buf=ex

  !compute the electric field ex at (x_{i+1/2}, y_{j})
  do i=1,geomx%nx-1
     do j=1,geomx%ny
        ex(i,j)=0.5_8*(buf(i+1,j)+buf(i,j)) 
     enddo
  enddo
  do j=1,geomx%ny
     ex(geomx%nx,j)=0.5_8*(buf(1,j)+buf(geomx%nx,j)) 
  enddo

  buf=ey
  !compute the electric field ey at (x_{i}, y_{j+1/2})
  
  do j=1,geomx%ny-1
     do i=1,geomx%nx
        ey(i,j)=0.5_8*(buf(i,j+1)+buf(i,j)) 
     enddo
  enddo
  do i=1,geomx%nx
     ey(i,geomx%ny)=0.5_8*(buf(i,1)+buf(i,geomx%ny)) 
  enddo

end subroutine average


subroutine average_current(vlas2d,jx,jy,jx1,jy1)
  type(vlasov2d) :: vlas2d
  sll_real64, dimension(:,:),pointer :: jx,jy,jx1,jy1
  integer :: i,j
  
  do j=1,vlas2d%geomx%ny
     do i=1,vlas2d%geomx%nx-1
        jx(i,j)=0.5_8*(jx1(i+1,j)+jx1(i,j))
     enddo
     jx(vlas2d%geomx%nx,j)=0.5_8*(jx1(1,j)+jx1(vlas2d%geomx%nx,j))
  enddo
  
  do i=1,vlas2d%geomx%nx
     do j=1,vlas2d%geomx%ny-1
        jy(i,j)=0.5_8*(jy1(i,j+1)+jy1(i,j))
     enddo
     jy(i,vlas2d%geomx%ny)=0.5_8*(jy1(i,1)+jy1(i,vlas2d%geomx%ny))
  enddo

end subroutine average_current


subroutine average_int(vlas2d,ex,ey,ex2,ey2)
  type(vlasov2d) :: vlas2d
  sll_real64, dimension(:,:),pointer :: ex,ey,ex2,ey2
  integer :: i,j
  
  
  ex2=0._8;ey2=0._8
  do j=1,vlas2d%geomx%ny
     do i=2,vlas2d%geomx%nx
        ex2(i,j)=0.5_8*(ex(i,j)+ex(i-1,j))
     enddo
     ex2(1,j)=0.5_8*(ex(1,j)+ex(vlas2d%geomx%nx,j))
  enddo
  
  do i=1,vlas2d%geomx%nx
     do j=2,vlas2d%geomx%ny
        ey2(i,j)=0.5_8*(ey(i,j)+ey(i,j-1))
     enddo
     ey2(i,1)=0.5_8*(ey(i,1)+ey(i,vlas2d%geomx%ny))
  enddo
  
end subroutine average_int


subroutine verif_charge(vlas2d,rho,ex,ey)
  type(vlasov2d) :: vlas2d
  sll_real64, dimension(:,:),pointer :: rho,ex,ey
  integer :: i,j
  sll_real64 :: mass

  mass=0._8
  do i=2,vlas2d%geomx%nx
     do j=2,vlas2d%geomx%ny
        mass=mass+abs((ex(i,j)-ex(i-1,j))/vlas2d%geomx%dx &
             +(ey(i,j)-ey(i,j-1))/vlas2d%geomx%dy-(rho(i,j)-1._8))
     enddo
  enddo
  !i=1, for all j=2,ny
!  do j=2,vlas2d%geomx%ny
!     mass=mass+abs((ex(1,j)-ex(vlas2d%geomx%nx,j))/vlas2d%geomx%dx &
!          +(ey(1,j)-ey(1,j-1))/vlas2d%geomx%dy-(rho(1,j)-1._8))
!  enddo
  !j=1, for all i=2,nx
!  do i=2,vlas2d%geomx%nx
!     mass=mass+abs((ex(i,1)-ex(i-1,1))/vlas2d%geomx%dx &
!          +(ey(i,1)-ey(i,vlas2d%geomx%ny))/vlas2d%geomx%dy-(rho(i,1)-1._8))
!  enddo
  !i=j=1
  mass=mass+abs((ex(1,1)-ex(vlas2d%geomx%nx,1))/vlas2d%geomx%dx &
       +(ey(1,1)-ey(1,vlas2d%geomx%ny))/vlas2d%geomx%dy-(rho(1,1)-1._8))


  if (mass>1.e-11) then 
     print *,'verif Poisson',mass*vlas2d%geomx%dx*vlas2d%geomx%dy
     print *,' '
  endif


end subroutine verif_charge



subroutine verif_rho(vlas2d,rho_pred,rho_tn,rho_ref,Jx,Jy)
  type(vlasov2d) :: vlas2d
  sll_real64, dimension(1:vlas2d%geomx%nx,1:vlas2d%geomx%ny), intent(in) :: rho_tn
  sll_real64, dimension(1:vlas2d%geomx%nx,1:vlas2d%geomx%ny),intent(out) :: rho_pred
  sll_real64, dimension(:,:),pointer :: rho_ref,Jx,Jy
  integer :: i,j
  
  do i=2,vlas2d%geomx%nx
     do j=2,vlas2d%geomx%ny
        rho_pred(i,j)=rho_tn(i,j)-(dt/vlas2d%geomx%dx)*(Jx(i,j)-Jx(i-1,j)) & 
             -(dt/vlas2d%geomx%dy)*(Jy(i,j)-Jy(i,j-1))
     enddo
  enddo
  do j=2,vlas2d%geomx%ny
     rho_pred(1,j)=rho_tn(1,j)-(dt/vlas2d%geomx%dx)*(Jx(1,j)-Jx(vlas2d%geomx%nx,j)) & 
          -(dt/vlas2d%geomx%dy)*(Jy(1,j)-Jy(1,j-1))
  enddo
  do i=2,vlas2d%geomx%nx
     rho_pred(i,1)=rho_tn(i,1)-(dt/vlas2d%geomx%dx)*(Jx(i,1)-Jx(i-1,1)) & 
          -(dt/vlas2d%geomx%dy)*(Jy(i,1)-Jy(i,vlas2d%geomx%ny))
  enddo
  rho_pred(1,1)=rho_tn(1,1)-(dt/vlas2d%geomx%dx)*(Jx(1,1)-Jx(vlas2d%geomx%nx,1)) & 
       -(dt/vlas2d%geomx%dy)*(Jy(1,1)-Jy(1,vlas2d%geomx%ny))
  
  if (sum(abs(rho_ref-rho_pred))>1.e-11) then 
     print *,'verif rho',sum(abs(rho_ref-rho_pred))*vlas2d%geomx%dx*vlas2d%geomx%dy
     print *,' '
  endif
  
end subroutine verif_rho

subroutine compute_rho_for_charge(vlas2d,rho,ex,ey)
  type(vlasov2d) :: vlas2d
  sll_real64, dimension(:,:),pointer :: rho,ex,ey
  integer :: i,j
  
  !compute rho such that discrete Poisson is ensured
  ![ex(x_{i+1/2}, y_j)-ex(x_{i-1/2}, y_j)]/dx + [ey(x_i, y_{j+1/2})-ey(x_i, y_{j-1/2})]/dy 
  != rho(x_i,y_j)-1
  
  rho=0._8
  do i=2,vlas2d%geomx%nx
     do j=2,vlas2d%geomx%ny
        rho(i,j)=(ex(i,j)-ex(i-1,j))/vlas2d%geomx%dx+(ey(i,j)-ey(i,j-1))/vlas2d%geomx%dy+1._8
     enddo
  enddo
  !i=1, for all j=2,ny
  do j=2,vlas2d%geomx%ny
     rho(1,j)=(ex(1,j)-ex(vlas2d%geomx%nx,j))/vlas2d%geomx%dx+(ey(1,j)-ey(1,j-1))/vlas2d%geomx%dy+1._8
  enddo
  do i=2,vlas2d%geomx%nx
     rho(i,1)=(ex(i,1)-ex(i-1,1))/vlas2d%geomx%dx+(ey(i,1)-ey(i,vlas2d%geomx%ny))/vlas2d%geomx%dy+1._8
  enddo
  rho(1,1)=(ex(1,1)-ex(vlas2d%geomx%nx,1))/vlas2d%geomx%dx+(ey(1,1)-ey(1,vlas2d%geomx%ny))/vlas2d%geomx%dy+1._8
  
end subroutine compute_rho_for_charge

subroutine solve_poisson_mg(this)

   type(mudpack_2d) :: this
   sll_real64 :: mass, avg

   mass=sum(rho(1:geomx%nx,1:geomx%ny))/(real(geomx%nx,8)*real(geomx%ny,8))

   jxp(1:geomx%nx,1:geomx%ny) = rho     -mass
   jxp(geomx%nx+1,1:geomx%ny) = rho(1,:)-mass
   jxp(1:geomx%nx,geomx%ny+1) = rho(:,1)-mass
   jxp(geomx%nx+1,geomx%ny+1) = rho(1,1)-mass
   
   call solve_mudpack_cartesian(this, phi, jxp)

   avg = sum(phi)
   phi = phi - avg / (geomx%nx*geomx%ny)
   
   do j = 1, geomx%ny
      do i = 1, geomx%nx
         ex(i,j) = (phi(i+1,j)-phi(i,j)) / geomx%dx
         ey(i,j) = (phi(i,j+1)-phi(i,j)) / geomx%dy
      end do
   end do
   
   nrj=0._wp
   do i=1,geomx%nx
      do j=1,geomx%ny
         nrj=nrj+ex(i,j)*ex(i,j)+ey(i,j)*ey(i,j)          
      enddo
   enddo
      
   nrj=nrj*geomx%dx*geomx%dy
   if (nrj>1.e-30) then 
      nrj=0.5_wp*log(nrj)
   else
      nrj=-10**9
   endif

end subroutine solve_poisson_mg

end program vm4d_transpose
