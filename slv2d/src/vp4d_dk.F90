!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!> @author
!> Michel Mehrenberger
!> Nicolas Crouseilles
!> @brief 
!> first drift kinetic simulation in slv2d
!> from the basis of vp4d_transpose
!> and CG_polar
!> should be obsolete and replaced/translated in simulation_class

program vp4d_dk

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"
!#include "sll_poisson_solvers.h"

use used_precision  
use geometry_module
use diagnostiques_module
use vlasov4d_plot
use poisson2dpp_seq
use vlasov2d_dk_module
use splinepx_class
use splinepy_class
use polar_operators
use polar_advection

implicit none

type (geometry)    :: geomx 
type (geometry)    :: geomv 
type (poisson2dpp) :: poisson 
type (vlasov2d)    :: vlas2d 
type (splinepx)    :: splx 
type (splinepy)    :: sply
type(sll_plan_poisson_polar), pointer :: plan_poisson
type(sll_SL_polar), pointer :: plan_sl
sll_real64, dimension(:,:,:,:), pointer :: f4d,f4d_old,f_equil
sll_real64, dimension(:,:),     pointer :: rho
sll_real64, dimension(:,:,:),     pointer :: rho_dk
sll_real64, dimension(:,:,:),     pointer :: phi
sll_real64, dimension(:,:,:,:),     pointer :: adv_field
sll_real64, dimension(:,:),     pointer :: e_x
sll_real64, dimension(:,:),     pointer :: e_y 
sll_real64, dimension(:,:),     pointer :: profile 

sll_int32  :: i,j,k,nbiter,scheme_case=2,qns_case=1
!scheme_case=1 EULER WITH 2DT
!scheme_case=2 !PREDICTOR CORRECTOR
!qns_case=0 !NO QUASI NEUTRAL
!qns_case=1 !NO ZONAL FLOW
!qns_case=2 !ZONAL FLOW
sll_real64 :: dt
sll_real64 :: R0     
sll_int32  :: fdiag, fthdiag,adv_case(4) 
sll_int32  :: iter 
sll_int32  :: jstartx, jendx, jstartv, jendv
sll_real64 :: nrj,tmp,tmp_tab(2)
sll_real64 :: tcpu1, tcpu2

sll_int32 :: my_num, num_threads, comm
sll_int32 :: nr,ntheta,err!,bc(2)
sll_real64 :: rmin,dr,r,y,z

sll_int32 :: kmin(2),kmax(2)

sll_real64 ::rpeak,deltarn,kappan,deltarTi,eps!,R0

sll_int32 :: ierr

adv_case(1) = 1   !1: normal 2: delta-f
adv_case(2) = 1   !1: BSL2D, 2: CSL2D
adv_case(3:4) = (/4,8/)   ! for CSL2D (3,0): for ppm0; (3,1) for ppm1; (3,2) for ppm2 (4,d): for FD(2d+1)



call sll_boot_collective()

my_num = sll_get_collective_rank(sll_world_collective)
num_threads = sll_get_collective_size(sll_world_collective)
comm   = sll_world_collective%comm


!kmin=(/9,4/)
!kmax=(/9,4/)

kmin=(/5,1/)
kmax=(/5,6/)



! initialisation global
tcpu1 = MPI_WTIME()
if (my_num == MPI_MASTER) then
   print*,'#MPI Version of slv2d running on ',num_threads, ' processors'
end if

!call initglobal_dk(geomx,geomv,dt,nbiter,fdiag,fthdiag,R0)



call init_dk(geomx,geomv,nbiter,fdiag,fthdiag,jstartv,jendv,jstartx,jendx, &
               f4d,f4d_old,f_equil,rho_dk,phi,adv_field,e_x,e_y,profile,vlas2d,plan_poisson,plan_sl,splx,sply)

   iter=0
   if (mod(iter,fdiag) == 0) then
     !print *,'#param are',iter/fdiag,geomx%x0,geomv%x0,jstartx,jendx,jstartv,jendv,VXVY_VIEW
     if(iter/fdiag==0)then
       call plot_mesh4d(geomx, geomv, jstartx,jendx,jstartv,jendv)
     endif
     call plot_df(f4d-f_equil,iter/fdiag,geomx,geomv,jstartx,jendx,jstartv,jendv,XY_VIEW)
   endif




dt=plan_sl%adv%dt
if (my_num == MPI_MASTER) then
  print *,'#solution with rsync'
  print *,'#Nr=',geomx%nx-1
  print *,'#Ntheta=',geomx%ny
  print *,'#Nphi=',geomv%nx
  print *,'#Nvpar=',geomv%ny
  print *,'#rmin=',geomx%x0
  print *,'#rmax=',geomx%x1
  print *,'#vmin=',geomv%y0
  print *,'#vmax=',geomv%y1
  print *,'#thetamin=',geomx%y0
  print *,'#thetamax=',geomx%y1
  print *,'#phimin=',geomv%x0
  print *,'#phimax=',geomv%x1
  print *,'#dt=',dt
  print *,'#nbiter=',nbiter
  print *,'#fdiag=',fdiag
  print *,'#fthdiag=',fthdiag  
endif
  

iter=0
   call transposexv(vlas2d,f4d)
   
   call densite_charge_dk(vlas2d,rho_dk)
!   do i=1,geomx%nx
!     tmp = sum(rho_dk(i,1:geomx%ny,1:geomv%nx))/real(geomx%ny*geomv%nx,f64)
!     rho_dk(i,:,:) = (rho_dk(i,:,:)-tmp)/profile(1,i)
!   enddo  
   do i=1,geomx%nx
     do j=1,geomx%ny+1
       !tmp = sum(rho_dk(i,j,1:geomv%nx))/real(geomv%nx,f64)
       !rho_dk(i,j,:) = (rho_dk(i,j,:)-tmp)/profile(1,i)
       
       !rho_dk(i,j,:) = (rho_dk(i,j,:))/profile(1,i)-1._f64
       
       !rho_dk(i,j,:) = (rho_dk(i,j,:))-1._f64
     enddo
     !print *,i,rho_dk(i,1,1),profile(1,i)
   enddo  
   
   !print *,maxval(abs(rho_dk))
   
   !stop
   
   
   !check first the poisson solver
   
   rpeak=7.2500725_f64!6.525_f64
   deltarn=2.899971_f64!2.03_f64
   kappan=5.51729655227586810E-002_f64!133.98_f64
   R0=1._f64!239.25_f64
   deltarTi=1.4499855_f64!1.015_f64
   eps=1.e-2_f64
   !deltarn=deltarn/deltarTi
   
!   do i=1,geomx%nx
!     r = geomx%x0+real(i-1,f64)*geomx%dx
!     tmp = -kappan/R0*cosh((r-rpeak)/deltarn)**(-2)
!     do j=1,geomx%ny+1
!       y=geomx%y0+real(j-1,f64)*geomx%dy
!       
!       phi(i,j,:)=cos(kmin(1)*y)*exp(-(r-rpeak)**2/(deltarn/deltarTi))*&
!       &(1._f64/profile(3,i)+2._f64/(deltarn/deltarTi)-4._f64*(r-rpeak)**2/(deltarn/deltarTi)**2&
!       &+2._f64*(r-rpeak)/(deltarn/deltarTi)*(1._f64/r+tmp)+kmin(1)**2/r**2)
!       !&-rho_dk(i,j,:)
!       
!       !phi(i,j,:)=exp(-(r-rpeak)**2/(deltarn/deltarTi))
!     enddo
!     print *,r,maxval((phi(i,:,:))),minval((phi(i,:,:))),&
!     &maxval((rho_dk(i,:,:))), minval((rho_dk(i,:,:)))
!   enddo  
!   
!   !stop
!   

   call sll_gnuplot_write(profile(1,:),'n0_r_init_slv2d',ierr)
   call sll_gnuplot_write(profile(2,:),'Ti_r_init_slv2d',ierr)
   call sll_gnuplot_write(profile(3,:),'Te_r_init_slv2d',ierr)
   
   call sll_gnuplot_write(rho_dk(:,1,1)/profile(1,:)-1._f64,'rho_slv2d_0',ierr)
 
   call solve_quasi_neutral(vlas2d,plan_poisson,rho_dk,phi,profile,qns_case)
   
   print *,'#delta=',deltarn/deltarTi
   !stop
   !compute error
   do i=1,geomx%nx
     r = geomx%x0+real(i-1,f64)*geomx%dx
     phi(i,1:geomx%ny,geomv%nx+1)=phi(i,1:geomx%ny,1)
     phi(i,geomx%ny+1,:)=phi(i,1,:)     
     do j=1,geomx%ny+1
       y=geomx%y0+real(j-1,f64)*geomx%dy
       do k=1,geomv%nx+1
         z=geomv%x0+real(k-1,f64)*geomv%dx
         rho_dk(i,j,k)=phi(i,j,k)-eps*exp(-(r-rpeak)**2/(deltarn/deltarTi))*&
         &cos(kmin(1)*y+kmin(2)*2._f64*sll_pi/(geomv%x1-geomv%x0)*z)
       enddo  
     enddo
     !print *,r,maxval((phi(i,:,:))),minval((phi(i,:,:)))!,minval((phi(i,1:geomx%ny,1:geomv%nx))),exp(-(r-rpeak)**2/deltarn)
   enddo
   
   print *,"#error for phi initial",maxval(abs(rho_dk))
   
   !stop
   
   call transposevx(vlas2d,f4d)
   call thdiag(vlas2d,f4d,phi,real(iter,f64)*dt,jstartv,kmin,kmax)
   

do iter=1,nbiter


   f4d_old=f4d
   call transposexv(vlas2d,f4d)
   call compute_field_dk(vlas2d,plan_sl%grad,phi,adv_field)

   ! prediction step
   call advection_x3_dk(vlas2d,0.5_f64*dt)
   call advection_x4_dk(vlas2d,adv_field(3,:,:,:),0.5_f64*dt)
   call transposevx(vlas2d,f4d)
   plan_sl%adv%dt = 0.5_f64*dt   
   call advection_x_dk(vlas2d,plan_sl%adv,f4d,adv_field,adv_case)
   call transposexv(vlas2d,f4d)
   
   
   !compute field at time t_{n+1/2}
   call densite_charge_dk(vlas2d,rho_dk)
   call solve_quasi_neutral(vlas2d,plan_poisson,rho_dk,phi,profile,qns_case)

   if(scheme_case==2)then

   call compute_field_dk(vlas2d,plan_sl%grad,phi,adv_field)
   
   !advection from tn to t_{n+1}
   call transposexv(vlas2d,f4d_old)
   call advection_x3_dk(vlas2d,0.5_f64*dt)
   call advection_x4_dk(vlas2d,adv_field(3,:,:,:),0.5_f64*dt)
   call transposevx(vlas2d,f4d)
   plan_sl%adv%dt = dt   
   call advection_x_dk(vlas2d,plan_sl%adv,f4d,adv_field,adv_case)
   call transposexv(vlas2d,f4d)
   call advection_x4_dk(vlas2d,adv_field(3,:,:,:),0.5_f64*dt)
   call advection_x3_dk(vlas2d,0.5_f64*dt)
   

   call densite_charge_dk(vlas2d,rho_dk)
   call solve_quasi_neutral(vlas2d,plan_poisson,rho_dk,phi,profile,qns_case)



   endif
   
   call transposevx(vlas2d,f4d)

   call thdiag(vlas2d,f4d,phi,real(iter,f64)*dt,jstartv,kmin,kmax)
!if (my_num==MPI_MASTER) then
!   print *,my_num,iter,dt,real(iter,f64)*dt
!end if

   if (mod(iter,fdiag) == 0) then
     !print *,'#param are',iter/fdiag,geomx%x0,geomv%x0,jstartx,jendx,jstartv,jendv,VXVY_VIEW
     if(iter/fdiag==1)then
       call plot_mesh4d(geomx, geomv, jstartx,jendx,jstartv,jendv)
     endif
     call plot_df(f4d-f_equil,iter/fdiag,geomx,geomv,jstartx,jendx,jstartv,jendv,XY_VIEW)
   endif

enddo   


   
   if (my_num == MPI_MASTER) then
   
     !write(20,*),maxval(abs(phi(:,:,:))),minval((phi(:,:,:)))
     do i=1,geomx%nx
       write (40,*),geomx%x0+real(i-1,f64)*geomx%dx,&
       &sum(phi(i,1:geomx%ny,1:geomv%nx))*geomx%dy*geomv%dx,&
       &sum(rho_dk(i,1:geomx%ny,1:geomv%nx))*geomx%dy*geomv%dx/(geomx%y1*geomv%x1),&
       &sum(rho_dk(i,1:geomx%ny,1:geomv%nx))/(geomx%ny*geomv%nx),&
       profile(1,i),profile(2,i),profile(3,i)
     enddo 
   
   endif

!tcpu2 = MPI_WTIME()
!if (my_num == MPI_MASTER) &
!   write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

call sll_halt_collective()

print*,'#PASSED'



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

  call new(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")
  call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")

 end subroutine initglobal

 subroutine initglobal_dk(geomx,geomv,dt,nbiter,fdiag,fthdiag,R0)

  type(geometry)  :: geomx       ! geometrie globale du probleme
  type(geometry)  :: geomv       ! geometrie globale du probleme
  sll_int32       :: iflag,ierr  ! indicateur d'erreur
  
  !MESH
  sll_int32       :: Nr,Ntheta,Nphi,Nvpar,nbiter,fdiag,fthdiag
  sll_real64      :: a,rhomin,rhomax,Lz,aspect_ratio,vmax,dt,R0
  
  sll_int32 :: my_num, num_threads, comm
  namelist /MESH/ Nr,Ntheta,Nphi,Nvpar,a,rhomin,rhomax,Lz,aspect_ratio,vmax,nbiter,dt,&
           &fdiag,fthdiag    
!  &MESH
!  Nr           = 127
!  Ntheta       = 128
!  Nphi         = 32
!  Nvpar        = 47
!  a            = 32.
!  rhomin       = 0.2
!  rhomax       = 0.8
!  Ltheta       = 6.283185307179586476925286766559005768394
!  Lz         = 6.283185307179586476925286766559005768394
!  aspect_ratio = 3.
!  vmax      = 6.
!  nbiter      = 30
!  dt      = 16.
!  fdiag      = 30
!  fthdiag      = 30

  my_num = sll_get_collective_rank(sll_world_collective)
  num_threads = sll_get_collective_size(sll_world_collective)
  comm   = sll_world_collective%comm
  
  if (my_num == MPI_MASTER) then
    call fichinit()
    read(idata,NML=MESH)
  end if

  !MESH  
  call mpi_bcast(Nr,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Ntheta,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Nphi,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Nvpar,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(a,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(rhomin,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(rhomax,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(Lz,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(aspect_ratio,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(vmax,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(nbiter,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(dt,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(fdiag,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fthdiag,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)

  

  call new(geomx,a*rhomin,0._f64,a*rhomax,2._f64*sll_pi,Nr+1,Ntheta,iflag,"pery")  
  !call new(geomv,0._f64,-vmax,Lz,vmax,Nphi,Nvpar,iflag,"perxy")
  call new(geomv,0._f64,-vmax,Lz,vmax,Nphi,Nvpar,iflag,"perx")
  R0 = a*aspect_ratio
  
 end subroutine initglobal_dk



 subroutine init_dk(geomx,geomv,nbiter,fdiag,fthdiag,jstartv,jendv,jstartx,jendx, &
                      f,f_old,f_equil,rho,phi,adv_field,e_x,e_y,profile,vlas2d,plan_poisson,plan_sl,splx,sply)

  type(geometry) :: geomx
  type(geometry) :: geomv
  sll_int32      :: jstartv,jendv
  sll_int32      :: jstartx,jendx

  sll_real64, dimension(:,:,:)    ,pointer :: rho
  sll_real64, dimension(:,:,:)    ,pointer :: phi
  sll_real64, dimension(:,:,:,:)    ,pointer :: adv_field
  sll_real64, dimension(:,:)    ,pointer :: e_x
  sll_real64, dimension(:,:)    ,pointer :: e_y
  sll_real64, dimension(:,:)    ,pointer :: profile
  sll_real64, dimension(:,:,:,:),pointer :: f,f_old,f_equil
  type(vlasov2d)    :: vlas2d
  type(splinepx)    :: splx
  type(splinepy)    :: sply
  !type(poisson2dpp) :: poisson
  type(sll_plan_poisson_polar), pointer :: plan_poisson
  type(sll_SL_polar), pointer :: plan_sl
  sll_int32  :: ipiece_size_v
  sll_int32  :: ipiece_size_x

  sll_real64 :: xi,vx,vy,v2,x,y,eps,ky,kvx,tmp_mode
  sll_int32  :: i,j,iv,jv,iflag,ierr,m,n
  sll_int32  :: my_num, num_threads, comm
  
  sll_real64, dimension(:), allocatable :: dlog_density,inv_Te


  
  !MESH
  sll_int32       :: Nr,Ntheta,Nphi,Nvpar,nbiter,fdiag,fthdiag
  sll_real64      :: a,rhomin,rhomax,Lz,aspect_ratio,vmax,dt,R0

  !EQUIL
  sll_int32       :: modethmin,modethmax,modezmin,modezmax
  logical         :: zonal_flow
  sll_real64      :: rpeak,kappan,kappaTi,kappaTe,deltarn,deltarTi,deltarTe,epsilon,tmp_tab(4)
  
  sll_int32  :: grad,carac,bc(2)     
  
  




  
  namelist /MESH/ Nr,Ntheta,Nphi,Nvpar,a,rhomin,rhomax,Lz,aspect_ratio,vmax,nbiter,dt,&
           &fdiag,fthdiag    


!  &MESH
!  Nr           = 127
!  Ntheta       = 128
!  Nphi         = 32
!  Nvpar        = 47
!  a            = 32.
!  rhomin       = 0.2
!  rhomax       = 0.8
!  Ltheta       = 6.283185307179586476925286766559005768394
!  Lz         = 6.283185307179586476925286766559005768394
!  aspect_ratio = 3.
!  vmax      = 6.
!  nbiter      = 30
!  dt      = 16.
!  fdiag      = 30
!  fthdiag      = 30



  namelist /EQUIL/ rpeak,kappan,kappaTi,kappaTe,deltarn,deltarTi,deltarTe,&
           &epsilon,modethmin,modethmax,modezmin,modezmax,zonal_flow



!  &EQUIL
!  rpeak          = 0.5     
!  kappan         = 4.
!  kappaTi        = 27.
!  kappaTe        = 27.
!  deltarn        = 0.08
!  deltarTi       = 0.08
!  deltarTe       = 0.08
!  epsilon        = .001
!  modethmin       = 1
!  modethmax       = 16
!  modezmin        = 1
!  modezmax        = 8
!  zonal_flow     = .true.





  my_num = sll_get_collective_rank(sll_world_collective)
  num_threads = sll_get_collective_size(sll_world_collective)
  comm   = sll_world_collective%comm
  
  if (my_num == MPI_MASTER) then
    call fichinit()
    read(idata,NML=MESH)
  end if

  !MESH  
  call mpi_bcast(Nr,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Ntheta,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Nphi,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(Nvpar,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(a,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(rhomin,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(rhomax,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(Lz,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(aspect_ratio,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(vmax,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(nbiter,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(dt,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(fdiag,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fthdiag,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)

  

  call new(geomx,a*rhomin,0._f64,a*rhomax,2._f64*sll_pi,Nr+1,Ntheta,iflag,"pery")  
  !call new(geomv,0._f64,-vmax,Lz,vmax,Nphi,Nvpar,iflag,"perxy")
  call new(geomv,0._f64,-vmax,Lz,vmax,Nphi,Nvpar,iflag,"perx")
  
  R0 = a*aspect_ratio


  !print *,'eh ben alors ?',my_num,geomx%nx



  !here we choice the boundary conditions, the <ay to derive phi and the way to compute
  !the characteristics
  
  bc = (/SLL_DIRICHLET,SLL_DIRICHLET/)
  grad = 3  !3: splines 1: finite diff order 2 
  carac = 4 !4: Verlet, 1: explicit Euler
   


  !my_num      = sll_get_collective_rank(sll_world_collective)
  !num_threads = sll_get_collective_size(sll_world_collective)
  !comm        = sll_world_collective%comm

  print *,'#xmin, xmax=',geomx%x0,geomx%x1
  
  !print *,'xmax=',geomx%x0+(real(geomx%nx,f64)-1)*geomx%dx
  
  
  !stop

  if (my_num == MPI_MASTER) then
   
    call fichinit()
    read(idata,NML=EQUIL)
    
    rpeak = geomx%x0+rpeak*(geomx%x1-geomx%x0)
    
    !kappan = kappan/(geomx%x1-geomx%x0)   !new before *
    !kappaTi = kappaTi/(geomx%x1-geomx%x0) !new before *
    !kappaTe = kappaTe/(geomx%x1-geomx%x0) !new before *
    !deltarn = deltarn*(geomx%x1-geomx%x0) 
    !deltarTi = deltarTi*(geomx%x1-geomx%x0)
    !deltarTe = deltarTe*(geomx%x1-geomx%x0)

  end if

  !EQUIL
  call mpi_bcast(rpeak,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(kappan,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(kappaTi,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(kappaTe,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(deltarn,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(deltarTi,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(deltarTe,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(epsilon,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
  call mpi_bcast(modethmin,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(modethmax,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(modezmin,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(modezmax,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(zonal_flow,   1,MPI_LOGICAL ,MPI_MASTER,comm,ierr)




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
  !print*,'#init zone ',my_num,jstartx,jendx,ipiece_size_x, &
  !                           jstartv,jendv,ipiece_size_v
  
  
  
  SLL_ALLOCATE(f(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)
  SLL_ALLOCATE(f_old(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)
  SLL_ALLOCATE(f_equil(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)

  SLL_ALLOCATE(rho(geomx%nx,geomx%ny+1,geomv%nx+1),iflag)
  SLL_ALLOCATE(phi(geomx%nx,geomx%ny+1,geomv%nx+1),iflag)
  SLL_ALLOCATE(adv_field(3,geomx%nx,geomx%ny+1,geomv%nx+1),iflag)
  SLL_ALLOCATE(e_x(geomx%nx,geomx%ny),iflag)
  SLL_ALLOCATE(e_y(geomx%nx,geomx%ny),iflag)
  
  SLL_ALLOCATE(profile(3,geomx%nx),iflag)
  
  !call compute_profile(vlas2d,profile(1,1:geomx%nx),geomx,rpeak,deltarn,kappan,R0)
  !call compute_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak,deltarTi,kappaTi,R0)
  !call normalize_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak)

  print *,'#rpeak=',rpeak
  print *,'#deltarn=',deltarn
  print *,'#kappan=',kappan
  print *,'#R0=',R0
  print *,'#deltarTi=',deltarTi
  print *,'#epsilon=',epsilon
  
  !epsilon=._f64
  
  R0=1._f64

  !redefinition of profile
  
  call compute_profile(vlas2d,profile(1,1:geomx%nx),geomx,rpeak,deltarn,kappan,R0)
  call compute_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak,deltarTi,kappaTi,R0)
  call compute_profile(vlas2d,profile(3,1:geomx%nx),geomx,rpeak,deltarTe,kappaTe,R0)
  call normalize_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak)
  call normalize_profile(vlas2d,profile(3,1:geomx%nx),geomx,rpeak)

  !initialization of physical constants in vlas2d
  vlas2d%rpeak=rpeak
  vlas2d%kappan=kappan
  vlas2d%kappaTi=kappaTi
  vlas2d%kappaTe=kappaTe
  vlas2d%deltarn=deltarn
  vlas2d%deltarTi=deltarTi
  vlas2d%deltarTe=deltarTe
  vlas2d%n0rp=n0rp_compute(vlas2d,profile(1,1:geomx%nx),geomx,rpeak)
  
  
  print *,'#n0rp=',vlas2d%n0rp,profile(1,geomx%nx/2),profile(1,geomx%nx/2+1)



  SLL_ALLOCATE(dlog_density(geomx%nx),iflag)
  SLL_ALLOCATE(inv_Te(geomx%nx),iflag)
  
  do i=1,geomx%nx
    inv_Te(i) = 1._f64/profile(3,i)
    x = geomx%x0+(i-1)*geomx%dx
    dlog_density(i) = -kappan/R0*cosh((x-rpeak)/deltarn)**(-2)
  enddo
  
  
  if (my_num == MPI_MASTER) then
    !print *,'#the param are:',kappaTi,deltarTi,rpeak,geomx%x0,geomx%x1,geomx%x0+(geomx%nx*geomx%dx)
    !stop
   ky  = 2._f64*pi/(real(geomx%ny,f64)*geomx%dy)
   
   !kappaTi=4._f64/(geomx%x1-geomx%x0)
   !deltarTi=0.1_f64*(geomx%x1-geomx%x0)
   !kappan=0.8_f64/(geomx%x1-geomx%x0)
   !deltarn=0.2_f64*(geomx%x1-geomx%x0)
   
   tmp=0._f64
   do i=1,geomx%nx
      x = geomx%x0+(i-1)*geomx%dx
      tmp=tmp+profile(1,i)*x*geomx%dx
   enddo
   print *,'#int no=',tmp
   tmp=0._f64
   do i=1,geomx%nx
      x = geomx%x0+(i-1)*geomx%dx
      tmp=tmp+x*geomx%dx
   enddo
   print *,'#int r=',tmp
   
   !profile(1,:)=profile(1,:)/108._f64*15._f64
   !stop
    do i=1,geomx%nx
      x = geomx%x0+(i-1)*geomx%dx
      !print *,'x=',x,i
              tmp_tab(1)=2._f64/(deltarn/deltarTi)
              tmp_tab(2)=tmp_tab(1)-tmp_tab(1)**2*(x-rpeak)**2
              tmp_tab(2)=tmp_tab(2)+(1._f64/x+dlog_density(i))*tmp_tab(1)*(x-rpeak)
              tmp_tab(2)=tmp_tab(2)+(real(modethmin,f64)*ky/x)**2+inv_Te(i)
              tmp_tab(2)=tmp_tab(2)*exp(-(x-rpeak)**2/(deltarn/deltarTi))
      !print *,'tmp_tab(2)=',tmp_tab(2)
      write(20,*),x,profile(1,i),profile(2,i),profile(3,i),&
       exp(-kappaTi*deltarTi*tanh((x-rpeak)/deltarTi)),&
      exp(-kappan*deltarn*tanh((x-rpeak)/deltarn)),&
      exp(-0.4*tanh((x-rpeak)/(deltarTi*(geomx%x1-geomx%x0)))),&
      exp(-kappaTi/R0*deltarTi*tanh((x-rpeak)/deltarTi)),&
      exp(-kappaTe/R0*deltarTe*tanh((x-rpeak)/deltarTe)),&
      exp(-kappan/R0*deltarn*tanh((x-rpeak)/deltarn)),&
      -kappaTi/R0*cosh((x-rpeak)/deltarTe)**(-2)/dlog_density(i),&
       tmp_tab(2)

      !dlog_density(i),-kappaTi/R0*cosh((x-rpeak)/deltarTe)**(-2),&
      !&-kappaTi/R0*cosh((x-rpeak)/deltarTe)**(-2)/dlog_density(i)
    enddo
    print *,'#kappaTi/R0=',kappaTi/R0
    print *,'#kappaTe/R0=',kappaTe/R0
    print *,'#kappan/R0=',kappan/R0
    print *,'#deltarTi=',deltarTi
    print *,'#deltarTe=',deltarTe
    print *,'#deltarn=',deltarn
    
    !stop

  endif



  !warning redefining profile for initialization of distribution function  
  !profile(1,:)=1._f64
  !profile(2,:)=1._f64
  
  profile(1,:)=profile(1,:)/sqrt(2._f64*sll_pi*profile(2,:))
  profile(2,:)=-0.5_f64/profile(2,:)
  
  !print *,'#m=',modethmin,modethmax
  
  !stop
  
  
  
  ky  = 2._f64*pi/(real(geomx%ny,f64)*geomx%dy)
  kvx  = 2._f64*pi/(real(geomv%nx,f64)*geomv%dx)
  
  !ky=1._f64
  
  do jv=jstartv,jendv
     vy = geomv%y0+real(jv-1,f64)*geomv%dy
     v2 = vy*vy
     do iv=1,geomv%nx
        vx = geomv%x0+real(iv-1,f64)*geomv%dx
        do j=1,geomx%ny
           y=geomx%y0+real(j-1,f64)*geomx%dy           
           tmp_mode=0._f64
           do m=modethmin,modethmax
             do n=modezmin,modezmax
               tmp_mode=tmp_mode+cos(real(n,f64)*kvx*vx+real(m,f64)*ky*y)
             enddo
           enddo  
!           do i=1,geomx%nx
!    x=geomx%x0+(i-1)*geomx%dx
!    tmp_mode = exp(-(x-rpeak)**2/(deltarn/deltarTi))
!    profile(1,i)=profile(1,i)*tmp_mode
!  enddo
           !tmp_mode=1._f64+tmp_mode*epsilon
           do i=1,geomx%nx
              x=geomx%x0+(i-1)*geomx%dx              
              !f(i,j,iv,jv)=tmp_mode*profile(1,i)*exp(v2*profile(2,i))
              
              tmp_tab(1)=2._f64/(deltarn/deltarTi)
              tmp_tab(2)=tmp_tab(1)-tmp_tab(1)**2*(x-rpeak)**2
              tmp_tab(2)=tmp_tab(2)+(1._f64/x+dlog_density(i))*tmp_tab(1)*(x-rpeak)
              tmp_tab(2)=tmp_tab(2)+(real(modethmin,f64)*ky/x)**2+inv_Te(i)
              !tmp_tab(2)=1._f64
              tmp_tab(2)=tmp_tab(2)*exp(-(x-rpeak)**2/(deltarn/deltarTi))
              tmp_tab(2)=exp(-(x-rpeak)**2/(4._f64*deltarn/deltarTi))
              !tmp_tab(2)=exp(-(x-rpeak)**2/(2._f64*deltarn/deltarTi))
              
              !if(i==1) print *,'tmp_tab=',tmp_tab(2),4._f64*deltarn/deltarTi
              f_equil(i,j,iv,jv)=compute_equil(vlas2d,x,vy)
              !profile(1,i)*exp(v2*profile(2,i))
              f(i,j,iv,jv)=(1._f64+tmp_mode*epsilon*tmp_tab(2))*&
                &f_equil(i,j,iv,jv)
                
           end do
        end do
     end do
  end do
  
  print *,'#diff=',maxval(abs(f-f_equil)),epsilon
  
  !stop
  
  

  
  
  
  
  
  !,geomx%nx*geomx%ny*geomv%nx*geomv%ny,&
  !geomx%nx*geomx%ny*geomv%nx*geomv%ny-mass*num_threads!*geomx%dx*geomx%dy*geomv%dx*geomv%dy,2*pi*2*pi*12*0.6*32


  !redefinition of profile
  
  call compute_profile_analytic(vlas2d,geomx,profile)
  
  !call compute_profile(vlas2d,profile(1,1:geomx%nx),geomx,rpeak,deltarn,kappan,R0)
  !call compute_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak,deltarTi,kappaTi,R0)
  !call compute_profile(vlas2d,profile(3,1:geomx%nx),geomx,rpeak,deltarTe,kappaTe,R0)
  !call normalize_profile(vlas2d,profile(2,1:geomx%nx),geomx,rpeak)
  !call normalize_profile(vlas2d,profile(3,1:geomx%nx),geomx,rpeak)

   !profile(1,:)=profile(1,:)/108._f64*15._f64
  
  
  
  plan_sl => new_SL(geomx%x0,geomx%x1,geomx%dx,geomx%dy,dt,geomx%nx-1,geomx%ny,grad,carac,bc)
  
  
  
  
  plan_poisson => new_plan_poisson_polar(geomx%dx,geomx%x0,geomx%nx-1,geomx%ny,bc,&
    &dlog_density,inv_Te)
  
  
  
  
  
  SLL_DEALLOCATE_ARRAY(dlog_density,iflag)
  SLL_DEALLOCATE_ARRAY(inv_Te,iflag)
  
  
  
  !call new(poisson, rho, geomx, iflag)
  call new(vlas2d, geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  !call new(splx,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  !call new(sply,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  
  
 end subroutine init_dk







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

  xi  = 0.90_f64
  eps = 0.05_f64
  kx  = 2._f64*pi/((geomx%nx)*geomx%dx)
  ky  = 2._f64*pi/((geomx%ny)*geomx%dy)
  do jv=jstartv,jendv
     vy = geomv%y0+(jv-1)*geomv%dy
     do iv=1,geomv%nx
        vx = geomv%x0+(iv-1)*geomv%dx
        v2 = vx*vx+vy*vy
        do j=1,geomx%ny
           y=geomx%y0+(j-1)*geomx%dy
           do i=1,geomx%nx
              x=geomx%x0+(i-1)*geomx%dx
              f(i,j,iv,jv)=(1+eps*cos(kx*x))*1._f64/(2._f64*pi)*exp(-0.5_f64*v2)
           end do
        end do
     end do
  end do

  call new(poisson, rho, geomx, iflag)
  call new(vlas2d, geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(splx,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)
  call new(sply,   geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)

 end subroutine initlocal





end program vp4d_dk
