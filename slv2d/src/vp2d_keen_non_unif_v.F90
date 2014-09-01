program vp2d_keen_non_unif_v

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use used_precision  
  use geometry_module
  use diagnostiques_module
  use vlasov4d_plot
#ifdef _FFTW
  use poisson2d_periodic
#else
  use poisson2dpp_seq
#endif
!  use sll_vlasov2d
  use vlasov2d_non_unif_v_module
  use splinepp_class
  use splinenn_class
  use interp_non_unif_pp_class

  implicit none

  type :: drive_param
     sll_real64 :: t0, twL, twR, tstart, tflat, tL, tR
     sll_real64 :: adr
     logical    :: turn_drive_off
     sll_real64 :: kdr
     sll_real64 :: omegadr
     sll_real64 :: Edrmax
  end type drive_param

  type(geometry)    :: geomx 
  type(geometry)    :: geomv 
  type(vlasov2d)    :: vlas2d 
  type(poisson2dpp) :: poisson 
  type(drive_param) :: dr_param
  ! old version of splines
  !type (splinepp)    :: splx 
  !type (splinenn)    :: sply

  sll_real64, dimension(:,:,:,:), pointer :: f4d
  sll_real64, dimension(:,:),     pointer :: rho
  sll_real64, dimension(:,:),     pointer :: e_x
  sll_real64, dimension(:,:),     pointer :: e_y 

  sll_int32  :: nbiter, iter 
  sll_int32  :: fdiag, fthdiag  
  sll_int32  :: jstartx, jendx, jstartv, jendv

  sll_real64 :: dt     
  sll_real64 :: nrj, tcpu1, tcpu2
  sll_real64 :: xx, hypergaussian

  sll_int32  :: i,j 

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

  call initglobal(geomx,geomv,dr_param,dt,nbiter,fdiag,fthdiag)

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
     write(*,*) 'dt, nbiter, fdiag, fthdiag:'
     write(*,"(g13.3,1x,3i6)") dt,nbiter,fdiag,fthdiag
  endif

  call initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
       f4d,rho,e_x,e_y,vlas2d,poisson)

  call plot_mesh4d(geomx,geomv,jstartx,jendx,jstartv,jendv)

  !call advection_x1(vlas2d,f4d,0.5*dt)
  !call advection_x2(vlas2d,f4d,0.5*dt)
  call advection_x(vlas2d,f4d,.5*dt)

  do iter=1,nbiter

     call transposexv(vlas2d,f4d)

     ! compute density (q is not applied)
     call densite_charge(vlas2d,rho)

     ! here rho is the density and e_x, e_y is the force (q=-1 compensate)
     call solve(poisson,e_x,e_y,rho,nrj)

     ! add applied field (potential is of form a(t) * cos(k*x-w*t)
     ! because of q=-1 this yields f_x_app = a*k*sin(k*x-w*t), f_y_app = 0
     call PFenvelope(dr_param%adr, iter*dt, dr_param%tflat, dr_param%tL, &
          dr_param%tR, dr_param%twL, dr_param%twR, &
          dr_param%t0, dr_param%turn_drive_off)
     do j = 1, geomx%ny
        do i = 1, geomx%nx
           ! use normalized hypergaussian envelope in x
           xx = geomx%x0 + (i-1)*geomx%dx
           hypergaussian = 0._f64*exp(-0.5_8 * xx * xx) ** 3 / 1.4472025  
           e_y(i,j) = e_y(i,j) +  hypergaussian * dr_param%Edrmax * &
                dr_param%adr * dr_param%kdr * &
                sin(dr_param%kdr * (j-1) * geomx%dy - dr_param%omegadr*iter*dt)
        enddo
     end do

     !call advection_x3(vlas2d,e_x,dt)
     !call advection_x4(vlas2d,e_y,dt)
     call advection_v(vlas2d,e_x,e_y,dt)
     
     call transposevx(vlas2d,f4d)

     if (mod(iter,fthdiag).eq.0) then
        !call advection_x1(vlas2d,f4d,0.5*dt)
        !call advection_x2(vlas2d,f4d,0.5*dt)
        call advection_x(vlas2d,f4d,.5*dt)
        call thdiag(vlas2d,f4d,nrj,iter*dt,jstartv)  
        if (mod(iter,fdiag) == 0) then 
           call plot_df(f4d,iter/fdiag,geomx,geomv,jstartx,jendx, &
                jstartv,jendv, XY_VIEW)
           call plot_df(f4d,iter/fdiag,geomx,geomv,jstartx,jendx, &
                jstartv,jendv, XVX_VIEW)
           call plot_df(f4d,iter/fdiag,geomx,geomv,jstartx,jendx, &
                jstartv,jendv, YVY_VIEW)
           call plot_df(f4d,iter/fdiag,geomx,geomv,jstartx,jendx, &
                jstartv,jendv, VXVY_VIEW)
        end if
        !call advection_x2(vlas2d,f4d,0.5*dt)
        !call advection_x1(vlas2d,f4d,0.5*dt)
        call advection_x(vlas2d,f4d,.5*dt)
     else
        !call advection_x1(vlas2d,f4d,1.0*dt)
        !call advection_x2(vlas2d,f4d,1.0*dt)
        call advection_x(vlas2d,f4d, dt)
     end if

  end do

  tcpu2 = MPI_WTIME()
  if (my_num == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*num_threads

  call sll_halt_collective()

  print*,'PASSED'


contains


  subroutine initglobal(geomx,geomv,dr_param,dt,nbiter,fdiag,fthdiag)

    type(geometry)  :: geomx       ! geometrie globale du probleme
    type(geometry)  :: geomv       ! geometrie globale du probleme
    type(drive_param) :: dr_param  ! drive parameters
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
    sll_real64 :: t0, twL, twR, tstart, tflat, tL, tR
    sll_real64 :: adr
    logical    :: turn_drive_off
    sll_real64 :: kdr
    sll_real64 :: omegadr
    sll_real64 :: Edrmax
    sll_int32       :: iflag,ierr  ! indicateur d'erreur

    namelist /time/ dt, nbiter
    namelist /diag/ fdiag, fthdiag
    namelist /phys_space/ x0,x1,y0,y1,nx,ny
    namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
    namelist / drive / t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, &
         Edrmax, omegadr, kdr
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
       read(idata,NML=drive)
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

    call mpi_bcast(t0,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(twL,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(twR,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(tstart,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(tflat,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(tL,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(tR,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(turn_drive_off, 1,MPI_LOGICAL,MPI_MASTER,comm,ierr)
    call mpi_bcast(Edrmax,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(omegadr,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(kdr,     1,MPI_REAL8,MPI_MASTER,comm,ierr)

    dr_param%t0 = t0
    dr_param%twL = twL 
    dr_param%twR = twR
    dr_param%tstart = tstart 
    dr_param%tflat = tflat
    dr_param%tL = tL
    dr_param%tR = tR
    dr_param%turn_drive_off = turn_drive_off
    dr_param%Edrmax = Edrmax
    dr_param%omegadr = omegadr
    dr_param%kdr = kdr

    call initialize(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")
    call initialize(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"perxy")

!print*, 'bcast ',my_num, t0, twL, twR, tstart, tflat, tL, tR, turn_drive_off, Edrmax, omegadr, kdr
  end subroutine initglobal

  subroutine initlocal(geomx,geomv,jstartv,jendv,jstartx,jendx, &
       f,rho,e_x,e_y,vlas2d,poisson)

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

    type(poisson2dpp) :: poisson

    sll_int32  :: ipiece_size_v
    sll_int32  :: ipiece_size_x

    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,iv,jv,iflag
    sll_int32  :: my_num, num_threads, comm
    sll_real64 :: xi, eps, kx, ky

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

    kx  = 2_f64*sll_pi/((geomx%nx)*geomx%dx)
    ky  = 2_f64*sll_pi/((geomx%ny)*geomx%dy)
    eps=0.05_f64
    do jv=jstartv,jendv
       !vy = geomv%y0+(jv-1)*geomv%dy
       vy = geomv%y0+vlas2d%interpv%node_pos_y(jv)*(geomv%y1-geomv%y0)
       do iv=1,geomv%nx
          !vx = geomv%x0+(iv-1)*geomv%dx
          vx = geomv%x0+vlas2d%interpv%node_pos_x(iv)*(geomv%x1-geomv%x0)
          v2 = vx*vx+vy*vy
          do j=1,geomx%ny
             y=geomx%y0+(j-1)*geomx%dy
             do i=1,geomx%nx
                x=geomx%x0+(i-1)*geomx%dx
                !f(i,j,iv,jv)= 1/(2*sll_pi)*exp(-.5*v2)
                !f(i,j,iv,jv)= 1._f64/(2._f64*sll_pi)*exp(-.5_f64*v2)*(1._f64+0.01_f64*cos(kx*x)*cos(ky*y))
                f(i,j,iv,jv)=(1._f64+eps*cos(kx*x))*1._f64/(2._f64*sll_pi)*exp(-.5_f64*v2)
             end do
          end do
       end do
    end do

#ifdef _FFTW
    call new(poisson, e_x, e_y,   geomx, iflag)
#else
    call new(poisson, rho,   geomx, iflag)
#endif

    call new(vlas2d, geomx, geomv, iflag, jstartx, jendx, jstartv, jendv)

  end subroutine initlocal

  subroutine PFenvelope(S, t, tflat, tL, tR, twL, twR, t0, &
       turn_drive_off)

    ! DESCRIPTION
    ! -----------
    ! S: the wave form at a given point in time. This wave form is 
    !    not scaled (its maximum value is 1).
    ! t: the time at which the envelope is being evaluated
    ! tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
    !    envelope, defined in the main portion of this program.
    ! turn_drive_off: 1 if the drive should be turned off after a time
    !    tflat, and 0 otherwise

    sll_real64, intent(in) :: t, tflat, tL, tR, twL, twR, t0
    sll_real64, intent(out) :: S
    logical, intent(in) :: turn_drive_off
    ! local variables
    sll_int32 :: i 
    sll_real64 :: epsilon

    ! The envelope function is defined such that it is zero at t0,
    ! rises to 1 smoothly, stay constant for tflat, and returns
    ! smoothly to zero.
    if(turn_drive_off) then
       epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
       S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
       S = S / (1-epsilon)
    else
       epsilon = 0.5*(tanh((t0-tL)/twL) + 1)
       S = 0.5*(tanh((t-tL)/twL) + 1) - epsilon
       S = S / (1-epsilon)
    endif
    if(S<0) then
       S = 0.
    endif
    return
  end subroutine PFenvelope

end program vp2d_keen_non_unif_v
