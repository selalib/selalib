program vm4d

#include "selalib.h"
  use used_precision  
  use geometry_module
  use diagnostiques_module
  use sll_vlasov4d
  use sll_cubic_spline_interpolator_1d
  use maxwell2dfdtd_module
  use poisson2d_periodic
  use remapper

  implicit none

  type(geometry)       :: geomx 
  type(geometry)       :: geomv 
  type(vlasov4d)       :: vlas4d 
  type (maxwell2dfdtd) :: maxw2dfdtd 

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
  type(cubic_spline_1d_interpolator), target :: spl_x3
  type(cubic_spline_1d_interpolator), target :: spl_x4

  sll_real64, dimension(:,:),     pointer :: bz
  sll_real64, dimension(:,:),     pointer :: ex
  sll_real64, dimension(:,:),     pointer :: ey 
  sll_real64, dimension(:,:),     pointer :: jx
  sll_real64, dimension(:,:),     pointer :: jy 

  sll_int32  :: nbiter, iter 
  sll_int32  :: fdiag, fthdiag  

  sll_real64 :: dt     
  sll_real64 :: nrj, tcpu1, tcpu2

  sll_int32  :: prank, comm
  sll_int64  :: psize

  sll_int32                                   :: loc_sz_i
  sll_int32                                   :: loc_sz_j
  sll_int32                                   :: loc_sz_k
  sll_int32                                   :: loc_sz_l

  sll_int32                                   :: jstartx 
  sll_int32                                   :: jendx 
  sll_int32                                   :: jstartv
  sll_int32                                   :: jendv   

  call sll_boot_collective()
  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm   = sll_world_collective%comm

  tcpu1 = MPI_WTIME()
  if (.not. is_power_of_two(psize)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if
  if (prank == MPI_MASTER) then
     print*,'MPI Version of slv2d running on ',psize, ' processors'
  end if

  call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)

  if (prank == MPI_MASTER) then
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

  call initlocal(jstartx,jendx,jstartv,jendv)

  call plot_mesh4d(geomx,geomv,jstartx,jendx,jstartv,jendv)

  call transposevx(vlas4d)
  call advection_x1(vlas4d,0.5*dt)
  call advection_x2(vlas4d,0.5*dt)

  do iter=1,nbiter

     if (mod(iter,fdiag) == 0) then 
        call plot_df(vlas4d%f,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, XY_VIEW)
        call plot_df(vlas4d%f,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, XVX_VIEW)
        call plot_df(vlas4d%f,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, YVY_VIEW)
        call plot_df(vlas4d%f,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, VXVY_VIEW)
     end if

     call transposexv(vlas4d)

     call densite_courant(vlas4d,jx,jy)
     
     if (iter == 1) then
     !   call solve_faraday(maxw2dfdtd,ex,ey,bz,0.5_f64*dt)
         call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,0.5_f64*dt)
     else
     !   call solve_faraday(maxw2dfdtd,ex,ey,bz,dt)
         call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)
     end if
     
     call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)

     call advection_x3(vlas4d,dt,ex,bz)
     call advection_x4(vlas4d,dt,ey,bz)

     call transposevx(vlas4d)

     call advection_x1(vlas4d,dt)
     call advection_x2(vlas4d,dt)

     if (mod(iter,fthdiag).eq.0) then 
      call thdiag(vlas4d,nrj,iter*dt)
     endif

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize


  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

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

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm   = sll_world_collective%comm

    if (prank == MPI_MASTER) then

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
    call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"perxy")

  end subroutine initglobal

  subroutine initlocal(jstartx,jendx,jstartv,jendv)

    sll_int32  :: jstartx 
    sll_int32  :: jendx 
    sll_int32  :: jstartv
    sll_int32  :: jendv   
    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,k,l,iflag
    sll_real64 :: xi, eps, kx, ky
    sll_int32  :: gi, gj, gk, gl
    sll_int32, dimension(4) :: global_indices

    type(poisson2dpp) :: poisson 
    sll_real64, dimension(:,:), pointer :: rho 
    sll_int32 :: psize

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    call spl_x1%initialize(geomx%nx, geomx%x0, geomx%x1, PERIODIC_SPLINE)
    call spl_x2%initialize(geomx%ny, geomx%y0, geomx%y1, PERIODIC_SPLINE)
    call spl_x3%initialize(geomv%nx, geomv%x0, geomv%x1, PERIODIC_SPLINE)
    call spl_x4%initialize(geomv%ny, geomv%y0, geomv%y1, PERIODIC_SPLINE)

    call new(vlas4d,geomx,geomv,spl_x1,spl_x2,spl_x3,spl_x4,iflag)

    call compute_local_sizes_4d(vlas4d%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

    xi  = 0.90_f64
    eps = 0.05_f64
    kx  = 2_f64*sll_pi/((geomx%nx)*geomx%dx)
    ky  = 2_f64*sll_pi/((geomx%ny)*geomx%dy)

    do l=1,loc_sz_l 
    do k=1,loc_sz_k
    do j=1,loc_sz_j
    do i=1,loc_sz_i

       global_indices = local_to_global_4D(vlas4d%layout_x,(/i,j,k,l/)) 
       gi = global_indices(1)
       gj = global_indices(2)
       gk = global_indices(3)
       gl = global_indices(4)

       x  = geomx%x0+(gi-1)*geomx%dx
       y  = geomx%y0+(gj-1)*geomx%dy
       vx = geomv%x0+(gk-1)*geomv%dx
       vy = geomv%y0+(gl-1)*geomv%dy

       v2 = vx*vx+vy*vy

       vlas4d%f(i,j,k,l)=(1+eps*cos(kx*x))*1/(2*sll_pi)*exp(-.5*v2)

    end do
    end do
    end do
    end do

    SLL_ALLOCATE(ex(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(ey(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(jx(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(jy(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(bz(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(rho(geomx%nx,geomx%ny),iflag)

    call new(poisson, ex, ey, geomx, iflag)
    call transposexv(vlas4d)
    call densite_charge(vlas4d,rho)
    call solve(poisson,ex,ey,rho,nrj)
    write(*,*) "NRJ =",nrj
    call free(poisson)

    jstartx = get_layout_4D_j_min( vlas4d%layout_v, prank )
    jendx   = get_layout_4D_j_max( vlas4d%layout_v, prank )
    jstartv = get_layout_4D_l_min( vlas4d%layout_x, prank )
    jendv   = get_layout_4D_l_max( vlas4d%layout_x, prank )

    call new(maxw2dfdtd, geomx, iflag, jstartx, jendx)

  end subroutine initlocal

end program vm4d
