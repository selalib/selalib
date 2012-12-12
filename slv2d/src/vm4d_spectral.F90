program vm4d_spectral

#include "selalib.h"

  use used_precision  
  use geometry_module
  use diagnostiques_module
  use sll_vlasov4d_spectral
  use sll_cubic_spline_interpolator_1d
  use sll_cubic_spline_interpolator_2d
  use sll_maxwell
  use sll_maxwell_2d_pstd
  use poisson2d_periodic
  use remapper

  implicit none

  type(geometry)     :: geomx 
  type(geometry)     :: geomv 
  type(vlasov4d)     :: vlas4d 
  type(maxwell_pstd) :: maxwell
  type(poisson2dpp)  :: poisson 

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
  type(cubic_spline_1d_interpolator), target :: spl_x3
  type(cubic_spline_1d_interpolator), target :: spl_x4
  type(cubic_spline_2d_interpolator), target :: spl_x3x4

  sll_int32  :: nbiter, iter , fdiag, fthdiag  
  sll_real64 :: dt, nrj, tcpu1, tcpu2

  sll_int32  :: prank, comm
  sll_int64  :: psize

  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
  sll_int32  :: jstartx, jendx, jstartv, jendv   

  call sll_boot_collective()
  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  tcpu1 = MPI_WTIME()
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

  call transposexv(vlas4d)
  call densite_charge(vlas4d)
  call solve(poisson,vlas4d%ex,vlas4d%ey,vlas4d%rho,nrj)

  call transposevx(vlas4d)
  call advection_x1(vlas4d,0.5*dt)
  call advection_x2(vlas4d,0.5*dt)

  do iter=1,nbiter

     if (iter ==1 .or. mod(iter,fdiag) == 0) then 
        call write_xmf_file(vlas4d,iter/fdiag)
     end if

     !call densite_charge(vlas4d,rho)
     call transposexv(vlas4d)
     !call densite_courant(vlas4d)
     !call ampere_te(maxwell,vlas4d%ex,vlas4d%ey,vlas4d%bz,dt,vlas4d%jx,vlas4d%jy) 
     call densite_charge(vlas4d)
     call solve(poisson,vlas4d%ex,vlas4d%ey,vlas4d%rho,nrj)

     call advection_x3(vlas4d,dt)
     call advection_x4(vlas4d,dt)
     !call advection_x3x4(vlas4d,dt)

     call transposevx(vlas4d)

     call advection_x1(vlas4d,dt)
     call advection_x2(vlas4d,dt)

     if (mod(iter,fthdiag).eq.0) then 
        nrj=sum(vlas4d%ex*vlas4d%ex+vlas4d%ey*vlas4d%ey)*(vlas4d%geomx%dx)*(vlas4d%geomx%dy)
        nrj=0.5_wp*log(nrj)
        call thdiag(vlas4d,nrj,iter*dt)
     endif

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize


  call free(poisson)
  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains

!####################################################################################

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
    sll_int32       :: error,ierr  ! indicateur d'erreur

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

    call new(geomx,x0,y0,x1,y1,nx,ny,error,"perxy")
    call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,error,"perxy")

  end subroutine initglobal

  subroutine initlocal(jstartx,jendx,jstartv,jendv)

    sll_int32  :: jstartx 
    sll_int32  :: jendx 
    sll_int32  :: jstartv
    sll_int32  :: jendv   
    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,k,l,error
    sll_real64 :: xi, eps, kx, ky
    sll_int32  :: gi, gj, gk, gl
    sll_int32, dimension(4) :: global_indices

    sll_int32 :: psize

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    call spl_x1%initialize(geomx%nx, geomx%x0, geomx%x1, PERIODIC_SPLINE)
    call spl_x2%initialize(geomx%nx, geomx%y0, geomx%y1, PERIODIC_SPLINE)
    call spl_x3%initialize(geomv%nx, geomv%x0, geomv%x1, PERIODIC_SPLINE)
    call spl_x4%initialize(geomv%ny, geomv%y0, geomv%y1, PERIODIC_SPLINE)

    call spl_x3x4%initialize(geomv%nx, geomv%ny,                        &
    &                        geomv%x0, geomv%x1, geomv%y0, geomv%y1,    &
    &                        PERIODIC_SPLINE, PERIODIC_SPLINE)

    call new(vlas4d,geomx,geomv,spl_x1,spl_x2,spl_x3,spl_x4,spl_x3x4,error)

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

    call initialize(maxwell, geomx%x0, geomx%x1, geomx%nx, &
                             geomx%y0, geomx%y1, geomx%ny, TE_POLARIZATION)

    call new(poisson, vlas4d%ex, vlas4d%ey, geomx, error)

    jstartx = get_layout_4D_j_min( vlas4d%layout_v, prank )
    jendx   = get_layout_4D_j_max( vlas4d%layout_v, prank )
    jstartv = get_layout_4D_l_min( vlas4d%layout_x, prank )
    jendv   = get_layout_4D_l_max( vlas4d%layout_x, prank )

  end subroutine initlocal

end program vm4d_spectral
