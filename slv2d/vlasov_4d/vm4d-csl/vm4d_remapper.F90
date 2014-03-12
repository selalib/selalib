program vm4d

#include "selalib-mpi.h"

  use geometry_module
  use diagnostiques_module
  use sll_vlasov4d_base
  use sll_vlasov4d_maxwell

  implicit none

  type(geometry)            :: geomx 
  type(geometry)            :: geomv 
  type(vlasov4d_maxwell)    :: vlasov4d 
  type(maxwell_2d_pstd)     :: maxwell
  type(poisson_2d_periodic) :: poisson 

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
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
     write(*,"(g13.3,1x,3i5)") dt,nbiter,fdiag,fthdiag
  endif

  call initlocal(jstartx,jendx,jstartv,jendv)

  call transposexv(vlasov4d)
  call compute_charge(vlasov4d)
  call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
  !call faraday(maxwell, vlasov4d%ex, vlasov4d%ey, vlasov4d%bz, 0.5*dt)   
  call transposevx(vlasov4d)
  call advection_x1(vlasov4d,0.5*dt)
  call advection_x2(vlasov4d,0.5*dt)

  do iter=1,nbiter

     if (iter ==1 .or. mod(iter,fdiag) == 0) then 
        call write_xmf_file(vlasov4d,iter/fdiag)
     end if

     call transposexv(vlasov4d)
     call compute_current(vlasov4d)
     call ampere(maxwell,vlasov4d%ex,vlasov4d%ey,vlasov4d%bz,dt,vlasov4d%jx,vlasov4d%jy) 
     !call faraday(maxwell, vlasov4d%ex, vlasov4d%ey, vlasov4d%bz, 0.5*dt)   
     call advection_x3x4(vlasov4d,dt)
     !call faraday(maxwell, vlasov4d%ex, vlasov4d%ey, vlasov4d%bz, 0.5*dt)   
     call transposevx(vlasov4d)
     call advection_x1(vlasov4d,dt)
     call advection_x2(vlasov4d,dt)

     if (mod(iter,fthdiag).eq.0) then 
        nrj=sum(vlasov4d%ex*vlasov4d%ex+vlasov4d%ey*vlasov4d%ey)*(vlasov4d%geomx%dx)*(vlasov4d%geomx%dy)
        nrj=0.5_wp*log(nrj)
        call thdiag(vlasov4d,nrj,iter*dt)
     endif

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call delete(poisson)
  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains

!####################################################################################

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

    call spl_x1%initialize(geomx%nx, geomx%x0, geomx%x1, SLL_PERIODIC)
    call spl_x2%initialize(geomx%ny, geomx%y0, geomx%y1, SLL_PERIODIC)

    call spl_x3x4%initialize(geomv%nx, geomv%ny,                        &
    &                        geomv%x0, geomv%x1, geomv%y0, geomv%y1,    &
    &                        SLL_PERIODIC, SLL_PERIODIC)

    call new(vlasov4d,geomx,geomv,spl_x1,spl_x2,spl_x3x4,error)

    call compute_local_sizes_4d(vlasov4d%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

    xi  = 0.90_f64
    eps = 0.05_f64
    kx  = 2_f64*sll_pi/((geomx%nx)*geomx%dx)
    ky  = 2_f64*sll_pi/((geomx%ny)*geomx%dy)

    do l=1,loc_sz_l 
    do k=1,loc_sz_k
    do j=1,loc_sz_j
    do i=1,loc_sz_i

       global_indices = local_to_global_4D(vlasov4d%layout_x,(/i,j,k,l/)) 
       gi = global_indices(1)
       gj = global_indices(2)
       gk = global_indices(3)
       gl = global_indices(4)

       x  = geomx%x0+(gi-1)*geomx%dx
       y  = geomx%y0+(gj-1)*geomx%dy
       vx = geomv%x0+(gk-1)*geomv%dx
       vy = geomv%y0+(gl-1)*geomv%dy

       v2 = vx*vx+vy*vy
       vlasov4d%f(i,j,k,l)=(1._f64+eps*cos(kx*x))*1._f64/(2._f64*sll_pi)*exp(-0.5_f64*v2)

    end do
    end do
    end do
    end do

    call initialize(maxwell, geomx%x0, geomx%x1, geomx%nx, &
                             geomx%y0, geomx%y1, geomx%ny, TE_POLARIZATION)

    call initialize(poisson, geomx%x0, geomx%x1, geomx%nx, &
                             geomx%y0, geomx%y1, geomx%ny, error)

    jstartx = get_layout_4D_j_min( vlasov4d%layout_v, prank )
    jendx   = get_layout_4D_j_max( vlasov4d%layout_v, prank )
    jstartv = get_layout_4D_l_min( vlasov4d%layout_x, prank )
    jendv   = get_layout_4D_l_max( vlasov4d%layout_x, prank )

  end subroutine initlocal

end program vm4d
