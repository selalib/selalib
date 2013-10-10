program vp4d

#include "selalib-mpi.h"

  use sll_vlasov4d_base
  use sll_vlasov4d_poisson

  implicit none

  type(sll_logical_mesh_2d), pointer :: geomx 
  type(sll_logical_mesh_2d), pointer :: geomv 
  type(vlasov4d_poisson)    :: vlasov4d 
  type(poisson_2d_periodic) :: poisson 

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
  type(cubic_spline_1d_interpolator), target :: spl_x3
  type(cubic_spline_1d_interpolator), target :: spl_x4

  sll_int32  :: nbiter, iter, fdiag, fthdiag  
  sll_int32  :: prank, comm
  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
  sll_int32  :: jstartx, jendx, jstartv, jendv   
  sll_int64  :: psize
  sll_real64 :: dt, nrj, tcpu1, tcpu2

  call sll_boot_collective()

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

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

  call initlocal(jstartx,jendx,jstartv,jendv)

  call advection_x1(vlasov4d,0.5*dt)
  call advection_x2(vlasov4d,0.5*dt)

  do iter=1,nbiter

     if (iter == 1 .or. mod(iter,fdiag) == 0) then 
        call write_xmf_file(vlasov4d,iter/fdiag)
     end if

     call transposexv(vlasov4d)

     call compute_charge(vlasov4d)
     call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)

     call advection_x3(vlasov4d,dt)
     call advection_x4(vlasov4d,dt)

     call transposevx(vlasov4d)

     call advection_x1(vlasov4d,dt)
     call advection_x2(vlasov4d,dt)

     if (mod(iter,fthdiag).eq.0) then 
        call write_energy(vlasov4d,iter*dt)
     end if

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call delete(poisson)

  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains


  subroutine initlocal(jstartx,jendx,jstartv,jendv)

    use sll_xdmf_parallel

    sll_int32               :: jstartx, jendx, jstartv, jendv   
    sll_int32               :: i,j,k,l,error
    sll_real64              :: vx,vy,v2,x,y
    sll_real64              :: xi, eps, kx, ky
    sll_int32               :: gi, gj, gk, gl
    sll_int32, dimension(4) :: global_indices
    sll_real64              :: dimx, dimy, factor
    character(len=4)        :: prefix = "vp4d"
    integer(HID_T)          :: file_id
    integer(HSSIZE_T)       :: offset(2)
    integer(HSIZE_T)        :: global_dims(2)
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2

    prank = sll_get_collective_rank(sll_world_collective)
    comm  = sll_world_collective%comm


    call spl_x1%initialize(nc_eta1, eta1_min, eta1_max, SLL_PERIODIC)
    call spl_x2%initialize(nc_eta2, eta2_min, eta2_max, SLL_PERIODIC)
    call spl_x3%initialize(nc_eta3, eta3_min, eta3_max, SLL_PERIODIC)
    call spl_x4%initialize(nc_eta4, eta4_min, eta4_max, SLL_PERIODIC)

    call initialize(vlasov4d,geomx,geomv,spl_x1,spl_x2,spl_x3,spl_x4,error)

    call compute_local_sizes_4d(vlasov4d%layout_x, &
         loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

    xi  = 0.90_f64
    eps = 0.05_f64
    dimx  = 4*sll_pi !geomx%nx*geomx%dx
    dimy  = 4*sll_pi !geomx%ny*geomx%dy
    kx  = 2_f64*sll_pi/dimx
    ky  = 2_f64*sll_pi/dimy

    factor = dimx*dimy*(2*sll_pi+eps*(sin(kx*eta1_max)-sin(kx*eta1_min)))
    factor = 2.*sll_pi
    
    do l=1,loc_sz_l 
    do k=1,loc_sz_k
    do j=1,loc_sz_j
    do i=1,loc_sz_i

       global_indices = local_to_global_4D(vlasov4d%layout_x,(/i,j,k,l/)) 
       gi = global_indices(1)
       gj = global_indices(2)
       gk = global_indices(3)
       gl = global_indices(4)

       x  = eta1_min+(gi-1)*delta_eta1
       y  = eta2_min+(gj-1)*delta_eta2
       vx = eta3_min+(gk-1)*delta_eta3
       vy = eta4_min+(gl-1)*delta_eta4

       v2 = vx*vx+vy*vy

       vlasov4d%f(i,j,k,l)=(1+eps*cos(kx*x))*exp(-.5*v2)/factor

    end do
    end do
    end do
    end do

    call initialize(poisson, eta1_min, eta1_max, nc_eta1, &
                             eta2_min, eta2_max, nc_eta2, error)

    jstartx = get_layout_4D_j_min( vlasov4d%layout_v, prank )
    jendx   = get_layout_4D_j_max( vlasov4d%layout_v, prank )
    jstartv = get_layout_4D_l_min( vlasov4d%layout_x, prank )
    jendv   = get_layout_4D_l_max( vlasov4d%layout_x, prank )

    call compute_charge(vlasov4d)
    call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)


    SLL_ALLOCATE(x1(nc_eta1,nc_eta2),error)
    SLL_ALLOCATE(x2(nc_eta1,nc_eta2),error)
   
    do j=1,nc_eta2
    do i=1,nc_eta1
      x1(i,j) = eta1_min + (i-1) * delta_eta1
      x2(i,j) = eta2_min + (j-1) * delta_eta2
    end do
    end do

    global_dims(1) = nc_eta1
    global_dims(2) = nc_eta2
    offset = 0

    call sll_xdmf_open(prank,"fields.xmf",prefix,nc_eta1,nc_eta2,file_id,error)
    call sll_xdmf_write_array(prefix,global_dims,offset,x1,'x1',error)
    call sll_xdmf_write_array(prefix,global_dims,offset,x2,'x2',error)
    call sll_xdmf_write_array(prefix,global_dims,offset, &
                              vlasov4d%rho,"rho",error,file_id,"Node")
    call sll_xdmf_write_array(prefix,global_dims,offset, &
                              vlasov4d%ex ,"ex" ,error,file_id,"Node")
    call sll_xdmf_close(file_id,error)

end subroutine initlocal



end program vp4d
