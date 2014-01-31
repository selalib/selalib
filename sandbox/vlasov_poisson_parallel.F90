program vlasov_poisson_parallel

#include "selalib-mpi.h"

  implicit none

  sll_int32,  parameter :: nstep    = 200
  sll_real64, parameter :: delta_t  = 0.05
  sll_real64, parameter :: eta1_min = 0.0_f64
  sll_real64, parameter :: eta1_max = 4*sll_pi
  sll_real64, parameter :: eta2_min = 0.0_f64
  sll_real64, parameter :: eta2_max = 2*sll_pi
  sll_real64, parameter :: eta3_min = -6.0_f64
  sll_real64, parameter :: eta3_max = +6.0_f64
  sll_real64, parameter :: eta4_min = -6.0_f64
  sll_real64, parameter :: eta4_max = +6.0_f64

  sll_int32, parameter :: nc_eta1 = 32
  sll_int32, parameter :: nc_eta2 = 32
  sll_int32, parameter :: nc_eta3 = 32
  sll_int32, parameter :: nc_eta4 = 32

  sll_real64 :: delta_eta1 = (eta1_max-eta1_min)/nc_eta1
  sll_real64 :: delta_eta2 = (eta2_max-eta2_min)/nc_eta2
  sll_real64 :: delta_eta3 = (eta3_max-eta3_min)/nc_eta3
  sll_real64 :: delta_eta4 = (eta4_max-eta4_min)/nc_eta4

  sll_real64 :: eps
  

  type(poisson_2d_periodic) :: poisson 

  sll_real64, dimension(:,:,:,:),  pointer   :: f
  sll_real64, dimension(:,:,:,:),  pointer   :: ft
  type(layout_4D), pointer                   :: layout_x
  type(layout_4D), pointer                   :: layout_v
  type(remap_plan_4D_real64), pointer        :: x_to_v 
  type(remap_plan_4D_real64), pointer        :: v_to_x
  sll_real64, dimension(nc_eta1+1,nc_eta2+1) :: ex
  sll_real64, dimension(nc_eta1+1,nc_eta2+1) :: ey
  sll_real64, dimension(nc_eta1+1,nc_eta2+1) :: rho

  class(sll_interpolator_1d_base), pointer   :: interp_eta1
  class(sll_interpolator_1d_base), pointer   :: interp_eta2
  class(sll_interpolator_1d_base), pointer   :: interp_eta3
  class(sll_interpolator_1d_base), pointer   :: interp_eta4

  type(cubic_spline_1d_interpolator), target :: spl_eta1
  type(cubic_spline_1d_interpolator), target :: spl_eta2
  type(cubic_spline_1d_interpolator), target :: spl_eta3
  type(cubic_spline_1d_interpolator), target :: spl_eta4

  sll_int32  :: istep
  sll_int32  :: prank, comm
  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
  sll_int64  :: psize
  sll_real64 :: tcpu1, tcpu2
  sll_int32  :: i, j, k, l, gi, gj, gk, gl
  sll_int32  :: global_indices(4)
  sll_real64 :: x, y, vx, vy, v2, kx
  sll_int32  :: error

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

  call spl_eta1%initialize(nc_eta1, eta1_min, eta1_max, SLL_PERIODIC)
  call spl_eta2%initialize(nc_eta2, eta2_min, eta2_max, SLL_PERIODIC)
  call spl_eta3%initialize(nc_eta3, eta3_min, eta3_max, SLL_HERMITE)
  call spl_eta4%initialize(nc_eta4, eta4_min, eta4_max, SLL_HERMITE)

  interp_eta1 => spl_eta1
  interp_eta2 => spl_eta2
  interp_eta3 => spl_eta3
  interp_eta4 => spl_eta4

  layout_x => new_layout_4D( sll_world_collective )        
  call initialize_layout_with_distributed_4D_array( &
             nc_eta1+1, nc_eta2+1, nc_eta3+1, nc_eta4+1,    &
             1,1,1,int(psize,4),layout_x)

  if ( prank == MPI_MASTER ) call sll_view_lims_4D( layout_x )
  call flush(6)

  call compute_local_sizes_4d(layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(f(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),error)

  layout_v => new_layout_4D( sll_world_collective )
  call initialize_layout_with_distributed_4D_array( &
              nc_eta1+1, nc_eta2+1, nc_eta3+1, nc_eta4+1,    &
              1,int(psize,4),1,1,layout_v)

  if ( prank == MPI_MASTER ) call sll_view_lims_4D( layout_v )
  call flush(6)

  call compute_local_sizes_4d(layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(ft(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),error)

  x_to_v => new_remap_plan( layout_x, layout_v, f)     
  v_to_x => new_remap_plan( layout_v, layout_x, ft)     
  
  eps  = 0.05_f64
  kx   = 0.5_f64

  do l=1,loc_sz_l 
  do k=1,loc_sz_k
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(layout_v,(/i,j,k,l/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     gk = global_indices(3)
     gl = global_indices(4)

     x  = eta1_min+(gi-1)*delta_eta1
     y  = eta2_min+(gj-1)*delta_eta2
     vx = eta3_min+(gk-1)*delta_eta3
     vy = eta4_min+(gl-1)*delta_eta4

     v2 = vx*vx+vy*vy

     ft(i,j,k,l)=(1+eps*cos(kx*x))*exp(-.5*v2)/(2.*sll_pi)

  end do
  end do
  end do
  end do

  call initialize(poisson, eta1_min, eta1_max, nc_eta1, &
                           eta2_min, eta2_max, nc_eta2, error)

  call compute_charge()
  call solve(poisson,ex,ey,rho)

  call advection_eta1(0.5*delta_t)
  call advection_eta2(0.5*delta_t)

  do istep=1, nstep

     call apply_remap_4D( x_to_v, f, ft )

     call compute_charge()
     call solve(poisson,ex,ey,rho)

     call advection_eta3(delta_t)
     call advection_eta4(delta_t)

     call apply_remap_4D( v_to_x, ft, f )

     call advection_eta1(delta_t)
     call advection_eta2(delta_t)

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains

 subroutine compute_charge()

   sll_int32                          :: error
   sll_int32                          :: comm
   sll_real64                         :: dvxvy

   sll_real64, dimension(nc_eta1+1,nc_eta2+1) :: locrho

   dvxvy = delta_eta3*delta_eta4

   call compute_local_sizes_4d(layout_v, &
        loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
   
   locrho = 0.
   do j=1,loc_sz_j
   do i=1,loc_sz_i
      global_indices = local_to_global_4D(layout_v,(/i,j,1,1/)) 
      gi = global_indices(1)
      gj = global_indices(2)
      locrho(gi,gj) = sum(ft(i,j,:,:))*dvxvy 
   end do
   end do
   rho = 0.
   comm  = sll_world_collective%comm
   call mpi_barrier(comm,error)
   call mpi_allreduce(locrho,rho,(nc_eta1+1)*(nc_eta2+1),MPI_REAL8,MPI_SUM,comm,error)

 end subroutine compute_charge


 subroutine advection_eta1(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_4d(layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  do l=1,loc_sz_l
  do k=1,loc_sz_k
     global_indices = local_to_global_4D(layout_x,(/1,1,k,l/)) 
     gk = global_indices(3)
     alpha = (eta3_min +(gk-1)*delta_eta3)*dt
     do j=1,loc_sz_j
        f(:,j,k,l) = &
           interp_eta1%interpolate_array_disp(loc_sz_i,f(:,j,k,l),alpha)
     end do
  end do
  end do

 end subroutine advection_eta1

 subroutine advection_eta2(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_4d(layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l

    global_indices = local_to_global_4D(layout_x,(/1,1,1,l/)) 
    gl = global_indices(4)
    alpha = (eta4_min +(gl-1)*delta_eta4)*dt

    do k=1,loc_sz_k
    do i=1,loc_sz_i

       f(i,:,k,l) = interp_eta2%interpolate_array_disp(loc_sz_j,f(i,:,k,l),alpha)

    end do
    end do

  end do

 end subroutine advection_eta2

 subroutine advection_eta3(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_4d(layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  do l=1,loc_sz_l
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(layout_v,(/i,j,1,l/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = ex(gi,gj)*dt

     ft(i,j,:,l) = interp_eta3%interpolate_array_disp(loc_sz_k,ft(i,j,:,l),alpha)

  end do
  end do
  end do

 end subroutine advection_eta3

 subroutine advection_eta4(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_4d(layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  do k=1,loc_sz_k
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(layout_v,(/i,j,1,1/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     alpha = ey(gi,gj)*dt
     ft(i,j,k,:) = interp_eta4%interpolate_array_disp(loc_sz_l,ft(i,j,k,:),alpha)

  end do
  end do
  end do

 end subroutine advection_eta4

end program vlasov_poisson_parallel
