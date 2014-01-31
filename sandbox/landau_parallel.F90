program landau_parallel

#include "selalib-mpi.h"

  implicit none

  class(sll_interpolator_1d_base), pointer   :: interp_x
  class(sll_interpolator_1d_base), pointer   :: interp_v
  type(cubic_spline_1d_interpolator), target :: spl_x
  type(cubic_spline_1d_interpolator), target :: spl_v
  sll_real64, dimension(:,:),  pointer       :: f
  sll_real64, dimension(:,:),  pointer       :: ft
  type(layout_2D), pointer                   :: layout_x
  type(layout_2D), pointer                   :: layout_v
  type(remap_plan_2D_real64), pointer        :: x_to_v 
  type(remap_plan_2D_real64), pointer        :: v_to_x
  sll_real64, dimension(:),   pointer        :: efield
  sll_real64, dimension(:),   pointer        :: rho
  type(poisson_1d_periodic)                  :: poisson

  sll_int32  :: i_step, j_step
  sll_int32  :: prank, comm
  sll_int32  :: loc_sz_i, loc_sz_j
  sll_int64  :: psize
  sll_real64 :: tcpu1
  sll_real64 :: tcpu2
  sll_int32  :: gi
  sll_int32  :: gj
  sll_int32  :: global_indices(2)
  sll_int32  :: error
  sll_real64 :: x
  sll_real64 :: v
  sll_real64 :: eps
  sll_real64 :: kx
  sll_real64 :: vsq

!###########################################################
  !Simulation parameters and geometry sizes                !
                                                           !
  sll_int32,  parameter :: nc_x = 256                      !
  sll_int32,  parameter :: nc_v = 256                      !
  sll_real64, parameter :: x_min = 0.0_f64                 !
  sll_real64, parameter :: x_max = 4.*sll_pi               !
  sll_real64, parameter :: v_min = - 6.0_f64               !
  sll_real64, parameter :: v_max = + 6.0_f64               !
  sll_real64, parameter :: delta_x = (x_max-x_min)/nc_x    !
  sll_real64, parameter :: delta_v = (v_max-v_min)/nc_v    !
  sll_real64, parameter :: delta_t = 0.05_f64              !
  sll_int32,  parameter :: n_step  = 1000                  !
  sll_real64, dimension(n_step) :: nrj                     !
                                                           !
!###########################################################

  sll_int32 :: i, j

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

  call spl_x%initialize(nc_x+1,x_min,x_max,SLL_PERIODIC)

  call spl_v%initialize(nc_v+1,v_min,v_max,SLL_HERMITE)

  interp_x => spl_x
  interp_v => spl_v

  layout_x => new_layout_2D( sll_world_collective )        

  call initialize_layout_with_distributed_2D_array( &
             nc_x+1, nc_v+1, 1,int(psize,4),layout_x)

  if ( prank == MPI_MASTER ) call sll_view_lims_2D( layout_x )
  call flush(6)

  call compute_local_sizes_2d(layout_x,loc_sz_i,loc_sz_j)        
  SLL_CLEAR_ALLOCATE(f(1:loc_sz_i,1:loc_sz_j),error)

  layout_v => new_layout_2D( sll_world_collective )

  call initialize_layout_with_distributed_2D_array( &
              nc_x+1, nc_v+1, int(psize,4),1,layout_v)

  if ( prank == MPI_MASTER ) call sll_view_lims_2D( layout_v )
  call flush(6)

  call compute_local_sizes_2d(layout_v,loc_sz_i,loc_sz_j)        

  SLL_CLEAR_ALLOCATE(ft(1:loc_sz_i,1:loc_sz_j),error)

  SLL_CLEAR_ALLOCATE(rho(1:nc_x+1),error)
  SLL_CLEAR_ALLOCATE(efield(1:nc_x+1),error)

  x_to_v => new_remap_plan( layout_x, layout_v, f)     
  v_to_x => new_remap_plan( layout_v, layout_x, ft)     
  
  call compute_local_sizes_2d(layout_v,loc_sz_i,loc_sz_j)        

  eps  = 0.05_f64
  kx   = 0.50_f64

  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_2D(layout_v,(/i,j/)) 
     gi = global_indices(1)
     gj = global_indices(2)

     x  = x_min+(gi-1)*delta_x
     v  = v_min+(gj-1)*delta_v

     vsq = v*v

     ft(i,j)=(1.0_f64+eps*cos(kx*x))*exp(-.5*vsq)/(2.*sll_pi)

  end do
  end do

  call initialize(poisson, x_min, x_max, nc_x, error)
  call compute_charge()
  call solve(poisson,efield,rho)

  call apply_remap_2D( v_to_x, ft, f )
  call advection_x(0.5*delta_t)

  do i_step=1, n_step

      call apply_remap_2D( x_to_v, f, ft )

      call compute_charge()
      call solve(poisson,efield,rho)

      call advection_v(delta_t)

      call apply_remap_2D( v_to_x, ft, f )

      call advection_x(delta_t)

      if (prank == MPI_MASTER) then
         nrj(i_step) = 0.5_f64*log(sum(efield*efield)*delta_x)
         write(*,100) .0,n_step*delta_t,-29.5,0.5
         do j_step = 1, i_step
            print'(2e15.3)', (j_step-1)*delta_t, nrj(j_step)
         end do
         print'(a)','e'
      end if

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call delete_layout_2D(layout_x)
  call delete_layout_2D(layout_v)
  SLL_DEALLOCATE_ARRAY(f, error)
  SLL_DEALLOCATE_ARRAY(ft, error)

  call sll_halt_collective()

100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

contains

 subroutine compute_charge()

   sll_int32                     :: comm
   sll_real64, dimension(nc_x+1) :: locrho

   call compute_local_sizes_2d(layout_v,loc_sz_i,loc_sz_j)        
   
   locrho = 0.0_f64
   do i=1,loc_sz_i
      global_indices = local_to_global_2D(layout_v,(/i,1/)) 
      gi = global_indices(1)
      gj = global_indices(2)
      locrho(gi) = sum(ft(i,:))*delta_v 
   end do
   rho = 0.0_f64
   comm  = sll_world_collective%comm
   call mpi_barrier(comm,error)
   call mpi_allreduce(locrho,rho,nc_x+1,MPI_REAL8,MPI_SUM,comm,error)

 end subroutine compute_charge

 subroutine advection_x(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_2d(layout_x,loc_sz_i,loc_sz_j)

  do j=1,loc_sz_j

     global_indices = local_to_global_2D(layout_x,(/1,j/)) 
     gj = global_indices(2)
     alpha = (v_min +(gj-1)*delta_v)*dt
     f(:,j) = interp_x%interpolate_array_disp(loc_sz_i,f(:,j),alpha)

  end do

 end subroutine advection_x

 subroutine advection_v(dt)

  sll_real64, intent(in) :: dt
  sll_real64 :: alpha

  call compute_local_sizes_2d(layout_v,loc_sz_i,loc_sz_j)        

  do i=1,loc_sz_i

     global_indices = local_to_global_2D(layout_v,(/i,1/)) 
     gi = global_indices(1)
     alpha = efield(gi)*dt
     ft(i,:) = interp_v%interpolate_array_disp(loc_sz_j,ft(i,:),alpha)

  end do

 end subroutine advection_v

end program landau_parallel

