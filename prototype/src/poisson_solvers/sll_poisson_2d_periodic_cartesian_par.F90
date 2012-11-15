!> @brief 
!> Selalib periodic 2D poisson solver for cartesian coordinates.
!> Start date: Sep. 4, 2012
!> Last modification: Sep. 4, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!**************************************************************************

module sll_poisson_2d_periodic_cartesian_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
  use remapper
  use sll_fft
  use numeric_constants
  use sll_collective

  implicit none

  ! While the solution region is 2D, the following solver is meant to be used 
  ! with 4D arrays with size 1 in the last two dimensions. The reason for this
  ! is that the input data, the rho array, is the result of a reduction
  ! process in the last 2 dimensions. The computation of the potential
  ! which may be an in-place operation, needs to be treated as a 4D array
  ! in order to remap it appropriately with the 4D distribution function.
  ! There might be ways around this, like using 'reshape'...

  type poisson_2d_periodic_plan_cartesian_par
     sll_int32                                 :: ncx   ! number of cells  
     sll_int32                                 :: ncy   ! number of cells
     sll_real64                                :: Lx    ! domain length 
     sll_real64                                :: Ly    ! domain length
     type(sll_fft_plan), pointer               :: px
     type(sll_fft_plan), pointer               :: py
     type(sll_fft_plan), pointer               :: px_inv
     type(sll_fft_plan), pointer               :: py_inv
     type(layout_2D),  pointer                 :: layout_seq_x1
     type(layout_2D),  pointer                 :: layout_seq_x2
     sll_int32                                 :: seq_x1_local_sz_x1
     sll_int32                                 :: seq_x1_local_sz_x2
     sll_int32                                 :: seq_x2_local_sz_x1
     sll_int32                                 :: seq_x2_local_sz_x2
     sll_comp64, dimension(:,:), pointer       :: fft_x_array
     sll_comp64, dimension(:,:), pointer       :: fft_y_array
     type(remap_plan_2D_comp64), pointer       :: rmp_xy
     type(remap_plan_2D_comp64), pointer       :: rmp_yx
  end type poisson_2d_periodic_plan_cartesian_par

contains

  ! Presently, this function receives the geometric information as 
  ! individual arguments. We should consider passing the 'simple geometry'
  ! object that we have for the cartesian cases.
  function new_poisson_2d_periodic_plan_cartesian_par( &
    start_layout, &
    ncx, &
    ncy, &
    Lx, &
    Ly ) result(plan)

    type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan
    type(layout_2D), pointer                     :: start_layout
    sll_int32                                    :: ncx
    sll_int32                                    :: ncy
    sll_real64                                   :: Lx
    sll_real64                                   :: Ly
    sll_int64                                    :: colsz ! collective size
    type(sll_collective_t), pointer              :: collective
    ! number of processors
    sll_int64                                    :: nprocx1
    sll_int64                                    :: nprocx2
    sll_int32                                    :: ierr 
    sll_int32                                    :: loc_sz_x1
    sll_int32                                    :: loc_sz_x2
    sll_int32                                    :: seq_x1_local_sz_x1
    sll_int32                                    :: seq_x1_local_sz_x2
    sll_int32                                    :: seq_x2_local_sz_x1
    sll_int32                                    :: seq_x2_local_sz_x2

    ! The collective to be used is the one that comes with the given layout.
    collective => get_layout_collective( start_layout )
    colsz      = sll_get_collective_size( collective )

    if ( (.not.is_power_of_two(int(ncx,i64))) .and. &
         (.not.is_power_of_two(int(ncy,i64))) ) then     
       print *, 'This test needs to run with numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan, ierr)

    ! Geometry
    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_seq_x1 => start_layout
    call compute_local_sizes_2d( &
         plan%layout_seq_x1, &
         loc_sz_x1, &
         loc_sz_x2 )

    plan%seq_x1_local_sz_x1 = loc_sz_x1
    plan%seq_x1_local_sz_x2 = loc_sz_x2

    SLL_ALLOCATE( plan%fft_x_array(loc_sz_x1,loc_sz_x2),ierr)

    ! For FFTs (in x-direction)
    plan%px => fft_new_plan( &
         loc_sz_x1, &
         loc_sz_x2, &
         plan%fft_x_array, &
         plan%fft_x_array, &
         FFT_FORWARD, &
         FFT_ONLY_FIRST_DIRECTION)!+FFT_NORMALIZE )

    plan%px_inv => fft_new_plan( &
         loc_sz_x1, &
         loc_sz_x2, &
         plan%fft_x_array, &
         plan%fft_x_array, &
         FFT_INVERSE, &
         FFT_ONLY_FIRST_DIRECTION) !+FFT_NORMALIZE )

    ! Layout and local sizes for FFTs in y-direction (x2)
    plan%layout_seq_x2 => new_layout_2D( collective )
    nprocx1 = colsz
    nprocx2 = 1

    call initialize_layout_with_distributed_2D_array( &
         ncx, &
         ncy, &
         int(nprocx1,i32), &
         int(nprocx2,i32), &
         plan%layout_seq_x2 )

    call compute_local_sizes_2d( &
         plan%layout_seq_x2, &
         loc_sz_x1, &
         loc_sz_x2 )
print *, 'inside new plan: ', loc_sz_x1, loc_sz_x2
    plan%seq_x2_local_sz_x1 = loc_sz_x1
    plan%seq_x2_local_sz_x2 = loc_sz_x2
    SLL_ALLOCATE( plan%fft_y_array(loc_sz_x1,loc_sz_x2), ierr )

    ! For FFTs (in y-direction)

    plan%py => fft_new_plan( &
         loc_sz_x1, &
         loc_sz_x2, &
         plan%fft_y_array, &
         plan%fft_y_array, &
         FFT_FORWARD, &
         FFT_ONLY_SECOND_DIRECTION)! + FFT_NORMALIZE )

    plan%py_inv => fft_new_plan( &
         loc_sz_x1, &
         loc_sz_x2, &
         plan%fft_y_array, &
         plan%fft_y_array, &
         FFT_INVERSE, &
         FFT_ONLY_SECOND_DIRECTION)! + FFT_NORMALIZE )

    plan%rmp_xy => &
     new_remap_plan(plan%layout_seq_x1, plan%layout_seq_x2, plan%fft_x_array)
    plan%rmp_yx => &
     new_remap_plan(plan%layout_seq_x2, plan%layout_seq_x1, plan%fft_y_array)
  end function new_poisson_2d_periodic_plan_cartesian_par

  subroutine solve_poisson_2d_periodic_cartesian_par(plan, rho, phi)
    type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan
    sll_real64, dimension(:,:)                     :: rho
    sll_real64, dimension(:,:)                     :: phi
    ! global sizes
    sll_int32                                      :: ncx, ncy
    sll_int32                                      :: npx_loc, npy_loc
    sll_int32                                      :: i, j
    sll_int32                                      :: ierr
    ! Reciprocals of domain lengths.
    sll_real64                                     :: r_Lx, r_Ly
    sll_real64                                     :: kx, ky
    sll_comp64                                     :: val
    sll_real64                                     :: normalization
    sll_int32                                      :: myrank
    sll_int64                                      :: colsz ! collective size
    type(layout_2D), pointer                       :: layout_x
    type(layout_2D), pointer                       :: layout_y
    sll_int32, dimension(1:2)                      :: global
    sll_int32                                      :: gi, gj

    ncx  = plan%ncx
    ncy  = plan%ncy
    r_Lx = 1.0_f64/plan%Lx
    r_Ly = 1.0_f64/plan%Ly
    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_seq_x1
    layout_y => plan%layout_seq_x2
    call verify_argument_sizes_par(layout_x, rho, phi)

    ! FFTs in x-direction
    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 

    ! The input is handled internally as complex arrays
    plan%fft_x_array = -cmplx(rho, 0_f64, kind=f64)

    call fft_apply_plan(plan%px, plan%fft_x_array, plan%fft_x_array)

    ! FFTs in y-direction
    npx_loc = plan%seq_x2_local_sz_x1
    npy_loc = plan%seq_x2_local_sz_x2
 
    call apply_remap_2D( plan%rmp_xy, plan%fft_x_array, plan%fft_y_array )
 
    call fft_apply_plan(plan%py, plan%fft_y_array, plan%fft_y_array) 
print *, 'inside solve, ffty array = ', plan%fft_y_array
    ! This should be inside the FFT plan...
    normalization = 1.0_f64/(ncx*ncy)

    ! Apply the kernel 
    do j=1, npy_loc
       do i=1, npx_loc
          ! Make sure that the first mode is set to zero so that we get an
          ! answer with zero mean value. This step assumes that the (1,1) point
          ! will always be at the border of any splitting of the domains. This 
          ! seems safe in this case.
          
          global = local_to_global_2D( layout_y, (/i, j/))
          gi = global(1)
          gj = global(2)
          
          if( (gi == 1) .and. (gj == 1) ) then
             ! NOTE: WE HAVE A POTENTIALLY DANGEROUS SITUATION HERE. THE 2D
             ! FFT WE HAVE CARRIED OUT WAS THE RESULT OF A SEQUENCE OF TWO
             ! INDEPENDENT FFTS. HENCE WHEN WE CALL fft_set_mode_complx_2d()
             ! BASED ON ONE OF THE FFT PLANS, IT WILL BELIEVE THAT ONLY ONE
             ! DIRECTION WAS TRANSFORMED. IN CASE THAT THE UNDERLYING ARRAYS
             ! HAVE BEEN SCRAMBLED IN ANY WAY, THIS WOULD YIELD THE WRONG
             ! RESULT. WE MAY BE GETTING LUCKY HERE IF EVERYTHING IS IN 
             ! NATURAL ORDER, BUT THIS NEEDS TO BE FIXED. CAN THIS BE FIXED BY
             ! USING ONLY 1D FFTs? HARDLY, BECAUSE OF THE SAME SITUATION. IT 
             ! SEEMS THAT THE 2d FFT DONE IN A PARALLEL WAY, NEEDS TO BE 
             ! IMPLEMENTED INDEPENDENTLY AND WITH KNOWLEDGE OF LAYOUTS, REMAP, 
             ! ETC. ...
             call fft_set_mode_complx_2d( &
                  plan%py,&
                  plan%fft_y_array,&
                  (0.0_f64,0.0_f64),&
                  1,&
                  1)
          else
             kx  = real(gi-1,f64)
             ky  = real(gj-1,f64)
             ! Crucial step: While the results of the FFT are basically a
             ! list of the modes in the 0..n-1 sense, the algorithm that we
             ! are using uses a decomposition with negative frequencies,
             ! i.e.: -n/2..n/2-1. Therefore a step is needed in which we are
             ! effectively reinterpreting the FFT results in this light.
             ! More specifically, we take the second half of the transformed
             ! array and apply a modified factor, proper of a negative 
             ! frequency.
             if( kx .ge. ncx/2 ) then
                kx = kx - ncx
             end if

             if( ky .ge. ncy/2 ) then
                ky = ky - ncy
             end if

              val = -plan%fft_y_array(i,j)*normalization / &
                  ( ( (kx*r_Lx)**2 + (ky*r_Ly)**2)*4.0_f64*sll_pi**2)
              call fft_set_mode_complx_2d(plan%py,plan%fft_y_array,val,i,j)
          end if
       enddo
    enddo

    ! Inverse FFTs in y-direction
    call fft_apply_plan(plan%py_inv, plan%fft_y_array, plan%fft_y_array) 

    ! Prepare to take inverse FFTs in x-direction
    call apply_remap_2D( plan%rmp_yx, plan%fft_y_array, plan%fft_x_array )

    npx_loc = plan%seq_x1_local_sz_x1 
    npy_loc = plan%seq_x1_local_sz_x2 
    call fft_apply_plan(plan%px_inv, plan%fft_x_array, plan%fft_x_array)
    phi = real(plan%fft_x_array, f64)
  end subroutine solve_poisson_2d_periodic_cartesian_par


  subroutine delete_poisson_2d_periodic_plan_cartesian_par(plan)
    type (poisson_2d_periodic_plan_cartesian_par), pointer :: plan
    sll_int32                                              :: ierr

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_poisson_3d_periodic_plan_par(): ', &
            'passed plan is not associated.'
       STOP
    end if

    call fft_delete_plan(plan%px)
    call fft_delete_plan(plan%py)
    call fft_delete_plan(plan%px_inv)
    call fft_delete_plan(plan%py_inv)

!    call delete( plan%layout_x ) ! can't delete this, the plan does not own it
    call delete( plan%layout_seq_x1 )
    call delete( plan%layout_seq_x2 )
    SLL_DEALLOCATE_ARRAY(plan%fft_x_array, ierr)
    SLL_DEALLOCATE_ARRAY(plan%fft_y_array, ierr)
    call delete( plan%rmp_xy )
    call delete( plan%rmp_yx )
    SLL_DEALLOCATE(plan, ierr)
  end subroutine delete_poisson_2d_periodic_plan_cartesian_par

  subroutine verify_argument_sizes_par(layout, rho, phi)
    type(layout_2D), pointer       :: layout
    sll_real64, dimension(:,:)     :: rho
    sll_real64, dimension(:,:)     :: phi
    sll_int32,  dimension(2)       :: n ! nx_loc, ny_loc
    sll_int32                      :: i

    ! Note that this checks for strict sizes, not an array being bigger
    ! than a certain size, but exactly a desired size... This may be a bit
    ! too stringent.
    call compute_local_sizes_2d( layout, n(1), n(2) )

    do i=1,2
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
          print*, 'ERROR: solve_poisson_2d_periodic_cartesian_par()', &
               'size of either rho or phi does not match expected size. '
          if (i==1) then
             print*, 'solve_poisson_3d_periodic_cartesian_par(): ', &
                  'mismatch in direction x'
          elseif (i==2) then
             print*, 'solve_poisson_3d_periodic_cartesian_par(): ', &
                  'mismatch in direction y'
          endif
          print *, 'Exiting...'
          stop
       endif
    enddo
  end subroutine verify_argument_sizes_par


end module sll_poisson_2d_periodic_cartesian_par
