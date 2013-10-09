!> @brief 
!> Selalib periodic 3D poisson solver
!> Start date: Feb. 08, 2012
!> Last modification: Sep 11, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!**************************************************************************

module sll_poisson_3d_periodic_par

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
  use sll_fft
  use sll_constants
  use sll_collective
  use sll_remapper
  implicit none

  !> Structure to solve Poisson equation on 3d mesh with periodic boundary
  !> conditions. Solver is parallel and numerical method is based on fft 
  !> transform.  Number of cells, which in this periodic case is equal to 
  !> the number of points.
  type poisson_3d_periodic_plan_par
     sll_int32                             :: ncx       !< number of cells in x
     sll_int32                             :: ncy       !< number of cells in y
     sll_int32                             :: ncz       !< number of cells in z
     sll_real64                            :: Lx        !< x domain length
     sll_real64                            :: Ly        !< y domain length
     sll_real64                            :: Lz        !< z domain length
     type(sll_fft_plan), pointer           :: px        !< fft plan in x
     type(sll_fft_plan), pointer           :: py        !< fft plan in y
     type(sll_fft_plan), pointer           :: pz        !< fft plan in z
     type(sll_fft_plan), pointer           :: px_inv    !< inverse fft in x
     type(sll_fft_plan), pointer           :: py_inv    !< inverse fft in y
     type(sll_fft_plan), pointer           :: pz_inv    !< inverse fft in z
     type(layout_3D),  pointer             :: layout_x  !< x layout for remap
     type(layout_3D),  pointer             :: layout_y  !< y layout for remap
     type(layout_3D),  pointer             :: layout_z  !< z layout for remap
     sll_int32, dimension(3,3)             :: loc_sizes !< local sizes
     sll_comp64, dimension(:,:,:), pointer :: array_x   !< x array component
     sll_comp64, dimension(:,:,:), pointer :: array_y   !< y array component
     sll_comp64, dimension(:,:,:), pointer :: array_z   !< z array component
     type(remap_plan_3D_comp64), pointer   :: rmp3_xy   !< transpose from x to y
     type(remap_plan_3D_comp64), pointer   :: rmp3_yz   !< transpose from y to z
     type(remap_plan_3D_comp64), pointer   :: rmp3_zy   !< transpose from z to y
     type(remap_plan_3D_comp64), pointer   :: rmp3_yx   !< transpose from y to x
  end type poisson_3d_periodic_plan_par

contains


  !> Allocate the structure for the 3d parallel Poisson solver
  function new_poisson_3d_periodic_plan_par( &
    start_layout, &
    ncx, &
    ncy, &
    ncz, &
    Lx, &
    Ly, &
    Lz) result(plan)

    type(layout_3D),  pointer                    :: start_layout !< intiial layout
    sll_int32                                    :: ncx !< number of cells in x
    sll_int32                                    :: ncy !< number of cells in y
    sll_int32                                    :: ncz !< number of cells in z
    sll_real64                                   :: Lx  !< x length
    sll_real64                                   :: Ly  !< y length
    sll_real64                                   :: Lz  !< z length
    sll_comp64,                   dimension(ncx) :: x   !< 1d array in x
    sll_comp64,                   dimension(ncy) :: y   !< 1d array in y
    sll_comp64,                   dimension(ncz) :: z   !< 1d array in z
    sll_int64                                    :: colsz ! collective size
    type(sll_collective_t), pointer              :: collective
    ! npx, npy, npz are the numbers of processors in directions x, y, z
    sll_int32                                    :: npx, npy, npz
    sll_int32                                    :: e
    sll_int32                                    :: ierr
    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_int32, dimension(3,3)                    :: loc_sizes

    collective => get_layout_collective(start_layout)
    colsz      = sll_get_collective_size(collective)

    if ( colsz > min(ncx,ncy,ncz) ) then     
       print *, 'This test needs to run in a number of processes which',  &
                ' is less than or equal', min(ncx,ncy,ncz), ' in order to ', &
                'be able properly split the arrays.'
       print *, 'Exiting...'
       stop
    end if
    if ( (.not.is_power_of_two(int(ncx,i64))) .and. &
         (.not.is_power_of_two(int(ncy,i64))) .and. &
         (.not.is_power_of_two(int(ncz,i64))) ) then     
       print *, 'This test needs to run on numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan, ierr)
    plan%ncx = ncx
    plan%ncy = ncy
    plan%ncz = ncz
    plan%Lx  = Lx
    plan%Ly  = Ly
    plan%Lz  = Lz

    ! For FFTs (in each direction)
    plan%px => fft_new_plan( ncx, x, x, FFT_FORWARD )
    plan%py => fft_new_plan( ncy, y, y, FFT_FORWARD )
    plan%pz => fft_new_plan( ncz, z, z, FFT_FORWARD )

    ! For inverse FFTs (in each direction)
    plan%px_inv => fft_new_plan( ncx, x, x, FFT_INVERSE )
    plan%py_inv => fft_new_plan( ncy, y, y, FFT_INVERSE )
    plan%pz_inv => fft_new_plan( ncz, z, z, FFT_INVERSE )

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_x => start_layout
    call compute_local_sizes_3d( &
         plan%layout_x, &
         loc_sizes(1,1), &
         loc_sizes(1,2), &
         loc_sizes(1,3)  )

    ! Layout and local sizes for FFTs in y-direction
    plan%layout_y => new_layout_3D( collective )
    e = int(log(real(colsz))/log(2.))
    npx = 2**(e/2)
    npy = 1
    npz = 2**(e-e/2)  ! int(colsz)/npx
    call initialize_layout_with_distributed_3D_array( &
         ncx, &
         ncy, &
         ncz, &
         npx, &
         npy, &
         npz, &
         plan%layout_y )

    call compute_local_sizes( &
         plan%layout_y,  &
         loc_sizes(2,1), &
         loc_sizes(2,2), &
         loc_sizes(2,3) )

    ! Layout and local sizes for FFTs in z-direction
    plan%layout_z => new_layout_3D( collective )
    ! npx remains the same. Exchange npy and npz.
    npy = npz
    npz = 1
    call initialize_layout_with_distributed_3D_array( &
         ncx, &
         ncy, &
         ncz, &
         npx, &
         npy, &
         npz, &
         plan%layout_z )

    call compute_local_sizes( &
         plan%layout_z, &
         loc_sizes(3,1), &
         loc_sizes(3,2), &
         loc_sizes(3,3) )

    plan%loc_sizes = loc_sizes

    SLL_ALLOCATE(plan%array_x(loc_sizes(1,1),loc_sizes(1,2),loc_sizes(1,3)),ierr)
    SLL_ALLOCATE(plan%array_y(loc_sizes(2,1),loc_sizes(2,2),loc_sizes(2,3)),ierr)
    SLL_ALLOCATE(plan%array_z(loc_sizes(3,1),loc_sizes(3,2),loc_sizes(3,3)),ierr)

    plan%rmp3_xy => NEW_REMAP_PLAN(plan%layout_x, plan%layout_y, plan%array_x)
    plan%rmp3_yz => NEW_REMAP_PLAN(plan%layout_y, plan%layout_z, plan%array_y)
    plan%rmp3_zy => NEW_REMAP_PLAN(plan%layout_z, plan%layout_y, plan%array_z)
    plan%rmp3_yx => NEW_REMAP_PLAN(plan%layout_y, plan%layout_x, plan%array_y)
  end function new_poisson_3d_periodic_plan_par

  !> Compute the 3d potential from the Poisson equation with periodic
  !> boundary conditions.
  subroutine solve_poisson_3d_periodic_par(plan, rho, phi)
    type (poisson_3d_periodic_plan_par), pointer :: plan !< Solver structure
    sll_real64, dimension(:,:,:)                 :: rho  !< Charge density
    sll_real64, dimension(:,:,:)                 :: phi  !< Electric potential
    sll_int32                                    :: nx, ny, nz
    ! nx, ny, nz are the numbers of points - 1 in directions x, y ,z
    sll_int32                                    :: nx_loc, ny_loc, nz_loc
    sll_int32                                    :: i, j, k
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z
    type(layout_3D), pointer                     :: layout_x
    type(layout_3D), pointer                     :: layout_y
    type(layout_3D), pointer                     :: layout_z
    sll_int32, dimension(1:3)                    :: global
    sll_int32                                    :: gi, gj, gk

    ! Get geometry information
    nx = plan%ncx
    ny = plan%ncy
    nz = plan%ncz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_x
    layout_y => plan%layout_y
    layout_z => plan%layout_z

    call verify_argument_sizes_par(layout_x, rho, phi)

    ! FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    plan%array_x = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz_loc
       do j=1,ny_loc
          call fft_apply_plan(plan%px, plan%array_x(:,j,k), plan%array_x(:,j,k))
       enddo
    enddo

    ! FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    call apply_remap_3D( plan%rmp3_xy, plan%array_x, plan%array_y ) 
    do k=1,nz_loc
       do i=1,nx_loc
          call fft_apply_plan(plan%py, plan%array_y(i,:,k), plan%array_y(i,:,k))
       enddo
    enddo

    ! FFTs in z-direction
    nx_loc = plan%loc_sizes(3,1) 
    ny_loc = plan%loc_sizes(3,2) 
    nz_loc = plan%loc_sizes(3,3)
    call apply_remap_3D( plan%rmp3_yz, plan%array_y, plan%array_z ) 
    do j=1,ny_loc
       do i=1,nx_loc
          call fft_apply_plan(plan%pz, plan%array_z(i,j,:), plan%array_z(i,j,:))
       enddo
    enddo

    ! move this normalization elsewhere where the modes are modified
    plan%array_z = plan%array_z/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    do k=1,nz_loc
       do j=1,ny_loc
          do i=1,nx_loc
             global = local_to_global_3D( layout_z, (/i, j, k/))
             gi = global(1)
             gj = global(2)
             gk = global(3)
             if( (gi==1) .and. (gj==1) .and. (gk==1) ) then
                call fft_set_mode_complx_1d( &
                     plan%pz, &
                     plan%array_z(1,1,:), &
                     (0.0_f64,0.0_f64), &
                     1)
             else
                if (gi<=nx/2) then
                   ind_x = real(gi-1,f64)
                else
                   ind_x = real(nx-(gi-1),f64)
                endif
                if (gj<=ny/2) then
                   ind_y = real(gj-1,f64)
                else
                   ind_y = real(ny-(gj-1),f64)
                endif
                if (gk<=nz/2) then
                   ind_z = real(gk-1,f64)
                else
                   ind_z = real(nz-(gk-1),f64)
                endif
!!$           if ( (ind_x==0) .and. (ind_y==0) .and. (ind_z==0) ) then
!!$               if ( rho(i,j,k) /= 0.d0 ) then     
!!$                  print *,'3D: periodic poisson cannot be solved without', &
!!$                                                      ' global_rho(1,1,1)=0'
!!$                  print *, 'Exiting...'
!!$                  stop
!!$              endif
!!$              plan%array_z(i,j,k) = 0.d0
!!$           else
                plan%array_z(i,j,k) = plan%array_z(i,j,k)/(4*sll_pi**2 * &
                            ((ind_x/Lx)**2 + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny_loc
       do i=1,nx_loc
          call fft_apply_plan( &
               plan%pz_inv, &
               plan%array_z(i,j,:), &
               plan%array_z(i,j,:))
       enddo
    enddo

    ! Inverse FFTs in y-direction
    nx_loc = plan%loc_sizes(2,1) 
    ny_loc = plan%loc_sizes(2,2) 
    nz_loc = plan%loc_sizes(2,3)
    call apply_remap_3D( plan%rmp3_zy, plan%array_z, plan%array_y )
    do k=1,nz_loc
       do i=1,nx_loc
          call fft_apply_plan( &
               plan%py_inv, &
               plan%array_y(i,:,k), &
               plan%array_y(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in x-direction
    nx_loc = plan%loc_sizes(1,1) 
    ny_loc = plan%loc_sizes(1,2) 
    nz_loc = plan%loc_sizes(1,3)
    call apply_remap_3D( plan%rmp3_yx, plan%array_y, plan%array_x ) 
    do k=1,nz_loc
       do j=1,ny_loc
          call fft_apply_plan( &
               plan%px_inv, &
               plan%array_x(:,j,k), &
               plan%array_x(:,j,k) )
       enddo
    enddo

    phi = real(plan%array_x, f64)

  end subroutine solve_poisson_3d_periodic_par


  !> Delete the solver structure
  subroutine delete_poisson_3d_periodic_plan_par(plan)
    type (poisson_3d_periodic_plan_par), pointer :: plan
    sll_int32                                    :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call fft_delete_plan(plan%px)
    call fft_delete_plan(plan%py)
    call fft_delete_plan(plan%pz)

    call fft_delete_plan(plan%px_inv)
    call fft_delete_plan(plan%py_inv)
    call fft_delete_plan(plan%pz_inv)

    call delete_layout_3D( plan%layout_x )
    call delete_layout_3D( plan%layout_y )
    call delete_layout_3D( plan%layout_z )

    SLL_DEALLOCATE_ARRAY(plan%array_x, ierr)
    SLL_DEALLOCATE_ARRAY(plan%array_y, ierr)
    SLL_DEALLOCATE_ARRAY(plan%array_z, ierr)
    SLL_DEALLOCATE(plan, ierr)
  end subroutine delete_poisson_3d_periodic_plan_par

  !> Check sizes of arrays in input
  subroutine verify_argument_sizes_par(layout, rho, phi)
    type(layout_3D), pointer   :: layout
    sll_real64, dimension(:,:,:) :: rho
    sll_real64, dimension(:,:,:) :: phi
    sll_int32,  dimension(3)     :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                    :: i

    call compute_local_sizes( layout, n(1), n(2), n(3) )

    do i=1,3
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
          print*, 'Input sizes passed to solve_poisson_3d_periodic_par do not match'
          if (i==1) then
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction x'
          elseif (i==2) then
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction y'
          else
             print*, 'Input sizes passed to "solve_poisson_3d_periodic_par" ', &
                  'do not match in direction z'
          endif
          print *, 'Exiting...'
          stop
       endif
    enddo

  end subroutine verify_argument_sizes_par


end module sll_poisson_3d_periodic_par
