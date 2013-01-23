!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Pierre Navaro 
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


#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)  \
sz_array = int((array_size/2+1),C_SIZE_T);                \
p_array = fftw_alloc_complex(sz_array);                   \
call c_f_pointer(p_array, array, [array_size/2+1])        \

#define D_DX(field)                                           \
call fftw_execute_dft_r2c(plan%fwx, field, plan%fft_x_array);       \
plan%fft_x_array = -cmplx(0.0_f64,plan%kx,kind=f64)*plan%fft_x_array;     \
call fftw_execute_dft_c2r(plan%bwx, plan%fft_x_array, plan%d_dx);   \
plan%d_dx = plan%d_dx / nx

#define D_DY(field)                                           \
call fftw_execute_dft_r2c(plan%fwy, field, plan%fft_y_array);       \
plan%fft_y_array = -cmplx(0.0_f64,plan%ky,kind=f64)*plan%fft_y_array;     \
call fftw_execute_dft_c2r(plan%bwy, plan%fft_y_array, plan%d_dy);    \
plan%d_dy = plan%d_dy / ny

#define MPI_MASTER 0

!> @brief 
!> Selalib periodic 2D maxwell solver for cartesian coordinates.
!   
!> @author                    
!> Pierre Navaro 
module sll_maxwell_2d_periodic_cartesian_par
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"

  use, intrinsic :: iso_c_binding
  use remapper
  use sll_fft
  use numeric_constants
  use sll_collective
  use sll_maxwell

  implicit none

  !> solve subroutine parameter to compute electric field
  sll_int32, public, parameter :: FARADAY = 0  
  !> solve subroutine parameter to compute magnetic field
  sll_int32, public, parameter :: AMPERE  = 1

  !> Maxwell solver object
  type maxwell_2d_periodic_plan_cartesian_par
     sll_int32                           :: ncx         !< number of cells x
     sll_int32                           :: ncy         !< number of cells y
     sll_real64                          :: Lx          !< domain length 
     sll_real64                          :: Ly          !< domain length
     type(C_PTR), pointer                :: fwx         !< fftw plan forward x
     type(C_PTR), pointer                :: fwy         !< fftw plan forward y
     type(C_PTR), pointer                :: bwx         !< fftw plan backward x
     type(C_PTR), pointer                :: bwy         !< fftw plan backward y
     type(layout_2D),  pointer           :: layout_x    !< layout sequential in x
     type(layout_2D),  pointer           :: layout_y    !< layout sequential in y
     sll_int32                           :: npx_loc     !< local points number x
     sll_int32                           :: npy_loc     !< local points number y
     sll_real64, dimension(:), pointer   :: d_dx        !< x derivative
     sll_real64, dimension(:), pointer   :: d_dy        !< y derivative
     sll_comp64, dimension(:), pointer   :: fft_x_array !< fft x transform
     sll_comp64, dimension(:), pointer   :: fft_y_array !< fft y transform
     type(remap_plan_2D_comp64), pointer :: rmp_xy      !< remap x->y pointer
     type(remap_plan_2D_comp64), pointer :: rmp_yx      !< remap y->x pointer
     sll_real64                          :: mu_0        !< magnetic permeability
     type(C_PTR)                         :: p_x_array   !< x pointer to memory
     type(C_PTR)                         :: p_y_array   !< y pointer to memory
  end type maxwell_2d_periodic_plan_cartesian_par

include 'fftw3.f03'

contains

  !> Presently, this function receives the geometric information as 
  !> individual arguments. We should consider passing the 'simple geometry'
  !> object that we have for the cartesian cases.
  function new_maxwell_2d_periodic_plan_cartesian_par( &
    start_layout, &
    ncx, &
    ncy, &
    Lx, &
    Ly ) result(plan)

    type(maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object
    type(layout_2D), pointer                     :: start_layout  !< layout
    sll_int32                                    :: ncx           !< x cell number
    sll_int32                                    :: ncy           !< y cell number
    sll_real64                                   :: Lx            !< domain x size
    sll_real64                                   :: Ly            !< domain y size
    sll_int64                                    :: prank         !< processor rank
    sll_int64                                    :: psize         !< processor size
    type(sll_collective_t), pointer              :: collective    !< mpi object
    sll_int64                                    :: nprocx        !< procs number x
    sll_int64                                    :: nprocy        !< procs number y
    sll_int32                                    :: ierr          !< error code
    sll_int32                                    :: npx_loc       !< x local points
    sll_int32                                    :: npy_loc

    integer(C_SIZE_T) :: sz_x_array
    integer(C_SIZE_T) :: sz_y_array

    ! The collective to be used is the one that comes with the given layout.
    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

    if ( (.not.is_power_of_two(int(ncx,i64))) .and. &
         (.not.is_power_of_two(int(ncy,i64))) ) then     
       print *, 'This test needs to run with numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan%d_dx(ncx), ierr)
    SLL_ALLOCATE(plan%d_dy(ncy), ierr)
    SLL_ALLOCATE(plan, ierr)

    ! Geometry
    plan%ncx = ncx
    plan%ncy = ncy
    plan%Lx  = Lx
    plan%Ly  = Ly

    FFTW_ALLOCATE(plan%fft_x_array,ncx/2+1,sz_x_array,plan%p_x_array)
    FFTW_ALLOCATE(plan%fft_y_array,ncy/2+1,sz_y_array,plan%p_y_array)

    plan%fwx = fftw_plan_dft_r2c_1d(ncx, plan%d_dx,  plan%fft_x_array, FFTW_ESTIMATE)
    plan%bwx = fftw_plan_dft_c2r_1d(ncx, plan%fft_x_array, plan%d_dx,  FFTW_ESTIMATE)
    plan%fwy = fftw_plan_dft_r2c_1d(ncy, plan%d_dy,  plan%fft_y_array, FFTW_ESTIMATE)
    plan%bwy = fftw_plan_dft_c2r_1d(ncy, plan%fft_y_array, plan%d_dy,  FFTW_ESTIMATE)

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_x => start_layout
    call compute_local_sizes_2d( plan%layout_x, npx_loc, npy_loc )
    plan%npy_loc = npy_loc

    ! Layout and local sizes for FFTs in y-direction
    plan%layout_y => new_layout_2D( collective )
    nprocx = psize
    nprocy = 1

    call initialize_layout_with_distributed_2D_array( &
         ncx, &
         ncy, &
         int(nprocx,i32), &
         int(nprocy,i32), &
         plan%layout_y )

    call compute_local_sizes_2d(plan%layout_y,npx_loc,npy_loc )
    plan%npx_loc = npx_loc

    !plan%rmp_xy => new_remap_plan(plan%layout_x, plan%layout_y, plan%fft_x_array)
    !plan%rmp_yx => new_remap_plan(plan%layout_y, plan%layout_x, plan%fft_y_array)

  end function new_maxwell_2d_periodic_plan_cartesian_par

!********************************************************************************

  !> Solve maxwell equations
  subroutine solve_maxwell_2d_periodic_cartesian_par(plan, dt, fx, fy, fz, equation)

    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    sll_real64, dimension(:,:), target            :: fx !< Ex or Bx
    sll_real64, dimension(:,:), target            :: fy !< Ey or By
    sll_real64, dimension(:,:), target            :: fz !< Bz or Ez
    sll_real64, dimension(:,:), pointer           :: ex
    sll_real64, dimension(:,:), pointer           :: ey
    sll_real64, dimension(:,:), pointer           :: bz
    sll_int32, intent(in)                         :: equation !< ampere-maxwell or faraday
    ! global sizes
    sll_int32                                     :: ncx     !< global x cell number
    sll_int32                                     :: ncy     !< global y cell number
    sll_int32                                     :: npx_loc !< local  x cell number
    sll_int32                                     :: npy_loc !< local  y cell number
    sll_int32                                     :: i
    sll_int32                                     :: j
    sll_int32                                     :: ierr
    ! Reciprocals of domain lengths.
    sll_real64                                    :: r_Lx, r_Ly
    sll_real64                                    :: kx, ky
    sll_comp64                                    :: val
    sll_real64                                    :: normalization
    sll_int32                                     :: prank
    sll_int64                                     :: psize
    type(layout_2D), pointer                      :: layout_x
    type(layout_2D), pointer                      :: layout_y
    sll_int32, dimension(1:2)                     :: global
    sll_int32                                     :: gi, gj
    sll_real64, intent(in)                        :: dt     !< time step
    sll_real64                                    :: dt_mu

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

    ex => fx
    ey => fy
    bz => fz

    ncx  = plan%ncx
    ncy  = plan%ncy
    r_Lx = 1.0_f64/plan%Lx
    r_Ly = 1.0_f64/plan%Ly
    ! Get layouts to compute FFTs (in each direction)
    layout_x => plan%layout_x
    layout_y => plan%layout_y
    call verify_argument_sizes_par(layout_x, ex, ey, bz)

    ! FFTs in x-direction
    npx_loc = plan%npx_loc 
    npy_loc = plan%npy_loc 

    dt_mu = dt / plan%mu_0 

    do i = 1, npx_loc
      global = local_to_global_2D( layout_y, (/i, j/))
      gi = global(1)
      gj = global(2)
      !D_DY(ex(i,1:ny))
      !bz(i,1:ny) = hz(i,1:ny) + dt * plan%d_dy
    end do

    do j = 1, npy_loc
      !D_DX(ey(1:nx,j))
      !bz(1:nx,j) = bz(1:nx,j) - dt * plan%d_dx
    end do

    !call apply_remap_2D( plan%rmp_xy, plan%fft_x_array, plan%fft_y_array )

    !call apply_remap_2D( plan%rmp_yx, plan%fft_y_array, plan%fft_x_array )

  end subroutine solve_maxwell_2d_periodic_cartesian_par


  !> Delete maxwell solver object
  subroutine delete_maxwell_2d_periodic_plan_cartesian_par(plan)
    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan
    sll_int32                                              :: ierr

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_maxwell_3d_periodic_plan_par(): ', &
            'passed plan is not associated.'
       STOP
    end if

    if (c_associated(plan%p_x_array)) call fftw_free(plan%p_x_array)
    if (c_associated(plan%p_y_array)) call fftw_free(plan%p_y_array)
    call dfftw_destroy_plan(plan%fwx)
    call dfftw_destroy_plan(plan%fwy)
    call dfftw_destroy_plan(plan%bwx)
    call dfftw_destroy_plan(plan%bwy)

    call delete( plan%layout_x )
    call delete( plan%layout_y )
    SLL_DEALLOCATE_ARRAY(plan%fft_x_array, ierr)
    SLL_DEALLOCATE_ARRAY(plan%fft_y_array, ierr)
    call delete( plan%rmp_xy )
    call delete( plan%rmp_yx )
    SLL_DEALLOCATE(plan, ierr)

  end subroutine delete_maxwell_2d_periodic_plan_cartesian_par

  !> Check array size.
  subroutine verify_argument_sizes_par(layout, fx, fy, fz)
    type(layout_2D), pointer       :: layout
    sll_real64, dimension(:,:)     :: fx
    sll_real64, dimension(:,:)     :: fy
    sll_real64, dimension(:,:)     :: fz
    sll_int32,  dimension(2)       :: n ! nx_loc, ny_loc
    sll_int32                      :: i

    call compute_local_sizes_2d( layout, n(1), n(2) )

    do i=1,2
       if (n(i)/=size(fx,i) .or. n(i)/=size(fy,i) .or. n(i)/=size(fz,i)) then
          print*, 'ERROR: solve_maxwell_2d_periodic_cartesian_par()', &
               'size of either ex,ey or bz does not match expected size. '
          if (i==1) then
             print*, 'solve_maxwell_2d_periodic_cartesian_par(): ', &
                  'mismatch in direction x'
          elseif (i==2) then
             print*, 'solve_maxwell_2d_periodic_cartesian_par(): ', &
                  'mismatch in direction y'
          endif
          print *, 'Exiting...'
          stop
       endif
    enddo
  end subroutine verify_argument_sizes_par

end module sll_maxwell_2d_periodic_cartesian_par
