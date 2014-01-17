!**************************************************************
!  Copyright INRIA, CNRS
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

#include "sll_fftw.h"

#define D_DX(field)                                                   \
call fftw_execute_dft_r2c(plan%fwx, field, plan%fft_x_array);         \
plan%fft_x_array = -cmplx(0.0_f64,plan%kx,kind=f64)*plan%fft_x_array; \
call fftw_execute_dft_c2r(plan%bwx, plan%fft_x_array, plan%d_dx);     \
plan%d_dx = plan%d_dx / plan%ncx

#define D_DY(field)                                                   \
call fftw_execute_dft_r2c(plan%fwy, field, plan%fft_y_array);         \
plan%fft_y_array = -cmplx(0.0_f64,plan%ky,kind=f64)*plan%fft_y_array; \
call fftw_execute_dft_c2r(plan%bwy, plan%fft_y_array, plan%d_dy);     \
plan%d_dy = plan%d_dy / plan%ncy

#define MPI_MASTER 0

!> @brief 
!> Selalib periodic 2D maxwell solver for cartesian coordinates.
!>   
module sll_maxwell_2d_periodic_cartesian_par
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_constants.h"

use sll_remapper
use sll_collective
use fftw3

implicit none

!> Maxwell solver 2D object, PSTD scheme
type maxwell_2d_periodic_plan_cartesian_par
   sll_int32                           :: ncx         !< number of cells x
   sll_int32                           :: ncy         !< number of cells y
   fftw_plan                           :: fwx         !< fftw plan forward x
   fftw_plan                           :: fwy         !< fftw plan forward y
   fftw_plan                           :: bwx         !< fftw plan backward x
   fftw_plan                           :: bwy         !< fftw plan backward y
   sll_real64, dimension(:), pointer   :: d_dx        !< x derivative
   sll_real64, dimension(:), pointer   :: d_dy        !< y derivative
   fftw_comp,  dimension(:), pointer   :: fft_x_array !< fft x transform
   fftw_comp,  dimension(:), pointer   :: fft_y_array !< fft y transform
   sll_real64                          :: e_0         !< electric conductibility
   sll_real64                          :: mu_0        !< magnetic permeability
   fftw_plan                           :: p_x_array   !< x pointer to memory
   fftw_plan                           :: p_y_array   !< y pointer to memory
   type(layout_2D),  pointer           :: layout_x    !< layout sequential in x
   type(layout_2D),  pointer           :: layout_y    !< layout sequential in y
   type(remap_plan_2D_real64), pointer :: rmp_xy      !< remap x->y pointer
   type(remap_plan_2D_real64), pointer :: rmp_yx      !< remap y->x pointer
   sll_real64, dimension(:,:), pointer :: fz_x        !< array sequential in x
   sll_real64, dimension(:,:), pointer :: fz_y        !< array sequential in y
   sll_real64, dimension(:), pointer   :: kx          !< x wave number
   sll_real64, dimension(:), pointer   :: ky          !< y wave number
end type maxwell_2d_periodic_plan_cartesian_par

sll_int32, private :: i
sll_int32, private :: j

contains

  !> Presently, this function receives the geometric information as 
  !> individual arguments. We should consider passing the 'simple geometry'
  !> object that we have for the cartesian cases.
  !> @return
  function new_maxwell_2d_periodic_plan_cartesian_par( &
    layout_x, layout_y, ncx, ncy, Lx, Ly ) result(plan)

    !> maxwell object
    type(maxwell_2d_periodic_plan_cartesian_par), pointer :: plan 
    type(layout_2D), pointer :: layout_x !< sequential in x direction
    type(layout_2D), pointer :: layout_y !< sequential in y direction
    sll_int32                :: ncx      !< x cell number
    sll_int32                :: ncy      !< y cell number
    sll_real64               :: Lx       !< Domain x length
    sll_real64               :: Ly       !< Domain y length
    sll_int64                :: prank    !< processor rank
    sll_int64                :: psize    !< processor size
    sll_int32                :: error    !< error code
    sll_int32                :: nx_loc   !< x local points
    sll_int32                :: ny_loc   !< y local points
    fftw_int                 :: sz_x_array
    fftw_int                 :: sz_y_array
    sll_real64               :: kx0
    sll_real64               :: ky0

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

    if ( (.not.is_power_of_two(int(ncx,i64))) .and. &
         (.not.is_power_of_two(int(ncy,i64))) ) then     
       print *, 'This test needs to run with numbers of cells which are',  &
                'powers of 2.'
       print *, 'Exiting...'
       stop
    end if

    SLL_ALLOCATE(plan,                 error)
    SLL_CLEAR_ALLOCATE(plan%d_dx(1:ncx), error)
    SLL_CLEAR_ALLOCATE(plan%d_dy(1:ncy), error)

    plan%ncx = ncx
    plan%ncy = ncy

    FFTW_ALLOCATE(plan%fft_x_array,ncx/2+1,sz_x_array,plan%p_x_array)
    FFTW_ALLOCATE(plan%fft_y_array,ncy/2+1,sz_y_array,plan%p_y_array)

    NEW_FFTW_PLAN_R2C_1D(plan%fwx,ncx,plan%d_dx,plan%fft_x_array)
    NEW_FFTW_PLAN_C2R_1D(plan%bwx,ncx,plan%fft_x_array,plan%d_dx)
    NEW_FFTW_PLAN_R2C_1D(plan%fwy,ncy,plan%d_dy,plan%fft_y_array)
    NEW_FFTW_PLAN_C2R_1D(plan%bwy,ncy,plan%fft_y_array,plan%d_dy)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)

    ! Layout and local sizes for FFTs in x-direction
    plan%layout_x => layout_x
    call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
    SLL_CLEAR_ALLOCATE(plan%fz_x(1:nx_loc,1:ny_loc),error)

    ! Layout and local sizes for FFTs in y-direction
    plan%layout_y => layout_y
    call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
    SLL_CLEAR_ALLOCATE(plan%fz_y(1:nx_loc,1:ny_loc),error)

    plan%rmp_xy => new_remap_plan(plan%layout_x, plan%layout_y, plan%fz_x)
    plan%rmp_yx => new_remap_plan(plan%layout_y, plan%layout_x, plan%fz_y)

    SLL_ALLOCATE(plan%kx(ncx/2+1), error)
    SLL_ALLOCATE(plan%ky(ncy/2+1), error)
   
    kx0 = 2._f64*sll_pi/Lx
    ky0 = 2._f64*sll_pi/Ly

    do i=2,ncx/2+1
       plan%kx(i) = (i-1)*kx0
    end do
    plan%kx(1) = 1.0_f64
    do j=2,ncy/2+1
       plan%ky(j) = (j-1)*ky0
    end do
    plan%ky(1) = 1.0_f64

  end function new_maxwell_2d_periodic_plan_cartesian_par

!********************************************************************************

  !> Solve faraday equation (TE mode)
  subroutine faraday_te(plan,dt,ex,ey)

    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    sll_real64, dimension(:,:) :: ex       !< x electric field
    sll_real64, dimension(:,:) :: ey       !< y electric field
    sll_int32                  :: ncx      !< global x cell number
    sll_int32                  :: ncy      !< global y cell number
    sll_int32                  :: nx_loc   !< local  x cell number
    sll_int32                  :: ny_loc   !< local  y cell number
    sll_int32                  :: prank
    sll_int64                  :: psize
    sll_real64, intent(in)     :: dt       !< time step
    sll_real64                 :: dt_mu

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

    ncx  = plan%ncx
    ncy  = plan%ncy

#ifdef DEBUG
    call verify_argument_sizes_par(plan%layout_x, ey)
    call verify_argument_sizes_par(plan%layout_y, ex)
#endif

    dt_mu = dt / plan%mu_0 

    call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
    do j = 1, ny_loc
       D_DX(ey(:,j))
       plan%fz_x(:,j) = plan%fz_x(:,j) - dt_mu * plan%d_dx
    end do

    call apply_remap_2D( plan%rmp_xy,plan%fz_x,plan%fz_y)

    call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
    do i = 1, nx_loc
       D_DY(ex(i,:))
       plan%fz_y(i,:) = plan%fz_y(i,:) + dt_mu * plan%d_dy
    end do
      
  end subroutine faraday_te

  !> Solve Ampere-Maxwell equation (TE mode)
  subroutine ampere_te(plan,dt,ex,ey,jx,jy)

    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    sll_real64, dimension(:,:), intent(inout) :: ex     !< Ex field
    sll_real64, dimension(:,:), intent(inout) :: ey     !< Ey field
    sll_real64, dimension(:,:), optional      :: jx     !< Jx field
    sll_real64, dimension(:,:), optional      :: jy     !< Jy field
    sll_real64, intent(in)                    :: dt     !< time step
    sll_int32                                 :: nx_loc !< local  x cell number
    sll_int32                                 :: ny_loc !< local  y cell number
    sll_int32                                 :: prank
    sll_int64                                 :: psize
    sll_real64                                :: dt_e

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

#ifdef DEBUG
    call verify_argument_sizes_par(plan%layout_x, ey)
    call verify_argument_sizes_par(plan%layout_y, ex)
#endif

    dt_e = dt / plan%e_0

    call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
    do i = 1, nx_loc
       D_DY(plan%fz_y(i,:))
       ex(i,:) = ex(i,:) + dt_e * plan%d_dy
    end do

    if (present(jx)) then
#ifdef DEBUG
       call verify_argument_sizes_par(plan%layout_y, jx)
#endif
       ex(:,:) = ex(:,:) - dt_e * jx(:,:) 
    end if

    call apply_remap_2D( plan%rmp_xy,plan%fz_x,plan%fz_y)

    call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
    do j = 1, ny_loc
       D_DX(plan%fz_x(:,j))
       ey(:,j) = ey(:,j) - dt_e * plan%d_dx
    end do

    if (present(jy)) then
#ifdef DEBUG
       call verify_argument_sizes_par(plan%layout_x, jy)
#endif
       ey(:,:) = ey(:,:) - dt * jy(:,:) / plan%e_0
    endif


  end subroutine ampere_te

  !> Solve Ampere-Maxwell equation (TE mode)
  subroutine ampere_tm(plan,dt,bx,by,jz)

    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    sll_real64, dimension(:,:) :: bx           !< electric field sequential in y
    sll_real64, dimension(:,:) :: by           !< electric field sequential in x
    sll_real64, dimension(:,:), optional :: jz !< current density sequential in 	
    sll_int32                  :: ncx          !< global x cell number
    sll_int32                  :: ncy          !< global y cell number
    sll_int32                  :: nx_loc       !< local  x cell number
    sll_int32                  :: ny_loc       !< local  y cell number
    sll_int32                  :: prank
    sll_int64                  :: psize
    sll_real64, intent(in)     :: dt           !< time step
    sll_real64                 :: dt_e

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

    ncx  = plan%ncx
    ncy  = plan%ncy

#ifdef DEBUG
    call verify_argument_sizes_par(plan%layout_x, by)
    call verify_argument_sizes_par(plan%layout_y, bx)
#endif

    dt_e = dt / plan%e_0

    call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
    do j = 1, ny_loc
      D_DX(by(:,j))
      plan%fz_x(:,j) = plan%fz_x(:,j) + dt_e * plan%d_dx
    end do

    call apply_remap_2D( plan%rmp_xy,plan%fz_x,plan%fz_y)

    call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
    do i = 1, nx_loc
      D_DY(bx(i,:))
      plan%fz_y(i,:) = plan%fz_y(i,:) - dt_e * plan%d_dy
    end do

    if (present(jz)) then
#ifdef DEBUG
    call verify_argument_sizes_par(plan%layout_y, jz)
#endif
      plan%fz_y = plan%fz_y - dt_e * jz 
    end if
      
  end subroutine ampere_tm

  !> Solve Ampere-Maxwell equation (TE mode)
  subroutine faraday_tm(plan,dt,bx,by)

    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    sll_real64, dimension(:,:), intent(inout) :: bx     !< Bx field
    sll_real64, dimension(:,:), intent(inout) :: by     !< By field
    sll_real64, intent(in)                    :: dt     !< time step
    sll_int32                                 :: nx_loc !< local  x cell number
    sll_int32                                 :: ny_loc !< local  y cell number
    sll_int32                                 :: prank
    sll_int64                                 :: psize
    sll_real64                                :: dt_mu

    prank = sll_get_collective_rank( sll_world_collective )
    psize = sll_get_collective_size( sll_world_collective )

#ifdef DEBUG
    call verify_argument_sizes_par(plan%layout_x, by)
    call verify_argument_sizes_par(plan%layout_y, bx)
#endif

    dt_mu = dt / plan%mu_0

    call compute_local_sizes_2d(plan%layout_y,nx_loc,ny_loc)
    do i = 1, nx_loc
       D_DY(plan%fz_y(i,:))
       bx(i,:) = bx(i,:) - dt_mu * plan%d_dy
    end do

    call apply_remap_2D( plan%rmp_xy,plan%fz_x,plan%fz_y)

    call compute_local_sizes_2d(plan%layout_x,nx_loc,ny_loc)
    do j = 1, ny_loc
       D_DX(plan%fz_x(:,j))
       by(:,j) = by(:,j) + dt_mu * plan%d_dx
    end do

  end subroutine faraday_tm

  !> Delete maxwell solver object

  !> Delete maxwell solver object
  subroutine delete_maxwell_2d_periodic_plan_cartesian_par(plan)
    type (maxwell_2d_periodic_plan_cartesian_par), pointer :: plan !< maxwell object

    if( .not. associated(plan) ) then
       print *, 'ERROR, delete_maxwell_3d_periodic_plan_par(): ', &
            'passed plan is not associated.'
       STOP
    end if

#ifdef FFTW_F2003
    if (c_associated(plan%p_x_array)) call fftw_free(plan%p_x_array)
    if (c_associated(plan%p_y_array)) call fftw_free(plan%p_y_array)
#endif

    call fftw_destroy_plan(plan%fwx)
    call fftw_destroy_plan(plan%fwy)
    call fftw_destroy_plan(plan%bwx)
    call fftw_destroy_plan(plan%bwy)

    call delete(plan%rmp_xy)
    call delete(plan%rmp_yx)
    call delete(plan%layout_x)
    call delete(plan%layout_y)

  end subroutine delete_maxwell_2d_periodic_plan_cartesian_par

  !> Check array size.
  subroutine verify_argument_sizes_par(layout, array)
    type(layout_2D), pointer       :: layout !< layout
    sll_real64, dimension(:,:)     :: array  !< array
    sll_int32,  dimension(2)       :: n      !< array dimension
    sll_int32                      :: i

    call compute_local_sizes_2d( layout, n(1), n(2) )

    do i=1,2
       if (n(i)/=size(array,i)) then
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
