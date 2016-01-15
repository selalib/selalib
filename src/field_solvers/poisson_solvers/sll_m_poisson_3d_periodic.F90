!> @ingroup poisson_solvers
!> @brief 
!> 3D poisson solver with periodic boundary conditions

module sll_m_poisson_3d_periodic

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_fft, only: &
    sll_s_fft_apply_plan_c2c_1d, &
    sll_p_fft_backward, &
    sll_s_fft_delete_plan, &
    sll_p_fft_forward, &
    sll_f_fft_new_plan_c2c_1d, &
    sll_t_fft_plan

  implicit none

  public :: &
    sll_s_delete_poisson_3d_periodic, &
    sll_f_new_poisson_3d_periodic, &
    sll_t_poisson_3d_periodic, &
    sll_s_solve_poisson_3d_periodic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Structure to solve Poisson equation on 3d domain. Mesh is cartesian and
  !> all boundary conditions are periodic. Numerical method is FFT based.
  type sll_t_poisson_3d_periodic
     sll_int32                   :: nx     !< Number of points-1 in x-direction
     sll_int32                   :: ny     !< Number of points-1 in y-direction
     sll_int32                   :: nz     !< Number of points-1 in z-direction
     sll_real64                  :: Lx     !< x length of domain
     sll_real64                  :: Ly     !< y length of domain
     sll_real64                  :: Lz     !< z length of domain
     type(sll_t_fft_plan), pointer :: px     !< forward fft plan along x
     type(sll_t_fft_plan), pointer :: py     !< forward fft plan along y
     type(sll_t_fft_plan), pointer :: pz     !< forward fft plan along z
     type(sll_t_fft_plan), pointer :: px_inv !< backward fft plan along x
     type(sll_t_fft_plan), pointer :: py_inv !< backward fft plan along y
     type(sll_t_fft_plan), pointer :: pz_inv !< backward fft plan along z
     sll_comp64, dimension(:,:,:), pointer :: hat_rho !< fft of RHS
     sll_comp64, dimension(:,:,:), pointer :: hat_phi !< fft of potential
  end type sll_t_poisson_3d_periodic

contains


  !> Allocate a structure to solve Poisson equation on 3d cartesian mesh
  !> with periodic boundary conditions
  !> @return
  function sll_f_new_poisson_3d_periodic(nx ,ny ,nz, Lx, Ly, Lz) &
                                                         result(plan)

    sll_int32                                    :: nx   !< number of points in x
    sll_int32                                    :: ny   !< number of points in y
    sll_int32                                    :: nz   !< number of points in z
    sll_comp64,                    dimension(nx) :: x
    sll_comp64,                    dimension(ny) :: y
    sll_comp64,                    dimension(nz) :: z
    sll_int32                                    :: ierr
    sll_real64                                   :: Lx   !< Length along x
    sll_real64                                   :: Ly   !< Length along y
    sll_real64                                   :: Lz   !< Length along z
    type (sll_t_poisson_3d_periodic), pointer :: plan !< Poisson solver type

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%hat_rho(nx,ny,nz), ierr)
    SLL_ALLOCATE(plan%hat_phi(nx,ny,nz), ierr)

    ! Geometry informations
    plan%nx = nx
    plan%ny = ny
    plan%nz = nz
    plan%Lx = Lx
    plan%Ly = Ly
    plan%Lz = Lz

    ! For FFTs (in each direction)
    plan%px => sll_f_fft_new_plan_c2c_1d( nx, x, x, sll_p_fft_forward )
    plan%py => sll_f_fft_new_plan_c2c_1d( ny, y, y, sll_p_fft_forward )
    plan%pz => sll_f_fft_new_plan_c2c_1d( nz, z, z, sll_p_fft_forward )

    ! For inverse FFTs (in each direction)
    plan%px_inv => sll_f_fft_new_plan_c2c_1d( nx, x, x, sll_p_fft_backward )
    plan%py_inv => sll_f_fft_new_plan_c2c_1d( ny, y, y, sll_p_fft_backward )
    plan%pz_inv => sll_f_fft_new_plan_c2c_1d( nz, z, z, sll_p_fft_backward )

  end function sll_f_new_poisson_3d_periodic

  !> Compute the potential from 3d Poisson solver
  subroutine sll_s_solve_poisson_3d_periodic(plan, rho, phi)

    type (sll_t_poisson_3d_periodic), pointer :: plan !< Solver structure
    sll_real64, dimension(:,:,:)                 :: rho  !< charge density
    sll_real64, dimension(:,:,:)                 :: phi  !< Electric potential
    sll_int32                                    :: nx, ny, nz
    sll_int32                                    :: i, j, k
    sll_real64                                   :: Lx, Ly, Lz
    sll_real64                                   :: ind_x, ind_y, ind_z
    logical, save                                :: flag = .true.

    ! Checking input arguments consistency
    if (flag) then
       call verify_argument_sizes(plan, rho, phi)
       flag = .false.
    end if

    nx = plan%nx
    ny = plan%ny
    nz = plan%nz
    Lx = plan%Lx
    Ly = plan%Ly
    Lz = plan%Lz

    ! FFTs in x-direction
    plan%hat_rho = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz
       do j=1,ny
          call sll_s_fft_apply_plan_c2c_1d( plan%px, plan%hat_rho(:,j,k), plan%hat_rho(:,j,k) )
       enddo
    enddo

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call sll_s_fft_apply_plan_c2c_1d( plan%py, plan%hat_rho(i,:,k), plan%hat_rho(i,:,k) )
       enddo
    enddo

    ! FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call sll_s_fft_apply_plan_c2c_1d( plan%pz, plan%hat_rho(i,j,:), plan%hat_rho(i,j,:) )
       enddo
    enddo

    plan%hat_rho = plan%hat_rho/(nx*ny*nz)

    ! Compute hat_phi, phi = inv_fft(hat_phi)
    do k=1,nz
       do j=1,ny
          do i=1,nx
             if (i<=nx/2) then
                ind_x = real(i-1,f64)
             else
                ind_x = real(nx-(i-1),f64)
             endif
             if (j<=ny/2) then
                ind_y = real(j-1,f64)
             else
                ind_y = real(ny-(j-1),f64)
             endif
             if (k<=nz/2) then
                ind_z = real(k-1,f64)
             else
                ind_z = real(nz-(k-1),f64)
             endif
             if ( (ind_x==0) .and. (ind_y==0) .and. (ind_z==0) ) then
                plan%hat_phi(i,j,k) = (0._f64,0._f64)
             else
                plan%hat_phi(i,j,k) = &
                    plan%hat_rho(i,j,k)/(4*sll_p_pi**2*((ind_x/Lx)**2 &
                    + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call sll_s_fft_apply_plan_c2c_1d( plan%pz_inv, plan%hat_phi(i,j,:), plan%hat_phi(i,j,:) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call sll_s_fft_apply_plan_c2c_1d( plan%py_inv, plan%hat_phi(i,:,k), plan%hat_phi(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call sll_s_fft_apply_plan_c2c_1d( plan%px_inv, plan%hat_phi(:,j,k), plan%hat_phi(:,j,k) )
       enddo
    enddo

    phi = real(plan%hat_phi, f64)

  end subroutine sll_s_solve_poisson_3d_periodic


  !> Delete the 3d poisson solver object
  subroutine sll_s_delete_poisson_3d_periodic(plan)

    type(sll_t_poisson_3d_periodic), pointer :: plan !< Poisson solver object
    sll_int32                                    :: ierr

    ! Fixme: some error checking, whether the poisson pointer is associated
    ! for instance
    SLL_ASSERT( associated(plan) )

    call sll_s_fft_delete_plan(plan%px)
    call sll_s_fft_delete_plan(plan%py)
    call sll_s_fft_delete_plan(plan%pz)

    call sll_s_fft_delete_plan(plan%px_inv)
    call sll_s_fft_delete_plan(plan%py_inv)
    call sll_s_fft_delete_plan(plan%pz_inv)

    SLL_DEALLOCATE_ARRAY(plan%hat_rho, ierr)
    SLL_DEALLOCATE_ARRAY(plan%hat_phi, ierr)
    SLL_DEALLOCATE(plan, ierr)

  end subroutine sll_s_delete_poisson_3d_periodic

  !> Subroutine to check that arrays are compatible with the solver
  subroutine verify_argument_sizes(plan, rho, phi)

    type (sll_t_poisson_3d_periodic), pointer :: plan
    sll_real64, dimension(:,:,:)                 :: rho
    sll_real64, dimension(:,:,:)                 :: phi
    sll_int32,  dimension(3)                     :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                                    :: i

    n(1) = plan%nx
    n(2) = plan%ny
    n(3) = plan%nz

    do i=1,3
       if ( (n(i)/=size(rho,i)) .or. (n(i)/=size(phi,i))  ) then
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

  end subroutine verify_argument_sizes

end module sll_m_poisson_3d_periodic