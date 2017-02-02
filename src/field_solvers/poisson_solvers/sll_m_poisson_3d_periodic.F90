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
    sll_s_fft_exec_c2c_1d, &
    sll_p_fft_backward, &
    sll_s_fft_free, &
    sll_p_fft_forward, &
    sll_s_fft_init_c2c_1d, &
    sll_t_fft

  implicit none

  public :: &
    sll_s_poisson_3d_periodic_init, &
    sll_s_poisson_3d_periodic_free, &
    sll_t_poisson_3d_periodic, &
    sll_s_poisson_3d_periodic_solve

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
     type(sll_t_fft)        :: px     !< forward fft plan along x
     type(sll_t_fft)        :: py     !< forward fft plan along y
     type(sll_t_fft)        :: pz     !< forward fft plan along z
     type(sll_t_fft)        :: px_inv !< backward fft plan along x
     type(sll_t_fft)        :: py_inv !< backward fft plan along y
     type(sll_t_fft)        :: pz_inv !< backward fft plan along z
     sll_comp64, dimension(:,:,:), pointer :: hat_rho !< fft of RHS
     sll_comp64, dimension(:,:,:), pointer :: hat_phi !< fft of potential
  end type sll_t_poisson_3d_periodic

contains


  !> Allocate a structure to solve Poisson equation on 3d cartesian mesh
  !> with periodic boundary conditions
  !> @return
  subroutine sll_s_poisson_3d_periodic_init(self, nx ,ny ,nz, Lx, Ly, Lz) 

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
    type (sll_t_poisson_3d_periodic)             :: self !< Poisson solver type

    SLL_ALLOCATE(self%hat_rho(nx,ny,nz), ierr)
    SLL_ALLOCATE(self%hat_phi(nx,ny,nz), ierr)

    ! Geometry informations
    self%nx = nx
    self%ny = ny
    self%nz = nz
    self%Lx = Lx
    self%Ly = Ly
    self%Lz = Lz

    ! For FFTs (in each direction)
    call sll_s_fft_init_c2c_1d( self%px, nx, x, x, sll_p_fft_forward )
    call sll_s_fft_init_c2c_1d( self%py, ny, y, y, sll_p_fft_forward )
    call sll_s_fft_init_c2c_1d( self%pz, nz, z, z, sll_p_fft_forward )

    ! For inverse FFTs (in each direction)
    call sll_s_fft_init_c2c_1d( self%px_inv, nx, x, x, sll_p_fft_backward )
    call sll_s_fft_init_c2c_1d( self%py_inv, ny, y, y, sll_p_fft_backward )
    call sll_s_fft_init_c2c_1d( self%pz_inv, nz, z, z, sll_p_fft_backward )

  end subroutine sll_s_poisson_3d_periodic_init

  !> Compute the potential from 3d Poisson solver
  subroutine sll_s_poisson_3d_periodic_solve(self, rho, phi)

    type (sll_t_poisson_3d_periodic) :: self !< Solver structure
    sll_real64, dimension(:,:,:)     :: rho  !< charge density
    sll_real64, dimension(:,:,:)     :: phi  !< Electric potential
    sll_int32                        :: nx, ny, nz
    sll_int32                        :: i, j, k
    sll_real64                       :: Lx, Ly, Lz
    sll_real64                       :: ind_x, ind_y, ind_z
    logical, save                    :: flag = .true.

    ! Checking input arguments consistency
    if (flag) then
       call verify_argument_sizes(self, rho, phi)
       flag = .false.
    end if

    nx = self%nx
    ny = self%ny
    nz = self%nz
    Lx = self%Lx
    Ly = self%Ly
    Lz = self%Lz

    ! FFTs in x-direction
    self%hat_rho = cmplx(rho, 0_f64, kind=f64)
    do k=1,nz
       do j=1,ny
          call sll_s_fft_exec_c2c_1d( self%px, self%hat_rho(:,j,k), self%hat_rho(:,j,k) )
       enddo
    enddo

    ! FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call sll_s_fft_exec_c2c_1d( self%py, self%hat_rho(i,:,k), self%hat_rho(i,:,k) )
       enddo
    enddo

    ! FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call sll_s_fft_exec_c2c_1d( self%pz, self%hat_rho(i,j,:), self%hat_rho(i,j,:) )
       enddo
    enddo

    self%hat_rho = self%hat_rho/(nx*ny*nz)

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
                self%hat_phi(i,j,k) = (0._f64,0._f64)
             else
                self%hat_phi(i,j,k) = &
                    self%hat_rho(i,j,k)/(4.0_f64*sll_p_pi**2*((ind_x/Lx)**2 &
                    + (ind_y/Ly)**2+(ind_z/Lz)**2))
             endif
          enddo
       enddo
    enddo

    ! Inverse FFTs in z-direction
    do j=1,ny
       do i=1,nx
          call sll_s_fft_exec_c2c_1d( self%pz_inv, self%hat_phi(i,j,:), self%hat_phi(i,j,:) )
       enddo
    enddo

    ! Inverse FFTs in y-direction
    do k=1,nz
       do i=1,nx
          call sll_s_fft_exec_c2c_1d( self%py_inv, self%hat_phi(i,:,k), self%hat_phi(i,:,k) )
       enddo
    enddo

    ! Inverse FFTs in x-direction
    do k=1,nz
       do j=1,ny
          call sll_s_fft_exec_c2c_1d( self%px_inv, self%hat_phi(:,j,k), self%hat_phi(:,j,k) )
       enddo
    enddo

    phi = real(self%hat_phi, f64)

  end subroutine sll_s_poisson_3d_periodic_solve


  !> Delete the 3d poisson solver object
  subroutine sll_s_poisson_3d_periodic_free(self)

    type(sll_t_poisson_3d_periodic) :: self !< Poisson solver object
    sll_int32                       :: ierr

    call sll_s_fft_free(self%px)
    call sll_s_fft_free(self%py)
    call sll_s_fft_free(self%pz)

    call sll_s_fft_free(self%px_inv)
    call sll_s_fft_free(self%py_inv)
    call sll_s_fft_free(self%pz_inv)

    SLL_DEALLOCATE_ARRAY(self%hat_rho, ierr)
    SLL_DEALLOCATE_ARRAY(self%hat_phi, ierr)

  end subroutine sll_s_poisson_3d_periodic_free

  !> Subroutine to check that arrays are compatible with the solver
  subroutine verify_argument_sizes(self, rho, phi)

    type (sll_t_poisson_3d_periodic) :: self
    sll_real64, dimension(:,:,:)     :: rho
    sll_real64, dimension(:,:,:)     :: phi
    sll_int32,  dimension(3)         :: n ! nx_loc, ny_loc, nz_loc
    sll_int32                        :: i

    n(1) = self%nx
    n(2) = self%ny
    n(3) = self%nz

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
