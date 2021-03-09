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
      sll_s_fft_exec_c2c_3d, &
      sll_p_fft_backward, &
      sll_s_fft_free, &
      sll_p_fft_forward, &
      sll_s_fft_init_c2c_3d, &
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
      sll_int32                             :: nx      !< Number of points in x-direction
      sll_int32                             :: ny      !< Number of points in y-direction
      sll_int32                             :: nz      !< Number of points in z-direction
      sll_real64                            :: Lx      !< x length of domain
      sll_real64                            :: Ly      !< y length of domain
      sll_real64                            :: Lz      !< z length of domain
      type(sll_t_fft)                       :: fw      !< forward fft
      type(sll_t_fft)                       :: bw      !< backward fft plan along x
      sll_comp64, dimension(:, :, :), pointer :: hat_rho !< fft of RHS
      sll_comp64, dimension(:, :, :), pointer :: hat_phi !< fft of potential
   end type sll_t_poisson_3d_periodic

contains

   !> Allocate a structure to solve Poisson equation on 3d cartesian mesh
   !> with periodic boundary conditions
   !> @return
   subroutine sll_s_poisson_3d_periodic_init(self, nx, ny, nz, Lx, Ly, Lz)

      sll_int32                              :: nx   !< number of points in x
      sll_int32                              :: ny   !< number of points in y
      sll_int32                              :: nz   !< number of points in z
      sll_int32                              :: ierr
      sll_real64                             :: Lx   !< Length along x
      sll_real64                             :: Ly   !< Length along y
      sll_real64                             :: Lz   !< Length along z
      type(sll_t_poisson_3d_periodic)       :: self !< Poisson solver type

      SLL_ALLOCATE(self%hat_rho(1:nx, 1:ny, 1:nz), ierr)
      SLL_ALLOCATE(self%hat_phi(1:nx, 1:ny, 1:nz), ierr)

      ! Geometry informations
      self%nx = nx
      self%ny = ny
      self%nz = nz
      self%Lx = Lx
      self%Ly = Ly
      self%Lz = Lz

      call sll_s_fft_init_c2c_3d(self%fw, nx, ny, nz, self%hat_rho, self%hat_phi, sll_p_fft_forward)
      call sll_s_fft_init_c2c_3d(self%bw, nx, ny, nz, self%hat_phi, self%hat_rho, sll_p_fft_backward)

   end subroutine sll_s_poisson_3d_periodic_init

   !> Compute the potential from 3d Poisson solver
   subroutine sll_s_poisson_3d_periodic_solve(self, rho, phi)

      type(sll_t_poisson_3d_periodic) :: self !< Solver structure
      sll_real64, dimension(:, :, :)     :: rho  !< charge density
      sll_real64, dimension(:, :, :)     :: phi  !< Electric potential
      sll_int32                        :: nx, ny, nz
      sll_int32                        :: i, j, k
      sll_real64                       :: Lx, Ly, Lz
      sll_real64                       :: ind_x, ind_y, ind_z

      nx = self%nx
      ny = self%ny
      nz = self%nz
      Lx = self%Lx
      Ly = self%Ly
      Lz = self%Lz

      SLL_ASSERT(nx == size(rho, 1))
      SLL_ASSERT(ny == size(rho, 2))
      SLL_ASSERT(nz == size(rho, 3))
      SLL_ASSERT(nx == size(phi, 1))
      SLL_ASSERT(ny == size(phi, 2))
      SLL_ASSERT(nz == size(phi, 3))

      self%hat_phi = cmplx(rho, 0_f64, kind=f64)
      call sll_s_fft_exec_c2c_3d(self%fw, self%hat_phi, self%hat_rho)
      self%hat_rho = self%hat_rho/cmplx(nx*ny*nz, 0.0, f64)

      ! Compute hat_phi, phi = inv_fft(hat_phi)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (i <= nx/2) then
                  ind_x = real(i - 1, f64)
               else
                  ind_x = real(nx - (i - 1), f64)
               end if
               if (j <= ny/2) then
                  ind_y = real(j - 1, f64)
               else
                  ind_y = real(ny - (j - 1), f64)
               end if
               if (k <= nz/2) then
                  ind_z = real(k - 1, f64)
               else
                  ind_z = real(nz - (k - 1), f64)
               end if
               if ((ind_x == 0) .and. (ind_y == 0) .and. (ind_z == 0)) then
                  self%hat_rho(i, j, k) = (0._f64, 0._f64)
               else
                  self%hat_rho(i, j, k) = self%hat_rho(i, j, k)/ &
                                          (4.0_f64*sll_p_pi**2* &
                                           ((real(ind_x, f64)/Lx)**2 &
                                            + (real(ind_y, f64)/Ly)**2 &
                                            + (real(ind_z, f64)/Lz)**2))
               end if
            end do
         end do
      end do

      call sll_s_fft_exec_c2c_3d(self%bw, self%hat_rho, self%hat_phi)

      phi = real(self%hat_phi, f64)

   end subroutine sll_s_poisson_3d_periodic_solve

   !> Delete the 3d poisson solver object
   subroutine sll_s_poisson_3d_periodic_free(self)

      type(sll_t_poisson_3d_periodic) :: self !< Poisson solver object
      sll_int32                       :: ierr

      call sll_s_fft_free(self%fw)
      call sll_s_fft_free(self%bw)

      SLL_DEALLOCATE_ARRAY(self%hat_rho, ierr)
      SLL_DEALLOCATE_ARRAY(self%hat_phi, ierr)

   end subroutine sll_s_poisson_3d_periodic_free

end module sll_m_poisson_3d_periodic
