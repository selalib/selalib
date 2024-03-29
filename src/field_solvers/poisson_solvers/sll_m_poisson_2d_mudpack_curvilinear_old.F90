#ifndef DOXYGEN_SHOULD_SKIP_THIS
!**************************************************************
!  Copyright INRIA
!  Authors :
!     CALVI project team
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
!> @author
!> Adnane Hamiaz (hamiaz@math.unistra.fr)
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!**************************************************************

!solves \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!       +\sum_{i=1}^2B_i\partial_i phi
!       +C \phi = rho
!in polar coordinates
!self leads when A_{1,2}=A_{2,1}=0 and B_2 = 0
! A_11\partial_{1,1}\hat{phi}+B_1\partial_{1}\hat{phi}+(C+A_{2,2}k^2)\hat{phi} = \hat{rho}

!> @ingroup poisson_solvers
module sll_m_poisson_2d_mudpack_curvilinear_old
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_mudpack, only: &
!   mud2, &
!   mud24, &
!   mud24cr, &
!   mud2cr

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet

   use sll_m_cartesian_meshes, only: &
      sll_t_cartesian_mesh_2d

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base

   use sll_m_cubic_spline_interpolator_2d, only: &
      sll_f_new_cubic_spline_interpolator_2d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   use sll_m_mudpack_curvilinear, only: &
      sll_p_non_separable_with_cross_terms, &
      sll_p_non_separable_without_cross_terms

   use sll_m_poisson_2d_base, only: &
      sll_c_poisson_2d_base, &
      sll_i_function_of_position

   implicit none

   public :: &
      sll_f_new_poisson_2d_mudpack_curvilinear_old

   private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_poisson_2d_base) :: poisson_2d_mudpack_curvilinear_old

      !type(sll_plan_poisson_polar), pointer                   :: poiss
      sll_real64, dimension(:, :), pointer :: cxx_2d
      sll_real64, dimension(:, :), pointer :: cxy_2d
      sll_real64, dimension(:, :), pointer :: cyy_2d
      sll_real64, dimension(:, :), pointer :: cx_2d
      sll_real64, dimension(:, :), pointer :: cy_2d
      sll_real64, dimension(:, :), pointer :: ce_2d
      sll_real64, dimension(:, :), pointer :: rho
      class(sll_c_interpolator_2d), pointer   :: cxx_2d_interp
      class(sll_c_interpolator_2d), pointer   :: cyy_2d_interp
      class(sll_c_interpolator_2d), pointer   :: cxy_2d_interp
      class(sll_c_interpolator_2d), pointer   :: cx_2d_interp
      class(sll_c_interpolator_2d), pointer   :: cy_2d_interp
      class(sll_c_interpolator_2d), pointer   :: ce_2d_interp
      class(sll_c_interpolator_2d), pointer   :: a11_interp
      class(sll_c_interpolator_2d), pointer   :: a22_interp
      class(sll_c_interpolator_2d), pointer   :: a12_interp
      class(sll_c_interpolator_2d), pointer   :: a21_interp
      class(sll_c_coordinate_transformation_2d_base), pointer :: transformation
      sll_int32  :: mudpack_curvilinear_case
      sll_real64, dimension(:), pointer :: work !< array for tmp data
      sll_int32  :: mgopt(4) !< Option to control multigrid
      sll_int32  :: iprm(16) !< Indices to control grid sizes
      sll_real64 :: fprm(6)  !< Real to set boundary conditions
      sll_int32  :: iguess   !< Initial solution or loop over time

   contains
      procedure, pass(poisson) :: initialize => &
         initialize_poisson_2d_mudpack_curvilinear_old
      procedure, pass(poisson) :: compute_phi_from_rho => &
         compute_phi_from_rho_2d_mudpack_curvilinear
      procedure, pass(poisson) :: compute_E_from_rho => &
         compute_E_from_rho_2d_mudpack_curvilinear
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar

      !> Compute the squarred L_2 for given coefficients
      procedure :: &
         l2norm_squared => l2norm_squarred_2d_mudpack_curvilinear
      !> Compute the right hand side from a given function
      procedure :: &
         compute_rhs_from_function => compute_rhs_from_function_2d_mudpack_curvilinear
      !> Destructor
      procedure :: free => delete_2d_mudpack_curvilinear_solver

   end type poisson_2d_mudpack_curvilinear_old

   class(poisson_2d_mudpack_curvilinear_old), pointer   :: mudpack_curvilinear_wrapper => null()

contains
   function sll_f_new_poisson_2d_mudpack_curvilinear_old( &
      transf, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      bc_eta1_left, &
      bc_eta1_right, &
      bc_eta2_left, &
      bc_eta2_right, &
      bc_interp2d_eta1, &
      bc_interp2d_eta2, &
      b11, &
      b12, &
      b21, &
      b22, &
      c, &
      mudpack_curvilinear_case) &
      result(poisson)

      type(poisson_2d_mudpack_curvilinear_old), pointer :: poisson
      sll_real64, intent(in) :: eta1_min
      sll_real64, intent(in) :: eta1_max
      sll_int32, intent(in) :: nc_eta1
      sll_real64, intent(in) :: eta2_min
      sll_real64, intent(in) :: eta2_max
      sll_int32, intent(in) :: nc_eta2
      sll_int32, intent(in) :: bc_eta1_left
      sll_int32, intent(in) :: bc_eta1_right
      sll_int32, intent(in) :: bc_eta2_left
      sll_int32, intent(in) :: bc_eta2_right
      sll_int32, intent(in) :: bc_interp2d_eta1
      sll_int32, intent(in) :: bc_interp2d_eta2
      sll_int32, intent(in), optional :: mudpack_curvilinear_case
      sll_real64, dimension(:, :), intent(in) :: b11
      sll_real64, dimension(:, :), intent(in) :: b12
      sll_real64, dimension(:, :), intent(in) :: b21
      sll_real64, dimension(:, :), intent(in) :: b22
      sll_real64, dimension(:, :), intent(in) :: c
      class(sll_c_coordinate_transformation_2d_base), pointer, intent(in) :: transf

      sll_int32 :: ierr

      SLL_ALLOCATE(poisson, ierr)

      call initialize_poisson_2d_mudpack_curvilinear_old( &
         poisson, &
         transf, &
         eta1_min, &
         eta1_max, &
         nc_eta1, &
         eta2_min, &
         eta2_max, &
         nc_eta2, &
         bc_eta1_left, &
         bc_eta1_right, &
         bc_eta2_left, &
         bc_eta2_right, &
         bc_interp2d_eta1, &
         bc_interp2d_eta2, &
         b11, &
         b12, &
         b21, &
         b22, &
         c, &
         mudpack_curvilinear_case)

   end function sll_f_new_poisson_2d_mudpack_curvilinear_old

   subroutine initialize_poisson_2d_mudpack_curvilinear_old( &
      poisson, &
      transf, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      bc_eta1_left, &
      bc_eta1_right, &
      bc_eta2_left, &
      bc_eta2_right, &
      bc_interp2d_eta1, &
      bc_interp2d_eta2, &
      b11, &
      b12, &
      b21, &
      b22, &
      c, &
      mudpack_curvilinear_case)
      class(poisson_2d_mudpack_curvilinear_old), target :: poisson
      sll_real64, intent(in) :: eta1_min
      sll_real64, intent(in) :: eta1_max
      sll_int32, intent(in) :: nc_eta1
      sll_real64, intent(in) :: eta2_min
      sll_real64, intent(in) :: eta2_max
      sll_int32, intent(in) :: nc_eta2
      sll_int32, intent(in) :: bc_eta1_left
      sll_int32, intent(in) :: bc_eta1_right
      sll_int32, intent(in) :: bc_eta2_left
      sll_int32, intent(in) :: bc_eta2_right
      sll_int32, intent(in) :: bc_interp2d_eta1
      sll_int32, intent(in) :: bc_interp2d_eta2
      sll_int32, intent(in), optional :: mudpack_curvilinear_case
      sll_real64, dimension(:, :), intent(in) :: b11
      sll_real64, dimension(:, :), intent(in) :: b12
      sll_real64, dimension(:, :), intent(in) :: b21
      sll_real64, dimension(:, :), intent(in) :: b22
      sll_real64, dimension(:, :), intent(in) :: c
      class(sll_c_coordinate_transformation_2d_base), pointer :: transf
      sll_real64, dimension(:, :), allocatable :: a12_array
      sll_real64, dimension(:, :), allocatable :: a21_array
      sll_int32 :: ierr
    !!!! begin variables for mudpack_curvilinear
      sll_int32, parameter   :: iixp = 2, jjyq = 2
      sll_int32               :: icall, iiex, jjey, llwork
      sll_real64, pointer :: phi(:) !< electric potential
      sll_real64, pointer :: rhs(:) !< charge density
      !put integer and floating point argument names in contiguous
      !storeage for labelling in vectors iprm,fprm
      sll_int32  :: iprm(16)
      sll_real64 :: fprm(6)
      !sll_int32  :: i
      sll_int32 :: error
      sll_int32  :: intl, nxa, nxb, nyc, nyd, ixp, jyq, iex, jey, nx, ny
      sll_int32  :: iguess, maxcy, method, nwork, lwrkqd, itero
      common/itmud2sp/intl, nxa, nxb, nyc, nyd, ixp, jyq, iex, jey, nx, ny, &
         iguess, maxcy, method, nwork, lwrkqd, itero
      sll_real64 :: xa, xb, yc, yd, tolmax, relmax
      common/ftmud2sp/xa, xb, yc, yd, tolmax, relmax
      equivalence(intl, iprm)
      equivalence(xa, fprm)

    !!!! end variables for mudpack_curvilinear
      sll_real64 :: delta1, delta2

      nx = nc_eta1 + 1
      ny = nc_eta2 + 1

      allocate (phi(nx*ny))
      allocate (rhs(nx*ny))

      delta1 = (eta1_max - eta1_min)/real(nc_eta1, f64)
      delta2 = (eta2_max - eta2_min)/real(nc_eta2, f64)
      ! set minimum required work space
      llwork = (7*(nx + 2)*(ny + 2) + 44*nx*ny)/3

      allocate (poisson%work(llwork))
      icall = 0
      iiex = ceiling(log((nx - 1.)/iixp)/log(2.)) + 1
      jjey = ceiling(log((ny - 1.)/jjyq)/log(2.)) + 1

      !set input integer arguments
      intl = 0

      !set boundary condition flags
      nxa = bc_eta1_left
      nxb = bc_eta1_right
      nyc = bc_eta2_left
      nyd = bc_eta2_right

      !set grid sizes from parameter statements
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey

      nx = ixp*(2**(iex - 1)) + 1
      ny = jyq*(2**(jey - 1)) + 1

      if (nx /= nc_eta1 + 1 .or. ny /= nc_eta2 + 1) then
         print *, "nx,nc_eta1+1=", nx, nc_eta1 + 1
         print *, "ny,nc_eta2+1=", ny, nc_eta2 + 1
         stop ' nx or ny different in sll_mudpack_curvilinear_cartesian '
      end if

      !set multigrid arguments (w(2,1) cycling with fully weighted
      !residual restriction and cubic prolongation)
      poisson%mgopt(1) = 2
      poisson%mgopt(2) = 2
      poisson%mgopt(3) = 1
      poisson%mgopt(4) = 3

      !set for three cycles to ensure second-order approximation is computed
      maxcy = 3

      !set no initial guess forcing full multigrid cycling
      poisson%iguess = 0
      iguess = poisson%iguess

      !set work space length approximation from parameter statement
      nwork = llwork

      !set point relaxation
      method = 0

      !set end points of solution rectangle in (x,y) space
      xa = eta1_min
      xb = eta1_max
      yc = eta2_min
      yd = eta2_max

      !set for no error control flag
      tolmax = 0.0_f64

!    write(*,101) (iprm(i),i=1,15)
!    write(*,102) (poisson%mgopt(i),i=1,4)
!    write(*,103) xa,xb,yc,yd,tolmax
!    write(*,104) intl

!call mud2sp(iprm,fprm,self%work,cofx,cofy,bndsp,rhs,phi,self%mgopt,error)
      if (present(mudpack_curvilinear_case)) then
         poisson%mudpack_curvilinear_case = mudpack_curvilinear_case
      else
         poisson%mudpack_curvilinear_case = sll_p_non_separable_with_cross_terms
      end if
      poisson%transformation => transf
      poisson%cxx_2d_interp => null()
      poisson%cyy_2d_interp => null()
      poisson%cx_2d_interp => null()
      poisson%cy_2d_interp => null()
      poisson%ce_2d_interp => null()
      poisson%a12_interp => null()
      poisson%a21_interp => null()
      SLL_ALLOCATE(poisson%rho(nc_eta1 + 1, nc_eta2 + 1), ierr)

      SLL_ALLOCATE(poisson%cxx_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)
      !SLL_ALLOCATE(poisson%cxy_2d(nc_eta1+1,nc_eta2+1),ierr)
      SLL_ALLOCATE(poisson%cyy_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)
      SLL_ALLOCATE(poisson%cx_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)
      SLL_ALLOCATE(poisson%cy_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)
      SLL_ALLOCATE(poisson%ce_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)

      poisson%a12_interp => sll_f_new_cubic_spline_interpolator_2d( &
                            nx, &
                            ny, &
                            eta1_min, &
                            eta1_max, &
                            eta2_min, &
                            eta2_max, &
                            bc_interp2d_eta1, &
                            bc_interp2d_eta2)
      poisson%a21_interp => sll_f_new_cubic_spline_interpolator_2d( &
                            nx, &
                            ny, &
                            eta1_min, &
                            eta1_max, &
                            eta2_min, &
                            eta2_max, &
                            bc_interp2d_eta1, &
                            bc_interp2d_eta2)
      SLL_ALLOCATE(a12_array(nx, ny), ierr)
      SLL_ALLOCATE(a21_array(nx, ny), ierr)
      call a12_a21_array(b11, b12, b21, b22, transf, eta1_min, &
                         eta2_min, delta1, delta2, nx, ny, a12_array, a21_array)
      call poisson%a12_interp%compute_interpolants(a12_array)
      call poisson%a21_interp%compute_interpolants(a21_array)

      poisson%cxx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                               nx, &
                               ny, &
                               eta1_min, &
                               eta1_max, &
                               eta2_min, &
                               eta2_max, &
                               bc_interp2d_eta1, &
                               bc_interp2d_eta2)
      poisson%cyy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                               nx, &
                               ny, &
                               eta1_min, &
                               eta1_max, &
                               eta2_min, &
                               eta2_max, &
                               bc_interp2d_eta1, &
                               bc_interp2d_eta2)
      call coefxxyy_array(b11, b12, b21, b22, transf, eta1_min, eta2_min, &
                          delta1, delta2, nx, ny, poisson%cxx_2d, poisson%cyy_2d)
      call poisson%cxx_2d_interp%compute_interpolants(poisson%cxx_2d)
      call poisson%cyy_2d_interp%compute_interpolants(poisson%cyy_2d)

      poisson%cx_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                              nx, &
                              ny, &
                              eta1_min, &
                              eta1_max, &
                              eta2_min, &
                              eta2_max, &
                              bc_interp2d_eta1, &
                              bc_interp2d_eta2)
      call coefx_array(eta1_min, eta2_min, delta1, delta2, nx, ny, &
                       poisson%cxx_2d_interp, poisson%a21_interp, poisson%cx_2d)
      call poisson%cx_2d_interp%compute_interpolants(poisson%cx_2d)

      poisson%cy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                              nx, &
                              ny, &
                              eta1_min, &
                              eta1_max, &
                              eta2_min, &
                              eta2_max, &
                              bc_interp2d_eta1, &
                              bc_interp2d_eta2)
      call coefy_array(eta1_min, eta2_min, delta1, delta2, nx, ny, &
                       poisson%cyy_2d_interp, poisson%a12_interp, poisson%cy_2d)
      call poisson%cy_2d_interp%compute_interpolants(poisson%cy_2d)

      poisson%ce_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                              nx, &
                              ny, &
                              eta1_min, &
                              eta1_max, &
                              eta2_min, &
                              eta2_max, &
                              bc_interp2d_eta1, &
                              bc_interp2d_eta2)
      poisson%ce_2d = -c
      call poisson%ce_2d_interp%compute_interpolants(poisson%ce_2d)

      !******sll_p_non_separable_with_cross_terms)
      select case (poisson%mudpack_curvilinear_case)

      case (sll_p_non_separable_without_cross_terms)
         if (associated(mudpack_curvilinear_wrapper)) then
            print *, '#Problem mudpack_curvilinear_wrapper is not null()'
            stop
         end if
         call associate_poisson(poisson)
         call mud2(iprm, fprm, poisson%work, &
                   mudpack_curvilinear_cof, &
                   mudpack_curvilinear_bndcr, &
                   rhs, &
                   phi, &
                   poisson%mgopt, &
                   error)
         mudpack_curvilinear_wrapper => null()

      case (sll_p_non_separable_with_cross_terms)
         SLL_ALLOCATE(poisson%cxy_2d(nc_eta1 + 1, nc_eta2 + 1), ierr)
         poisson%cxy_2d_interp => sll_f_new_cubic_spline_interpolator_2d( &
                                  nx, &
                                  ny, &
                                  eta1_min, &
                                  eta1_max, &
                                  eta2_min, &
                                  eta2_max, &
                                  bc_interp2d_eta1, &
                                  bc_interp2d_eta2)
         call coefxy_array(b11, b12, b21, b22, transf, eta1_min, eta2_min, &
                           delta1, delta2, nx, ny, poisson%cxy_2d)
         call poisson%cxy_2d_interp%compute_interpolants(poisson%cxy_2d)
         if (associated(mudpack_curvilinear_wrapper)) then
            print *, '#Problem mudpack_curvilinear_wrapper is not null()'
            stop
         end if
         call associate_poisson(poisson)
         call mud2cr(iprm, fprm, poisson%work, &
                     mudpack_curvilinear_cofcr, &
                     mudpack_curvilinear_bndcr, &
                     rhs, &
                     phi, &
                     poisson%mgopt, &
                     error)
         mudpack_curvilinear_wrapper => null()
      case default
         print *, '#bad mudpack_curvilinear_case', poisson%mudpack_curvilinear_case
         print *, '#in subroutine initialize_poisson_2d_mudpack_curvilinear_old'
         stop
      end select

   end subroutine initialize_poisson_2d_mudpack_curvilinear_old

   ! solves -\Delta phi = rho in 2d
   subroutine compute_phi_from_rho_2d_mudpack_curvilinear(poisson, phi, rho)
      class(poisson_2d_mudpack_curvilinear_old) :: poisson
      sll_real64, dimension(:, :), intent(in)     :: rho
      sll_real64, dimension(:, :), intent(out)    :: phi
      sll_int32 :: Nc_eta1
      sll_int32 :: Nc_eta2
      sll_real64 :: eta1_min
      sll_real64 :: eta2_min
      sll_real64 :: delta_eta1
      sll_real64 :: delta_eta2
      sll_real64 :: eta1
      sll_real64 :: eta2
      sll_int32 :: i1
      sll_int32 :: i2
      !sll_real64        :: phi(:,:)  !< Electric potential
      !sll_real64        :: rhs(:,:)  !< Charge density
      !put integer and floating point argument names in contiguous
      !storeage for labelling in vectors iprm,fprm
      sll_int32  :: iprm(16)
      sll_real64 :: fprm(6)
      sll_int32  :: error
      sll_int32  :: intl, nxa, nxb, nyc, nyd, ixp, jyq, iex, jey, nx, ny
      sll_int32  :: iguess, maxcy, method, nwork, lwrkqd, itero

      common/itmud2sp/intl, nxa, nxb, nyc, nyd, ixp, jyq, iex, jey, nx, ny, &
         iguess, maxcy, method, nwork, lwrkqd, itero
      sll_real64 :: xa, xb, yc, yd, tolmax, relmax
      common/ftmud2sp/xa, xb, yc, yd, tolmax, relmax

      equivalence(intl, iprm)
      equivalence(xa, fprm)

      !set initial guess because solve should be called every time step in a
      !time dependent problem and the elliptic operator does not depend on time.
      iguess = poisson%iguess

      !attempt solution
      intl = 1
      !write(*,106) intl,method,iguess

      associate (mesh => poisson%transformation%mesh)

         Nc_eta1 = mesh%num_cells1
         Nc_eta2 = mesh%num_cells2
         eta1_min = mesh%eta1_min
         eta2_min = mesh%eta2_min
         delta_eta1 = mesh%delta_eta1
         delta_eta2 = mesh%delta_eta2

      end associate

      poisson%rho(1:Nc_eta1 + 1, 1:Nc_eta2 + 1) = 0._f64
      do i2 = 1, Nc_eta2 + 1
         eta2 = eta2_min + real(i2 - 1, f64)*delta_eta2
         do i1 = 1, Nc_eta1 + 1
            eta1 = eta1_min + real(i1 - 1, f64)*delta_eta1
            poisson%rho(i1, i2) = -rho(i1, i2)*poisson%transformation%jacobian(eta1, eta2)
         end do
      end do

      if (nxa == sll_p_dirichlet) then
         do i2 = 1, Nc_eta2 + 1
            phi(1, i2) = 0._f64
         end do
      end if
      if (nxb == sll_p_dirichlet) then
         do i2 = 1, Nc_eta2 + 1
            phi(Nc_eta1 + 1, i2) = 0._f64
         end do
      end if
      if (nyc == sll_p_dirichlet) then
         do i1 = 1, Nc_eta1 + 1
            phi(i1, 1) = 0._f64
         end do
      end if
      if (nyd == sll_p_dirichlet) then
         do i1 = 1, Nc_eta1 + 1
            phi(i1, Nc_eta2 + 1) = 0._f64
         end do
      end if

      select case (poisson%mudpack_curvilinear_case)

      case (sll_p_non_separable_without_cross_terms)
         if (associated(mudpack_curvilinear_wrapper)) then
            print *, '#Problem mudpack_curvilinear_wrapper is not null()'
            stop
         end if
         call associate_poisson(poisson)
         call mud2(iprm, &
                   fprm, &
                   poisson%work, &
                   mudpack_curvilinear_cof, &
                   mudpack_curvilinear_bndcr, &
                   poisson%rho, &
                   phi, &
                   poisson%mgopt, &
                   error)
         !write(*,107) error
         if (error > 0) stop 0
         ! attempt to improve approximation to fourth order
         ! seems not to work for the moment
         call mud24(poisson%work, phi, error)
         !write (*,108) error
         if (error > 0) stop 0
         mudpack_curvilinear_wrapper => null()

      case (sll_p_non_separable_with_cross_terms)
         if (associated(mudpack_curvilinear_wrapper)) then
            print *, '#Problem mudpack_curvilinear_wrapper is not null()'
            stop
         end if
         call associate_poisson(poisson)
         call mud2cr(iprm, &
                     fprm, &
                     poisson%work, &
                     mudpack_curvilinear_cofcr, &
                     mudpack_curvilinear_bndcr, &
                     poisson%rho, &
                     phi, &
                     poisson%mgopt, &
                     error)
         !write(*,107) error
         if (error > 0) stop 0
         ! attempt to improve approximation to fourth order
         ! seems not to work for the moment
         call mud24cr(poisson%work, &
                      mudpack_curvilinear_cofcr, &
                      mudpack_curvilinear_bndcr, &
                      phi, &
                      error)
         !write (*,108) error
         if (error > 0) stop 0
         mudpack_curvilinear_wrapper => null()

      case default
         print *, '#bad mudpack_curvilinear_case', poisson%mudpack_curvilinear_case
         print *, '#in subroutine initialize_poisson_2d_mudpack_curvilinear_old'
         stop
      end select

   end subroutine compute_phi_from_rho_2d_mudpack_curvilinear

   ! solves E = -\nabla Phi with -\Delta phi = rho in 2d
   subroutine compute_E_from_rho_2d_mudpack_curvilinear(poisson, E1, E2, rho)
      class(poisson_2d_mudpack_curvilinear_old) :: poisson
      sll_real64, dimension(:, :), intent(in) :: rho
      sll_real64, dimension(:, :), intent(out) :: E1
      sll_real64, dimension(:, :), intent(out) :: E2

      print *, '#compute_E_from_rho_2d_mudpack_curvilinear'
      print *, '#not implemented for the moment'
      E1 = 0._f64
      E2 = 0._f64
      print *, maxval(rho)

      if (.not. (associated(poisson%cxx_2d))) then
         print *, '#poisson%cxx_2d is not associated'
      end if

      stop

      !call solve( poisson%poiss, E1, E2, rho)

   end subroutine compute_E_from_rho_2d_mudpack_curvilinear

   function l2norm_squarred_2d_mudpack_curvilinear(poisson, coefs_dofs) result(r)
      class(poisson_2d_mudpack_curvilinear_old), intent(in)  :: poisson !< Poisson solver object.
      sll_real64, intent(in)                                  :: coefs_dofs(:, :) !< Values of the coefficient vectors for each DoF
      sll_real64                                     :: r

      print *, 'l2norm_squared not implemented for poisson_2d_mudpack_curvilinear_solver.'
      r = 0.0_f64

   end function l2norm_squarred_2d_mudpack_curvilinear

   subroutine compute_rhs_from_function_2d_mudpack_curvilinear(poisson, func, coefs_dofs)
      class(poisson_2d_mudpack_curvilinear_old)  :: poisson !< Poisson solver object.
      procedure(sll_i_function_of_position)          :: func !< Function to be projected.
      sll_real64, intent(out)                        :: coefs_dofs(:) !< Coefficients of the projection.

      print *, 'compute_rhs_from_function not implemented for poisson_2d_mudpack_curvilinear_solver.'

   end subroutine compute_rhs_from_function_2d_mudpack_curvilinear

   subroutine delete_2d_mudpack_curvilinear_solver(poisson)
      class(poisson_2d_mudpack_curvilinear_old)  :: poisson !< Poisson solver object.

   end subroutine delete_2d_mudpack_curvilinear_solver

   subroutine coefxxyy_array(b11, b12, b21, b22, transf, eta1_min, eta2_min, &
                             delta1, delta2, nx, ny, cxx_array, cyy_array)
      implicit none
      sll_real64                :: eta1, eta1_min, eta2_min
      sll_real64                :: eta2, delta1, delta2
      sll_int32                 :: i, j, nx, ny
      sll_real64, dimension(:, :):: cxx_array, cyy_array
      sll_real64, dimension(1:2, 1:2) :: jac_m
      class(sll_c_coordinate_transformation_2d_base), pointer :: transf
      sll_real64, dimension(:, :) :: b11
      sll_real64, dimension(:, :) :: b12
      sll_real64, dimension(:, :) :: b21
      sll_real64, dimension(:, :) :: b22

      do j = 1, ny
         eta2 = eta2_min + real(j - 1, f64)*delta2
         do i = 1, nx
            eta1 = eta1_min + real(i - 1, f64)*delta1
            jac_m = transf%jacobian_matrix(eta1, eta2)
            cxx_array(i, j) = (b11(i, j)*(jac_m(1, 2)*jac_m(1, 2) + jac_m(2, 2)*jac_m(2, 2)) - &
                            & b12(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))) &
                            /transf%jacobian(eta1, eta2)
            cyy_array(i, j) = (b22(i, j)*(jac_m(2, 1)*jac_m(2, 1) + jac_m(1, 1)*jac_m(1, 1)) - &
                            & b21(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))) &
                            /transf%jacobian(eta1, eta2)
         end do
      end do
   end subroutine coefxxyy_array

   subroutine coefxy_array(b11, b12, b21, b22, transf, eta1_min, eta2_min, &
                           delta1, delta2, nx, ny, cxy_array)
      implicit none
      sll_real64                :: eta1, eta1_min, eta2_min
      sll_real64                :: eta2, delta1, delta2
      sll_real64                :: a12, a21
      sll_int32                 :: i, j, nx, ny
      sll_real64, dimension(:, :):: cxy_array
      sll_real64, dimension(1:2, 1:2) :: jac_m
      class(sll_c_coordinate_transformation_2d_base), pointer :: transf
      sll_real64, dimension(:, :) :: b11
      sll_real64, dimension(:, :) :: b12
      sll_real64, dimension(:, :) :: b21
      sll_real64, dimension(:, :) :: b22

      do j = 1, ny
         eta2 = eta2_min + real(j - 1, f64)*delta2
         do i = 1, nx
            eta1 = eta1_min + real(i - 1, f64)*delta1
            jac_m = transf%jacobian_matrix(eta1, eta2)
            a12 = b12(i, j)*(jac_m(2, 1)*jac_m(2, 1) + jac_m(1, 1)*jac_m(1, 1)) - &
                            & b11(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))

            a21 = b21(i, j)*(jac_m(1, 2)*jac_m(1, 2) + jac_m(2, 2)*jac_m(2, 2)) - &
                            & b22(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))
            cxy_array(i, j) = (a12 + a21)/transf%jacobian(eta1, eta2)
            !write(100,*) eta1,eta2, cxy_array(i,j), transf%jacobian(eta1,eta2)
         end do
      end do
   end subroutine coefxy_array

   subroutine a12_a21_array(b11, b12, b21, b22, transf, eta1_min, eta2_min, delta1, delta2, nx, ny, a12_array, a21_array)
      implicit none
      sll_real64                :: eta1, eta1_min, eta2_min
      sll_real64                :: eta2, delta1, delta2
      sll_real64                :: a12, a21
      sll_int32                 :: i, j, nx, ny
      sll_real64, dimension(:, :):: a12_array
      sll_real64, dimension(:, :):: a21_array
      class(sll_c_coordinate_transformation_2d_base), pointer :: transf
      sll_real64, dimension(:, :) :: b11
      sll_real64, dimension(:, :) :: b12
      sll_real64, dimension(:, :) :: b21
      sll_real64, dimension(:, :) :: b22
      sll_real64, dimension(1:2, 1:2) :: jac_m

      do j = 1, ny
         eta2 = eta2_min + real(j - 1, f64)*delta2
         do i = 1, nx
            eta1 = eta1_min + real(i - 1, f64)*delta1
            jac_m = transf%jacobian_matrix(eta1, eta2)
            a12 = b12(i, j)*(jac_m(2, 1)*jac_m(2, 1) + jac_m(1, 1)*jac_m(1, 1)) - &
                            & b11(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))
            a12_array(i, j) = a12/transf%jacobian(eta1, eta2)
            a21 = b21(i, j)*(jac_m(1, 2)*jac_m(1, 2) + jac_m(2, 2)*jac_m(2, 2)) - &
                            & b22(i, j)*(jac_m(2, 1)*jac_m(2, 2) + jac_m(1, 1)*jac_m(1, 2))
            a21_array(i, j) = a21/transf%jacobian(eta1, eta2)
         end do
      end do
   end subroutine a12_a21_array

   subroutine coefx_array(eta1_min, eta2_min, &
                          delta1, delta2, nx, ny, cxx_2d_interp, a21_interp, cx_array)
      implicit none
      sll_real64                :: eta1, eta1_min, eta2_min
      sll_real64                :: eta2, delta1, delta2
      sll_int32                 :: i, j, nx, ny
      sll_real64, dimension(:, :):: cx_array
      class(sll_c_interpolator_2d), pointer   :: cxx_2d_interp
      class(sll_c_interpolator_2d), pointer   :: a21_interp

      do j = 1, ny
         eta2 = eta2_min + real(j - 1, f64)*delta2
         do i = 1, nx
            eta1 = eta1_min + real(i - 1, f64)*delta1
            cx_array(i, j) = cxx_2d_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2) + &
                             a21_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2)
         end do
      end do
   end subroutine coefx_array

   subroutine coefy_array(eta1_min, eta2_min, &
                          delta1, delta2, nx, ny, cyy_2d_interp, a12_interp, cy_array)
      implicit none
      sll_real64                :: eta1, eta1_min, eta2_min
      sll_real64                :: eta2, delta1, delta2
      sll_int32                 :: i, j, nx, ny
      sll_real64, dimension(:, :):: cy_array
      class(sll_c_interpolator_2d), pointer   :: cyy_2d_interp
      class(sll_c_interpolator_2d), pointer   :: a12_interp

      do j = 1, ny
         eta2 = eta2_min + real(j - 1, f64)*delta2
         do i = 1, nx
            eta1 = eta1_min + real(i - 1, f64)*delta1
            cy_array(i, j) = cyy_2d_interp%interpolate_from_interpolant_derivative_eta2(eta1, eta2) + &
                             a12_interp%interpolate_from_interpolant_derivative_eta1(eta1, eta2)
         end do
      end do
   end subroutine coefy_array

   subroutine mudpack_curvilinear_cof(x, y, cxx, cyy, cx, cy, ce)
      real(8)  :: x, cxx, cx
      real(8)  :: y, cyy, cy, ce
      cxx = mudpack_curvilinear_wrapper%cxx_2d_interp%interpolate_from_interpolant_value(x, y)
      cyy = mudpack_curvilinear_wrapper%cyy_2d_interp%interpolate_from_interpolant_value(x, y)
      cx = mudpack_curvilinear_wrapper%cx_2d_interp%interpolate_from_interpolant_value(x, y)
      cy = mudpack_curvilinear_wrapper%cy_2d_interp%interpolate_from_interpolant_value(x, y)
      ce = mudpack_curvilinear_wrapper%ce_2d_interp%interpolate_from_interpolant_value(x, y)
      return
   end subroutine

   subroutine mudpack_curvilinear_cofcr(x, y, cxx, cxy, cyy, cx, cy, ce)
      real(8)  :: x, cxx, cx, cxy
      real(8)  :: y, cyy, cy, ce
      cxx = mudpack_curvilinear_wrapper%cxx_2d_interp%interpolate_from_interpolant_value(x, y)
      cxy = mudpack_curvilinear_wrapper%cxy_2d_interp%interpolate_from_interpolant_value(x, y)
      cyy = mudpack_curvilinear_wrapper%cyy_2d_interp%interpolate_from_interpolant_value(x, y)
      cx = mudpack_curvilinear_wrapper%cx_2d_interp%interpolate_from_interpolant_value(x, y)
      cy = mudpack_curvilinear_wrapper%cy_2d_interp%interpolate_from_interpolant_value(x, y)
      ce = mudpack_curvilinear_wrapper%ce_2d_interp%interpolate_from_interpolant_value(x, y)

      return
   end subroutine

!> input mixed derivative b.c. to mud2sp
   subroutine mudpack_curvilinear_bndcr(kbdy, xory, alfa, gbdy)
      integer  :: kbdy
      real(8)  :: xory, alfa, gbdy, x, y, pe, px, py
      real(8)  :: xa, xb, yc, yd, tolmax, relmax
      common/ftmud2sp/xa, xb, yc, yd, tolmax, relmax

!subroutine not used in periodic case
      if (kbdy == 1) then  ! x=xa boundary
         y = xory
         x = xa
         alfa = -1._f64
         gbdy = px + alfa*pe
         return
      end if

      if (kbdy == 4) then  ! y=yd boundary
         y = yd
         x = xory
         alfa = 1._f64
         gbdy = py + alfa*pe
         return
      end if
   end subroutine

   subroutine associate_poisson(poisson)
      type(poisson_2d_mudpack_curvilinear_old), target :: poisson

      mudpack_curvilinear_wrapper => poisson

   end subroutine associate_poisson

end module sll_m_poisson_2d_mudpack_curvilinear_old

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

