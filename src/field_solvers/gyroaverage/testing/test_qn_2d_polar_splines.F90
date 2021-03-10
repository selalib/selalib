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

program test_qn_2d_polar_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: sll_p_pi
   use sll_m_gyroaverage_utilities, only: sll_s_compute_init_f_polar
   use sll_m_qn_2d_base, only: sll_c_qn_2d_base
   use sll_m_qn_2d_polar, only: sll_s_compute_gamma0, &
                                sll_s_qn_2d_polar_test_solve

   use sll_m_qn_2d_polar_splines_solver, only: &
      sll_t_qn_2d_polar_splines_solver

   implicit none

   class(sll_c_qn_2d_base), pointer :: qn
   type(sll_t_qn_2d_polar_splines_solver), target  :: qn_splines

   sll_real64                              :: mu_max
   sll_real64                              :: eta_min(2)
   sll_real64                              :: eta_max(2)
   sll_int32                               :: Nc(2)
   sll_int32                               :: mode(2)
   sll_int32                               :: N_points
   sll_int32                               :: N_mu
   sll_int32                               :: mu_case
   sll_real64, dimension(:), allocatable :: mu_points
   sll_real64, dimension(:), allocatable :: mu_weights
   sll_real64, dimension(:), allocatable :: lambda
   sll_real64, dimension(:), allocatable :: T_i
   sll_real64                              :: val
   sll_real64, dimension(:, :), allocatable :: phi
   sll_real64, dimension(:, :), allocatable :: phi_restr
   sll_real64, dimension(:, :), allocatable :: phi_init
   sll_real64, dimension(:, :), allocatable :: phi_qn
   sll_int32                               :: ierr
   sll_int32                               :: i

   eta_min(1) = 2._f64 ! 0.1_f64
   eta_max(1) = 18._f64 ! 0.9_f64
   eta_min(2) = 0._f64
   eta_max(2) = 2._f64*sll_p_pi

   mode(1) = 1
   mode(2) = 1
   Nc(1) = 32
   Nc(2) = 32 ! Ntheta doit etre pair
   mu_case = 1  ! 1 : Rectangles ! 2 : Lin-Lee ! 3 : Gauss-Laguerre
   N_mu = 20
   mu_max = 3.0_f64
   N_points = 32

   val = 1._f64
   call sll_s_compute_gamma0(mode, eta_min, eta_max, val)
   print *, "gamma0 = ", val

   SLL_ALLOCATE(phi(1:Nc(1) + 1, 1:Nc(2) + 1), ierr)
   SLL_ALLOCATE(phi_init(1:Nc(1) + 1, 1:Nc(2) + 1), ierr)
   SLL_ALLOCATE(phi_qn(1:Nc(1) + 1, 1:Nc(2) + 1), ierr)
   SLL_ALLOCATE(phi_restr(1:Nc(1) + 1, 1:Nc(2)), ierr)
   SLL_ALLOCATE(mu_points(0:N_mu - 1), ierr)
   SLL_ALLOCATE(mu_weights(0:N_mu - 1), ierr)
   SLL_ALLOCATE(lambda(1:Nc(1) + 1), ierr)
   SLL_ALLOCATE(T_i(1:Nc(1) + 1), ierr)

! mu integration
   if (mu_case == 1) then

      if (N_mu == 1) then
         mu_points(0) = 1._f64
         mu_weights(0) = 1._f64
      else
         do i = 1, N_mu
            mu_points(i - 1) = real(i, f64)*mu_max/real(N_mu, f64)
            mu_weights(i - 1) = mu_max/real(N_mu, f64)
         end do
      end if

   elseif (mu_case == 2) then

      if (N_mu == 2) then
         mu_points(0) = 0.4167845_f64
         mu_weights(0) = 0.7194_f64
         mu_points(1) = 2.495154605_f64
         mu_weights(1) = 0.2806_f64
      elseif (N_mu == 3) then
         mu_points(0) = 0.148131245_f64
         mu_weights(0) = 0.3583_f64
         mu_points(1) = 1._f64
         mu_weights(1) = 0.5004_f64
         mu_points(2) = 3.15959522_f64
         mu_weights(2) = 0.1413_f64
      end if

   elseif (mu_case == 3) then

      mu_points = [real(8) :: 111.751398098_8, &
                   0.0444893658333_8, &
                   0.23452610952_8, &
                   0.576884629302_8, &
                   1.07244875382_8, &
                   1.72240877644_8, &
                   2.52833670643_8, &
                   3.49221327302_8, &
                   4.61645676975_8, &
                   5.90395850417_8, &
                   7.35812673319_8, &
                   8.98294092421_8, &
                   10.7830186325_8, &
                   12.7636979867_8, &
                   14.9311397555_8, &
                   17.2924543367_8, &
                   19.8558609403_8, &
                   22.6308890132_8, &
                   25.6286360225_8, &
                   28.8621018163_8, &
                   32.346629154_8, &
                   36.1004948058_8, &
                   40.1457197715_8, &
                   44.5092079958_8, &
                   49.2243949873_8, &
                   54.3337213334_8, &
                   59.8925091621_8, &
                   65.9753772879_8, &
                   72.6876280907_8, &
                   80.1874469779_8, &
                   98.8295428683_8, &
                   88.7353404179_8]

      mu_weights = [real(8) :: 8.22732055833d-33, &
                    0.109218341952_8, &
                    0.210443107939_8, &
                    0.23521322967_8, &
                    0.195903335973_8, &
                    0.129983786286_8, &
                    0.0705786238657_8, &
                    0.0317609125092_8, &
                    0.0119182148348_8, &
                    0.00373881629461_8, &
                    0.000980803306615_8, &
                    0.000214864918801_8, &
                    3.92034196799d-05, &
                    5.93454161287d-06, &
                    7.41640457867d-07, &
                    7.60456787912d-08, &
                    6.35060222663d-09, &
                    4.28138297104d-10, &
                    2.30589949189d-11, &
                    9.79937928873d-13, &
                    3.23780165773d-14, &
                    8.17182344342d-16, &
                    1.54213383339d-17, &
                    2.11979229016d-19, &
                    2.05442967379d-21, &
                    1.34698258664d-23, &
                    5.6612941304d-26, &
                    1.41856054546d-28, &
                    1.91337549445d-31, &
                    1.1922487601d-34, &
                    1.33861694201d-42, &
                    2.67151121924d-38]

   end if

   print *, "mu_points =", mu_points
   print *, "mu_weights =", mu_weights

   call sll_s_compute_init_f_polar(phi, mode, Nc, eta_min, eta_max)
   phi_init = phi
   phi_restr(1:Nc(1) + 1, 1:Nc(2)) = phi(1:Nc(1) + 1, 1:Nc(2))

   do i = 1, Nc(1) + 1
      lambda(i) = 1.0_f64 + real(i, 8)*0.1_f64
      T_i(i) = 1.7_f64 + real(i, 8)*0.2_f64
   end do

   call qn_splines%init( &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      lambda, &
      T_i)

   qn => qn_splines

   call qn%precompute_qn(mu_points(0:N_mu - 1), mu_weights(0:N_mu - 1), N_mu)
   call qn%solve_qn(phi_restr)

   phi_qn(1:Nc(1) + 1, 1:Nc(2)) = phi_restr(1:Nc(1) + 1, 1:Nc(2))
   phi_qn(1:Nc(1) + 1, Nc(2) + 1) = phi_restr(1:Nc(1) + 1, 1)

   call sll_s_qn_2d_polar_test_solve( &
      Nc, &
      eta_min, &
      eta_max, &
      mu_points(0:N_mu - 1), &
      mu_weights(0:N_mu - 1), &
      N_mu, mode, &
      lambda, &
      T_i, &
      phi_init, &
      phi_qn)

   print *, '#PASSED'

end program test_qn_2d_polar_splines
