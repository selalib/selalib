module sll_m_reduction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_f_compute_integral_trapezoid_1d, &
      sll_s_compute_reduction_2d_to_0d, &
      sll_s_compute_reduction_4d_to_2d_direction34, &
      sll_s_compute_reduction_4d_to_3d_direction4, &
      sll_s_compute_reduction_diag_4d_to_2d_direction12

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   abstract interface
      function sll_integration_discrete_1d(data, Npts, delta, params) result(res)
         use sll_m_working_precision
         sll_real64                                  :: res
         sll_real64, dimension(:), intent(in)         :: data
         sll_int32, intent(in)                       :: Npts
         sll_real64, intent(in)                       :: delta
         sll_real64, dimension(:), intent(in), optional :: params
      end function sll_integration_discrete_1d
   end interface

contains

   !--------------------------------------------------
   ! Generic function which can be used for computing charge density
   ! by default we suppose a uniform mesh node based
   ! and trapezoid quadrature formula is used
   ! as an option we can use a way to integrate
   ! by providing a function whose signature is sll_integration_discrete_1d
   !---------------------------------------------------

   subroutine sll_s_compute_reduction_4d_to_3d_direction4( &
      data_4d, &
      data_3d, &
      Npts1, &
      Npts2, &
      Npts3, &
      Npts4, &
      delta4, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :, :, :), intent(in)    :: data_4d
      sll_real64, dimension(:, :, :), intent(out) :: data_3d
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_int32, intent(in)  :: Npts3
      sll_int32, intent(in)  :: Npts4
      sll_real64, intent(in) :: delta4
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2, i3, i4
      sll_real64 :: tmp

      if (Npts1 > size(data_4d, 1)) then
         print *, '#Problem for size1 in sll_s_compute_reduction_4d_to_3d_direction4'
         print *, '#Npts=', (/Npts1, Npts2, Npts3, Npts4/)
         print *, '#size(data_4d)=', (/ &
            size(data_4d, 1), &
            size(data_4d, 2), &
            size(data_4d, 3), &
            size(data_4d, 4)/)
         stop
      end if
      if (Npts2 > size(data_4d, 2)) then
         print *, '#Problem for size2 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (Npts3 > size(data_4d, 3)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (Npts4 > size(data_4d, 4)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (.not. (present(integration_func))) then
         do i3 = 1, Npts3
            do i2 = 1, Npts2
               do i1 = 1, Npts1
                  tmp = 0.5_f64*(data_4d(i1, i2, i3, 1) &
                                 + data_4d(i1, i2, i3, Npts4))
                  do i4 = 2, Npts4 - 1
                     tmp = tmp + data_4d(i1, i2, i3, i4)
                  end do
                  data_3d(i1, i2, i3) = tmp*delta4
                  !data_3d(i1,i2,i3) = delta4*sum(data_4d(i1,i2,i3,1:Npts4-1))
               end do
            end do
         end do
      else
         do i3 = 1, Npts3
            do i2 = 1, Npts2
               do i1 = 1, Npts1
                  data_3d(i1, i2, i3) = integration_func( &
                                        data_4d(i1, i2, i3, :), &
                                        Npts4, &
                                        delta4, &
                                        integration_func_params)
               end do
            end do
         end do
      end if
   end subroutine sll_s_compute_reduction_4d_to_3d_direction4

   subroutine compute_reduction_4d_to_3d_direction4_accumulate( &
      data_4d, &
      data_3d, &
      Npts1, &
      Npts2, &
      Npts3, &
      Npts4, &
      delta4, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :, :, :), intent(in)    :: data_4d
      sll_real64, dimension(:, :, :), intent(inout) :: data_3d
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_int32, intent(in)  :: Npts3
      sll_int32, intent(in)  :: Npts4
      sll_real64, intent(in) :: delta4
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2, i3, i4
      sll_real64 :: tmp

      if (Npts1 > size(data_4d, 1)) then
         print *, '#Problem for size1 in sll_s_compute_reduction_4d_to_3d_direction4'
         print *, '#Npts=', (/Npts1, Npts2, Npts3, Npts4/)
         print *, '#size(data_4d)=', (/ &
            size(data_4d, 1), &
            size(data_4d, 2), &
            size(data_4d, 3), &
            size(data_4d, 4)/)
         stop
      end if
      if (Npts2 > size(data_4d, 2)) then
         print *, '#Problem for size2 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (Npts3 > size(data_4d, 3)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (Npts4 > size(data_4d, 4)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_3d_direction4'
         stop
      end if
      if (.not. (present(integration_func))) then
         do i3 = 1, Npts3
            do i2 = 1, Npts2
               do i1 = 1, Npts1
                  tmp = 0.5_f64*(data_4d(i1, i2, i3, 1) &
                                 + data_4d(i1, i2, i3, Npts4))
                  do i4 = 2, Npts4 - 1
                     tmp = tmp + data_4d(i1, i2, i3, i4)
                  end do
                  data_3d(i1, i2, i3) = data_3d(i1, i2, i3) + tmp*delta4
                  !data_3d(i1,i2,i3) = delta4*sum(data_4d(i1,i2,i3,1:Npts4-1))
               end do
            end do
         end do
      else
         do i3 = 1, Npts3
            do i2 = 1, Npts2
               do i1 = 1, Npts1
                  data_3d(i1, i2, i3) = data_3d(i1, i2, i3) + integration_func( &
                                        data_4d(i1, i2, i3, :), &
                                        Npts4, &
                                        delta4, &
                                        integration_func_params)
               end do
            end do
         end do
      end if
   end subroutine compute_reduction_4d_to_3d_direction4_accumulate

   subroutine sll_s_compute_reduction_4d_to_2d_direction34( &
      data_4d, &
      data_2d, &
      Npts1, &
      Npts2, &
      Npts3, &
      Npts4, &
      delta3, &
      delta4, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :, :, :), intent(in)    :: data_4d
      sll_real64, dimension(:, :), intent(out) :: data_2d
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_int32, intent(in)  :: Npts3
      sll_int32, intent(in)  :: Npts4
      sll_real64, intent(in) :: delta3
      sll_real64, intent(in) :: delta4
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2, i3, i4
      sll_real64 :: tmp

      if (Npts1 > size(data_4d, 1)) then
         print *, '#Problem for size1 in sll_s_compute_reduction_4d_to_2d_direction34'
         print *, 'Npts1=', Npts1, 'size(data_4d,1)', size(data_4d, 1)
         stop
      end if
      if (Npts2 > size(data_4d, 2)) then
         print *, '#Problem for size2 in sll_s_compute_reduction_4d_to_2d_direction34'
         stop
      end if
      if (Npts3 > size(data_4d, 3)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_2d_direction34'
         stop
      end if
      if (Npts4 > size(data_4d, 4)) then
         print *, '#Problem for size3 in sll_s_compute_reduction_4d_to_2d_direction34'
         stop
      end if
      if (.not. (present(integration_func))) then
         do i2 = 1, Npts2
            do i1 = 1, Npts1
               data_2d(i1, i2) = 0._f64
               i3 = 1
               tmp = 0.5_f64*(data_4d(i1, i2, i3, 1) &
                              + data_4d(i1, i2, i3, Npts4))
               do i4 = 2, Npts4 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
               end do
               tmp = tmp*delta4
               data_2d(i1, i2) = data_2d(i1, i2) + 0.5_f64*tmp

               do i3 = 2, Npts3 - 1
                  tmp = 0.5_f64*(data_4d(i1, i2, i3, 1) &
                                 + data_4d(i1, i2, i3, Npts4))
                  do i4 = 2, Npts4 - 1
                     tmp = tmp + data_4d(i1, i2, i3, i4)
                  end do
                  tmp = tmp*delta4
                  data_2d(i1, i2) = data_2d(i1, i2) + tmp
               end do

               i3 = Npts3
               tmp = 0.5_f64*(data_4d(i1, i2, i3, 1) &
                              + data_4d(i1, i2, i3, Npts4))
               do i4 = 2, Npts4 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
               end do
               tmp = tmp*delta4
               data_2d(i1, i2) = data_2d(i1, i2) + 0.5_f64*tmp

               data_2d(i1, i2) = data_2d(i1, i2)*delta3
            end do
         end do
      else
         if (present(integration_func_params)) then
            print *, 'integration_func_params is present'
         end if
         print *, '#not implemented yet'
         print *, '#in sll_s_compute_reduction_4d_to_2d_direction34'
         stop
      end if
   end subroutine sll_s_compute_reduction_4d_to_2d_direction34

   subroutine compute_reduction_4d_to_2d_direction12( &
      data_4d, &
      data_2d, &
      Npts1, &
      Npts2, &
      Npts3, &
      Npts4, &
      delta1, &
      delta2, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :, :, :), intent(in)    :: data_4d
      sll_real64, dimension(:, :), intent(out) :: data_2d
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_int32, intent(in)  :: Npts3
      sll_int32, intent(in)  :: Npts4
      sll_real64, intent(in) :: delta1
      sll_real64, intent(in) :: delta2
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2, i3, i4
      sll_real64 :: tmp

      if (Npts1 > size(data_4d, 1)) then
         print *, '#Problem for size1 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts1=', Npts1, 'size(data_4d,1)', size(data_4d, 1)
         stop
      end if
      if (Npts2 > size(data_4d, 2)) then
         print *, '#Problem for size2 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts2=', Npts2, 'size(data_4d,2)', size(data_4d, 2)
         stop
      end if
      if (Npts3 > size(data_4d, 3)) then
         print *, '#Problem for size3 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts3=', Npts3, 'size(data_4d,3)', size(data_4d, 3)
         stop
      end if
      if (Npts4 > size(data_4d, 4)) then
         print *, '#Problem for size3 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts4=', Npts4, 'size(data_4d,4)', size(data_4d, 4)
         stop
      end if
      if (.not. (present(integration_func))) then
         do i4 = 1, Npts4
            do i3 = 1, Npts3
               data_2d(i3, i4) = 0._f64
               i1 = 1
               tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                              + data_4d(i1, Npts2, i3, i4))
               do i2 = 2, Npts2 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
               end do
               tmp = tmp*delta2
               data_2d(i3, i4) = data_2d(i3, i4) + 0.5_f64*tmp

               do i1 = 2, Npts1 - 1
                  tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                                 + data_4d(i1, Npts2, i3, i4))
                  do i2 = 2, Npts2 - 1
                     tmp = tmp + data_4d(i1, i2, i3, i4)
                  end do
                  tmp = tmp*delta2
                  data_2d(i3, i4) = data_2d(i3, i4) + tmp
               end do

               i1 = Npts1
               tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                              + data_4d(i1, Npts2, i3, i4))
               do i2 = 2, Npts2 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
               end do
               tmp = tmp*delta2
               data_2d(i3, i4) = data_2d(i3, i4) + 0.5_f64*tmp

               data_2d(i3, i4) = data_2d(i3, i4)*delta1
            end do
         end do
      else
         if (present(integration_func_params)) then
            print *, 'integration_func_params is present'
         end if
         print *, '#not implemented yet'
         print *, '#in compute_reduction_4d_to_2d_direction12'
         stop
      end if
   end subroutine compute_reduction_4d_to_2d_direction12

   subroutine sll_s_compute_reduction_diag_4d_to_2d_direction12( &
      data_4d, &
      data_diag_2d, &
      Npts1, &
      Npts2, &
      Npts3, &
      Npts4, &
      delta1, &
      delta2, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :, :, :), intent(in)    :: data_4d
      sll_real64, dimension(:, :, :), intent(out) :: data_diag_2d
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_int32, intent(in)  :: Npts3
      sll_int32, intent(in)  :: Npts4
      sll_real64, intent(in) :: delta1
      sll_real64, intent(in) :: delta2
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2, i3, i4
      sll_real64 :: tmp
      sll_real64 :: tmp_l1
      sll_real64 :: tmp_l2

      if (Npts1 > size(data_4d, 1)) then
         print *, '#Problem for size1 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts1=', Npts1, 'size(data_4d,1)', size(data_4d, 1)
         stop
      end if
      if (Npts2 > size(data_4d, 2)) then
         print *, '#Problem for size2 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts2=', Npts2, 'size(data_4d,2)', size(data_4d, 2)
         stop
      end if
      if (Npts3 > size(data_4d, 3)) then
         print *, '#Problem for size3 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts3=', Npts3, 'size(data_4d,3)', size(data_4d, 3)
         stop
      end if
      if (Npts4 > size(data_4d, 4)) then
         print *, '#Problem for size3 in compute_reduction_4d_to_2d_direction12'
         print *, 'Npts4=', Npts4, 'size(data_4d,4)', size(data_4d, 4)
         stop
      end if
      if (.not. (present(integration_func))) then
         do i4 = 1, Npts4
            do i3 = 1, Npts3
               data_diag_2d(i3, i4, 1:3) = 0._f64
               i1 = 1
               tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                              + data_4d(i1, Npts2, i3, i4))
               tmp_l1 = 0.5_f64*(abs(data_4d(i1, 1, i3, i4)) &
                                 + abs(data_4d(i1, Npts2, i3, i4)))
               tmp_l2 = 0.5_f64*((data_4d(i1, 1, i3, i4))**2 &
                                 + (data_4d(i1, Npts2, i3, i4))**2)
               do i2 = 2, Npts2 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
                  tmp_l1 = tmp_l1 + abs(data_4d(i1, i2, i3, i4))
                  tmp_l2 = tmp_l2 + data_4d(i1, i2, i3, i4)**2
               end do
               tmp = tmp*delta2
               tmp_l1 = tmp_l1*delta2
               tmp_l2 = tmp_l2*delta2
               data_diag_2d(i3, i4, 1) = data_diag_2d(i3, i4, 1) + 0.5_f64*tmp
               data_diag_2d(i3, i4, 2) = data_diag_2d(i3, i4, 2) + 0.5_f64*tmp_l1
               data_diag_2d(i3, i4, 3) = data_diag_2d(i3, i4, 3) + 0.5_f64*tmp_l2

               do i1 = 2, Npts1 - 1
                  tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                                 + data_4d(i1, Npts2, i3, i4))
                  tmp_l1 = 0.5_f64*(abs(data_4d(i1, 1, i3, i4)) &
                                    + abs(data_4d(i1, Npts2, i3, i4)))
                  tmp_l2 = 0.5_f64*((data_4d(i1, 1, i3, i4))**2 &
                                    + (data_4d(i1, Npts2, i3, i4))**2)
                  do i2 = 2, Npts2 - 1
                     tmp = tmp + data_4d(i1, i2, i3, i4)
                     tmp_l1 = tmp_l1 + abs(data_4d(i1, i2, i3, i4))
                     tmp_l2 = tmp_l2 + data_4d(i1, i2, i3, i4)**2
                  end do
                  tmp = tmp*delta2
                  tmp_l1 = tmp_l1*delta2
                  tmp_l2 = tmp_l2*delta2
                  data_diag_2d(i3, i4, 1) = data_diag_2d(i3, i4, 1) + tmp
                  data_diag_2d(i3, i4, 2) = data_diag_2d(i3, i4, 2) + tmp_l1
                  data_diag_2d(i3, i4, 3) = data_diag_2d(i3, i4, 3) + tmp_l2
               end do

               i1 = Npts1
               tmp = 0.5_f64*(data_4d(i1, 1, i3, i4) &
                              + data_4d(i1, Npts2, i3, i4))
               tmp_l1 = 0.5_f64*(abs(data_4d(i1, 1, i3, i4)) &
                                 + abs(data_4d(i1, Npts2, i3, i4)))
               tmp_l2 = 0.5_f64*((data_4d(i1, 1, i3, i4))**2 &
                                 + (data_4d(i1, Npts2, i3, i4))**2)
               do i2 = 2, Npts2 - 1
                  tmp = tmp + data_4d(i1, i2, i3, i4)
                  tmp_l1 = tmp_l1 + abs(data_4d(i1, i2, i3, i4))
                  tmp_l2 = tmp_l2 + data_4d(i1, i2, i3, i4)**2
               end do
               tmp = tmp*delta2
               tmp_l1 = tmp_l1*delta2
               tmp_l2 = tmp_l2*delta2
               data_diag_2d(i3, i4, 1) = data_diag_2d(i3, i4, 1) + 0.5_f64*tmp
               data_diag_2d(i3, i4, 2) = data_diag_2d(i3, i4, 2) + 0.5_f64*tmp_l1
               data_diag_2d(i3, i4, 3) = data_diag_2d(i3, i4, 3) + 0.5_f64*tmp_l2

               data_diag_2d(i3, i4, 1:3) = data_diag_2d(i3, i4, 1:3)*delta1
            end do
         end do
      else
         if (present(integration_func_params)) then
            print *, 'integration_func_params is present'
         end if
         print *, '#not implemented yet'
         print *, '#in compute_reduction_4d_to_2d_direction12'
         stop
      end if
   end subroutine sll_s_compute_reduction_diag_4d_to_2d_direction12

   subroutine sll_s_compute_reduction_2d_to_0d( &
      data_2d, &
      res, &
      Npts1, &
      Npts2, &
      delta1, &
      delta2, &
      integration_func, &
      integration_func_params)

      sll_real64, dimension(:, :), intent(in)    :: data_2d
      sll_real64, intent(out) :: res
      sll_int32, intent(in)  :: Npts1
      sll_int32, intent(in)  :: Npts2
      sll_real64, intent(in) :: delta1
      sll_real64, intent(in) :: delta2
      procedure(sll_integration_discrete_1d), optional :: integration_func
      sll_real64, dimension(:), optional :: integration_func_params
      sll_int32  :: i1, i2
      sll_real64 :: tmp

      if (Npts1 > size(data_2d, 1)) then
         print *, '#Problem for size1 in sll_s_compute_reduction_2d_to_0d'
         print *, 'Npts1=', Npts1, 'size(data_2d,1)', size(data_2d, 1)
         stop
      end if
      if (Npts2 > size(data_2d, 2)) then
         print *, '#Problem for size2 in sll_s_compute_reduction_2d_to_0d'
         print *, 'Npts2=', Npts2, 'size(data_2d,2)', size(data_2d, 2)
         stop
      end if
      if (.not. (present(integration_func))) then
         res = 0._f64
         i1 = 1
         tmp = 0.5_f64*(data_2d(i1, 1) &
                        + data_2d(i1, Npts2))
         do i2 = 2, Npts2 - 1
            tmp = tmp + data_2d(i1, i2)
         end do
         tmp = tmp*delta2
         res = res + 0.5_f64*tmp

         do i1 = 2, Npts1 - 1
            tmp = 0.5_f64*(data_2d(i1, 1) &
                           + data_2d(i1, Npts2))
            do i2 = 2, Npts2 - 1
               tmp = tmp + data_2d(i1, i2)
            end do
            tmp = tmp*delta2
            res = res + tmp
         end do

         i1 = Npts1
         tmp = 0.5_f64*(data_2d(i1, 1) &
                        + data_2d(i1, Npts2))
         do i2 = 2, Npts2 - 1
            tmp = tmp + data_2d(i1, i2)
         end do
         tmp = tmp*delta2
         res = res + 0.5_f64*tmp

         res = res*delta1
      else
         if (present(integration_func_params)) then
            print *, 'integration_func_params is present'
         end if
         print *, '#not implemented yet'
         print *, '#in sll_s_compute_reduction_2d_to_0d'
         stop
      end if
   end subroutine sll_s_compute_reduction_2d_to_0d

   function sll_f_compute_integral_trapezoid_1d(data, Npts, delta, func_params) result(res)
      sll_real64, dimension(:), intent(in)    :: data
      sll_int32, intent(in) :: Npts
      sll_real64, intent(in) :: delta
      sll_real64, dimension(:), intent(in), optional :: func_params
      sll_real64 :: res
      sll_int32 :: i

      if (present(func_params)) then
         if (size(func_params) > 100) then
            print *, '#size of func_params is >100'
         end if
      end if

      res = 0.5*(data(1) + data(Npts))
      do i = 2, Npts - 1
         res = res + data(i)
      end do
      res = res*delta
   end function sll_f_compute_integral_trapezoid_1d

   function compute_integral_conservative_1d(data, Npts, node_positions) result(res)
      sll_real64, dimension(:), intent(in)    :: data
      sll_int32, intent(in) :: Npts
      sll_real64, dimension(:), intent(in) :: node_positions
      sll_real64 :: res
      sll_int32 :: i

      res = 0._f64
      do i = 1, Npts - 1
         res = res + data(i)*(node_positions(i + 1) - node_positions(i))
      end do
   end function compute_integral_conservative_1d

end module sll_m_reduction
