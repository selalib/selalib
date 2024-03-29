module sll_m_fcisl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_1d, &
      sll_t_cartesian_mesh_1d

   implicit none

   public :: &
      sll_s_compute_at_aligned, &
      sll_s_compute_derivative_periodic, &
      sll_s_compute_iota_from_shift, &
      sll_s_compute_oblic_shift, &
      sll_s_compute_spaghetti, &
      sll_s_compute_spaghetti_and_shift_from_guess, &
      sll_s_compute_spaghetti_size_from_shift, &
      sll_s_compute_w_hermite, &
      sll_s_load_spaghetti, &
      sll_f_new_oblic_derivative, &
      sll_t_oblic_derivative, &
      sll_s_unload_spaghetti

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! a direction is chosen for interpolation
! derivative computation and advection

   type :: sll_t_oblic_derivative
      sll_int32 :: degree
      sll_real64, dimension(:, :), pointer :: buf
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x1
      type(sll_t_cartesian_mesh_1d), pointer :: mesh_x2
      class(sll_c_advector_1d), pointer :: adv
      sll_real64, dimension(:), pointer :: w
   end type sll_t_oblic_derivative

contains

   subroutine sll_s_compute_oblic_shift(iota, Nc_x1, shift, iota_modif)
      sll_real64, intent(in) :: iota
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(out) :: shift
      sll_real64, intent(out) :: iota_modif

      shift = floor(iota*Nc_x1 + 0.5)
      iota_modif = real(shift, f64)/real(Nc_x1, f64)

   end subroutine sll_s_compute_oblic_shift

   subroutine sll_s_compute_iota_from_shift(Nc_x1, shift, iota_modif)
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: shift
      sll_real64, intent(out) :: iota_modif

      iota_modif = real(shift, f64)/real(Nc_x1, f64)

   end subroutine sll_s_compute_iota_from_shift

   !warning: normalization to correct/choose
   subroutine sll_s_compute_at_aligned( &
      f_input, &
      f_output, &
      Nc_x1, &
      Nc_x2, &
      adv, &
      x1_min, &
      x1_max, &
      iota)
      sll_real64, dimension(:, :), intent(in) :: f_input
      sll_real64, dimension(:, :), intent(out) :: f_output
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: Nc_x2
      class(sll_c_advector_1d), pointer :: adv
      sll_real64, intent(in) :: x1_min
      sll_real64, intent(in) :: x1_max
      sll_real64, intent(in) :: iota
      !local variables
      sll_int32 :: i
      sll_real64 :: A
      sll_real64 :: delta_x1
      A = iota*real(Nc_x1, f64)/real(Nc_x2, f64)
      delta_x1 = (x1_max - x1_min)/real(Nc_x1, f64)

      do i = 1, Nc_x2 + 1
         call adv%advect_1d_constant( &
            A, &
            real(i - 1, f64)*delta_x1, &
            f_input(1:Nc_x1 + 1, i), &
            f_output(1:Nc_x1 + 1, i))
      end do
   end subroutine sll_s_compute_at_aligned

   subroutine sll_s_compute_spaghetti_size_from_shift( &
      Nc_x1, &
      shift, &
      spaghetti_size)
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: shift
      sll_int32, intent(out) :: spaghetti_size
      !local variables
      sll_int32 :: i
      sll_int32 :: i_loc

      spaghetti_size = 1
      do i = 1, Nc_x1
         i_loc = modulo(i*shift, Nc_x1)
         if (i_loc .ne. 0) then
            spaghetti_size = spaghetti_size + 1
         else
            exit
         end if
      end do
   end subroutine sll_s_compute_spaghetti_size_from_shift

!< compute shift and spaghetti_size
!< from iota_guess and spaghetti_size_guess
!< first search the existing spaghetti_size nearest of spaghetti_size_guess
!< for a shift between -Nc_x1 to Nc_x1
!< then looks for a shift that is near shift guess that has the same spaghetti_size
   subroutine sll_s_compute_spaghetti_and_shift_from_guess( &
      Nc_x1, &
      Nc_x2, &
      iota_guess, &
      spaghetti_size_guess, &
      shift, &
      spaghetti_size)
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: Nc_x2
      !sll_int32, intent(in) :: shift_guess
      sll_real64, intent(in) :: iota_guess
      sll_int32, intent(in) :: spaghetti_size_guess
      sll_int32, intent(out) :: shift
      sll_int32, intent(out) :: spaghetti_size
      !local variables
      sll_int32 :: shift_guess
      sll_int32 :: shift_plus
      sll_int32 :: shift_minus
      sll_int32 :: i
      sll_int32 :: s
      sll_int32 :: i_val
      sll_int32 :: val
      sll_int32 :: spaghetti_size_new

      SLL_ASSERT(Nc_x1 > 0)
      SLL_ASSERT(Nc_x2 > 0)

      s = 0
      val = 2*Nc_x1
      do i = -Nc_x1, Nc_x1
         call sll_s_compute_spaghetti_size_from_shift( &
            Nc_x1, &
            i, &
            spaghetti_size)
         if (abs(spaghetti_size_guess - spaghetti_size) < abs(val)) then
            val = spaghetti_size_guess - spaghetti_size
            i_val = i
            if (val == 0) then
               exit
            end if
         end if
      end do

      call sll_s_compute_spaghetti_size_from_shift( &
         Nc_x1, &
         i_val, &
         spaghetti_size)

      shift_guess = floor(iota_guess*Nc_x1)

      shift_plus = shift_guess
      shift_minus = shift_guess

      do i = 0, Nc_x1
         call sll_s_compute_spaghetti_size_from_shift( &
            Nc_x1, &
            shift_guess + i, &
            spaghetti_size_new)
         if (spaghetti_size_new .eq. spaghetti_size) then
            shift_plus = shift_guess + i
            exit
         end if
      end do
      do i = 0, Nc_x1
         call sll_s_compute_spaghetti_size_from_shift( &
            Nc_x1, &
            shift_guess - i, &
            spaghetti_size_new)
         if (spaghetti_size_new .eq. spaghetti_size) then
            shift_minus = shift_guess - i
            exit
         end if
      end do

      if (abs(shift_minus*Nc_x1 - iota_guess) < abs(shift_plus*Nc_x1 - iota_guess)) then
         shift = shift_minus
      else
         shift = shift_plus
      end if

      call sll_s_compute_spaghetti_size_from_shift( &
         Nc_x1, &
         shift, &
         spaghetti_size)

!    call sll_s_compute_spaghetti_size_from_shift( &
!      Nc_x1, &
!      shift_guess, &
!      spaghetti_size)

   end subroutine sll_s_compute_spaghetti_and_shift_from_guess

   subroutine sll_s_compute_spaghetti( &
      Nc_x1, &
      Nc_x2, &
      shift, &
      spaghetti_index, &
      spaghetti_size)
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: Nc_x2
      sll_int32, intent(in) :: shift
      sll_int32, dimension(:), intent(out) :: spaghetti_index
      sll_int32, intent(out) :: spaghetti_size
      !local variables
      sll_int32 :: i
      sll_int32 :: i_loc
      sll_int32 :: ii
      sll_int32 :: j
      sll_int32 :: q
      sll_int32, dimension(:), allocatable ::  check
      sll_int32 :: ierr

      SLL_ASSERT(Nc_x1 > 0)
      SLL_ASSERT(Nc_x2 > 0)
      SLL_ALLOCATE(check(Nc_x1 + 1), ierr)
      !0, shift,2*shift, Nc_x1*shift
      !while(modulo(k*shift,Nc_x1) .ne. 0)

      spaghetti_index = 0
      spaghetti_index(1) = 1
      spaghetti_size = 1
      do i = 1, Nc_x1
         i_loc = modulo(i*shift, Nc_x1)
         if (i_loc .ne. 0) then
            spaghetti_size = spaghetti_size + 1
            spaghetti_index(i + 1) = i_loc + 1
         else
            exit
         end if
      end do
      !print *,'#spaghetti_size=',shift,spaghetti_size
      !print *,'#i,i_loc=',i,i_loc
      q = Nc_x1/spaghetti_size
      do j = 1, q - 1
         do ii = 1, spaghetti_size
            spaghetti_index(ii + j*spaghetti_size) = modulo(spaghetti_index(ii) + j - 1, Nc_x1) + 1
         end do
      end do
      spaghetti_index(Nc_x1 + 1) = spaghetti_index(1)
      check = 0
      do i = 1, Nc_x1
         if (spaghetti_index(i) < 1) then
            print *, '#problem with spaghetti_index'
            stop
         end if
         if (spaghetti_index(i) > Nc_x1 + 1) then
            print *, '#problem with spaghetti_index'
            stop
         end if
         check(spaghetti_index(i)) = check(spaghetti_index(i)) + 1
         !print *,i,spaghetti_index(i)
      end do
      if (maxval(check(1:Nc_x1)) .ne. 1) then
         print *, '#problem checking spaghetti_index'
         print *, '#maxval=', maxval(check(1:Nc_x1))
         stop
      end if
      if (minval(check(1:Nc_x1)) .ne. 1) then
         print *, '#problem checking spaghetti_index'
         print *, '#minval=', minval(check(1:Nc_x1))
         stop
      end if

   end subroutine sll_s_compute_spaghetti

   subroutine sll_s_load_spaghetti( &
      input, &
      output, &
      spaghetti_index, &
      !    spaghetti_size, &
      Npts_x1, &
      Npts_x2)
      sll_real64, dimension(:, :), intent(in) :: input
      sll_real64, dimension(:), intent(out) :: output
      sll_int32, dimension(:), intent(in) :: spaghetti_index
!    sll_int32, intent(in) :: spaghetti_size
      sll_int32, intent(in) :: Npts_x1
      sll_int32, intent(in) :: Npts_x2
      !local variables
      sll_int32 :: i
      sll_int32 :: j
      sll_int32 :: s
      s = 0
      do i = 1, Npts_x1
         do j = 1, Npts_x2
            s = s + 1
            output(s) = input(spaghetti_index(i), j)
         end do
      end do
   end subroutine sll_s_load_spaghetti

   subroutine sll_s_unload_spaghetti( &
      input, &
      output, &
      spaghetti_index, &
      !    spaghetti_size, &
      Npts_x1, &
      Npts_x2)
      sll_real64, dimension(:), intent(in) :: input
      sll_real64, dimension(:, :), intent(out) :: output
      sll_int32, dimension(:), intent(in) :: spaghetti_index
!    sll_int32, intent(in) :: spaghetti_size
      sll_int32, intent(in) :: Npts_x1
      sll_int32, intent(in) :: Npts_x2
      !local variables
      sll_int32 :: i
      sll_int32 :: j
      sll_int32 :: s
      s = 0
      do i = 1, Npts_x1
         do j = 1, Npts_x2
            s = s + 1
            output(spaghetti_index(i), j) = input(s)
         end do
      end do
   end subroutine sll_s_unload_spaghetti

   function sll_f_new_oblic_derivative( &
      degree, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      Nc_x1, &
      Nc_x2, &
      adv &
      ) result(res)
      type(sll_t_oblic_derivative), pointer :: res
      sll_int32, intent(in) :: degree
      sll_real64, intent(in) :: x1_min
      sll_real64, intent(in) :: x1_max
      sll_real64, intent(in) :: x2_min
      sll_real64, intent(in) :: x2_max
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: Nc_x2
      class(sll_c_advector_1d), pointer :: adv
      !local variables
      sll_int32 :: ierr
      SLL_ALLOCATE(res, ierr)
      call initialize_oblic_derivative( &
         res, &
         degree, &
         x1_min, &
         x1_max, &
         x2_min, &
         x2_max, &
         Nc_x1, &
         Nc_x2, &
         adv)

   end function sll_f_new_oblic_derivative

   subroutine initialize_oblic_derivative( &
      deriv, &
      degree, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      Nc_x1, &
      Nc_x2, &
      adv)
      type(sll_t_oblic_derivative) :: deriv
      sll_int32, intent(in) :: degree
      sll_real64, intent(in) :: x1_min
      sll_real64, intent(in) :: x1_max
      sll_real64, intent(in) :: x2_min
      sll_real64, intent(in) :: x2_max
      sll_int32, intent(in) :: Nc_x1
      sll_int32, intent(in) :: Nc_x2
      class(sll_c_advector_1d), pointer :: adv
      sll_int32 :: ierr

      deriv%degree = degree
      deriv%mesh_x1 => sll_f_new_cartesian_mesh_1d(Nc_x1, eta_min=x1_min, eta_max=x1_max)
      deriv%mesh_x2 => sll_f_new_cartesian_mesh_1d(Nc_x2, eta_min=x2_min, eta_max=x2_max)
      deriv%adv => adv

      SLL_ALLOCATE(deriv%buf(1:Nc_x2 + 2*degree + 1, 1:Nc_x1 + 1), ierr)
      SLL_ALLOCATE(deriv%w(-degree:Nc_x2 + degree), ierr)

      call compute_finite_difference_init(deriv%w, 2*degree)

   end subroutine initialize_oblic_derivative

   subroutine sll_s_compute_w_hermite(w, r, s)
      sll_int32, intent(in)::r, s
      sll_real64, dimension(r:s), intent(out)::w
      sll_int32 ::i, j
      sll_real64::tmp

      !maple code for generation of w
      !for k from r to -1 do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
      !od:
      !for k from 1 to s do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
      !od:
      !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):

      do i = r, -1
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do

      do i = 1, s
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do
      tmp = 0._f64
      do i = r, -1
         tmp = tmp + w(i)
      end do
      do i = 1, s
         tmp = tmp + w(i)
      end do
      w(0) = -tmp

      !print *,'w',w
      !do ii=r,s
      !  print *,ii,w(r+s-ii)
      !enddo
      !

   end subroutine sll_s_compute_w_hermite

   subroutine sll_s_compute_derivative_periodic( &
      input, &
      output, &
      Nc, &
      w, &
      r, &
      s, &
      length)
      sll_real64, dimension(:), intent(in) :: input
      sll_real64, dimension(:), intent(out) :: output
      sll_int32, intent(in) :: Nc
      sll_int32, intent(in) :: r
      sll_int32, intent(in) :: s
      sll_real64, dimension(r:s), intent(in) :: w
      sll_real64, intent(in) :: length
      !local variables
      sll_int32 :: i
      sll_int32 :: ii
      sll_real64 :: tmp
      sll_int32 :: ind
      sll_real64 :: dx
      dx = length/real(Nc, f64)
      do i = 1, Nc + 1
         tmp = 0._f64
         do ii = r, s
            ind = modulo(i + ii - 1 + Nc, Nc) + 1
            tmp = tmp + w(ii)*input(ind)
         end do
         output(i) = tmp/dx
      end do

   end subroutine sll_s_compute_derivative_periodic

   subroutine compute_finite_difference_init(w, p)
      integer, intent(in)::p
      real(f64), dimension(-p/2:(p + 1)/2), intent(out)::w
      integer::r, s, i, j
      real(f64)::tmp

      r = -p/2
      s = (p + 1)/2

!    if(modulo(p,2)==0)then
!      r=-p/2
!      s=p/2
!    else
!      r=-(p-1)/2
!      s=(p+1)/2
!    endif

      !maple code for generation of w
      !for k from r to -1 do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
      !od:
      !for k from 1 to s do
      !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
      !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
      !od:
      !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):

      do i = r, -1
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do

      do i = 1, s
         tmp = 1._f64
         do j = r, i - 1
            tmp = tmp*real(i - j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(i - j, f64)
         end do
         tmp = 1._f64/tmp
         do j = r, -1
            tmp = tmp*real(-j, f64)
         end do
         do j = 1, i - 1
            tmp = tmp*real(-j, f64)
         end do
         do j = i + 1, s
            tmp = tmp*real(-j, f64)
         end do
         w(i) = tmp
      end do
      tmp = 0._f64
      do i = r, -1
         tmp = tmp + w(i)
      end do
      do i = 1, s
         tmp = tmp + w(i)
      end do
      w(0) = -tmp

   end subroutine compute_finite_difference_init

!  subroutine compute_field_from_phi_cartesian_1d(phi,mesh,A,interp)
!    sll_real64, dimension(:), intent(in) :: phi
!    sll_real64, dimension(:), intent(out) :: A
!    type(sll_t_cartesian_mesh_1d), pointer :: mesh
!    class(sll_c_interpolator_1d), pointer   :: interp
!    sll_int32 :: Nc_x1
!    sll_real64 :: x1_min
!    sll_real64 :: delta_x1
!    sll_real64 :: x1
!    sll_int32 :: i1
!
!    Nc_x1 = mesh%num_cells
!    x1_min = mesh%eta_min
!    delta_x1 = mesh%delta_eta
!
!    call interp%compute_interpolants(phi)
!
!    do i1=1,Nc_x1+1
!      x1=x1_min+real(i1-1,f64)*delta_x1
!      A(i1)=interp%interpolate_from_interpolant_derivative_eta1(x1)
!    end do
!  end subroutine compute_field_from_phi_cartesian_1d

end module sll_m_fcisl
