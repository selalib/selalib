!> @ingroup splines
!> @brief
!> Provides capabilities for data interpolation with cubic B-splines on
!> non-uniform meshes

module sll_m_cubic_non_uniform_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic

   use sll_m_tridiagonal, only: &
      sll_s_setup_cyclic_tridiag, &
      sll_o_solve_cyclic_tridiag

   implicit none

   public :: &
      sll_s_compute_spline_nonunif, &
      sll_s_compute_spline_nonunif_1d_periodic_aux2, &
      sll_t_cubic_nonunif_spline_1d, &
      sll_s_interpolate_array_value_nonunif, &
      sll_s_interpolate_array_value_nonunif_aux, &
      sll_f_new_cubic_nonunif_spline_1d, &
      sll_s_setup_spline_nonunif_1d_periodic_aux, &
      sll_o_delete

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Spline object
   type :: sll_t_cubic_nonunif_spline_1d
      sll_int32                         :: n_cells        !< number of cells
      sll_real64, dimension(:), pointer :: node_positions !< the non uniform
      !< mesh normalized on (0,1)
      sll_real64, dimension(:), pointer :: buf            !< memory buffer
      sll_int32                         :: size_buf       !< size of the buffer
      sll_int32, dimension(:), pointer ::  ibuf           !< memory buffer of integers
      sll_int32                         :: size_ibuf      !< integers buffer size
      sll_real64, dimension(:), pointer :: coeffs         !< spline coefficients
      sll_real64                        :: xmin           !< left corner
      sll_real64                        :: xmax           !< right corner
      sll_int32                         :: bc_type        !< periodic or hermite
      LOGICAL                           :: is_setup       !< check splines setup
      sll_real64                        :: slope_L        !< left slope, for Hermite
      sll_real64                        :: slope_R        !< right slope, for Hermite
   contains

      procedure :: init => sll_s_cubic_nonunif_spline_1d_init

   end type sll_t_cubic_nonunif_spline_1d

   interface sll_o_delete
      module procedure delete_cubic_nonunif_spline_1d
   end interface sll_o_delete

contains  ! ****************************************************************

   !> create new spline object
   function sll_f_new_cubic_nonunif_spline_1d(n_cells, bc_type)

      type(sll_t_cubic_nonunif_spline_1d), pointer  :: sll_f_new_cubic_nonunif_spline_1d
      sll_int32, intent(in)                        :: n_cells
      sll_int32, intent(in)                        :: bc_type
      sll_int32                                     :: ierr

      SLL_ALLOCATE(sll_f_new_cubic_nonunif_spline_1d, ierr)
      call sll_f_new_cubic_nonunif_spline_1d%init(n_cells, bc_type)

   end function sll_f_new_cubic_nonunif_spline_1d

   !> Initialize spline object
   subroutine sll_s_cubic_nonunif_spline_1d_init(self, n_cells, bc_type)

      class(sll_t_cubic_nonunif_spline_1d) :: self
      sll_int32, intent(in)               :: n_cells
      sll_int32, intent(in)               :: bc_type
      sll_int32                            :: ierr, size_buf, size_ibuf

      self%n_cells = n_cells
      self%bc_type = bc_type

      SLL_ALLOCATE(self%node_positions(-2:n_cells + 2), ierr)
      size_buf = 10*(n_cells + 1)
      size_ibuf = n_cells + 1
      self%size_buf = size_buf
      SLL_ALLOCATE(self%buf(size_buf), ierr)
      self%size_ibuf = size_ibuf
      SLL_ALLOCATE(self%ibuf(size_ibuf), ierr)
      SLL_ALLOCATE(self%coeffs(-1:n_cells + 1), ierr)
      self%is_setup = .FALSE.
      self%slope_L = 0.0_f64
      self%slope_R = 0.0_f64

   end subroutine sll_s_cubic_nonunif_spline_1d_init

   !> delete spline object
   subroutine delete_cubic_nonunif_spline_1D(spline, ierr)
      type(sll_t_cubic_nonunif_spline_1d), pointer :: spline
      sll_int32                    :: ierr
      SLL_ASSERT(associated(spline))
      SLL_DEALLOCATE(spline%node_positions, ierr)
      SLL_DEALLOCATE(spline%buf, ierr)
      SLL_DEALLOCATE(spline%coeffs, ierr)
      SLL_DEALLOCATE(spline%ibuf, ierr)
      SLL_DEALLOCATE(spline, ierr)
   end subroutine delete_cubic_nonunif_spline_1D

   !> compute splines coefficients
   subroutine sll_s_compute_spline_nonunif(f, spline, node_positions, sl, sr)
      sll_real64, dimension(:), intent(in) :: f    ! data to be fit
      type(sll_t_cubic_nonunif_spline_1d), pointer         :: spline
      sll_real64, dimension(:), intent(in), optional :: node_positions    ! non uniform mesh
      sll_real64, intent(in), optional     :: sl
      sll_real64, intent(in), optional     :: sr
      sll_int32  ::n_cells, bc_type
      sll_real64, dimension(:), pointer :: node_pos
      sll_real64 :: length

      !check that the spline is initialized
      SLL_ASSERT(associated(spline))
      n_cells = spline%n_cells
      bc_type = spline%bc_type
      !copy on non_uniform mesh in spline object

      if (present(node_positions)) then
         node_pos => spline%node_positions
         spline%xmin = node_positions(1)
         spline%xmax = node_positions(n_cells + 1)
         length = spline%xmax - spline%xmin
         node_pos(1:n_cells - 1) = (node_positions(2:n_cells) - node_positions(1))/length
         node_pos(0) = 0.0_f64
         node_pos(n_cells) = 1.0_f64
         select case (bc_type)
         case (sll_p_periodic)
            call sll_s_setup_spline_nonunif_1d_periodic_aux(node_pos, n_cells, spline%buf, spline%ibuf)
         case (sll_p_hermite)
            call setup_spline_nonunif_1D_hermite_aux(node_pos, n_cells, spline%buf, spline%ibuf)
         case default
            print *, 'ERROR: compute_spline_nonunif_1D(): not recognized boundary condition'
            STOP
         end select
         spline%is_setup = .TRUE.
      end if

      if (.not. (spline%is_setup)) then
         ! FIXME: THROW ERROR
         print *, 'ERROR: compute_spline_nonunif_1D_periodic(): '
         print *, 'node_positions are to be given first'
         STOP
      end if

      select case (bc_type)
      case (sll_p_periodic)
         call compute_spline_nonunif_1D_periodic(f, spline)
      case (sll_p_hermite)
         if (present(sl)) then
            spline%slope_L = sl
         end if
         if (present(sr)) then
            spline%slope_R = sr
         end if
         call compute_spline_nonunif_1D_hermite(f, spline)
      case default
         print *, 'ERROR: compute_spline_nonunif_1D(): not recognized boundary condition'
         STOP
      end select
   end subroutine sll_s_compute_spline_nonunif

   !>
   subroutine compute_spline_nonunif_1D_periodic(f, spline)
      sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
      type(sll_t_cubic_nonunif_spline_1d), pointer         :: spline
      sll_real64, dimension(:), pointer :: buf, fp, coeffs
      sll_int32, dimension(:), pointer :: ibuf
      sll_int32 :: nc
      if (.not. associated(spline)) then
         ! FIXME: THROW ERROR
         print *, 'ERROR: compute_spline_nonunif_1D_periodic(): ', &
            'uninitialized spline object passed as argument. Exiting... '
         STOP
      end if
      if (.not. (size(f) .ge. spline%n_cells + 1)) then
         ! FIXME: THROW ERROR
         print *, 'ERROR: compute_spline_nonunif_1D_periodic(): '
         write (*, '(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_cells + 1, ' . Passed size: ', size(f)
         STOP
      end if
      nc = spline%n_cells

      fp => f
      buf => spline%buf
      ibuf => spline%ibuf
      coeffs => spline%coeffs
      call compute_spline_nonunif_1D_periodic_aux(fp, nc, buf, ibuf, coeffs)
   end subroutine compute_spline_nonunif_1D_periodic

   !>
   subroutine compute_spline_nonunif_1D_hermite(f, spline)
      sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
      type(sll_t_cubic_nonunif_spline_1d), pointer         :: spline
      sll_real64, dimension(:), pointer :: buf, fp, coeffs
      sll_int32, dimension(:), pointer :: ibuf
      sll_int32 :: nc
      sll_real64 :: lift(4, 2)
      if (.not. associated(spline)) then
         ! FIXME: THROW ERROR
         print *, 'ERROR: compute_spline_nonunif_1D_hermite(): ', &
            'uninitialized spline object passed as argument. Exiting... '
         STOP
      end if
      if (.not. (size(f) .ge. spline%n_cells + 1)) then
         ! FIXME: THROW ERROR
         print *, 'ERROR: compute_spline_nonunif_1D_hermite(): '
         write (*, '(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_cells + 1, ' . Passed size: ', size(f)
         STOP
      end if
      nc = spline%n_cells

      fp => f
      buf => spline%buf
      ibuf => spline%ibuf
      coeffs => spline%coeffs

      !lift(1,1) = 1._f64/3._f64*(spline%node_positions(1)-spline%node_positions(0))*spline%slope_L
      !lift(1,2) = -1._f64/3._f64*(spline%node_positions(nc)-spline%node_positions(nc-1))*spline%slope_R
      !lift(2,1) = -2._f64/3._f64*(spline%node_positions(1)+spline%node_positions(2)-2.0_f64*spline%node_positions(0))*spline%slope_L
      !lift(2,2) = 2._f64/3._f64*(2.0_f64*spline%node_positions(nc)-spline%node_positions(nc-1)-spline%node_positions(nc-2))*spline%slope_R

      lift(1, 1) = 1._f64/3._f64*(spline%node_positions(1) - spline%node_positions(0))*spline%slope_L
      lift(1, 2) = -1._f64/3._f64*(spline%node_positions(nc) - spline%node_positions(nc - 1))*spline%slope_R
      lift(2, 1) = -1._f64/3._f64*(spline%node_positions(1) - spline%node_positions(-1))&
       &*(spline%node_positions(1) - spline%node_positions(-2))/(spline%node_positions(1) - spline%node_positions(0))*spline%slope_L
      lift(2, 2) = 1._f64/3._f64*(spline%node_positions(nc + 1) - spline%node_positions(nc - 1))&
  &*(spline%node_positions(nc+2)-spline%node_positions(nc-1))/(spline%node_positions(nc)-spline%node_positions(nc-1))*spline%slope_R

      lift(1:2, 1:2) = lift(1:2, 1:2)*(spline%xmax - spline%xmin)

      lift(3, 1) = (spline%node_positions(2) - spline%node_positions(0))*(spline%node_positions(1) - spline%node_positions(0))
     lift(3,1) = lift(3,1)-(spline%node_positions(0)-spline%node_positions(-1))*(spline%node_positions(0)-spline%node_positions(-2))
    lift(3,1) = lift(3,1)/((spline%node_positions(2)-spline%node_positions(-1))*(spline%node_positions(1)-spline%node_positions(0)))

      lift(3,2) = (spline%node_positions(nc+1)-spline%node_positions(nc))*(spline%node_positions(nc+2)-spline%node_positions(nc-1))
    lift(3,2) = lift(3,2)/((spline%node_positions(nc+1)-spline%node_positions(nc-2))*(spline%node_positions(nc)-spline%node_positions(nc-1)))

      lift(4, 1) = (spline%node_positions(0) - spline%node_positions(-1))*(spline%node_positions(1) - spline%node_positions(-2))
    lift(4,1) = lift(4,1)/((spline%node_positions(2)-spline%node_positions(-1))*(spline%node_positions(1)-spline%node_positions(0)))

lift(4, 2) = (spline%node_positions(nc) - spline%node_positions(nc - 2))*(spline%node_positions(nc) - spline%node_positions(nc - 1))
    lift(4,2) = lift(4,2)-(spline%node_positions(nc+2)-spline%node_positions(nc))*(spline%node_positions(nc+1)-spline%node_positions(nc))
    lift(4,2) = lift(4,2)/((spline%node_positions(nc+1)-spline%node_positions(nc-2))*(spline%node_positions(nc)-spline%node_positions(nc-1)))

      call compute_spline_nonunif_1D_hermite_aux(fp, nc, buf, ibuf, coeffs, lift)
   end subroutine compute_spline_nonunif_1D_hermite

!0.95165885066842626      -0.66689806247771388
!      0.16666666666666666       0.59333333333333371       0.23999999999999971
!    3.07692307692324215E-003  0.20571428571428574       0.79120879120879106
!    3.57142857142822359E-002  0.79761904761903923       0.16666666666667851
!
!      0.16666666666667851       0.79761904761903923       3.57142857142822359E-002
!      0.79120879120879106       0.20571428571428568       3.07692307692324129E-003
!      0.23999999999999971       0.59333333333333360       0.16666666666666666

   !>
   subroutine sll_s_setup_spline_nonunif_1d_periodic_aux(node_pos, N, buf, ibuf)
      sll_real64, dimension(:), pointer :: node_pos, buf
      sll_int32, intent(in) :: N
      sll_real64, dimension(:), pointer :: a, cts
      sll_int32, dimension(:), pointer  :: ipiv, ibuf
      sll_int32 :: i
      a => buf(1:3*N)
      cts => buf(3*N + 1:10*N)
      ipiv => ibuf(1:N)

      node_pos(-1) = node_pos(N - 1) - (node_pos(N) - node_pos(0))
      node_pos(-2) = node_pos(N - 2) - (node_pos(N) - node_pos(0))
      node_pos(N + 1) = node_pos(1) + (node_pos(N) - node_pos(0))
      node_pos(N + 2) = node_pos(2) + (node_pos(N) - node_pos(0))

      !fill a with mesh information
      do i = 0, N - 1
         !subdiagonal terms
      a(3*i+1)=(node_pos(i+1)-node_pos(i))*(node_pos(i+1)-node_pos(i))/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+1)-node_pos(i-2)))
         !a(3*i+1)=(node_pos(i+1)-node_pos(i))/(node_pos(i+1)-node_pos(i-1))*((node_pos(i+1)-node_pos(i))/(node_pos(i+1)-node_pos(i-2)))
         !overdiagonal terms
         !a(3*i+3)=(node_pos(i)-node_pos(i-1))*(node_pos(i)-node_pos(i-1))/((node_pos(i+1)-node_pos(i-1))*(node_pos(i+2)-node_pos(i-1)))
      a(3*i+3)=(node_pos(i)-node_pos(i-1))/(node_pos(i+1)-node_pos(i-1))*((node_pos(i)-node_pos(i-1))/(node_pos(i+2)-node_pos(i-1)))
         !diagonal terms
         a(3*i + 2) = 1.0_f64 - a(3*i + 1) - a(3*i + 3)
      end do
      !print *,'#a=',a
      !initialize the tridiagonal solver
      call sll_s_setup_cyclic_tridiag(a, N, cts, ipiv)
      !print *,'#cts=',maxval(cts),minval(cts)
      !print *,'#ipiv=',ipiv
   end subroutine sll_s_setup_spline_nonunif_1d_periodic_aux

   subroutine setup_spline_nonunif_1D_hermite_aux(node_pos, N, buf, ibuf)
      sll_real64, dimension(:), pointer :: node_pos, buf
      sll_int32, intent(in) :: N
      sll_real64, dimension(:), pointer :: a, cts
      sll_int32, dimension(:), pointer  :: ipiv, ibuf
      sll_int32 :: i, Np
      Np = N + 1
      a => buf(1:3*Np)
      cts => buf(3*Np + 1:10*Np)
      ipiv => ibuf(1:Np)

      !symmetric case
      node_pos(-1) = 2._f64*node_pos(0) - node_pos(1)
      node_pos(-2) = 2._f64*node_pos(0) - node_pos(2)
      node_pos(N + 1) = 2.0_f64*node_pos(N) - node_pos(N - 1)
      node_pos(N + 2) = 2.0_f64*node_pos(N) - node_pos(N - 2)

      !triple point case
      node_pos(-1) = node_pos(0)
      node_pos(-2) = node_pos(0)
      node_pos(N + 1) = node_pos(N)
      node_pos(N + 2) = node_pos(N)

      !fill a with mesh information
      !a(1)=0.0_f64
      !a(2)=(node_pos(2)-node_pos(0))/(node_pos(1)+node_pos(2)-2._f64*node_pos(0))
      !a(3)=1.0_f64-a(2)
      a(1) = 0.0_f64
      a(2) = (node_pos(2) - node_pos(0))/(node_pos(2) - node_pos(-1))
      a(3) = 1.0_f64 - a(2)
      do i = 1, N - 1
         !subdiagonal terms
         a(3*i + 1) = (node_pos(i + 1) - node_pos(i))**2/((node_pos(i + 1) - node_pos(i - 1))*(node_pos(i + 1) - node_pos(i - 2)))
         !overdiagonal terms
         a(3*i + 3) = (node_pos(i) - node_pos(i - 1))**2/((node_pos(i + 1) - node_pos(i - 1))*(node_pos(i + 2) - node_pos(i - 1)))
         !diagonal terms
         a(3*i + 2) = 1.0_f64 - a(3*i + 1) - a(3*i + 3)
      end do
      !a(3*N+2)=(node_pos(N)-node_pos(N-2))/(2._f64*node_pos(N)-node_pos(N-1)-node_pos(N-2))
      !a(3*N+1)=1.0_f64-a(3*N+2)
      !a(3*N+3)=0.0_f64
      a(3*N + 2) = (node_pos(N) - node_pos(N - 2))/(node_pos(N + 1) - node_pos(N - 2))
      a(3*N + 1) = 1.0_f64 - a(3*N + 2)
      a(3*N + 3) = 0.0_f64

      !initialize the tridiagonal solver
      call sll_s_setup_cyclic_tridiag(a, Np, cts, ipiv)

   end subroutine setup_spline_nonunif_1D_hermite_aux

   !>
   subroutine compute_spline_nonunif_1D_periodic_aux(f, N, buf, ibuf, coeffs)
      sll_real64, dimension(:), pointer :: f, buf, coeffs
      sll_int32, intent(in) :: N
      sll_real64, dimension(:), pointer :: cts!, a
      sll_int32, dimension(:), pointer  :: ipiv, ibuf
      !sll_real64 :: linf_err,tmp
      !a    => buf(1:3*N)
      cts => buf(3*N + 1:10*N)
      ipiv => ibuf(1:N)

      !compute the spline coefficients
      call sll_o_solve_cyclic_tridiag(cts, ipiv, f, N, coeffs(0:N - 1))
      coeffs(-1) = coeffs(N - 1)
      coeffs(N) = coeffs(0)
      coeffs(N + 1) = coeffs(1)

      !linf_err=0._f64
      !do i=1,N
      !  tmp=a(3*(i-1)+1)*coeffs(i-2)+a(3*(i-1)+2)*coeffs(i-1)+a(3*(i-1)+3)*coeffs(i)-f(i)
      !  if(abs(tmp)>linf_err)then
      !    linf_err=abs(tmp)
      !  endif
      !enddo
      !print *,'error of compute_spline=',linf_err
   end subroutine compute_spline_nonunif_1D_periodic_aux

   !>
   subroutine sll_s_compute_spline_nonunif_1d_periodic_aux2(f, N, buf, ibuf, coeffs)
      sll_real64, dimension(:) :: f
      sll_real64, dimension(:), pointer :: buf
      sll_real64, dimension(:), pointer :: coeffs
      sll_int32, intent(in) :: N
      sll_real64, dimension(:), pointer :: cts!, a
      sll_int32, dimension(:), pointer  :: ipiv, ibuf
      sll_real64, dimension(:), pointer :: a
      sll_real64, dimension(:), allocatable :: b
      sll_real64, dimension(:), pointer :: x
      sll_int32 :: ierr

      !sll_real64 :: linf_err,tmp
      !a    => buf(1:3*N)
      cts => buf(3*N + 1:10*N)
      ipiv => ibuf(1:N)

      a => buf(1:3*N)
      !b => f(1:N)
      x => coeffs(0:N - 1)

      SLL_ALLOCATE(b(N), ierr)
      b(1:N) = f(1:N)

      !compute the spline coefficients
      !call sll_o_solve_cyclic_tridiag( cts, ipiv, f, N, coeffs(0:N-1) )
      coeffs(0:N - 1) = 0._f64
      !do i=0,N-1
      !  if(coeffs(i) .ne. 0._f64)then
      !    print *,'#coeffs=',coeffs
      !    stop
      !  endif
      !enddo

      call sll_o_solve_cyclic_tridiag(cts, ipiv, f(1:N), N, coeffs(0:N - 1))

      coeffs(-1) = coeffs(N - 1)
      coeffs(N) = coeffs(0)
      coeffs(N + 1) = coeffs(1)

      !linf_err=0._f64
      !do i=1,N
      !  tmp=a(3*(i-1)+1)*coeffs(i-2)+a(3*(i-1)+2)*coeffs(i-1)+a(3*(i-1)+3)*coeffs(i)-f(i)
      !  if(abs(tmp)>linf_err)then
      !    linf_err=abs(tmp)
      !  endif
      !enddo
      !print *,'error of compute_spline=',linf_err
   end subroutine sll_s_compute_spline_nonunif_1d_periodic_aux2

   subroutine compute_spline_nonunif_1D_hermite_aux(f, N, buf, ibuf, coeffs, lift)
      sll_real64, dimension(:), pointer :: f, buf, coeffs
      sll_int32, intent(in) :: N
      sll_real64, intent(in) :: lift(4, 2)
      sll_real64, dimension(:), pointer :: cts
      sll_int32, dimension(:), pointer  :: ipiv, ibuf
      sll_int32 :: Np
      Np = N + 1
      cts => buf(3*Np + 1:10*Np)
      ipiv => ibuf(1:Np)

      !compute the spline coefficients with inplace tridiagonal solver
      coeffs(1:N - 1) = f(2:N)
      coeffs(0) = f(1) + lift(1, 1)
      coeffs(N) = f(N + 1) + lift(1, 2)
      call sll_o_solve_cyclic_tridiag(cts, ipiv, coeffs(0:N), Np, coeffs(0:N))
      !coeffs(-1) = coeffs(1)+lift(2,1)
      !coeffs(N+1) = coeffs(N-1)+lift(2,2)

      coeffs(-1) = lift(3, 1)*coeffs(0) + lift(4, 1)*coeffs(1) + lift(2, 1)
      coeffs(N + 1) = lift(3, 2)*coeffs(N - 1) + lift(4, 2)*coeffs(N) + lift(2, 2)

   end subroutine compute_spline_nonunif_1D_hermite_aux

   !> get spline interpolate at point x
   function interpolate_value_nonunif(x, spline)
      !warning not tested, should not work->t o be adapted like array
      sll_real64                                     :: interpolate_value_nonunif
      sll_real64, intent(in)                         :: x
      type(sll_t_cubic_nonunif_spline_1d), pointer      :: spline
      sll_int32 :: j
      sll_real64 ::xx
      sll_real64, dimension(:), pointer :: Xj, coef
      sll_real64 :: w(0:3)
      xx = (x - spline%xmin)/(spline%xmax - spline%xmin)
      if (.not. ((xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64))) then
         print *, 'bad_value of x=', x, 'xmin=', spline%xmin, 'xmax=', spline%xmax, xx, xx - 1.0_f64
         print *, 'in function interpolate_value_nonunif()'
         STOP
      end if
      !SLL_ASSERT( (xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64 ))

      !localization of xx
      if (xx == 1.0_f64) then
         j = spline%n_cells
      else
         j = 0
         do while (spline%node_positions(j) .le. xx)
            j = j + 1
         end do
      end if
      Xj => spline%node_positions(j:)

      SLL_ASSERT((xx .ge. Xj(0)) .and. (xx .le. Xj(1)))
      SLL_ASSERT((Xj(0) == spline%node_positions(j - 1)) .and. (Xj(1) == spline%node_positions(j)))

      !compute weights
      w(0) = (Xj(1) - xx)*(Xj(1) - xx)*(Xj(1) - xx)/((Xj(1) - Xj(0))*(Xj(1) - Xj(-1))*(Xj(1) - Xj(-2))); 
      w(1) = (Xj(1) - xx)*(Xj(1) - xx)*(xx - Xj(-2))/((Xj(1) - Xj(0))*(Xj(1) - Xj(-1))*(Xj(1) - Xj(-2))); 
      w(1) = w(1) + (Xj(2) - xx)*(Xj(1) - xx)*(xx - Xj(-1))/((Xj(1) - Xj(0))*(Xj(1) - Xj(-1))*(Xj(2) - Xj(-1))); 
      w(1) = w(1) + (Xj(2) - xx)*(Xj(2) - xx)*(xx - Xj(0))/((Xj(1) - Xj(0))*(Xj(2) - Xj(0))*(Xj(2) - Xj(-1))); 
      w(2) = (Xj(1) - xx)*(xx - Xj(-1))*(xx - Xj(-1))/((Xj(1) - Xj(0))*(Xj(1) - Xj(-1))*(Xj(2) - Xj(-1))); 
      w(2) = w(2) + (Xj(2) - xx)*(xx - Xj(-1))*(xx - Xj(0))/((Xj(1) - Xj(0))*(Xj(2) - Xj(0))*(Xj(2) - Xj(-1))); 
      w(2) = w(2) + (Xj(3) - xx)*(xx - Xj(0))*(xx - Xj(0))/((Xj(1) - Xj(0))*(Xj(2) - Xj(0))*(Xj(3) - Xj(0))); 
      w(3) = (xx - Xj(0))*(xx - Xj(0))*(xx - Xj(0))/((Xj(1) - Xj(0))*(Xj(2) - Xj(0))*(Xj(3) - Xj(0))); 
      coef => spline%coeffs(j - 1:)

      interpolate_value_nonunif = w(0)*coef(0) + w(1)*coef(1) + w(2)*coef(2) + w(3)*coef(3)
   end function interpolate_value_nonunif

   !> get spline interpolate on an array of points
   subroutine sll_s_interpolate_array_value_nonunif(a_in, a_out, n, spline)
      sll_int32, intent(in) :: n
      sll_real64, dimension(1:n), intent(in)  :: a_in
      sll_real64, dimension(1:n), intent(out) :: a_out
      sll_real64                         :: x
      type(sll_t_cubic_nonunif_spline_1d), pointer      :: spline
      sll_int32 :: i, j, shift = 3
      sll_real64 ::xx
      sll_real64, dimension(:), pointer :: Xj
      sll_real64, dimension(:), pointer :: coef
      sll_real64 :: w(4)

      x = a_in(1)
      xx = (x - spline%xmin)/(spline%xmax - spline%xmin)
      if (.not. ((xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64))) then
         print *, 'bad_value of x=', x, 'xmin=', spline%xmin, 'xmax=', spline%xmax, xx, xx - 1.0_f64
         print *, 'in subroutine sll_s_interpolate_array_value_nonunif()'
         STOP
      end if
      !SLL_ASSERT( (xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64 ))

      !localization of xx
      j = 0
      if (xx == 1.0_f64) then
         j = spline%n_cells
      else
         do while (spline%node_positions(j) .le. xx)
            j = j + 1
         end do
      end if

      do i = 1, n

         x = (a_in(i) - spline%xmin)/(spline%xmax - spline%xmin)
         if (.not. ((x .ge. 0.0_f64) .and. (x .le. 1.0_f64))) then
            print *, 'bad_value of a_in(', i, ')=', a_in(i), 'xmin=', spline%xmin, 'xmax=', spline%xmax
            print *, 'in subroutine sll_s_interpolate_array_value_nonunif()'
            STOP
         end if
         if (x == 1.0_f64) then
            j = spline%n_cells
         else
            if (x .ge. xx) then
               do while (spline%node_positions(j) .le. x)
                  j = j + 1
               end do
            else
               do while (spline%node_positions(j) .gt. x)
                  j = j - 1
               end do
               j = j + 1
            end if
         end if
         xx = x
         Xj => spline%node_positions(j - shift:)

         !print *,i,Xj(0+shift),xx,Xj(1+shift)
         !stop
         if (.not. ((xx .ge. Xj(0 + shift)) .and. (xx .lt. Xj(1 + shift)))) then
            if (xx .ne. 1.0_f64) then
               print *, Xj(0 + shift), xx, Xj(1 + shift)

               stop
            end if
         end if
         SLL_ASSERT((xx .ge. Xj(0 + shift)) .and. (xx .le. Xj(1 + shift)))
         SLL_ASSERT((Xj(0 + shift) == spline%node_positions(j - 1)) .and. (Xj(1 + shift) == spline%node_positions(j)))

         !compute weights
      w(1)=(Xj(shift+1)-xx)*(Xj(shift+1)-xx)*(Xj(shift+1)-xx)/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+1)-Xj(shift-2))); 
      w(2)=(Xj(shift+1)-xx)*(Xj(shift+1)-xx)*(xx-Xj(shift-2))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+1)-Xj(shift-2))); 
      w(2)=w(2)+(Xj(shift+2)-xx)*(Xj(shift+1)-xx)*(xx-Xj(shift-1))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+2)-Xj(shift-1))); 
      w(2)=w(2)+(Xj(shift+2)-xx)*(Xj(shift+2)-xx)*(xx-Xj(shift+0))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+2)-Xj(shift-1))); 
      w(3)=(Xj(shift+1)-xx)*(xx-Xj(shift-1))*(xx-Xj(shift-1))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+1)-Xj(shift-1))*(Xj(shift+2)-Xj(shift-1))); 
      w(3)=w(3)+(Xj(shift+2)-xx)*(xx-Xj(shift-1))*(xx-Xj(shift+0))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+2)-Xj(shift-1))); 
      w(3)=w(3)+(Xj(shift+3)-xx)*(xx-Xj(shift+0))*(xx-Xj(shift+0))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+3)-Xj(shift+0))); 
      w(4)=(xx-Xj(shift+0))*(xx-Xj(shift+0))*(xx-Xj(shift+0))/((Xj(shift+1)-Xj(shift+0))*(Xj(shift+2)-Xj(shift+0))*(Xj(shift+3)-Xj(shift+0))); 
         !coef => spline%coeffs(j-1:)
         !print *,i,xx,j,w(0),w(1),w(2),w(3),Xj(-2:3)
         !a_out(i) = w(0)*coef(0)+w(1)*coef(1)+w(2)*coef(2)+w(3)*coef(3)

         coef => spline%coeffs(j - 2:)
         !print *,i,xx,j,w(0),w(1),w(2),w(3),Xj(-2:3)
         a_out(i) = w(1)*coef(1) + w(2)*coef(2) + w(3)*coef(3) + w(4)*coef(4)

      end do
   end subroutine sll_s_interpolate_array_value_nonunif

   !>
   subroutine sll_s_interpolate_array_value_nonunif_aux(a_in, a_out, n, node_pos, coeffs, n_cells)
      sll_int32, intent(in) :: n, n_cells
      sll_real64, dimension(1:n), intent(in)  :: a_in
      !sll_real64, dimension(-1:n+1), intent(in)  :: node_pos
      sll_real64, dimension(1:n), intent(out) :: a_out
      sll_real64                         :: x
      !type(sll_t_cubic_nonunif_spline_1d), pointer      :: spline
      sll_int32 :: i, j, shift = 3
      sll_real64 ::xx
      sll_real64, dimension(:), pointer :: Xj
      sll_real64, dimension(:), pointer :: coef, coeffs, node_pos
      sll_real64 :: w(4)

      !do i=-1,n_cells+1
      !  print *,i,node_pos(i)
      !enddo

      xx = a_in(1)
      !xx = (x-spline%xmin)/(spline%xmax-spline%xmin)
      if (.not. ((xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64))) then
         print *, 'bad_value of x=', x!, 'xmin=', spline%xmin, 'xmax=', spline%xmax, xx,xx-1.0_f64
         print *, 'in subroutine sll_s_interpolate_array_value_nonunif()'
         STOP
      end if
      !SLL_ASSERT( (xx .ge. 0.0_f64) .and. (xx .le. 1.0_f64 ))

      !localization of xx
      j = 0
      if (xx == 1.0_f64) then
         j = n_cells
      else
         do while (node_pos(j) .le. xx)
            j = j + 1
         end do
      end if

      do i = 1, n

         x = a_in(i)
         if (.not. ((x .ge. 0.0_f64) .and. (x .le. 1.0_f64))) then
            print *, 'bad_value of a_in(', i, ')=', a_in(i)!, 'xmin=', spline%xmin, 'xmax=', spline%xmax
            print *, 'in subroutine sll_s_interpolate_array_value_nonunif()'
            STOP
         end if
         if (x == 1.0_f64) then
            j = n_cells!spline%n_cells
         else
            if (x .ge. xx) then
               do while (node_pos(j) .le. x)
                  j = j + 1
               end do
            else
               do while (node_pos(j) .gt. x)
                  j = j - 1
               end do
               j = j + 1
            end if
         end if
         xx = x
         Xj => node_pos(j - shift:)

         !print *,i,Xj(0+shift),xx,Xj(1+shift)
         !print *,i,Xj(-2+shift:3+shift)
         !stop
         if (.not. ((xx .ge. Xj(0 + shift)) .and. (xx .lt. Xj(1 + shift)))) then
            if (xx .ne. 1.0_f64) then
               print *, Xj(0 + shift), xx, Xj(1 + shift)

               stop
            end if
         end if
         !SLL_ASSERT( (xx .ge. Xj(0+shift)) .and. (xx .le. Xj(1+shift) ))
         !SLL_ASSERT( (Xj(0+shift)==spline%node_positions(j-1)) .and. (Xj(1+shift)==spline%node_positions(j)))

         !compute weights
         w(1) = (Xj(shift + 1) - xx)*(Xj(shift + 1) - xx)*(Xj(shift + 1) - xx)&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 1) - Xj(shift - 1))*(Xj(shift + 1) - Xj(shift - 2))); 
         w(2) = (Xj(shift + 1) - xx)*(Xj(shift + 1) - xx)*(xx - Xj(shift - 2))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 1) - Xj(shift - 1))*(Xj(shift + 1) - Xj(shift - 2))); 
         w(2) = w(2) + (Xj(shift + 2) - xx)*(Xj(shift + 1) - xx)*(xx - Xj(shift - 1))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 1) - Xj(shift - 1))*(Xj(shift + 2) - Xj(shift - 1))); 
         w(2) = w(2) + (Xj(shift + 2) - xx)*(Xj(shift + 2) - xx)*(xx - Xj(shift + 0))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift - 1))); 
         w(3) = (Xj(shift + 1) - xx)*(xx - Xj(shift - 1))*(xx - Xj(shift - 1))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 1) - Xj(shift - 1))*(Xj(shift + 2) - Xj(shift - 1))); 
         w(3) = w(3) + (Xj(shift + 2) - xx)*(xx - Xj(shift - 1))*(xx - Xj(shift + 0))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift - 1))); 
         w(3) = w(3) + (Xj(shift + 3) - xx)*(xx - Xj(shift + 0))*(xx - Xj(shift + 0))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift + 0))*(Xj(shift + 3) - Xj(shift + 0))); 
         w(4) = (xx - Xj(shift + 0))*(xx - Xj(shift + 0))*(xx - Xj(shift + 0))&
         &/((Xj(shift + 1) - Xj(shift + 0))*(Xj(shift + 2) - Xj(shift + 0))*(Xj(shift + 3) - Xj(shift + 0))); 
         !coef => spline%coeffs(j-1:)
         !print *,i,xx,j,w(1),w(2),w(3),w(4)!,Xj(-2:3)
         !a_out(i) = w(0)*coef(0)+w(1)*coef(1)+w(2)*coef(2)+w(3)*coef(3)

         coef => coeffs(j - 2:)
         a_out(i) = w(1)*coef(1) + w(2)*coef(2) + w(3)*coef(3) + w(4)*coef(4)
         !print *,i,xx,j,w(1),w(2),w(3),w(4),coef(1:4),a_out(i)
         !print *,Xj(shift-2:shift+3)
      end do
      !do i=1,n
      !  print *,i,a_in(i),a_out(i)
      !enddo
   end subroutine sll_s_interpolate_array_value_nonunif_aux

end module sll_m_cubic_non_uniform_splines

