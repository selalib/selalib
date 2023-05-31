!> @ingroup advection
!> @brief
!> Module implementing spline advection for the setting of a domain decomposition in 6d with extra buffers for the halo cells. The spline calculations are localized by an algorithm from signal processing.
!> @author
!> Katharina Kormann
!> Klaus Reuter, Max Planck Computing and Data Facility (MPCDF)

!  Note:  This module was (partly) rewritten such that Fortran range operators ":"
!  are avoided wherever possible.  These operators slow down the code a lot.

#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif

! ---
! Cache blocking greatly improves access to the large 6D arrays but heavily complicates the code.%
! For the interpolators below, it gives a performance boost of at least 2x.
! Performance needs to be measured with all the cores busy on a machine (MPI/OMP/hybrid),
! not only with a single thread!
! ---

module sll_m_advection_6d_spline_dd_slim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_constants, only : &
       sll_p_twopi

  use sll_m_interpolators_1d_base, only: &
       sll_c_interpolator_1d

  use sll_m_lagrange_interpolator_1d, only: &
       sll_p_lagrange_fixed, &
       sll_p_lagrange_centered

  use sll_m_decomposition, only: &
       sll_t_decomposition_slim_6d, &
       sll_t_cartesian_topology_6d, &
       sll_s_apply_halo_exchange_slim_6d_real64, &
       sll_f_apply_halo_exchange_slim_6d_real64, &
       sll_s_allocate_bc_buffers_6d, &
       sll_s_allocate_bc_buffers_6d_part, &
       sll_s_deallocate_bc_buffers, &
       sll_s_apply_bc_exchange_slim_6d_real64

  use sll_m_cubic_spline_halo_1d, only : &
    sll_s_cubic_spline_halo_1d_prepare_exchange, &
    sll_s_cubic_spline_halo_1d_finish_boundary_conditions, &
    sll_s_cubic_spline_halo_1d_compute_interpolant, &
    sll_s_cubic_spline_halo_1d_eval_disp

#ifdef _OPENMP
  use omp_lib
! TODO : confirm safety of collapse(2)
#define OMP_COLLAPSE collapse(2)
#define OMP_SCHEDULE schedule(static)
#endif

  implicit none

  ! --- spline routines ---
  public :: &
       sll_t_advection_6d_spline_dd_slim, &
       sll_s_advection_6d_spline_dd_slim_init, &
       sll_s_advection_6d_spline_dd_slim_free, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta1, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta2, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta3, &
       sll_s_advection_6d_spline_dd_slim_advect_eta4, &
       sll_s_advection_6d_spline_dd_slim_advect_eta5, &
       sll_s_advection_6d_spline_dd_slim_advect_eta6, &
       sll_s_advection_6d_spline_dd_slim_advect_eta1_dispeta45, &
       sll_s_advection_6d_spline_dd_slim_advect_eta2_dispeta45

  private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_t_advection_6d_spline_dd_slim
    sll_real64, allocatable :: displacement_eta1(:)
    sll_real64, allocatable :: displacement_eta2(:)
    sll_real64, allocatable :: displacement_eta3(:)

    sll_int32 , allocatable :: halo_blocks_eta1(:,:,:)
    sll_int32 , allocatable :: halo_blocks_eta2(:,:,:)
    sll_int32 , allocatable :: halo_blocks_eta3(:,:,:)

    sll_int32 , allocatable :: halo_blocks5d_eta2(:,:,:)
    sll_int32 , allocatable :: halo_blocks5d_eta3(:,:,:)

    sll_int32 , allocatable :: halo_width_eta1(:,:)
    sll_int32 , allocatable :: halo_width_eta2(:,:)
    sll_int32 , allocatable :: halo_width_eta3(:,:)

    sll_int32, allocatable :: idisplacement_eta1(:)
    sll_int32, allocatable :: idisplacement_eta2(:)
    sll_int32, allocatable :: idisplacement_eta3(:)

    sll_int32 :: n_halo_blocks(3)
  end type sll_t_advection_6d_spline_dd_slim

contains


  subroutine sll_s_advection_6d_spline_dd_slim_free(self)
    class(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !<  advector object

    if (allocated(self%displacement_eta1)) deallocate(self%displacement_eta1)
    if (allocated(self%displacement_eta2)) deallocate(self%displacement_eta2)
    if (allocated(self%displacement_eta3)) deallocate(self%displacement_eta3)

    if (allocated(self%halo_blocks_eta1)) deallocate(self%halo_blocks_eta1)
    if (allocated(self%halo_blocks_eta2)) deallocate(self%halo_blocks_eta2)
    if (allocated(self%halo_blocks_eta3)) deallocate(self%halo_blocks_eta3)

    if (allocated(self%halo_blocks5d_eta2)) deallocate(self%halo_blocks5d_eta2)
    if (allocated(self%halo_blocks5d_eta3)) deallocate(self%halo_blocks5d_eta3)

    if (allocated(self%halo_width_eta1)) deallocate(self%halo_width_eta1)
    if (allocated(self%halo_width_eta2)) deallocate(self%halo_width_eta2)
    if (allocated(self%halo_width_eta3)) deallocate(self%halo_width_eta3)

    if (allocated(self%idisplacement_eta1)) deallocate(self%idisplacement_eta1)
    if (allocated(self%idisplacement_eta2)) deallocate(self%idisplacement_eta2)
    if (allocated(self%idisplacement_eta3)) deallocate(self%idisplacement_eta3)
  end subroutine sll_s_advection_6d_spline_dd_slim_free


  subroutine sll_s_advection_6d_spline_dd_slim_init(self, decomposition, displacement_eta1, &
                                                          displacement_eta2, displacement_eta3)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !<  advector object
    type(sll_t_decomposition_slim_6d),  intent(in) :: decomposition

    sll_real64, intent(in) :: displacement_eta1(:)
    sll_real64, intent(in) :: displacement_eta2(:)
    sll_real64, intent(in) :: displacement_eta3(:)

    allocate(self%displacement_eta1(decomposition%local%mn(4):decomposition%local%mx(4)))
    self%displacement_eta1 = displacement_eta1
    allocate(self%displacement_eta2(decomposition%local%mn(5):decomposition%local%mx(5)))
    self%displacement_eta2 = displacement_eta2
    allocate(self%displacement_eta3(decomposition%local%mn(6):decomposition%local%mx(6)))
    self%displacement_eta3 = displacement_eta3

    call make_blocks_spline(4, decomposition, &
         self%displacement_eta1, self%idisplacement_eta1, self%halo_blocks_eta1, &
         self%halo_width_eta1, self%n_halo_blocks(1) )
    call make_blocks_spline(5, decomposition, &
         self%displacement_eta2, self%idisplacement_eta2, self%halo_blocks_eta2, &
         self%halo_width_eta2, self%n_halo_blocks(2) )
    call make_blocks_spline(6, decomposition, &
         self%displacement_eta3, self%idisplacement_eta3, self%halo_blocks_eta3, &
         self%halo_width_eta3, self%n_halo_blocks(3) )

    allocate(self%halo_blocks5d_eta2(5, 2, self%n_halo_blocks(2)))
    allocate(self%halo_blocks5d_eta3(5, 2, self%n_halo_blocks(3)))

    self%halo_blocks5d_eta2(1,:,:) = self%halo_blocks_eta2(1,:,:)
    self%halo_blocks5d_eta2(2:5,:,:) = self%halo_blocks_eta2(3:6,:,:)
    self%halo_blocks5d_eta3(1:2,:,:) = self%halo_blocks_eta3(1:2,:,:)
    self%halo_blocks5d_eta3(3:5,:,:) = self%halo_blocks_eta3(4:6,:,:)
#ifdef DISABLE_CACHE_BLOCKING
    write(*,*) "sll_m_advection_6d_spline_dd_slim :: cache blocking disabled"
#endif
  end subroutine sll_s_advection_6d_spline_dd_slim_init


  ! --- calculate depth of the buffers used for cache blocking ---
  function get_wx(decomposition, id_in)
    type(sll_t_decomposition_slim_6d), intent(in) :: decomposition
    sll_int32, intent(in), optional :: id_in
    sll_int32 :: get_wx, wx, id
#ifndef DISABLE_CACHE_BLOCKING
    ! --- by default we want to block in the first dimension
    if (present(id_in)) then
      id = id_in
    else
      id = 1
    endif
!    ! --- calculate width of prefetch cache buffer
!    wx = 8  ! SandyBridge cache line is 64 Bytes long
!    do while(wx > 1)
!      if (mod(decomposition%local%nw(id),wx) == 0) then
!        exit
!      else
!        wx = wx / 2
!      endif
!    enddo
    wx = decomposition%local%nw(id)
#else
    wx = 1
#endif
    get_wx = wx
  end function get_wx


  !> Helper function to calculate the communication blocks for the spline interpolation.
  subroutine make_blocks_spline( ind, decomposition,  disp, disp_int, halo_blocks , halo_width, n_halo_blocks)
    sll_int32,  intent( in    ) :: ind
    type(sll_t_decomposition_slim_6d), intent( in    ) :: decomposition
    sll_real64, intent( inout ) :: disp(decomposition%local%mn(ind):decomposition%local%mx(ind))
    sll_int32,  intent(   out ), allocatable :: disp_int(:)
    sll_int32,  intent(   out ), allocatable :: halo_blocks(:,:,:)
    sll_int32,  intent(   out ), allocatable :: halo_width(:,:)
    sll_int32,  intent(   out ) :: n_halo_blocks

    sll_int32 :: index_range(2)
    sll_int32 :: blocks, j, bl
    sll_int32 :: box1, box2
    sll_int32 :: si

    index_range = [ decomposition%local%mn(ind), decomposition%local%mx(ind) ]

    ! However, not all blocks might be present, so let us check first how many blocks there are or rather the interval of boxes that are there
    ! Since we assume monotonic displacement, we only need to identify the box of the first and last value
    bl = index_range(1)
    if ( abs(disp (bl) ) == 0.0_f64) bl = bl+1
    box1 = floor( disp(bl) )
    bl = index_range(2)
    if ( abs(disp(bl) ) == 0.0_f64) bl = bl-1
    box2 = floor( disp(bl) )

    ! Compute number of blocks
    blocks = abs(box2-box1)+1

    ! Now that we know the number of blocks, we can allocate the array holding the block information
    allocate( halo_blocks(6, 2, blocks) )
    allocate( halo_width (2, blocks) )
    allocate( disp_int( blocks ) )
    do j=1,blocks
       halo_blocks(:, 1, j) = decomposition%local%mn
       halo_blocks(:, 2, j) = decomposition%local%mx
    end do

    ! We have to distinguish increasing and decreasing displacement
    if (box1 > box2 ) then
       j = index_range(1)
       do bl = 1,blocks!box1, box2
          if ( abs(disp(j)) == 0.0_f64 ) j = j+1
          si = box1-bl+1
          halo_width(1, bl ) = -si!max( - si, 0 )!stencil/2 - (box1-bl+1)-1
          halo_width(2, bl ) = si+1!max( si+1, 0 )!stencil/2 + (box1-bl+1)
          disp_int( bl ) = si
          halo_blocks(ind, 1, bl ) = j
          do while( disp(j) > box1-bl+1 )
             j = j+1
             if ( j > index_range(2)) exit
          end do
          if ( abs(disp(j-1)) == 0.0_f64 ) then
             halo_blocks(ind, 2, bl ) = j-2
          else
             halo_blocks(ind, 2, bl ) = j-1
          end if
       end do
    else
       j = index_range(1)
       do bl=box1, box2
          if ( abs(disp(j)) == 0.0_f64 ) j = j+1
          halo_width(1, bl-box1+1) = -bl!max( - bl, 0 )!stencil/2 - bl-1
          halo_width(2, bl-box1+1) = bl+1!max( bl+1, 0 )!stencil/2 + bl
          disp_int( bl-box1+1 ) = bl
          halo_blocks(ind, 1, bl-box1+1) = j
          do while ( (disp(j) < bl+1) )
             j = j+1
             if ( j > index_range(2) ) exit
          end do
          if ( abs(disp(j-1)) == 0.0_f64 ) then
             halo_blocks(ind, 2, bl-box1+1) = j-2
          else
             halo_blocks(ind, 2, bl-box1+1) = j-1
          end if
       end do
    end if

    n_halo_blocks = blocks

    ! Compute normalized displacement
    do j= index_range(1), index_range(2)
       disp(j) = disp(j) - real( floor( disp(j)), f64)
       SLL_ASSERT( disp(j) >= 0.0_f64 )
       SLL_ASSERT( disp(j) <= 1.0_f64 )
    end do
  end subroutine make_blocks_spline


  !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta1(self, topology, decomposition,  f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx,bl
    sll_int32 :: num_points
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    num_points = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! NOTE: Buffering is done here for the second dim "j".
    !       This implementation uses 0-based indexing in w!
    wx = get_wx(decomposition, 2)

    do bl=1, self%n_halo_blocks(id)
       ! cache values to avoid multiple dereferencing inside the loops
       l_mn = c_mn-self%halo_width_eta1(1, bl)
       l_mx = c_mn-1
       r_mn = c_mx+1
       r_mx = c_mx+self%halo_width_eta1(2, bl)

       !call sll_s_allocate_bc_buffers_6d(decomposition, id)
       call sll_s_allocate_bc_buffers_6d_part(decomposition, id, &
            self%halo_blocks_eta1(2:6,1,bl), self%halo_blocks_eta1(2:6,2,bl))
!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2)
       allocate(d0(0:wx-1))
       allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       ! Compute the remote part of the boundary conditions
       do n=loop_mn(6), loop_mx(6)
          do m=loop_mn(5), loop_mx(5)
             do l=self%halo_blocks_eta1(id+3, 1, bl), self%halo_blocks_eta1(id+3, 2, bl)
                do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                   do j=loop_mn(2), loop_mx(2), wx
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                               f6d(:,j+w,k,l,m,n), self%idisplacement_eta1(bl), &
                               num_points, d0(w), c_np2(w))
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_left_send(j+w,k,l,m,n) = c_np2(w)
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_right_send(j+w,k,l,m,n) = d0(w)
                      enddo
                   end do
#else
                   do j=loop_mn(2), loop_mx(2)
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                             f6d(:,j,k,l,m,n), self%idisplacement_eta1(bl), &
                             num_points, d0(0), c_np2(0) )
                      decomposition%local%bc_left_send(j,k,l,m,n) = c_np2(0)
                      decomposition%local%bc_right_send(j,k,l,m,n) = d0(0)
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
      deallocate(c_np2)
      deallocate(d0)
!$omp single
       ! Exchange boundary conditions
       call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
       ! Exchange data for the neighboring cells
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta1(1, bl), &
            self%halo_width_eta1(2, bl), &
            self%halo_blocks_eta1(:,:,bl))
       l_halo => decomposition%local%halo_left%buf
       r_halo => decomposition%local%halo_right%buf
!$omp end single
       allocate(buf_i(l_mn-1:r_mx+1,0:wx-1))
       allocate(buf_o(l_mn-1:r_mx,0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       do n=loop_mn(6), loop_mx(6)
          do m=loop_mn(5), loop_mx(5)
             do l=self%halo_blocks_eta1(id+3, 1, bl), self%halo_blocks_eta1(id+3, 2, bl)
                do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                   do j=loop_mn(2), loop_mx(2), wx
!                      buf_i(c_mn:c_mx) = f6d(c_mn:c_mx,j,k,l,m,n)
                      do w=0,wx-1
                        do i=c_mn,c_mx
                          buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                        enddo
                      enddo
                      ! add local contributions to boundary conditions for the spline
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                               buf_i(c_mn:c_mx,w), &
                               self%idisplacement_eta1(bl), num_points, &
                               decomposition%local%bc_left(j+w,k,l,m,n), &
                               decomposition%local%bc_right(j+w,k,l,m,n))
                      enddo

                      ! fill input buffer piecewise
                      do w=0,wx-1
                        buf_i(l_mn-1,w) = decomposition%local%bc_left(j+w,k,l,m,n)!d_0
                      enddo

                      ! For splines, we always have halos and precisely one side

                      if ( l_mn <= l_mx ) then
!                       buf_i(l_mn:l_mx) = l_halo(:,j,k,l,m,n)
                        do w=0,wx-1
                          do i=l_mn,l_mx
                            buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                          enddo
                        enddo
                      else
!                       buf_i(r_mn:r_mx) = r_halo(:,j,k,l,m,n)
                        do w=0,wx-1
                          do i=r_mn,r_mx
                            buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                          enddo
                        enddo
                      end if
                      do w=0,wx-1
                        buf_i(r_mx+1,w) = decomposition%local%bc_right(j+w,k,l,m,n)!c_np2
                      enddo

                      ! compute interpolant
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                                buf_i(:,w), num_points, buf_o(:,w), buf_i(:,w) )
                      enddo

                      ! perform interpolation
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_eval_disp( &
                             buf_i(:,w), self%displacement_eta1(l), num_points, buf_o(:,w))
                      enddo

                      ! copy-back interpolated values
!                      f6d(:,j,k,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do w=0,wx-1
                        do i=0,num_points-1
                          f6d(c_mn+i,j+w,k,l,m,n) = buf_o(l_mn-1+i,w)
                        enddo
                      enddo
                   end do
#else
                   do j=loop_mn(2), loop_mx(2)
!                      buf_i(c_mn:c_mx) = f6d(c_mn:c_mx,j,k,l,m,n)
                      do i=c_mn,c_mx
                        buf_i(i,0) = f6d(i,j,k,l,m,n)
                      enddo
                      ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                           buf_i(c_mn:c_mx,0), &
                           self%idisplacement_eta1( bl ), num_points, &
                           decomposition%local%bc_left(j,k,l,m,n), &
                           decomposition%local%bc_right(j,k,l,m,n))
                      ! fill input buffer piecewise
                      buf_i(l_mn-1,0) = decomposition%local%bc_left(j,k,l,m,n)!d_0
                      ! For splines, we always have halos and precisely one side
                      if ( l_mn <= l_mx ) then
!                       buf_i(l_mn:l_mx) = l_halo(:,j,k,l,m,n)
                        do i=l_mn,l_mx
                          buf_i(i,0) = l_halo(i,j,k,l,m,n)
                        enddo
                      else
!                       buf_i(r_mn:r_mx) = r_halo(:,j,k,l,m,n)
                        do i=r_mn,r_mx
                          buf_i(i,0) = r_halo(i,j,k,l,m,n)
                        enddo
                      end if
                      buf_i(r_mx+1,0) = decomposition%local%bc_right(j,k,l,m,n)!c_np2
                      ! compute interpolant
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                             buf_i(:,0), num_points, buf_o(:,0), buf_i(:,0) )
                      ! perform interpolation
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_eval_disp( &
                             buf_i(:,0), self%displacement_eta1(l), num_points, buf_o(:,0))
                      ! copy-back interpolated values
!                      f6d(:,j,k,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do i=0,num_points-1
                        f6d(c_mn+i,j,k,l,m,n) = buf_o(l_mn-1+i,0)
                      enddo
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
      deallocate(buf_i)
      deallocate(buf_o)
!$omp end parallel

      call sll_s_deallocate_bc_buffers(decomposition)
    end do
  end subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta1


  !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta2(self, topology, decomposition,  f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 2  ! eta1
    sll_int32 :: i,j,k,l,m,n,bl,w
    sll_int32 :: num_points
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    num_points = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)

    do bl=1, self%n_halo_blocks(id)
       !call sll_s_allocate_bc_buffers_6d(decomposition, id)
       call sll_s_allocate_bc_buffers_6d_part(decomposition, id, &
            self%halo_blocks5d_eta2(:,1,bl), self%halo_blocks5d_eta2(:,2,bl))
       l_mn = c_mn-self%halo_width_eta2(1, bl)
       l_mx = c_mn-1
       r_mn = c_mx+1
       r_mx = c_mx+self%halo_width_eta2(2, bl)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2)
       allocate(buf_i(l_mn-1:r_mx+1,0:wx-1))
       allocate(buf_o(l_mn-1:r_mx,0:wx-1))
       allocate(d0(0:wx-1))
       allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       ! Compute the remote part of the boundary conditions
       do n=loop_mn(6), loop_mx(6)
          do m=self%halo_blocks_eta2(id+3, 1, bl), self%halo_blocks_eta2(id+3, 2, bl)
             do l=loop_mn(4), loop_mx(4)
                do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                   do i=loop_mn(1), loop_mx(1), wx
!                      buf_i(c_mn:c_mx) = f6d(i,:,k,l,m,n)
                      do j=c_mn,c_mx
                        do w=0,wx-1
                          buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                        enddo
                      enddo
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                                buf_i(c_mn:c_mx,w), self%idisplacement_eta2(bl), &
                                num_points, d0(w), c_np2(w) )
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_left_send(i+w,k,l,m,n) = c_np2(w)
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_right_send(i+w,k,l,m,n) = d0(w)
                      enddo
                   end do
#else
                   do i=loop_mn(1), loop_mx(1)
!                      buf_i(c_mn:c_mx) = f6d(i,:,k,l,m,n)
                      do j=c_mn,c_mx
                        buf_i(j,0) = f6d(i,j,k,l,m,n)
                      enddo
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                           buf_i(c_mn:c_mx,0), self%idisplacement_eta2( bl), &
                           num_points, d0(0), c_np2(0))
                      decomposition%local%bc_left_send(i,k,l,m,n) = c_np2(0)
                      decomposition%local%bc_right_send(i,k,l,m,n) = d0(0)
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
       deallocate(c_np2)
       deallocate(d0)
!$omp single
       ! Exchange boundary conditions
       call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
       ! Exchange data for the neighboring cells
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta2(1, bl), &
            self%halo_width_eta2(2, bl), &
            self%halo_blocks_eta2(:,:,bl))
       l_halo => decomposition%local%halo_left%buf
       r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       do n=loop_mn(6), loop_mx(6)
          do m=self%halo_blocks_eta2(id+3, 1, bl), self%halo_blocks_eta2(id+3, 2, bl)
             do l=loop_mn(4), loop_mx(4)
                do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                   do i=loop_mn(1), loop_mx(1), wx
!                      buf_i(c_mn:c_mx) = f6d(i,:,k,l,m,n)
                      do j=c_mn,c_mx
                        do w=0,wx-1
                          buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                        enddo
                      enddo
                      ! Add local contributions to boundary conditions for the spline
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                                buf_i(c_mn:c_mx,w), self%idisplacement_eta2( bl ), num_points, &
                                decomposition%local%bc_left(i+w,k,l,m,n), &
                                decomposition%local%bc_right(i+w,k,l,m,n) )
                      enddo
                      ! fill input buffer piecewise
                      do w=0,wx-1
                        buf_i(l_mn-1,w) = decomposition%local%bc_left(i+w,k,l,m,n)!d_0
                      enddo
                      ! For splines, we always have halos and precisely one side
                      if ( l_mn <= l_mx ) then
!                         buf_i(l_mn:l_mx) = l_halo(i,:,k,l,m,n)
                        do j=l_mn,l_mx
                          do w=0,wx-1
                            buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                          enddo
                        enddo
                      else
!                         buf_i(r_mn:r_mx) = r_halo(i,:,k,l,m,n)
                        do j=r_mn,r_mx
                          do w=0,wx-1
                            buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                          enddo
                        enddo
                      end if
                      do w=0,wx-1
                        buf_i(r_mx+1,w) = decomposition%local%bc_right(i+w,k,l,m,n)!c_np2
                      enddo
                      ! compute interpolant
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                             buf_i(:,w), num_points, buf_o(:,w), buf_i(:,w) )
                      enddo
                      ! perform interpolation
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_eval_disp( &
                              buf_i(:,w), self%displacement_eta2(m), num_points, buf_o(:,w))
                      enddo
                      ! copy-back interpolated values
!                      f6d(i,:,k,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do j=0,num_points-1
                        do w=0,wx-1
                          f6d(i+w,c_mn+j,k,l,m,n) = buf_o(l_mn-1+j,w)
                        enddo
                      enddo
                   end do
#else
                   do i=loop_mn(1), loop_mx(1)
!                      buf_i(c_mn:c_mx) = f6d(i,:,k,l,m,n)
                      do j=c_mn,c_mx
                        buf_i(j,0) = f6d(i,j,k,l,m,n)
                      enddo
                      ! Add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                           buf_i(c_mn:c_mx,0), self%idisplacement_eta2( bl ), num_points, &
                           decomposition%local%bc_left(i,k,l,m,n), &
                           decomposition%local%bc_right(i,k,l,m,n) )
                      ! fill input buffer piecewise
                      buf_i(l_mn-1,0) = decomposition%local%bc_left(i,k,l,m,n)!d_0
                      ! For splines, we always have halos and precisely one side
                      if ( l_mn <= l_mx ) then
!                         buf_i(l_mn:l_mx) = l_halo(i,:,k,l,m,n)
                        do j=l_mn,l_mx
                          buf_i(j,0) = l_halo(i,j,k,l,m,n)
                        enddo
                      else
!                         buf_i(r_mn:r_mx) = r_halo(i,:,k,l,m,n)
                        do j=r_mn,r_mx
                          buf_i(j,0) = r_halo(i,j,k,l,m,n)
                        enddo
                      end if
                      buf_i(r_mx+1,0) = decomposition%local%bc_right(i,k,l,m,n)!c_np2
                      ! compute interpolant
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                           buf_i(:,0), num_points, buf_o(:,0), buf_i(:,0) )
                      ! perform interpolation
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_eval_disp( &
                           buf_i(:,0), self%displacement_eta2(m), num_points, buf_o(:,0))
                      ! copy-back interpolated values
!                      f6d(i,:,k,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do j=0,num_points-1
                        f6d(i,c_mn+j,k,l,m,n) = buf_o(l_mn-1+j,0)
                      enddo
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
       deallocate(buf_i)
       deallocate(buf_o)
!$omp end parallel
       call sll_s_deallocate_bc_buffers(decomposition)
    end do
  end subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta2


  !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta3(self, topology, decomposition,  f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 3  ! eta1
    sll_int32 :: i,j,k,l,m,n,bl,w
    sll_int32 :: num_points
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    num_points = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)

    do bl=1, self%n_halo_blocks(id)
       !call sll_s_allocate_bc_buffers_6d(decomposition, id)
       call sll_s_allocate_bc_buffers_6d_part(decomposition, id, &
            self%halo_blocks5d_eta3(:,1,bl), self%halo_blocks5d_eta3(:,2,bl))
       l_mn = c_mn-self%halo_width_eta3(1, bl)
       l_mx = c_mn-1
       r_mn = c_mx+1
       r_mx = c_mx+self%halo_width_eta3(2, bl)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2)
       allocate(buf_i(l_mn-1:r_mx+1,0:wx-1))
       allocate(buf_o(l_mn-1:r_mx,0:wx-1))
       allocate(d0(0:wx-1))
       allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       ! Compute the remote part of the boundary conditions
       do n=self%halo_blocks_eta3(id+3, 1, bl), self%halo_blocks_eta3(id+3, 2, bl)
          do m=loop_mn(5), loop_mx(5)
             do l=loop_mn(4), loop_mx(4)
                do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                   do i=loop_mn(1), loop_mx(1), wx
!                      buf_i(c_mn:c_mx) = f6d(i,j,:,l,m,n)
                      do k=c_mn,c_mx
                        do w=0,wx-1
                          buf_i(k,w) = f6d(i+w,j,k,l,m,n)
                        enddo
                      enddo
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                                buf_i(c_mn:c_mx,w), self%idisplacement_eta3(bl), &
                                num_points, d0(w), c_np2(w))
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_left_send(i+w,j,l,m,n) = c_np2(w)
                      enddo
                      do w=0,wx-1
                        decomposition%local%bc_right_send(i+w,j,l,m,n) = d0(w)
                      enddo
                   end do
#else
                   do i=loop_mn(1), loop_mx(1)
!                      buf_i(c_mn:c_mx) = f6d(i,j,:,l,m,n)
                      do k=c_mn,c_mx
                        buf_i(k,0) = f6d(i,j,k,l,m,n)
                      enddo
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                              buf_i(c_mn:c_mx,0), self%idisplacement_eta3(bl), &
                              num_points, d0(0), c_np2(0) )
                      decomposition%local%bc_left_send(i,j,l,m,n) = c_np2(0)
                      decomposition%local%bc_right_send(i,j,l,m,n) = d0(0)
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
      deallocate(c_np2)
      deallocate(d0)
!$omp single
       ! Exchange boundary conditions
       call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
       ! Exchange data for the neighboring cells
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta3(1, bl), &
            self%halo_width_eta3(2, bl), &
            self%halo_blocks_eta3(:,:,bl))
       l_halo => decomposition%local%halo_left%buf
       r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
       do n=self%halo_blocks_eta3(id+3, 1, bl), self%halo_blocks_eta3(id+3, 2, bl)
          do m=loop_mn(5), loop_mx(5)
             do l=loop_mn(4), loop_mx(4)
                do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                   do i=loop_mn(1), loop_mx(1), wx
!                      buf_i(c_mn:c_mx) = f6d(i,j,:,l,m,n)
                      do k=c_mn,c_mx
                        do w=0,wx-1
                          buf_i(k,w) = f6d(i+w,j,k,l,m,n)
                        enddo
                      enddo
                      ! Add local contributions to boundary conditions for the spline
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                                buf_i(c_mn:c_mx,w), self%idisplacement_eta3( bl ), num_points, &
                                decomposition%local%bc_left(i+w,j,l,m,n), &
                                decomposition%local%bc_right(i+w,j,l,m,n) )
                      enddo
                      ! fill input buffer piecewise
                      do w=0,wx-1
                        buf_i(l_mn-1,w) = decomposition%local%bc_left(i+w,j,l,m,n)!d_0
                      enddo
                      ! For splines, we always have halos and precisely one side
                      if ( l_mn <= l_mx ) then
!                         buf_i(l_mn:l_mx) = l_halo(i,j,:,l,m,n)
                        do k=l_mn,l_mx
                          do w=0,wx-1
                            buf_i(k,w) = l_halo(i+w,j,k,l,m,n)
                          enddo
                        enddo
                      else
!                         buf_i(r_mn:r_mx) = r_halo(i,j,:,l,m,n)
                        do k=r_mn,r_mx
                          do w=0,wx-1
                            buf_i(k,w) = r_halo(i+w,j,k,l,m,n)
                          enddo
                        enddo
                      end if
                      do w=0,wx-1
                        buf_i(r_mx+1,w) = decomposition%local%bc_right(i+w,j,l,m,n)!c_np2
                      enddo
                      ! compute interpolant
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                                buf_i(:,w), num_points, buf_o(:,w), buf_i(:,w))
                      enddo

                      ! perform interpolation
                      do w=0,wx-1
!DIR$ FORCEINLINE
                        call sll_s_cubic_spline_halo_1d_eval_disp( &
                                buf_i(:,w), self%displacement_eta3(n), num_points, buf_o(:,w))
                      enddo
                      ! copy-back interpolated values
!                      f6d(i,j,:,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do k=0,num_points-1
                        do w=0,wx-1
                          f6d(i+w,j,c_mn+k,l,m,n) = buf_o(l_mn-1+k,w)
                        enddo
                      enddo
                   end do
#else
                   do i=loop_mn(1), loop_mx(1)
!                      buf_i(c_mn:c_mx) = f6d(i,j,:,l,m,n)
                      do k=c_mn,c_mx
                        buf_i(k,0) = f6d(i,j,k,l,m,n)
                      enddo
                      ! Add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                              buf_i(c_mn:c_mx,0), self%idisplacement_eta3( bl ), num_points, &
                              decomposition%local%bc_left(i,j,l,m,n), &
                              decomposition%local%bc_right(i,j,l,m,n) )
                      ! fill input buffer piecewise
                      buf_i(l_mn-1,0) = decomposition%local%bc_left(i,j,l,m,n)!d_0
                      ! For splines, we always have halos and precisely one side
                      if ( l_mn <= l_mx ) then
!                         buf_i(l_mn:l_mx) = l_halo(i,j,:,l,m,n)
                        do k=l_mn,l_mx
                          buf_i(k,0) = l_halo(i,j,k,l,m,n)
                        enddo
                      else
!                         buf_i(r_mn:r_mx) = r_halo(i,j,:,l,m,n)
                        do k=r_mn,r_mx
                          buf_i(k,0) = r_halo(i,j,k,l,m,n)
                        enddo
                      end if
                      buf_i(r_mx+1,0) = decomposition%local%bc_right(i,j,l,m,n)!c_np2
                      ! compute interpolant
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                              buf_i(:,0), num_points, buf_o(:,0), buf_i(:,0) )

                      ! perform interpolation
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_eval_disp( &
                              buf_i(:,0), self%displacement_eta3(n), num_points, buf_o(:,0))

                      ! copy-back interpolated values
!                      f6d(i,j,:,l,m,n) = buf_o(l_mn-1:r_mx-2)!(c_mn:c_mx) !TODO: First until last -3
                      do k=0,num_points-1
                        f6d(i,j,c_mn+k,l,m,n) = buf_o(l_mn-1+k,0)
                      enddo
                   end do
#endif
                end do
             end do
          end do
       end do
!$omp end do
       deallocate(buf_i)
       deallocate(buf_o)
!$omp end parallel
       call sll_s_deallocate_bc_buffers(decomposition)
    end do
  end subroutine sll_s_advection_6d_spline_dd_slim_fadvect_eta3


  !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_spline_dd_slim_advect_eta4(self, topology, decomposition, displacement, f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 4  ! eta1
    sll_int32 :: i,j,k,l,m,n,w
    sll_int32 :: c_mn
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: n_loc, ind_max
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32, allocatable :: si(:)
    sll_int32, allocatable :: indm(:)
    sll_real64, allocatable :: alpha(:)
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    call sll_s_allocate_bc_buffers_6d(decomposition, id)

    ! cache values to avoid multiple dereferencing inside the loops
    c_mn = decomposition%local%mn(id)
    n_loc = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ind_max = n_loc+3
    wx = get_wx(decomposition)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2,si,alpha,indm)
    ! WARNING: We use zero based indexing for the cache blocking to avoid numerous "+1"
    !          in indexing operations, different from the other advector routines!!!
    allocate(buf_i(1:ind_max, 0:wx-1))
    allocate(buf_o(1:ind_max-1, 0:wx-1))
    allocate(si(0:wx-1))
    allocate(indm(0:wx-1))
    allocate(alpha(0:wx-1))
    allocate(d0(0:wx-1))
    allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    ! Compute the remote part of the boundary conditions
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
                   do l=0,n_loc-1
                     do w=0,wx-1
                       buf_i(1+l,w) = f6d(i+w,j,k,c_mn+l,m,n)
                     enddo
                   enddo
                   do w=0,wx-1
                      si(w) = floor(displacement(i+w,j,k))
                   enddo
                   do w=0,wx-1
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                              buf_i(1:n_loc,w), si(w), n_loc, d0(w), c_np2(w) )
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_left_send(i+w,j,k,m,n) = c_np2(w)
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_right_send(i+w,j,k,m,n) = d0(w)
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
                   do l=0,n_loc-1
                     buf_i(1+l,0) = f6d(i,j,k,c_mn+l,m,n)
                   enddo
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                          buf_i(1:n_loc,0), si(0), n_loc, d0(0), c_np2(0) )
                   decomposition%local%bc_left_send(i,j,k,m,n) = c_np2(0)
                   decomposition%local%bc_right_send(i,j,k,m,n) = d0(0)
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(c_np2)
    deallocate(d0)
!$omp single
    ! Exchange boundary conditions
    ! TODO: restrict exchange to the actual range needed by the blocks
    call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
    ! Exchange data for the neighboring cells
    call sll_f_apply_halo_exchange_slim_6d_real64(topology, &
         decomposition, &
         f6d, &
         id, &
         1, &
         1)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   do w=0,wx-1
                     si(w) = floor( displacement(i+w,j,k) )
                   enddo
                   do w=0,wx-1
                     alpha(w) = displacement(i+w,j,k) - real(si(w), f64)
                   enddo
                   do w=0,wx-1
                     if ( si(w) == 0 ) then
                        indm(w) = 1
                        buf_i(n_loc+2,w) = r_halo(i+w,j,k,loop_mx(4)+1,m,n)
                     else ! si = -1
                        indm(w) = 2
                        buf_i(indm(w),w) = l_halo(i+w,j,k,c_mn-1,m,n)
                     end if
                   enddo
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do l=0,n_loc-1
                     do w=0,wx-1
                       buf_i(indm(w)+1+l,w) = f6d(i+w,j,k,c_mn+l,m,n)
                     enddo
                   enddo
                   ! add local contributions to boundary conditions for the spline
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                            buf_i(indm(w)+1:n_loc+indm(w),w), si(w), n_loc,  &
                            decomposition%local%bc_left(i+w,j,k,m,n), &
                            decomposition%local%bc_right(i+w,j,k,m,n) )
                   enddo
                   do w=0,wx-1
                     buf_i(1, w) = decomposition%local%bc_left(i+w,j,k,m,n)
                   enddo
                   do w=0,wx-1
                     buf_i(ind_max, w) = decomposition%local%bc_right(i+w,j,k,m,n)
                   enddo
                   ! compute interpolant
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                            buf_i(:,w), n_loc, buf_o(:,w), buf_i(:,w) )
                   enddo
                   ! perform interpolation
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_eval_disp( &
                            buf_i(:,w), alpha(w), n_loc, buf_o(:,w))
                   enddo
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do l=0,n_loc-1
                     do w=0,wx-1
                       f6d(i+w,j,k,c_mn+l,m,n) = buf_o(1+l,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
                   alpha(0) = displacement(i,j,k) - real(si(0), f64)
                   if ( si(0) == 0 ) then
                      indm(0) = 1
                      buf_i(n_loc+2,0) = r_halo(i,j,k,loop_mx(4)+1,m,n)
                   else ! si = -1
                      indm(0) = 2
                      buf_i(indm(0),0) = l_halo(i,j,k,c_mn-1,m,n)
                   end if
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do l=0,n_loc-1
                     buf_i(indm(0)+1+l,0) = f6d(i,j,k,c_mn+l,m,n)
                   enddo
                   ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                        buf_i(indm(0)+1:n_loc+indm(0),0), si(0), n_loc,  &
                        decomposition%local%bc_left(i,j,k,m,n), &
                        decomposition%local%bc_right(i,j,k,m,n) )
                   buf_i( 1, 0 ) = decomposition%local%bc_left(i,j,k,m,n)
                   buf_i( ind_max, 0 ) = decomposition%local%bc_right(i,j,k,m,n)
                   ! compute interpolant
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,0), n_loc, buf_o(:,0), buf_i(:,0) )
                   ! perform interpolation
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_eval_disp( &
                          buf_i(:,0), alpha(0), n_loc, buf_o(:,0))
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do l=0,n_loc-1
                     f6d(i,j,k,c_mn+l,m,n) = buf_o(1+l,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(buf_i)
    deallocate(buf_o)
    deallocate(si)
    deallocate(indm)
    deallocate(alpha)
!$omp end parallel
    call sll_s_deallocate_bc_buffers(decomposition)
  end subroutine sll_s_advection_6d_spline_dd_slim_advect_eta4


!> Advection along eta4 with displacement dependent on eta1-eta3
  subroutine sll_s_advection_6d_spline_dd_slim_advect_eta5(self, topology, decomposition, displacement, f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 5  ! eta1
    sll_int32 :: i,j,k,l,m,n,w
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: n_loc, ind_max
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32, allocatable :: si(:)
    sll_int32, allocatable :: indm(:)
    sll_real64, allocatable :: alpha(:)
    sll_int32 :: c_mn
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    call sll_s_allocate_bc_buffers_6d(decomposition, id)

    ! cache values to avoid multiple dereferencing inside the loops
    c_mn = decomposition%local%mn(id)
    n_loc = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ind_max = n_loc+3
    wx = get_wx(decomposition)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2,si,alpha,indm)
    allocate(buf_i(1:ind_max, 0:wx-1))
    allocate(buf_o(1:ind_max-1, 0:wx-1))
    allocate(si(0:wx-1))
    allocate(indm(0:wx-1))
    allocate(alpha(0:wx-1))
    allocate(d0(0:wx-1))
    allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    ! Compute the remote part of the boundary conditions
    do n=loop_mn(6), loop_mx(6)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   do w=0,wx-1
                     si(w) = floor( displacement(i+w,j,k) )
                   enddo
!                   buf_i(1:n_loc) = f6d(i,j,k,l,:,n)
                   do m=0,n_loc-1
                     do w=0,wx-1
                       buf_i(1+m,w) = f6d(i+w,j,k,l,c_mn+m,n)
                     enddo
                   enddo
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                            buf_i(1:n_loc,w), si(w), n_loc, d0(w), c_np2(w) )
                   enddo
                   do w=0,wx-1
                     decomposition%local%bc_left_send(i+w,j,k,l,n) = c_np2(w)
                   enddo
                   do w=0,wx-1
                     decomposition%local%bc_right_send(i+w,j,k,l,n) = d0(w)
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
!                   buf_i(1:n_loc) = f6d(i,j,k,l,:,n)
                   do m=0,n_loc-1
                     buf_i(1+m,0) = f6d(i,j,k,l,c_mn+m,n)
                   enddo
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                          buf_i(1:n_loc,0), si(0), n_loc, d0(0), c_np2(0) )
                   decomposition%local%bc_left_send(i,j,k,l,n) = c_np2(0)
                   decomposition%local%bc_right_send(i,j,k,l,n) = d0(0)
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(c_np2)
    deallocate(d0)
!$omp single
    ! Exchange boundary conditions
    ! TODO: restrict exchange to the actual range needed by the blocks
    call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
    ! Exchange data for the neighboring cells
    call sll_f_apply_halo_exchange_slim_6d_real64(topology, &
         decomposition, &
         f6d, &
         id, &
         1, &
         1)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   do w=0,wx-1
                     si(w) = floor( displacement(i+w,j,k) )
                   enddo
                   do w=0,wx-1
                     alpha(w) = displacement(i+w,j,k) - real(si(w), f64)
                   enddo
                   do w=0,wx-1
                     if ( si(w) == 0 ) then
                        indm(w) = 1
                        buf_i(n_loc+2,w) = r_halo(i+w,j,k,l,loop_mx(5)+1,n)
                     else ! si = -1
                        indm(w) = 2
                        buf_i(indm(w),w) = l_halo(i+w,j,k,l,c_mn-1,n)
                     end if
                   enddo
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,l,:,n)
                   do m=0,n_loc-1
                     do w=0,wx-1
                       buf_i(indm(w)+1+m,w) = f6d(i+w,j,k,l,c_mn+m,n)
                     enddo
                   enddo

                   ! add local contributions to boundary conditions for the spline
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                            buf_i(indm(w)+1:n_loc+indm(w),w), si(w), n_loc,  &
                            decomposition%local%bc_left(i+w,j,k,l,n), &
                            decomposition%local%bc_right(i+w,j,k,l,n) )
                   enddo
                   do w=0,wx-1
                     buf_i(1,w) = decomposition%local%bc_left(i+w,j,k,l,n)
                   enddo
                   do w=0,wx-1
                     buf_i(ind_max,w) = decomposition%local%bc_right(i+w,j,k,l,n)
                   enddo
                   ! compute interpolant
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                            buf_i(:,w), n_loc, buf_o(:,w), buf_i(:,w) )
                   enddo
                   ! perform interpolation
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_eval_disp( &
                            buf_i(:,w), alpha(w), n_loc, buf_o(:,w))
                   enddo
                   ! copy-back interpolated values
!                   f6d(i,j,k,l,:,n) = buf_o(1:n_loc)!First until last -3
                   do m=0,n_loc-1
                     do w=0,wx-1
                       f6d(i+w,j,k,l,c_mn+m,n) = buf_o(1+m,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
                   alpha(0) = displacement(i,j,k) - real(si(0), f64)
                   if ( si(0) == 0 ) then
                      indm(0) = 1
                      buf_i(n_loc+2,0) = r_halo(i,j,k,l,loop_mx(5)+1,n)
                   else ! si = -1
                      indm(0) = 2
                      buf_i(indm(0),0) = l_halo(i,j,k,l,c_mn-1,n)
                   end if

!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,l,:,n)
                   do m=0,n_loc-1
                     buf_i(indm(0)+1+m,0) = f6d(i,j,k,l,c_mn+m,n)
                   enddo

                   ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                          buf_i(indm(0)+1:n_loc+indm(0),0), si(0), n_loc,  &
                          decomposition%local%bc_left(i,j,k,l,n), &
                          decomposition%local%bc_right(i,j,k,l,n) )
                   buf_i(1,0) = decomposition%local%bc_left(i,j,k,l,n)
                   buf_i(ind_max,0) = decomposition%local%bc_right(i,j,k,l,n)
                   ! compute interpolant
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,0), n_loc, buf_o(:,0), buf_i(:,0) )
                   ! perform interpolation
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_eval_disp( &
                          buf_i(:,0), alpha(0), n_loc, buf_o(:,0))
                   ! copy-back interpolated values
!                   f6d(i,j,k,l,:,n) = buf_o(1:n_loc)!First until last -3
                   do m=0,n_loc-1
                     f6d(i,j,k,l,c_mn+m,n) = buf_o(1+m,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(buf_i)
    deallocate(buf_o)
    deallocate(si)
    deallocate(indm)
    deallocate(alpha)
!$omp end parallel
    call sll_s_deallocate_bc_buffers(decomposition)
  end subroutine sll_s_advection_6d_spline_dd_slim_advect_eta5


!> Advection along eta5 with displacement dependent on eta1-eta3
  subroutine sll_s_advection_6d_spline_dd_slim_advect_eta6(self, topology, decomposition, displacement,  f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 6  ! eta6
    sll_int32 :: i,j,k,l,m,n,w
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: n_loc, ind_max
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32, allocatable :: si(:)
    sll_int32, allocatable :: indm(:)
    sll_real64, allocatable :: alpha(:)
    sll_int32 :: c_mn
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    call sll_s_allocate_bc_buffers_6d(decomposition, id)

    ! cache values to avoid multiple dereferencing inside the loops
    c_mn = decomposition%local%mn(id)
    n_loc = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ind_max = n_loc+3
    wx = get_wx(decomposition)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2,si,alpha,indm)
    allocate(buf_i(1:ind_max, 0:wx-1))
    allocate(buf_o(1:ind_max-1, 0:wx-1))
    allocate(si(0:wx-1))
    allocate(indm(0:wx-1))
    allocate(alpha(0:wx-1))
    allocate(d0(0:wx-1))
    allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    ! Compute the remote part of the boundary conditions
    do m=loop_mn(5), loop_mx(5)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   do w=0,wx-1
                     si(w) = floor(displacement(i+w,j,k))
                   enddo
!                   buf_i(1:n_loc) = f6d(i,j,k,l,m,:)
                   do n=0,n_loc-1
                     do w=0,wx-1
                       buf_i(1+n,w) = f6d(i+w,j,k,l,m,c_mn+n)
                     enddo
                   enddo
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                            buf_i(1:n_loc,w), si(w), n_loc, d0(w), c_np2(w) )
                   enddo
                   do w=0,wx-1
                     decomposition%local%bc_left_send(i+w,j,k,l,m) = c_np2(w)
                   enddo
                   do w=0,wx-1
                     decomposition%local%bc_right_send(i+w,j,k,l,m) = d0(w)
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
!                   buf_i(1:n_loc) = f6d(i,j,k,l,m,:)
                   do n=0,n_loc-1
                     buf_i(1+n,0) = f6d(i,j,k,l,m,c_mn+n)
                   enddo
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                          buf_i(1:n_loc,0), si(0), n_loc, d0(0), c_np2(0) )
                   decomposition%local%bc_left_send(i,j,k,l,m) = c_np2(0)
                   decomposition%local%bc_right_send(i,j,k,l,m) = d0(0)
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(c_np2)
    deallocate(d0)
!$omp single
    ! Exchange boundary conditions
    ! TODO: restrict exchange to the actual range needed by the blocks
    call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
    ! Exchange data for the neighboring cells
    call sll_f_apply_halo_exchange_slim_6d_real64(topology, &
         decomposition, &
         f6d, &
         id, &
         1, &
         1)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do m=loop_mn(5), loop_mx(5)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   do w=0,wx-1
                     si(w) = floor( displacement(i+w,j,k) )
                   enddo
                   do w=0,wx-1
                     alpha(w) = displacement(i+w,j,k) - real(si(w), f64)
                   enddo
                   do w=0,wx-1
                     if ( si(w) == 0 ) then
                        indm(w) = 1
                        buf_i(n_loc+2,w) = r_halo(i+w,j,k,l,m,loop_mx(6)+1)
                     else ! si = -1
                        indm(w) = 2
                        buf_i(indm(w),w) = l_halo(i+w,j,k,l,m,c_mn-1)
                     end if
                   enddo

!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,l,m,:)
                   do n=0,n_loc-1
                     do w=0,wx-1
                       buf_i(indm(w)+1+n,w) = f6d(i+w,j,k,l,m,c_mn+n)
                     enddo
                   enddo

                   ! add local contributions to boundary conditions for the spline
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                            buf_i(indm(w)+1:n_loc+indm(w),w), si(w), n_loc,  &
                            decomposition%local%bc_left(i+w,j,k,l,m), &
                            decomposition%local%bc_right(i+w,j,k,l,m) )
                   enddo
                   do w=0,wx-1
                     buf_i(1,w) = decomposition%local%bc_left(i+w,j,k,l,m)
                   enddo
                   do w=0,wx-1
                     buf_i(ind_max,w) = decomposition%local%bc_right(i+w,j,k,l,m)
                   enddo

                   ! compute interpolant
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                            buf_i(:,w), n_loc, buf_o(:,w), buf_i(:,w))
                   enddo

                   ! perform interpolation
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_eval_disp( &
                            buf_i(:,w), alpha(w), n_loc, buf_o(:,w))
                   enddo

                   ! copy-back interpolated values
!                   f6d(i,j,k,l,m,:) = buf_o(1:n_loc)!First until last -3
                   do n=0,n_loc-1
                     do w=0,wx-1
                       f6d(i+w,j,k,l,m,c_mn+n) = buf_o(1+n,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si(0) = floor( displacement(i,j,k) )
                   alpha(0) = displacement(i,j,k) - real(si(0), f64)
                   if ( si(0) == 0 ) then
                      indm(0) = 1
                      buf_i(n_loc+2,0) = r_halo(i,j,k,l,m,loop_mx(6)+1)
                   else ! si = -1
                      indm(0) = 2
                      buf_i(indm(0),0) = l_halo(i,j,k,l,m,c_mn-1)
                   end if

!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,l,m,:)
                   do n=0,n_loc-1
                     buf_i(indm(0)+1+n,0) = f6d(i,j,k,l,m,c_mn+n)
                   enddo

                   ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                          buf_i(indm(0)+1:n_loc+indm(0),0), si(0), n_loc,  &
                          decomposition%local%bc_left(i,j,k,l,m), &
                          decomposition%local%bc_right(i,j,k,l,m) )

                   buf_i(1,0) = decomposition%local%bc_left(i,j,k,l,m)
                   buf_i(ind_max,0) = decomposition%local%bc_right(i,j,k,l,m)

                   ! compute interpolant
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,0), n_loc, buf_o(:,0), buf_i(:,0))

                   ! perform interpolation
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_eval_disp( &
                        buf_i(:,0), alpha(0), n_loc, buf_o(:,0))

                   ! copy-back interpolated values
!                   f6d(i,j,k,l,m,:) = buf_o(1:n_loc)!First until last -3
                   do n=0,n_loc-1
                     f6d(i,j,k,l,m,c_mn+n) = buf_o(1+n,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(buf_i)
    deallocate(buf_o)
    deallocate(si)
    deallocate(indm)
    deallocate(alpha)
!$omp end parallel
    call sll_s_deallocate_bc_buffers(decomposition)
  end subroutine sll_s_advection_6d_spline_dd_slim_advect_eta6


  
 !> Advection along eta2 with displacement dependent on eta4 and eta5 
  subroutine sll_s_advection_6d_spline_dd_slim_advect_eta2_dispeta45(self, topology, decomposition, displacement, f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w
    sll_int32 :: c_mn
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: n_loc, ind_max
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32 :: si
    sll_int32 :: indm
    sll_real64 :: alpha
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    call sll_s_allocate_bc_buffers_6d(decomposition, id)

    ! cache values to avoid multiple dereferencing inside the loops
    c_mn = decomposition%local%mn(id)
    n_loc = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ind_max = n_loc+3
    wx = get_wx(decomposition)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2,si,alpha,indm)
    ! WARNING: We use zero based indexing for the cache blocking to avoid numerous "+1"
    !          in indexing operations, different from the other advector routines!!!
    allocate(buf_i(1:ind_max, 0:wx-1))
    allocate(buf_o(1:ind_max-1, 0:wx-1))
    allocate(d0(0:wx-1))
    allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    ! Compute the remote part of the boundary conditions
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
             !do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
                   do j=0,n_loc-1
                     do w=0,wx-1
                       buf_i(1+j,w) = f6d(i+w,c_mn+j,k,l,m,n)
                     enddo
                   enddo
                   si = floor(displacement(l,m))
                   do w=0,wx-1
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                              buf_i(1:n_loc,w), si, n_loc, d0(w), c_np2(w) )
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_left_send(i+w,k,l,m,n) = c_np2(w)
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_right_send(i+w,k,l,m,n) = d0(w)
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   si = floor( displacement(l,m) )
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
                   do j=0,n_loc-1
                     buf_i(1+j,0) = f6d(i,c_mn+j,k,l,m,n)
                   enddo
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                          buf_i(1:n_loc,0), si, n_loc, d0(0), c_np2(0) )
                   decomposition%local%bc_left_send(i,k,l,m,n) = c_np2(0)
                   decomposition%local%bc_right_send(i,k,l,m,n) = d0(0)
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(c_np2)
    deallocate(d0)
!$omp single
    ! Exchange boundary conditions
    ! TODO: restrict exchange to the actual range needed by the blocks
    call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
    ! Exchange data for the neighboring cells
    call sll_f_apply_halo_exchange_slim_6d_real64(topology, &
         decomposition, &
         f6d, &
         id, &
         1, &
         1)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
             !do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   si = floor( displacement(l,m) )
                   if ( si == 0 ) then
                      indm = 1
                   else ! si = -1
                      indm = 2
                   end if
                   if ( si> 0 ) print*, '######## si', si
                   if (si<-1) print*, '####### si', si
                   alpha = displacement(l,m) - real(si, f64)
                   do w=0,wx-1
                     if ( si == 0 ) then
                        buf_i(n_loc+2,w) = r_halo(i+w,loop_mx(2)+1,k,l,m,n)
                     else ! si = -1
                        buf_i(indm,w) = l_halo(i+w,c_mn-1,k,l,m,n)
                     end if
                   enddo
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do j=0,n_loc-1
                     do w=0,wx-1
                       buf_i(indm+1+j,w) = f6d(i+w,c_mn+j,k,l,m,n)
                     enddo
                   enddo
                   ! add local contributions to boundary conditions for the spline
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                            buf_i(indm+1:n_loc+indm,w), si, n_loc,  &
                            decomposition%local%bc_left(i+w,k,l,m,n), &
                            decomposition%local%bc_right(i+w,k,l,m,n) )
                   enddo
                   do w=0,wx-1
                     buf_i(1, w) = decomposition%local%bc_left(i+w,k,l,m,n)
                   enddo
                   do w=0,wx-1
                     buf_i(ind_max, w) = decomposition%local%bc_right(i+w,k,l,m,n)
                   enddo
                   ! compute interpolant
                   do w=0,wx-1
                     ! print*, i,k,l,m,n
                     ! write(25,*) buf_i(:,w)
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,w), n_loc, buf_o(:,w), buf_i(:,w) )
                     !write(26,*) buf_i(:,w)
                     !stop
                  enddo
                   ! perform interpolation
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_eval_disp( &
                            buf_i(:,w), alpha, n_loc, buf_o(:,w))
                   enddo
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do j=0,n_loc-1
                     do w=0,wx-1
                       f6d(i+w,c_mn+j,k,l,m,n) = buf_o(1+j,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   
                   si = floor( displacement(l,m) )
                   if ( si == 0 ) then
                      indm = 1
                   else ! si = -1
                      indm = 2
                   end if
                   alpha = displacement(l,m) - real(si, f64)
                   if ( si == 0 ) then
                      buf_i(n_loc+2,0) = r_halo(i,loop_mx(2)+1,k,l,m,n)
                   else ! si = -1
                      buf_i(indm,0) = l_halo(i,c_mn-1,k,l,m,n)
                   end if
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do j=0,n_loc-1
                     buf_i(indm+1+j,0) = f6d(i,c_mn+j,k,l,m,n)
                   enddo
                   ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                        buf_i(indm+1:n_loc+indm,0), si, n_loc,  &
                        decomposition%local%bc_left(i,k,l,m,n), &
                        decomposition%local%bc_right(i,k,l,m,n) )
                   buf_i( 1, 0 ) = decomposition%local%bc_left(i,k,l,m,n)
                   buf_i( ind_max, 0 ) = decomposition%local%bc_right(i,k,l,m,n)
                   ! compute interpolant
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,0), n_loc, buf_o(:,0), buf_i(:,0) )
                   ! perform interpolation
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_eval_disp( &
                          buf_i(:,0), alpha, n_loc, buf_o(:,0))
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do j=0,n_loc-1
                     f6d(i,c_mn+j,k,l,m,n) = buf_o(1+j,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(buf_i)
    deallocate(buf_o)
!$omp end parallel
    call sll_s_deallocate_bc_buffers(decomposition)
  end subroutine sll_s_advection_6d_spline_dd_slim_advect_eta2_dispeta45

  
 !> Advection along eta1 with displacement dependent on eta4 and eta5 
  subroutine sll_s_advection_6d_spline_dd_slim_advect_eta1_dispeta45(self, topology, decomposition, displacement, f6d)
    type(sll_t_advection_6d_spline_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 1  ! eta2
    sll_int32 :: i,j,k,l,m,n,w
    sll_int32 :: c_mn
    sll_real64, allocatable :: buf_i(:,:)
    sll_real64, allocatable :: buf_o(:,:)
    sll_int32 :: n_loc, ind_max
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_int32 :: si
    sll_int32 :: indm
    sll_real64 :: alpha
    sll_int32 :: wx
    HALO_DTYPE, allocatable :: d0(:), c_np2(:)
    sll_int32, pointer :: loop_mn(:), loop_mx(:)

    call sll_s_allocate_bc_buffers_6d(decomposition, id)

    ! cache values to avoid multiple dereferencing inside the loops
    c_mn = decomposition%local%mn(id)
    n_loc = decomposition%local%nw(id)
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ind_max = n_loc+3
    wx = get_wx(decomposition)

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,d0,c_np2,si,alpha,indm)
    ! WARNING: We use zero based indexing for the cache blocking to avoid numerous "+1"
    !          in indexing operations, different from the other advector routines!!!
    allocate(buf_i(1:ind_max, 0:wx-1))
    allocate(buf_o(1:ind_max-1, 0:wx-1))
    allocate(d0(0:wx-1))
    allocate(c_np2(0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    ! Compute the remote part of the boundary conditions
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
                  ! do i=0,n_loc-1
                  !   do w=0,wx-1
                  !     buf_i(1+i,w) = f6d(c_mn+i,j+w,k,l,m,n)
                  !   enddo
                  ! enddo
                   si = floor(displacement(l,m))
                   if ( si> 0 ) print*, '######## si', si
                   if (si<-1) print*, '####### si', si
                   do w=0,wx-1
!DIR$ FORCEINLINE
                      call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                           f6d(:,j+w,k,l,m,n),  si, n_loc, d0(w), c_np2(w) )
                              !buf_i(1:n_loc,w), si, n_loc, d0(w), c_np2(w) )
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_left_send(j+w,k,l,m,n) = c_np2(w)
                   enddo
                   do w=0,wx-1
                      decomposition%local%bc_right_send(j+w,k,l,m,n) = d0(w)
                   enddo
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   si = floor( displacement(l,m) )
                   if ( si> 0 ) print*, '######## si', si
                   if (si<-1) print*, '####### si', si
!                   buf_i(1:n_loc) = f6d(i,j,k,:,m,n)
               !    do i=0,n_loc-1
               !      buf_i(1+i,0) = f6d(c_mn+i,j,k,l,m,n)
               !    enddo
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_prepare_exchange( &
                          f6d(:,j,k,l,m,n), si, n_loc, d0(0), c_np2(0) )
                         ! buf_i(1:n_loc,0), si, n_loc, d0(0), c_np2(0) )
                   decomposition%local%bc_left_send(j,k,l,m,n) = c_np2(0)
                   decomposition%local%bc_right_send(j,k,l,m,n) = d0(0)
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(c_np2)
    deallocate(d0)
!$omp single
    ! Exchange boundary conditions
    ! TODO: restrict exchange to the actual range needed by the blocks
    call sll_s_apply_bc_exchange_slim_6d_real64(topology, decomposition, id)
    ! Exchange data for the neighboring cells
    call sll_f_apply_halo_exchange_slim_6d_real64(topology, &
         decomposition, &
         f6d, &
         id, &
         1, &
         1)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
!$omp end single
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
             !do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
                   si = floor( displacement(l,m) )
                   if ( si == 0 ) then
                      indm = 1
                   else ! si = -1
                      indm = 2
                   end if
                   alpha = displacement(l,m) - real(si, f64)
                   do w=0,wx-1
                     if ( si == 0 ) then
                        buf_i(n_loc+2,w) = r_halo(loop_mx(1)+1,j+w,k,l,m,n)
                     else ! si = -1
                        buf_i(indm,w) = l_halo(c_mn-1,j+w,k,l,m,n)
                     end if
                   enddo
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do w=0,wx-1
                      do i=0,n_loc-1
                       buf_i(indm+1+i,w) = f6d(c_mn+i,j+w,k,l,m,n)
                     enddo
                   enddo
                   ! add local contributions to boundary conditions for the spline
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                            buf_i(indm+1:n_loc+indm,w), si, n_loc,  &
                            decomposition%local%bc_left(j+w,k,l,m,n), &
                            decomposition%local%bc_right(j+w,k,l,m,n) )
                   enddo
                   do w=0,wx-1
                     buf_i(1, w) = decomposition%local%bc_left(j+w,k,l,m,n)
                   enddo
                   do w=0,wx-1
                     buf_i(ind_max, w) = decomposition%local%bc_right(j+w,k,l,m,n)
                   enddo
                   ! compute interpolant
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                            buf_i(:,w), n_loc, buf_o(:,w), buf_i(:,w) )
                   enddo
                   ! perform interpolation
                   do w=0,wx-1
!DIR$ FORCEINLINE
                     call sll_s_cubic_spline_halo_1d_eval_disp( &
                            buf_i(:,w), alpha, n_loc, buf_o(:,w))
                   enddo
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do w=0,wx-1
                      do i=0,n_loc-1
                         f6d(c_mn+i,j+w,k,l,m,n) = buf_o(1+i,w)
                      enddo
                   enddo
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   si = floor( displacement(l,m) )
                   if ( si == 0 ) then
                      indm = 1
                   else ! si = -1
                      indm = 2
                   end if
                   alpha = displacement(l,m) - real(si, f64)
                   if ( si == 0 ) then
                      buf_i(n_loc+2,0) = r_halo(loop_mx(1)+1,j,k,l,m,n)
                   else ! si = -1
                      buf_i(indm,0) = l_halo(c_mn-1,j,k,l,m,n)
                   end if
!                   buf_i(indm+1:n_loc+indm) = f6d(i,j,k,:,m,n)
                   do i=0,n_loc-1
                     buf_i(indm+1+i,0) = f6d(c_mn+i,j,k,l,m,n)
                   enddo
                   ! add local contributions to boundary conditions for the spline
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_finish_boundary_conditions( &
                        buf_i(indm+1:n_loc+indm,0), si, n_loc,  &
                        decomposition%local%bc_left(j,k,l,m,n), &
                        decomposition%local%bc_right(j,k,l,m,n) )
                   buf_i( 1, 0 ) = decomposition%local%bc_left(j,k,l,m,n)
                   buf_i( ind_max, 0 ) = decomposition%local%bc_right(j,k,l,m,n)
                   ! compute interpolant
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_compute_interpolant( &
                          buf_i(:,0), n_loc, buf_o(:,0), buf_i(:,0) )
                   ! perform interpolation
!DIR$ FORCEINLINE
                   call sll_s_cubic_spline_halo_1d_eval_disp( &
                          buf_i(:,0), alpha, n_loc, buf_o(:,0))
                   ! copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(1:n_loc)!First until last -3
                   do i=0,n_loc-1
                     f6d(c_cm+i,j,k,l,m,n) = buf_o(1+i,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
    deallocate(buf_i)
    deallocate(buf_o)
!$omp end parallel
    call sll_s_deallocate_bc_buffers(decomposition)
  end subroutine sll_s_advection_6d_spline_dd_slim_advect_eta1_dispeta45

  
end module sll_m_advection_6d_spline_dd_slim
