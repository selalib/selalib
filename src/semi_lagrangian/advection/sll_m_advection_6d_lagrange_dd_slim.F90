#ifdef USE_HALO_REAL32
#define HALO_DTYPE sll_real32
#else
#define HALO_DTYPE sll_real64
#endif

! --- confirmed 20160825: cache blocking ON gives the same result as OFF ---

module sll_m_advection_6d_lagrange_dd_slim
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
  
  use sll_m_distribution_function_initializer_6d, only: &
       sll_t_array

  use sll_m_lagrange_interpolation_1d_fast, only : &
       sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells, &
       sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells, &
       sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells, &
       sll_s_lagrange_interpolation_1d_fast_haloc_cells

#ifdef _OPENMP
  use omp_lib
! TODO : confirm safety of collapse(2)
#define OMP_COLLAPSE collapse(2)
#define OMP_SCHEDULE schedule(static)
#endif

#ifdef USE_FMEMPOOL
  use fmempool
#endif

  implicit none

  ! --- Lagrange routines ---
  public :: &
       sll_t_advection_6d_lagrange_dd_slim, &
       sll_s_advection_6d_lagrange_dd_slim_init, &
       sll_s_advection_6d_lagrange_dd_slim_free, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta4, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta5, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta6, &
       sll_s_advection_6d_lagrange_dd_slim_set_eta123, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta1, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta2, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1_dispeta45, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2_dispeta45, &
       sll_s_advection_6d_clagrange_dd_slim_advect_eta1_d45, &
       sll_s_advection_6d_clagrange_dd_slim_advect_eta2_d45, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1_givenv, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2_givenv, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta3_givenv

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_t_advection_6d_lagrange_dd_slim
     sll_int32 :: lagrange_width(2)

    sll_real64, allocatable :: displacement_eta1(:)
    sll_real64, allocatable :: displacement_eta2(:)
    sll_real64, allocatable :: displacement_eta3(:)

    sll_int32 , allocatable :: halo_blocks_eta1(:,:,:)
    sll_int32 , allocatable :: halo_blocks_eta2(:,:,:)
    sll_int32 , allocatable :: halo_blocks_eta3(:,:,:)

    sll_int32 , allocatable :: halo_width_eta1(:,:)
    sll_int32 , allocatable :: halo_width_eta2(:,:)
    sll_int32 , allocatable :: halo_width_eta3(:,:)

    sll_int32 :: n_halo_blocks(3)
  end type sll_t_advection_6d_lagrange_dd_slim


contains


  subroutine sll_s_advection_6d_lagrange_dd_slim_free(self)
    class(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !<  advector object

    if (allocated(self%displacement_eta1)) deallocate(self%displacement_eta1)
    if (allocated(self%displacement_eta2)) deallocate(self%displacement_eta2)
    if (allocated(self%displacement_eta3)) deallocate(self%displacement_eta3)

    if (allocated(self%halo_blocks_eta1)) deallocate(self%halo_blocks_eta1)
    if (allocated(self%halo_blocks_eta2)) deallocate(self%halo_blocks_eta2)
    if (allocated(self%halo_blocks_eta3)) deallocate(self%halo_blocks_eta3)

    if (allocated(self%halo_width_eta1)) deallocate(self%halo_width_eta1)
    if (allocated(self%halo_width_eta2)) deallocate(self%halo_width_eta2)
    if (allocated(self%halo_width_eta3)) deallocate(self%halo_width_eta3)
  end subroutine sll_s_advection_6d_lagrange_dd_slim_free


  subroutine sll_s_advection_6d_lagrange_dd_slim_init(self, lagrange_width )
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !<  advector object
    sll_int32, intent( in ) :: lagrange_width(2)
    self%lagrange_width = lagrange_width
#ifdef DISABLE_CACHE_BLOCKING
    write(*,*) "sll_m_advection_6d_lagrange_dd_slim :: cache blocking disabled"
#endif
  end subroutine sll_s_advection_6d_lagrange_dd_slim_init


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
    ! --- calculate width of prefetch cache buffer
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



! --- SETUP ROUTINES FOR LAGRANGE INTERPOLATION WITH VARIABLE DISPLACEMENT BELOW ---

  subroutine sll_s_advection_6d_lagrange_dd_slim_set_eta123(self, decomposition, displacement_eta1, &
                                                            displacement_eta2, displacement_eta3)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !<  advector object
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

    call make_blocks_lagrange(4, self%lagrange_width(1), decomposition, &
         self%displacement_eta1, self%halo_blocks_eta1, &
         self%halo_width_eta1, self%n_halo_blocks(1))
    call make_blocks_lagrange(5, self%lagrange_width(1), decomposition, &
         self%displacement_eta2, self%halo_blocks_eta2, &
         self%halo_width_eta2, self%n_halo_blocks(2))
    call make_blocks_lagrange(6, self%lagrange_width(1), decomposition, &
         self%displacement_eta3, self%halo_blocks_eta3, &
         self%halo_width_eta3, self%n_halo_blocks(3))
  end subroutine sll_s_advection_6d_lagrange_dd_slim_set_eta123


  !> Helper function to calculate the communication blocks for the Lagrange interpolation.
  subroutine make_blocks_lagrange(ind, stencil, decomposition, disp, halo_blocks , halo_width, n_halo_blocks)
    sll_int32,  intent( in    ) :: ind
    sll_int32,  intent( in    ) :: stencil
    type(sll_t_decomposition_slim_6d), intent( in    ) :: decomposition
    sll_real64, intent( inout ) :: disp(decomposition%local%mn(ind):decomposition%local%mx(ind))
    sll_int32,  intent(   out ), allocatable :: halo_blocks(:,:,:)
    sll_int32,  intent(   out ), allocatable :: halo_width(:,:)
    sll_int32,  intent(   out ) :: n_halo_blocks

    sll_int32 :: index_range(2)
    sll_int32 :: blocks, j, bl
    sll_int32 :: box1, box2

    index_range = [ decomposition%local%mn(ind), decomposition%local%mx(ind) ]

    ! Determine the number of blocks: We want to allow for displacements of maximally
    ! stencil/2 cells (such that the total number of halo cells can be stencil-1 in each
    ! of the blocks.
    ! For a stencil of width "stencil" (even) the maximum number of blocks is
    ! stencil
    blocks = stencil
    ! However, not all blocks might be present, so let us check first how many blocks there are or rather the interval of boxes that are there
    ! Since we assume monotonic displacement, we only need to identify the box of the first and last value
    bl = index_range(1)
    if ( abs(disp (bl) ) == 0.0_f64) bl = bl+1
    box1 = floor( disp(bl) )
    bl = index_range(2)
    if ( abs(disp(bl) ) == 0.0_f64) bl = bl-1
    box2 = floor( disp(bl) )

    ! The allowed boxes in this numbering are -stencil/2, ...., stencil/2-1 (otherwise the displacement was too large in modulus)
    SLL_ASSERT( box1 >= -stencil/2 )
    SLL_ASSERT( box1 <   stencil/2 )
    SLL_ASSERT( box2 >= -stencil/2 )
    SLL_ASSERT( box2 <   stencil/2 )
    ! Compute number of blocks
    blocks = abs(box2-box1)+1

    ! Now that we know the number of blocks, we can allocate the array holding the block information
    allocate( halo_blocks(6, 2, blocks) )
    allocate( halo_width (2, blocks) )
    do j=1,blocks
       halo_blocks(:, 1, j) = decomposition%local%mn
       halo_blocks(:, 2, j) = decomposition%local%mx
    end do

    ! We have to distinguish increasing and decreasing displacement
    if (box1 > box2 ) then
       j = index_range(1)
       do bl = 1,blocks!box1, box2
          if ( abs(disp(j)) == 0.0_f64 ) j = j+1
          halo_width(1, bl ) = stencil/2 - (box1-bl+1)-1
          halo_width(2, bl ) = stencil/2 + (box1-bl+1)
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
          halo_width(1, bl-box1+1) = stencil/2 - bl-1
          halo_width(2, bl-box1+1) = stencil/2 + bl
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
  end subroutine make_blocks_lagrange



  !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta1(self, topology, decomposition,  f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
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
    sll_int32 :: j

    do j=1, self%n_halo_blocks(id)
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta1(1, j), &
            self%halo_width_eta1(2, j), &
            self%halo_blocks_eta1(:,:,j))

       call sll_s_advection_6d_lagrange_dd_slim_advect_eta1_core(&
            self, &
            decomposition, &
            self%halo_blocks_eta1(id+3, :, j), &
            self%displacement_eta1, &
            f6d)
    end do
  end subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta1

 !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_core(self, decomposition, eta4_cut, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_int32,  intent(in) :: eta4_cut(2)
    sll_real64, intent(in) :: displacement(decomposition%local%mn(4):decomposition%local%mx(4)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! It turns out that blocking over the second dimension is not significantly faster,
    ! even slower that no blocking, so be disable it by default!
#ifdef CACHE_BLOCKING_J
    ! NOTE: Buffering is done here for the second dim "j".
    !       This implementation uses 0-based indexing in w!
    wx = get_wx(decomposition, 2)
#else
    wx = 1
#endif
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=eta4_cut(1), eta4_cut(2)
             do k=loop_mn(3), loop_mx(3)
#ifdef CACHE_BLOCKING_J
                do j=loop_mn(2), loop_mx(2), wx
                   do w=0,wx-1
                     do i=l_mn, l_mx
                       buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                     enddo
                   enddo
                   do w=0,wx-1
                     do i=c_mn, c_mx
                       buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                     enddo
                   enddo
                   do w=0,wx-1
                     do i=r_mn, r_mx
                        buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                     enddo
                   enddo
                   ! perform interpolation
                   do w=0,wx-1
                     call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                            buf_i(:,w), &
                            buf_o(:,w), &
                            displacement(l), &
                            lw)
                   enddo
                   ! copy-back interpolated values
!                   f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx)
                   do w=0,wx-1
                     do i=c_mn, c_mx
                       f6d(i,j+w,k,l,m,n) = buf_o(i,w)
                     enddo
                   enddo
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   ! (1) fill input buffer piecewise
!                   if ( l_mn <= l_mx ) then
!                      buf_i(l_mn:l_mx) = l_halo(:,j,k,l,m,n)
!                   end if
!                   buf_i(c_mn:c_mx) = f6d(:,j,k,l,m,n)
!                   if (r_mn <= r_mx ) then
!                      buf_i(r_mn:r_mx) = r_halo(:,j,k,l,m,n)
!                   end if
                   do i=l_mn, l_mx
                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do i=c_mn, c_mx
                     buf_i(i,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do i=r_mn, r_mx
                      buf_i(i,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(l), &
                          lw)
                   ! (3) copy-back interpolated values
!                   f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx)
                   do i=c_mn, c_mx
                     f6d(i,j,k,l,m,n) = buf_o(i,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_core


  !> Advection along eta2 with displacement dependent on eta5
  subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta2(self, topology, decomposition, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: j

    do j=1, self%n_halo_blocks(id)
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta2(1, j), &
            self%halo_width_eta2(2, j), &
            self%halo_blocks_eta2(:,:,j))

       call sll_s_advection_6d_lagrange_dd_slim_advect_eta2_core(&
            self, &
            decomposition, &
            self%halo_blocks_eta2(id+3, :, j), &
            self%displacement_eta2, &
            f6d)
    end do
  end subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta2

  !> Advection along eta2 with displacement dependent on eta5
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_core(self, decomposition, eta5_cut, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_int32,  intent(in) :: eta5_cut(2)
    sll_real64, intent(in) :: displacement(decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=eta5_cut(1), eta5_cut(2)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   ! (1) fill input buffer piecewise
                   do j=l_mn,l_mx
                     do w=0,wx-1
                       buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   do j=c_mn,c_mx
                     do w=0,wx-1
                       buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   do j=r_mn,r_mx
                     do w=0,wx-1
                       buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   ! (2) perform interpolation
                   do w=0,wx-1
                     call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                            buf_i(:,w), &
                            buf_o(:,w), &
                            displacement(m), &
                            lw)
                   enddo
                   ! (3) copy-back interpolated values
!                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx)
                   do j=c_mn,c_mx
                     do w=0,wx-1
                       f6d(i+w,j,k,l,m,n) = buf_o(j,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
                   do j=l_mn,l_mx
                     buf_i(j,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do j=c_mn,c_mx
                     buf_i(j,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do j=r_mn,r_mx
                     buf_i(j,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(m), &
                          lw)
                   ! (3) copy-back interpolated values
!                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx)
                   do j=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(j,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_core


  !> Advection along eta3 with displacement dependent on eta6
  subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3(self, topology, decomposition, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_cartesian_topology_6d), intent( in ) :: topology
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output

    sll_int32, parameter :: id = 3  ! eta3
    sll_int32 :: j

    do j=1, self%n_halo_blocks(id)
       call sll_s_apply_halo_exchange_slim_6d_real64(topology, &
            decomposition, &
            f6d, &
            id, &
            self%halo_width_eta3(1, j), &
            self%halo_width_eta3(2, j), &
            self%halo_blocks_eta3(:,:,j))

       call sll_s_advection_6d_lagrange_dd_slim_advect_eta3_core(&
            self, &
            decomposition, &
            self%halo_blocks_eta3(id+3, :, j), &
            self%displacement_eta3, &
            f6d)
    end do
  end subroutine sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3

    !> Advection along eta3 with displacement dependent on eta6
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3_core(self, decomposition, eta6_cut, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_int32, intent(in) :: eta6_cut(2)
    sll_real64, intent(in) :: displacement(decomposition%local%mn(6):decomposition%local%mx(6)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 3  ! eta3
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=eta6_cut(1), eta6_cut(2)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                   ! (1) fill input buffer piecewise
                   do k=l_mn,l_mx
                     do w=0,wx-1
                       buf_i(k,w) = l_halo(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   do k=c_mn,c_mx
                     do w=0,wx-1
                       buf_i(k,w) = f6d(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   do k=r_mn,r_mx
                     do w=0,wx-1
                       buf_i(k,w) = r_halo(i+w,j,k,l,m,n)
                     enddo
                   enddo
                   ! (2) perform interpolation
                   do w=0,wx-1
                     call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                            buf_i(:,w), &
                            buf_o(:,w), &
                            displacement(n), &
                            lw)
                   enddo
                   ! (3) copy-back interpolated values
!                   f6d(i,j,:,l,m,n) = buf_o(c_mn:c_mx)
                   do k=c_mn,c_mx
                     do w=0,wx-1
                       f6d(i+w,j,k,l,m,n) = buf_o(k,w)
                     enddo
                   enddo
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
                   do k=l_mn,l_mx
                     buf_i(k,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do k=c_mn,c_mx
                     buf_i(k,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do k=r_mn,r_mx
                     buf_i(k,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_centered_halo_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(n), &
                          lw)
                   ! (3) copy-back interpolated values
!                   f6d(i,j,:,l,m,n) = buf_o(c_mn:c_mx)
                   do k=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(k,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3_core


! --- LAGRANGE INTERPOLATION ROUTINES WITH *FIXED* DISCPLACEMENT BELOW ---
!     (previously sll_m_advection_6d_lagrange_dd_slim.BAK.F90)


 !> Advection along eta1 with displacement dependent on eta4
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(4):decomposition%local%mx(4)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! NOTE: Buffering is done here for the second dim "j".
    !       Seems to be marginally beneficial.
    wx = get_wx(decomposition, 2)
    lw = self%lagrange_width(1)

! !$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
!     allocate(buf_i(l_mn:r_mx))
!     allocate(buf_o(l_mn:r_mx))
! !$omp do OMP_COLLAPSE OMP_SCHEDULE
!     do n=decomposition%local%mn(6), decomposition%local%mx(6)
!        do m=decomposition%local%mn(5), decomposition%local%mx(5)
!           do l=decomposition%local%mn(4), decomposition%local%mx(4)
!              do k=decomposition%local%mn(3), decomposition%local%mx(3)
!                 do j=decomposition%local%mn(2), decomposition%local%mx(2)
!                    ! (1) fill input buffer piecewise
!                    buf_i(l_mn:l_mx) = l_halo(:,j,k,l,m,n)
!                    buf_i(c_mn:c_mx) = f6d(:,j,k,l,m,n)
!                    buf_i(r_mn:r_mx) = r_halo(:,j,k,l,m,n)
!                    ! (2) perform interpolation
!                    call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                           buf_i, &
!                           buf_o, &
!                           displacement(l), &
!                           self%lagrange_width(1))
!                    ! (3) copy-back interpolated values
!                    f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx)
!                 end do
!              end do
!           end do
!        end do
!     end do
! !$omp end do
!     deallocate(buf_i)
!     deallocate(buf_o)
! !$omp end parallel

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
                  do w=0,wx-1
                    do i=l_mn,l_mx
                      buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=r_mn,r_mx
                      buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,w), &
                          buf_o(:,w), &
                          displacement(l), &
                          lw)
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      f6d(i,j+w,k,l,m,n) = buf_o(i,w)
                    end do
                  end do
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(:,j,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(:,j,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(:,j,k,l,m,n)
                   do i=l_mn,l_mx
                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do i=c_mn,c_mx
                     buf_i(i,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do i=r_mn,r_mx
                     buf_i(i,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(l), &
                          lw)
                   ! (3) copy-back interpolated values
!                   f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx,1)
                   do i=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(i,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1

 !> Advection along eta1 with displacement dependent on eta45
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_dispeta45(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(4):decomposition%local%mx(4), decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! NOTE: Buffering is done here for the second dim "j".
    !       Seems to be marginally beneficial.
    wx = get_wx(decomposition, 2)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
                  do w=0,wx-1
                    do i=l_mn,l_mx
                      buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=r_mn,r_mx
                      buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,w), &
                          buf_o(:,w), &
                          displacement(l,m), &
                          lw)
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      f6d(i,j+w,k,l,m,n) = buf_o(i,w)
                    end do
                  end do
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(:,j,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(:,j,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(:,j,k,l,m,n)
                   do i=l_mn,l_mx
                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do i=c_mn,c_mx
                     buf_i(i,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do i=r_mn,r_mx
                     buf_i(i,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(l,m), &
                          lw)
                   ! (3) copy-back interpolated values
!                   f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx,1)
                   do i=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(i,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_dispeta45


  !> Advection along eta2 with displacement dependent on eta5
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

! !$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
!     allocate(buf_i(l_mn:r_mx))
!     allocate(buf_o(l_mn:r_mx))
! !$omp do OMP_COLLAPSE OMP_SCHEDULE
!     do n=decomposition%local%mn(6), decomposition%local%mx(6)
!        do m=decomposition%local%mn(5), decomposition%local%mx(5)
!           do l=decomposition%local%mn(4), decomposition%local%mx(4)
!              do k=decomposition%local%mn(3), decomposition%local%mx(3)
!                 do i=decomposition%local%mn(1), decomposition%local%mx(1)
! !                   buf_i = f6d(i,:,k,l,m,n)
!                    ! (1) fill input buffer piecewise
!                    buf_i(l_mn:l_mx) = l_halo(i,:,k,l,m,n)
!                    buf_i(c_mn:c_mx) = f6d(i,:,k,l,m,n)
!                    buf_i(r_mn:r_mx) = r_halo(i,:,k,l,m,n)
!                    ! (2) perform interpolation
!                    call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                           buf_i, &
!                           buf_o, &
!                           displacement(m), &
!                           self%lagrange_width(1))
!                    ! (3) copy-back interpolated values
!                    f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx)
! !                   f6d(i,:,k,l,m,n) = buf_o
!                 end do
!              end do
!           end do
!        end do
!     end do
! !$omp end do
!     deallocate(buf_i)
!     deallocate(buf_o)
! !$omp end parallel

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
            do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
              do i=loop_mn(1), loop_mx(1), wx
                do j=l_mn,l_mx
                  do w=0,wx-1
                    buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                  end do
                end do
                do j=r_mn,r_mx
                  do w=0,wx-1
                    buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do w=0,wx-1
                  call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                         buf_i(:,w), &
                         buf_o(:,w), &
                         displacement(m), &
                         lw)
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    f6d(i+w,j,k,l,m,n) = buf_o(j,w)
                  end do
                end do
              end do
#else
              do i=loop_mn(1), loop_mx(1)
                 ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(i,:,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(i,:,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(i,:,k,l,m,n)
                 do j=l_mn,l_mx
                   buf_i(j,0) = l_halo(i,j,k,l,m,n)
                 enddo
                 do j=c_mn,c_mx
                   buf_i(j,0) = f6d(i,j,k,l,m,n)
                 enddo
                 do j=r_mn,r_mx
                   buf_i(j,0) = r_halo(i,j,k,l,m,n)
                 enddo
                 ! (2) perform interpolation
                 call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                        buf_i(:,0), &
                        buf_o(:,0), &
                        displacement(m), &
                        lw)
!                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx,1)
                 do j=c_mn,c_mx
                   f6d(i,j,k,l,m,n) = buf_o(j,0)
                 enddo
              end do
#endif
            end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2


!> Advection along eta2 with displacement dependent on eta45
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_dispeta45(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(&
         decomposition%local%mn(4):decomposition%local%mx(4),&
         decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
            do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
              do i=loop_mn(1), loop_mx(1), wx
                do j=l_mn,l_mx
                  do w=0,wx-1
                    buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                  end do
                end do
                do j=r_mn,r_mx
                  do w=0,wx-1
                    buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do w=0,wx-1
                  call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                         buf_i(:,w), &
                         buf_o(:,w), &
                         displacement(l,m), &
                         lw)
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    f6d(i+w,j,k,l,m,n) = buf_o(j,w)
                  end do
                end do
              end do
#else
              do i=loop_mn(1), loop_mx(1)
                 ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(i,:,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(i,:,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(i,:,k,l,m,n)
                 do j=l_mn,l_mx
                   buf_i(j,0) = l_halo(i,j,k,l,m,n)
                 enddo
                 do j=c_mn,c_mx
                   buf_i(j,0) = f6d(i,j,k,l,m,n)
                 enddo
                 do j=r_mn,r_mx
                   buf_i(j,0) = r_halo(i,j,k,l,m,n)
                 enddo
                 ! (2) perform interpolation
                 call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                        buf_i(:,0), &
                        buf_o(:,0), &
                        displacement(l,m), &
                        lw)
!                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx,1)
                 do j=c_mn,c_mx
                   f6d(i,j,k,l,m,n) = buf_o(j,0)
                 enddo
              end do
#endif
            end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_dispeta45





  !> Advection along eta3 with displacement dependent on eta6
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(6):decomposition%local%mx(6)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 3  ! eta3
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do j=loop_mn(2), loop_mx(2)
!                do i=loop_mn(1), loop_mx(1)
!                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx) = l_halo(i,j,:,l,m,n)
!                   buf_i(c_mn:c_mx) = f6d(i,j,:,l,m,n)
!                   buf_i(r_mn:r_mx) = r_halo(i,j,:,l,m,n)
!                   ! (2) perform interpolation
!                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                          buf_i, &
!                          buf_o, &
!                          displacement(n), &
!                          self%lagrange_width(1))
!                   ! (3) copy-back interpolated values
!                   f6d(i,j,:,l,m,n) = buf_o(c_mn:c_mx)
!                end do
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                  do k=l_mn,l_mx
                    do w=0,wx-1
                      buf_i(k,w) = l_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do k=c_mn,c_mx
                    do w=0,wx-1
                      buf_i(k,w) = f6d(i+w,j,k,l,m,n)
                    end do
                  end do
                  do k=r_mn,r_mx
                    do w=0,wx-1
                      buf_i(k,w) = r_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,w), &
                          buf_o(:,w), &
                          displacement(n), &
                          lw)
                  end do
                  do k=c_mn,c_mx
                    do w=0,wx-1
                      f6d(i+w,j,k,l,m,n) = buf_o(k,w)
                    end do
                  end do
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(i,j,:,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(i,j,:,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(i,j,:,l,m,n)
                   do k=l_mn,l_mx
                     buf_i(k,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do k=c_mn,c_mx
                     buf_i(k,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do k=r_mn,r_mx
                     buf_i(k,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(n), &
                          lw)
!                   f6d(i,j,:,l,m,n) = buf_o(c_mn:c_mx,1)
                   do k=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(k,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3


  !> Advection along eta4 with displacement dependent on eta1-3
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta4(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
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
    sll_int32, parameter :: id = 4  ! eta4
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(2)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(2)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
      do m=loop_mn(5), loop_mx(5)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
!                do i=loop_mn(1), loop_mx(1)
!                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx) = l_halo(i,j,k,:,m,n)
!                   buf_i(c_mn:c_mx) = f6d(i,j,k,:,m,n)
!                   buf_i(r_mn:r_mx) = r_halo(i,j,k,:,m,n)
!                   ! (2) perform interpolation
!                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                          buf_i, &
!                          buf_o, &
!                          displacement(i,j,k), &
!                          self%lagrange_width(2))
!                   ! (3) copy-back interpolated values
!                   f6d(i,j,k,:,m,n) = buf_o(c_mn:c_mx)
!                end do
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                  do l=l_mn,l_mx
                    do w=0,wx-1
                      buf_i(l,w) = l_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do l=c_mn,c_mx
                    do w=0,wx-1
                      buf_i(l,w) = f6d(i+w,j,k,l,m,n)
                    end do
                  end do
                  do l=r_mn,r_mx
                    do w=0,wx-1
                      buf_i(l,w) = r_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                           buf_i(:,w), &
                           buf_o(:,w), &
                           displacement(i+w,j,k), &
                           lw)
                  end do
                  do l=c_mn,c_mx
                    do w=0,wx-1
                      f6d(i+w,j,k,l,m,n) = buf_o(l,w)
                    end do
                  end do
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx) = l_halo(i,j,k,:,m,n)
!                   buf_i(c_mn:c_mx) = f6d(i,j,k,:,m,n)
!                   buf_i(r_mn:r_mx) = r_halo(i,j,k,:,m,n)
                   do l=l_mn,l_mx
                     buf_i(l,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do l=c_mn,c_mx
                     buf_i(l,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do l=r_mn,r_mx
                     buf_i(l,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(i,j,k), &
                          lw)
                   do l=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(l,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta4


  !> Advection along eta5 with displacement dependent on eta1-3
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta5(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
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
    sll_int32, parameter :: id = 5  ! eta5
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(2)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(2)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)
!                do i=loop_mn(1), loop_mx(1)
!                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx) = l_halo(i,j,k,l,:,n)
!                   buf_i(c_mn:c_mx) = f6d(i,j,k,l,:,n)
!                   buf_i(r_mn:r_mx) = r_halo(i,j,k,l,:,n)
!                   ! (2) perform interpolation
!                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                          buf_i, &
!                          buf_o, &
!                          displacement(i,j,k), &
!                          self%lagrange_width(2))
!                   ! (3) copy-back interpolated values
!                   f6d(i,j,k,l,:,n) = buf_o(c_mn:c_mx)
!                end do
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                  do m=l_mn,l_mx
                    do w=0,wx-1
                      buf_i(m,w) = l_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do m=c_mn,c_mx
                    do w=0,wx-1
                      buf_i(m,w) = f6d(i+w,j,k,l,m,n)
                    end do
                  end do
                  do m=r_mn,r_mx
                    do w=0,wx-1
                      buf_i(m,w) = r_halo(i+w,j,k,l,m,n)
                    end do
                 end do
                 
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,w), &
                          buf_o(:,w), &
                          displacement(i+w,j,k), &
                          lw)
                  end do
                  do m=c_mn,c_mx
                    do w=0,wx-1
                      f6d(i+w,j,k,l,m,n) = buf_o(m,w)
                    end do
                  end do
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
                   do m=l_mn,l_mx
                     buf_i(m,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do m=c_mn,c_mx
                     buf_i(m,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do m=r_mn,r_mx
                     buf_i(m,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(i,j,k), &
                          lw)
                   do m=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(m,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta5


  !> Advection along eta6 with displacement dependent on eta1-3
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta6(self, decomposition, displacement, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
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
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(2)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(2)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do m=loop_mn(5), loop_mx(5)
       do l=loop_mn(4), loop_mx(4)
          do k=loop_mn(3), loop_mx(3)
             do j=loop_mn(2), loop_mx(2)

!                do i=loop_mn(1), loop_mx(1)
!                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx) = l_halo(i,j,k,l,m,:)
!                   buf_i(c_mn:c_mx) = f6d(i,j,k,l,m,:)
!                   buf_i(r_mn:r_mx) = r_halo(i,j,k,l,m,:)
!                   ! (2) perform interpolation
!                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
!                          buf_i, &
!                          buf_o, &
!                          displacement(i,j,k), &
!                          self%lagrange_width(2))
!                   ! (3) copy-back interpolated values
!                   f6d(i,j,k,l,m,:) = buf_o(c_mn:c_mx)
!                end do
#ifndef DISABLE_CACHE_BLOCKING
                do i=loop_mn(1), loop_mx(1), wx
                  do n=l_mn,l_mx
                    do w=0,wx-1
                      buf_i(n,w) = l_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do n=c_mn,c_mx
                    do w=0,wx-1
                      buf_i(n,w) = f6d(i+w,j,k,l,m,n)
                    end do
                  end do
                  do n=r_mn,r_mx
                    do w=0,wx-1
                      buf_i(n,w) = r_halo(i+w,j,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,w), &
                          buf_o(:,w), &
                          displacement(i+w,j,k), &
                          lw)
                  end do
                  do n=c_mn,c_mx
                    do w=0,wx-1
                      f6d(i+w,j,k,l,m,n) = buf_o(n,w)
                    end do
                  end do
                end do
#else
                do i=loop_mn(1), loop_mx(1)
                   ! (1) fill input buffer piecewise
                   do n=l_mn,l_mx
                     buf_i(n,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do n=c_mn,c_mx
                     buf_i(n,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do n=r_mn,r_mx
                     buf_i(n,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells( &
                          buf_i(:,0), &
                          buf_o(:,0), &
                          displacement(i,j,k), &
                          lw)
!                   f6d(i,j,k,:,m,n) = buf_o(c_mn:c_mx)
                   do n=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(n,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta6



 !> Advection along eta1 with displacement dependent on eta4-5
  subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta1_d45(self, decomposition, displacement, f6d)

    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(4):decomposition%local%mx(4),decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    ! ---
    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx, n_omp_threads, deg
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx, l_mn_l, r_mx_l
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: si
    sll_real64 :: alpha

    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)  !decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)  !decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)  !decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)  !decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! NOTE: Buffering is done here for the second dim "j".
    !       Seems to be marginally beneficial.
    wx = get_wx(decomposition, 2)
    deg = self%lagrange_width(1)
#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,si,alpha,l_mn_l,r_mx_l)
!    allocate(buf_i(l_mn:r_mx,0:wx-1))
!    allocate(buf_o(l_mn:r_mx,0:wx-1))
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)

             si = floor(displacement(l,m))
             alpha = displacement(l,m) - real(si,f64)
             l_mn_l = -deg/2+1+si+c_mn
             r_mx_l = c_mx+deg/2+si

             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
                  do w=0,wx-1
                    do i=l_mn_l,l_mx
                      buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=r_mn,r_mx_l
                      buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                    end do
                 end do
                  do w=0,wx-1
                   call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
                          buf_i(l_mn_l:r_mx_l,w), &
                          buf_o(l_mn_l:r_mx_l,w), &
                          alpha, &
                          deg)
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      f6d(i,j+w,k,l,m,n) = buf_o(i+si,w)
                    end do
                  end do
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   ! (1) fill input buffer piecewise
                   do i=l_mn_l,l_mx
                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do i=c_mn,c_mx
                     buf_i(i,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do i=r_mn,r_mx_l
                     buf_i(i,0) = r_halo(i,j,k,l,m,n)
                   enddo
                   ! (2) perform interpolation
                   call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
                          buf_i(l_mn_l:r_mx_l,0), &
                          buf_o(l_mn_l:r_mx_l,0), &
                          alpha, &
                          deg)
                   ! (3) copy-back interpolated values
                   do i=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(i+si,0)
                   enddo
                end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta1_d45


  !> Advection along eta2 with displacement dependent on eta4-5
  subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta2_d45(self, decomposition, displacement, f6d)

    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    sll_real64, intent(in) :: displacement(decomposition%local%mn(4):decomposition%local%mx(4),decomposition%local%mn(5):decomposition%local%mx(5)) !< displacement vector
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    ! ---
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w,wx, n_omp_threads, deg
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx, l_mn_l, r_mx_l
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2)
#endif
    sll_int32 :: si
    sll_real64 :: alpha

    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)  !decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)  !decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)  !decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)  !decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition, 1)
    deg = self%lagrange_width(1)
#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,si,alpha,l_mn_l,r_mx_l)
!    allocate(buf_i(l_mn:r_mx,0:wx-1))
!    allocate(buf_o(l_mn:r_mx,0:wx-1))
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mn, mx)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(l_mn:r_mx, 0:wx-1))
#endif
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)

             si = floor(displacement(l,m))
             alpha = displacement(l,m) - real(si,f64)
             l_mn_l = -deg/2+1+si+c_mn
             r_mx_l = c_mx+deg/2+si

            do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
              do i=loop_mn(1), loop_mx(1), wx
                do j=l_mn_l,l_mx
                  do w=0,wx-1
                    buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                  end do
                end do
                do j=r_mn,r_mx_l
                  do w=0,wx-1
                    buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do w=0,wx-1
                  call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
                         buf_i(l_mn_l:r_mx_l,w), &
                         buf_o(l_mn_l:r_mx_l,w), &
                         alpha, &
                         deg)
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    f6d(i+w,j,k,l,m,n) = buf_o(j+si,w)
                  end do
                end do
              end do
#else
              do i=loop_mn(1), loop_mx(1)
                 do j=l_mn_l,l_mx
                   buf_i(j,0) = l_halo(i,j,k,l,m,n)
                 enddo
                 do j=c_mn,c_mx
                   buf_i(j,0) = f6d(i,j,k,l,m,n)
                 enddo
                 do j=r_mn,r_mx_l
                   buf_i(j,0) = r_halo(i,j,k,l,m,n)
                 enddo
                 call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
                        buf_i(l_mn_l:r_mx_l,0), &
                        buf_o(l_mn_l:r_mx_l,0), &
                        alpha, &
                        deg)
                 do j=c_mn,c_mx
                   f6d(i,j,k,l,m,n) = buf_o(j+si,0)
                 enddo
              end do
#endif
            end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
!$omp end parallel
  end subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta2_d45

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!  From here the block of the advection routines for x with x dependent velocities
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Advection along eta1 with displacement dependent on eta45 (fixed interval Lagrange interpolation )
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_givenv(self, decomposition, rdelta_x, velocity_grid, velocity_shift, rotation_matrix, propagation_matrix, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    type(sll_t_array), intent(in) :: velocity_grid(6)
    sll_real64, intent(in) :: velocity_shift(decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3),3)
    sll_real64, intent(in) :: rotation_matrix(3,3)
    sll_real64, intent(in) :: propagation_matrix(3,3)
    sll_real64, intent(in) :: rdelta_x !< 1/delta_x(1)
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 1  ! eta1
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
    sll_real64, pointer :: displacement(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2), mno(2), mxo(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_real64 :: v_coord(3)

    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    ! NOTE: Buffering is done here for the second dim "j".
    !       Seems to be marginally beneficial.
    wx = get_wx(decomposition, 2)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
    mno = [c_mn, 0]
    mxo = [c_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,displacement, v_coord)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mno, mxo)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(c_mn:c_mx, 0:wx-1))
#endif
    allocate(displacement(c_mn:c_mx, 0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
             do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
                do j=loop_mn(2), loop_mx(2), wx
                  do w=0,wx-1
                    do i=l_mn,l_mx
                      buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      buf_i(i,w) = f6d(i,j+w,k,l,m,n)
                    end do
                  end do
                  do w=0,wx-1
                    do i=r_mn,r_mx
                      buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
                    end do
                 end do
                 

                 ! Compute the displacement vector
                 
                 do w=0,wx-1
                    do i=c_mn, c_mx

                       v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j+w,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j+w,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j+w,k,3)

                       displacement(i,w) = rdelta_x * ( &
                            propagation_matrix(1,1) * v_coord(1)+&
                            propagation_matrix(1,2) * v_coord(2)+&
                            propagation_matrix(1,3) * v_coord(3) )
                       if (abs(displacement(i,w))> 1.0_f64 )  then
                          print*, i,j+w,k,l,m,n, displacement(i,w)
                          stop
                       end if
                       !print*, v_coord
                    end do
                 end do
                 !print*, displacement
                 
                 do w=0,wx-1                  
                    ! (2) perform interpolation
                    call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                         buf_i(:,w), &
                         buf_o(:,w), &
                         displacement(:,w), &
                         lw, &
                         0)
                  end do
                  do w=0,wx-1
                    do i=c_mn,c_mx
                      f6d(i,j+w,k,l,m,n) = buf_o(i,w)
                    end do
                  end do
                end do
#else
                do j=loop_mn(2), loop_mx(2)
                   
                   ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(:,j,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(:,j,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(:,j,k,l,m,n)
                   do i=l_mn,l_mx
                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
                   enddo
                   do i=c_mn,c_mx
                     buf_i(i,0) = f6d(i,j,k,l,m,n)
                   enddo
                   do i=r_mn,r_mx
                     buf_i(i,0) = r_halo(i,j,k,l,m,n)
                  enddo

                  ! Compute the displacement vector
                  do i=c_mn, c_mx
                       v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,3)

                       displacement(i,0) = rdelta_x * ( &
                            propagation_matrix(1,1) * v_coord(1)+&
                            propagation_matrix(1,2) * v_coord(2)+&
                            propagation_matrix(1,3) * v_coord(3) )
                  end do
  
                  ! (2) perform interpolation
                  call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                       buf_i(:,0), &
                       buf_o(:,0), &
                       displacement, &
                       lw, &
                       0)
               
                  ! (3) copy-back interpolated values
                  !                   f6d(:,j,k,l,m,n) = buf_o(c_mn:c_mx,1)
                  do i=c_mn,c_mx
                     f6d(i,j,k,l,m,n) = buf_o(i,0)
                  enddo
               end do
#endif
             end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
    deallocate(displacement)
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta1_givenv

!> Advection along eta2 with displacement dependent on eta45
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_givenv(self, decomposition, rdelta_x, velocity_grid, velocity_shift, rotation_matrix, propagation_matrix, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    type(sll_t_array), intent(in) :: velocity_grid(6)
    sll_real64, intent(in) :: velocity_shift(decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3),3)
    sll_real64, intent(in) :: rotation_matrix(3,3)
    sll_real64, intent(in) :: propagation_matrix(3,3)
    sll_real64, intent(in) :: rdelta_x !< 1/delta_x(1)
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    sll_int32, parameter :: id = 2  ! eta2
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
    sll_real64, pointer :: displacement(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2), mno(2), mxo(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_real64 :: v_coord(3)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
    mno = [c_mn, 0]
    mxo = [c_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,displacement, v_coord)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mno, mxo)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(c_mn:c_mx, 0:wx-1))
#endif
    allocate(displacement(c_mn:c_mx, 0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
            do k=loop_mn(3), loop_mx(3)
#ifndef DISABLE_CACHE_BLOCKING
              do i=loop_mn(1), loop_mx(1), wx
                do j=l_mn,l_mx
                  do w=0,wx-1
                    buf_i(j,w) = l_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do j=c_mn,c_mx
                  do w=0,wx-1
                    buf_i(j,w) = f6d(i+w,j,k,l,m,n)
                  end do
                end do
                do j=r_mn,r_mx
                  do w=0,wx-1
                    buf_i(j,w) = r_halo(i+w,j,k,l,m,n)
                  end do
               end do
               ! Compute the displacement vector
               do j=c_mn, c_mx
                  do w=0,wx-1
                     v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,3)
                     
                     displacement(j,w) = rdelta_x * ( &
                          propagation_matrix(2,1) * v_coord(1)+&
                          propagation_matrix(2,2) * v_coord(2)+&
                          propagation_matrix(2,3) * v_coord(3) )
                     if (abs(displacement(j,w))> 1.0_f64 )  then
                        print*, i+w,j,k,l,m,n, displacement(j,w)
                        stop
                     end if
                     !displacement(j,w) = 0.0_f64!- delta_tx * velocity_grid(i+w,j,k,2)
                  end do
               end do
                do w=0,wx-1
                    ! (2) perform interpolation
                    call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                         buf_i(:,w), &
                         buf_o(:,w), &
                         displacement(:,w), &
                         lw, &
                         0)
                end do
                do j=c_mn,c_mx
                   do w=0,wx-1
                    f6d(i+w,j,k,l,m,n) = buf_o(j,w)
                  end do
                end do
             end do
#else
              do i=loop_mn(1), loop_mx(1)
                 ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(i,:,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(i,:,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(i,:,k,l,m,n)
                 do j=l_mn,l_mx
                   buf_i(j,0) = l_halo(i,j,k,l,m,n)
                 enddo
                 do j=c_mn,c_mx
                   buf_i(j,0) = f6d(i,j,k,l,m,n)
                 enddo
                 do j=r_mn,r_mx
                   buf_i(j,0) = r_halo(i,j,k,l,m,n)
                enddo
                ! Compute the displacement vector
                do j=c_mn, c_mx
                   v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,3)
                     
                     displacement(j,0) = rdelta_x * ( &
                          propagation_matrix(2,1) * v_coord(1)+&
                          propagation_matrix(2,2) * v_coord(2)+&
                          propagation_matrix(2,3) * v_coord(3) )
                    !displacement(j,0) = - delta_tx * velocity_grid(i,j,k,2)
                 end do
                 ! (2) perform interpolation
                 call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                      buf_i(:,0), &
                      buf_o(:,0), &
                      displacement(:,0), &
                      lw, &
                      0)
                 !                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx,1)
                 do j=c_mn,c_mx
                   f6d(i,j,k,l,m,n) = buf_o(j,0)
                 enddo
              end do
#endif
            end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
    deallocate(displacement)
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta2_givenv
  
  !> Advection along eta3 with displacement dependent on eta45
  subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3_givenv(self, decomposition, rdelta_x, velocity_grid, velocity_shift, rotation_matrix, propagation_matrix, f6d)
    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
    type(sll_t_array), intent(in) :: velocity_grid(6)
    sll_real64, intent(in) :: velocity_shift(decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3),3)
    sll_real64, intent(in) :: rotation_matrix(3,3)
    sll_real64, intent(in) :: propagation_matrix(3,3)
    sll_real64, intent(in) :: rdelta_x !< 1/delta_x(1)
    sll_real64, intent(inout) :: f6d(&
         decomposition%local%mn(1):decomposition%local%mx(1), &
         decomposition%local%mn(2):decomposition%local%mx(2), &
         decomposition%local%mn(3):decomposition%local%mx(3), &
         decomposition%local%mn(4):decomposition%local%mx(4), &
         decomposition%local%mn(5):decomposition%local%mx(5), &
         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
    
    sll_int32, parameter :: id = 3  ! eta3
    sll_int32 :: i,j,k,l,m,n,w,wx, lw
    sll_real64, pointer :: buf_i(:,:)
    sll_real64, pointer :: buf_o(:,:)
    sll_real64, pointer :: displacement(:,:)
#ifdef USE_FMEMPOOL
    sll_int32 :: mn(2), mx(2), mno(2), mxo(2)
#endif
    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx
    sll_int32, pointer :: loop_mn(:), loop_mx(:)
    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
    sll_real64 :: v_coord(3)

    ! NOTE: Current restriction for the memory-saving `slim` halos,
    !       to be lowered later for the real one-sided halo implementation.
    !       One-sided halos require a more clever index handling (though straight-forward).
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == decomposition%local%halo_right%nw(id))
    SLL_ASSERT_ALWAYS(decomposition%local%halo_left%nw(id) == (self%lagrange_width(1)-1)/2)
    ! cache values to avoid multiple dereferencing inside the loops
    l_mn = decomposition%local%halo_left%mn(id)
    l_mx = decomposition%local%halo_left%mx(id)
    c_mn = decomposition%local%mn(id)
    c_mx = decomposition%local%mx(id)
    r_mn = decomposition%local%halo_right%mn(id)
    r_mx = decomposition%local%halo_right%mx(id)
    l_halo => decomposition%local%halo_left%buf
    r_halo => decomposition%local%halo_right%buf
    loop_mn => decomposition%local%mn
    loop_mx => decomposition%local%mx
    wx = get_wx(decomposition)
    lw = self%lagrange_width(1)

#ifdef USE_FMEMPOOL
    mn = [l_mn, 0]
    mx = [r_mx, wx-1]
    mno = [c_mn, 0]
    mxo = [c_mx, wx-1]
#endif

!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,displacement,v_coord)
#ifdef USE_FMEMPOOL
    call mp_acquire(buf_i, mn, mx)
    call mp_acquire(buf_o, mno, mxo)
#else
    allocate(buf_i(l_mn:r_mx, 0:wx-1))
    allocate(buf_o(c_mn:c_mx, 0:wx-1))
#endif
    allocate(displacement(c_mn:c_mx, 0:wx-1))
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=loop_mn(6), loop_mx(6)
       do m=loop_mn(5), loop_mx(5)
          do l=loop_mn(4), loop_mx(4)
            do j=loop_mn(2), loop_mx(2)
#ifndef DISABLE_CACHE_BLOCKING
              do i=loop_mn(1), loop_mx(1), wx
                do k=l_mn,l_mx
                  do w=0,wx-1
                    buf_i(k,w) = l_halo(i+w,j,k,l,m,n)
                  end do
                end do
                do k=c_mn,c_mx
                  do w=0,wx-1
                    buf_i(k,w) = f6d(i+w,j,k,l,m,n)
                  end do
                end do
                do k=r_mn,r_mx
                  do w=0,wx-1
                    buf_i(k,w) = r_halo(i+w,j,k,l,m,n)
                  end do
               end do
               ! Compute the displacement vector
               do k=c_mn, c_mx
                  do w=0,wx-1
                     v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i+w,j,k,3)
                     
                     displacement(k,w) = rdelta_x * ( &
                          propagation_matrix(3,1) * v_coord(1)+&
                          propagation_matrix(3,2) * v_coord(2)+&
                          propagation_matrix(3,3) * v_coord(3) )
                  end do
               end do
                do w=0,wx-1
                    ! (2) perform interpolation
                    call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                         buf_i(:,w), &
                         buf_o(:,w), &
                         displacement(:,w), &
                         lw, &
                         0)
                end do
                do k=c_mn,c_mx
                  do w=0,wx-1
                    f6d(i+w,j,k,l,m,n) = buf_o(k,w)
                  end do
                end do
              end do
#else
              do i=loop_mn(1), loop_mx(1)
                 ! (1) fill input buffer piecewise
!                   buf_i(l_mn:l_mx,1) = l_halo(i,:,k,l,m,n)
!                   buf_i(c_mn:c_mx,1) = f6d(i,:,k,l,m,n)
!                   buf_i(r_mn:r_mx,1) = r_halo(i,:,k,l,m,n)
                 do k=l_mn,l_mx
                   buf_i(k,0) = l_halo(i,j,k,l,m,n)
                 enddo
                 do k=c_mn,c_mx
                   buf_i(k,0) = f6d(i,j,k,l,m,n)
                 enddo
                 do k=r_mn,r_mx
                   buf_i(k,0) = r_halo(i,j,k,l,m,n)
                 enddo               ! Compute the displacement vector
                 do k=c_mn, c_mx
                    !displacement(k,0) = - delta_tx * velocity_grid(i,j,k,3)
                    v_coord(1) = rotation_matrix(1,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(1,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(1,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,1)
                     v_coord(2) = rotation_matrix(2,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(2,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(2,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,2)
                     v_coord(3) = rotation_matrix(3,1) * velocity_grid(4)%vals(l-loop_mn(4)+1) + &
                          rotation_matrix(3,2) * velocity_grid(5)%vals(m-loop_mn(5)+1) + &
                          rotation_matrix(3,3) * velocity_grid(6)%vals(n-loop_mn(6)+1) + &
                          velocity_shift(i,j,k,3)
                     
                     displacement(k,0) = rdelta_x * ( &
                          propagation_matrix(3,1) * v_coord(1)+&
                          propagation_matrix(3,2) * v_coord(2)+&
                          propagation_matrix(3,3) * v_coord(3) )
                 end do
                 ! (2) perform interpolation
                 call sll_s_lagrange_interpolation_1d_fast_haloc_cells( &
                      buf_i(:,0), &
                      buf_o(:,0), &
                      displacement(:,0), &
                      lw, &
                      0)
                 !                   f6d(i,:,k,l,m,n) = buf_o(c_mn:c_mx,1)
                 do k=c_mn,c_mx
                   f6d(i,j,k,l,m,n) = buf_o(k,0)
                 enddo
              end do
#endif
            end do
          end do
       end do
    end do
!$omp end do
#ifdef USE_FMEMPOOL
    call mp_release(buf_i)
    call mp_release(buf_o)
#else
    deallocate(buf_i)
    deallocate(buf_o)
#endif
    deallocate(displacement)
!$omp end parallel
  end subroutine sll_s_advection_6d_lagrange_dd_slim_advect_eta3_givenv
  

!!$  !> Advection along eta1 with displacement with centered lagrange interpolation
!!$  subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta1_givenv(self, decomposition, velocity_grid, delta_t, f6d)
!!$
!!$    type(sll_t_advection_6d_lagrange_dd_slim), intent(inout) :: self !< advector object
!!$    type(sll_t_decomposition_slim_6d), target, intent(inout) :: decomposition
!!$    sll_real64, intent(in) :: velocity_grid(decomposition%local%mn(1):decomposition%local%mx(1),decomposition%local%mn(2):decomposition%local%mx(2),decomposition%local%mn(3):decomposition%local%mx(3),3) !< displacement vector
!!$    sll_real64, intent(in) :: delta_t
!!$    sll_real64, intent(inout) :: f6d(&
!!$         decomposition%local%mn(1):decomposition%local%mx(1), &
!!$         decomposition%local%mn(2):decomposition%local%mx(2), &
!!$         decomposition%local%mn(3):decomposition%local%mx(3), &
!!$         decomposition%local%mn(4):decomposition%local%mx(4), &
!!$         decomposition%local%mn(5):decomposition%local%mx(5), &
!!$         decomposition%local%mn(6):decomposition%local%mx(6)) !< value of the function on input and advected function values on output
!!$    ! ---
!!$    sll_int32, parameter :: id = 1  ! eta1
!!$    sll_int32 :: i,j,k,l,m,n,w,wx, n_omp_threads, deg
!!$    sll_real64, pointer :: buf_i(:,:)
!!$    sll_real64, pointer :: buf_o(:,:)
!!$    sll_int32 :: l_mn, l_mx, c_mn, c_mx, r_mn, r_mx, l_mn_l, r_mx_l
!!$    sll_int32, pointer :: loop_mn(:), loop_mx(:)
!!$    HALO_DTYPE, pointer :: l_halo(:,:,:,:,:,:)
!!$    HALO_DTYPE, pointer :: r_halo(:,:,:,:,:,:)
!!$#ifdef USE_FMEMPOOL
!!$    sll_int32 :: mn(2), mx(2)
!!$#endif
!!$    sll_int32 :: si
!!$    sll_real64 :: alpha
!!$    sll_real64 :: displacement
!!$
!!$    ! cache values to avoid multiple dereferencing inside the loops
!!$    l_mn = decomposition%local%halo_left%mn(id)  !decomposition%local%halo_left%mn(id)
!!$    l_mx = decomposition%local%halo_left%mx(id)  !decomposition%local%halo_left%mx(id)
!!$    c_mn = decomposition%local%mn(id)
!!$    c_mx = decomposition%local%mx(id)
!!$    r_mn = decomposition%local%halo_right%mn(id)  !decomposition%local%halo_right%mn(id)
!!$    r_mx = decomposition%local%halo_right%mx(id)  !decomposition%local%halo_right%mx(id)
!!$    l_halo => decomposition%local%halo_left%buf
!!$    r_halo => decomposition%local%halo_right%buf
!!$    loop_mn => decomposition%local%mn
!!$    loop_mx => decomposition%local%mx
!!$    ! NOTE: Buffering is done here for the second dim "j".
!!$    !       Seems to be marginally beneficial.
!!$    wx = get_wx(decomposition, 2)
!!$    deg = self%lagrange_width(1)
!!$#ifdef USE_FMEMPOOL
!!$    mn = [l_mn, 0]
!!$    mx = [r_mx, wx-1]
!!$#endif
!!$
!!$!$omp parallel default(shared) private(i,j,k,l,m,n,w,buf_i,buf_o,si,alpha,l_mn_l,r_mx_l)
!!$!    allocate(buf_i(l_mn:r_mx,0:wx-1))
!!$!    allocate(buf_o(l_mn:r_mx,0:wx-1))
!!$#ifdef USE_FMEMPOOL
!!$    call mp_acquire(buf_i, mn, mx)
!!$    call mp_acquire(buf_o, mn, mx)
!!$#else
!!$    allocate(buf_i(l_mn:r_mx, 0:wx-1))
!!$    allocate(buf_o(l_mn:r_mx, 0:wx-1))
!!$#endif
!!$!$omp do OMP_COLLAPSE OMP_SCHEDULE
!!$    do n=loop_mn(6), loop_mx(6)
!!$       do m=loop_mn(5), loop_mx(5)
!!$          do l=loop_mn(4), loop_mx(4)
!!$
!!$             do k=loop_mn(3), loop_mx(3)
!!$#ifndef DISABLE_CACHE_BLOCKING
!!$                do j=loop_mn(2), loop_mx(2), wx
!!$                  do w=0,wx-1
!!$                    do i=l_mn_l,l_mx
!!$                      buf_i(i,w) = l_halo(i,j+w,k,l,m,n)
!!$                    end do
!!$                  end do
!!$                  do w=0,wx-1
!!$                    do i=c_mn,c_mx
!!$                      buf_i(i,w) = f6d(i,j+w,k,l,m,n)
!!$                    end do
!!$                  end do
!!$                  do w=0,wx-1
!!$                    do i=r_mn,r_mx_l
!!$                      buf_i(i,w) = r_halo(i,j+w,k,l,m,n)
!!$                    end do
!!$                 end do
!!$                  do w=0,wx-1
!!$                   call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
!!$                          buf_i(l_mn_l:r_mx_l,w), &
!!$                          buf_o(l_mn_l:r_mx_l,w), &
!!$                          alpha, &
!!$                          deg)
!!$                  end do
!!$                  do w=0,wx-1
!!$                    do i=c_mn,c_mx
!!$                      f6d(i,j+w,k,l,m,n) = buf_o(i+si,w)
!!$                    end do
!!$                  end do
!!$                end do
!!$#else
!!$                do j=loop_mn(2), loop_mx(2)
!!$
!!$                   displacement = -delta_t * velocity_grid(c_mn:c_mx,j,k,1) 
!!$
!!$                   si = floor(displacement(l,m))
!!$                   alpha = displacement(l,m) - real(si,f64)
!!$                   l_mn_l = -deg/2+1+si+c_mn
!!$                   r_mx_l = c_mx+deg/2+si
!!$                   
!!$                   ! (1) fill input buffer piecewise
!!$                   do i=l_mn_l,l_mx
!!$                     buf_i(i,0) = l_halo(i,j,k,l,m,n)
!!$                   enddo
!!$                   do i=c_mn,c_mx
!!$                     buf_i(i,0) = f6d(i,j,k,l,m,n)
!!$                   enddo
!!$                   do i=r_mn,r_mx_l
!!$                     buf_i(i,0) = r_halo(i,j,k,l,m,n)
!!$                   enddo
!!$                   ! (2) perform interpolation
!!$                   call sll_s_lagrange_interpolation_1d_fast_disp_even_halo_cells( &
!!$                          buf_i(l_mn_l:r_mx_l,0), &
!!$                          buf_o(l_mn_l:r_mx_l,0), &
!!$                          alpha, &
!!$                          deg)
!!$                   ! (3) copy-back interpolated values
!!$                   do i=c_mn,c_mx
!!$                     f6d(i,j,k,l,m,n) = buf_o(i+si,0)
!!$                   enddo
!!$                end do
!!$#endif
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$!$omp end do
!!$#ifdef USE_FMEMPOOL
!!$    call mp_release(buf_i)
!!$    call mp_release(buf_o)
!!$#else
!!$    deallocate(buf_i)
!!$    deallocate(buf_o)
!!$#endif
!!$!$omp end parallel
!!$  end subroutine sll_s_advection_6d_clagrange_dd_slim_advect_eta1_d45
  
  

end module sll_m_advection_6d_lagrange_dd_slim
