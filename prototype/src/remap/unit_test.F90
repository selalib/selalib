program remap_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
 use sll_collective
#include "sll_remap.h"

  implicit none

  ! THIS IS A VERY BASIC TEST AND SOMETHING TRULY SUPERIOR AND ROBUST NEEDS
  ! TO BE DONE BEFORE THE REMAPPER IS DEPLOYED.


  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  integer, dimension(:,:,:), allocatable :: a3
  integer, dimension(:,:,:), allocatable :: b3
  sll_real64, dimension(:,:,:), allocatable :: ad3
  sll_real64, dimension(:,:,:), allocatable :: bd3
  ! Take a 3D array of dimensions 8X8X1
  integer, parameter                 :: total_sz_i = 8
  integer, parameter                 :: total_sz_j = 8
  integer, parameter                 :: total_sz_k = 1
  ! the process mesh
  integer, parameter                 :: pi = 4
  integer, parameter                 :: pj = 4
  integer, parameter                 :: pk = 1

  ! Split it in  16 processes, each with a local chunk 2X2X1
  integer                            :: local_sz_i 
  integer                            :: local_sz_j 
  integer                            :: local_sz_k 
  integer                            :: ierr
  integer                            :: myrank
  integer                            :: colsz        ! collective size
  integer                            :: i,j,k
  integer                            :: i_min, i_max
  integer                            :: j_min, j_max
  integer                            :: k_min, k_max
  integer                            :: node
  integer, dimension(1:3)            :: gcoords
  ! Remap stuff
  type(layout_3D_t), pointer         :: conf3_init
  type(layout_3D_t), pointer         :: conf3_final
  type(layout_3D_t), pointer         :: random_layout1

  type(remap_plan_3D_t), pointer     :: rmp3

!  integer, dimension(:), pointer :: limits

  print *, ' '
  print *, '--------------- REMAP test ---------------------'
  print *, ' '

  call flush()
  call sll_boot_collective()
!  SLL_ALLOCATE( a(local_sz), ierr )

  local_sz_i = total_sz_i/pi
  local_sz_j = total_sz_j/pj
  local_sz_k = total_sz_k/pk
  SLL_ALLOCATE( a3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  SLL_ALLOCATE( b3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  SLL_ALLOCATE( ad3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  SLL_ALLOCATE( bd3(1:local_sz_i,1:local_sz_j,1:local_sz_k), ierr )
  myrank    = sll_get_collective_rank(sll_world_collective)
  colsz     = sll_get_collective_size(sll_world_collective)

  conf3_init  => new_layout_3D( sll_world_collective )
  conf3_final => new_layout_3D( sll_world_collective )
  random_layout1 => new_layout_3D( sll_world_collective )

  do k=0, pk-1
     do j=0, pj-1
        do i=0, pi-1
           node = i+pi*(j+pj*k) ! linear index of node
           i_min = i*local_sz_i + 1
           i_max = i*local_sz_i + local_sz_i
           j_min = j*local_sz_j + 1
           j_max = j*local_sz_j + local_sz_j
           k_min = k*local_sz_k + 1
           k_max = k*local_sz_k + local_sz_k
           call set_layout_i_min( conf3_init, node, i_min )
           call set_layout_i_max( conf3_init, node, i_max )
           call set_layout_j_min( conf3_init, node, j_min )
           call set_layout_j_max( conf3_init, node, j_max )
           call set_layout_k_min( conf3_init, node, k_min )
           call set_layout_k_max( conf3_init, node, k_max )
        end do
     end do
  end do
      


!  call sll_view_lims_3D( conf3_init )

  ! Initialize the data. We use the information in the layout.
  do k=1, local_sz_k
     do j=1, local_sz_j
        do i=1, local_sz_i
           gcoords =  local_to_global_3D( conf3_init, (/i,j,k/) )
!           write (*,'(a,i4)') 'gcoords in rank: ', myrank
!           print *, gcoords(:)
!           call flush()
           a3(i,j,k) = gcoords(1) + &
                total_sz_i*((gcoords(2)-1) + total_sz_j*(gcoords(3)-1))
           ad3(i,j,k) = real(gcoords(1) + &
                total_sz_i*((gcoords(2)-1) + total_sz_j*(gcoords(3)-1)),f64)
        end do
     end do
  end do

  write (*,'(a,i4)') 'From rank: ', myrank
  print *, a3(:,:,:)
  print *, ad3(:,:,:)
  call flush()

  ! Initialize the final layout, in this case, just a transposition
  do k=0, pk-1
     do j=0, pj-1
        do i=0, pi-1
           node = i+pi*(j+pj*k) ! linear index of node
           i_min = i*local_sz_i + 1
           i_max = i*local_sz_i + local_sz_i
           j_min = j*local_sz_j + 1
           j_max = j*local_sz_j + local_sz_j
           k_min = k*local_sz_k + 1
           k_max = k*local_sz_k + local_sz_k
           call set_layout_i_min( conf3_final, node, j_min )
           call set_layout_i_max( conf3_final, node, j_max )
           call set_layout_j_min( conf3_final, node, i_min )
           call set_layout_j_max( conf3_final, node, i_max )
           call set_layout_k_min( conf3_final, node, k_min )
           call set_layout_k_max( conf3_final, node, k_max )
        end do
     end do
  end do
  
!  call sll_view_lims_3D( conf3_final )      
  

  rmp3 => NEW_REMAPPER_PLAN_3D( conf3_init, conf3_final, ad3)!INT32_SIZEOF(a3(1,1,1)) )
  call apply_remap_3D( rmp3, ad3, bd3 )
  print *, 'Remap operation completed.'
  write (*,'(a, i4)') 'the output data in rank: ', myrank
  print *, bd3(:,:,:)
  call flush()
  
  call delete_layout_3D( conf3_init )
  call delete_layout_3D( conf3_final )

  call sll_collective_barrier(sll_world_collective)
  call sll_halt_collective()
  print *, 'TEST COMPLETE'
  call flush()
!do i=1,30
!   call random_number(fpn)
!print *,gaussian_dev()
!end do
  !SLL_ALLOCATE(limits(pi), ierr)
#if 0
  call split_interval_randomly(1000,4,limits)
  print *, limits(:)
  call flush()
#endif


contains

  subroutine split_interval_randomly(n, num_halvings, ans)
    integer, intent(in) :: n
    integer, intent(in) :: num_halvings
    integer, dimension(:), pointer :: ans
    integer :: load_i = 1
    integer :: ierr
    SLL_ALLOCATE(ans(2*(2**num_halvings)),ierr)
    call split_interval_randomly_aux( 1, n, 0, num_halvings, load_i, ans )
  end subroutine split_interval_randomly

  recursive subroutine split_interval_randomly_aux( lo, hi, gen, lim, loadi, ans )
    intrinsic :: random_number
    integer, intent(in)                :: lo
    integer, intent(in)                :: hi
    integer, intent(in)                :: gen ! splitting generation
    integer, intent(in)                :: lim ! maximum number of splittings
    integer :: mid
!    real :: rand
    integer, parameter                 :: spread = 25
    ! dangerous... the following should be found in the enclosing scope...
    integer, intent(inout) :: loadi
    integer, dimension(:), pointer :: ans
!    integer, dimension(:), intent(out) :: ans ! answer
    if( (hi-lo).eq. 1 ) then
       ans(loadi) = lo
       ans(loadi+1) = hi
       loadi = loadi + 2
 !      write (*,'(i8,i8)') lo, hi
    else if (gen .eq. lim) then
       ans(loadi) = lo
       ans(loadi+1) = hi
       loadi = loadi + 2
 !      write (*,'(i8,i8)') lo, hi
    else
!       call random_number(rand)
       ! decided to use the gaussian because the uniform deviate can give
       ! a number very close to the 0 or one and such a small interval is
       ! hard to subdivide. In such case one may end up with less intervals
       ! and this is a problem since the number of processors is fixed from
       ! the beginning.
       mid = (hi+lo)/2 + gaussian_dev()*spread
!       mid = int(rand*(hi-lo))+lo
       call split_interval_randomly_aux( lo,   mid, gen+1, lim, loadi, ans )
       call split_interval_randomly_aux( mid+1, hi, gen+1, lim, loadi, ans )
    end if
  end subroutine split_interval_randomly_aux

  function gaussian_dev()
    intrinsic :: random_number
    real :: gaussian_dev
    real :: v1, v2
    real :: ran1, ran2
    real :: rsq
    real :: fac
    do
       call random_number(ran1)
       call random_number(ran2)
       v1 = 2.0*ran1-1.0
       v2 = 2.0*ran2-1.0
       rsq = v1*v1 + v2*v2
       if( (rsq .lt. 1.0) .and. (rsq .gt. 0.0) ) exit
    end do
    fac = sqrt(-2.0*log(rsq)/rsq)
    gaussian_dev = v1*fac
  end function gaussian_dev

end program remap_test
