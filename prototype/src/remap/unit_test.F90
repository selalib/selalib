program remap_test
  use sll_collective
#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
  implicit none
  
  ! THIS IS A VERY BASIC TEST AND SOMETHING TRULY SUPERIOR AND ROBUST NEEDS
  ! TO BE DONE BEFORE THE REMAPPER IS DEPLOYED.
  
  
  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array1, local_array2
  ! Take a 3D array of dimensions global_sz_i*global_sz_j*globall_sz_k
  integer , parameter                       :: global_sz_i = 16
  integer , parameter                       :: global_sz_j = 16
  integer , parameter                       :: global_sz_k = 16
  ! Local sizes
  integer                                   :: local_sz_i
  integer                                   :: local_sz_j
  integer                                   :: local_sz_k
  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  ! Split it in  16 processes, each with a local chunk 
  ! local_sz_i*local_sz_j*local_sz_k
  integer                                   :: ierr
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size
  ! Remap stuff
  type(layout_3D_t), pointer                :: layout1
  type(layout_3D_t), pointer                :: layout2
  
  type(remap_plan_3D_t), pointer            :: rmp3
  
  !  integer, dimension(:), pointer :: limits
  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 10
  integer                                   :: i_test
  integer                                   :: i, j, k
  sll_int32, dimension(1:3)                 :: global_index
  sll_real32   , dimension(1)               :: sum4test
  integer                                   :: ok

  call flush()
  call sll_boot_collective()

  
  colsz = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- REMAP test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
  end if

  
!  if (is_power_of_two(colsz)) then     
!  do, i_test=1, nbtest
  layout1  => new_layout_3D( sll_world_collective )        
  call two_power_rand_factorization(colsz, npi, npj, npk)
  
  if( myrank .eq. 0) then
     print *, 'for the initial configuration, '
     print *, 'a total of ', colsz, ' processors is being split in a ', &
          'processor mesh of dimensions: ', npi, npj, npk
  end if
  
  call initialize_layout_with_distributed_3D_array( &
       global_sz_i, &
       global_sz_j, &
       global_sz_k, &
       npi, &
       npj, &
       npk, &
       layout1 )
  if( myrank .eq. 0 ) then   
     print *, 'Printing layout1: '
     call sll_view_lims_3D( layout1 )
  end if
  
  call get_local_sz( layout1, local_sz_i, local_sz_j, local_sz_k )        
  SLL_ALLOCATE( local_array1(local_sz_i, local_sz_j, local_sz_k), ierr )
  
  do k=1,local_sz_k
     do j=1,local_sz_j 
        do i=1,local_sz_i
           global_index =  local_to_global_3D( layout1, (/i, j, k/) )
           local_array1(i,j,k) = cos(1.*global_index(1)) * &
                sin(1.*global_index(2)) * global_index(3)
        enddo
     enddo
  enddo
  
  layout2  => new_layout_3D( sll_world_collective )
  call two_power_rand_factorization(colsz, npi, npj, npk)

  if( myrank .eq. 0) then
     print *, 'for the target configuration, '
     print *, 'a total of ', colsz, ' processors is being split in a ', &
          'processor mesh of dimensions: ', npi, npj, npk
  end if

  call initialize_layout_with_distributed_3D_array( &
       global_sz_i, &
       global_sz_j, &
       global_sz_k, &
       npi, &
       npj, &
       npk, &
       layout2 )
        
  if( myrank .eq. 0 ) then   
     print *, 'Printing layout2: '
     call sll_view_lims_3D( layout2 )
  end if
  
  call get_local_sz( layout2, local_sz_i, local_sz_j, local_sz_k )
  SLL_ALLOCATE( local_array2(local_sz_i, local_sz_j, local_sz_k), ierr )
  
  rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, local_array1)
  call apply_remap_3D( rmp3, local_array1, local_array2 )


  ! Compute the difference between data contained in local_array2 and
  ! expected data. This difference must be null in each point
  do k=1,local_sz_k
     do j=1,local_sz_j 
        do i=1,local_sz_i
           global_index =  local_to_global_3D( layout2, (/i, j, k/) )
           local_array2(i,j,k) = local_array2(i,j,k) - cos(1.*global_index(1)) &
                                 * sin(1.*global_index(2)) * global_index(3)
        enddo
     enddo
  enddo

  ! As the difference described above is null in each point, the sum of the corresponding 
  ! absolute values must be null. Each processor compute a local sum and all local sums are
  ! finally added and the result is sent to processor 0 which will check if equal 0 to 
  ! validate the test. (*)
  call sll_collective_reduce_real( sll_world_collective, (/ real(sum(abs(local_array2))) /), &
                                   1, MPI_SUM, 0, sum4test )

  print *, 'Remap operation completed.'
  call flush()
  
  call delete_layout_3D( layout1 )
  call delete_layout_3D( layout2 )
  
  call sll_collective_barrier(sll_world_collective)
  
  SLL_DEALLOCATE_ARRAY(local_array1, ierr)
  SLL_DEALLOCATE_ARRAY(local_array2, ierr)
  
  ok = 1
  if (myrank==0) then
     if (sum4test(1)/=0.) then ! Refer to (*)
        ok = 0 
     endif
  endif
  
!enddo
 
  if (myrank==0) then
     if (ok==1) then
        print*, '"remap" unit test: PASS'
     else
        print*, '"remap" unit test: NOT PASS'
     endif
  endif

! else
 !    print*, 'The number of processors must be a power of 2'
 ! endif
  
  call sll_halt_collective()
  !print *, 'Rank ', myrank, ': TEST COMPLETE'
  !call flush()

#if 0
  call split_interval_randomly(1000,4,limits)
  print *, limits(:)
  call flush()
#endif
  
contains

  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2, n3
    integer   :: expo, expo1, expo2, expo3
    if (is_power_of_two(colsz)) then   
       expo = int(log(real(n))/log(2.))  
       call random_number(rand_real)
       expo1 = int(rand_real*expo)
       call random_number(rand_real)
       expo2 = int(rand_real*(expo-expo1))
       expo3 = expo-(expo1+expo2)
       n1 = 2**expo1
       n2 = 2**expo2
       n3 = 2**expo3
    else
       print*, 'The number of processors must be a power of 2'
    endif
  end subroutine two_power_rand_factorization
  
  subroutine get_local_sz( layout, local_sz_i, local_sz_j, local_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: local_sz_i
    sll_int32 :: local_sz_j
    sll_int32 :: local_sz_k
    sll_int32 :: my_rank
    if( .not. associated(layout) ) then
       print *, 'not-associated layout passed to new_distributed_mesh_3D'
       print *, 'Exiting...'
       STOP
    end if
    my_rank = sll_get_collective_rank(get_layout_3D_collective(layout))
    i_min = get_layout_3D_i_min( layout, my_rank )
    i_max = get_layout_3D_i_max( layout, my_rank )
    j_min = get_layout_3D_j_min( layout, my_rank )
    j_max = get_layout_3D_j_max( layout, my_rank )
    k_min = get_layout_3D_k_min( layout, my_rank )
    k_max = get_layout_3D_k_max( layout, my_rank )
    local_sz_i = i_max - i_min + 1
    local_sz_j = j_max - j_min + 1
    local_sz_k = k_max - k_min + 1
    end subroutine get_local_sz
      
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
    
