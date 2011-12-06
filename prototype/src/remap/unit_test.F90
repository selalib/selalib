program remap_test
  use sll_collective
#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
  implicit none
  
  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array1, local_array2
  ! Take a 3D array of dimensions global_sz_i*global_sz_j*global_sz_k
  integer , parameter                       :: global_sz_i = 1024
  integer , parameter                       :: global_sz_j = 1024
  integer , parameter                       :: global_sz_k = 1
  !  logical                                   :: end_result
  ! Local sizes
  integer                                   :: local_sz_i_init
  integer                                   :: local_sz_j_init
  integer                                   :: local_sz_k_init
  integer                                   :: local_sz_i_final
  integer                                   :: local_sz_j_final
  integer                                   :: local_sz_k_final

  ! the process mesh
  integer                                   :: npi
  integer                                   :: npj
  integer                                   :: npk
  sll_int32                                 :: gi, gj, gk
  integer                                   :: ierr
  integer                                   :: myrank
  sll_int64                                 :: colsz        ! collective size
  ! Remap stuff
  type(layout_3D_t), pointer                :: layout1
  type(layout_3D_t), pointer                :: layout2
  type(remap_plan_3D_t), pointer            :: rmp3

  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 10
  integer                                   :: i_test
  integer                                   :: i, j, k
  sll_int32, dimension(1:3)                 :: global_index
  sll_real32   , dimension(1)               :: sum4test
  integer                                   :: ok

  ! Boot parallel environment
  call sll_boot_collective()
!  end_result = .true.
  colsz  = sll_get_collective_size(sll_world_collective)
  myrank = sll_get_collective_rank(sll_world_collective)

  if( myrank .eq. 0) then
     print *, ' '
     print *, '--------------- REMAP test ---------------------'
     print *, ' '
     print *, 'Running a test on ', colsz, 'processes'
     call flush()
  end if

  if (.not. is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if

  do, i_test=1, nbtest
 !    call sll_collective_barrier(sll_world_collective)
     if( myrank .eq. 0 ) then
        print *, i_test, myrank, 'Iteration ', i_test, ' of ', nbtest
        call flush()
     end if
     layout1  => new_layout_3D( sll_world_collective )        
!     call two_power_rand_factorization(colsz, npi, npj, npk)
!     call factorize_in_random_2powers( colsz, npi, npj )
npi = 8
npj = 4
     npk = 1 

     if( myrank .eq. 0) then
        print *, i_test, myrank, 'for the initial configuration, '
        print *, i_test, myrank, 'a total of ', colsz, &
             ' processors is being split in a ', &
             'processor mesh of dimensions: ', npi, npj, npk
        call flush()
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
        print *, i_test, myrank, 'Printing layout1: '
        call sll_view_lims_3D( layout1 )
     end if
     
     call compute_local_sizes( &
          layout1, &
          local_sz_i_init, &
          local_sz_j_init, &
          local_sz_k_init )        
 !    print *, 'proceeding to allocate dimensions: ', local_sz_i, local_sz_j, &
 !         local_sz_k
     SLL_ALLOCATE( local_array1(local_sz_i_init, local_sz_j_init, local_sz_k_init), ierr )
 
     ! initialize the local data    
     do k=1,local_sz_k_init
        do j=1,local_sz_j_init 
           do i=1,local_sz_i_init
              global_index =  local_to_global_3D( layout1, (/i, j, k/) )
              gi = global_index(1)
              gj = global_index(2)
              gk = global_index(3)
              local_array1(i,j,k) = cos(1.0*gi) * sin(1.0*gj) * gk
           enddo
        enddo
     enddo
     
     layout2  => new_layout_3D( sll_world_collective )
!     call two_power_rand_factorization(colsz, npi, npj, npk)
!     call factorize_in_random_2powers(colsz, npi, npj)
npi = 2
npj = 16
     npk = 1    
     if( myrank .eq. 0) then
        print *, i_test, myrank, 'for the target configuration, '
        print *, i_test, myrank, 'a total of ', colsz, &
             ' processors is being split in a ', &
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
        print *, i_test, myrank, 'Printing layout2: '
        call sll_view_lims_3D( layout2 )
     end if
     
     call compute_local_sizes( &
          layout2, &
          local_sz_i_final, &
          local_sz_j_final, &
          local_sz_k_final )
     SLL_ALLOCATE( local_array2(local_sz_i_final, local_sz_j_final, local_sz_k_final), ierr )
  
     ! Do the remap
     print *, i_test, myrank, 'from rank: ', myrank, &
          'Proceeding to build remap plan'
     call flush()     
     rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, local_array1)
     print *, i_test, myrank, 'from rank: ', myrank, &
     'Proceeding to apply remap'
     call flush()     

     call apply_remap_3D( rmp3, local_array1, local_array2 )
     print *, i_test, myrank, 'from rank: ', myrank, 'Completed remap'
     call flush()     

     ! compare results with expected data
     do k=1,local_sz_k_final
        do j=1,local_sz_j_final 
           do i=1,local_sz_i_final
              global_index =  local_to_global_3D( layout2, (/i, j, k/) )
              gi = global_index(1)
              gj = global_index(2)
              gk = global_index(3)
              local_array2(i,j,k) = local_array2(i,j,k) - &
              cos(1.0*gi)*sin(1.0*gj)*gk
           end do
        end do
     end do
     
     !print *, 'from rank ', myrank, 'differences = ' ,local_array2(:,:,:)
     call sll_collective_reduce_real( &
          sll_world_collective, &
          (/ real(sum(abs(local_array2))) /), &
          1, &
          MPI_SUM, &
          0, &
          sum4test )
     if( myrank .eq. 0 ) then
        print *, i_test, myrank, 'result of reduction = ', sum4test
     end if

     print *, i_test, myrank, 'rank: ', myrank, 'Remap operation completed.'
     call flush()
  
     if (myrank==0) then
        if (sum4test(1)==0.) then
           print*, ' '
           print*, i_test, myrank, ' "remap" unit test: PASS'
           print*, ' '
        else
           print*, ' '
           print*, i_test, myrank, '"remap" unit test: FAIL'
           print*, ' '
        endif
     endif
     call delete_layout_3D( layout1 )
     call delete_layout_3D( layout2 )
     SLL_DEALLOCATE_ARRAY(local_array1, ierr)
     SLL_DEALLOCATE_ARRAY(local_array2, ierr)

!     call sll_collective_barrier(sll_world_collective)
     print *, i_test, myrank, 'Process ', myrank, ' finishing loop.'
     call flush()
!     call sll_collective_barrier(sll_world_collective)
     if( myrank .eq. 0) then
        print *, ' '
        print *, '-------------------------------------------'
        print *, ' '
        call flush()
     end if
     call flush() 
     call sll_collective_barrier(sll_world_collective)
 
  enddo
 

  call sll_collective_barrier(sll_world_collective)

  ! As the difference described above is null in each point, the sum of the 
  ! corresponding absolute values must be null. Each processor compute a local 
  ! sum and all local sums are finally added and the result is sent to 
  ! processor 0 which will check if equal 0 to validate the test. (*)
  call sll_collective_reduce_real( &
       sll_world_collective, &
       (/ real(sum(abs(local_array2))) /), &
       1, &
       MPI_SUM, &
       0, &
       sum4test )

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
       stop
    endif
  end subroutine two_power_rand_factorization

  subroutine factorize_in_random_2powers(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  :: n1, n2
    integer   :: expo, expo1, expo2
    if (is_power_of_two(colsz)) then   
       expo = int(log(real(n))/log(2.))  
       call random_number(rand_real)
       expo1 = int(rand_real*expo)
       call random_number(rand_real)
       expo2 = expo -expo1
       n1 = 2**expo1
       n2 = 2**expo2
    else
       print*, 'The number of processors must be a power of 2'
       stop
    endif
  end subroutine factorize_in_random_2powers

  
  subroutine compute_local_sizes( layout, local_sz_i, local_sz_j, local_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(out) :: local_sz_i
    sll_int32, intent(out) :: local_sz_j
    sll_int32, intent(out) :: local_sz_k
    sll_int32 :: i_min
    sll_int32 :: i_max
    sll_int32 :: j_min
    sll_int32 :: j_max
    sll_int32 :: k_min
    sll_int32 :: k_max
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
  end subroutine compute_local_sizes
      
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
    
