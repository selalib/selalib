program remap_test
  use sll_collective
#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "misc_utils.h"
  implicit none
  
  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  integer, dimension(:,:,:), allocatable    :: local_array1, local_array2, arrays_diff
  ! Take a 3D array of dimensions ni*nj*nk
  ! ni, nj, nk: global sizes
  integer , parameter                       :: ni = 512
  integer , parameter                       :: nj = 512
  integer , parameter                       :: nk = 256
  ! Local sizes
  integer                                   :: loc_sz_i_init
  integer                                   :: loc_sz_j_init
  integer                                   :: loc_sz_k_init
  integer                                   :: loc_sz_i_final
  integer                                   :: loc_sz_j_final
  integer                                   :: loc_sz_k_final

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
  integer, parameter                        :: nbtest = 1
  integer                                   :: i_test
  integer                                   :: i, j, k
  sll_int32, dimension(3)                   :: global_indices, g
  sll_real32   , dimension(1)               :: prod4test
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

  ok = 1
  do, i_test=1, nbtest
     if( myrank .eq. 0 ) then
        print *, 'Iteration ', i_test, ' of ', nbtest
     end if
     layout1  => new_layout_3D( sll_world_collective )        
     call two_power_rand_factorization(colsz, npi, npj, npk)
!     call factorize_in_random_2powers( colsz, npi, npj )
!     npk = 1 
     if( myrank .eq. 0 ) then
        print *, 'source configuration: ', npi, npj, npk
     end if

     call initialize_layout_with_distributed_3D_array( &
          ni, &
          nj, &
          nk, &
          npi, &
          npj, &
          npk, &
          layout1 )
     
     call compute_local_sizes( &
          layout1, &
          loc_sz_i_init, &
          loc_sz_j_init, &
          loc_sz_k_init )        

     SLL_ALLOCATE(local_array1(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init),ierr)
 
     ! initialize the local data    
     do k=1,loc_sz_k_init
        do j=1,loc_sz_j_init 
           do i=1,loc_sz_i_init
              global_indices =  local_to_global_3D( layout1, (/i, j, k/) )
              gi = global_indices(1)
              gj = global_indices(2)
              gk = global_indices(3)
              local_array1(i,j,k) = gi + (gj-1)*ni + (gk-1)*ni*nj
           enddo
        enddo
     enddo
     
     layout2  => new_layout_3D( sll_world_collective )
     call two_power_rand_factorization(colsz, npi, npj, npk)
 !    call factorize_in_random_2powers(colsz, npi, npj)
 !    npk = 1    
     if( myrank .eq. 0 ) then
        print *, 'target configuration: ', npi, npj, npk
     end if

     call initialize_layout_with_distributed_3D_array( &
          ni, &
          nj, &
          nk, &
          npi, &
          npj, &
          npk, &
          layout2 )

     call reorganize_randomly(layout2)
     
     call compute_local_sizes( &
          layout2, &
          loc_sz_i_final, &
          loc_sz_j_final, &
          loc_sz_k_final )

     SLL_ALLOCATE( local_array2(loc_sz_i_final, loc_sz_j_final, loc_sz_k_final), ierr )
    
     rmp3 => NEW_REMAPPER_PLAN_3D( layout1, layout2, local_array1)     

     call apply_remap_3D( rmp3, local_array1, local_array2 ) 

     SLL_ALLOCATE(arrays_diff(loc_sz_i_final,loc_sz_j_final,loc_sz_k_final),ierr ) 

     ! compare results with expected data
     do k=1,loc_sz_k_final
        do j=1,loc_sz_j_final 
           do i=1,loc_sz_i_final
              global_indices =  local_to_global_3D( layout2, (/i, j, k/) )
              gi = global_indices(1)
              gj = global_indices(2)
              gk = global_indices(3)
              arrays_diff(i,j,k) = local_array2(i,j,k) - &
                   (gi + (gj-1)*ni + (gk-1)*ni*nj)
              if (arrays_diff(i,j,k)/=0) then
                 ok = 0
                 print*, i_test, myrank, '"remap" unit test: FAIL'
                 print *, 'local indices: ', '(', i, j, k, ')'
                 print*, 'in global indices: (',gi, ',', gj, ',', gk, ')'
                 print*, 'local array1(',gi, ',', gj, ',', gk, ')=', &
                      local_array1(i,j,k)  
                 print*, 'local array2(',gi, ',', gj, ',', gk, ')=', &
                      local_array2(i,j,k)
                 g = theoretical_global_3D_indices(local_array2(i,j,k), ni, nj)
                 print*, 'Theoretical indices: (',g(1), ',', g(2),',',g(3), ')'
                 if(myrank .eq. 1) then
                    print *, local_array2(:,:,:)
                 end if

                 print *, i_test, myrank, 'Printing layout1: '
                 call sll_view_lims_3D( layout1 )
                 print *, i_test, myrank, 'Printing layout2: '
                 call sll_view_lims_3D( layout2 )

                 print*, 'program stopped'
                 stop
              end if
           end do
        end do
     end do

     ! As the difference described above is null in each point, the sum of the 
     ! corresponding absolute values must be null. Each processor compute a local 
     ! sum and all local sums are finally added and the result is sent to 
     ! processor 0 which will check if equal 0 to validate the test. (*)
     call sll_collective_reduce_real( &
          sll_world_collective, &
          (/ real(ok) /), &
          1, &
          MPI_PROD, &
          0, &
          prod4test )

     if( myrank .eq. 0) then
        print *, ' '
        print *, '-------------------------------------------'
        print *, ' '
        call flush()
     end if
     call flush() 
       
     call sll_collective_barrier(sll_world_collective)
  
     call delete_layout_3D( layout1 )
     call delete_layout_3D( layout2 )
     SLL_DEALLOCATE_ARRAY(local_array1, ierr)
     SLL_DEALLOCATE_ARRAY(local_array2, ierr)
     SLL_DEALLOCATE_ARRAY(arrays_diff, ierr)
  enddo

  if (myrank==0) then
     if (prod4test(1)==1.) then
        print*, 'TEST PASSED'
     endif
  endif
  
  call sll_halt_collective()

#if 0
  call split_interval_randomly(1000,4,limits)
  print *, limits(:)
  call flush()
#endif
  
contains

  subroutine reorganize_randomly(layout_3D)
    implicit none
    type(layout_3D_t), pointer :: layout_3D
    integer                    :: i, colsz, proc_n, proc_p
    real                       :: rand_real
    colsz = sll_get_num_nodes(layout_3D)
    do i=1, colsz
       call random_number(rand_real)
       proc_n = int(rand_real*(colsz-1))
       call random_number(rand_real)
       proc_p = int(rand_real*(colsz-1))
       call swap_box_3D(proc_n, proc_p, layout_3D)
    enddo    
  end subroutine reorganize_randomly

  subroutine swap_box_3D(proc_n, proc_p, layout_3D)
    implicit none
    integer                    :: proc_n, proc_p
    type(layout_3D_t), pointer :: layout_3D
    integer                    :: i_min_n, i_max_n, j_min_n, j_max_n, &
    k_min_n, k_max_n, i_min_p, i_max_p, j_min_p, j_max_p, k_min_p, k_max_p    
    ! Get proc_n contents from layout_3D
    i_min_n = get_layout_3D_i_min( layout_3D, proc_n )
    i_max_n = get_layout_3D_i_max( layout_3D, proc_n )
    j_min_n = get_layout_3D_j_min( layout_3D, proc_n )
    j_max_n = get_layout_3D_j_max( layout_3D, proc_n )
    k_min_n = get_layout_3D_k_min( layout_3D, proc_n )
    k_max_n = get_layout_3D_k_max( layout_3D, proc_n )
    ! Get proc_p contents from layout_3D
    i_min_p = get_layout_3D_i_min( layout_3D, proc_p )
    i_max_p = get_layout_3D_i_max( layout_3D, proc_p )
    j_min_p = get_layout_3D_j_min( layout_3D, proc_p )
    j_max_p = get_layout_3D_j_max( layout_3D, proc_p )
    k_min_p = get_layout_3D_k_min( layout_3D, proc_p )
    k_max_p = get_layout_3D_k_max( layout_3D, proc_p )
    ! Set proc_n contents in layout_3D
    call set_layout_3D_i_min( layout_3D, proc_n, i_min_p )
    call set_layout_3D_i_max( layout_3D, proc_n, i_max_p)
    call set_layout_3D_j_min( layout_3D, proc_n, j_min_p )
    call set_layout_3D_j_max( layout_3D, proc_n, j_max_p )
    call set_layout_3D_k_min( layout_3D, proc_n, k_min_p )
    call set_layout_3D_k_max( layout_3D, proc_n, k_max_p )
    ! Set proc_p contents in layout_3D
    call set_layout_3D_i_min( layout_3D, proc_p, i_min_n )
    call set_layout_3D_i_max( layout_3D, proc_p, i_max_n )
    call set_layout_3D_j_min( layout_3D, proc_p, j_min_n )
    call set_layout_3D_j_max( layout_3D, proc_p, j_max_n )
    call set_layout_3D_k_min( layout_3D, proc_p, k_min_n )
    call set_layout_3D_k_max( layout_3D, proc_p, k_max_n )   
  end subroutine swap_box_3D
  
  function theoretical_global_3D_indices(d, ni, nj)
    integer, dimension(1:3) :: theoretical_global_3D_indices
    integer, intent(in)      :: d, ni, nj
    integer                  :: q
#if 0
    integer                  :: val
    val = d/(ni*nj)
    theoretical_global_3D_indices(3) = val
    theoretical_global_3D_indices(2) = 
    theoretical_global_3D_indices(1) = 
#endif
#if 1
    if(mod(d,ni) /= 0) then
       theoretical_global_3D_indices(1) = mod(d,ni)
    else
       theoretical_global_3D_indices(1) = ni
    end if
    q = d/ni
    theoretical_global_3D_indices(2) = mod(q,nj) + 1
    theoretical_global_3D_indices(3) = q/nj + 1
#endif
  end function theoretical_global_3D_indices
  
  subroutine two_power_rand_factorization(n, n1, n2, n3)
    sll_int64, intent(in) :: n
    integer, intent(out) ::n1, n2, n3
    integer   :: expo, expo1, expo2, expo3
    if (.not.is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    expo3 = expo - (expo1+expo2)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
  end subroutine two_power_rand_factorization

  subroutine factorize_in_random_2powers(n, n1, n2)
    sll_int64, intent(in) :: n
    integer, intent(out)  :: n1, n2
    integer   :: expo, expo1, expo2
    if (.not.is_power_of_two(colsz)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = expo - expo1
    n1 = 2**expo1
    n2 = 2**expo2
  end subroutine factorize_in_random_2powers

  
  subroutine compute_local_sizes( layout, loc_sz_i, loc_sz_j, loc_sz_k )
    type(layout_3D_t), pointer :: layout
    sll_int32, intent(out) :: loc_sz_i
    sll_int32, intent(out) :: loc_sz_j
    sll_int32, intent(out) :: loc_sz_k
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
    loc_sz_i = i_max - i_min + 1
    loc_sz_j = j_max - j_min + 1
    loc_sz_k = k_max - k_min + 1
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
        integer, parameter                 :: spread = 25
        ! dangerous... the following should be found in the enclosing scope...
        integer, intent(inout) :: loadi
        integer, dimension(:), pointer :: ans
        if( (hi-lo).eq. 1 ) then
           ans(loadi) = lo
           ans(loadi+1) = hi
           loadi = loadi + 2
        else if (gen .eq. lim) then
           ans(loadi) = lo
           ans(loadi+1) = hi
           loadi = loadi + 2
        else
           ! call random_number(rand)
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
    
