
!****************************************************************************
!
! Selalib      
! Module: unit_test.F90
!
!> @brief 
!> remapper unit test
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!****************************************************************************

program remap_test
  use sll_collective, only: sll_boot_collective, sll_halt_collective
  use sll_remapper
!#include "sll_remap.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
  implicit none

  ! Test of the 3D remapper takes a 3D array whose global size Nx*Ny*Nz,
  ! distributed among NPi*NPj*NPk processors.
  sll_real64, dimension(:,:,:), allocatable :: local_array1
  sll_real64, dimension(:,:,:), allocatable :: local_array2
  sll_real64, dimension(:,:,:), allocatable ::  arrays_diff
  ! Take a 3D array of dimensions ni*nj*nk
  ! ni, nj, nk: global sizes
  integer , parameter                       :: ni = 64
  integer , parameter                       :: nj = 64
  integer , parameter                       :: nk = 32
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
  type(layout_3D), pointer                  :: layout1
  type(layout_3D), pointer                  :: layout2
  type(remap_plan_3D_real64), pointer       :: rmp3

  sll_real64                                :: rand_real
  integer, parameter                        :: nbtest = 25
  integer                                   :: i_test
  integer                                   :: i, j, k
  sll_int32, dimension(3)                   :: global_indices, g
  sll_real32   , dimension(1)               :: prod4test
  integer                                   :: ok
  sll_int32, dimension(3)                   :: tmpa

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
     call flush(6)
  end if

  if (.not. is_power_of_two(colsz)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if

  ok = 1
  do, i_test=1, nbtest
     call flush(6)
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
              tmpa(:) = (/i, j, k/)
              global_indices =  local_to_global_3D( layout1, tmpa )
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
    
     rmp3 => NEW_REMAP_PLAN( layout1, layout2, local_array1)     

     call apply_remap_3D( rmp3, local_array1, local_array2 ) 

     SLL_ALLOCATE(arrays_diff(loc_sz_i_final,loc_sz_j_final,loc_sz_k_final),ierr ) 
#if 0
     if( myrank .eq. 0 ) then
        print *, i_test, myrank, 'Printing layout1: '
        call sll_view_lims_3D( layout1 )
        print *, i_test, myrank, 'Printing layout2: '
        call sll_view_lims_3D( layout2 )
     end if
#endif
     ! compare results with expected data
     do k=1,loc_sz_k_final
        do j=1,loc_sz_j_final 
           do i=1,loc_sz_i_final
              tmpa(:) = (/i, j, k/)
              global_indices =  local_to_global_3D( layout2, tmpa )
              gi = global_indices(1)
              gj = global_indices(2)
              gk = global_indices(3)
              arrays_diff(i,j,k) = local_array2(i,j,k) - &
                   (gi + (gj-1)*ni + (gk-1)*ni*nj)
              if (arrays_diff(i,j,k)/=0) then
                 ok = 0
                 print*, i_test, myrank, '"remap" unit test: FAIL'
                 print *, i_test, myrank, 'local indices: ', '(', i, j, k, ')'
                 print *, 'global indices wrt target layout'
                 print*, myrank, 'in global indices: (',gi,',', gj,',', gk,')'

                 print*, myrank, 'local array1(',gi, ',', gj, ',', gk, ')=', &
                      local_array1(i,j,k)  
                 print*, myrank, 'local array2(',gi, ',', gj, ',', gk, ')=', &
                      local_array2(i,j,k)
                 g = theoretical_global_3D_indices(&
                      int(local_array2(i,j,k),i32),&
                      ni, &
                      nj)
                 print*, 'Theoretical indices: (',g(1), ',', g(2),',',g(3), ')'
            !     if(myrank .eq. 1) then
                    print *, local_array2(:,:,:)
             !    end if

                 print *, i_test, myrank, 'Printing layout1: '
                 call sll_view_lims_3D( layout1 )
                 print *, i_test, myrank, 'Printing layout2: '
                 call sll_view_lims_3D( layout2 )

                 print*, 'program stopped by failure'
                 stop
              end if
              call flush(6)
           end do
        end do
     end do

     ! As the difference described above is null in each point, the sum of the 
     ! corresponding absolute values must be null. Each processor compute a local 
     ! sum and all local sums are finally added and the result is sent to 
     ! processor 0 which will check if equal 0 to validate the test. (*)
     call sll_collective_reduce( &
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
        call flush(6)
     end if
     call flush(6) 
       
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
  
contains

  subroutine reorganize_randomly(layout)
    implicit none
    type(layout_3D), pointer   :: layout
    integer                    :: i, colsz, proc_n, proc_p
    real                       :: rand_real
    colsz = sll_get_num_nodes(layout)
    do i=0, colsz-1
       call random_number(rand_real)
       proc_n = int(rand_real*(colsz-1))
       call random_number(rand_real)
       proc_p = int(rand_real*(colsz-1))
       call swap_box_3D(proc_n, proc_p, layout)
    enddo    
  end subroutine reorganize_randomly

  subroutine swap_box_3D(proc_n, proc_p, layout)
    implicit none
    integer                  :: proc_n, proc_p
    type(layout_3D), pointer :: layout
    integer                  :: i_min_n, i_max_n, j_min_n, j_max_n, &
    k_min_n, k_max_n, i_min_p, i_max_p, j_min_p, j_max_p, k_min_p, k_max_p    
    ! Get proc_n contents from layout
    i_min_n = get_layout_3D_i_min( layout, proc_n )
    i_max_n = get_layout_3D_i_max( layout, proc_n )
    j_min_n = get_layout_3D_j_min( layout, proc_n )
    j_max_n = get_layout_3D_j_max( layout, proc_n )
    k_min_n = get_layout_3D_k_min( layout, proc_n )
    k_max_n = get_layout_3D_k_max( layout, proc_n )
    ! Get proc_p contents from layout
    i_min_p = get_layout_3D_i_min( layout, proc_p )
    i_max_p = get_layout_3D_i_max( layout, proc_p )
    j_min_p = get_layout_3D_j_min( layout, proc_p )
    j_max_p = get_layout_3D_j_max( layout, proc_p )
    k_min_p = get_layout_3D_k_min( layout, proc_p )
    k_max_p = get_layout_3D_k_max( layout, proc_p )
    ! Set proc_n contents in layout
    call set_layout_3D_i_min( layout, proc_n, i_min_p )
    call set_layout_3D_i_max( layout, proc_n, i_max_p)
    call set_layout_3D_j_min( layout, proc_n, j_min_p )
    call set_layout_3D_j_max( layout, proc_n, j_max_p )
    call set_layout_3D_k_min( layout, proc_n, k_min_p )
    call set_layout_3D_k_max( layout, proc_n, k_max_p )
    ! Set proc_p contents in layout
    call set_layout_3D_i_min( layout, proc_p, i_min_n )
    call set_layout_3D_i_max( layout, proc_p, i_max_n )
    call set_layout_3D_j_min( layout, proc_p, j_min_n )
    call set_layout_3D_j_max( layout, proc_p, j_max_n )
    call set_layout_3D_k_min( layout, proc_p, k_min_n )
    call set_layout_3D_k_max( layout, proc_p, k_max_n )   
  end subroutine swap_box_3D
  
  function theoretical_global_3D_indices(d, ni, nj) !result(ind)
    !integer, dimension(1:3) :: ind
    integer, dimension(1:3) :: theoretical_global_3D_indices
    integer, intent(in)     :: d, ni, nj
    integer                 :: i, j, k
    integer                 :: q
    sll_real64              :: tmp
    ! We know that the linear index relates with the other
    ! indices by:
    !
    !          linear = i + ni*(j-1) + ni*nj*(k-1)
    !
    ! where i,j,k are the indices sought. We start by working backwards.
    tmp = real(d,f64)/real(ni*nj,f64)
    k   = ceiling(tmp)
    ! reduce the size of the number we are working with
    q   = d - (k-1)*ni*nj
    tmp = real(q,f64)/real(ni,f64)
    j   = ceiling(tmp)
    ! reduce again
    q   = q - (j-1)*ni
    i   = q
    theoretical_global_3D_indices(1:3) = (/i,j,k/) 
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
           mid = (hi+lo)/2 + int(gaussian_dev())*spread
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
    
