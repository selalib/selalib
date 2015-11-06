program test_vector_space

#include "sll_working_precision.h"
  use sll_m_vector_space_base, only: sll_vector_space_base
  use sll_m_vector_space_real_arrays, only: &
    sll_vector_space_real_1d, &
    sll_vector_space_real_2d
 
  implicit none

  sll_real64,             target, allocatable :: a1(:), a2(:,:)
  type( sll_vector_space_real_1d )            :: v1
  type( sll_vector_space_real_2d )            :: v2
  class( sll_vector_space_base ), allocatable :: v3
  sll_int32                                   :: i, j, n, m

  interface test_association
    procedure test_association_1d
    procedure test_association_2d
    procedure test_association_3d
  end interface

  ! Allocate arrays
  n = 3
  m = 5
  allocate( a1(m)   )
  allocate( a2(n,m) )

  !----------------------------------------------------------------------------
  ! Fill in arrays and print them to standard output
  do i=1,m
    a1(i) = real(i,f64)
  end do

  do i=1,n
    do j=1,m
      a2(i,j) = real(i*10 + j,f64)
    end do
  end do

  write(*,'(*(g12.5))') a1
  write(*,*)
  do i = 1,size( a2,1 )
    write(*,'(*(g12.5))') a2(i,:)
  end do
  write(*,*)

  !----------------------------------------------------------------------------
  ! Link vectors to arrays, and print them to std out
  call v1%attach( a1 )
  call v2%attach( a2 )

  write(*,'(*(g12.5))') v1%array
  write(*,*)
  do i = 1,size( v2%array,1 )
    write(*,'(*(g12.5))') v2%array(i,:)
  end do
  write(*,*)

  !----------------------------------------------------------------------------
  ! Allocate a new generic vector by "sourcing" from one of the others
  call v2%source( v3 )
  call v3%mult( 100.0_f64, v2 )
  call v3%incr_mult( 2.0_f64, v2 )

!  v3 = v1
!  call v3%add( v1, v2 )

  select type( v3 ); type is( sll_vector_space_real_2d )

    do i = 1,size( v3%array,1 )
      write(*,'(*(g12.5))') v3%array(i,:)
    end do
    call test_association( v3%array,'v3' )
    call v3%delete()
    call test_association( v3%array,'v3' )

  end select

!==============================================================================
contains
!==============================================================================
  
  !----------------------------------------------------------------------------
  subroutine test_association_1d( v, vname )
    sll_real64, pointer, intent( in ) :: v(:)
    character( len=* ) , intent( in ) :: vname

    if( associated( v ) ) then
      write(*,*) vname, " is associated"
    else
      write(*,*) vname, " is not associated"
    endif

  end subroutine test_association_1d

  !----------------------------------------------------------------------------
  subroutine test_association_2d( v, vname )
    sll_real64, pointer, intent( in ) :: v(:,:)
    character( len=* ) , intent( in ) :: vname

    if( associated( v ) ) then
      write(*,*) vname, " is associated"
    else
      write(*,*) vname, " is not associated"
    endif

  end subroutine test_association_2d

  !----------------------------------------------------------------------------
  subroutine test_association_3d( v, vname )
    sll_real64, pointer, intent( in ) :: v(:,:,:)
    character( len=* ) , intent( in ) :: vname

    if( associated( v ) ) then
      write(*,*) vname, " is associated"
    else
      write(*,*) vname, " is not associated"
    endif

  end subroutine test_association_3d

!==============================================================================
end program test_vector_space
