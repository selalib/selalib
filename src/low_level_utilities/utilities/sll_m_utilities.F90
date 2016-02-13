!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup utilities
!> Some common numerical utilities
module sll_m_utilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_p_byte_size, &
    sll_s_compute_bloc, &
    sll_s_compute_mesh_from_bloc, &
    sll_s_display_matrix_2d_integer, &
    sll_s_int2string, &
    sll_f_is_even, &
    sll_f_is_power_of_two, &
    sll_s_mpe_decomp1d, &
    sll_s_pfenvelope, &
    sll_o_display, &
    sll_o_factorial, &
    sll_s_new_file_id

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  intrinsic :: selected_int_kind ! this line gives an error, why?

  ! Tentative implementation of a standard-compliant way to get the
  ! memory footprint of a variable. This is our yardstick...
  !
  ! selected_int_kind(r) returns the default integer scalar that is the kind 
  ! type parameter value for an integer data type able to represent all 
  ! integer values n in the range -10^(-r) < n < 10^r, where r is a scalar 
  ! integer. If more than one is available, a kind with least decimal exponent 
  ! range is chosen (and least kind value if several have least decimal 
  ! exponent range). If no corresponding kind is availalble, the result is -1. 
  ! (Metcalf & Reid. F90/95 2nd ed. p. 176).
  !
  ! We are, maybe dangerously, relying on the common practice of many compilers
  ! of using the kind values to indicate the number of bytes of storage
  ! occupied by a value. But this is not mandated by the standard. For 
  ! instance a compiler that only has available 4-byte integers may still
  ! support kind values of 1, 2 and 4 to 'ease portability' from other 
  ! platforms. The size of k1 will become our basic 'yardstick' to measure
  ! the size of a memory footprint of a variable. When we ask for the size
  ! of 'var', the answer will be given in terms of how many 'yardsticks'
  ! are needed to represent 'var'. The implications of this assumption
  ! need to be checked further.

  !> kind of integers
  integer, parameter :: sll_p_byte_size = selected_int_kind(0)

  !> Functions to display on screen matrix or vector
  interface sll_o_display
     module procedure sll_s_display_matrix_2d_integer
     module procedure display_matrix_2d_real
     module procedure display_vector_integer
     module procedure display_vector_real
  end interface sll_o_display

  !> @param logical variable used to print time history
  logical :: flag = .true.

  !> Return factorial
  interface sll_o_factorial
     module procedure factorial_int32, factorial_int64
  end interface sll_o_factorial

contains

  !> Check if an integer is equal to \f[2^n\f]
  function sll_f_is_power_of_two( n )
    sll_int64, intent(in) :: n
    logical               :: sll_f_is_power_of_two

    intrinsic :: not, iand

    if( (n>0) .and. (0 .eq. (iand(n,(n-1)))) ) then
       sll_f_is_power_of_two = .true.
    else
       sll_f_is_power_of_two = .false.
    end if
  end function sll_f_is_power_of_two


  !> Check if an integer is even
  function sll_f_is_even( n )
    sll_int32, intent(in) :: n
    logical               :: sll_f_is_even

    intrinsic :: modulo

    if( modulo(n,2) .eq. 0 ) then
       sll_f_is_even = .true.
    else
       sll_f_is_even = .false.
    end if
  end function sll_f_is_even


  !> It would have been nice to declare the next functions as 'pure' functions,
  !> but it is safer to be able to indicate when their arguments have fallen
  !> out of range as this number is so limited anyway.
  function factorial_int32(n) result(fac)
    sll_int32, intent(in) :: n
    sll_int64             :: fac

    sll_int64 :: acc
    sll_int64 :: i

    if(n < 0) then
       print *, 'ERROR, factorial_int32(): n < 0 or n > 20, if latter, ', &
            'function will overflow a 64-bit integer. n =  ', n
    end if

    acc = 1
    if( n >= 1 ) then
       do i=int(n,kind=i64),1,-1
          acc = acc*i
       end do
    end if
    ! case n == 0 is already taken care of. No protection for negative input.
    fac = acc
  end function factorial_int32


  !> It would have been nice to declare the next functions as 'pure' functions,
  !> but it is safer to be able to indicate when their arguments have fallen
  !> out of range as this number is so limited anyway.
  function factorial_int64(n) result(fac)
    sll_int64, intent(in) :: n
    sll_int64             :: fac

    sll_int64 :: acc
    sll_int64 :: i

    if( (n < 0) .or. (n > 20) ) then
       print *, 'ERROR, factorial_int64(): either a negative n was passed: ', &
            'or n > 20, which will overflow a 64-bit integer. n = ', n
    end if

    acc = 1
    if( n >= 1 ) then
       do i=n,1,-1
          acc = acc*i
       end do
    end if
    ! case n == 0 is already taken care of. 
    fac = acc
  end function factorial_int64


  !> Convert an integer < 9999 to a 4 characters string
  subroutine sll_s_int2string( istep, cstep )
    integer         , intent(in ) :: istep   !< input integer
    character(len=4), intent(out) :: cstep   !< output string

    character(len=1) :: aa, bb, cc, dd
    integer          :: kk1, kk2, kk3, kk4

    if ( istep >= 0 .and. istep < 10000) then
       kk1 = istep/1000
       aa  = char(kk1 + 48)
       kk2 = (istep - kk1*1000)/100
       bb  = char(kk2 + 48)
       kk3 = (istep - (kk1*1000) - (kk2*100))/10
       cc  = char(kk3 + 48)
       kk4 = (istep - (kk1*1000) - (kk2*100) - (kk3*10))/1
       dd  = char(kk4 + 48)
       cstep = aa//bb//cc//dd
    else
       SLL_WARNING( 'sll_s_int2string', 'index is negative or greater than 9999' )
       print*, 'index =', istep
       cstep = 'xxxx'
    end if

  end subroutine sll_s_int2string


  !> Get a file unit number free before creating a file
  subroutine sll_s_new_file_id( file_id, error )
    sll_int32, intent(out) :: file_id   !< file unit number
    sll_int32, intent(out) :: error     !< error code

    logical :: lopen
      
    error=1

    do 100 file_id=20,99
  
       inquire(unit=file_id,opened=lopen)
       if(lopen) then
          cycle
       else
          open(file_id,status='SCRATCH',err=100)
          close(file_id,status='DELETE',err=100)
          error=0
          exit
       end if
 
    100 continue

    !SLL_ASSERT(error == 0)
   
  end subroutine sll_s_new_file_id


  !> Display a vector to screen
  subroutine display_vector_real( array, real_format )
   sll_real64, dimension(:), intent(in) :: array
   character(len=*)        , intent(in) :: real_format

   character(len=20) :: display_format
   sll_int32         :: n, i

   n = size(array,1)

   write(display_format, "('(''|''',i4,a,''' |'')')") n, real_format

   write(*,*)
   write(*,display_format) (array(i), i = 1, n)
   write(*,*)

  end subroutine display_vector_real


  !> Display a vector to screen
  subroutine display_vector_integer( array, integer_format )
   sll_int32, dimension(:), intent(in) :: array
   character(len=*)       , intent(in) :: integer_format

   character(len=20) :: display_format
   sll_int32         :: n, i

   n = size(array)

   write(display_format, "('(''|''',i4,a,''' |'')')") n, integer_format

   write(*,*)
   write(*,display_format) (array(i), i = 1, n)
   write(*,*)

  end subroutine display_vector_integer


  !> Display matrix to screen
  subroutine display_matrix_2d_real( array, real_format )
   sll_real64, dimension(:,:), intent(in) :: array
   character(len=*)          , intent(in) :: real_format

   character(len=20) :: display_format
   sll_int32         :: n1, n2, i, j

   n1 = size(array,1)
   n2 = size(array,2)

   write(display_format, "('(''|''',i4,a,''' |'')')") n2, real_format

   write(*,*)
   do i = 1, n1
      write(*,display_format) (array(i,j), j = 1, n2)
   end do
   write(*,*)

  end subroutine display_matrix_2d_real


  !> Display matrix to screen
  subroutine sll_s_display_matrix_2d_integer( array, integer_format )
   sll_int32, dimension(:,:), intent(in) :: array
   character(len=*)         , intent(in) :: integer_format

   character(len=20) :: display_format
   sll_int32         :: n1, n2, i, j

   n1 = size(array,1)
   n2 = size(array,2)

   write(display_format, "('(''|''',i4,a,''' |'')')") n2, integer_format

   write(*,*)
   do i = 1, n1
      write(*,display_format) (array(i,j), j = 1, n2)
   end do
   write(*,*)

  end subroutine sll_s_display_matrix_2d_integer


!> Subroutine to open data file for slv2d and
!> create the thf.dat file to write results
subroutine initialize_file( data_file_id, thf_file_id )
  sll_int32, intent(out) :: data_file_id !< namelist file for slv2d
  sll_int32, intent(out) :: thf_file_id  !< thf file for energy plot

  character(len=*), parameter :: this_sub_name = "initialize_file"
  character(len=72)           :: filename
  integer                     :: IO_stat
  sll_int32                   :: error

  call get_command_argument( 1, filename)

  call sll_s_new_file_id(data_file_id, error)
  open(data_file_id,file=trim(filename),IOStat=IO_stat)
  if (IO_stat/=0) SLL_ERROR( this_sub_name, "Miss argument file.nml" )

  call sll_s_new_file_id(thf_file_id, error)
  open(thf_file_id,file="thf.dat",IOStat=IO_stat, position='append')
  if (IO_stat/=0) SLL_ERROR( this_sub_name, "Cannot open file thf.dat" )

  rewind(thf_file_id)
  close(thf_file_id)

end subroutine initialize_file
 

!> Routine from slv2d to write diagnostics
subroutine time_history( file_id, desc, fformat, array )
   sll_int32               , intent(in) :: file_id !< file unit number
   character(3)            , intent(in) :: desc    !< name of the diagnostics
   character(14)           , intent(in) :: fformat !< fortran output format
   sll_real64, dimension(:), intent(in) :: array   !< data array
    
   if (desc(1:3)=="thf") then
      open(file_id,file="thf.dat",position='append')
      if (flag) then
         rewind(file_id)
         flag = .false.
      end if
      write(file_id,fformat) array
      close(file_id)
   else
      write(*,*) desc," not recognized"
   endif
    
end subroutine time_history

!------------------------------------------------------------------------
!>  @brief
!>  From the MPE library
!>  @details
!>  This file contains a routine for producing a decomposition of a 1-d 
!>  array when given a number of processors.  It may be used in "direct" 
!>  product decomposition.  The values returned assume a "global" domain 
!>  in [1:n]
!------------------------------------------------------------------------
subroutine sll_s_mpe_decomp1d( n, numprocs, myid, s, e)
   sll_int32, intent(in)  :: n
   sll_int32, intent(in)  :: numprocs
   sll_int32, intent(in)  :: myid
   sll_int32, intent(out) :: s
   sll_int32, intent(out) :: e

   sll_int32 :: nlocal
   sll_int32 :: deficit

   nlocal  = n / numprocs
   s       = myid * nlocal + 1
   deficit = mod(n,numprocs)
   s       = s + min(myid,deficit)
   if (myid  < deficit) then
       nlocal = nlocal + 1
   endif
   e = s + nlocal - 1
   if (e  >  n .or. myid == numprocs-1) e = n

end subroutine sll_s_mpe_decomp1d


!> S: the wave form at a given point in time. This wave form is 
!>    not scaled (its maximum value is 1).
!> t: the time at which the envelope is being evaluated
!> tflat, tL, tR, twL, twR, tstart, t0: the parameters defining the
!>    envelope, defined in the main portion of this program.
!> turn_drive_off: 1 if the drive should be turned off after a time
!>    tflat, and 0 otherwise
subroutine sll_s_pfenvelope(S,               &
                      t,               &
                      tflat,           &
                      tL,              &
                      tR,              &
                      twL,             &
                      twR,             &
                      t0,              &
                      turn_drive_off)

  sll_real64, intent(out) :: S
  sll_real64, intent(in)  :: t
  sll_real64, intent(in)  :: tflat
  sll_real64, intent(in)  :: tL
  sll_real64, intent(in)  :: tR
  sll_real64, intent(in)  :: twL
  sll_real64, intent(in)  :: twR
  sll_real64, intent(in)  :: t0
  logical,    intent(in)  :: turn_drive_off

  sll_real64 :: epsilon

  ! The envelope function is defined such that it is zero at t0,
  ! rises to 1 smoothly, stay constant for tflat, and returns
  ! smoothly to zero.
  if (turn_drive_off) then
     epsilon = 0.5*(tanh((t0-tL)/twL) - tanh((t0-tR)/twR))
     S = 0.5*(tanh((t-tL)/twL) - tanh((t-tR)/twR)) - epsilon
     S = S / (1-epsilon)
  else
     epsilon = 0.5*(tanh((t0-tL)/twL) + 1.0_f64)
     S = 0.5*(tanh((t-tL)/twL) + 1.0_f64) - epsilon
     S = S / (1.0_f64-epsilon)
  endif
  if (S<0) then
     S = 0.0_f64
  endif
  S = S + 0.*tflat ! for use of unused
  return

end subroutine sll_s_pfenvelope


!> - Input: 
!>  + a=bloc_coord(1) b=bloc_coord(2)
!>  + (a,b) subset (0,1) is the refine zone
!>  + bloc_index(1) = density of points in (0,a) 
!>  + bloc_index(2) = density of points in (a,b) 
!>  + bloc_index(3) = density of points in (b,1)
!>
!> - Output:
!>  + 0<=i1<i1+N_fine<=N and x(i1)=a, x(i1+N_fine)=b (approx), x(0)=0, x(N)=1
!>  + bloc_coord(1) = x(i1)
!>  + bloc_coord(2) = x(i1+N_fine)
!>  + bloc_index(1) = i1 
!>  + bloc_index(2) = N_fine 
!>  + bloc_index(3) = N-i1-N_fine
subroutine sll_s_compute_bloc( bloc_coord, bloc_index, N )

  sll_real64, intent(inout)  :: bloc_coord(2)
  sll_int32,  intent(inout)  :: bloc_index(3)
  sll_int32,  intent(in)     :: N

  sll_real64 :: a,b
  sll_int32  :: i1,i2,N_coarse,N_local,N_fine
  
  a=bloc_coord(1)
  b=bloc_coord(2)
  
  !case of uniform mesh with refined zone
  !we have a coarse mesh with N_coarse
  !N=i1+N_local*(i2-i1)+N_coarse-i2
  !N_fine=N_local*(i2-i1)
  !x(i1)=i1/N_coarse x(i1+N_fine)=i2/N_coarse

  if ((bloc_index(1)==1).and.(bloc_index(3)==1)) then      

    N_local = bloc_index(2)
    N_coarse = floor(real(N,f64)/(1._f64+(b-a)*(real(N_local,f64)-1._f64)))
    if (N_local/=1) then
      i2 = (N-N_coarse)/(N_local-1)
    else
      i2 = floor((b-a)*N_coarse)  
    endif   
    N_coarse      = N-i2*(N_local-1)
    i1            = floor(a*N_coarse)
    i2            = i2+i1
    bloc_index(1) = i1
    N_fine        = N_local*(i2-i1)
    bloc_index(2) = N_fine
    bloc_index(3) = N-i1-N_fine
    bloc_coord(1) = real(i1,f64)/real(N_coarse,f64)
    bloc_coord(2) = real(i2,f64)/real(N_coarse,f64)
         
    print *,'#uniform fine mesh would be:',N_coarse*N_local
    print *,'#N_coarse=',N_coarse
    print *,'#saving:',real(N,f64)/real(N_coarse*N_local,f64)
    print *,'#new x(i1),x(i1+N_fine)=',bloc_coord(1),bloc_coord(2)
    print *,'#error for x(i1),x(i1+N_fine)=',bloc_coord(1)-a,bloc_coord(2)-b
    print *,'#i1,i1+N_fine,N_fine,N=',i1,i1+N_fine,N_fine,N

  else

    print*, 'case in sll_s_compute_bloc not implemented yet'

  endif
  
end subroutine sll_s_compute_bloc


!> - Input:   
!>   + x1=bloc_coord(1),x2=bloc_coord(2)
!>   + with 0<i1<i2<N i1=bloc_index(1), i2=i1+bloc_index(2)
!>   + N=bloc_index(1)+bloc_index(2)+bloc_index(3)
!> 
!> Output:  
!>   + node_positions(1:N+1)
!>   + with constraints node_positions(i1+1)=x1,node_positions(i2+1)=x2
!>   + node_positions(1)=0, node_positions(N+1)=1
subroutine sll_s_compute_mesh_from_bloc( bloc_coord, bloc_index, node_positions )

  sll_real64,               intent(in)  :: bloc_coord(2)
  sll_int32,                intent(in)  :: bloc_index(3)
  sll_real64, dimension(:), intent(out) :: node_positions

  sll_int32  :: i, i1, i2, N
  sll_real64 :: dx
  
  N                     = bloc_index(1)+bloc_index(2)+bloc_index(3)
  i1                    = bloc_index(1)
  i2                    = i1+bloc_index(2)
  node_positions(1:N+1) = -1._f64
  node_positions(1)     = 0._f64
  node_positions(i1+1)  = bloc_coord(1)
  node_positions(i2+1)  = bloc_coord(2)
  node_positions(N+1)   = 1._f64
  
  !piecewise linear mapping (maybe enhanced like in complete mesh)
  if(bloc_index(1).ne. 0)then
    dx=bloc_coord(1)/real(bloc_index(1),f64)
    do i=2,bloc_index(1)
      node_positions(i) = (real(i,f64)-1._f64)*dx
    enddo
  endif
  if(bloc_index(2).ne.0)then  
    dx=(bloc_coord(2)-bloc_coord(1))/real(bloc_index(2),f64)
    do i=2,bloc_index(2)
      node_positions(i+i1)=bloc_coord(1)+(real(i,f64)-1._f64)*dx
    enddo
  endif
  if(bloc_index(3).ne.0)then  
    dx=(1._f64-bloc_coord(2))/real(bloc_index(3),f64)
    do i=2,bloc_index(3)
      node_positions(i+i2)=bloc_coord(2)+(real(i,f64)-1._f64)*dx
    enddo
  endif      
end subroutine sll_s_compute_mesh_from_bloc


end module sll_m_utilities
