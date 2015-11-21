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
!> @author
!> Michel Mehrenberger
!> @brief 
!> initialization of 1d cartesian mesh from namelist
!> @details
!> <br>
!> Partly generated from mesh_1d_cart.gnml file
!>\code
!> character(len=256) choice "landau"
!>\endcode
!> and num4.clone file
!>\code
!>_1
!>_2
!>_3
!>_4
!>\endcode
!> Default parameters correspond to namelist
!>\code
!> &mesh_1d_cart
!>   choice = "landau"
!> /
!>\endcode
!> and clones as
!>\code
!> &mesh_1d_cart_1
!>   choice_1 = "landau"
!> /
!>\endcode
!>...
!>\code
!> &mesh_1d_cart_4
!>   choice_4 = "landau"
!> /
!>\endcode

!> Examples of calls of interface (generic)
!>\code
!> !print namelist info
!> call sll_s_nml_mesh_1d_cart( &
!>  filename, &
!>  proc_id=sll_get_collective_rank(sll_world_collective))
!>\endcode

  !-----------------------------------------------------------------
  !  SPECIFIC DOCUMENTATION (BEGIN)
  !-----------------------------------------------------------------

!> <br>
!> Possibilities for namelists variables
!>\code
!> choice = "unif" -> sll_m_nml_mesh_1d_unif_cart
!> choice = "landau" -> sll_m_nml_mesh_1d_landau_cart
!> choice = "two_grid" -> sll_m_nml_mesh_1d_two_grid_cart
!>\endcode
!> Examples of calls of interface (specific):
!> <br>
!> 1. Allocation and initialization of 
!> <code>sll_real64, pointer :: array(:)</code>
!> according to <code>choice</code> read from namelist 
!> <code>mesh_1d_cart</code> in <code>filename</code>
!>\code
!> call sll_s_nml_mesh_1d_cart(filename,array)
!>\endcode
!> 2. Same, but read this time, from namelist <code>mesh_1d_cart_1</code>.
!>\code
!> call sll_s_nml_mesh_1d_cart(filename,array,clone="_1")
!>\endcode

  !-----------------------------------------------------------------
  !  SPECIFIC DOCUMENTATION (END)
  !-----------------------------------------------------------------



module sll_m_nml_mesh_1d_cart
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"

  use sll_m_utilities, only : &
    sll_new_file_id

  !-----------------------------------------------------------------
  !  SPECIFIC INCLUDE (BEGIN)
  !-----------------------------------------------------------------

  use sll_m_nml_mesh_1d_landau_cart, only : &
    sll_s_nml_mesh_1d_landau_cart
  use sll_m_nml_mesh_1d_unif_cart, only : &
    sll_s_nml_mesh_1d_unif_cart
  use sll_m_nml_mesh_1d_two_grid_cart, only : &
    sll_s_nml_mesh_1d_two_grid_cart

  !-----------------------------------------------------------------
  !  SPECIFIC INCLUDE (END)
  !-----------------------------------------------------------------

  implicit none

  public :: &
    sll_s_nml_mesh_1d_cart, &
    sll_t_nml_mesh_1d_cart

  private

  type :: sll_t_nml_mesh_1d_cart
    character(len=256) :: choice
    character(len=256) :: label
  contains
    procedure, pass(self) :: init=>init_nml_mesh_1d_cart
    procedure, pass(self) :: init_1=>init_nml_mesh_1d_cart_1
    procedure, pass(self) :: init_2=>init_nml_mesh_1d_cart_2
    procedure, pass(self) :: init_3=>init_nml_mesh_1d_cart_3
    procedure, pass(self) :: init_4=>init_nml_mesh_1d_cart_4
    procedure, pass(self) :: init_clone=>init_clone_nml_mesh_1d_cart
  end type sll_t_nml_mesh_1d_cart

  !-----------------------------------------------------------------
  !  SPECIFIC DECLARATION (BEGIN)
  !-----------------------------------------------------------------

  interface sll_s_nml_mesh_1d_cart
      module procedure &
        s_nml_mesh_1d_cart_array, &
        s_nml_mesh_1d_cart_print
  end interface sll_s_nml_mesh_1d_cart

  !-----------------------------------------------------------------
  !  SPECIFIC DECLARATION (END)
  !-----------------------------------------------------------------

contains

  !-----------------------------------------------------------------
  !  SPECIFIC SUBROUTINES (BEGIN)
  !-----------------------------------------------------------------
  
  !> @brief create 1d array from namelist  
  subroutine s_nml_mesh_1d_cart_array( &
    filename, &
    array, &
    clone, & 
    proc_id)

    character(len=*), intent(in)    :: filename !< namelist file input
    sll_real64, pointer, intent(out) :: array(:) !< output array
    character(len=*), intent(in), optional :: clone !< optional choice of clone
    sll_int32, intent(in), optional :: proc_id !< optional id of proc
    
    type(sll_t_nml_mesh_1d_cart) :: self
    character(len=256) :: err_msg
    character(len=256) :: caller
    
    if(present(clone))then
      call self%init_clone(clone,filename,proc_id)
    else
      call self%init(filename,proc_id)          
    endif


    caller = 's_nml_mesh_1d_cart_array'
    select case (self%choice)
      case("unif")
        call sll_s_nml_mesh_1d_unif_cart( &
          filename, &
          array, &
          clone, &
          proc_id)
      case("landau")
        call sll_s_nml_mesh_1d_landau_cart( &
          filename, &
          array, &
          clone, &
          proc_id)
      case("two_grid")
        call sll_s_nml_mesh_1d_two_grid_cart( &
          filename, &
          array, &
          clone, &
          proc_id)
      case default
        err_msg = 'bad value for self%choice'
        SLL_ERROR(trim(caller),trim(err_msg))
    end select
  
  end subroutine s_nml_mesh_1d_cart_array

  !-----------------------------------------------------------------
  !  SPECIFIC SUBROUTINES (END)
  !-----------------------------------------------------------------

  !> @brief print namelist info  
  subroutine s_nml_mesh_1d_cart_print( &
    filename, &
    clone, &
    proc_id)

    character(len=*), intent(in) :: filename !< namelist file input
    character(len=*), intent(in), optional :: clone !< optional choice of clone
    sll_int32, intent(in), optional :: proc_id !< optional id of proc

    type(sll_t_nml_mesh_1d_cart) :: self
    sll_int32 :: proc_id_loc

    if(present(clone))then
      call self%init_clone(clone,filename,proc_id)
    else
      call self%init(filename,proc_id)
    endif

    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    if(proc_id_loc==0)then
      print *,'#nml_mesh_1d_cart:'

      print *,'#label=',trim(self%label)
      print *,'#choice=',self%choice
    endif

  end subroutine s_nml_mesh_1d_cart_print

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  subroutine init_clone_nml_mesh_1d_cart( &
    self, &
    clone, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in) :: clone
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    character(len=256) :: err_msg
    character(len=256) :: caller

    caller = 'init_clone_nml_mesh_1d_cart'
    select case (clone)
      case("_1")
        call self%init_1(filename,proc_id)
      case("_2")
        call self%init_2(filename,proc_id)
      case("_3")
        call self%init_3(filename,proc_id)
      case("_4")
        call self%init_4(filename,proc_id)
      case default
        err_msg = 'bad value for clone'
        SLL_ERROR(trim(caller),trim(err_msg))
    end select

  end subroutine init_clone_nml_mesh_1d_cart

  subroutine init_nml_mesh_1d_cart( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    character(len=256) :: choice
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_cart/ &
      choice
    caller = 'init_nml_mesh_1d_cart'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      choice)

    call sll_new_file_id(namelist_id, ierr)
    open( &
      unit = namelist_id, &
      file=trim(filename)//'.nml', &
      IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = &
         'failed to open first file '//trim(filename)//'.nml'
       SLL_ERROR(trim(caller),trim(err_msg))
    end if

    read(namelist_id, mesh_1d_cart)
    self%label="no_label"
    self%choice = choice
    close(namelist_id)

  end subroutine init_nml_mesh_1d_cart

  subroutine init_nml_mesh_1d_cart_1( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    character(len=256) :: choice_1
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_cart_1/ &
      choice_1

    caller = 'init_nml_mesh_1d_cart_1'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      choice_1)

    call sll_new_file_id(namelist_id, ierr)
    open( &
      unit = namelist_id, &
      file=trim(filename)//'.nml', &
      IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = &
         'failed to open first file '//trim(filename)//'.nml'
       SLL_ERROR(trim(caller),trim(err_msg))
    end if

    read(namelist_id, mesh_1d_cart_1)
    self%label="_1"
    self%choice = choice_1
    close(namelist_id)

  end subroutine init_nml_mesh_1d_cart_1

  subroutine init_nml_mesh_1d_cart_2( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    character(len=256) :: choice_2
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_cart_2/ &
      choice_2

    caller = 'init_nml_mesh_1d_cart_2'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      choice_2)

    call sll_new_file_id(namelist_id, ierr)
    open( &
      unit = namelist_id, &
      file=trim(filename)//'.nml', &
      IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = &
         'failed to open first file '//trim(filename)//'.nml'
       SLL_ERROR(trim(caller),trim(err_msg))
    end if

    read(namelist_id, mesh_1d_cart_2)
    self%label="_2"
    self%choice = choice_2
    close(namelist_id)

  end subroutine init_nml_mesh_1d_cart_2

  subroutine init_nml_mesh_1d_cart_3( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    character(len=256) :: choice_3
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_cart_3/ &
      choice_3

    caller = 'init_nml_mesh_1d_cart_3'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      choice_3)

    call sll_new_file_id(namelist_id, ierr)
    open( &
      unit = namelist_id, &
      file=trim(filename)//'.nml', &
      IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = &
         'failed to open first file '//trim(filename)//'.nml'
       SLL_ERROR(trim(caller),trim(err_msg))
    end if

    read(namelist_id, mesh_1d_cart_3)
    self%label="_3"
    self%choice = choice_3
    close(namelist_id)

  end subroutine init_nml_mesh_1d_cart_3

  subroutine init_nml_mesh_1d_cart_4( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    character(len=256) :: choice_4
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_cart_4/ &
      choice_4

    caller = 'init_nml_mesh_1d_cart_4'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      choice_4)

    call sll_new_file_id(namelist_id, ierr)
    open( &
      unit = namelist_id, &
      file=trim(filename)//'.nml', &
      IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       err_msg = &
         'failed to open first file '//trim(filename)//'.nml'
       SLL_ERROR(trim(caller),trim(err_msg))
    end if

    read(namelist_id, mesh_1d_cart_4)
    self%label="_4"
    self%choice = choice_4
    close(namelist_id)

  end subroutine init_nml_mesh_1d_cart_4

  subroutine set_default_values( &
    choice)
    character(len=256), intent(inout) :: choice

    choice = "landau"

  end subroutine set_default_values

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_nml_mesh_1d_cart
