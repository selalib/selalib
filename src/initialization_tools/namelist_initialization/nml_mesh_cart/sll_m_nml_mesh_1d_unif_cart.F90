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
!> initialization of 1d uniform cartesian mesh from namelist
!> @details
!> <br>
!> Partly generated from mesh_1d_unif_cart.gnml file
!>\code
!>sll_int32 num_cells 32
!>sll_real64 eta_min 0._f64
!>sll_real64 eta_max 1._f64
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
!> &mesh_1d_unif_cart
!>   num_cells = 32
!>   eta_min = 0.
!>   eta_max = 1.
!> /
!>\endcode
!> and clones as
!>\code
!> &mesh_1d_unif_cart_1
!>   num_cells_1 = 32
!>   eta_min_1 = 0.
!>   eta_max_1 = 1.
!> /
!>\endcode
!>...
!>\code
!> &mesh_1d_unif_cart_4
!>   num_cells_4 = 32
!>   eta_min_4 = 0.
!>   eta_max_4 = 1.
!> /
!>\endcode

!> Examples of calls of interface (generic)
!>\code
!> !print namelist info
!> call sll_s_nml_mesh_1d_unif_cart( &
!>  filename, &
!>  proc_id=sll_get_collective_rank(sll_world_collective))
!>\endcode

  !-----------------------------------------------------------------
  !  SPECIFIC DOCUMENTATION (BEGIN)
  !-----------------------------------------------------------------

!> <br>
!> Possibilities for namelists variables
!>\code
!> num_cells : integer >= 1
!> eta_min : real 
!> eta_max : real > eta_min
!>\endcode
!> The 1d mesh is <code>[eta_min,eta_max]</code> discretized in <code>num_cells</code> cells
!> <br>
!> Examples of calls of interface (specific):
!> <br>
!> 1. Allocation and initialization of 
!> <code>sll_real64, pointer :: array(:)</code>
!> according to <code>choice</code> read from namelist 
!> <code>mesh_1d_unif_cart</code> in <code>filename</code>
!>\code
!> call sll_s_nml_mesh_1d_unif_cart(filename,array)
!>\endcode
!> The output <code>array</code> is of size <code>num_cells+1</code>
!> with  
!>\code
!> array(i) = eta_min+(i-1)*(eta_max-eta_min), i=1,num_cells+1
!>\endcode
!> 2. Same, but read this time, from namelist <code>mesh_1d_unif_cart_1</code>.
!>\code
!> call sll_s_nml_mesh_1d_unif_cart(filename,array,clone="_1")
!>\endcode
!> 3. Allocation and initialization of 
!> <code>type(sll_cartesian_mesh_1d), pointer :: mesh</code>
!> according to <code>choice</code> read from namelist 
!> <code>mesh_1d_unif_cart</code> in <code>filename</code>
!>\code
!> call sll_s_nml_mesh_1d_unif_cart(filename,mesh)
!>\endcode
!> 4. Same, but read this time, from namelist <code>mesh_1d_unif_cart_1</code>.
!>\code
!> call sll_s_nml_mesh_1d_unif_cart(filename,mesh,clone="_1")
!>\endcode

  !-----------------------------------------------------------------
  !  SPECIFIC DOCUMENTATION (END)
  !-----------------------------------------------------------------

module sll_m_nml_mesh_1d_unif_cart
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"

  use sll_m_utilities, only : &
    sll_new_file_id

  !-----------------------------------------------------------------
  !  SPECIFIC INCLUDE (BEGIN)
  !-----------------------------------------------------------------

  use sll_m_cartesian_meshes, only : &
    sll_cartesian_mesh_1d, &
    new_cartesian_mesh_1d  

  !-----------------------------------------------------------------
  !  SPECIFIC INCLUDE (END)
  !-----------------------------------------------------------------

  implicit none

  public :: &
    sll_s_nml_mesh_1d_unif_cart, &
    sll_t_nml_mesh_1d_unif_cart

  private

  type :: sll_t_nml_mesh_1d_unif_cart
    sll_int32 :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    character(len=256) :: label
  contains
    procedure, pass(self) :: init=>init_nml_mesh_1d_unif_cart
    procedure, pass(self) :: init_1=>init_nml_mesh_1d_unif_cart_1
    procedure, pass(self) :: init_2=>init_nml_mesh_1d_unif_cart_2
    procedure, pass(self) :: init_3=>init_nml_mesh_1d_unif_cart_3
    procedure, pass(self) :: init_4=>init_nml_mesh_1d_unif_cart_4
    procedure, pass(self) :: init_clone=>init_clone_nml_mesh_1d_unif_cart
  end type sll_t_nml_mesh_1d_unif_cart

  !-----------------------------------------------------------------
  !  SPECIFIC DECLARATION (BEGIN)
  !-----------------------------------------------------------------

  interface sll_s_nml_mesh_1d_unif_cart
      module procedure &
        s_nml_mesh_1d_unif_cart_array, &
        s_nml_mesh_1d_unif_cart_mesh, &
        s_nml_mesh_1d_unif_cart_print
  end interface sll_s_nml_mesh_1d_unif_cart

  !-----------------------------------------------------------------
  !  SPECIFIC DECLARATION (END)
  !-----------------------------------------------------------------

contains

  !-----------------------------------------------------------------
  !  SPECIFIC SUBROUTINES (BEGIN)
  !-----------------------------------------------------------------

  !> @brief create 1d array from namelist
  subroutine s_nml_mesh_1d_unif_cart_array( &
    filename, & !< namelist file input
    array, & !< output array
    clone, & !< optional choice of clone
    proc_id) !< optional id of proc

    character(len=*), intent(in)    :: filename
    sll_real64, pointer, intent(out) :: array(:)
    character(len=*), intent(in), optional :: clone
    sll_int32, intent(in), optional :: proc_id
    
    type(sll_t_nml_mesh_1d_unif_cart) :: self
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    sll_real64 :: delta_eta
    
    if(present(clone))then
      call self%init_clone(clone,filename,proc_id)
    else
      call self%init(filename,proc_id)          
    endif
    
    num_cells = self%num_cells
    eta_min = self%eta_min
    eta_max = self%eta_max
    delta_eta = (eta_max-eta_min)/real(num_cells,f64)
    
    SLL_ALLOCATE(array(num_cells+1),ierr)
    do i=1,num_cells+1
      array(i) = eta_min+real(i-1,f64)*delta_eta
    enddo
    
  
  end subroutine s_nml_mesh_1d_unif_cart_array

  !> @brief create 1d (uniform) cartesian mesh from namelist
  subroutine s_nml_mesh_1d_unif_cart_mesh( &
    filename, & !< namelist file input
    mesh, & !< output mesh
    clone, & !< optional choice of clone
    proc_id) !< optional id of proc

    character(len=*), intent(in)    :: filename
    type(sll_cartesian_mesh_1d), pointer, intent(out) :: mesh
    character(len=*), intent(in), optional :: clone
    sll_int32, intent(in), optional :: proc_id
    
    type(sll_t_nml_mesh_1d_unif_cart) :: self
    sll_int32 :: ierr
    sll_int32 :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    
    if(present(clone))then
      call self%init_clone(clone,filename,proc_id)
    else
      call self%init(filename,proc_id)          
    endif
    
    num_cells = self%num_cells
    eta_min = self%eta_min
    eta_max = self%eta_max

    mesh => new_cartesian_mesh_1d( &
      num_cells, &
      eta_min=eta_min, &
      eta_max=eta_max)
          
  end subroutine s_nml_mesh_1d_unif_cart_mesh

  !-----------------------------------------------------------------
  !  SPECIFIC SUBROUTINES (END)
  !-----------------------------------------------------------------

  !> @brief print namelist info  
  subroutine s_nml_mesh_1d_unif_cart_print( &
    filename, & !< namelist file input
    clone, & !< optional choice of clone
    proc_id & !< optional id of proc
  )
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: clone
    sll_int32, intent(in), optional :: proc_id

    type(sll_t_nml_mesh_1d_unif_cart) :: self
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
      print *,'#nml_mesh_1d_unif_cart:'

      print *,'#label=',trim(self%label)
      print *,'#num_cells=',self%num_cells
      print *,'#eta_min=',self%eta_min
      print *,'#eta_max=',self%eta_max
    endif

  end subroutine s_nml_mesh_1d_unif_cart_print

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  subroutine init_clone_nml_mesh_1d_unif_cart( &
    self, &
    clone, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in) :: clone
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    character(len=256) :: err_msg
    character(len=256) :: caller

    caller = 'init_clone_nml_mesh_1d_unif_cart'
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

  end subroutine init_clone_nml_mesh_1d_unif_cart

  subroutine init_nml_mesh_1d_unif_cart( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    sll_int32 :: num_cells
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_unif_cart/ &
      num_cells, &
      eta_min, &
      eta_max
    caller = 'init_nml_mesh_1d_unif_cart'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      num_cells, &
      eta_min, &
      eta_max)

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

    read(namelist_id, mesh_1d_unif_cart)
    self%label="no_label"
    self%num_cells = num_cells
    self%eta_min = eta_min
    self%eta_max = eta_max
    close(namelist_id)

  end subroutine init_nml_mesh_1d_unif_cart

  subroutine init_nml_mesh_1d_unif_cart_1( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    sll_int32 :: num_cells_1
    sll_real64 :: eta_min_1
    sll_real64 :: eta_max_1
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_unif_cart_1/ &
      num_cells_1, &
      eta_min_1, &
      eta_max_1

    caller = 'init_nml_mesh_1d_unif_cart_1'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      num_cells_1, &
      eta_min_1, &
      eta_max_1)

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

    read(namelist_id, mesh_1d_unif_cart_1)
    self%label="_1"
    self%num_cells = num_cells_1
    self%eta_min = eta_min_1
    self%eta_max = eta_max_1
    close(namelist_id)

  end subroutine init_nml_mesh_1d_unif_cart_1

  subroutine init_nml_mesh_1d_unif_cart_2( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    sll_int32 :: num_cells_2
    sll_real64 :: eta_min_2
    sll_real64 :: eta_max_2
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_unif_cart_2/ &
      num_cells_2, &
      eta_min_2, &
      eta_max_2

    caller = 'init_nml_mesh_1d_unif_cart_2'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      num_cells_2, &
      eta_min_2, &
      eta_max_2)

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

    read(namelist_id, mesh_1d_unif_cart_2)
    self%label="_2"
    self%num_cells = num_cells_2
    self%eta_min = eta_min_2
    self%eta_max = eta_max_2
    close(namelist_id)

  end subroutine init_nml_mesh_1d_unif_cart_2

  subroutine init_nml_mesh_1d_unif_cart_3( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    sll_int32 :: num_cells_3
    sll_real64 :: eta_min_3
    sll_real64 :: eta_max_3
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_unif_cart_3/ &
      num_cells_3, &
      eta_min_3, &
      eta_max_3

    caller = 'init_nml_mesh_1d_unif_cart_3'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      num_cells_3, &
      eta_min_3, &
      eta_max_3)

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

    read(namelist_id, mesh_1d_unif_cart_3)
    self%label="_3"
    self%num_cells = num_cells_3
    self%eta_min = eta_min_3
    self%eta_max = eta_max_3
    close(namelist_id)

  end subroutine init_nml_mesh_1d_unif_cart_3

  subroutine init_nml_mesh_1d_unif_cart_4( &
    self, &
    filename, &
    proc_id)
    class(sll_t_nml_mesh_1d_unif_cart), intent(inout) :: self
    character(len=*), intent(in)    :: filename
    sll_int32, intent(in), optional :: proc_id

    sll_int32 :: namelist_id
    sll_int32 :: ierr
    sll_int32 :: IO_stat
    character(len=256) :: err_msg
    character(len=256) :: caller
    sll_int32 :: num_cells_4
    sll_real64 :: eta_min_4
    sll_real64 :: eta_max_4
    sll_int32 :: proc_id_loc

    namelist /mesh_1d_unif_cart_4/ &
      num_cells_4, &
      eta_min_4, &
      eta_max_4

    caller = 'init_nml_mesh_1d_unif_cart_4'
    if(present(proc_id))then
      proc_id_loc = proc_id
    else
      proc_id_loc = 0
    endif

    call set_default_values( &
      num_cells_4, &
      eta_min_4, &
      eta_max_4)

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

    read(namelist_id, mesh_1d_unif_cart_4)
    self%label="_4"
    self%num_cells = num_cells_4
    self%eta_min = eta_min_4
    self%eta_max = eta_max_4
    close(namelist_id)

  end subroutine init_nml_mesh_1d_unif_cart_4

  subroutine set_default_values( &
    num_cells, &
    eta_min, &
    eta_max)
    sll_int32, intent(inout) :: num_cells
    sll_real64, intent(inout) :: eta_min
    sll_real64, intent(inout) :: eta_max

    num_cells = 32
    eta_min = 0._f64
    eta_max = 1._f64

  end subroutine set_default_values

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_nml_mesh_1d_unif_cart
