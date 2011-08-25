module sll_timer
  use iso_c_binding
  implicit none

  interface 
     type(c_ptr) function set_time_mark() bind(c, name='set_time_mark_C')
       use iso_c_binding
     end function set_time_mark
  end interface

  ! Note that the 'value' attribute is essential for this interface to work.
  ! If not present, Fortran will take the address of the pointers passed and
  ! a mess will be unleashed...
  interface
     type(c_ptr) function reset_time_mark( mark ) &
          bind(c, name='reset_time_mark_C')
       use iso_c_binding
       type(c_ptr), value :: mark
     end function reset_time_mark
  end interface

  interface
     real(c_double) function time_elapsed_since( t0 ) &
          bind( c, name='time_elapsed_since_C')
       use iso_c_binding
       type(c_ptr), intent(in), value :: t0
     end function time_elapsed_since
  end interface

end module sll_timer
