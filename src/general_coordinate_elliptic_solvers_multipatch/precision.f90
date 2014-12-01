Module used_precision
  Implicit None
 
  Integer, Parameter :: sp=Kind(1.0), dp=Kind(1.0d0)

  Integer, Parameter :: wp = dp
  
  Real(wp),Parameter :: zero = 0.0_wp
  Real(wp),Parameter :: one = 1.0_wp
  Real(wp),Parameter :: two = 2.0_wp
  Real(wp),Parameter :: pi = 3.141592653589793_wp
  Real(wp),Parameter :: epsilon = 1.e-14
  logical ,Parameter :: per = .true. 
  logical ,Parameter :: nat = .false. 
contains
	! IT SEEMS THAT gfortran DOES NOT LIKE TO EXPORT A MODULE WITHOUT A SUBROUTINE
  subroutine print_one ()
    implicit none
    print*,'1'
  end subroutine print_one
End Module used_precision

