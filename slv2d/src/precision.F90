MODULE used_precision
  !=========================================
  !    
  !    File:          mabo2d.F90
  !    Project:       vlasov
  !    Author(s):     Eric Sonnendrucker
  !    Creation:      03.02.1998
  !    Last modified: 10.03.1999 for Cray T3E
  !    
  !=========================================

  implicit none

  Integer, Parameter :: wp = 8

  REAL(wp),PARAMETER :: zero = 0.0_wp
  REAL(wp),PARAMETER :: one = 1.0_wp
  REAL(wp),PARAMETER :: two = 2.0_wp
  REAL(wp),PARAMETER :: pi = 3.141592653589793_wp
  REAL(wp),PARAMETER :: epsilon = 1.e-14

END MODULE used_precision
