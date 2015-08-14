!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations in 1D
!> @details
!> Contains the abstract class to create a Maxwell solver in 1D.

module sll_m_maxwell_1d_base
#include "sll_working_precision.h"
#include "sll_utilities.h"

  implicit none
  private
  public :: sll_plot_two_fields_1d
  
  type, public, abstract :: sll_maxwell_1d_base

   contains
     procedure(compute_field1_from_field2), deferred :: &
          compute_E_from_B !< Solve E and B part of Amperes law with B constant in time
     procedure(compute_field1_from_field2), deferred :: &
          compute_B_from_E !< Solve Faraday equation with E constant in time
     procedure(signature_compute_E_from_rho_1d), deferred :: &
          compute_E_from_rho !< Solve E from rho using Poisson
     !procedure(signature_solve), deferred :: &
     !     solve !< Solve Amperes law and Faraday equation
  end type sll_maxwell_1d_base

  abstract interface 
     subroutine compute_field1_from_field2(this, delta_t, field_in, field_out)
     use sll_working_precision
     import sll_maxwell_1d_base
     
     class(sll_maxwell_1d_base) :: this
     sll_real64, intent(in)     :: delta_t
     sll_real64, intent(in)     :: field_in(:)
     sll_real64, intent(inout)  :: field_out(:)
   end subroutine compute_field1_from_field2
  end interface

  abstract interface    
    subroutine signature_compute_E_from_rho_1d(this, E, rho )
      use sll_working_precision
      import sll_maxwell_1d_base       
      class(sll_maxwell_1d_base) :: this
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_1d
  end interface

contains
  !> write files to visualize 1d fields with gnuplot
  subroutine sll_plot_two_fields_1d(fname, n1, f1, f2, iplot, time )
    character(len=*),             intent(in) :: fname !< output file name
    sll_int32,                    intent(in) :: n1    !< size of f1 and f2 
    sll_real64, dimension(n1),    intent(in) :: f1    !< first field 2d
    sll_real64, dimension(n1),    intent(in) :: f2    !< second field 2d
    sll_int32,                    intent(in) :: iplot !< plot counter
    sll_real64,                   intent(in) :: time  !< step time

    integer          :: i, j
    character(len=4) :: cplot

    call int2string(iplot, cplot)

    !write domains
    open( 80, file = fname//cplot//".dat" )
    do i=1,n1
          write(80,*) i, sngl(f1(i)), sngl(f2(i))
    end do
    close(80)

    open( 90, file = fname//'plots.gnu', position="append" )
    if ( iplot == 1 ) then
       rewind(90)
       !write(90,*)"set xr[-0.1:1.1]"
       !write(90,*)"set yr[-0.1:1.1]"
       !write(90,*)"set zr[-1.1:1.1]"
       !write(90,*)"set cbrange[-1:1]"
       !write(90,*)"set pm3d"
       !write(90,*)"set surf"
       write(90,*)"set term x11"
    end if
    write(90,*)"set title 'Time = ",time,"'"
    write(90,"(a)",advance='no')"plot '"//fname//cplot//".dat' w lines"
    write(90,"(a)")",'"//fname//cplot//".dat' u 1:3 w lines"
    write(90,*)"pause -1"
    close(90)

  end subroutine sll_plot_two_fields_1d

end module sll_m_maxwell_1d_base
