!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann
!> @brief Binomial filter for smooting of fields
!> @details ...

module sll_m_binomial_filter

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d
  
  implicit none

  public :: &
       sll_t_binomial_filter

  private
  
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type, extends(sll_c_filter_base_1d) :: sll_t_binomial_filter

     sll_real64, allocatable :: scratch(:)

     logical :: odd
     sll_int32 :: double_iters
     
   contains

     procedure :: init => init_binomial
     procedure :: apply => apply_binomial
     procedure :: apply_inplace => apply_inplace_binomial

  end type sll_t_binomial_filter
  
contains

  subroutine init_binomial( self, iterations, n_dofs, mode  )
    class( sll_t_binomial_filter ), intent( out ) :: self
    sll_int32, intent( in ) :: iterations
    sll_int32, intent( in ) :: n_dofs
    sll_int32, optional, intent( in ) :: mode
    
    self%iterations = iterations
    self%n_dofs = n_dofs
    allocate( self%scratch(n_dofs) )

    if (mod(iterations, 2) .eq. 1 ) then
       self%odd = .true.
       self%double_iters = iterations/2
    else
       self%odd = .false.
       self%double_iters = iterations/2
    end if

    !print*, 'Binomial filter', self%double_iters, self%odd
    
  end subroutine init_binomial

  subroutine apply_binomial( self, field_in, field_out ) 
    class( sll_t_binomial_filter ), intent( inout ) :: self
    sll_real64,                            intent(in)  :: field_in(:) !< array for the coefficients of the fields 
    sll_real64,                            intent(out)  :: field_out(:) !< array for the coefficients of the fields

    sll_int32 :: iter

    if ( self%odd .eqv. .true. ) then
       call sll_s_binomial_filter( field_in, field_out, self%n_dofs )
    else
       field_out = field_in
    end if
       
    
    do iter=1,self%double_iters
       call sll_s_binomial_filter( field_out, self%scratch, self%n_dofs )
       call sll_s_binomial_filter( self%scratch, field_out, self%n_dofs )
    end do
  end subroutine apply_binomial
    
    
  subroutine apply_inplace_binomial( self, field )
    class( sll_t_binomial_filter ), intent( inout ) :: self
    sll_real64,                            intent(inout)  :: field(:) !< array for the coefficients of the efields 

    sll_int32 :: iter

    do iter=1,self%double_iters
       call sll_s_binomial_filter( field, self%scratch, self%n_dofs )
       call sll_s_binomial_filter( self%scratch, field, self%n_dofs )
    end do

    if ( self%odd .eqv. .true. ) then
       call sll_s_binomial_filter( field, self%scratch, self%n_dofs )
       field = self%scratch
    end if
    
  end subroutine apply_inplace_binomial
  
  subroutine sll_s_binomial_filter( field_in, field_out, n_dofs )
    sll_real64,                            intent(in)  :: field_in(:) !< array for the coefficients of the efields 
    sll_real64,                            intent(out)  :: field_out(:) !< array for the coefficients of the efields
    sll_int32, intent(in) :: n_dofs

    sll_int32 :: j
    
    ! Filter
    field_out(1) = 0.25_f64*( field_in(1)*2.0_f64 + &
               field_in(n_dofs)+field_in(2))
    do j=2,n_dofs-1
       field_out(j) = 0.25_f64*( field_in(j)*2.0_f64 + &
            field_in(j-1)+field_in(j+1))
    end do
    field_out(n_dofs) = 0.25_f64*( field_in(n_dofs)*2.0_f64 + &
         field_in(n_dofs-1)+field_in(1))

  end subroutine sll_s_binomial_filter

  

  

end module sll_m_binomial_filter
