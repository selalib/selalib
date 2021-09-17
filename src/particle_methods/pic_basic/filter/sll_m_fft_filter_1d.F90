module sll_m_fft_filter_1d
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d
  
  use sll_m_fft

  implicit none
  private

  public :: sll_t_fft_filter_1d



  type, extends(sll_c_filter_base_1d) :: sll_t_fft_filter_1d

     sll_int32 :: k_min
     sll_int32 :: k_max
     type(sll_t_fft) :: fft
     type(sll_t_fft) :: ifft


   contains

     procedure :: init => init_fft_1d 
     procedure :: apply => apply_fft_1d
     procedure :: apply_inplace => apply_inplace_fft_1d

  end type sll_t_fft_filter_1d

contains

  subroutine init_fft_1d ( self, iterations, n_dofs, mode )
    class( sll_t_fft_filter_1d), intent( out ) :: self
    sll_int32, intent( in ) :: iterations
    sll_int32, intent( in ) :: n_dofs
    sll_int32, optional, intent( in ) :: mode
    !local variables
    sll_real64:: array1(n_dofs),  array2(n_dofs)

    self%iterations = iterations
    self%n_dofs = n_dofs
    self%k_min = max(1, mode-1 )
    self%k_max = min(n_dofs/2-1, mode+1)
    
    
    call sll_s_fft_init_r2r_1d( self%fft, n_dofs, array1, array2, sll_p_fft_forward, normalized=.false.)
    call sll_s_fft_init_r2r_1d( self%ifft, n_dofs, array1, array1, &
         sll_p_fft_backward, normalized=.false.)

  end subroutine init_fft_1d 

  subroutine apply_fft_1d (self, field_in, field_out )
    class( sll_t_fft_filter_1d), intent( inout ) :: self
    sll_real64,                  intent( in    ) :: field_in(:) !< array for the coefficients of the fields 
    sll_real64,                  intent(   out ) :: field_out(:) !< array for the coefficients of the fields
    !local variables
    sll_int32 :: i

    field_out = field_in
    do i = 1, self%iterations
       call sll_s_fft_filter( field_in, field_out, self%n_dofs, self%fft, self%ifft, self%k_min, self%k_max )
    end do

  end subroutine apply_fft_1d 
  
  subroutine apply_inplace_fft_1d (self, field )
    class( sll_t_fft_filter_1d), intent( inout ) :: self
    sll_real64,                  intent( inout ) :: field(:)
    !local variables
    sll_int32 :: i

    do i = 1, self%iterations
       call sll_s_fft_filter( field, field, self%n_dofs, self%fft, self%ifft, self%k_min, self%k_max )
    end do
    
  end subroutine apply_inplace_fft_1d


  subroutine sll_s_fft_filter( field_in, field_out, n_dofs, fft, ifft, k_min, k_max )
    sll_real64,                  intent( in    ) :: field_in(:)
    sll_real64,                  intent(   out ) :: field_out(:)
    sll_int32, intent( in ) :: n_dofs
    type(sll_t_fft), intent( inout ) :: fft
    type(sll_t_fft), intent( inout ) :: ifft
    sll_int32, intent( in ) :: k_min
    sll_int32, intent( in ) :: k_max
    !local variables
    sll_int32 :: i
    sll_real64:: array1(n_dofs), array2(n_dofs)
    
    array1 = field_in
    
    ! Fourier transform
    call sll_s_fft_exec_r2r_1d ( fft, array1, array2 )
    
    ! Filter
    array1 = 0._f64
    !0 mode
    array1(1) = array2(1)
    array1(n_dofs/2+1) = array2(n_dofs/2+1)
    ! 1 mode
    do i = k_min, k_max
       array1(i+1) = array2(i+1)
       array1(n_dofs+1-i) = array2(n_dofs+1-i)
    end do

    ! Inverse Fourier transform
    call sll_s_fft_exec_r2r_1d( ifft, array1 , field_out )
    ! normalize
    field_out = field_out/ real(n_dofs, f64)
    
  end subroutine sll_s_fft_filter
  
end module sll_m_fft_filter_1d
