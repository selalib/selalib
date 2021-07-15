module sll_m_fft_filter_3d
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d
  
   use sll_m_fft

  implicit none
  private

  public :: sll_t_fft_filter_3d
  


  type, extends(sll_c_filter_base_3d):: sll_t_fft_filter_3d

     sll_int32 :: k_min(3)
     sll_int32 :: k_max(3)

     type(sll_t_fft) :: fft(3)
     type(sll_t_fft) :: ifft(3)

   contains

     procedure :: init => init_fft_3d 
     procedure :: apply => apply_fft_3d
     procedure :: apply_inplace => apply_inplace_fft_3d
     
  end type sll_t_fft_filter_3d

contains

  subroutine init_fft_3d( self, iterations, n_dofs, mode )
    class( sll_t_fft_filter_3d), intent( out ) :: self
    sll_int32, intent( in ) :: iterations
    sll_int32, intent( in ) :: n_dofs(3)
    sll_int32, optional, intent( in ) :: mode(3)
    !local variables
    sll_comp64 :: array1d_x(n_dofs(1)),  array1d_y(n_dofs(2)),  array1d_z(n_dofs(3))

    self%iterations = iterations
    self%n_dofs = n_dofs
    self%k_min = max(1, mode-1 )
    self%k_max = min(n_dofs/2-1 ,mode+1)

    call sll_s_fft_init_c2c_1d( self%fft(1), n_dofs(1), array1d_x, array1d_x, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft(2), n_dofs(2), array1d_y, array1d_y, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%fft(3), n_dofs(3), array1d_z, array1d_z, sll_p_fft_forward)
    call sll_s_fft_init_c2c_1d( self%ifft(1), n_dofs(1), array1d_x, array1d_x, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft(2), n_dofs(2), array1d_y, array1d_y, &
         sll_p_fft_backward, normalized=.true.)
    call sll_s_fft_init_c2c_1d( self%ifft(3), n_dofs(3), array1d_z, array1d_z, &
         sll_p_fft_backward, normalized=.true.)
    
  end subroutine init_fft_3d

  subroutine apply_fft_3d(self, field_in, field_out )
    class( sll_t_fft_filter_3d), intent( inout ) :: self
    sll_real64,                  intent( in    ) :: field_in(:) !< array for the coefficients of the fields 
    sll_real64,                  intent(   out ) :: field_out(:) !< array for the coefficients of the fields
    !local variables
    sll_int32 :: i

    field_out = field_in
    do i = 1, self%iterations
       call sll_s_fft_filter( field_in, field_out, self%n_dofs, self%fft, self%ifft, self%k_min, self%k_max ) 
    end do
    
  end subroutine apply_fft_3d

   subroutine apply_inplace_fft_3d(self, field )
    class( sll_t_fft_filter_3d), intent( inout ) :: self
    sll_real64,                  intent( inout ) :: field(:) !< array for the coefficients of the fields 
     !local variables
    sll_int32 :: i

    do i = 1, self%iterations
       call sll_s_fft_filter( field, field, self%n_dofs, self%fft, self%ifft, self%k_min, self%k_max )
    end do
    
  end subroutine apply_inplace_fft_3d


  subroutine sll_s_fft_filter ( field_in, field_out, n_dofs, fft, ifft, k_min, k_max )
    sll_real64,                  intent( in    ) :: field_in(:)
    sll_real64,                  intent(   out ) :: field_out(:)
    sll_int32, intent( in ) :: n_dofs(3)
    type(sll_t_fft), intent( inout ) :: fft(3)
    type(sll_t_fft), intent( inout ) :: ifft(3)
    sll_int32, intent( in ) :: k_min(3)
    sll_int32, intent( in ) :: k_max(3)
    !local variables
    sll_int32 :: ind, i, j, k
    sll_comp64 :: scratch(n_dofs(1), n_dofs(2), n_dofs(3)), scratch2(n_dofs(1), n_dofs(2), n_dofs(3))
    sll_comp64 :: array1d_x(n_dofs(1)),  array1d_y(n_dofs(2)),  array1d_z(n_dofs(3))

    ! Fourier transform
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind = ind+1
             array1d_x(i) = cmplx( field_in(ind), 0_f64, kind=f64)
          end do
          call sll_s_fft_exec_c2c_1d( fft(1), array1d_x, array1d_x)
          do i=1,n_dofs(1)
             scratch(i,j,k) = array1d_x(i)
          end do
       end do
    end do
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( fft(2), array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch(i,j,k) = array1d_y(j)
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( fft(3), array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch(i,j,k) = array1d_z(k)
          end do
       end do
    end do

    !filtering
    scratch2 = 0._f64

    !filter kx
!!$    ! copy 0-mode
!!$    scratch2(1,:,:) = scratch(1,:,:)
!!$    scratch2(n_dofs(1)/2+1,:,:) = scratch(n_dofs(1)/2+1,:,:)
!!$    
!!$    ! copy from k_min to k_max 
!!$    do i= k_min(1), k_max(1)
!!$       scratch2(1+i,:,:) = scratch(1+i,:,:)
!!$       scratch2(n_dofs(1)+1-i,:,:) = scratch(n_dofs(1)+1-i,:,:)
!!$    end do
!!$  
!!$
!!$    !filter kx and ky
!!$    scratch2(1,1,:) = scratch(1,1,:)
!!$    scratch2(n_dofs(1)/2+1,1,:) = scratch(n_dofs(1)/2+1,1,:)
!!$    do j= k_min(1), k_max(1)
!!$       scratch2(1+j,1,:) = scratch(1+j,1,:)
!!$       scratch2(n_dofs(1)+1-j,1,:) = scratch(n_dofs(1)+1-j,1,:)
!!$    end do
!!$
!!$    scratch2(1,n_dofs(2)/2+1,:) = scratch(1,n_dofs(2)/2+1,:)
!!$    scratch2(n_dofs(1)/2+1,n_dofs(2)/2+1,:) = scratch(n_dofs(1)/2+1,n_dofs(2)/2+1,:)
!!$    do j= k_min(1), k_max(1)
!!$       scratch2(1+j,n_dofs(2)/2+1,:) = scratch(1+j,n_dofs(2)/2+1,:)
!!$       scratch2(n_dofs(1)+1-j,n_dofs(2)/2+1,:) = scratch(n_dofs(1)+1-j,n_dofs(2)/2+1,:)
!!$    end do
!!$
!!$    do k= k_min(2), k_max(2)
!!$       scratch2(1,1+k,:) = scratch(1,1+k,:)
!!$       scratch2(n_dofs(1)/2+1,1+k,:) = scratch(n_dofs(1)/2+1,1+k,:)
!!$       do j= k_min(1), k_max(1)
!!$          scratch2(1+j,1+k,:) = scratch(1+j,1+k,:)
!!$          scratch2(n_dofs(1)+1-j,1+k,:) = scratch(n_dofs(1)+1-j,1+k,:)
!!$       end do
!!$
!!$       scratch2(1,n_dofs(2)+1-k,:) = scratch(1,n_dofs(2)+1-k,:)
!!$       scratch2(n_dofs(1)/2+1,n_dofs(2)+1-k,:) = scratch(n_dofs(1)/2+1,n_dofs(2)+1-k,:)
!!$       do j= k_min(1), k_max(1)
!!$          scratch2(1+j,n_dofs(2)+1-k,:) = scratch(1+j,n_dofs(2)+1-k,:)
!!$          scratch2(n_dofs(1)+1-j,n_dofs(2)+1-k,:) = scratch(n_dofs(1)+1-j,n_dofs(2)+1-k,:)
!!$       end do
!!$    end do
!!$
!!$    
    !filter ky and kz
    scratch2(:,1,1) = scratch(:,1,1)
    scratch2(:,n_dofs(2)/2+1,1) = scratch(:,n_dofs(2)/2+1,1)
    do j= k_min(2), k_max(2)
       scratch2(:,1+j,1) = scratch(:,1+j,1)
       scratch2(:,n_dofs(2)+1-j,1) = scratch(:,n_dofs(2)+1-j,1)
    end do

    scratch2(:,1,n_dofs(3)/2+1) = scratch(:,1,n_dofs(3)/2+1)
    scratch2(:,n_dofs(2)/2+1,n_dofs(3)/2+1) = scratch(:,n_dofs(2)/2+1,n_dofs(3)/2+1)
    do j= k_min(2), k_max(2)
       scratch2(:,1+j,n_dofs(3)/2+1) = scratch(:,1+j,n_dofs(3)/2+1)
       scratch2(:,n_dofs(2)+1-j,n_dofs(3)/2+1) = scratch(:,n_dofs(2)+1-j,n_dofs(3)/2+1)
    end do

    do k= k_min(3), k_max(3)
       scratch2(:,1,1+k) = scratch(:,1,1+k)
       scratch2(:,n_dofs(2)/2+1,1+k) = scratch(:,n_dofs(2)/2+1,1+k)
       do j= k_min(2), k_max(2)
          scratch2(:,1+j,1+k) = scratch(:,1+j,1+k)
          scratch2(:,n_dofs(2)+1-j,1+k) = scratch(:,n_dofs(2)+1-j,1+k)
       end do

       scratch2(:,1,n_dofs(3)+1-k) = scratch(:,1,n_dofs(3)+1-k)
       scratch2(:,n_dofs(2)/2+1,n_dofs(3)+1-k) = scratch(:,n_dofs(2)/2+1,n_dofs(3)+1-k)
       do j= k_min(2), k_max(2)
          scratch2(:,1+j,n_dofs(3)+1-k) = scratch(:,1+j,n_dofs(3)+1-k)
          scratch2(:,n_dofs(2)+1-j,n_dofs(3)+1-k) = scratch(:,n_dofs(2)+1-j,n_dofs(3)+1-k)
       end do
    end do
    
    
    ! Inverse Fourier transform
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          do k=1,n_dofs(3)
             array1d_z(k) = scratch2(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( ifft(3), array1d_z, array1d_z)
          do k=1,n_dofs(3)
             scratch2(i,j,k) = array1d_z(k)
          end do
       end do
    end do

    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          do j=1,n_dofs(2)
             array1d_y(j) = scratch2(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( ifft(2), array1d_y, array1d_y)
          do j=1,n_dofs(2)
             scratch2(i,j,k) = array1d_y(j)
          end do
       end do
    end do
    
    ind=0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             array1d_x(i) = scratch2(i,j,k)
          end do
          call sll_s_fft_exec_c2c_1d( ifft(1), array1d_x, array1d_x)
          
          do i=1,n_dofs(1)
             ind = ind+1
             field_out(ind) = real( array1d_x(i), kind=f64 )
          end do
       end do
    end do
  end subroutine sll_s_fft_filter

end module sll_m_fft_filter_3d
