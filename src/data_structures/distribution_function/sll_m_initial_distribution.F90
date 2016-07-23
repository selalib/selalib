module sll_m_initial_distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only : &
       sll_p_twopi

  implicit none
  
  public :: sll_p_sumcos_onegaussian, &
       sll_p_cossum_onegaussian, &
       sll_c_distribution_params, &
       sll_t_cos_gaussian, &
       sll_s_initial_distribution_new
  
  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Descriptors for various distributions
  sll_int32, parameter :: sll_p_sumcos_onegaussian = 0
  sll_int32, parameter :: sll_p_cossum_onegaussian = 1
  sll_int32, parameter :: sll_p_sumcos_twogaussian = 2
  sll_int32, parameter :: sll_p_cossum_twogaussian = 3

  !> Abstract data type for parameters of initial distribution
  type, abstract :: sll_c_distribution_params
     sll_int32 :: dims(2) !< Number of spatial and velocity dimensions

   contains
     procedure( signature_eval ), deferred :: eval
     procedure( signature_evalx), deferred :: evalx
     procedure( signature_evalv), deferred :: evalv
     procedure( signature_empty), deferred :: free

  end type sll_c_distribution_params

  !> 
  type, extends(sll_c_distribution_params) :: sll_t_cos_gaussian
     sll_real64, allocatable :: kx(:,:)
     sll_real64, allocatable :: alpha(:)
     sll_real64, allocatable :: v_thermal(:,:)
     sll_real64, allocatable :: v_mean(:,:)
     sll_real64, allocatable :: delta(:)
     sll_real64, allocatable :: normal(:)
     sll_int32               :: n_gaussians
     sll_int32               :: n_cos

   contains
     procedure :: eval => sll_f_cos_gaussian
     procedure :: free => free_cos_gaussian
     procedure :: evalx => sll_f_cos
     procedure :: evalv => sll_f_gaussian
     procedure :: init => cos_gaussian_init
     
  end type sll_t_cos_gaussian

  abstract interface
     subroutine signature_empty( self )
       import sll_c_distribution_params
       class( sll_c_distribution_params ), intent(inout) :: self
       
     end subroutine signature_empty
  end interface
  
  abstract interface
     function signature_eval( self, x, v ) result( fval )
       use sll_m_working_precision
       import sll_c_distribution_params
       class( sll_c_distribution_params ) :: self
       sll_real64 :: x(:)
       sll_real64 :: v(:)
       sll_real64 :: fval
       
     end function signature_eval
  end interface
  
  abstract interface
     function signature_evalx( self, x ) result( fval )
       use sll_m_working_precision
       import sll_c_distribution_params
       class( sll_c_distribution_params ) :: self
       sll_real64 :: x(:)
       sll_real64 :: fval
       
     end function signature_evalx
  end interface

  abstract interface
     function signature_evalv( self, v ) result( fval )
       use sll_m_working_precision
       import sll_c_distribution_params
       class( sll_c_distribution_params ) :: self
       sll_real64 :: v(:)
       sll_real64 :: fval
       
     end function signature_evalv
  end interface
contains

  function sll_f_cos_gaussian( self, x, v ) result( fval )
    class( sll_t_cos_gaussian ) :: self
    sll_real64 :: x(:)
    sll_real64 :: v(:)
    sll_real64 :: fval

    sll_real64 :: fexp
    sll_int32  :: j

    fval = 1.0_f64
    do j=1,self%n_cos
       fval = fval + self%alpha(j) * cos( sum(self%kx(:,j) * x) )
    end do

    fexp = 0.0_f64
    do j=1,self%n_gaussians
       fexp = fexp + self%normal(j)*self%delta(j)* &
            exp( -0.5_f64 * sum( ((v-self%v_mean(:,j))/self%v_thermal(:,j))**2 ) )
    end do

    fval = fval*fexp

  end function sll_f_cos_gaussian

  function sll_f_cos( self, x ) result( fval )
    class( sll_t_cos_gaussian ) :: self
    sll_real64 :: x(:)
    sll_real64 :: fval

    sll_int32  :: j

    fval = 1.0_f64
    do j=1,self%n_cos
       fval = fval + self%alpha(j) * cos( sum(self%kx(:,j) * x) )
    end do

  end function sll_f_cos
  
  function sll_f_gaussian( self, v ) result( fval )
    class( sll_t_cos_gaussian ) :: self
    sll_real64 :: v(:)
    sll_real64 :: fval

    sll_int32  :: j

    fval = 0.0_f64
    do j=1,self%n_gaussians
       fval = fval + self%normal(j)*self%delta(j)* &
            exp( -0.5_f64 * sum( ((v-self%v_mean(:,j))/self%v_thermal(:,j))**2 ) )
    end do

  end function sll_f_gaussian
  
  subroutine free_cos_gaussian( self )
    class( sll_t_cos_gaussian ), intent( inout ) :: self


    if (allocated(self%kx)) deallocate(self%kx)
    if (allocated(self%alpha)) deallocate(self%alpha)   
    if (allocated(self%v_thermal)) deallocate(self%v_thermal)
    if (allocated(self%v_mean)) deallocate(self%v_mean)
    if (allocated(self%normal)) deallocate(self%normal)
    if (allocated(self%delta)) deallocate(self%delta)

  end subroutine free_cos_gaussian


  subroutine cos_gaussian_init( self, descriptor, dims, file_id )
    class( sll_t_cos_gaussian ), intent( out ) :: self
    sll_int32, intent( in    ) :: descriptor !< descriptor of the test case
    sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
    sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format

    self%dims = dims
    select case( descriptor )
    case( sll_p_sumcos_onegaussian )
       call sumcos_onegaussian_init( file_id, self )
       
    case( sll_p_cossum_onegaussian )
       call cossum_onegaussian_init( file_id, self )
       
    case( sll_p_cossum_twogaussian )
       call cossum_twogaussian_init( file_id, self )
       
    case( sll_p_sumcos_twogaussian )
       call sumcos_twogaussian_init( file_id, self )
       
    end select
    

  end subroutine cos_gaussian_init


  subroutine  sll_s_initial_distribution_new( distribution, dims, file_id, params )
    sll_int32, intent( in    ) :: distribution !< descriptor of the test case
    sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
    sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
    class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.
    
    
    params%dims = dims
    select case( distribution )
    case( sll_p_sumcos_onegaussian )
       allocate( sll_t_cos_gaussian :: params )
       select type( params )
       type is( sll_t_cos_gaussian )
          call sumcos_onegaussian_init( file_id, params )
       end select
    case( sll_p_cossum_onegaussian )
       allocate( sll_t_cos_gaussian :: params )
       select type( params )
       type is( sll_t_cos_gaussian )
          call cossum_onegaussian_init( file_id, params )
       end select
    case( sll_p_cossum_twogaussian )
       allocate( sll_t_cos_gaussian :: params )
       select type( params )
       type is( sll_t_cos_gaussian )
          call cossum_twogaussian_init( file_id, params )
       end select
    case( sll_p_sumcos_twogaussian )
       allocate( sll_t_cos_gaussian :: params )
       select type( params )
       type is( sll_t_cos_gaussian )
          call sumcos_twogaussian_init( file_id, params )
       end select
    end select
    
    
  end subroutine sll_s_initial_distribution_new


  subroutine sumcos_onegaussian_init( file_id, params )
    sll_int32, intent( in ) :: file_id
    type( sll_t_cos_gaussian ), intent( inout ) :: params

    ! local variables
    sll_real64 :: kx(params%dims(1))
    sll_real64 :: alpha(params%dims(1))
    sll_real64 :: v_thermal(params%dims(2))
    sll_real64 :: v_mean(params%dims(2))
    sll_int32  :: j
    
    namelist /sumcos_onegaussian/ kx, alpha, v_thermal, v_mean
    
    read(file_id, sumcos_onegaussian)

    allocate( params%kx(params%dims(1),params%dims(1)) )
    allocate( params%alpha(params%dims(1)) )
    allocate( params%v_thermal(params%dims(2),1) )
    allocate( params%v_mean(params%dims(2),1) )
    allocate( params%normal(1) )
    allocate( params%delta(1) )

    params%n_cos = params%dims(1)
    params%kx = 0.0_f64
    do j=1, params%dims(1)
       params%kx(:,j) = kx(j)
    end do
    params%alpha = alpha
    params%v_thermal(:,1) = v_thermal
    params%v_mean(:,1) = v_mean

    
    params%normal = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
         product(params%v_thermal(:,1)))

    params%delta(1) = 0.0_f64
    
  end subroutine sumcos_onegaussian_init



  subroutine cossum_onegaussian_init( file_id, params )
    sll_int32, intent( in ) :: file_id
    type( sll_t_cos_gaussian ), intent( inout ) :: params

    ! local variables
    sll_real64 :: kx(params%dims(1))
    sll_real64 :: alpha
    sll_real64 :: v_thermal(params%dims(2))
    sll_real64 :: v_mean(params%dims(2))
    
    namelist /cossum_onegaussian/ kx, alpha, v_thermal, v_mean
    
    read(file_id, cossum_onegaussian)
    
    allocate( params%kx(params%dims(1),1) )
    allocate( params%alpha( 1 ) )
    allocate( params%v_thermal(params%dims(2),1) )
    allocate( params%v_mean(params%dims(2),1) )
    allocate( params%normal(1) )
    allocate( params%delta(1) )

    params%n_cos = 1
    params%kx(:,1) = kx
    params%alpha = alpha
    params%v_thermal(:,1) = v_thermal
    params%v_mean (:,1) = v_mean

    params%normal = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
         product(params%v_thermal(:,1)))
    params%delta(1) = 0.0_f64
    
  end subroutine cossum_onegaussian_init
  
  
  subroutine cossum_twogaussian_init( file_id, params )
    sll_int32, intent( in ) :: file_id
    type( sll_t_cos_gaussian ), intent( inout ) :: params

    ! local variables
    sll_real64 :: kx(params%dims(1))
    sll_real64 :: alpha
    sll_real64 :: v_thermal_1(params%dims(2))
    sll_real64 :: v_mean_1(params%dims(2))
    sll_real64 :: v_thermal_2(params%dims(2))
    sll_real64 :: v_mean_2(params%dims(2))
    sll_real64 :: delta

    sll_int32  :: j
    
    namelist /cossum_twogaussian/ kx, alpha, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta
    
    read(file_id, cossum_twogaussian)

    allocate( params%kx(params%dims(1),1) )
    allocate( params%v_thermal(params%dims(2),2) )
    allocate( params%v_mean(params%dims(2),2) )
    allocate( params%normal(2) )
    allocate( params%delta(2) )

    params%n_cos = 1
    params%kx(:,1) = kx
    params%alpha = alpha
    params%v_thermal(:,1) = v_thermal_1
    params%v_mean(:,1) = v_mean_1
    params%v_thermal(:,2) = v_thermal_2
    params%v_mean(:,2) = v_mean_2
    params%delta(1) = delta
    params%delta(2) = 1.0_f64 - delta

    do j=1,2
       params%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
            product(params%v_thermal(:,j)))

    end do
    
  end subroutine cossum_twogaussian_init

   subroutine sumcos_twogaussian_init( file_id, params )
    sll_int32, intent( in ) :: file_id
    type( sll_t_cos_gaussian ), intent( inout ) :: params

    ! local variables
    sll_real64 :: kx(params%dims(1))
    sll_real64 :: alpha(params%dims(1))
    sll_real64 :: v_thermal_1(params%dims(2))
    sll_real64 :: v_mean_1(params%dims(2))
    sll_real64 :: v_thermal_2(params%dims(2))
    sll_real64 :: v_mean_2(params%dims(2))
    sll_real64 :: delta

    sll_int32  :: j
    
    namelist /cossum_twogaussian/ kx, alpha, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta
    
    read(file_id, cossum_twogaussian)

    allocate( params%kx(params%dims(1),params%dims(1)) )
    allocate( params%v_thermal(params%dims(2),2) )
    allocate( params%v_mean(params%dims(2),2) )
    allocate( params%normal(2) )
    allocate( params%delta(2) )

    params%n_cos = params%dims(1)
    params%kx = 0.0_f64
    do j=1, params%dims(1)
       params%kx(:,j) = kx(j)
    end do
    params%alpha = alpha
    params%v_thermal(:,1) = v_thermal_1
    params%v_mean(:,1) = v_mean_1
    params%v_thermal(:,2) = v_thermal_2
    params%v_mean(:,2) = v_mean_2
    params%delta(1) = delta
    params%delta(2) = 1.0_f64 - delta

    do j=1,2
       params%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
            product(params%v_thermal(:,j)))

    end do
    
  end subroutine sumcos_twogaussian_init
  

end module sll_m_initial_distribution
