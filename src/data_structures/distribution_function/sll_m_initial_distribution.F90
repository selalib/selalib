!> @ingroup distribution_function
!> @author Katharina Kormann
!> @brief Parameters to define common initial distributions
!> @details ...
module sll_m_initial_distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_constants, only : &
       sll_p_twopi

  implicit none
  
  public :: sll_p_sumcos_onegaussian, &
       sll_p_cossum_onegaussian, &
       sll_p_sumcos_twogaussian, &
       sll_p_cossum_twogaussian, &
       sll_c_distribution_params, &
       sll_t_params_cos_gaussian, &
       sll_s_initial_distribution_new, &
       sll_s_initial_distribution_new_descriptor
  
  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Descriptors for various distributions
  sll_int32, parameter :: sll_p_sumcos_onegaussian = 0 !< Descriptor for (1+\sum cos( kx * x_i))*exp(-0.5(v-v_mean)**2/v_thermal**2)
  sll_int32, parameter :: sll_p_cossum_onegaussian = 1 !< Descriptor for (1+cos( \sum kx_i * x_i))*exp(-0.5(v-v_mean)**2/v_thermal**2)
  sll_int32, parameter :: sll_p_sumcos_twogaussian = 2 !< as sll_p_sumcos_onegaussian but with sum of two Gaussians
  sll_int32, parameter :: sll_p_cossum_twogaussian = 3 !< as sll_p_sumcos_onegaussian but with sum of two Gaussians

  !> Abstract data type for parameters of initial distribution
  type, abstract :: sll_c_distribution_params
     sll_int32 :: dims(2) !< Number of spatial and velocity dimensions

   contains
     procedure( signature_eval ), deferred :: eval_xv_density   !< Evaluate the distribution function
     procedure( signature_evalx), deferred :: eval_x_density  !< Evaluate the charge density
     procedure( signature_evalv), deferred :: eval_v_density  !< Evaluate the v-dependence (integrated over x)
     procedure( signature_empty), deferred :: free   !< Destructor

  end type sll_c_distribution_params

  !> Data type for distribution function with (multiple) Gaussians in v and one plus cosine perturbations in x.
  type, extends(sll_c_distribution_params) :: sll_t_params_cos_gaussian
     sll_real64, allocatable :: kx(:,:)  !< values of the wave numbers (first index dimension, second index for multiple cosines)
     sll_real64, allocatable :: alpha(:) !< strength of perturbations
     sll_real64, allocatable :: v_thermal(:,:) !< variance of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
     sll_real64, allocatable :: v_mean(:,:)    !< mean value of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
     sll_real64, allocatable :: delta(:) !< Portion of each Gaussian
     sll_real64, allocatable :: normal(:) !< Normalization constant of each Gaussian
     sll_int32               :: n_gaussians !< Number of Gaussians
     sll_int32               :: n_cos !< Number of cosines

   contains
     procedure :: eval_xv_density => sll_f_cos_gaussian !< Evaluate the distribution function
     procedure :: free => free_cos_gaussian  !< Descructor
     procedure :: eval_x_density => sll_f_cos         !< Evaluate the charge density
     procedure :: eval_v_density => sll_f_gaussian    !< Evaluate the v-dependence (integrated over x)
     procedure :: init => cos_gaussian_init  !< Initialization
     
  end type sll_t_params_cos_gaussian

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

!------------------------------------------------------------------------------- 
! Define the procedures of the type sll_t_params_cos_gaussian
!-------------------------------------------------------------------------------
  
  function sll_f_cos_gaussian( self, x, v ) result( fval )
    class( sll_t_params_cos_gaussian ) :: self
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
    class( sll_t_params_cos_gaussian ) :: self
    sll_real64 :: x(:)
    sll_real64 :: fval

    sll_int32  :: j

    fval = 1.0_f64
    do j=1,self%n_cos
       fval = fval + self%alpha(j) * cos( sum(self%kx(:,j) * x) )
    end do

  end function sll_f_cos
  
  function sll_f_gaussian( self, v ) result( fval )
    class( sll_t_params_cos_gaussian ) :: self
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
    class( sll_t_params_cos_gaussian ), intent( inout ) :: self


    if (allocated(self%kx)) deallocate(self%kx)
    if (allocated(self%alpha)) deallocate(self%alpha)   
    if (allocated(self%v_thermal)) deallocate(self%v_thermal)
    if (allocated(self%v_mean)) deallocate(self%v_mean)
    if (allocated(self%normal)) deallocate(self%normal)
    if (allocated(self%delta)) deallocate(self%delta)

  end subroutine free_cos_gaussian


  subroutine cos_gaussian_init( self, descriptor, dims, file_id )
    class( sll_t_params_cos_gaussian ), intent( out ) :: self
    sll_int32, intent( in    ) :: descriptor !< descriptor of the test case
    sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
    sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format

    self%dims = dims
    select case( descriptor )
    case( sll_p_sumcos_onegaussian )
       select case( self%dims(1) )
       case (1)
          select case (self%dims(2) )
          case(1)
             call sumcos_onegaussian_init_1d1v( file_id, self )
          case(2)
             call sumcos_onegaussian_init_1d2v( file_id, self )
          end select
       case(2)
          call sumcos_onegaussian_init_2d2v( file_id, self )
       case(3)
          call sumcos_onegaussian_init_3d3v( file_id, self )
       end select
    case( sll_p_cossum_onegaussian )
       select case( self%dims(1) )
       case (1)
          select case (self%dims(2) )
          case(1)
             call cossum_onegaussian_init_1d1v( file_id, self )
          case(2)
             call cossum_onegaussian_init_1d2v( file_id, self )
          end select
       case(2)
          call cossum_onegaussian_init_2d2v( file_id, self )
       case(3)
          call cossum_onegaussian_init_3d3v( file_id, self )
       end select
       
    case( sll_p_cossum_twogaussian )
       select case( self%dims(1) )
       case (1)
          select case (self%dims(2) )
          case(1)
             call cossum_twogaussian_init_1d1v( file_id, self )
          case(2)
             call cossum_twogaussian_init_1d2v( file_id, self )
          end select
       case(2)
          call cossum_twogaussian_init_2d2v( file_id, self )
       case(3)
          call cossum_twogaussian_init_3d3v( file_id, self )
       end select
    case( sll_p_sumcos_twogaussian )
       select case( self%dims(1) )
       case (1)
          select case (self%dims(2) )
          case(1)
             call sumcos_twogaussian_init_1d1v( file_id, self )
          case(2)
             call sumcos_twogaussian_init_1d2v( file_id, self )
          end select
       case(2)
          call sumcos_twogaussian_init_2d2v( file_id, self )
       case(3)
          call sumcos_twogaussian_init_3d3v( file_id, self )
       end select
    end select
    

  end subroutine cos_gaussian_init


  
  !------------------------------------------------------------------------------- 
  ! Factory function with specific functions called depending on chosen distribution type
  !-------------------------------------------------------------------------------
  !> Factory function for sll_c_distribution_params, parameters read form input file
  subroutine  sll_s_initial_distribution_new( distribution, dims, file_id, params )
    character(len=*), intent( in    ) :: distribution !< descriptor of the test case
    sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
    sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
    class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.

    sll_int32 :: descriptor
    
    select case( distribution )
    case( "sumcos_onegaussian" )
       allocate( sll_t_params_cos_gaussian :: params )
       descriptor = sll_p_sumcos_onegaussian
    case( "cossum_onegaussian" )
       allocate( sll_t_params_cos_gaussian :: params )
       descriptor = sll_p_cossum_onegaussian
    case( "cossum_twogaussian" )
       allocate( sll_t_params_cos_gaussian :: params )
       descriptor = sll_p_cossum_twogaussian
    case( "sumcos_twogaussian" )
       allocate( sll_t_params_cos_gaussian :: params )
       descriptor = sll_p_sumcos_twogaussian
    case default
       SLL_ERROR('Initial distribution not implemented.','sll_s_initial_distribution_new')
    end select

    select type( params )
    type is( sll_t_params_cos_gaussian )
       call params%init( descriptor, dims, file_id )
    end select
    
    
  end subroutine sll_s_initial_distribution_new

       
  !> Factory function for sll_c_distribution_params, parameters read form input file. Version build upon descriptors
  subroutine  sll_s_initial_distribution_new_descriptor( distribution, dims, file_id, params )
    sll_int32, intent( in    ) :: distribution !< descriptor of the test case
    sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
    sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
    class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.
    
    
    select case( distribution )
    case( sll_p_sumcos_onegaussian )
       allocate( sll_t_params_cos_gaussian :: params )
    case( sll_p_cossum_onegaussian )
       allocate( sll_t_params_cos_gaussian :: params )
    case( sll_p_cossum_twogaussian )
       allocate( sll_t_params_cos_gaussian :: params )
    case( sll_p_sumcos_twogaussian )
       allocate( sll_t_params_cos_gaussian :: params )
    case default
       SLL_ERROR('Initial distribution not implemented.','sll_s_initial_distribution_new')
    end select
    
    select type( params )
    type is( sll_t_params_cos_gaussian )
       call params%init( distribution, dims, file_id )
    end select
    
  end subroutine sll_s_initial_distribution_new_descriptor


! Since assumed shape arrays are not allowed in namelists, we need to define a separate function for each combination of dimensions in x and v. We use a macro to avoid code dublication.
  
#define MAKE_COS_ONEGAUSSIAN_INIT( fname, dimx, dimv, dimalpha )\
  subroutine fname( file_id, params );\
    sll_int32, intent( in ) :: file_id;\
    type( sll_t_params_cos_gaussian ), intent( inout ) :: params; \
    sll_real64 :: kx(dimx); \
    sll_real64 :: alpha(dimalpha); \
    sll_real64 :: v_thermal(dimv); \
    sll_real64 :: v_mean(dimv); \
    sll_int32  :: j; \
    namelist /cos_onegaussian/ kx, alpha, v_thermal, v_mean; \
    params%n_cos = dimalpha;
!-----------------------------------------------------------------

MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"   
    
MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"

MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_2d2v, 2, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90" 

MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_3d3v, 3, 3, 1 )
#include "sll_k_make_cos_onegaussian_init.F90" 

MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"   
    
MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"

MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_2d2v, 2, 2, 2 )
#include "sll_k_make_cos_onegaussian_init.F90" 

MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_3d3v, 3, 3, 3 )
#include "sll_k_make_cos_onegaussian_init.F90"


#define MAKE_COS_TWOGAUSSIAN_INIT( fname, dimx, dimv, dimalpha ) \
  subroutine fname( file_id, params ); \
    sll_int32, intent( in ) :: file_id; \
    type( sll_t_params_cos_gaussian ), intent( inout ) :: params; \
    sll_real64 :: kx(dimx); \
    sll_real64 :: alpha(dimalpha); \
    sll_real64 :: v_thermal_1(dimv); \
    sll_real64 :: v_mean_1(dimv); \
    sll_real64 :: v_thermal_2(dimv); \
    sll_real64 :: v_mean_2(dimv); \
    sll_real64 :: delta; \
    sll_int32  :: j; \
    namelist /cos_twogaussian/ kx, alpha, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta; \
    params%n_cos = dimalpha;
    !-------------------------------
  
MAKE_COS_TWOGAUSSIAN_INIT( cossum_twogaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_twogaussian_init.F90"   
    
MAKE_COS_TWOGAUSSIAN_INIT( cossum_twogaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_twogaussian_init.F90"

MAKE_COS_TWOGAUSSIAN_INIT( cossum_twogaussian_init_2d2v, 2, 2, 1 )
#include "sll_k_make_cos_twogaussian_init.F90" 

MAKE_COS_TWOGAUSSIAN_INIT( cossum_twogaussian_init_3d3v, 3, 3, 1 )
#include "sll_k_make_cos_twogaussian_init.F90" 

MAKE_COS_TWOGAUSSIAN_INIT( sumcos_twogaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_twogaussian_init.F90"   
    
MAKE_COS_TWOGAUSSIAN_INIT( sumcos_twogaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_twogaussian_init.F90"

MAKE_COS_TWOGAUSSIAN_INIT( sumcos_twogaussian_init_2d2v, 2, 2, 2 )
#include "sll_k_make_cos_twogaussian_init.F90" 

MAKE_COS_TWOGAUSSIAN_INIT( sumcos_twogaussian_init_3d3v, 3, 3, 3 )
#include "sll_k_make_cos_twogaussian_init.F90"


end module sll_m_initial_distribution
