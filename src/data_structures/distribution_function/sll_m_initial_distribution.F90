!> @ingroup distribution_function
!> @author Katharina Kormann
!> @brief Parameters to define common initial distributions
!> @details ...
module sll_m_initial_distribution
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_constants, only : &
       sll_p_pi, sll_p_twopi
  
  use sll_m_prob, only: &
       sll_s_normal_cdf_inv

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  implicit none

  public :: sll_p_sumcos_onegaussian, &
       sll_p_cossum_onegaussian, &
       sll_p_sumcos_twogaussian, &
       sll_p_cossum_twogaussian, &
       sll_p_cossum_multigaussian1, &
       sll_c_distribution_params, &
       sll_t_params_cos_gaussian, &
       sll_s_initial_distribution_new, &
       sll_t_params_cos_gaussian_screwpinch, &
       sll_t_params_noise_gaussian, &
       sll_s_initial_distribution_file_new, &
       sll_s_initial_distribution_new_descriptor

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Descriptors for various distributions
  sll_int32, parameter :: sll_p_sumcos_onegaussian = 0 !< Descriptor for (1+\sum cos( kx * x_i))*exp(-0.5(v-v_mean)**2/v_thermal**2)
  sll_int32, parameter :: sll_p_cossum_onegaussian = 1 !< Descriptor for (1+cos( \sum kx_i * x_i))*exp(-0.5(v-v_mean)**2/v_thermal**2)
  sll_int32, parameter :: sll_p_sumcos_twogaussian = 2 !< as sll_p_sumcos_onegaussian but with sum of two Gaussians
  sll_int32, parameter :: sll_p_cossum_twogaussian = 3 !< as sll_p_sumcos_onegaussian but with sum of two Gaussians
  sll_int32, parameter :: sll_p_cossum_multigaussian1 = 4
  sll_int32, parameter :: sll_p_cossum_multigaussian2 = 5
  sll_int32, parameter :: sll_p_cossum_multigaussian11 = 14
  sll_int32, parameter :: sll_p_noise_multigaussian1 = 15
  sll_int32, parameter :: sll_p_noise_multigaussian11 = 16

  !> Abstract data type for parameters of initial distribution
  type, abstract :: sll_c_distribution_params
     sll_int32 :: dims(2) !< Number of spatial and velocity dimensions
     sll_int32               :: n_gaussians !< Number of Gaussians
     sll_real64, allocatable :: v_thermal(:,:) !< variance of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
     sll_real64, allocatable :: v_mean(:,:)    !< mean value of the Gaussian ( first index velocity dimension, second index multiple Gaussians)
     sll_real64, allocatable :: delta(:) !< Portion of each Gaussian
     sll_real64, allocatable :: normal(:) !< Normalization constant of each Gaussian

   contains
     procedure( signature_eval ), deferred :: eval_xv_density !< Evaluate the distribution function
     procedure( signature_evalx), deferred :: eval_x_density  !< Evaluate the charge density
     procedure( signature_evalv), deferred :: eval_v_density  !< Evaluate the v-dependence (integrated over x)
     procedure( signature_empty), deferred :: free            !< Destructor

  end type sll_c_distribution_params

  !> Data type for distribution function with (multiple) Gaussians in v and one plus cosine perturbations in x.
  type, extends(sll_c_distribution_params) :: sll_t_params_cos_gaussian
     sll_real64, allocatable :: kx(:,:)  !< value of the modenumber of the wave which is exited (first index dimension, second index for multiple cosines)
  sll_real64, allocatable :: modnum(:,:)  !< value of the modenumber of the wave which is exited (first index dimension, second index for multiple cosines)
  sll_real64, allocatable :: alpha(:) !< strength of perturbations
  sll_real64, allocatable :: phase_shift(:) !< phase shift in the cosine
  sll_int32               :: n_cos !< Number of cosines

contains
  procedure :: eval_xv_density => sll_f_cos_gaussian !< Evaluate the distribution function
  procedure :: free => free_cos_gaussian  !< Descructor
  procedure :: eval_x_density => sll_f_cos         !< Evaluate the charge density
  procedure :: eval_v_density => sll_f_gaussian    !< Evaluate the v-dependence (integrated over x)
  procedure :: init => cos_gaussian_init  !< Initialization

end type sll_t_params_cos_gaussian

type, extends(sll_c_distribution_params) :: sll_t_params_noise_gaussian
   ! For the Gaussians in velocity

   ! For the white noise in space
   sll_real64 :: alpha !< strength of the noise
   sll_int32 :: n_boxes(3) !< number of boxes for randomization
   sll_real64 :: rdx(3) !< reciprocal of delta_x
   sll_real64, allocatable :: noise_vector(:,:,:)
   type(sll_t_profile_functions) :: profile
   
 contains
   procedure :: eval_xv_density => sll_f_noise_gaussian !< Evaluate the distribution function
   procedure :: free => free_noise_gaussian  !< Descructor
   procedure :: eval_x_density => sll_f_noise         !< Evaluate the charge density
   procedure :: eval_v_density => sll_f_gaussian_pnoise    !< Evaluate the v-dependence (integrated over x)
   procedure :: init => noise_gaussian_init  !< Initialization

end type sll_t_params_noise_gaussian

type, extends(sll_t_params_cos_gaussian) :: sll_t_params_cos_gaussian_screwpinch
   type(sll_t_profile_functions) :: profile
 contains
   procedure :: eval_xv_density => sll_f_cos_gaussian_screwpinch !< Evaluate the distribution function
   procedure :: eval_x_density => sll_f_cos_screwpinch         !< Evaluate the charge density
   procedure :: eval_v_density => sll_f_gaussian_screwpinch    !< Evaluate the v-dependence (integrated over x)
   !procedure :: free => free_cos_gaussian_screwpinch  !< Descructor
   !procedure :: init => cos_gaussian_screwpinch_init  !< Initialization

end type sll_t_params_cos_gaussian_screwpinch

abstract interface
  subroutine signature_empty( self )
    import sll_c_distribution_params
    class( sll_c_distribution_params ), intent(inout) :: self

  end subroutine signature_empty
end interface

abstract interface
  function signature_eval( self, x, v, m ) result( fval )
    use sll_m_working_precision
    import sll_c_distribution_params
    class( sll_c_distribution_params ) :: self
    sll_real64 :: x(:)
    sll_real64 :: v(:)
    sll_real64, optional :: m
    sll_real64 :: fval

  end function signature_eval
end interface

abstract interface
  function signature_evalx( self, x, v ) result( fval )
    use sll_m_working_precision
    import sll_c_distribution_params
    class( sll_c_distribution_params ) :: self
    sll_real64 :: x(:)
    sll_real64, optional :: v(:)
    sll_real64 :: fval

  end function signature_evalx
end interface

abstract interface
  function signature_evalv( self, v, x, m ) result( fval )
    use sll_m_working_precision
    import sll_c_distribution_params
    class( sll_c_distribution_params ) :: self
    sll_real64 :: v(:)
    sll_real64, optional :: x(:)
     sll_real64, optional :: m
    sll_real64 :: fval

  end function signature_evalv
end interface
contains

  !------------------------------------------------------------------------------- 
  ! Define the procedures of the type sll_t_params_cos_gaussian
  !-------------------------------------------------------------------------------

function sll_f_cos_gaussian( self, x, v, m ) result( fval )
 class( sll_t_params_cos_gaussian ) :: self
 sll_real64 :: x(:)
 sll_real64 :: v(:)
 sll_real64, optional :: m
 sll_real64 :: fval
 !local variables
 sll_real64 :: fexp

 fval = sll_f_cos( self, x )
 fexp = sll_f_gaussian( self, v, x, m )

 fval = fval*fexp

end function sll_f_cos_gaussian


function sll_f_cos( self, x, v ) result( fval )
 class( sll_t_params_cos_gaussian ) :: self
 sll_real64 :: x(:)
 sll_real64, optional :: v(:)
 sll_real64 :: fval
 !local variables
 sll_int32  :: j

 fval = 1.0_f64
!!$ do j = 1, self%n_cos
!!$    fval = fval + self%alpha(j) * (cos( self%kx(1,j) * (x(1) + self%modnum(1,j)*v(2))+  self%kx(2,j)*(x(2) - self%modnum(2,j)*v(1)) - self%phase_shift(j)*sll_p_pi ) - &
!!$         0.5_f64*cos( sum(self%kx(:,j) * x) )*exp(-0.5_f64*sum(self%kx(1:2,j)**2) )*sum(self%modnum(:,j)*0.5_f64) )
!!$ end do
!!$
!!$
 do j = 1, self%n_cos
    fval = fval + self%alpha(j) * cos( sum(self%kx(:,j) * x) - self%phase_shift(j)*sll_p_pi )
 end do

!!$ do j=1,self%n_cos
!!$    fval = fval + self%alpha(j) * product( cos( (self%kx(:,j) * x(:)) - self%phase_shift(j)*sll_p_pi ) )
!!$ end do

end function sll_f_cos


function sll_f_gaussian( self, v, x, m ) result( fval )
 class( sll_t_params_cos_gaussian ) :: self
 sll_real64 :: v(:)
 sll_real64, optional :: x(:)
 sll_real64, optional :: m
 sll_real64 :: fval
 !local variables
 sll_int32  :: j

 fval = 0.0_f64
 do j = 1, self%n_gaussians
    fval = fval + self%normal(j)*self%delta(j)* &
         exp( -0.5_f64 * sum( ((v-self%v_mean(:,j))/self%v_thermal(:,j))**2 ) )
 end do

end function sll_f_gaussian

function sll_f_gaussian_pnoise( self, v, x, m ) result( fval )
 class( sll_t_params_noise_gaussian ) :: self
 sll_real64 :: v(:)
 sll_real64, optional :: x(:)
 sll_real64, optional :: m
 sll_real64 :: fval
 !local variables
 sll_int32  :: j

 fval = 0.0_f64
 do j=1,self%n_gaussians
    fval = fval + self%normal(j)*self%delta(j)* &
         exp( -0.5_f64 * sum( ((v-self%v_mean(:,j))/self%v_thermal(:,j))**2 ) )
 end do

end function sll_f_gaussian_pnoise


function sll_f_noise( self, x, v ) result(fval)
  class( sll_t_params_noise_gaussian ) :: self
  sll_real64 :: x(:)
  sll_real64, optional :: v(:)
  sll_real64 :: fval
  !local variables
  sll_int32 :: box(3)
  sll_real64 :: xbox(3)

  xbox = x* self%rdx
  box = floor(xbox)+1

  if(self%dims(1) == 1 ) then
     xbox = xbox - real(box-1,f64)
     fval = (1.0_f64-xbox(1)) * self%noise_vector(box(1)+1,1,1) + xbox(1)* self%noise_vector(box(1),1,1)
  else if (self%dims(1) == 3 ) then
     fval = self%noise_vector(box(1), box(2), box(3) )
  else
     print*, 'wrong dimension in eval_x_noise'
  end if
    
end function sll_f_noise

function sll_f_noise_gaussian( self, x, v, m ) result( fval )
 class( sll_t_params_noise_gaussian ) :: self
 sll_real64 :: x(:)
 sll_real64 :: v(:)
 sll_real64, optional :: m
 sll_real64 :: fval
 
 fval = self%eval_x_density( x ) * self%eval_v_density( v )
 
end function sll_f_noise_gaussian

function sll_f_cos_gaussian_screwpinch( self, x, v, m ) result( fval )
 class( sll_t_params_cos_gaussian_screwpinch ) :: self
 sll_real64 :: x(:)
 sll_real64 :: v(:)
 sll_real64, optional :: m
 sll_real64 :: fval
 !local variables
 sll_real64 :: fexp

 fval = sll_f_cos_screwpinch( self, x, v )
 fexp = sll_f_gaussian_screwpinch( self, v, x, m )

 fval = fval*fexp!/self%profile%rho_0(x(1))

end function sll_f_cos_gaussian_screwpinch


function sll_f_cos_screwpinch( self, x, v ) result( fval )
 class( sll_t_params_cos_gaussian_screwpinch ) :: self
 sll_real64 :: x(:)
 sll_real64, optional :: v(:)
 sll_real64 :: fval
 !local variables
 sll_int32  :: j
 
 fval = 1.0_f64
 do j=1,self%n_cos
    fval = fval + self%alpha(j) * ( cos( sll_p_twopi*sum(self%modnum(:,j) * x) - self%kx(2,j)*v(1) - self%phase_shift(j)*sll_p_pi  ) - 0.5_f64*cos( sll_p_twopi*sum(self%modnum(:,j) * x) - self%phase_shift(j)*sll_p_pi )*exp(-0.5_f64*self%profile%T_i(x(1) + self%kx(1,j)*v(2))*self%kx(2,j)**2 ) ) * self%profile%radial_distrib(x(1))
 end do

 fval = fval! * self%profile%rho_0(x(1))

end function sll_f_cos_screwpinch

function sll_f_gaussian_screwpinch( self, v, x, m ) result( fexp )
 class( sll_t_params_cos_gaussian_screwpinch ) :: self
 sll_real64 :: v(:)
 sll_real64, optional :: x(:)
 sll_real64, optional :: m
 sll_real64 :: fexp
 !local variables
 sll_int32  :: j
 
 fexp = 0.0_f64
 do j=1,self%n_gaussians
    fexp = fexp + self%profile%rho_0(x(1))/sqrt(sll_p_twopi*self%profile%T_i(x(1))/m)**3 * &
         exp( -0.5_f64 * sum( (v-self%v_mean(:,j))**2/(self%profile%T_i(x(1))/m) ) )
 end do

end function sll_f_gaussian_screwpinch


subroutine free_cos_gaussian( self )
 class( sll_t_params_cos_gaussian ), intent( inout ) :: self

 if (allocated(self%kx)) deallocate(self%kx)
 if (allocated(self%modnum)) deallocate(self%modnum)
 if (allocated(self%alpha)) deallocate(self%alpha)
 if (allocated(self%phase_shift)) deallocate(self%phase_shift)   
 if (allocated(self%v_thermal)) deallocate(self%v_thermal)
 if (allocated(self%v_mean)) deallocate(self%v_mean)
 if (allocated(self%normal)) deallocate(self%normal)
 if (allocated(self%delta)) deallocate(self%delta)

end subroutine free_cos_gaussian

subroutine free_noise_gaussian( self )
 class( sll_t_params_noise_gaussian ), intent( inout ) :: self

 if (allocated(self%v_thermal)) deallocate(self%v_thermal)
 if (allocated(self%v_mean)) deallocate(self%v_mean)
 if (allocated(self%normal)) deallocate(self%normal)
 if (allocated(self%delta)) deallocate(self%delta)
 if (allocated(self%noise_vector)) deallocate(self%noise_vector)
 
end subroutine free_noise_gaussian

subroutine noise_gaussian_init( self, n_gaussians, dims, file_id, profile )
  class( sll_t_params_noise_gaussian ), intent( out ) :: self
  sll_int32, intent( in    ) :: n_gaussians !< descriptor of the test case
  sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
  sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
  type(sll_t_profile_functions), optional :: profile

  if( present(profile) ) then
     self%profile = profile
  end if
  self%dims = dims
  self%n_gaussians = n_gaussians
  select case( self%dims(1) )
  case (1)
     call noise_gaussian_init_1d2v( self, file_id )
  case(3)
     call noise_gaussian_init_3d3v( self, file_id )
  end select
end subroutine noise_gaussian_init

subroutine noise_gaussian_init_1d2v( self, file_id )
  class( sll_t_params_noise_gaussian ), intent( inout ) :: self
  sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
  sll_real64 :: alpha
  sll_real64 :: v_thermal_1(self%n_gaussians)
  sll_real64 :: v_mean_1(self%n_gaussians)
  sll_real64 :: v_thermal_2(self%n_gaussians)
  sll_real64 :: v_mean_2(self%n_gaussians)
  sll_real64 :: delta(self%n_gaussians)
  sll_int32 :: n_boxes
  sll_real64 :: domain   !< domain length (to set delta x)
  sll_int32 :: i, j
  sll_int32 :: rnd_seed_size
  sll_int32, allocatable :: rnd_seed(:)
  sll_real64 :: rnd
  
  namelist /noise_multigaussian/ alpha, n_boxes, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta, domain

  

  read(file_id, noise_multigaussian)

  allocate( self%noise_vector(1:n_boxes+1,1,1) )
  allocate( self%v_thermal(1:self%dims(2), 1:self%n_gaussians) )
  allocate( self%v_mean(1:self%dims(2), 1:self%n_gaussians) )
  allocate( self%normal(1:self%n_gaussians) )
  allocate( self%delta(1:self%n_gaussians) )

  
 
  self%alpha = alpha
  self%v_thermal(1,:) = v_thermal_1
  self%v_mean(1,:) = v_mean_1
  self%v_thermal(2,:) = v_thermal_2
  self%v_mean(2,:) = v_mean_2
  self%delta = delta
  self%n_boxes = n_boxes
  
  do j=1,self%n_gaussians
     self%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(self%dims(2),f64))*&
          product(self%v_thermal(:,j)))   
  end do

  ! Set random seed: Same on all processes
  call random_seed(size=rnd_seed_size)
  allocate( rnd_seed(rnd_seed_size) )
  do i=1, rnd_seed_size
     rnd_seed(i) = 15*i
  end do
  call random_seed(put=rnd_seed)
  do i=1,n_boxes
     call random_number( rnd )
     call sll_s_normal_cdf_inv( rnd, 1.0_f64, alpha, self%noise_vector(i,1,1) )
  end do

  ! Make sure the function integrates up to one
  rnd = sum(self%noise_vector(1:n_boxes,1,1))/real(n_boxes, f64)
  rnd = rnd - 1.0_f64
  self%noise_vector(1:n_boxes,1,1) = self%noise_vector(1:n_boxes,1,1) - rnd
  
  ! Periodic boundary conditions
  self%noise_vector(n_boxes+1,1,1) = self%noise_vector(1,1,1)

  self%rdx = real( self%n_boxes, f64 ) / domain

end subroutine noise_gaussian_init_1d2v

subroutine noise_gaussian_init_3d3v( self, file_id )
  class( sll_t_params_noise_gaussian ), intent( inout ) :: self
  sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
  sll_real64 :: alpha
  sll_real64 :: v_thermal(3,self%n_gaussians)
  sll_real64 :: v_mean(3,self%n_gaussians)
  sll_real64 :: delta(self%n_gaussians)
  sll_int32 :: n_boxes(3)
  sll_real64 :: domain(3)   !< domain length (to set delta x)
  sll_int32 :: i, j, k
  sll_int32 :: rnd_seed_size
  sll_int32, allocatable :: rnd_seed(:)
  sll_real64 :: rnd
  
  namelist /noise_multigaussian/ alpha, n_boxes, v_thermal, v_mean, delta, domain

  read(file_id, noise_multigaussian)

  allocate( self%noise_vector(1:n_boxes(1)+1, 1:n_boxes(2)+1, 1:n_boxes(3)+1) )
  allocate( self%v_thermal(1:self%dims(2), 1:self%n_gaussians) )
  allocate( self%v_mean(1:self%dims(2), 1:self%n_gaussians) )
  allocate( self%normal(1:self%n_gaussians) )
  allocate( self%delta(1:self%n_gaussians) )

  
 
  self%alpha = alpha
  self%v_thermal = v_thermal
  self%v_mean = v_mean
  self%delta = delta
  self%n_boxes = n_boxes
  
  do j=1,self%n_gaussians
     self%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(self%dims(2),f64))*&
          product(self%v_thermal(:,j)))   
  end do

  ! Set random seed: Same on all processes
  call random_seed(size=rnd_seed_size)
  allocate( rnd_seed(rnd_seed_size) )
  do i=1, rnd_seed_size
     rnd_seed(i) = 15*i
  end do
  call random_seed(put=rnd_seed)
  do i=1,n_boxes(1)
     do j = 1, n_boxes(2)
        do k = 1, n_boxes(3)
           call random_number( rnd )
           call sll_s_normal_cdf_inv( rnd, 1.0_f64, alpha, self%noise_vector(i,j,k) )
        end do
     end do
  end do

  ! Make sure the function integrates up to one
  do j = 1, n_boxes(2)
     do k = 1, n_boxes(3)
        rnd = sum(self%noise_vector(1:n_boxes(1),j,k))/real(n_boxes(1), f64)
        rnd = rnd - 1.0_f64
        self%noise_vector(1:n_boxes(1),j,k) = self%noise_vector(1:n_boxes(1),j,k) - rnd
     end do
  end do
  ! Periodic boundary conditions
  !self%noise_vector(n_boxes+1) = self%noise_vector(1)

  self%rdx = real( self%n_boxes, f64 ) / domain

end subroutine noise_gaussian_init_3d3v


subroutine cos_gaussian_init( self, descriptor, dims, file_id, profile )
 class( sll_t_params_cos_gaussian ), intent( out ) :: self
 sll_int32, intent( in    ) :: descriptor !< descriptor of the test case
 sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
 sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
 type(sll_t_profile_functions), optional :: profile

 select type( self)
 type is(sll_t_params_cos_gaussian_screwpinch)
    if( present(profile) ) then
       self%profile = profile
    else
       print*, 'Error: no profile given'
       STOP
    end if
 end select
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
       select case ( self%dims(2) )
       case ( 2)
          call sumcos_onegaussian_init_2d2v( file_id, self )
       case(3)
          call sumcos_onegaussian_init_2d3v( file_id, self )
       end select
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
       select case ( self%dims(2) )
       case ( 2)
          call cossum_onegaussian_init_2d2v( file_id, self )
       case (3)
          call cossum_onegaussian_init_2d3v( file_id, self )
       end select
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
 case ( sll_p_cossum_multigaussian1 )
    select case (self%dims(1) )
    case(1)
       select case( self%dims(2) )
       case(2)
          call cossum_multigaussian_init_1d2v( file_id, self, 1 )
       end select
    end select
 case ( sll_p_cossum_multigaussian2 )
    select case (self%dims(1) )
    case(1)
       select case( self%dims(2) )
       case(2)
          call cossum_multigaussian_init_1d2v( file_id, self, 2 )
       end select
    end select
 case ( sll_p_cossum_multigaussian11 )
    select case (self%dims(1) )
    case(1)
       select case( self%dims(2) )
       case(2)
          call cossum_multigaussian_init_1d2v( file_id, self, 11 )
       end select
    end select
 case default
    SLL_ERROR('Initial distribution not implemented.','cos_gaussian_init')
 end select


end subroutine cos_gaussian_init


!------------------------------------------------------------------------------- 
! Factory function with specific functions called depending on chosen distribution type
!-------------------------------------------------------------------------------
!> Factory function for sll_c_distribution_params, parameters read form input file
subroutine  sll_s_initial_distribution_new( distribution, dims, file_id, params, profile )
 character(len=*), intent( in    ) :: distribution !< descriptor of the test case
 sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
 sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
 class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.
 type(sll_t_profile_functions), optional :: profile

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
 case( "cossum_multigaussian1" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian1
 case( "cossum_multigaussian2" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian2
 case( "cossum_multigaussian11" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian11
 case ( "noise_multigaussian1" )
    allocate( sll_t_params_noise_gaussian :: params )
    descriptor = sll_p_noise_multigaussian1
 case ( "noise_multigaussian11" )
    allocate( sll_t_params_noise_gaussian :: params )
    descriptor = sll_p_noise_multigaussian11
 case( "cossum_onegaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_cossum_onegaussian
 case( "sumcos_onegaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_sumcos_onegaussian
 case( "cossum_twogaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_cossum_twogaussian
 case( "sumcos_twogaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_sumcos_twogaussian
 case default
    SLL_ERROR('Initial distribution not implemented.','sll_s_initial_distribution_new')
 end select

 select type( params )
 type is( sll_t_params_cos_gaussian )
    call params%init( descriptor, dims, file_id )
 type is( sll_t_params_noise_gaussian )
    select case ( descriptor )
    case (sll_p_noise_multigaussian1 )
       call params%init( 1, dims, file_id, profile )
    case (sll_p_noise_multigaussian11 )
       call params%init( 11, dims, file_id, profile )
    case default
       !call params%init( descriptor, dims, file_id, profile )
    end select
 type is( sll_t_params_cos_gaussian_screwpinch )
    call params%init( descriptor, dims, file_id, profile )
 end select


end subroutine sll_s_initial_distribution_new


!------------------------------------------------------------------------------- 
! Factory function with specific functions called depending on chosen distribution type
!-------------------------------------------------------------------------------
!> Factory function for sll_c_distribution_params, parameters read form input file
subroutine  sll_s_initial_distribution_file_new( dims, nml_file, params, profile )
 sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
 character(*), intent( in    ) :: nml_file    !< nml-file with parameters in unified format
 class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.
 type(sll_t_profile_functions), optional :: profile

 character(len=256) :: type !< descriptor of the test case
 sll_int32 :: descriptor
 sll_int32 :: file_id, ierr

 namelist /initial_distribution/ type

 ! Read what type
 open(newunit=file_id, file=trim(nml_file), iostat=ierr)
 if ( ierr /= 0 ) then       
    SLL_ERROR('NML file could not be opened.','sll_s_initial_distribution_file_new')
 end if
 read( file_id, initial_distribution )

 select case( type )
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
 case( "cossum_multigaussian1" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian1
 case( "cossum_multigaussian2" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian2
 case( "cossum_multigaussian11" )
    allocate( sll_t_params_cos_gaussian :: params )
    descriptor = sll_p_cossum_multigaussian11
 case ( "noise_multigaussian1" )
    allocate( sll_t_params_noise_gaussian :: params )
    descriptor = sll_p_noise_multigaussian1
 case ( "noise_multigaussian11" )
    allocate( sll_t_params_noise_gaussian :: params )
    descriptor = sll_p_noise_multigaussian11
 case( "cossum_onegaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_cossum_onegaussian
 case( "sumcos_onegaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_sumcos_onegaussian
 case( "cossum_twogaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_cossum_twogaussian
 case( "sumcos_twogaussian_screwpinch" )
    allocate( sll_t_params_cos_gaussian_screwpinch :: params )
    descriptor = sll_p_sumcos_twogaussian
 case default
    SLL_ERROR('Initial distribution not implemented.','sll_s_initial_distribution_file_new')
 end select
 
 select type( params )
 type is( sll_t_params_cos_gaussian )
    call params%init( descriptor, dims, file_id )
 type is( sll_t_params_noise_gaussian )
    select case ( descriptor )
    case (sll_p_noise_multigaussian1 )
       call params%init( 1, dims, file_id, profile )
    case (sll_p_noise_multigaussian11 )
       call params%init( 11, dims, file_id, profile )
    case default
       !call params%init( descriptor, dims, file_id, profile )
    end select
 type is( sll_t_params_cos_gaussian_screwpinch )
    call params%init( descriptor, dims, file_id, profile )
 end select

 close(file_id)

end subroutine sll_s_initial_distribution_file_new


!> Factory function for sll_c_distribution_params, parameters read form input file. Version build upon descriptors
subroutine  sll_s_initial_distribution_new_descriptor( distribution, dims, file_id, params, profile )
 sll_int32, intent( in    ) :: distribution !< descriptor of the test case
 sll_int32, intent( in    ) :: dims(2) !< number of spatial and velocity dimensions
 sll_int32, intent( in    ) :: file_id    !< nml-file with parameters in unified format
 class(sll_c_distribution_params), allocatable, intent(   out ) ::  params    !< real array specifying the parameters for the given test case in the predefined order.
 type(sll_t_profile_functions), optional :: profile


 select case( distribution )
 case( sll_p_sumcos_onegaussian )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_cossum_onegaussian )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_cossum_twogaussian )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_sumcos_twogaussian )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_cossum_multigaussian1 )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_cossum_multigaussian2 )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_cossum_multigaussian11 )
    allocate( sll_t_params_cos_gaussian :: params )
 case( sll_p_noise_multigaussian1 )
    allocate( sll_t_params_noise_gaussian :: params )
 case( sll_p_noise_multigaussian11 )
    allocate( sll_t_params_noise_gaussian :: params )
 case default
    SLL_ERROR('Initial distribution not implemented.','sll_s_initial_distribution_new')
 end select

 select type( params )
 type is( sll_t_params_cos_gaussian )
    call params%init( distribution, dims, file_id )
 type is( sll_t_params_noise_gaussian )
    select case ( distribution)
    case (sll_p_noise_multigaussian1 )
       call params%init( 1, dims, file_id, profile )
    case (sll_p_noise_multigaussian11 )
       call params%init( 11, dims, file_id, profile )
    end select
 type is( sll_t_params_cos_gaussian_screwpinch )
    call params%init( distribution, dims, file_id )
 end select

end subroutine sll_s_initial_distribution_new_descriptor


! Since assumed shape arrays are not allowed in namelists, we need to define a separate function for each combination of dimensions in x and v. We use a macro to avoid code dublication.

#define MAKE_COS_ONEGAUSSIAN_INIT( fname, dimx, dimv, dimalpha )\
subroutine fname( file_id, params );\
 sll_int32, intent( in ) :: file_id;\
 type( sll_t_params_cos_gaussian ), intent( inout ) :: params; \
 sll_real64 :: kx(dimx)= 0.0_f64; \
 sll_real64 :: modnum(dimx)= 0.0_f64; \
 sll_real64 :: alpha(dimalpha); \
 sll_real64 :: phase_shift(dimalpha) = 0.0_f64; \
 sll_real64 :: v_thermal(dimv); \
 sll_real64 :: v_mean(dimv)= 0.0_f64; \
 sll_int32  :: j; \
 namelist /cos_onegaussian/ kx, modnum, alpha, phase_shift, v_thermal, v_mean; \
 params%n_cos = dimalpha;
 !-----------------------------------------------------------------

 MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"   

 MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"

 MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_2d2v, 2, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90" 

 MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_2d3v, 2, 3, 1 )
#include "sll_k_make_cos_onegaussian_init.F90" 

 MAKE_COS_ONEGAUSSIAN_INIT( cossum_onegaussian_init_3d3v, 3, 3, 1 )
#include "sll_k_make_cos_onegaussian_init.F90" 

 MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_1d1v, 1, 1, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"   

 MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_1d2v, 1, 2, 1 )
#include "sll_k_make_cos_onegaussian_init.F90"

 MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_2d2v, 2, 2, 2 )
#include "sll_k_make_cos_onegaussian_init.F90" 

 MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_2d3v, 2, 3, 2 )
#include "sll_k_make_cos_onegaussian_init.F90" 

 MAKE_COS_ONEGAUSSIAN_INIT( sumcos_onegaussian_init_3d3v, 3, 3, 3 )
#include "sll_k_make_cos_onegaussian_init.F90"


#define MAKE_COS_TWOGAUSSIAN_INIT( fname, dimx, dimv, dimalpha ) \
 subroutine fname( file_id, params ); \
   sll_int32, intent( in ) :: file_id; \
   type( sll_t_params_cos_gaussian ), intent( inout ) :: params; \
   sll_real64 :: kx(dimx)= 0.0_f64; \
   sll_real64 :: modnum(dimx)= 0.0_f64; \
   sll_real64 :: alpha(dimalpha); \
   sll_real64 :: phase_shift(dimalpha) = 0.0_f64; \
   sll_real64 :: v_thermal_1(dimv); \
   sll_real64 :: v_mean_1(dimv); \
   sll_real64 :: v_thermal_2(dimv); \
   sll_real64 :: v_mean_2(dimv); \
   sll_real64 :: delta; \
   sll_int32  :: j; \
   namelist /cos_twogaussian/ kx, modnum, alpha, phase_shift, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta; \
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


   !> 1d2v subroutine for initialization of sum of arbitrary number of Gaussians.
   !> Note that v_thermal_1/2 refers her to the velocity dimension instead of the number of the Gaussian in the sum as for the twogaussian case.
   subroutine cossum_multigaussian_init_1d2v( file_id, params, n_gaussians )
     sll_int32, intent( in ) :: file_id
     type( sll_t_params_cos_gaussian ), intent( inout ) :: params
     sll_int32, intent( in ) :: n_gaussians
     sll_real64 :: kx(1)
     sll_real64 :: modnum(1)
     sll_real64 :: alpha(1)
     sll_real64 :: phase_shift(1) = 0.0_f64
     sll_real64 :: v_thermal_1(n_gaussians)
     sll_real64 :: v_mean_1(n_gaussians)
     sll_real64 :: v_thermal_2(n_gaussians)
     sll_real64 :: v_mean_2(n_gaussians)
     sll_real64 :: delta(n_gaussians)
     sll_int32  :: j
     namelist /cos_multigaussian/ kx, modnum, alpha, phase_shift, v_thermal_1, v_mean_1, v_thermal_2, v_mean_2, delta
     params%n_cos = 1
     params%n_gaussians = n_gaussians

     read(file_id, cos_multigaussian)

     allocate( params%alpha(1:params%n_cos) );
     allocate( params%phase_shift(1:params%n_cos) );
     allocate( params%kx(1:params%dims(1), 1:params%n_cos) )
     allocate( params%modnum(1:params%dims(1), 1:params%n_cos) )
     allocate( params%v_thermal(1:params%dims(2), 1:params%n_gaussians) )
     allocate( params%v_mean(1:params%dims(2), 1:params%n_gaussians) )
     allocate( params%normal(1:params%n_gaussians) )
     allocate( params%delta(1:params%n_gaussians) )

     if ( params%n_cos == 1 ) then
        params%kx(:,1) = kx
        params%modnum(:,1) = modnum
     else
        params%modnum = 0.0_f64
        do j=1, params%n_cos
           params%modnum(j,j) = modnum(j)
        end do
     end if
     params%alpha = alpha
     params%phase_shift=phase_shift
     params%v_thermal(1,:) = v_thermal_1
     params%v_mean(1,:) = v_mean_1
     params%v_thermal(2,:) = v_thermal_2
     params%v_mean(2,:) = v_mean_2
     params%delta = delta

     do j=1,params%n_gaussians
        params%normal(j) = 1.0_f64/(sll_p_twopi**(0.5_f64*real(params%dims(2),f64))*&
             product(params%v_thermal(:,j)))

     end do

   end subroutine cossum_multigaussian_init_1d2v


 end module sll_m_initial_distribution
