
module sll_m_distribution_function_initializer_6d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_collective, only : &
       sll_o_collective_allreduce, &
       sll_v_world_collective
  
  use sll_m_constants, only: &
       sll_p_twopi

  use sll_mpi, only : &
       mpi_sum
  
#ifdef _OPENMP
  use omp_lib
#define OMP_COLLAPSE collapse(1)
#define OMP_SCHEDULE schedule(static)
#endif

  implicit none

  public :: sll_s_set_local_grid, &
       sll_s_set_local_grid_en,&
       sll_t_array, &
       sll_c_distribution_params_6d, &
       sll_s_distribution_params_6d_new, &
       sll_s_distribution_initializer_6d, &
       sll_p_landau_prod, &
       sll_p_landau_sum, &
       sll_p_landau_diag, &
       sll_p_landau_sum_df, &
       sll_p_pslab, &
       sll_p_delta, &
       sll_p_pslab2, &
       sll_p_twogaussian_sum,&
       sll_s_compute_velocity_transformation,&
       sll_s_compute_velocity_transformation_en

  private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_landau_prod = 0
  sll_int32, parameter :: sll_p_landau_sum = 1
  sll_int32, parameter :: sll_p_landau_diag = 2
  sll_int32, parameter :: sll_p_pslab = 3
  sll_int32, parameter :: sll_p_twogaussian_prod = 4
  sll_int32, parameter :: sll_p_twogaussian_sum = 5
  sll_int32, parameter :: sll_p_twogaussian_diag = 6
  sll_int32, parameter :: sll_p_pslab2 = 7
  sll_int32, parameter :: sll_p_delta = 9
  sll_int32, parameter :: sll_p_landau_sum_df = 8

  !> Array type
  type :: sll_t_array
     sll_real64, allocatable ::  vals(:) !< Values of the array
  end type sll_t_array

  
  !> Abstract parameter type
  type, abstract :: sll_c_distribution_params_6d
     sll_int32 :: distrib_type

   contains
     procedure( signature_init   ), deferred  :: init
     procedure( signature_eval   ), deferred  :: eval
     procedure( signature_eval_v ), deferred  :: eval_v

  end type sll_c_distribution_params_6d

  !> Type to define parameters of Landau damping
  type, extends(sll_c_distribution_params_6d) :: sll_t_landau_sum_parameters_6d
     sll_real64 :: v_thermal(3)
     sll_real64 :: alpha(3)
     sll_real64 :: kx(3)
     sll_real64 :: phase(3)
     sll_real64 :: v_max(3)
     sll_real64 :: factor

   contains

     procedure :: init => init_landau_sum
     procedure :: eval => eval_landau_sum
     procedure :: eval_v => eval_landau_v_sum
     
  end type sll_t_landau_sum_parameters_6d
  

  !> Type to define parameters of Landau damping
  type, extends(sll_t_landau_sum_parameters_6d) :: sll_t_landau_sum_df_parameters_6d
     
   contains
     procedure :: eval => eval_landau_sum_df
     
  end type sll_t_landau_sum_df_parameters_6d
  
  !> Type to define parameters of Landau damping
  type, extends(sll_c_distribution_params_6d) :: sll_t_landau_prod_parameters_6d
     sll_real64 :: v_thermal(3)
     sll_real64 :: alpha
     sll_real64 :: kx(3)
     sll_real64 :: v_max(3)
     sll_real64 :: factor

   contains

     procedure :: init => init_landau_prod
     procedure :: eval => eval_landau_prod
     procedure :: eval_v => eval_landau_v
     
  end type sll_t_landau_prod_parameters_6d

  !> Type to define parameters of Landau damping
  type, extends(sll_t_landau_prod_parameters_6d) :: sll_t_landau_diag_parameters_6d
     
   contains
     procedure :: eval => eval_landau_diag
     
  end type sll_t_landau_diag_parameters_6d


  !> Type to specify parameter for double Gaussian (includes  bump-on-tail and TSI)
  type, extends(sll_c_distribution_params_6d) :: sll_t_twogaussian_parameters_6d
     sll_real64 :: v_thermal(2,3)
     sll_real64 :: v_mean(2,3)
     sll_real64 :: alpha
     sll_real64 :: kx(3)
     sll_real64 :: v_max(3)
     sll_real64 :: factor
     sll_real64 :: delta(2)

   contains

     procedure :: init => init_twogaussian
     procedure :: eval => eval_twogaussian_sum
     procedure :: eval_v => eval_twogaussian_v
     
  end type sll_t_twogaussian_parameters_6d

  !> Type
  type :: sll_t_itg_parameters_6d
     sll_real64 :: v_thermal(3)
     sll_real64 :: kappa_n0
     sll_real64 :: kappa_Ti
     sll_real64 :: kappa_Te
     sll_real64 :: delta_rn0
     sll_real64 :: delta_rTi
     sll_real64 :: delta_rTe
     sll_real64 :: kx(3)
     sll_real64 :: alpha
     sll_real64 :: B0
     sll_real64 :: C_Te
     sll_real64 :: rp
  end type sll_t_itg_parameters_6d

  type, extends(sll_c_distribution_params_6d) :: sll_t_pslab_parameters_6d
     sll_real64 :: kappa_Ti
     sll_real64 :: kx(3)
     sll_real64 :: alpha
     sll_real64 :: B0
     sll_real64 :: C_Te
     sll_real64 :: factor
     sll_real64 :: v_thermal(3)
     sll_real64 :: v_mean(3)
     
   contains
     procedure :: init => init_pslab
     procedure :: eval => eval_pslab
     procedure :: eval_v => eval_pslab_v
     
  end type sll_t_pslab_parameters_6d
  
  


 type, extends(sll_t_pslab_parameters_6d) :: sll_t_pslab2_parameters_6d
      
  sll_real64 :: randArray(509)= (/-0.2664, 0.014, -0.3229, -0.3171, -0.0291, -0.0062, -0.4745, 0.4168, &
                        0.1802, -0.3243, 0.3694, -0.1246, -0.1856, -0.0366, 0.0695, -0.2344, &
                        -0.4531, 0.0384, -0.0048, -0.4186, -0.4135, 0.0035, 0.473, -0.3343, &
                        0.4867, -0.0914, 0.2798, -0.0493, -0.4807, -0.275, 0.3846, -0.3377, &
                        -0.4438, -0.3822, -0.3808, -0.2062, 0.3435, -0.465, -0.4356, 0.1279, &
                        -0.4913, 0.2428, 0.0293, -0.0828, -0.3102, 0.3547, 0.4727, 0.1805, &
                        -0.4724, -0.1745, -0.3345, -0.2994, 0.0773, 0.0896, 0.2664, 0.299, &
                        -0.1944, -0.2085, 0.3635, 0.1152, -0.4352, -0.1202, -0.1496, 0.4634, &
                        0.4438, -0.3141, 0.3368, -0.3054, 0.2811, -0.2799, -0.453, -0.3705, &
                        -0.3277, 0.3093, 0.293, -0.3352, 0.2645, -0.2042, -0.0412, -0.1305, &
                        0.3631, 0.2099, 0.1429, 0.3308, 0.1216, -0.3194, -0.4364, -0.4331, &
                        0.2451, -0.1217, -0.0064, -0.4491, 0.474, 0.0612, 0.4629, 0.2188, &
                        -0.2955, -0.2415, 0.309, 0.3778, 0.1484, 0.1198, -0.0777, 0.3699, &
                        0.415, -0.2132, 0.01, -0.045, -0.4163, 0.3464, -0.1753, -0.2682, &
                        0.3636, -0.4797, 0.1802, -0.1344, -0.3197, 0.3575, -0.0722, 0.372, &
                        -0.2739, -0.3028, -0.1227, -0.0174, -0.4454, 0.3344, 0.3524, -0.481, &
                        0.0944, -0.325, -0.4205, 0.3596, -0.3042, -0.2806, -0.429, 0.2148, &
                        -0.1049, 0.1622, -0.0293, 0.2618, 0.3408, -0.3008, 0.061, 0.2094, &
                        0.1375, 0.1928, -0.0813, -0.484, 0.2006, 0.0031, -0.1039, -0.4471, &
                        0.4876, -0.007, -0.2788, 0.1259, 0.4924, -0.2759, 0.4467, 0.4321, &
                        -0.0065, 0.3155, -0.0429, 0.1107, 0.2936, -0.239, 0.382, 0.2367, &
                        -0.3248, -0.4305, -0.304, -0.1935, -0.1988, 0.2261, 0.1036, -0.1348, &
                        -0.3115, 0.363, 0.2836, 0.2853, 0.4983, 0.4221, 0.3994, 0.4019, &
                        0.2773, -0.2255, -0.0043, -0.0997, -0.2185, 0.1486, -0.4017, 0.3383, &
                        0.167, -0.3125, -0.1966, 0.2756, -0.0057, 0.3256, 0.1876, -0.3863, &
                        -0.0464, 0.2262, 0.2781, -0.246, -0.1694, -0.0606, -0.4085, 0.4251, &
                        0.2965, 0.3739, -0.4765, 0.2349, -0.2762, -0.3162, 0.1564, -0.058, &
                        0.264, -0.0621, 0.2979, -0.1459, 0.4982, -0.1117, -0.3846, -0.2493, &
                        -0.1303, -0.2645, 0.2048, -0.3789, 0.4384, -0.196, -0.2605, -0.3464, &
                        0.443, 0.3775, -0.093, 0.1543, 0.3666, 0.0055, 0.015, 0.124, &
                        0.3759, -0.0529, -0.4103, -0.1574, -0.0722, -0.116, 0.3529, 0.3186, &
                        -0.4562, -0.3968, -0.4112, 0.0024, 0.0191, 0.2281, -0.3646, -0.3362, &
                        0.4484, -0.1759, -0.0501, 0.4268, 0.1152, -0.4264, -0.1195, 0.2677, &
                        0.0351, -0.1557, -0.0494, 0.3777, -0.2088, -0.4571, -0.4257, 0.004, &
                        0.0283, -0.2383, 0.3535, 0.3623, -0.4942, 0.0078, 0.0773, -0.2971, &
                        -0.2633, 0.081, -0.2694, 0.2239, -0.2029, 0.1425, 0.0181, 0.2826, &
                        -0.3549, -0.3291, -0.3042, 0.0888, -0.161, 0.01, 0.4109, -0.026, &
                        0.4173, 0.1474, -0.377, -0.4248, 0.3494, -0.4673, -0.2663, 0.0858, &
                        0.2815, -0.0831, 0.1734, -0.1101, 0.3472, -0.4237, -0.499, -0.4278, &
                        0.0403, -0.2564, 0.106, -0.1339, 0.2183, -0.4593, 0.3897, 0.1814, &
                        0.3096, 0.3329, -0.2162, -0.1349, -0.3752, -0.1155, -0.4494, -0.0356, &
                        -0.1776, -0.2952, -0.4692, -0.2821, 0.1987, 0.1159, -0.1637, -0.3908, &
                        -0.4081, -0.311, 0.1897, -0.4715, 0.4399, -0.4585, -0.2052, -0.3093, &
                        -0.0484, -0.0195, 0.3651, -0.3401, -0.2006, 0.0011, 0.2373, -0.1851, &
                        -0.2044, 0.0305, -0.3288, -0.3535, -0.051, 0.2206, -0.0382, 0.325, &
                        0.4449, -0.0674, 0.1765, 0.2122, 0.2747, 0.364, 0.3926, 0.2361, &
                        -0.4706, -0.0562, 0.0346, -0.4615, 0.2632, 0.2095, 0.346, -0.0197, &
                        -0.3338, 0.4065, -0.3779, -0.0901, 0.0981, -0.1179, -0.1619, -0.2831, &
                        0.4709, 0.3433, -0.0723, 0.1744, -0.2182, -0.0365, 0.3851, -0.4263, &
                        0.2903, -0.0463, 0.4824, -0.349, -0.2941, 0.1386, 0.228, -0.026, &
                        0.1105, -0.3866, 0.417, -0.0021, 0.4718, -0.4834, -0.0377, 0.2265, &
                        0.4044, -0.4975, 0.0467, -0.4142, -0.1741, 0.2117, 0.4677, -0.1504, &
                        -0.1259, -0.4895, 0.1251, -0.45, 0.4062, -0.2606, 0.0356, 0.2806, &
                        -0.4247, 0.2262, -0.3132, -0.0803, -0.0144, 0.2835, -0.0234, -0.2817, &
                        -0.1099, -0.0675, -0.3592, 0.131, -0.0506, 0.4867, 0.0172, -0.1996, &
                        0.2327, -0.1522, 0.3792, 0.0062, -0.1455, 0.4794, -0.4917, 0.3584, &
                        0.387, -0.4079, 0.2842, -0.4815, 0.1635, 0.0122, -0.4488, -0.3681, &
                        -0.4737, 0.3576, 0.2158, 0.0163, 0.0637, -0.4268, -0.2131, 0.2677, &
                        -0.0808, 0.0681, 0.3314, 0.1626, -0.489, -0.2921, 0.1211, 0.3653, &
                        0.224, 0.427, -0.0387, 0.4869, -0.0176, -0.4717, 0.4876, -0.2194, &
                        -0.4679, -0.1669, 0.1821, 0.0641, -0.1187, 0.4273, -0.1504, 0.2218, &
                        0.1033, 0.2177, -0.4011, 0.0596, 0.0083, 0.4719, -0.3228, -0.1893, &
                        -0.3938, 0.233, -0.3116, 0.1987, -0.0399, -0.3335, 0.4565, 0.0106, &
                        0.1017, 0.1004, 0.3754, 0.4765, -0.1246 /) 

   contains
     procedure :: eval => eval_pslab2
     
 
     
  end type sll_t_pslab2_parameters_6d
 type, extends(sll_t_pslab2_parameters_6d) :: sll_t_delta_parameters_6d

   contains
     procedure :: init => init_delta
     procedure :: eval => eval_delta
     
  end type sll_t_delta_parameters_6d
! 

  
  abstract interface
     subroutine signature_init ( self, file_id )
       use sll_m_working_precision
       import sll_c_distribution_params_6d
       class( sll_c_distribution_params_6d ), intent( out ) :: self
       sll_int32 :: file_id
       
     end subroutine signature_init
  end interface

  
  abstract interface
     function signature_eval ( self, x ) result(res)
       use sll_m_working_precision
       import sll_c_distribution_params_6d
       class( sll_c_distribution_params_6d ), intent( in ) :: self
       sll_real64,                 intent(in)   :: x(6)
       sll_real64 :: res
       
       
     end function signature_eval
  end interface
  
  abstract interface
     function signature_eval_v ( self, x ) result(res)
       use sll_m_working_precision
       import sll_c_distribution_params_6d
       class( sll_c_distribution_params_6d ), intent( in ) :: self
       sll_real64,                 intent(in)   :: x(3)
       sll_real64 :: res
       
       
     end function signature_eval_v
  end interface

contains

  subroutine sll_s_set_local_grid(&
       local_sizes, &
       indices_min, &
       eta_min, &
       delta_eta, &
       tensor_grid)
    sll_int32,                        intent(in)   :: local_sizes(6)
    sll_int32,                        intent(in)   :: indices_min(6)
    sll_real64,                       intent(in)   :: eta_min(6)
    sll_real64,                       intent(in)   :: delta_eta(6)
    type(sll_t_array),                intent(out)  :: tensor_grid(6)

    sll_int32 :: j, k, ierr

    do j=1,6
       allocate(tensor_grid(j)%vals(local_sizes(j)), stat=ierr)
       SLL_ASSERT(ierr == 0 )
       do k=1,local_sizes(j)
          tensor_grid(j)%vals(k) = eta_min(j) + &
               delta_eta(j) * real(indices_min(j)-2+k, f64)
       end do
    end do

  end subroutine sll_s_set_local_grid

  
  
  subroutine sll_s_set_local_grid_en(&
       local_sizes, &
       indices_min, &
       eta_min, &
       delta_eta, &
       tensor_grid)
    sll_int32,                        intent(in)   :: local_sizes(6)
    sll_int32,                        intent(in)   :: indices_min(6)
    sll_real64,                       intent(in)   :: eta_min(6)
    sll_real64,                       intent(in)   :: delta_eta(6)
    type(sll_t_array),                intent(out)  :: tensor_grid(6)
    sll_real64 :: vtrans

    sll_int32 :: j, k, ierr

    do j=1,6
       allocate(tensor_grid(j)%vals(local_sizes(j)), stat=ierr)
       SLL_ASSERT(ierr == 0 )
       do k=1,local_sizes(j)
          if (j .gt. 3) then
            call sll_s_compute_velocity_transformation_en(eta_min(j) + &
               delta_eta(j) * real(indices_min(j)-2+k, f64), vtrans)
            tensor_grid(j)%vals(k) = vtrans
          else
            tensor_grid(j)%vals(k) = eta_min(j) + &
               delta_eta(j) * real(indices_min(j)-2+k, f64)
          endif
       end do
    end do

  end subroutine sll_s_set_local_grid_en


  !< Function to initialized the parameters (read in from file)
  subroutine sll_s_distribution_params_6d_new( params, distrib_type, file_id )
    class( sll_c_distribution_params_6d ), allocatable, intent( out ) :: params
    character(*), intent(in) :: distrib_type
    sll_int32 :: file_id

    select case ( distrib_type )
    case ( "landau_sum" )
       allocate( sll_t_landau_sum_parameters_6d :: params )
       params%distrib_type = sll_p_landau_sum
    case ( "landau_prod" )
       allocate( sll_t_landau_prod_parameters_6d :: params )
       params%distrib_type = sll_p_landau_prod
    case ( "landau_diag" )
       allocate( sll_t_landau_diag_parameters_6d :: params )
       params%distrib_type = sll_p_landau_diag
    case ( "landau_sum_df" )
       allocate( sll_t_landau_sum_df_parameters_6d :: params )
       params%distrib_type = sll_p_landau_sum_df
    case ( "pslab" )
       allocate( sll_t_pslab_parameters_6d :: params )
       params%distrib_type = sll_p_pslab
    case ( "delta" )
       allocate( sll_t_delta_parameters_6d :: params )
       params%distrib_type = sll_p_delta
    case ( "pslab2" )
       allocate( sll_t_pslab2_parameters_6d :: params )
       params%distrib_type = sll_p_pslab2
    case ( "twogaussian_sum" )
       allocate( sll_t_twogaussian_parameters_6d :: params )
       params%distrib_type = sll_p_twogaussian_sum
    case default
       SLL_ERROR('Type of initial distribution not implemented.', 'sll_s_distribution_params_6d_new')
    end select

    call params%init( file_id )
    
  end subroutine sll_s_distribution_params_6d_new
  
  
  !> Initialize distribution function with given distribution parameter
  subroutine sll_s_distribution_initializer_6d( &
       local_sizes, &
       data_indices_min, &
       params, &
       tensor_grid, &
       fdistrib)
    sll_int32,                         intent(in)   :: local_sizes(6)
    sll_int32,                         intent(in)   :: data_indices_min(6)
    class(sll_c_distribution_params_6d),intent(in)   :: params
    type(sll_t_array),                 intent(in)   :: tensor_grid(6)
    sll_real64,                        intent(out)  :: fdistrib(1:,1:,1:,1:,1:,1:)
    sll_int32 :: i,j,k,l,m,n


!$omp parallel default(shared) private(i,j,k,l,m,n)
!$omp do OMP_COLLAPSE OMP_SCHEDULE
    do n=1, local_sizes(6)
       do m=1, local_sizes(5)
          do l=1, local_sizes(4)
             do k=1, local_sizes(3)
                do j=1, local_sizes(2)
                   do i=1, local_sizes(1)
                      fdistrib(i,j,k,l,m,n) = &
                           params%eval( [tensor_grid(1)%vals(i), tensor_grid(2)%vals(j), tensor_grid(3)%vals(k),  tensor_grid(4)%vals(l), tensor_grid(5)%vals(m), tensor_grid(6)%vals(n)])
                   end do
                end do
             end do
          end do
       end do
    end do
!$omp end do
!$omp end parallel
  end subroutine sll_s_distribution_initializer_6d
  
  
  
  
    subroutine sll_s_compute_velocity_transformation(vin, vtrans)
      sll_real64,    intent(in) :: vin
      sll_real64, intent(out) :: vtrans
      
      sll_real64    :: vmax  = 6.
      sll_real64    :: E = exp(1.)
      sll_real64    :: PI=4.D0*DATAN(1.D0)

      
!      
!       vtrans = 1.6766955862167576*((10*vmax*Sin((11.242085967477433*vin)/vmax**2))/(7.*Pi) &
!                         - (5*vmax*Sin((22.484171934954865*vin)/vmax**2))/(21.*Pi) &
!                         + (10*vmax*Sin((33.726257902432295*vin)/vmax**2))/(693.*Pi) &
!                         + (5*vmax*Sin((44.96834386990973*vin)/vmax**2))/(6006.*Pi) &
!                         + (2*vmax*Sin((56.21042983738717*vin)/vmax**2))/(15015.*Pi))
      
      vtrans =         (vmax*(211629600.*Sin((Pi*vin)/vmax) &
      - 79361100.*Sin((2*Pi*vin)/vmax) &
      + 32558400.*Sin((3*Pi*vin)/vmax) &
      - 12209400.*Sin((4*Pi*vin)/vmax) &
      + 3907008.*Sin((5*Pi*vin)/vmax)  &
      - 4275.*(238.*Sin((6*Pi*vin)/vmax)&
      - 48.*Sin((7.*Pi*vin)/vmax)&
      + 7.*Sin((8.*Pi*vin)/vmax)) &
      + 2800.*Sin((9.*Pi*vin)/vmax)&
      - 126.*Sin((10.*Pi*vin)/vmax)))/(1.1639628e8*Pi)

!       1.4909655684059773*((5.*vmax*Sin((12.642515911142883*vin)/vmax**2))/(3.*Pi) &
!                                - (10.*vmax*Sin((25.285031822285767*vin)/vmax**2))/(21.*Pi)&
!                                + (5.*vmax*Sin((37.927547733428646*vin)/vmax**2))/(42.*Pi) &
!                                - (5.*vmax*Sin((50.57006364457153*vin)/vmax**2))/(252.*Pi)&
!                                + (vmax*Sin((63.21257955571442*vin)/vmax**2))/(630.*Pi))

  end subroutine
  
  
  subroutine sll_s_compute_velocity_transformation_en(vin, vtrans)
      sll_real64,    intent(in) :: vin
      sll_real64, intent(out) :: vtrans
      
      sll_real64    :: vmax  = 6.
      sll_real64    :: E = exp(1.)
      sll_real64    :: PI=4.D0*DATAN(1.D0)

      
     
      vtrans = vin!1.6766955862167576*((10*vmax*Sin((11.242085967477433*vin)/vmax**2))/(7.*Pi) &
!                         - (5*vmax*Sin((22.484171934954865*vin)/vmax**2))/(21.*Pi) &
!                         + (10*vmax*Sin((33.726257902432295*vin)/vmax**2))/(693.*Pi) &
!                         + (5*vmax*Sin((44.96834386990973*vin)/vmax**2))/(6006.*Pi) &
!                         + (2*vmax*Sin((56.21042983738717*vin)/vmax**2))/(15015.*Pi))
      
      
!       vtrans = Sqrt(vmax**2*(196163108820. - 62249523265.*Cos((Pi*vin)/vmax) &
!             + 16481327830.*Cos((2*Pi*vin)/vmax) - 4782807075.*Cos((3*Pi*vin)/vmax) &
!             + 1302114020.*Cos((4*Pi*vin)/vmax) - 304886885.*Cos((5*Pi*vin)/vmax) &
!             + 57244242.*Cos((6*Pi*vin)/vmax) - 7953631.*Cos((7*Pi*vin)/vmax) &
!             + 720496.*Cos((8*Pi*vin)/vmax) - 31752.*Cos((9*Pi*vin)/vmax))&
!             *Sin((Pi*vin)/(2.*vmax))**2)/(630.*Sqrt(92378.)*Pi)
!       
!       1.4909655684059773*((5.*vmax*Sin((12.642515911142883*vin)/vmax**2))/(3.*Pi) &
!                                - (10.*vmax*Sin((25.285031822285767*vin)/vmax**2))/(21.*Pi)&
!                                + (5.*vmax*Sin((37.927547733428646*vin)/vmax**2))/(42.*Pi) &
!                                - (5.*vmax*Sin((50.57006364457153*vin)/vmax**2))/(252.*Pi)&
!                                + (vmax*Sin((63.21257955571442*vin)/vmax**2))/(630.*Pi))

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now the section defining the various types of paramters that are currently available.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! First the init functions

  subroutine init_landau_sum( self, file_id )
    class( sll_t_landau_sum_parameters_6d ), intent( out ) :: self
    sll_int32 :: file_id

    sll_real64 :: v_thermal(3)
    sll_real64 :: alpha(3)
    sll_real64 :: kx(3)
    sll_real64 :: phase(3) = [0.0_f64,0.0_f64,0.0_f64]
    
    namelist /landau_params/ v_thermal, alpha, kx, phase

    read(file_id, landau_params)

    self%v_thermal = v_thermal
    self%kx = kx
    self%alpha = alpha
    self%phase = phase
    self%factor = 1.0_f64/(sqrt(sll_p_twopi)**3*product(self%v_thermal))

  end subroutine init_landau_sum

  subroutine init_landau_prod ( self, file_id )
    class( sll_t_landau_prod_parameters_6d ), intent( out ) :: self
    sll_int32 :: file_id

    sll_real64 :: v_thermal(3)
    sll_real64 :: alpha
    sll_real64 :: kx(3)
    
    namelist /landau_params/ v_thermal, alpha, kx
    
    read(file_id, landau_params)

    self%v_thermal = v_thermal
    self%kx = kx
    self%alpha = alpha
    self%factor = 1.0_f64/(sqrt(sll_p_twopi)**3*product(self%v_thermal))

  end subroutine init_landau_prod
  

  subroutine init_twogaussian ( self, file_id )
    class( sll_t_twogaussian_parameters_6d ), intent( out ) :: self
    sll_int32 :: file_id

    sll_real64 :: v_thermal1(3)
    sll_real64 :: v_thermal2(3)
    sll_real64 :: v_mean1(3)
    sll_real64 :: v_mean2(3)
    sll_real64 :: portion
    sll_real64 :: alpha
    sll_real64 :: kx(3)
    
    namelist /twogaussian_params/ v_thermal1, v_thermal2, v_mean1, v_mean2, portion, alpha, kx

    read(file_id, twogaussian_params)
    
    self%v_thermal(1,:) = v_thermal1
    self%v_thermal(2,:) = v_thermal2
    self%v_mean(1,:) = v_mean1
    self%v_mean(2,:) = v_mean2
    self%kx = kx
    self%alpha = alpha
!    self%delta = [1-portion, portion]
    self%delta = [(1-portion)/sqrt(product(v_thermal1)), &
         portion/sqrt(product(v_thermal2))]
    self%factor = 1.0_f64/(sqrt(sll_p_twopi)**3)
    
  end subroutine init_twogaussian


  subroutine init_pslab ( self, file_id )
    class( sll_t_pslab_parameters_6d ), intent( out ) :: self
    sll_int32 :: file_id

    sll_real64 :: kappa_ti
    sll_real64 :: c_te
    sll_real64 :: b0
    sll_real64 :: v_thermal(3) = [1.0_f64, 1.0_f64, 1.0_f64]
    sll_real64 :: v_mean(3) = [0.0_f64, 0.0_f64, 0.0_f64] 
    sll_real64 :: alpha
    sll_real64 :: kx(3)
    
    namelist /landau_params/ v_thermal, alpha, kx, v_mean
    namelist /pslab_params/ kappa_ti, b0, c_te

    read(file_id, landau_params)
    read(file_id, pslab_params)
    
    self%kappa_ti = kappa_ti
    self%c_te = c_te
    self%b0 = b0
    self%alpha = alpha
    self%v_thermal = v_thermal
    self%v_mean = v_mean
    self%kx = kx
    self%factor =  1.0_f64/(sqrt(sll_p_twopi)**3)
    
  end subroutine init_pslab
  
  

  subroutine init_delta ( self, file_id )
    class( sll_t_delta_parameters_6d ), intent( out ) :: self
    sll_int32 :: file_id

    sll_real64 :: kappa_ti
    sll_real64 :: c_te
    sll_real64 :: b0
    sll_real64 :: v_thermal(3) = [1.0_f64, 1.0_f64, 1.0_f64]
    sll_real64 :: v_mean(3) = [0.0_f64, 0.0_f64, 0.0_f64] 
    sll_real64 :: alpha
    sll_real64 :: kx(3)
    
    namelist /landau_params/ v_thermal, alpha, kx, v_mean
    namelist /pslab_params/ kappa_ti, b0, c_te

    read(file_id, landau_params)
    read(file_id, pslab_params)
    
    self%kappa_ti = kappa_ti
    self%c_te = c_te
    self%b0 = b0
    self%alpha = alpha
    self%v_thermal = v_thermal
    self%v_mean = v_mean
    self%kx = kx
    self%factor =  1.0_f64/(sqrt(sll_p_twopi)**3)
    
  end subroutine init_delta

  ! Now the eval functions
  function eval_landau_sum ( self, x ) result(res)
    class( sll_t_landau_sum_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res
    
    res =  self%factor*(1.0_f64 + &
         self%alpha(1)*cos(self%kx(1)*x(1))+ &
         self%alpha(2)*cos(self%kx(2)*x(2))+ &
         self%alpha(3)*cos(self%kx(3)*x(3))) *&
         exp(-0.5_f64*((x(4)/self%v_thermal(1))**2+ &
         (x(5)/self%v_thermal(2))**2+ &
         (x(6)/self%v_thermal(3))**2))

  end function eval_landau_sum

  
  function eval_landau_sum_df ( self, x ) result(res)
    class( sll_t_landau_sum_df_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res
    
    res =  self%factor*( &
         self%alpha(1)*cos(self%kx(1)*x(1))+ &
         self%alpha(2)*cos(self%kx(2)*x(2))+ &
         self%alpha(3)*cos(self%kx(3)*x(3))) *&
         exp(-0.5_f64*((x(4)/self%v_thermal(1))**2+ &
         (x(5)/self%v_thermal(2))**2+ &
         (x(6)/self%v_thermal(3))**2))

  end function eval_landau_sum_df
  
  function eval_landau_prod ( self, x ) result(res)
    class( sll_t_landau_prod_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res
    
    res =  self%factor*(1.0_f64 + self%alpha*( &
         cos(self%kx(1)*x(1))* &
         cos(self%kx(2)*x(2))* &
         cos(self%kx(3)*x(3)))) *&
         exp(-0.5_f64*((x(4)/self%v_thermal(1))**2+ &
         (x(5)/self%v_thermal(2))**2+ &
         (x(6)/self%v_thermal(3))**2))

  end function eval_landau_prod
  
  function eval_landau_v ( self, x ) result(res)
    class( sll_t_landau_prod_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(3)
    sll_real64 :: res
    
    res =  self%factor* &
         exp(-0.5_f64*((x(1)/self%v_thermal(1))**2+ &
         (x(2)/self%v_thermal(2))**2+ &
         (x(3)/self%v_thermal(3))**2))

  end function eval_landau_v
  
  function eval_landau_v_sum ( self, x ) result(res)
    class( sll_t_landau_sum_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(3)
    sll_real64 :: res
    
    res =  self%factor* &
         exp(-0.5_f64*((x(1)/self%v_thermal(1))**2+ &
         (x(2)/self%v_thermal(2))**2+ &
         (x(3)/self%v_thermal(3))**2))

  end function eval_landau_v_sum

  function eval_landau_diag ( self, x ) result(res)
    class( sll_t_landau_diag_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res

    res =  self%factor*(1.0_f64 + self%alpha*( &
         cos(self%kx(1)*x(1)+ &
         self%kx(2)*x(2)+ &
         self%kx(3)*x(3)))) *&
         exp(-0.5_f64*((x(4)/self%v_thermal(1))**2+ &
         (x(5)/self%v_thermal(2))**2+ &
         (x(6)/self%v_thermal(3))**2))
 
  end function eval_landau_diag



  function eval_twogaussian_sum ( self, x ) result(res)
    class( sll_t_twogaussian_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res
    
    res =  self%factor*(1.0_f64 + self%alpha*( &
         cos(self%kx(1)*x(1))+ &
         cos(self%kx(2)*x(2))+ &
         cos(self%kx(3)*x(3)))) *&
         ( self%delta(1) * exp(-0.5_f64*(&
         ((x(4)-self%v_mean(1,1))/self%v_thermal(1,1))**2+ &
         ((x(5)-self%v_mean(1,2))/self%v_thermal(1,2))**2+ &
         ((x(6)-self%v_mean(1,3))/self%v_thermal(1,3))**2)) + &
         self%delta(2) * exp(-0.5_f64*(&
         ((x(4)-self%v_mean(2,1))/self%v_thermal(2,1))**2+ &
         ((x(5)-self%v_mean(2,2))/self%v_thermal(2,2))**2+ &
         ((x(6)-self%v_mean(2,3))/self%v_thermal(2,3))**2)))

  end function eval_twogaussian_sum

  function eval_twogaussian_v ( self, x ) result(res)
    class( sll_t_twogaussian_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(3)
    sll_real64 :: res
    
    res =  self%factor *&
         ( self%delta(1) * exp(-0.5_f64*(&
         ((x(1)-self%v_mean(1,1))/self%v_thermal(1,1))**2+ &
         ((x(2)-self%v_mean(1,2))/self%v_thermal(1,2))**2+ &
         ((x(3)-self%v_mean(1,3))/self%v_thermal(1,3))**2)) + &
         self%delta(2) * exp(-0.5_f64*(&
         ((x(1)-self%v_mean(2,1))/self%v_thermal(2,1))**2+ &
         ((x(2)-self%v_mean(2,2))/self%v_thermal(2,2))**2+ &
         ((x(3)-self%v_mean(2,3))/self%v_thermal(2,3))**2)))

  end function eval_twogaussian_v

  function eval_pslab ( self, x ) result(res)
    class( sll_t_pslab_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res
    
!     res =  self%factor* &
!          (1.0_f64 + self%alpha* &
!          cos(self%kx(1)*x(1) + &
!          self%kx(2)*x(2) + &
!          self%kx(3)*x(3))* &
!          exp(-0.5_f64*((x(4)- self%v_mean(1))**2+ &
!          (x(5)-self%v_mean(2))**2+ &
!          (x(6)-self%v_mean(3))**2)))

    if (self%B0 .ne. 0)then
    
        res =  self%factor* &
            (1.0_f64 + self%alpha* &
            (cos(self%kx(1)*(x(1)- &
            x(5)/self%B0) + &
            self%kx(2)*(x(2) + &
            x(4)/self%B0) + &
            self%kx(3)*x(3)) - &
            0.5_f64*cos(self%kx(1)*x(1) + &
            self%kx(2)*x(2) + &
            self%kx(3)*x(3))* &
            exp(-0.5_f64*(self%kx(1)**2+self%kx(2)**2)/self%B0**2)))*&
            exp(-0.5_f64*((x(4)- self%v_mean(1))**2+ &
            (x(5)-self%v_mean(2))**2+ &
            (x(6)-self%v_mean(3))**2))
    else
        res =  self%factor* &
            (1.0_f64 + self%alpha* &
            sin(self%kx(1)*(x(1)) + &
            self%kx(2)*(x(2)) + &
            self%kx(3)*x(3))) *&
            exp(-0.5_f64*((x(4)- self%v_mean(1))**2+ &
            (x(5)-self%v_mean(2))**2+ &
            (x(6)-self%v_mean(3))**2))
    endif
    !          
         
!         res =  self%factor* &
!          (1.0_f64 + self%alpha* &
!          exp(-(x(1)-x(5)/self%B0-sll_p_twopi/self%kx(1)/2)**2 * self%kx(1)**2)*&
!          exp(-(x(2)+x(4)/self%B0-sll_p_twopi/self%kx(2)/2)**2 * self%kx(2)**2)*&
!          exp(-(x(3)-sll_p_twopi/self%kx(3)/2)**2 * self%kx(3)**2)- &
!         0.5/(exp(((3 + 4*self%kx(2)**2 + 4*self%kx(1)**2*(1 + self%kx(2)**2))*(sll_p_twopi/2.0)**2 + self%kx(2)**2*x(2)**2 + self%kx(3)**2*x(3)**2 +    2*self%kx(2)**2*self%kx(3)**2*x(3)**2 - 2*(sll_p_twopi/2.0)* (self%kx(1)*(x(1) + 2*self%kx(2)**2*x(1)) + self%kx(2)*x(2) + self%kx(3)*x(3) + 2*self%kx(2)**2*self%kx(3)*x(3) + 2*self%kx(1)**2*(self%kx(2)*x(2) + self%kx(3)*x(3) + 2*self%kx(2)**2*self%kx(3)*x(3))) +self%kx(1)**2*((1 + 2*self%kx(2)**2)*x(1)**2 + 2*self%kx(3)**2*x(3)**2 + 2*self%kx(2)**2*(x(2)**2 + 2*self%kx(3)**2*x(3)**2)))/((1 + 2*self%kx(1)**2)*(1 + 2*self%kx(2)**2)))*Sqrt(1 + 2*self%kx(1)**2)*Sqrt(1 + 2*self%kx(2)**2)))*&
!          exp(-0.5_f64*((x(4)- self%v_mean(1))**2+ &
!          (x(5)-self%v_mean(2))**2+ &
!          (x(6)-self%v_mean(3))**2))     
    
  end function eval_pslab
  
  function eval_delta ( self, x ) result(res)
    class( sll_t_delta_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res,r 
     
    res = 1.0
    
    r = self%randArray(modulo(int(x(1)*241+73)*int(x(2)*163+211)*int(x(3)*31+139),509)+1)
    
    
    res =   self%factor * (res+self%alpha*r)*exp(-0.5_f64*((x(4)- self%v_mean(1))**2+ &
                        (x(5)-self%v_mean(2))**2+ &
                        (x(6)-self%v_mean(3))**2))
    
    
  end function eval_delta


  function eval_pslab2 ( self, x ) result(res)
    class( sll_t_pslab2_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(6)
    sll_real64 :: res,r
    sll_int32 :: l,m,n

!     print*, self%factor, self%alpha, self%kx, x, self%B, 
!     
!     res =self%factor* &
!          (1.0_f64 + self%alpha* &
!          (cos(self%kx(1)*(x(1)- &
!          x(5)/self%B0) + &
!          self%kx(2)*(x(2)+ &
!          x(4)/self%B0) + &
!          self%kx(3)*x(3)) - &
!          0.5_f64*cos(self%kx(1)*(x(1)) + &
!          self%kx(2)*(x(2)) + &
!          self%kx(3)*x(3))* &
!          exp(-0.5_f64*(self%kx(1)**2+self%kx(2)**2)/self%B0**2)))*&
!          exp(-0.5_f64*((x(4)-self%v_mean(1))**2+ &
!          (x(5)-self%v_mean(2))**2+ &
!          (x(6)-self%v_mean(3))**2)) + &
!          0.1_f64*self%factor* &
!          self%alpha* &
!          (cos(self%kx(2)*(x(2) + &
!          x(4)/self%B0) + &
!          self%kx(3)*x(3)) - &
!          0.5_f64*cos(self%kx(2)*(x(2)) + &
!          self%kx(3)*x(3))* &
!          exp(-0.5_f64*(self%kx(2)**2)/self%B0**2))*&
!          exp(-0.5_f64*((x(4)-self%v_mean(1))**2+ &
!          (x(5)-self%v_mean(2))**2+ &
!          (x(6)-self%v_mean(3))**2))
!          
!          
        res = 0.0
         
    do m=-4,4
       do n=-4,4
           do l=-4,4
               if (n .ne. 0 .or. l .ne. 0 .or. m .ne. 0 ) then
                   r = self%randArray(modulo((m+8)*(n+8)*(l+8),512))
                   
                   res = res +  r* 0.1*exp(-0.5_f64*(n**2+l**2+m**2) * 0.05) *&
                       self%alpha* &
                       (cos(l*self%kx(1)*(x(1)- &
                       x(5)/self%B0) + &
                       m*self%kx(2)*(x(2) + &
                       x(4)/self%B0) + &
                       n*self%kx(3)*x(3))- &
                       0.5_f64*cos(l*self%kx(1)*(x(1)) + &
                       m*self%kx(2)*(x(2)) + &
                       n*self%kx(3)*x(3))* &
                       exp(-0.5_f64*(l**2*self%kx(1)**2+m**2*self%kx(2)**2)/self%B0**2)    )                   
                   end if
               end do
           end do 
       end do     
        
                        
    res = self%factor*(1.0_f64 + self%alpha*(res))*&
            exp(-0.5_f64*(  (x(4)-self%v_mean(1))**2+ &
                            (x(5)-self%v_mean(2))**2+ &
                            (x(6)-self%v_mean(3))**2)) 
         
         
         

  end function eval_pslab2

    function eval_pslab_v ( self, x ) result(res)
    class( sll_t_pslab_parameters_6d ), intent( in ) :: self
    sll_real64,                               intent(in)   :: x(3)
    sll_real64 :: res
    
    res =  self%factor* &
         exp(-0.5_f64*((x(1)- self%v_mean(1))**2+ &
         (x(2)-self%v_mean(2))**2+ &
         (x(3)-self%v_mean(3))**2))

  end function eval_pslab_v

  ! TODO: Old function, needs to be adapted.
!!$  !> Initialize distribution function with Landau initial value
!!$  subroutine sll_s_initialize_itg_6d( &
!!$       local_sizes, &
!!$       indices_min, &
!!$       eta_min, &
!!$       delta_eta, &
!!$       params, &
!!$       fdistrib, &
!!$       tensor_grid, &
!!$       Terdn0r)
!!$    sll_int32,                        intent(in)   :: local_sizes(6)
!!$    sll_int32,                        intent(in)   :: indices_min(6)
!!$    sll_real64,                       intent(in)   :: eta_min(6)
!!$    sll_real64,                       intent(in)   :: delta_eta(6)
!!$    type(sll_t_itg_parameters_6d), intent(in)   :: params
!!$    sll_real64,                       intent(out)  :: fdistrib(1:,1:,1:,1:,1:,1:)
!!$    type(sll_t_array),                intent(out)  :: tensor_grid(6)
!!$    sll_real64,                       intent(out)  :: Terdn0r(:)
!!$
!!$    sll_int32 :: i,j,k,l,m,n
!!$    sll_real64 :: factor, x, cn0(1), cn0_local(1)
!!$    sll_real64, allocatable :: n0r(:), Ter(:), Tir(:)
!!$
!!$    factor    = 1.0_f64/(sqrt(sll_p_twopi)**3*product(params%v_thermal))
!!$
!!$    call sll_s_set_local_grid(local_sizes, indices_min, eta_min, delta_eta, tensor_grid)
!!$
!!$
!!$    allocate(n0r(local_sizes(1)))
!!$    allocate(Ter(local_sizes(1)))
!!$    allocate(Tir(local_sizes(1)))
!!$
!!$
!!$    do i=1, local_sizes(1)
!!$       x = tensor_grid(1)%vals(i)
!!$       n0r(i) = exp(-params%kappa_n0*params%delta_rn0* &
!!$            tanh((x-params%rp)/params%delta_rn0))
!!$       Ter(i) = params%C_Te* exp(-params%kappa_Te*params%delta_rTe* &
!!$            tanh((x-params%rp)/params%delta_rTe))
!!$       Tir(i) = exp(-params%kappa_Ti*params%delta_rTi* &
!!$            tanh((x-params%rp)/params%delta_rTi))
!!$       Terdn0r(i) = Ter(i)/n0r(i)
!!$    end do
!!$    cn0_local = sum(n0r)
!!$    ! TODO: Exchange
!!$    call sll_o_collective_allreduce( sll_v_world_collective, cn0_local, 1, MPI_SUM, cn0 ) 
!!$
!!$    n0r = n0r/cn0
!!$    Terdn0r = Terdn0r*cn0
!!$
!!$
!!$    do n=1, local_sizes(6)
!!$       do m=1, local_sizes(5)
!!$          do l=1, local_sizes(4)
!!$             do k=1, local_sizes(3)
!!$                do j=1, local_sizes(2)
!!$                   do i=1, local_sizes(1)
!!$                      ! TODO: Update
!!$                      fdistrib(i,j,k,l,m,n) = &
!!$                           n0r(i)/(sll_p_twopi)**(1.5_f64)* &
!!$
!!$                           factor*(1.0_f64 + params%alpha* &
!!$                           cos(params%kx(1)*tensor_grid(1)%vals(i))* &
!!$                           cos(params%kx(2)*tensor_grid(2)%vals(j))* &
!!$                           cos(params%kx(3)*tensor_grid(3)%vals(k))) *&
!!$                           exp(-0.5_f64/Tir(i)*((tensor_grid(4)%vals(l)/params%v_thermal(1))**2+ &
!!$                           (tensor_grid(5)%vals(m)/params%v_thermal(2))**2+ &
!!$                           (tensor_grid(6)%vals(n)/params%v_thermal(3))**2))
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine sll_s_initialize_itg_6d

  
  
end module sll_m_distribution_function_initializer_6d
