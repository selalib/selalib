module sim_bsl_vp_2d2v_cart_multipatch_helper
#include "sll_working_precision.h"

contains
  
  function func_one( eta1, eta2, params ) result(res)
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = 1.0_8
  end function func_one
  
  function func_minus_one( eta1, eta2, params ) result(res)
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = -1.0_8
  end function func_minus_one
  
  function func_zero( eta1, eta2, params ) result(res)
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = 0.0_8
  end function func_zero
  
  
  function func_epsi( eta1, eta2, params ) result(res)
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = 0.00001_8
  end function func_epsi
  
  function electric_field_ext_1(x,y,params) result(res)
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = params(1)*x
  end function electric_field_ext_1
  
  
  function electric_field_ext_2(x,y,params) result(res)
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: res
    res = params(2)*y
  end function electric_field_ext_2
  
end module sim_bsl_vp_2d2v_cart_multipatch_helper


! External functions used as parameters in the above unit test:





  
! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with an arbitrary coordinate 
!   transformation in the x,y variables.
! - parallel
! - The coordinate transformation is defined by patches.

program sim_bsl_vp_2d2v_cart_multipatch
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_halt_collective

  use sll_m_common_array_initializers, only: &
    sll_f_gaussian_beam_initializer_4d

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_multipatch, only: &
    sll_t_coordinate_transformation_multipatch_2d

  use sll_m_general_coordinate_elliptic_solver, only: &
    sll_p_es_gauss_legendre

  use sll_m_sim_bsl_vp_2d2v_cart_multipatch, only: &
    sll_s_initialize_4d_qns_gen_mp, &
    sll_s_run_4d_qns_general_mp, &
    sll_t_simulation_4d_qns_general_multipatch

  use sim_bsl_vp_2d2v_cart_multipatch_helper

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_t_simulation_4d_qns_general_multipatch)     :: simulation
  type(sll_t_cartesian_mesh_2d), pointer      :: mv ! delete
  type(sll_t_coordinate_transformation_multipatch_2d) :: mp

  sll_real64, dimension(1:8) :: gaussian_beam_params
  sll_real64, dimension(1:2) :: elec_field_ext_params
  sll_real64, dimension(1) :: f_zero_params
  sll_real64, dimension(1) :: f_one_params
  sll_real64, dimension(1) :: f_minus_one_params
  sll_real64, dimension(1) :: f_epsi_params


  print *, 'Booting parallel environment...'
  call sll_s_boot_collective() ! Wrap this up somewhere else


  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)
  
  f_zero_params(:) = (/0.0_f64/)
  f_one_params(:) = (/0.0_f64/)
  f_minus_one_params(:) = (/0.0_f64/)
  f_epsi_params(:) = (/0.0_f64/)
  elec_field_ext_params(:) =  (/64.0_f64,64.0_f64/)!(/36.0_f64,36.0_f64/)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:


  
  call simulation%init_from_file(filename_local)
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

! hardwired, this should be consistent with whatever is read from a file
#define NPTS3 64
#define NPTS4 64
#define SPL_DEG_VX 3 
#define SPL_DEG_VY 3


  ! logical mesh for velocity coordinates
  mv => sll_f_new_cartesian_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-9.0_f64, eta1_max=9.0_f64, &
       eta2_min=-9.0_f64, eta2_max=9.0_f64)

  ! ---------------------------------------------------------------------
  ! coordinate transformation associated with space coordinates
  ! ---------------------------------------------------------------------

  call mp%read_from_file("circle_mp5_pts12")

  ! ---------------------------------------------------------------------
  ! define the values of the parameters for the sll_m_gaussian beam initializer
  ! ---------------------------------------------------------------------

  ! sll_f_gaussian_beam_initializer_4d parameters
  !TODO : fix parameteres

  gaussian_beam_params(1) = 1.0_f64!3.0_f64/2.0_f64!vth
  gaussian_beam_params(2) = 0.05_f64!1.0_f64!xth
  gaussian_beam_params(3) = 0.0_f64!vxc
  gaussian_beam_params(4) = 0.0_f64!vyc
  gaussian_beam_params(5) = 0.0_f64!xc
  gaussian_beam_params(6) = 0.0_f64!yc
  gaussian_beam_params(7) = sll_p_pi*15*8!sll_p_pi*15*18!n0
  gaussian_beam_params(8) = 6.0_f64 !radius
  

  ! ---------------------------------------------------------------------
  ! initialize simulation object with the above parameters
  ! ---------------------------------------------------------------------

  call simulation%initialize( &
       mv, &
       mp, &
       sll_f_gaussian_beam_initializer_4d, &
       gaussian_beam_params, &
       func_one,  &  ! a11
       f_one_params, &
       func_zero, &   !a12
       f_zero_params, &
       func_zero, &   !a21
       f_zero_params, &
       func_one,  &   !a22
       f_one_params, &
       func_zero, &   !b1
       f_zero_params, &
       func_zero, &   !der1 b1
       func_zero, &   !der2 b1
       func_zero, &   ! b2
       f_zero_params, &
       func_zero, &   !der1 b2
       func_zero, &   !der2 b2
       func_zero, &   ! c
       f_zero_params, &
       SPL_DEG_VX, & 
       SPL_DEG_VY, & 
       sll_p_es_gauss_legendre,&
       sll_p_es_gauss_legendre,&
       func_zero, &
       func_zero, &
       elec_field_ext_params,&
       100)


  print *, ' f initialized '

  call simulation%run( )
 ! call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_s_halt_collective()

  !call delete(simulation)
  !call delete(transformation_x)

end program sim_bsl_vp_2d2v_cart_multipatch

