!------------------------------------------------------------------------------!
! Implementation of sgxsgv tensor product sparse grid semi-Lagrangian solver for
! 4D Vlasov-Poisson
! Implemented test cases: Landau damping (test_case = SLL_LANDAU) and
! TSI (test_case = SLL_TSI) where TSI along x1v1 and Landau along x2v2.
! Author: Katharina Kormann, IPP
!------------------------------------------------------------------------------

module sll_m_sim_bsl_vp_2d2v_cart_sparsegrid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_collective_alltoallv_double, &
    sll_s_collective_alltoallv_int, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_2d_sparse_grid_fft, only: &
    sll_t_fft_derivative

  use sll_m_remapper, only: &
    sll_o_apply_remap_2d, &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_2d_real64

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_sparse_grid_2d, only: &
    sll_t_sparse_grid_interpolator_2d

  implicit none

  public :: &
    sll_t_sim_sl_vp_2d2v_cart_sparsegrid

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_LANDAU = 0
  sll_int32, parameter :: SLL_TSI = 1

  
  type, extends(sll_c_simulation_base_class) :: sll_t_sim_sl_vp_2d2v_cart_sparsegrid
     logical     :: is_mdeltaf
     sll_int32   :: test_case
     sll_int32   :: levelx
     sll_int32   :: levelv
     sll_int32   :: order_sg

     ! Domain
     sll_real64  :: eta_min(4), eta_max(4)

     !Time domain
     sll_int32  :: n_time_steps
     sll_real64 :: delta_t

     !Sparse grid parameters
     sll_int32 :: level, order
     sll_int32, dimension(2) :: levelsx,levelsv, dorder1, dorder2

     !Distribution function 6D
     sll_real64 :: eps, v2, v0
     sll_real64, dimension(:,:), allocatable :: f_x, f_v,ft_v,f2_x,ft2_v
     
     !Electric fields and charge density
     sll_real64, dimension(:), allocatable :: ex
     sll_real64, dimension(:), allocatable :: ey
     sll_real64, dimension(:), allocatable :: ez
     sll_real64, dimension(:), allocatable :: rho, rho_loc
     sll_real64, dimension(:), allocatable :: dv

     !Poisson solver
     type(sll_t_fft_derivative)  :: poisson
     
     ! Interpolator
     type(sll_t_sparse_grid_interpolator_2d)   :: interp_x, interp_v

     ! For parallelization
     type(sll_t_remap_plan_2d_real64), pointer :: remap_x2v, remap_v2x ! Remapper between layoutx and layoutv
     type(sll_t_layout_2d), pointer :: layoutx, layoutv ! Mesh discribution
     sll_int32 :: colsz ! no of processors
     sll_int32 :: myrank ! mpi rank
     sll_int32, dimension(2) :: local_size_x, local_size_v ! Local size for x and v layouts, respectively
     sll_int32, dimension(2) :: global_x, global_v ! Index shift (add to get global index from local index
     sll_int32, dimension(:), allocatable :: recv_displs, send_displs, recv2, send_cnt, recv_cnt

     !Diagnostics
     sll_real64, dimension(:), allocatable :: nrj

     contains
       procedure :: init_from_file => init_sparsegrid_2d2v
       procedure :: run => run_sparsegrid_2d2v

    end type sll_t_sim_sl_vp_2d2v_cart_sparsegrid

contains
!------------------------------------------------------------------------------!
  subroutine init_sparsegrid_2d2v (sim, filename)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    character(len=*), intent(in)                               :: filename

    ! Local variables
    sll_int32 :: io_stat
    logical   :: is_mdeltaf
    character(len=256) :: test_case
    sll_int32   :: levelx
    sll_int32   :: levelv
    sll_int32   :: order_sg
    sll_real64  :: delta_t
    sll_int32   :: n_time_steps

    sll_int32, parameter :: input_file = 99

    namelist /sim_params/ is_mdeltaf, test_case, levelx, levelv, order_sg, delta_t, n_time_steps

    ! Read parameters from file
    open( unit = input_file, file=trim(filename), IOStat= io_stat)
    if (io_stat /= 0) then
       print*, 'init_sparsegrid_2d2v() failed to open file', filename
       STOP
    end if
    read(input_file, sim_params)
    close(input_file)

    sim%is_mdeltaf = is_mdeltaf
    sim%levelx = levelx
    sim%levelv = levelv
    sim%order_sg = order_sg
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps

    select case(test_case)
    case("SLL_LANDAU")
       sim%test_case = SLL_LANDAU
    case("SLL_TSI")
       sim%test_case = SLL_TSI
    case default
       print*, '#test case', test_case, 'not implemented.'
       STOP
    end select

  end subroutine init_sparsegrid_2d2v

!------------------------------------------------------------------------------!
  subroutine run_sparsegrid_2d2v(sim)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim

    sll_int32  :: j, i1, i2
    sll_int32  :: error
    sll_real64 :: eta(4)
    sll_real64 :: kxy(3)
    sll_real64 :: time
    sll_int32  :: i_step

    ! MPI parameters
    sim%colsz = sll_f_get_collective_size(sll_v_world_collective)
    sim%myrank = sll_f_get_collective_rank(sll_v_world_collective)

    !Sparse grid
    sim%levelsx(1) = sim%levelx; sim%levelsx(2) = sim%levelsx(1);
    sim%levelsv(1) = sim%levelv; sim%levelsv(2) = sim%levelsv(1);
    
    ! Information on ordering for 1D interpolations
    sim%dorder1 = (/1, 2 /); 
    sim%dorder2 = (/2, 1 /);

    ! Set the domain for one of the two test cases
    if (sim%test_case == SLL_LANDAU ) then
       !x domain
       sim%eta_min(1) =  0.0_f64; sim%eta_max(1) =  4.0_f64 * sll_p_pi
       sim%eta_min(2) = sim%eta_min(1); sim%eta_max(2) = sim%eta_max(1)
       !v domain
       sim%eta_min(3) = -6.0_f64; sim%eta_max(3) = 6.0_f64
       sim%eta_min(4) = sim%eta_min(3); sim%eta_max(4) = sim%eta_max(3)
    elseif (sim%test_case  == SLL_TSI) then
       !x domain
       sim%eta_min(1) =  0.0_f64; sim%eta_max(1) =  10.0_f64 * sll_p_pi
       sim%eta_min(2) =  0.0_f64; sim%eta_max(2) =  4.0_f64 * sll_p_pi
       !v domain
       sim%eta_min(3) = -8.0_f64; sim%eta_max(3) = 8.0_f64
       sim%eta_min(4) = -6.0_f64; sim%eta_max(4) = 6.0_f64
    end if

    ! Set the perturbation
    if (sim%test_case == SLL_LANDAU) then
       sim%eps = 0.01_f64
    elseif (sim%test_case == SLL_TSI) then
       sim%eps = 0.001_f64
       sim%v0 = 2.4_f64
    end if

    ! Initialize the sparse grids in x and v
    call sim%interp_x%initialize(sim%levelsx,sim%order_sg, sim%order_sg+1,0, &
         sim%eta_min(1:2),sim%eta_max(1:2), 0, 0);
    call sim%interp_v%initialize(sim%levelsv,sim%order_sg,sim%order_sg+1,0, &
         sim%eta_min(3:4), sim%eta_max(3:4),0, 1);
    ! Initialize Poisson
    call sim%poisson%initialize(sim%interp_x);


    ! Now we compute the two parallel layouts
    sim%layoutx => sll_f_new_layout_2d(sll_v_world_collective);
    call sll_o_initialize_layout_with_distributed_array(&
         sim%interp_x%size_basis, sim%interp_v%size_basis, 1, sim%colsz, sim%layoutx);
    call sll_o_compute_local_sizes(sim%layoutx, sim%local_size_x(1), sim%local_size_x(2));
    SLL_ALLOCATE(sim%f_x(sim%local_size_x(1),sim%local_size_x(2)),error);
    sim%global_x = sll_o_local_to_global(sim%layoutx, (/1,1/));
    sim%global_x = sim%global_x - 1;
    
    sim%layoutv => sll_f_new_layout_2d(sll_v_world_collective);
    call sll_o_initialize_layout_with_distributed_array(&
         sim%interp_x%size_basis, sim%interp_v%size_basis, sim%colsz, 1, sim%layoutv);
    call sll_o_compute_local_sizes(sim%layoutv, sim%local_size_v(1), sim%local_size_v(2));
    SLL_ALLOCATE(sim%f_v(sim%local_size_v(1),sim%local_size_v(2)),error);
    sim%global_v = sll_o_local_to_global(sim%layoutv, (/1,1/));
    sim%global_v = sim%global_v - 1;
    
    sim%remap_x2v => sll_o_new_remap_plan(sim%layoutx, sim%layoutv, sim%f_x);
    sim%remap_v2x => sll_o_new_remap_plan(sim%layoutv, sim%layoutx, sim%f_v);
    
    ! Information for rho communication
    SLL_ALLOCATE(sim%recv_displs(sim%colsz),error);
    SLL_ALLOCATE(sim%send_displs(sim%colsz),error);
    SLL_ALLOCATE(sim%recv2(sim%colsz),error);
    SLL_ALLOCATE(sim%send_cnt(sim%colsz),error);
    SLL_ALLOCATE(sim%recv_cnt(sim%colsz),error);
    sim%send_displs = 0;
    sim%recv2 = (/ (j, j=0,sim%colsz-1) /);
    sim%recv_displs = 0
    sim%send_cnt = 1;
    sim%recv_cnt = 1;
    call sll_s_collective_alltoallv_int((/sim%global_v(1)/), sim%send_cnt, sim%send_displs,&
         sim%recv_displs, sim%recv_cnt, sim%recv2, sll_v_world_collective);
    sim%send_cnt = sim%local_size_v(1);
    do j=1,sim%colsz-1
       sim%recv_cnt(j) = sim%recv_displs(j+1)-sim%recv_displs(j);
    end do
    sim%recv_cnt(sim%colsz) = sim%interp_x%size_basis-sim%recv_displs(sim%colsz);
    
    ! Print MPI information
    if (sim%myrank == 0) then
       print*, 'Parallel layout'
       print*, 'No. of processors:', sim%colsz
    end if
    print*, 'Size and start of x layout  of processor', sim%myrank, ':', sim%local_size_x, sim%global_x
    print*, 'Size and start of v layout  of processor', sim%myrank, ':', sim%local_size_v, sim%global_v


    ! Allocate the arrays
    SLL_ALLOCATE(sim%f2_x(sim%local_size_x(1),sim%local_size_x(2)),error);
    SLL_ALLOCATE(sim%ft_v(sim%local_size_v(2),sim%local_size_v(1)),error);
    SLL_ALLOCATE(sim%ft2_v(sim%local_size_v(2),sim%local_size_v(1)),error);
    SLL_ALLOCATE(sim%rho(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%rho_loc(sim%local_size_v(1)),error)
    SLL_ALLOCATE(sim%ex(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%ey(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%ez(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%dv(sim%local_size_v(1)),error)
    
    SLL_ALLOCATE(sim%nrj(0:sim%n_time_steps), error)


  
    ! Set initial value
    do j=1,2
       kxy(j)  = 2.0_f64*sll_p_pi/(sim%eta_max(j)-sim%eta_min(j))
    end do

    do i1=1,sim%local_size_x(1)
       eta(1:2) =  sim%interp_x%hierarchy(sim%global_x(1)+i1)%coordinate
       do i2=1,sim%local_size_x(2)
          eta(3:4) =sim%interp_v%hierarchy(sim%global_x(2)+i2)%coordinate
          sim%f_x(i1,i2) = 1.0_f64/(2.0_f64*sll_p_pi)*&
               (1.0_f64+sim%eps*(cos(kxy(1)*eta(1))+&
               cos(kxy(2)*eta(2))))
          if (sim%test_case == SLL_TSI) then
             if (sim%is_mdeltaf .EQV. .FALSE.) then
                sim%f_x(i1,i2) = sim%f_x(i1,i2)*0.5_f64*&
                     (exp(-(eta(3)-sim%v0)**2*0.5_f64) + &
                     exp(-(eta(3)+sim%v0)**2*0.5_f64))* &
                     exp(-(eta(4)**2)*0.5_f64)
             else
                sim%f_x(i1,i2) = sim%f_x(i1,i2)*0.5_f64*&
                     (exp(-(eta(3)-sim%v0)**2*0.5_f64) + &
                     exp(-(eta(3)+sim%v0)**2*0.5_f64))
             end if
          elseif ((sim%test_case == SLL_LANDAU) .AND. (sim%is_mdeltaf .EQV. .FALSE.)) then
             sim%f_x(i1,i2) = sim%f_x(i1,i2)*exp(-.5_f64*( eta(3)**2+eta(4)**2))
          end if
       end do
    end do

    time = 0.0_f64
    ! Propagate v half a step
    call sll_o_apply_remap_2d(sim%remap_x2v, sim%f_x, sim%f_v);
    sim%ft_v = transpose(sim%f_v);
    
    if (sim%is_mdeltaf .EQV. .FALSE.) then
       call compute_rho(sim)
    else
       call compute_rho_df(sim)
    end if
    call sim%poisson%solve(sim%interp_x,sim%rho,sim%ex,sim%ey)

    if (sim%is_mdeltaf .EQV. .FALSE. ) then
       call advection_v(sim, sim%delta_t)
    else
       call advection_v_df(sim, sim%delta_t)
    end if
    sim%f_v = transpose(sim%ft_v);
    call sll_o_apply_remap_2d(sim%remap_v2x, sim%f_v, sim%f_x);
    call compute_energy(sim, sim%nrj(0));
       
   
    if (sim%myrank == 0) then
       open(11, file='sgxsgv4d_output.dat', position='append')
       rewind(11)
       write(11,*) time, sim%nrj(0)
    end if

    do i_step = 1, sim%n_time_steps !Loop over time

       if (sim%myrank == 0) then
          if(mod(i_step,100) == 0) then
             print*, i_step
          end if
       end if
       
       call advection_x(sim, sim%delta_t)
       
       call sll_o_apply_remap_2d(sim%remap_x2v, sim%f_x, sim%f_v);
       sim%ft_v = transpose(sim%f_v);

       if (sim%is_mdeltaf .EQV. .FALSE.) then
          call compute_rho(sim)
       else
          call compute_rho_df(sim)
       end if
       call sim%poisson%solve(sim%interp_x,sim%rho,sim%ex,sim%ey)

       if (sim%is_mdeltaf .EQV. .FALSE. ) then
          call advection_v(sim, sim%delta_t)
       else
          call advection_v_df(sim, sim%delta_t)
       end if
       sim%f_v = transpose(sim%ft_v);
       call sll_o_apply_remap_2d(sim%remap_v2x, sim%f_v, sim%f_x);
       call compute_energy(sim, sim%nrj(i_step));
       
       time  = time + sim%delta_t
       
       if (sim%myrank == 0) then
          open(11, file='sgxsgv4d_output.dat', position='append')
          write(11,*) time, sim%nrj(i_step)
          close(11)
       end if
       
    end do !next time step


  end subroutine run_sparsegrid_2d2v
  
!------------------------------------------------------------------------------!
  subroutine advection_x(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt
    sll_int32 :: i2
    
    do i2 = 1, sim%local_size_x(2)
       call sim%interp_x%compute_hierarchical_surplus(sim%f_x(:,i2));
       call sim%interp_x%interpolate_const_disp(sim%dorder1,&
            -dt*sim%interp_v%hierarchy(sim%global_x(2)+i2)%coordinate(1),&
            sim%f_x(:,i2), sim%f2_x(:,i2),.TRUE.)
       call sim%interp_x%interpolate_const_disp(sim%dorder2,&
            -dt*sim%interp_v%hierarchy(sim%global_x(2)+i2)%coordinate(2),&
            sim%f2_x(:,i2), sim%f_x(:,i2),.FALSE.)
    end do

  end subroutine advection_x
!------------------------------------------------------------------------------!


  subroutine advection_v(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt
    sll_int32 :: i1

    do i1 = 1, sim%local_size_v(1) 
       call sim%interp_v%compute_hierarchical_surplus(sim%ft_v(:,i1));
       call sim%interp_v%interpolate_const_disp(sim%dorder1,&
            dt*sim%ex(sim%global_v(1)+i1),sim%ft_v(:,i1), sim%ft2_v(:,i1),.TRUE.);
       call sim%interp_v%interpolate_const_disp(sim%dorder2,&
            dt*sim%ey(sim%global_v(1)+i1),sim%ft2_v(:,i1), sim%ft_v(:,i1),.FALSE.);
    end do
    
  end subroutine advection_v

  subroutine advection_v_df(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt
    sll_int32 :: i1, i2
    sll_real64 :: vv

    if (sim%test_case == SLL_LANDAU) then 
       do i1 = 1, sim%local_size_v(1) 
          call sim%interp_v%compute_hierarchical_surplus(sim%ft_v(:,i1));
          call sim%interp_v%interpolate_const_disp(sim%dorder1,&
               dt*sim%ex(sim%global_v(1)+i1),sim%ft_v(:,i1), sim%ft2_v(:,i1),.FALSE.);
       end do
       do i2 = 1, sim%local_size_v(2)
          vv = sim%interp_v%hierarchy(sim%global_v(2)+i2)%coordinate(1)
          sim%dv = sim%ex(sim%global_v(1)+1:sim%global_v(1)+sim%local_size_v(1))*dt
          sim%ft_v(i2,:) = sim%ft2_v(i2,:)*exp(-vv*sim%dv-sim%dv**2*0.5_f64)
       end do
       do i1 = 1, sim%local_size_v(1) 
          call sim%interp_v%compute_hierarchical_surplus(sim%ft_v(:,i1));
          call sim%interp_v%interpolate_const_disp(sim%dorder2,&
               dt*sim%ey(sim%global_v(1)+i1),sim%ft_v(:,i1), sim%ft2_v(:,i1),.FALSE.);
       end do
    elseif (sim%test_case == SLL_TSI) then
        do i1 = 1, sim%local_size_v(1) 
          call sim%interp_v%compute_hierarchical_surplus(sim%ft_v(:,i1));
          call sim%interp_v%interpolate_const_disp(sim%dorder1,&
               dt*sim%ex(sim%global_v(1)+i1),sim%ft_v(:,i1), sim%ft2_v(:,i1),.TRUE.);
       end do

       sim%ft_v = sim%ft2_v
       ! v2 interpolation
       do i1 = 1, sim%local_size_v(1)!interp_x%size_basis 
          !call sim%interp_v%compute_hierarchical_surplus(sim%ft_v(:,i1));
          call sim%interp_v%interpolate_const_disp(sim%dorder2,&
               dt*sim%ey(sim%global_v(1)+i1),sim%ft_v(:,i1), sim%ft2_v(:,i1),.FALSE.);
       end do
    end if
    ! Take care of the shift in the exponential    
    do i2 = 1, sim%local_size_v(2)
       vv = sim%interp_v%hierarchy(sim%global_v(2)+i2)%coordinate(2)
       sim%dv = sim%ey(sim%global_v(1)+1:sim%global_v(1)+sim%local_size_v(1))*dt
          sim%ft_v(i2,:) = sim%ft2_v(i2,:)*exp(-vv*sim%dv-sim%dv**2*0.5_f64)
    end do

  end subroutine advection_v_df

!------------------------------------------------------------------------------!

  subroutine compute_rho(sim)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    
    sll_int32 :: i1

    sim%ft2_v=sim%ft_v;
    do i1 = 1, sim%local_size_v(1)
       call sim%interp_v%compute_linear_hierarchical_surplus(sim%ft2_v(:,i1))
       call sim%interp_v%integrate_trapezoidal(sim%ft2_v(:,i1),&
            sim%rho_loc(i1))
       sim%rho_loc(i1) = 1.0_f64-sim%rho_loc(i1);
    end do
    
    call sll_s_collective_alltoallv_double(sim%rho_loc, sim%send_cnt, sim%send_displs,&
         sim%rho, sim%recv_cnt, sim%recv_displs, sll_v_world_collective);
    
  end subroutine compute_rho
  !------------------------------------------------------------------------------!

  subroutine compute_rho_df(sim)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    sll_real64 :: h, x0
    sll_real64, dimension(2) :: phi
    sll_int32 :: i1, i2

    sim%ft2_v=sim%ft_v;
    do i1 = 1, sim%local_size_v(1)
       call sim%interp_v%compute_linear_hierarchical_surplus(sim%ft2_v(:,i1))
       sim%rho_loc(i1) = 0.0_f64
       do i2 = 1, sim%local_size_v(2)
          if (sim%test_case == SLL_LANDAU) then
             if (sim%interp_v%hierarchy(sim%global_v(2)+i2)%level(1) == 0) then
                phi(1) = sqrt(2.0_f64*sll_p_pi)
             else
                h = sim%interp_v%length(1)/&
                     (real(2**sim%interp_v%hierarchy(sim%global_v(2)+i2)%level(1),f64));
                x0 = sim%interp_v%hierarchy(sim%global_v(2)+i2)%coordinate(1)
                phi(1) = int_f0_alpha(x0,h,0.5_f64);
             end if
          elseif (sim%test_case == SLL_TSI) then
             phi(1) = sim%interp_v%length(1)/&
                  real(max(2**(sim%interp_v%hierarchy(sim%global_v(2)+i2)%level(1)),2), f64)
          end if
          if (sim%interp_v%hierarchy(sim%global_v(2)+i2)%level(2) == 0) then
             phi(2) = sqrt(2.0_f64*sll_p_pi) 
          else
             h = sim%interp_v%length(2)&
                  /(real(2**sim%interp_v%hierarchy(sim%global_v(2)+i2)%level(2),f64));
             x0 = sim%interp_v%hierarchy(sim%global_v(2)+i2)%coordinate(2)
             phi(2) = int_f0_alpha(x0,h,0.5_f64);
          end if
          sim%rho_loc(i1) = sim%rho_loc(i1) + sim%ft2_v(i2,i1)*product(phi)
       end do
       sim%rho_loc(i1) = 1.0_f64-sim%rho_loc(i1);
    end do

    call sll_s_collective_alltoallv_double(sim%rho_loc, sim%send_cnt, sim%send_displs,&
         sim%rho, sim%recv_cnt, sim%recv_displs, sll_v_world_collective);
    
  end subroutine compute_rho_df
!------------------------------------------------------------------------------!

  subroutine compute_energy(sim, val)
    class(sll_t_sim_sl_vp_2d2v_cart_sparsegrid), intent(inout) :: sim
    sll_real64, intent(out) :: val
    
    sll_int32 :: i1

    do i1=1,sim%interp_x%size_basis
       sim%ex(i1) = sim%ex(i1)*sim%ex(i1) + sim%ey(i1)*sim%ey(i1)
    end do
    call sim%interp_x%compute_linear_hierarchical_surplus(sim%ex)
    call sim%interp_x%integrate_trapezoidal(sim%ex,&
         val);
  end subroutine compute_energy

!------------------------------------------------------------------------------!
function int_f0_alpha(x0,h,alpha) result(intf0)
  sll_real64, intent(in) :: x0, h, alpha
  sll_real64 :: intf0

  intf0 = sqrt(sll_p_pi/alpha)*0.5_f64*(&
       (1.0_f64-x0/h)*(erf(x0*sqrt(alpha))-erf((x0-h)*sqrt(alpha)))+&
       (1.0_f64+x0/h)*(erf((x0+h)*sqrt(alpha))-erf(x0*sqrt(alpha))))+&
       (exp(-(x0-h)**2*alpha)+exp(-(x0+h)**2*alpha)-2.0_f64*exp(-x0**2*alpha))/&
       (2.0_f64*alpha*h);

end function int_f0_alpha

!------------------------------------------------------------------------------!

end module sll_m_sim_bsl_vp_2d2v_cart_sparsegrid
