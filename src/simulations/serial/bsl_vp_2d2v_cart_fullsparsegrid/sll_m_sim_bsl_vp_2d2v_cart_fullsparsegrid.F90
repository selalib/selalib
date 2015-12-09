!------------------------------------------------------------------------------!
! Implementation of sgxv 4D sparse grid semi-Lagrangian solver for
! 4D Vlasov-Poisson
! Implemented test cases: Landau damping (test_case = SLL_LANDAU) and
! TSI (test_case = SLL_TSI) where TSI along x1v1 and Landau along x2v2.
!
!------------------------------------------------------------------------------

module sll_m_sim_bsl_vp_2d2d_cart_fullsparsegrid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_pi

  use sll_m_poisson_2d_sparse_grid_fft, only: &
    sll_fft_derivative

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_m_sparse_grid_2d, only: &
    sparse_grid_interpolator_2d

  use sll_m_sparse_grid_4d, only: &
    sparse_grid_interpolator_4d

  implicit none

  public :: &
    sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_LANDAU = 0
  sll_int32, parameter :: SLL_TSI = 1

  
  type, extends(sll_simulation_base_class) :: sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid
     logical     :: is_mdeltaf
     sll_int32   :: test_case
     sll_int32   :: levelxv
     sll_int32   :: order_sg

     ! Domain
     sll_real64  :: eta_min(4), eta_max(4)

     !Time domain
     sll_int32  :: n_time_steps
     sll_real64 :: delta_t

     !Sparse grid parameters
     sll_int32 :: order
     sll_int32, dimension(4) :: levels, dorder_x1, dorder_x2, dorder_v1, dorder_v2
     sll_int32 :: size_basis

     !Distribution function 6D
     sll_real64 :: eps, v2, v0
     sll_real64, dimension(:), allocatable :: f_xv, f2_xv
     
     !Electric fields and charge density
     sll_real64, dimension(:), allocatable :: ex
     sll_real64, dimension(:), allocatable :: ey
     sll_real64, dimension(:), allocatable :: ez
     sll_real64, dimension(:), allocatable :: rho

     !Poisson solver
     type(sll_fft_derivative)  :: poisson
     
     ! Interpolator
     type(sparse_grid_interpolator_2d)   :: interp_x
     type(sparse_grid_interpolator_4d)   :: interp_xv

     !Diagnostics
     sll_real64, dimension(:), allocatable :: nrj

     contains
       procedure :: init_from_file => init_sparsegrid_2d2v
       procedure :: run => run_sparsegrid_2d2v

    end type sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid

contains
!------------------------------------------------------------------------------!
  subroutine init_sparsegrid_2d2v (sim, filename)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    character(len=*), intent(in)                               :: filename

    ! Local variables
    sll_int32 :: io_stat
    logical   :: is_mdeltaf
    character(len=256) :: test_case
    sll_int32   :: levelxv
    sll_int32   :: order_sg
    sll_real64  :: delta_t
    sll_int32   :: n_time_steps

    sll_int32, parameter :: input_file = 99

    namelist /sim_params/ is_mdeltaf, test_case, levelxv, order_sg, delta_t, n_time_steps

    ! Read parameters from file
    open( unit = input_file, file=trim(filename), IOStat= io_stat)
    if (io_stat /= 0) then
       print*, 'init_sparsegrid_2d2v() failed to open file', filename
       STOP
    end if
    read(input_file, sim_params)
    close(input_file)

    sim%is_mdeltaf = is_mdeltaf
    sim%levelxv = levelxv
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
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim

    sll_int32  :: j, i1
    sll_int32  :: error
    sll_real64 :: eta(4)
    sll_real64 :: kxy(3)
    sll_real64 :: time
    sll_int32  :: i_step

    !Sparse grid
    sim%levels = sim%levelxv

    ! Information on ordering for 1D interpolations
    sim%dorder_x1 = [1, 3, 2, 4]
    sim%dorder_x2 = [2, 4, 1, 3]
    sim%dorder_v1 = [3, 1, 2, 4]
    sim%dorder_v2 = [4, 1, 2, 3]

    ! Set the domain for one of the two test cases
    if (sim%test_case == SLL_LANDAU ) then
       !x domain
       sim%eta_min(1) =  0.0_f64; sim%eta_max(1) =  4.0_f64 * sll_pi
       sim%eta_min(2) = sim%eta_min(1); sim%eta_max(2) = sim%eta_max(1)
       !v domain
       sim%eta_min(3) = -6.0_f64; sim%eta_max(3) = 6.0_f64
       sim%eta_min(4) = sim%eta_min(3); sim%eta_max(4) = sim%eta_max(3)
    elseif (sim%test_case  == SLL_TSI) then
       !x domain
       sim%eta_min(1) =  0.0_f64; sim%eta_max(1) =  10.0_f64 * sll_pi
       sim%eta_min(2) =  0.0_f64; sim%eta_max(2) =  4.0_f64 * sll_pi
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
    call sim%interp_x%initialize(sim%levels(1:2),sim%order_sg, sim%order_sg+1,0, &
         sim%eta_min(1:2),sim%eta_max(1:2), 0, 0);
    call sim%interp_xv%initialize(sim%levels, sim%order_sg, sim%order_sg+1, 0, &
         sim%eta_min, sim%eta_max)
    ! Initialize Poisson
    call sim%poisson%initialize(sim%interp_x);
    
    sim%size_basis = sim%interp_xv%size_basis

    ! Allocate the arrays
    SLL_ALLOCATE(sim%f_xv(sim%size_basis), error)
    SLL_ALLOCATE(sim%f2_xv(sim%size_basis), error)
    SLL_ALLOCATE(sim%rho(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%ex(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%ey(sim%interp_x%size_basis),error)
    SLL_ALLOCATE(sim%ez(sim%interp_x%size_basis),error)
    
    SLL_ALLOCATE(sim%nrj(0:sim%n_time_steps), error)


  
    ! Set initial value
    do j=1,2
       kxy(j)  = 2.0_f64*sll_pi/(sim%eta_max(j)-sim%eta_min(j))
    end do

    do i1=1,sim%size_basis
       eta = sim%interp_xv%hierarchy(i1)%coordinate
       sim%f_xv(i1) = 1.0_f64/(2.0_f64*sll_pi)*&
            (1.0_f64+sim%eps*(cos(kxy(1)*eta(1))+&
            cos(kxy(2)*eta(2))))
       if (sim%test_case == SLL_TSI) then
          if (sim%is_mdeltaf .EQV. .FALSE.) then
             sim%f_xv(i1) = sim%f_xv(i1)*0.5_f64*&
                  (exp(-(eta(3)-sim%v0)**2*0.5_f64) + &
                  exp(-(eta(3)+sim%v0)**2*0.5_f64))* &
                  exp(-(eta(4)**2)*0.5_f64)
          else
             sim%f_xv(i1) = sim%f_xv(i1)*0.5_f64*&
                  (exp(-(eta(3)-sim%v0)**2*0.5_f64) + &
                  exp(-(eta(3)+sim%v0)**2*0.5_f64))
          end if
       elseif ((sim%test_case == SLL_LANDAU) .AND. (sim%is_mdeltaf .EQV. .FALSE.)) then
          sim%f_xv(i1) = sim%f_xv(i1)*exp(-.5_f64*( eta(3)**2+eta(4)**2))
       end if
    end do

    time = 0.0_f64

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
    
    call compute_energy(sim, sim%nrj(0));
       
    open(11, file='sgxv4d_output.dat', position='append')
    rewind(11)
    write(11,*) time, sim%nrj(0)
    
    do i_step = 1, sim%n_time_steps !Loop over time

       if(mod(i_step,100) == 0) then
          print*, i_step
       end if
       
       call advection_x(sim, sim%delta_t)

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
       call compute_energy(sim, sim%nrj(i_step));
       
       time  = time + sim%delta_t
       
       open(11, file='sgxsgv4d_output.dat', position='append')
       write(11,*) time, sim%nrj(i_step)
       close(11)
       
    end do !next time step


  end subroutine run_sparsegrid_2d2v
  
!------------------------------------------------------------------------------!
  subroutine advection_x(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt

    ! Advection along x1
    call sim%interp_xv%compute_hierarchical_surplus(sim%f_xv)
    call sim%interp_xv%interpolate_disp_linnconst_in_1d( -dt,&
         sim%dorder_x1, sim%f_xv, sim%f2_xv)
   
    ! Advection along x2
    call sim%interp_xv%compute_hierarchical_surplus(sim%f2_xv)
    call sim%interp_xv%interpolate_disp_linnconst_in_1d( -dt,&
         sim%dorder_x2, sim%f2_xv, sim%f_xv)

  end subroutine advection_x
!------------------------------------------------------------------------------!


  subroutine advection_v(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt

    ! Advection along v1
    call sim%interp_xv%compute_hierarchical_surplus(sim%f_xv)
    call sim%interp_xv%interpolate_disp_nconst_in_2d(dt*sim%ex,&
         sim%dorder_v1, sim%f_xv, sim%f2_xv)

    ! Advection along v2
    call sim%interp_xv%compute_hierarchical_surplus(sim%f2_xv)
    call sim%interp_xv%interpolate_disp_nconst_in_2d(dt*sim%ey,&
         sim%dorder_v2, sim%f2_xv, sim%f_xv)
    
  end subroutine advection_v

  subroutine advection_v_df(sim, dt)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    sll_real64, intent(in) :: dt

    ! Advection along v1
    call sim%interp_xv%compute_hierarchical_surplus(sim%f_xv)
    call sim%interp_xv%interpolate_disp_nconst_in_2d(dt*sim%ex,&
         sim%dorder_v1, sim%f_xv, sim%f2_xv)
    ! Delta f along x1 only for Landau test case
    if (sim%test_case == SLL_LANDAU) then
       call scale_gaussian(sim, 3,dt,sim%ex,sim%f2_xv)
    end if

    ! Advection along v2
    call sim%interp_xv%compute_hierarchical_surplus(sim%f2_xv)
    call sim%interp_xv%interpolate_disp_nconst_in_2d(dt*sim%ey,&
         sim%dorder_v2, sim%f2_xv, sim%f_xv)
    call scale_gaussian(sim, 4,dt,sim%ey, sim%f_xv)

  end subroutine advection_v_df

!------------------------------------------------------------------------------!

  subroutine compute_rho(sim)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    
    sll_int32  :: counter, counter2, lev1, lev2, l1, l2, l3, l4, k1, k2, k, no 
    sll_real64 :: factor

    sim%f2_xv = sim%f_xv
    call sim%interp_xv%compute_linear_hierarchical_surplus(sim%f2_xv)
    counter = 1
    do lev1 = 0, max(sim%interp_xv%levels(1), sim%interp_xv%levels(2))
       do l1 = max(0, lev1-sim%interp_xv%levels(2)), min (lev1, sim%interp_xv%levels(1))
          l2 = lev1-l1
          counter2 = 0
          do k1 = 0, max(2**(l1-1),1)-1
             do k2 = 0, max(2**(l2-1),1)-1
                sim%rho(counter) = 0.0_f64
                do lev2 = 0, sim%interp_xv%max_level-lev1
                   factor = 1.0_f64/real(2**lev2,f64)
                   do l3 = max(0, lev2-sim%interp_xv%levels(4)), min(lev2, sim%interp_xv%levels(3))
                      l4 = lev2-l3
                      no = 2**(max(l3-1,0))*2**(max(l4-1,0))
                      do k=sim%interp_xv%index(l1,l2,l3,l4)+ no*counter2, &
                           sim%interp_xv%index(l1,l2,l3,l4)+ no*(counter2+1)-1
                         sim%rho(counter) = sim%rho(counter) + &
                              sim%f2_xv(k)*factor
                      end do
                   end do
                end do
                sim%rho(counter) = sim%rho(counter)* sim%interp_xv%length(3)*&
                     sim%interp_xv%length(4)
                counter = counter+1
                counter2 = counter2+1
             end do
          end do
       end do
    end do

    call dehira_landau(sim%interp_x, sim%rho)

    do k=1, sim%interp_x%size_basis
       sim%rho(k) = 1.0_f64 - sim%rho(k)
    end do

  end subroutine compute_rho
  !------------------------------------------------------------------------------!

  subroutine compute_rho_df(sim)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim


    sll_int32  :: counter, counter2, lev1, lev2, l1, l2, l3, l4, k1, k2, k, no 
    sll_real64 :: factor
    sll_real64 :: h, x0

    sim%f2_xv = sim%f_xv
    call sim%interp_xv%compute_linear_hierarchical_surplus(sim%f2_xv)
    counter = 1
    do lev1 = 0, sim%interp_xv%max_level
       do l1 = max(0, lev1-sim%interp_xv%levels(2)), min (lev1, sim%interp_xv%levels(1))
          l2 = lev1-l1
          counter2 = 0
          do k1 = 0, max(2**(l1-1),1)-1
             do k2 = 0, max(2**(l2-1),1)-1
                sim%rho(counter) = 0.0_f64
                do lev2 =  0, sim%interp_xv%max_level-lev1
                   do l3 = max(0, lev2-sim%interp_xv%levels(4)), min(lev2, sim%interp_xv%levels(3))
                      l4 = lev2-l3
                      no = 2**(max(l3-1,0))*2**(max(l4-1,0))
                      do k=sim%interp_xv%index(l1,l2,l3,l4)+ no*counter2, &
                           sim%interp_xv%index(l1,l2,l3,l4)+ no*(counter2+1)-1
                         if (sim%test_case == SLL_LANDAU) then
                            if(sim%interp_xv%hierarchy(k)%level(3) == 0) then
                               factor = sqrt(2.0_f64*sll_pi)
                            else
                               h = sim%interp_xv%length(3)/ &
                                    (real(2**sim%interp_xv%hierarchy(k)%level(3), f64))
                               x0 = sim%interp_xv%hierarchy(k)%coordinate(3)
                               factor = int_f0_alpha(x0,h,0.5_f64)
                            end if
                         elseif (sim%test_case == SLL_TSI) then
                            factor = sim%interp_xv%length(3)/&
                                 real(max(2**(sim%interp_xv%hierarchy(k)%level(3)),2), f64)
                         end if
                         if(sim%interp_xv%hierarchy(k)%level(4) == 0) then
                            factor = factor*sqrt(2.0_f64*sll_pi)
                         else
                            h = sim%interp_xv%length(4)/ &
                                 (real(2**sim%interp_xv%hierarchy(k)%level(4), f64))
                            x0 = sim%interp_xv%hierarchy(k)%coordinate(4)
                            factor = factor*int_f0_alpha(x0,h,0.5_f64)
                         end if
                         sim%rho(counter) = sim%rho(counter) + &
                              sim%f2_xv(k)*factor
                      end do
                   end do
                end do
                counter = counter+1
                counter2 = counter2+1
             end do
          end do
       end do
    end do

    call dehira_landau(sim%interp_x, sim%rho)

    do k=1, sim%interp_x%size_basis
       sim%rho(k) = 1.0_f64 - sim%rho(k)
    end do

  end subroutine compute_rho_df
!------------------------------------------------------------------------------!

  subroutine compute_energy(sim, val)
    class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
    sll_real64, intent(out) :: val
    

    sll_int32 :: i1

    do i1=1,sim%interp_x%size_basis
       sim%ex(i1) = sim%ex(i1)*sim%ex(i1) + sim%ey(i1)*sim%ey(i1)
    end do
    call sim%interp_x%compute_linear_hierarchical_surplus(sim%ex)
    call sim%interp_x%integrate_trapezoidal(sim%ex, val)

  end subroutine compute_energy

!------------------------------------------------------------------------------!
function int_f0_alpha(x0,h,alpha) result(intf0)
  sll_real64, intent(in) :: x0, h, alpha
  sll_real64 :: intf0

  intf0 = sqrt(sll_pi/alpha)*0.5_f64*(&
       (1.0_f64-x0/h)*(erf(x0*sqrt(alpha))-erf((x0-h)*sqrt(alpha)))+&
       (1.0_f64+x0/h)*(erf((x0+h)*sqrt(alpha))-erf(x0*sqrt(alpha))))+&
       (exp(-(x0-h)**2*alpha)+exp(-(x0+h)**2*alpha)-2.0_f64*exp(-x0**2*alpha))/&
       (2.0_f64*alpha*h);

end function int_f0_alpha

!------------------------------------------------------------------------------!

subroutine scale_gaussian(sim, dim,dt,efield,fscale)
  class(sll_t_sim_sl_vp_2d2v_cart_fullsparsegrid), intent(inout) :: sim
  sll_int32, intent(in) :: dim
  sll_real64, intent(in) :: dt
  sll_real64, dimension(:), intent(in) :: efield
  sll_real64, dimension(:), intent(inout) :: fscale

  sll_int32 :: l1,l2,l3,k1,k2,k, counter, counter4d,lev1, lev2
  sll_int32, dimension(4) :: lvec,novec,kvec  
  sll_int32 :: dimc
  sll_real64 :: vv, dv

     if (dim == 3) then
        dimc = 4;
     else
        dimc = 3;
     end if

     counter = 0;
     do lev1 = 0, sim%interp_xv%max_level
        do l1 = max(0, lev1-sim%interp_xv%levels(2)), min(lev1, sim%interp_xv%levels(1))
           novec(1) = max(2**(l1-1),1)
           lvec(1) = l1;
           l2 = lev1-l1
           novec(2) = max(2**(l2-1),1)
           lvec(2) = l2;
           do k1=0,novec(1)-1
              kvec(1) = k1
              do k2=0,novec(2)-1
                 kvec(2) = k2
                 counter = counter + 1;
                 dv = dt*efield(counter);
                 do lev2 = 0, sim%interp_xv%max_level-lev1
                    do l3 = max(0,lev2-sim%interp_xv%levels(dimc)), &
                         min(lev2,sim%interp_xv%levels(dim))
                       lvec(dim) = l3;
                       lvec(dimc) = lev2-l3;
                       novec(dim) = max(2**(l3-1),1);
                       novec(dimc) = max(2**(lvec(dimc)-1),1);
                       counter4d = sim%interp_xv%index(l1,l2,lvec(3),lvec(4))+&
                            kvec(1)*novec(2)*novec(3)*novec(4)+&
                            kvec(2)*novec(3)*novec(4)
                       do k = counter4d,counter4d+novec(3)*novec(4)-1
                          vv = sim%interp_xv%hierarchy(k)%coordinate(dim);
                          fscale(k) = fscale(k) *&
                               exp(-vv*dv-dv**2*0.5_f64);
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
   end subroutine scale_gaussian

!------------------------------------------------------------------------------!
 subroutine dehira_landau(interpolator,data)
    class(sparse_grid_interpolator_2d) :: interpolator
    sll_real64, dimension(:), intent(inout)   :: data
    sll_int32                                :: counter,level,l1,k1,k2
    sll_int32, dimension(2) :: l,no
    sll_real64 :: factor

    counter = 1
    do level = 0, interpolator%max_level
       do l1 = max(0,level-interpolator%levels(2)) , min(level,interpolator%levels(1))
          l(1) = l1
          no(1) = max(2**(l(1)-1),1)
          l(2) = level-l(1)
          no(2) = max(2**(l(2)-1),1)
          do k1=0,no(1)-1
             do k2=0,no(2)-1
                factor = -1.0_f64
                call dehira_landau_d_dimension&
                     (interpolator,data(counter),&
                     data,l,factor,counter,0)
                counter = counter + 1
             end do
          end do
       end do
    end do

  end subroutine dehira_landau


!------------------------------------------------------------------------------!
  recursive subroutine dehira_landau_d_dimension(interpolator,surplus,data_array,level,factor,index,d)
  class(sparse_grid_interpolator_2d) :: interpolator
  sll_real64, intent(inout) :: surplus
  sll_real64, dimension(:), intent(in) :: data_array
  sll_int32, dimension(:), intent(in) :: level
  sll_real64, intent(in) :: factor
  sll_int32, intent(in) :: index
  sll_int32, intent(in) :: d
  sll_int32 :: j

 ! if (index .NE. -1) then
     if(d>0) then
        surplus = surplus + data_array(index)*factor
     end if
     !write(16,*) d,factor,data_array(index),surplus
     do j=d+1,2
        if (level(j)>0) then
           call dehira_landau_d_dimension(interpolator,surplus,data_array,level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*j-1),j)
           call dehira_landau_d_dimension(interpolator,surplus,data_array,level,&
                -0.5_f64*factor,interpolator%hierarchy(index)%parent(2*j),j)
        end if
     end do
 ! end if

   end subroutine dehira_landau_d_dimension


end module sll_m_sim_bsl_vp_2d2d_cart_fullsparsegrid
