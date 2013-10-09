!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!> @author We will see
!> @brief 
!> Simulation class to solve slab drift kinetic equation in polar coordinates
!> (3d space (r,\theta,z) 1d velocity (v))
!> translation of slv2d/src/vp4d_dk.F90 program in simulation class
!> intended to be close to sll_simulation_4d_DK_hybrid_module
!> but specific use of polar coordinates
!> @details
!> Example of use in test program: see unit_test_4d_dk_polar.F90 file
!> 
!> \code
!>
!>  use sll_simulation_4d_drift_kinetic_polar
!>  type(sll_simulation_4d_vp_polar)    :: simulation
!>  call simulation%init_from_file(trim(filename))
!>  call simulation%run()
!>  call delete(simulation)
!> \endcode


module sll_simulation_4d_drift_kinetic_polar_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_fdistribu4D_DK
  use sll_logical_meshes
  implicit none

  
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_4d_drift_kinetic_polar

     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     ! Mesh parameters
     type(sll_logical_mesh_1d), pointer :: r_grid
     type(sll_logical_mesh_1d), pointer :: theta_grid
     type(sll_logical_mesh_1d), pointer :: phi_grid
     type(sll_logical_mesh_1d), pointer :: vpar_grid
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     sll_int32  :: nc_x3
     sll_int32  :: nc_x4
     !sll_real64 :: r_min
     !sll_real64 :: r_max
     !sll_real64 :: phi_min
     !sll_real64 :: phi_max
     !sll_real64 :: vpar_min
     !sll_real64 :: vpar_max
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     !sll_int32  :: spline_degree_eta1, spline_degree_eta2
     !sll_int32  :: spline_degree_eta3, spline_degree_eta4
     !--> Equilibrium
     sll_real64 :: tau0      !-> tau0 = Ti(rpeak)/Te(rpeak)
     sll_real64 :: rho_peak    
     sll_real64 :: kappan   
     sll_real64 :: deltarn  
     sll_real64 :: kappaTi  
     sll_real64 :: deltarTi 
     sll_real64 :: kappaTe  
     sll_real64 :: deltarTe
     sll_real64 :: n0_at_rpeak     
     !--> Pertubation
     sll_int32  :: perturb_choice
     sll_int32  :: mmode
     sll_int32  :: nmode
     sll_real64 :: eps_perturb   

     !--> 4D logical mesh (r,theta,phi,vpar)
     type(sll_logical_mesh_4d), pointer :: logical_mesh4d


     !--> Density and temperature profiles
     sll_real64, dimension(:)  , pointer :: n0_r
     sll_real64, dimension(:)  , pointer :: Ti_r
     sll_real64, dimension(:)  , pointer :: Te_r

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_x1x4


     !--> 4D distribution function 
     !----> sequential in (x1,x2) and parallel in (x3,x4)
     type(layout_4D), pointer :: layout4d_x1x2
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x1x2 
     !----> parallel in (x1,x2) and sequential in (x3,x4) 
     type(layout_4D), pointer :: layout4d_x3x4
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x3x4

     !--> 3D charge density and 3D electric potential
     !----> sequential in (x1,x2)
     type(layout_3D), pointer :: layout3d_x1x2
     sll_real64, dimension(:,:,:), pointer :: rho3d_x1x2 
     sll_real64, dimension(:,:,:), pointer :: phi3d_x1x2 
     !----> sequential in x3
     type(layout_3D), pointer :: layout3d_x3
     sll_real64, dimension(:,:,:), pointer :: rho3d_x3
     sll_real64, dimension(:,:,:), pointer :: phi3d_x3



   contains
     procedure, pass(sim) :: run => run_dk4d_polar
     procedure, pass(sim) :: init_from_file => init_dk4d_polar
  end type sll_simulation_4d_drift_kinetic_polar

  interface delete
     module procedure delete_dk4d_polar
  end interface delete

contains

  subroutine init_dk4d_polar( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99

    !--> Mesh
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: phi_min
    sll_real64 :: phi_max
    sll_real64 :: vpar_min
    sll_real64 :: vpar_max
    !--> Equilibrium
    sll_real64 :: tau0
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe     
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    !sll_int32  :: spline_degree

    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      r_min, r_max, phi_min, phi_max, &
      vpar_min, vpar_max
    namelist /equilibrium/ tau0, rho_peak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb
    namelist /sim_params/ dt, number_iterations!, spline_degree

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_dk4d_polar() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
    close(input_file)

    !--> Mesh
    sim%nc_x1    = num_cells_x1
    sim%nc_x2    = num_cells_x2
    sim%nc_x3    = num_cells_x3
    sim%nc_x4    = num_cells_x4
    !sim%r_min    = r_min
    !sim%r_max    = r_max
    !sim%phi_min  = phi_min  
    !sim%phi_max  = phi_max
    !sim%vpar_min = vpar_min
    !sim%vpar_min = vpar_max
    
    sim%r_grid => new_logical_mesh_1d(num_cells_x1,eta_min=r_min,eta_max=r_max)
    sim%theta_grid => new_logical_mesh_1d(num_cells_x2,&
      eta_min=0._f64,eta_max=2._f64*sll_pi)
    sim%phi_grid => new_logical_mesh_1d(num_cells_x3,eta_min=phi_min,eta_max=phi_max)
    sim%vpar_grid => new_logical_mesh_1d(num_cells_x4,eta_min=vpar_min,eta_max=vpar_max)
    
    !--> Equilibrium
    sim%tau0     = tau0
    sim%rho_peak = rho_peak 
    sim%kappan   = kappan
    sim%deltarn  = deltarn
    sim%kappaTi  = kappaTi
    sim%deltarTi = deltarTi
    sim%kappaTe  = kappaTe
    sim%deltarTe = deltarTe
    !--> Pertubation
    sim%perturb_choice = perturb_choice
    sim%mmode          = mmode
    sim%nmode          = nmode
    sim%eps_perturb    = eps_perturb
    !--> Algorithm
    sim%dt                 = dt
    sim%num_iterations     = number_iterations
    !sim%spline_degree_eta1 = spline_degree
    !sim%spline_degree_eta2 = spline_degree
    !sim%spline_degree_eta3 = spline_degree
    !sim%spline_degree_eta4 = spline_degree


    if(sll_get_collective_rank(sll_world_collective)==0)then
      print *,'##Mesh'
      print *,'#num_cells_x1=',num_cells_x1
      print *,'#num_cells_x2=',num_cells_x2
      print *,'#num_cells_x3=',num_cells_x3
      print *,'#num_cells_x4=',num_cells_x4
      print *,'#r_min=',r_min
      print *,'#r_max=',r_max
      print *,'#phi_min=',phi_min
      print *,'#phi_max=',phi_max
      print *,'#vpar_min=',vpar_min
      print *,'#vpar_max=',r_max
      print *,'##equilibrium'
      print *,'#tau0=',tau0
      print *,'#rho_peak=',rho_peak
      print *,'#kappan=',kappan
      print *,'#deltarn=',deltarn
      print *,'#kappaTi=',kappaTi
      print *,'#deltarTi=',deltarTi
      print *,'#kappaTe=',kappaTe
      print *,'#deltarTe=',deltarTe
      print *,'##perturbation'
      print *,'#perturb_choice=',perturb_choice
      print *,'#mmode=',mmode
      print *,'#nmode=',nmode
      print *,'#eps_perturb=',eps_perturb
      print *,'#dt=',dt
      print *,'#number_iterations=',number_iterations
    endif
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    
    call initialize_profiles_analytic(sim)    
    call allocate_fdistribu4d_DK(sim)
    call allocate_QN_DK( sim )


  end subroutine init_dk4d_polar

  subroutine run_dk4d_polar(sim)
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    !--> For initial profile HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"

    !*** Saving of the radial profiles in HDF5 file ***
    if (sll_get_collective_rank(sll_world_collective)==0) then
      call sll_hdf5_file_create(filename_prof,file_id,file_err)
      call sll_hdf5_write_array_1d(file_id,sim%n0_r,'n0_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti_r,'Ti_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te_r,'Te_r',file_err)
      call sll_hdf5_file_close(file_id,file_err)
    end if


    !*** Initialization of the distribution function ***
    !***  i.e f4d(t=t0)                              ***
    call initialize_fdistribu4d_DK(sim)
    
    call compute_charge_density(sim)


  end subroutine run_dk4d_polar


  subroutine delete_dk4d_polar( sim )
    class(sll_simulation_4d_drift_kinetic_polar) :: sim
    sll_int32 :: ierr
  end subroutine delete_dk4d_polar
  
  
  subroutine initialize_profiles_analytic(sim)
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    sll_int32 :: i,ierr,Nr
    sll_real64 :: r,dr,rp,tmp,r_min,r_max
    
    Nr = sim%r_grid%num_cells
    dr = sim%r_grid%delta_eta
    r_min = sim%r_grid%eta_min
    r_max = sim%r_grid%eta_max
    
    SLL_ALLOCATE(sim%n0_r(Nr+1),ierr)
    SLL_ALLOCATE(sim%Ti_r(Nr+1),ierr)
    SLL_ALLOCATE(sim%Te_r(Nr+1),ierr)
    
    rp = r_min+sim%rho_peak*(r_max-r_min)
    do i=1,Nr+1
      r = r_min+real(i-1,f64)*dr
      sim%n0_r(i) = exp(-sim%kappan*sim%deltarn*tanh((r-rp)/sim%deltarn))
      sim%Ti_r(i)=exp(-sim%kappaTi*sim%deltarTi*tanh((r-rp)/sim%deltarTi))    
      sim%Te_r(i)=exp(-sim%kappaTe*sim%deltarTe*tanh((r-rp)/sim%deltarTe))    
    enddo
    
    !we then change the normalization for n0_r
    tmp = 0.5_f64*(sim%n0_r(1)*r_min+sim%n0_r(Nr+1)*r_max)
    do i = 2,Nr
      r = r_min+real(i-1,f64)*dr
      tmp = tmp + sim%n0_r(i)*r
    enddo
    tmp = tmp/real(Nr,f64)
    sim%n0_r = sim%n0_r/tmp
    sim%n0_at_rpeak = 1._f64/tmp      
  
  end subroutine initialize_profiles_analytic
  



  !----------------------------------------------------
  ! Allocation of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine allocate_fdistribu4d_DK( sim )
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim

    sll_int32 :: ierr, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4

    ! layout for sequential operations in x3 and x4. 
    ! Make an even split for x1 and x2, or as close as 
    ! even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up 
    ! at initialization time.
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**(sim%power2/2)
       sim%nproc_x4 = 2**(sim%power2/2)
    else 
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**((sim%power2-1)/2)
       sim%nproc_x4 = 2**((sim%power2+1)/2)
    end if

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x1,x2) : sequential
    !-->  (x3,x4) : parallelized layout
    sim%layout4d_x1x2  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x1x2 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_x1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(sim%f4d_x1x2(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x1,x2) directions
    !-->  (x1,x2) : parallelized layout
    !-->  (x3,x4) : sequential
    ! switch x3 and x1:
    itemp        = sim%nproc_x1
    sim%nproc_x1 = sim%nproc_x3
    sim%nproc_x3 = itemp
    ! switch x4 and x2
    itemp        = sim%nproc_x2
    sim%nproc_x2 = sim%nproc_x4 
    sim%nproc_x4 = itemp

    sim%layout4d_x3x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nc_x4+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x3x4 )
        
    call compute_local_sizes_4d( sim%layout4d_x3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_x3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim)
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim

    sll_int32  :: ierr
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32, dimension(1:4) :: glob_ind
    sll_real64, dimension(:), pointer :: x1_node,x2_node,x3_node,x4_node
    sll_real64 :: rp,k_x2,k_x3
    sll_real64 :: tmp_mode,tmp
    sll_real64 :: r_min,r_max      

    !--> Initialization of the equilibrium distribution function
    SLL_ALLOCATE(sim%feq_x1x4(sim%nc_x1+1,sim%nc_x4+1),ierr)
    
    
    r_min = sim%r_grid%eta_min
    r_max = sim%r_grid%eta_max
        
    call initialize_x1_node_1d(sim%r_grid,x1_node)
    call initialize_x1_node_1d(sim%theta_grid,x2_node)
    call initialize_x1_node_1d(sim%phi_grid,x3_node)
    call initialize_x1_node_1d(sim%vpar_grid,x4_node)
    
    call init_fequilibrium( sim%nc_x1+1, &
      sim%nc_x4+1, &
      x1_node, &
      x4_node, &
      sim%n0_r, &
      sim%Ti_r, &
      sim%feq_x1x4 )

    !--> Initialization of the distribution function f4d_x3x4
    call compute_local_sizes_4d( sim%layout4d_x3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
   
   k_x2  = 2._f64*sll_pi/(sim%theta_grid%eta_max -sim%theta_grid%eta_min)
   k_x3  = 2._f64*sll_pi/(sim%phi_grid%eta_max-sim%phi_grid%eta_min)
      
   rp = r_min+sim%rho_peak*(r_max-r_min) 
    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind(:) = local_to_global_4D(sim%layout4d_x3x4, &
              (/iloc1,iloc2,iloc3,iloc4/))
            i1 = glob_ind(1)
            i2 = glob_ind(2)
            i3 = glob_ind(3)
            i4 = glob_ind(4)
            tmp_mode = cos(real(sim%nmode,f64)*k_x3*x3_node(i3)&
               +real(sim%mmode,f64)*k_x2*x2_node(i2))
            tmp = exp(-(x1_node(i1)-rp)**2/(4._f64*sim%deltarn/sim%deltarTi))   
            sim%f4d_x3x4(iloc1,iloc2,i3,i4) = &
              (1._f64+tmp_mode*sim%eps_perturb*tmp)*sim%feq_x1x4(i1,i4)            
          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(x1_node,ierr)
    SLL_DEALLOCATE(x2_node,ierr)
    SLL_DEALLOCATE(x3_node,ierr)
    SLL_DEALLOCATE(x4_node,ierr)
  end subroutine initialize_fdistribu4d_DK


  
  function compute_equil_analytic(sim,r,v)
    type(sll_simulation_4d_drift_kinetic_polar), intent(in) :: sim
    sll_real64,intent(in)::r,v
    sll_real64::compute_equil_analytic
    sll_real64:: tmp(2),rp,r_min,r_max
    r_min = sim%r_grid%eta_min
    r_max = sim%r_grid%eta_max
    
    rp = r_min+sim%rho_peak*(r_max-r_min)
    tmp(1) = sim%n0_at_rpeak*exp(-sim%kappan*sim%deltarn*tanh((r-rp)/sim%deltarn))
    tmp(2) = exp(-sim%kappaTi*sim%deltarTi*tanh((r-rp)/sim%deltarTi))  
    compute_equil_analytic = tmp(1)/sqrt(2._f64*sll_pi*tmp(2))*exp(-0.5_f64*v**2/tmp(2))
  
  end function compute_equil_analytic



  !----------------------------------------------------
  ! Allocation for QN solver
  !----------------------------------------------------
  subroutine allocate_QN_DK( sim )
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim

    type(sll_logical_mesh_2d), pointer :: logical_mesh2d
    sll_int32 :: ierr, itemp
    !sll_int32 :: i1, i2, i3, i4
    !sll_int32 :: iloc1, iloc2, iloc3, iloc4
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: nproc3d_x3


    ! layout for sequential operations in x3 
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**(sim%power2/2)
       sim%nproc_x4 = 2**(sim%power2/2)
    else 
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**((sim%power2-1)/2)
       sim%nproc_x4 = 2**((sim%power2+1)/2)
    end if

    !--> Initialization of rho3d_x1x2 and phi3d_x1x2
    !-->  (x1,x2) : sequential
    !-->  x3 : parallelized layout    
    sim%layout3d_x1x2  => new_layout_3D( sll_world_collective )
    nproc3d_x3 = sim%nproc_x3*sim%nproc_x4
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, & 
      sim%nc_x2+1, & 
      sim%nc_x3+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      nproc3d_x3, &
      sim%layout3d_x1x2 )
    call compute_local_sizes_3d( sim%layout3d_x1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)

    !--> Initialization of rho3d_x3 and phi3d_x3
    !-->  (x1,x2) : parallelized layout
    !-->  x3 : sequential
    ! switch x3 and x1:
    itemp        = sim%nproc_x1
    sim%nproc_x1 = sim%nproc_x3
    sim%nproc_x3 = itemp
    ! switch x4 and x2
    itemp        = sim%nproc_x2
    sim%nproc_x2 = sim%nproc_x4 
    sim%nproc_x4 = itemp
        
    sim%layout3d_x3  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%nc_x1+1, &
      sim%nc_x2+1, &
      sim%nc_x3+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%layout3d_x3 )
    call compute_local_sizes_3d( sim%layout3d_x3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    
!
!    !----->
!    sim%rho2d => new_scalar_field_2d_discrete_alt( &
!      "rho2d_x1x2", &
!      sim%interp_rho2d, &     
!      sim%transf_xy, &
!      sim%bc_left_eta1, &
!      sim%bc_right_eta1, &
!      sim%bc_left_eta2, &
!      sim%bc_right_eta2)
!
!    sim%phi2d => new_scalar_field_2d_discrete_alt( &
!      "phi2d_x1x2", &
!      sim%interp_phi2d, &     
!      sim%transf_xy, &
!      sim%bc_left_eta1, &
!      sim%bc_right_eta1, &
!      sim%bc_left_eta2, &
!      sim%bc_right_eta2)
  end subroutine allocate_QN_DK


  !-----------------------------------------------------------
  ! Computation of the charge density, i.e
  !  rho(eta1,eta2,eta3) = \int f(eta1,eta2,eta3,vpar) dvpar
  !  In : f4d_x3x4(x1 part,x2 part,x3=*,x4=*)
  !  Out: rho3d_x3(x1 part,x2 part,x3=*)
  !-----------------------------------------------------------
  subroutine compute_charge_density(sim)
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    
    call compute_reduction_x3x4_to_x3(sim%vpar_grid,sim%f4d_x3x4,sim%rho3d_x3)
  
  end subroutine compute_charge_density
  
  
  !--------------------------------------------------
  ! Generic function for computing charge density
  ! should be elsewhere
  ! we should also add a choice for the integration  
  !---------------------------------------------------  
  
  subroutine compute_reduction_x3x4_to_x3(logical_mesh_x4,data_4d_x3x4,data_3d_x3)
    type(sll_logical_mesh_1d)     , intent(in)    :: logical_mesh_x4
    sll_real64, dimension(:,:,:,:), intent(in)    :: data_4d_x3x4
    sll_real64, dimension(:,:,:)  , intent(out) :: data_3d_x3
    sll_int32  :: np_x1_loc
    sll_int32  :: np_x2_loc
    sll_int32  :: np_x3
    sll_int32  :: np_x4
    sll_int32  :: iloc1, iloc2, i3, i4
    sll_real64 :: delta_x4, tmp
    
    np_x1_loc  = size(data_4d_x3x4,1)
    np_x2_loc  = size(data_4d_x3x4,2)
    np_x3      = size(data_4d_x3x4,3)
    np_x4      = size(data_4d_x3x4,4)
    delta_x4 = logical_mesh_x4%delta_eta 

    !-> Computation of the charge density locally in (x1,x2) directions
    do i3 = 1,np_x3
      do iloc2 = 1,np_x2_loc
        do iloc1 = 1,np_x1_loc
          tmp = 0._f64
          do i4 = 1,np_x4
            tmp = tmp + data_4d_x3x4(iloc1,iloc2,i3,i4)*delta_x4
          end do
          data_3d_x3(iloc1,iloc2,i3) = tmp
        end do
      end do
    end do
  end subroutine compute_reduction_x3x4_to_x3
  
  
  subroutine solve_quasi_neutral(sim)
    type(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
  
  end subroutine solve_quasi_neutral

end module sll_simulation_4d_drift_kinetic_polar_module



