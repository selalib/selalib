module sll_simulation_4d_vp_eulerian_cartesian_finite_volume_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
  use sll_collective
  use sll_remapper
  use sll_poisson_2d_periodic_cartesian_par
  use sll_simulation_base
  use sll_parallel_array_initializer_module  
  use sll_logical_meshes
  use sll_mesh_calculus_2d_module
  use sll_gnuplot_parallel
  use sll_timer
  use sll_point_to_point_comms_module
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_4d_vp_eulerian_cartesian_finite_volume
  ! Parallel environment parameters
  sll_int32  :: buf_size
  sll_int32  :: world_size
  sll_int32  :: my_rank
  sll_int32  :: power2 ! 2^power2 = number of processes available
  ! Processor mesh sizes
  sll_int32  :: nproc_v1  
  sll_int32  :: nproc_v2
  sll_int32  :: nproc_x1
  sll_int32  :: nproc_x2
  ! Physics/numerical parameters
  sll_real64 :: dt
  sll_real64 :: cfl !cfl value 
  sll_real64 :: tmax
  sll_real64 :: Enorm
  sll_int32  :: num_iterations

  ! Mesh parameters
  sll_int32  :: nc_v1  ! velocity cells
  sll_int32  :: nc_v2
  sll_int32  :: nel_boundary
  sll_int32  :: nc_x1
  sll_int32  :: nc_x2

  sll_int32 :: degree ! polynomial degree
  sll_int32 :: test
  sll_int32 :: nsch

  sll_int32  :: np_v1  ! velocity nodes
  sll_int32  :: np_v2
  sll_int32 :: np_x1
  ! used to refer to port buffers
  sll_real64, dimension(:), pointer :: buf1,buf2,buf3,buf4
  !buf1 pour prendre les valeur de port 1 et buf2 pour prendre les valeur de port2 

  ! number of local nodes in each element
  sll_int32 :: np_loc
  sll_real64 :: eps=0.001_f64 !0 center flux, if not decentered flux

  ! for initializers
  type(sll_p2p_comm_real64), pointer :: comm
  type(sll_logical_mesh_2d), pointer    :: mesh2dx,mesh2dv
  class(sll_coordinate_transformation_2d_base),pointer     :: tx,tv
  type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan

  procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
  sll_real64, dimension(:), pointer :: params

  ! 1D interpolation points repartition
  sll_real64, dimension(:),pointer  :: interp_pts_1D 
  ! finite element approximation in the velocity space
  ! connectivity array
  sll_int32,dimension(:,:),pointer  :: connec,connecline
  ! velocity interpolation points
  sll_real64, dimension(:,:),pointer  :: vcoords

  ! skyline structure of the matrices
  ! profil and pointers to the columns
  sll_int32,dimension(:),pointer  :: prof
  sll_int32,dimension(:),pointer  :: mkld
  sll_int32 :: nsky
  ! volumes of the cells and surfaces of the faces
  sll_real64,dimension(:,:),allocatable :: volume
  sll_real64,dimension(:,:),allocatable :: surfx1,surfx2
  sll_real64, dimension(:),pointer  :: p
  sll_real64, dimension(:),pointer  :: M_diag,M_low,M_sup
  sll_real64, dimension(:),pointer  :: M1_diag,M1_low,M1_sup
  sll_real64, dimension(:),pointer  :: Av1_diag,Av1_low,Av1_sup
  sll_real64, dimension(:),pointer  :: Av2_diag,Av2_low,Av2_sup
  sll_real64, dimension(:),pointer  :: Bv1p_diag,Bv1p_low,Bv1p_sup
  sll_real64, dimension(:),pointer  :: Bv1m_diag,Bv1m_low,Bv1m_sup
  sll_real64, dimension(:),pointer  :: Bv2p_diag,Bv2p_low,Bv2p_sup
  sll_real64, dimension(:),pointer  :: Bv2m_diag,Bv2m_low,Bv2m_sup
  ! distribution functions at time steps n, star and n+1 
  ! communications are needed only in the x3 direction
  sll_real64, dimension(:,:,:,:), pointer     :: fn_v1v2
  !attention, it's a local variable
  sll_real64, dimension(:,:,:,:), pointer     :: fn_star_v1v2
  sll_real64, dimension(:,:,:,:), pointer     :: fnp1_v1v2
  sll_real64, dimension(:,:,:,:), pointer     :: dtfn_v1v2
  ! charge density
  sll_real64, dimension(:,:), allocatable     :: rho_x1
  sll_real64, dimension(:,:), allocatable     :: rho_split
  ! potential 
  sll_real64, dimension(:,:), allocatable     :: phi_x1
  sll_real64, dimension(:,:), allocatable     :: phi_split
  type(layout_4d),pointer :: sequential_v1v2_layout
  type(layout_2d),pointer :: phi_seq_x1_layout
  type(layout_2d),pointer :: split_rho_layout
  type(remap_plan_2D_real64), pointer :: split_to_seqx1
  type(remap_plan_2D_real64), pointer :: seqx1_to_split
contains
  procedure, pass(sim) :: run => run_vp_cart
  procedure, pass(sim) :: init_from_file => init_vp_cart

end type sll_simulation_4d_vp_eulerian_cartesian_finite_volume

interface sll_delete
  module procedure delete_vp_cart
end interface sll_delete

contains

subroutine init_vp_cart( sim, filename )
 intrinsic :: trim
 class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
 character(len=*), intent(in)                                :: filename
 sll_int32             :: IO_stat
 sll_real64            :: dt
 sll_int32             :: number_iterations
 sll_int32             :: num_cells_x1
 sll_int32             :: num_cells_x2
 sll_int32             :: num_cells_x3
 sll_int32             :: num_cells_x4
 sll_int32, parameter  :: input_file = 99

 namelist /sim_params/ dt, number_iterations
 namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3, num_cells_x4
 ! Try to add here other parameters to initialize the mesh values like
 ! xmin, xmax and also for the distribution function initializer.
 open(unit = input_file, file=trim(filename),IOStat=IO_stat)
 if( IO_stat /= 0 ) then
    print *, 'init_vp_cart() failed to open file ', filename
    STOP
 end if
 read(input_file, sim_params)
 read(input_file,grid_dims)
 close(input_file)

 sim%dt = dt
 sim%num_iterations = number_iterations
 ! In this particular simulation, since the system is periodic, the number
 ! of points is the same as the number of cells in all directions.
 !sim%nc_x1 = num_cells_x1
 !sim%nc_x2 = num_cells_x2
 !sim%nc_x3 = num_cells_x3
 !sim%nc_x4 = num_cells_x4
end subroutine init_vp_cart


subroutine initialize_vp4d( &
    sim, &
    mesh2dx, &
    mesh2dv, &
    tx,tv, &
    init_func, &
    params, &
    tmax )

 type(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout)     :: sim
 type(sll_logical_mesh_2d), pointer                    :: mesh2dx,mesh2dv
 class(sll_coordinate_transformation_2d_base),pointer          :: tx,tv
 procedure(sll_scalar_initializer_4d)                  :: init_func
 sll_real64, dimension(:), target                      :: params
 sll_real64 :: tmax
 sll_int32 :: ierr
!!$   sll_real64,dimension(:),pointer :: vx_mil
!!$   sll_real64,dimension(:),pointer :: vy_mil
!!$   sll_real64,dimension(:),pointer :: x_mil
!!$   sll_real64,dimension(:),pointer :: y_mil
!!$   sll_int32 :: ix
!!$   sll_int32 :: iy
!!$   sll_int32 :: ivx
!!$   sll_int32 :: ivy
!!$   sll_int32 :: mm
 sim%mesh2dx  => mesh2dx
 sim%mesh2dv  => mesh2dv
 sim%tx => tx
 sim%tv => tv
 sim%init_func => init_func
 sim%params    => params

 sim%tmax = tmax


end subroutine initialize_vp4d


! Note that the following function has no local variables, which is silly...
! This just happened since the guts of the unit test were transplanted here
! directly, but this should be cleaned up.
subroutine run_vp_cart(sim)
 class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
 logical :: exist
  type(sll_time_mark)  :: t0,t1
 sll_int32  :: file_id_1
 sll_int32  :: file_id_2
 sll_int32  :: file_id_3
 sll_int32  :: file_id_4
 sll_int32  :: loc_sz_v1
 sll_int32  :: loc_sz_v2
 sll_int32  :: loc_sz_x1
 sll_int32  :: loc_sz_x2
 sll_int32  :: i
 sll_int32  :: j
 sll_int32  :: k
 sll_int32  :: l
 sll_int32 :: count1, count2, count3, count4
 sll_real64 :: vmin,E2norm_ex
 sll_real64 :: vmax,time,time1
 sll_real64 :: dv
 sll_int32  :: ierr
 sll_int32  :: itime
 sll_int32  :: istat
 sll_int32  :: ic
 sll_int32  :: jc
 sll_real64 :: df,x
 sll_real64 ::erreurL2
 sll_real64 ::erreurL2_G
 sll_real64,dimension(:),pointer :: node,xmil 
 sll_real64,dimension(:,:),allocatable :: plotf2d_c1
 sll_real64,dimension(:,:),allocatable :: plotf2d_c2
 sll_real64,dimension(:,:),allocatable :: plotrho_split
 sll_real64,dimension(:,:),allocatable :: plotphi2d
 sll_real64,dimension(:,:),allocatable :: f_x_exact,f_vx_exact
 sll_real64,dimension(:,:),allocatable :: f_y_exact,f_vy_exact
 sll_real64,dimension(:,:),allocatable :: f_x_exact2,f_x_num
 sll_real64,dimension(:),pointer :: ww,w1
 sll_real64,dimension(:,:),pointer :: energ
 sll_real64,dimension(:),pointer :: vx_mil
 sll_real64,dimension(:),pointer :: vy_mil
 sll_real64,dimension(:),pointer :: x_mil
 sll_real64,dimension(:),pointer :: y_mil
 sll_real64,dimension(:,:),allocatable :: err
 sll_int32 :: ix
 sll_int32 :: iy
 sll_int32 :: ivx
 sll_int32 :: ivy
 sll_int32 :: ii,jj,mm
 sll_real64 :: t,v1,v2
 sll_real64 :: xref, yref,vxref,phi1
 sll_real64,dimension(2,2) :: jacob,invjacob
 sll_int32 :: iploc,ib1
 sll_int32,dimension(4)  :: global_indices
 sll_real64,dimension(1:2,1:2) :: jac_m,inv_jac
 sll_real64 :: det
 sll_real64 :: deltav,emax
 sll_real64 :: x1, x2,Ex,Ey
 sll_int32 :: icL,icR,jcL,jcR
#define BUFFER_SIZE 2
 sll_real64, dimension (BUFFER_SIZE) :: buffer
 sll_real64, dimension (BUFFER_SIZE) :: buffer_result
 sll_real64, dimension (BUFFER_SIZE) :: num_particles_local
 sll_real64, dimension (BUFFER_SIZE) :: num_particles_global
 sll_int32 :: buffer_counter
 sim%world_size = sll_get_collective_size(sll_world_collective)  
 sim%my_rank    = sll_get_collective_rank(sll_world_collective)  

 ! allocate the layouts...
 sim%sequential_v1v2_layout  => new_layout_4D( sll_world_collective )
 sim%phi_seq_x1_layout       => new_layout_2D( sll_world_collective )
 sim%split_rho_layout => new_layout_2D( sll_world_collective )

 sim%degree=sim%params(6)
 sim%nsch=sim%params(10)
 sim%test=sim%params(8)



 sim%nc_v1 = sim%mesh2dv%num_cells1
 sim%nc_v2 = sim%mesh2dv%num_cells2   
 sim%nc_x1 = sim%mesh2dx%num_cells1   
 sim%nc_x2 = sim%mesh2dx%num_cells2  

 sim%np_v1 = sim%degree * sim%nc_v1 + 1
 sim%np_v2 = sim%degree * sim%nc_v2 + 1

 sim%nproc_v1 = 1
 sim%nproc_v2 = 1
 if( (int(sqrt(real(sim%world_size))))**2/=sim%world_size ) then
    print *, 'sim%world_size must be a square number', sim%world_size
    STOP
 end if
 !write(*,*) 'test ', sqrt(real(9))
 !stop
 sim%nproc_x1 = int(sqrt(real(sim%world_size)))
 sim%nproc_x2 = int(sqrt(real(sim%world_size)))
 ! init the layout for the distribution function
 ! the mesh is split on the x1 et x2  direction
 call initialize_layout_with_distributed_4D_array( &
      sim%np_v1, &
      sim%np_v2, &
      sim%nc_x1, & !avance il n'y a pas +1
      sim%nc_x2, & !avance il n'y a pas +1
      sim%nproc_v1, &
      sim%nproc_v2, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%sequential_v1v2_layout)

 ! Allocate the array needed to store the local chunk of the distribution
 ! function data. First compute the local sizes.
 call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
      loc_sz_v1, &
      loc_sz_v2, &
      loc_sz_x1, &
      loc_sz_x2 )
 sim%buf_size=loc_sz_v2*loc_sz_v1*loc_sz_x1*loc_sz_x2
 !initialize for the comm
 sim%comm => new_comm_real64(sll_world_collective,4, sim%buf_size )
 call sll_configure_comm_real64_torus_2d( sim%comm,sim%nproc_x1,sim%nproc_x2 )
 !write(*,*) 'verify 2 loc_sz_x1',loc_sz_x1
 !write(*,*) 'verify 2 loc_sz_x2',loc_sz_x2
 !stop
 !write(*,*) 'sim%np_v1',sim%np_v1,loc_sz_v1
 !print *, 'just about to allocate dtfn_v1v2:', loc_sz_x1
 ! iz=0 and iz=loc_sz_x2+1 correspond to ghost cells.
 SLL_ALLOCATE(sim%fn_v1v2(loc_sz_v1,loc_sz_v2,loc_sz_x1,loc_sz_x2),ierr)
 !just essayer car c'est une variable local de RK4
 SLL_ALLOCATE(sim%fn_star_v1v2(loc_sz_v1,loc_sz_v2,loc_sz_x1,loc_sz_x2),ierr)
 SLL_ALLOCATE(sim%fnp1_v1v2(loc_sz_v1,loc_sz_v2,loc_sz_x1,loc_sz_x2),ierr)
 SLL_ALLOCATE(sim%dtfn_v1v2(loc_sz_v1,loc_sz_v2,loc_sz_x1,loc_sz_x2),ierr)

 SLL_ALLOCATE(ww(sim%np_v1*sim%np_v2),ierr)
 SLL_ALLOCATE(w1(sim%np_v1*sim%np_v2),ierr)
 SLL_ALLOCATE(energ(loc_sz_x1,loc_sz_x2),ierr)


 ! iz=0 corresponds to the mean values of rho and phi
 !write(*,*) 'verified loc_sz_x1',loc_sz_x1
 !write(*,*) 'verify loc_sz_x2',loc_sz_x2 
 SLL_ALLOCATE(sim%rho_split(loc_sz_x1,loc_sz_x2),ierr)
 !write(*,*) 'size rho_split',size(sim%rho_split(1,:))
 SLL_ALLOCATE(sim%phi_split(loc_sz_x1,loc_sz_x2),ierr)
 !write(*,*) sim%my_rank,'taille ',sim%nc_x1, sim%nc_x2, sim%nproc_x1, sim%nproc_x2
 call initialize_layout_with_distributed_2D_array( &
      sim%nc_x1, &
      sim%nc_x2, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%split_rho_layout)

 sim%nproc_x1 = 1
 sim%nproc_x2 = sim%world_size
!!$sim%nproc_x1 =1
!!$ sim%nproc_x2 = sim%world_size
 ! potential layout
 !write(*,*) ' sim%nproc_x1', sim%nproc_x1
 !write(*,*) ' sim%nproc_x2', sim%nproc_x2
 call initialize_layout_with_distributed_2D_array( &
      sim%nc_x1, &
      sim%nc_x2, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%phi_seq_x1_layout)
 call compute_local_sizes_2d( sim%phi_seq_x1_layout, loc_sz_x1, loc_sz_x2)

!!$write(*,*) sim%my_rank,'taille ',sim%nc_x1, sim%nc_x2, sim%nproc_x1, sim%nproc_x2
!!$write(*,*) sim%my_rank,'size ',loc_sz_x1,loc_sz_x2
!!$ call initialize_layout_with_distributed_2D_array( &
!!$      sim%nc_x1, &
!!$      sim%nc_x2, &
!!$      sim%nproc_x1, &
!!$      sim%nproc_x2, &
!!$      sim%split_rho_layout)
 !apres changer et faire inverse
 SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2),ierr)
 !write(*,*) 'size rho_x1',size(sim%rho_x1(1,:))
 SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
 !write(*,*) 'taille de rho_x1', size(sim%rho_x1(1,:)),size(sim%rho_x1(:,1))
 !write(*,*) 'arrive ici'
 sim%poisson_plan=>new_poisson_2d_periodic_plan_cartesian_par_alt( &
      sim%phi_seq_x1_layout, &
      sim%nc_x1, &
      sim%nc_x2, &
      sim%mesh2dx%eta1_max-sim%mesh2dx%eta1_min, &
      sim%mesh2dx%eta2_max-sim%mesh2dx%eta2_min )
 !write(*,*) 'arrive ici2'




!!!!!!!!!!
 !SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
 !print *, 'rank = ', sim%my_rank, 'local sizes of phi: ', loc_sz_x1, loc_sz_x2


 ! initialize here the distribution function

 ! the function is passed by the user when the init_vp subroutine is called.
 ! The routine sll_4d_parallel_array_initializer_cartesian is in 
 ! src/parallal_array_initializers/sll_parallel_array_initializer_module.F90
 ! the particular initializer is in
 ! parallel_array_initializers/sll_common_array_initializers_module.F90

!!$    call sll_4d_parallel_array_initializer_cartesian( &
!!$         sim%sequential_v1v2, &
!!$         sim%mesh4d, &
!!$         sim%fn_v1v2(:,:,:,1:loc_sz_x2), &
!!$         sim%init_func, &
!!$         sim%params)
 call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
      loc_sz_v1, &
      loc_sz_v2, &
      loc_sz_x1, &
      loc_sz_x2 )
 !call sll_view_lims_4D(sim%sequential_v1v2_layout)
!!$ write(*,*) 'taille fn_v1v2_1',size(sim%fn_v1v2(:,:,:,1:loc_sz_x2),1)
!!$ write(*,*) 'taille fn_v1v2_2',size(sim%fn_v1v2(:,:,:,1:loc_sz_x2),2)
!!$ write(*,*) 'taille fn_v1v2_3',size(sim%fn_v1v2(:,:,:,1:loc_sz_x2),3)
!!$ write(*,*) 'taille fn_v1v2_4',size(sim%fn_v1v2(:,:,:,1:loc_sz_x2),4)
!!$ write(*,*) 'loc_sz_x2',loc_sz_x2
 call sll_4d_parallel_array_initializer_finite_volume( &
      sim%sequential_v1v2_layout, &
      sim%mesh2dv, &
      sim%mesh2dx, &
      sim%fn_v1v2(:,:,:,1:loc_sz_x2), &
      sim%init_func , &
      sim%params, &
      sim%tv, &
      sim%tx, &
      sim%degree, &
      sim%degree)


 call velocity_mesh_connectivity(sim)


 !


!!$    do i=1,loc_sz_x1
!!$       do j=1,loc_sz_x2
!!$          sim%rho_x1(i,j)=0_f64
!!$       enddo
!!$    enddo



 ! With the distribution function initialized in at least one configuration,
 ! we can proceed to carry out the computation of the electric potential.
 ! First we need to compute the charge density. Some thoughts:
 !
 ! The computation of rho is a reduction process that takes as input a 4d
 ! array and that should return a 2d array (or alternatively, a 4d array
 ! of size 1 in the reduced directions). For example, a df of dimensions
 ! np1 X np2 X np3 X np4 might effectively end up as an array of dimensions
 ! np1 X np2 X np3 X 1  after a summation of all the values in x4. After a
 ! second summation along x3, the dimensions would be np1 X np2 X 1  X 1. 
 ! One simple-minded but inefficient way to prepare the data for the double
 ! reduction could be to have a layout in a process mesh NP1 X NP2 X 1 X 1
 ! where NP1xNP2 is the total number of processors available. The problem 
 ! here is that the end result would still need a layout change to be fed
 ! into the Poisson solver...
 !
 ! So can we do better? Let try working backwards. The desired input for
 ! the Poisson solver is a 2d array where sequential operations in x1 are
 ! possible. Hence, the last reduction operation
 !
 ! Let's start with the idea that we want to maintain the same number of
 ! processors busy when we launch the Poisson solver. This means that if the
 ! original process mesh for the 4d data is NP1xNP2xNP3x1 (we start with a
 ! layout that permits a reduction in x4), then the processor mesh for the
 ! Poisson step should be NP1'xNP2'x1x1 where NP1xNP2xNP3 = NP1'xNP2'


 ! plot the int function 



 ! mpi communications parameters
!!$    ranktop=mod(sim%my_rank+1,sim%world_size)
!!$    rankbottom=sim%my_rank-1
!!$    if (rankbottom.lt.0) rankbottom=sim%world_size-1
!!$    message_id=1
!!$    datasize=loc_sz_v1*loc_sz_v2*loc_sz_x1
!!$    datasizephi=loc_sz_x1
!!$    tagtop=sim%my_rank
!!$    tagbottom=ranktop


 ! init the volume and surface arrays
 SLL_ALLOCATE(sim%volume(loc_sz_x1,loc_sz_x2),ierr)
 ! vertical edges 
 SLL_ALLOCATE(sim%surfx1(loc_sz_x1+1,loc_sz_x2),ierr)
 ! horizontal edges 
 SLL_ALLOCATE(sim%surfx2(loc_sz_x1,loc_sz_x2+1),ierr)
!!$    write(*,*) 'sim%mesh2dx%eta1_min', sim%mesh2dx%eta1_min
!!$    write(*,*) 'sim%mesh2dx%eta1_max', sim%mesh2dx%eta1_max
!!$    write(*,*) 'sim%mesh2dx%num_cells1', sim%mesh2dx%num_cells1
!!$    write(*,*) 'sim%mesh2dx%delta_eta1', sim%mesh2dx%delta_eta1
!!$    stop
!!$    ! simple computation
!!$    ! to be generalized with generic transformation...
!!$    volume=sim%mesh2dx%delta_eta1 * &
!!$         sim%mesh2dx%delta_eta2 

 global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, (/1,1,1,1/) )
 ! cell volumes, check this...
 !print *, 'what is the size of loc_sz_x2??? ', loc_sz_x2
 do j=1,loc_sz_x2
    do i=1,loc_sz_x1 ! loc_sz_x1 is number of points, we're processing cells
       ic=i+global_indices(3)-1
       jc=j+global_indices(4)-1
       sim%volume(i,j)=cell_volume( sim%tx, ic, jc,3)          
    end do
 end do

 ! cell volumes, check this...
!!$ do j=1,loc_sz_x2
!!$    do i=1,loc_sz_x1-1 ! loc_sz_x1 is number of points, we're processing cells
!!$       global_indices(1:4) =  &
!!$            local_to_global_4D(sim%sequential_v1v2, (/1,1,i,j/) )
!!$       sim%volume(i,j) = &
!!$            cell_volume( sim%tx, global_indices(3),global_indices(4),3)
!!$    end do
!!$ end do

!!$    surface(1:loc_sz_x1*loc_sz_x2)= &
!!$         sim%mesh2dx%delta_eta1
!!$    surface(loc_sz_x1*loc_sz_x2+1:2*loc_sz_x1*loc_sz_x2)= &
!!$         sim%mesh2dx%delta_eta2

 do j=1,loc_sz_x2
    ic=1+global_indices(3)-1
    jc=j+global_indices(4)-1  
    sim%surfx1(1,j)=edge_length_eta1_minus( sim%tx, ic, jc,3)      
    do i=1,loc_sz_x1
       ic=i+global_indices(3)-1
       sim%surfx1(i+1,j)=edge_length_eta1_plus( sim%tx, ic, jc,3)  
    end do
 end do

 do i=1,loc_sz_x1
    ic=i+global_indices(3)-1
    jc=1+global_indices(4)-1  
    sim%surfx2(i,1)=edge_length_eta2_minus( sim%tx, ic, jc,3)      
    do j=1,loc_sz_x2
       jc=j+global_indices(4)-1
       sim%surfx2(i,j+1)=edge_length_eta2_plus( sim%tx, ic, jc,3)  
    end do
 end do

 !plan
 !write(*,*) 'planter ici 1'

 sim%split_to_seqx1 => &
      NEW_REMAP_PLAN(sim%split_rho_layout,sim%phi_seq_x1_layout,sim%rho_split)
 !write(*,*) 'planter ici 2'
 sim%seqx1_to_split => &
      NEW_REMAP_PLAN(sim%phi_seq_x1_layout,sim%split_rho_layout, sim%phi_x1)
 !write(*,*) 'planter ici 3'


 !write(*,*) 'vertical',sim%my_rank,sum(sim%surfx1(1,:))
 !write(*,*) 'horizontal',sim%my_rank,sum(sim%surfx2(1,:))
 !verify the calcul of rho
 !the initial electric field
 SLL_ALLOCATE(vx_mil(sim%np_v1),ierr)
 SLL_ALLOCATE(vy_mil(sim%np_v2),ierr)
 SLL_ALLOCATE(x_mil(sim%nc_x1),ierr)
 SLL_ALLOCATE(y_mil(sim%nc_x2),ierr)
 do ivx=1,sim%np_v1
    vx_mil(ivx)=sim%mesh2dv%eta1_min+(ivx-1)*sim%mesh2dv%delta_eta1/sim%degree
 enddo
 do ivy=1,sim%np_v2
    vy_mil(ivy)=sim%mesh2dv%eta2_min+(ivy-1)*sim%mesh2dv%delta_eta2/sim%degree
 enddo
 do ix=1,sim%nc_x1
    x_mil(ix)=sim%mesh2dx%eta1_min+ix*sim%mesh2dx%delta_eta1- &
         sim%mesh2dx%delta_eta1/2
 enddo
 do iy=1,sim%nc_x2
    y_mil(iy)=sim%mesh2dx%eta2_min+iY*sim%mesh2dx%delta_eta2- &
         sim%mesh2dx%delta_eta2/2
 enddo

 ! time loop
 t=0
 buffer_counter =1
 !compute the time step
 sim%cfl=sim%params(7)
!!$    write(*,*) 'Vxmax = ', sim%mesh2dv%eta1_max
!!$    write(*,*) 'deltav = ', sim%mesh2dx%delta_eta1
!!$    write(*,*) 'surf/volum =  ', 2*(sim%surfx1(1,1)+sim%surfx2(1,1))/ &
!!$         sim%volume(1,1)*sim%mesh2dx%delta_eta1
!!$    stop
 ! write(*,*) 'coucou'
 ! space cfl condition
 sim%dt = sim%cfl*sim%volume(1,1)/2/(sim%surfx1(1,1)+sim%surfx2(1,1))/ &
      max(sim%mesh2dv%eta1_max,sim%mesh2dv%eta2_max,abs(sim%mesh2dv%eta1_min), &
      abs(sim%mesh2dv%eta2_min))
 ! velocity cfl condition
 deltav=min(sim%mesh2dv%delta_eta1,sim%mesh2dv%delta_eta2)
 emax=sim%params(9)
 sim%dt=min(sim%dt,sim%cfl*deltav/emax/(sim%degree+1)) 

!!$    sim%dt=0.1
 write(*,*) 'dt = ', sim%dt
 !stop
 itime=0
 sim%params(11)=t
 do while(t.lt.sim%tmax)
    itime=itime+1
    sim%Enorm = 0.0_f64

    !compute rho_split c-a-d calculer avec parallel sur 2 directions

    sim%rho_split = 0.0_f64
    do i=1,loc_sz_x1
       do j=1,loc_sz_x2
          do ii=1,loc_sz_v1
             do jj=1,loc_sz_v2
                mm=loc_sz_v1*(jj-1)+ii
                sim%rho_split(i,j) = &
                     sim%rho_split(i,j)+sim%fn_v1v2(ii,jj,i,j)*sim%p(mm)
             enddo
          end do
       enddo
    enddo

!!$    !if (sim%my_rank==0) then
!!$       do i=1,loc_sz_x1
!!$          do j=1,loc_sz_x2
!!$             write(*,*) ' sim%rho_split',sim%my_rank,i, j,  sim%rho_split(i,j)
!!$          end do
!!$        end do
!!$    !end if
    !write(*,*) 'probleme ici 3'
    !allouer pour rho_x1 et phi_x1

    !ordonner rho_split sur x1 pour utiliser le Poisson solver.
    !call compute_local_sizes_2d( sim%phi_seq_x1_layout, loc_sz_x1, loc_sz_x2)
    !call sll_view_lims_2D(sim%split_rho_layout)
    !write(*,*) 'taille  sim%rho_x1',size( sim%rho_x1,1),size( sim%rho_x1,2)
    !write(*,*) sim%my_rank, 'taille  sim%rho_split',size( sim%rho_split,1),size( sim%rho_split,2)
    !call set_time_mark(t0)
    call apply_remap_2D(sim%split_to_seqx1, sim%rho_split,sim%rho_x1)
!!$if (sim%my_rank==0) then
!!$write(*,*) ' sim%rho_x1(1,1)=', sim%rho_x1(1,1)
!!$end if
    !stop
    !write(*,*) 'probleme ici 4'
!!$   do i=1,loc_sz_x1
!!$    write(*,*) 'rho_x1_temp (',i,',:) = ', rho_x1_temp(i,:)
!!$   end do
    !write()
!!$    do i=2,loc_sz_x1
!!$       do j=1,loc_sz_x2
!!$          sim%rho_x1(i,j) = (rho_x1_temp(i-1,j)+rho_x1_temp(i-1,j-1)+rho_x1_temp(i,j)+rho_x1_temp(i,j-1))/4
!!$       enddo
!!$    enddo
!!$    do j=1,loc_sz_x2+1
!!$       sim%rho_x1(1,j) = (rho_x1_temp(loc_sz_x1,j)+rho_x1_temp(loc_sz_x1,j-1)+rho_x1_temp(1,j)+rho_x1_temp(1,j-1))/4
!!$       sim%rho_x1(loc_sz_x1+1,j) =  sim%rho_x1(1,j)
!!$    enddo
!!$   do i=1,loc_sz_x1+1
!!$    write(*,*) 'rho (',i,',:) = ', sim%rho_x1(i,:)
!!$   end do

    ! solve the poisson equation
    !maillage logical
!!$    do i=1,loc_sz_x1
!!$       write (*,*) 'rho = ', sim%rho_x1(i,:)
!!$    enddo
    !stop
    sim%phi_x1=0.0_f64
!!$    write(*,*) 'taille de rho_x1', size(sim%rho_x1(1,:)),size(sim%rho_x1(:,1))
!!$    stop
    call compute_local_sizes_2d( sim%phi_seq_x1_layout, loc_sz_x1, loc_sz_x2)
    !write(*,*) 'planter poisson'

    call solve_poisson_2d_periodic_cartesian_par_alt(sim%poisson_plan, &
         sim%rho_x1, &
         sim%phi_x1(:,1:loc_sz_x2))

!!$    !if (sim%my_rank==0) then
!!$       do i=1,loc_sz_x1
!!$          do j=1,loc_sz_x2
!!$             write(*,*) ' sim%phi_x1',sim%my_rank,i, j,  sim%phi_x1(i,j)
!!$          end do
!!$        end do
!!$    !end if
    !sim%phi_x1(:,0)=sim%phi_x1(:,loc_sz_x2)
    !sim%phi_x1(:,loc_sz_x2+1)=sim%phi_x1(:,1)
!!$write(*,*) 'phi avance = ', sim%phi_x1(:,1)
!!$      sim%phi_x1=-sim%phi_x1
!!$write(*,*) 'phi apres = ', sim%phi_x1(:,1)
    !attention the sign of phi after the solver Poisson
    if((sim%test.ne.11).and.(sim%test.ne.12)) then
       sim%phi_x1=-sim%phi_x1
    end if
     ! sim%phi_x1=-sim%phi_x1
    !revient dans le split layout pour phi
    !write(*,*) sim%my_rank, 'here 1'
    call apply_remap_2D( sim%seqx1_to_split, sim%phi_x1, sim%phi_split)
    !time = time_elapsed_since(t0)
    !print *, 'time of solveur poisson is  : ',time
!!$    do i=1,loc_sz_x1
!!$       write (*,*)  'x',x_mil(i)-sim%mesh2dx%delta_eta1/2,'phi = ', sim%phi_x1(i,:)
!!$    enddo
    !stop
    !revenir dans le bon valeur du loc_sz
    !write(*,*) sim%my_rank, 'here 2'
    !write(*,*) 'taille of phi _split',size(sim%phi_split)
    call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
         loc_sz_v1, &
         loc_sz_v2, &
         loc_sz_x1, &
         loc_sz_x2 )
    t=t+sim%dt

 !call set_time_mark(t1)
    if (sim%nsch == 0) then
       call euler(sim)
    elseif (sim%nsch == 1) then
       call RK2(sim)
    elseif (sim%nsch == 2) then
       call RK4(sim)
    end if
  !time1 = time_elapsed_since(t1)
  !print *, 'time of schema is  : ',time1

  call sll_new_file_id(file_id_4,ierr)

    if((sim%test==1).or.(sim%test==9).or.(sim%test==5).or.(sim%test==10) &
         .or.(sim%test==11).or.(sim%test==12)) then
       sim%buf1 => get_buffer(sim%comm,1) 
       do k=1,loc_sz_x2
          sim%buf1(k)=sim%phi_split(1,k)
       end do
       call comm_send_real64(sim%comm,1,loc_sz_x2)

       sim%buf2 => get_buffer(sim%comm,2) 
       do k=1,loc_sz_x2
          sim%buf2(k)=sim%phi_split(loc_sz_x1,k)
       end do
       call comm_send_real64(sim%comm,2,loc_sz_x2)

       sim%buf3 => get_buffer(sim%comm,3) 
       do k=1,loc_sz_x1
          sim%buf3(k)=sim%phi_split(k,1)
       end do
       call comm_send_real64(sim%comm,3,loc_sz_x1)

       sim%buf4 => get_buffer(sim%comm,4) 
       do k=1,loc_sz_x1
          sim%buf4(k)=sim%phi_split(k,loc_sz_x2)
       end do
       call comm_send_real64(sim%comm,4,loc_sz_x1)
       do ic=2,loc_sz_x1-1
          do jc=2,loc_sz_x2-1
             icL=ic-1
             icR=ic+1
             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
                  (/1,1,1,1/) )
             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
                  sim%mesh2dx%delta_eta1
             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
                  sim%mesh2dx%delta_eta2
             jac_m=sim%tx%jacobian_matrix(x1,x2)
             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
             Ex=-(sim%phi_split(icR,jc)-sim%phi_split(icL,jc))/2/ &
                  sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(ic,jc+1)- &
                  sim%phi_split(ic,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
             Ey=-(sim%phi_split(ic,jc+1)-sim%phi_split(ic,jc-1))/2/ &
                  sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,jc)- &
                  sim%phi_split(icL,jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
             det=sim%tx%jacobian(x1,x2)
             if(sim%test==2) then
                Ex=1.0_f64
                Ey=0.0_f64
             endif
             if(sim%test==3) then
                Ex=0.0_f64
                Ey=1.0_f64
             endif
             sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
                  sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
             !write(*,*) 'pro mpi 1', sim%Enorm

          end do
       end do


       call comm_receive_real64(sim%comm,3,count3)
       call comm_receive_real64(sim%comm,4,count4)
       if ((count3.ne.loc_sz_x1).or.(count4.ne.loc_sz_x1)) then
          write(*,*) 'problem avec send de mpi'
          stop
       endif
       sim%buf3 => get_buffer(sim%comm,3)
       sim%buf4 => get_buffer(sim%comm,4)
       !pour jc=1 et Loc_sz_x2
       do ic=2,loc_sz_x1-1
          icL=ic-1
          icR=ic+1
          global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
               (/1,1,1,1/) )
          x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
               sim%mesh2dx%delta_eta1
          x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
               sim%mesh2dx%delta_eta2
          jac_m=sim%tx%jacobian_matrix(x1,x2)
          inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)

          Ex=-(sim%phi_split(icR,1)-sim%phi_split(icL,1))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(ic,2)- &
               sim%buf3(ic))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
          !write(*,*) 'Ex',ic, jc,  Ex
          Ey=-(sim%phi_split(ic,2)-sim%buf3(ic))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,1)- &
               sim%phi_split(icL,1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
          !write(*,*) 'pro mpi 2', sim%Enorm
          Ex=-(sim%phi_split(icR,loc_sz_x2)-sim%phi_split(icL,loc_sz_x2))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf4(ic)- &
               sim%phi_split(ic,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)

          Ey=-(sim%buf4(ic)-sim%phi_split(ic,loc_sz_x2-1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,loc_sz_x2)- &
               sim%phi_split(icL,loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
       end do
       call comm_receive_real64(sim%comm,1,count1)
       call comm_receive_real64(sim%comm,2,count2)
       if ((count1.ne.loc_sz_x2).or.(count2.ne.loc_sz_x2)) then
          write(*,*) 'problem avec send de mpi'
          stop
       endif
       sim%buf1 => get_buffer(sim%comm,1)
       sim%buf2 => get_buffer(sim%comm,2)
          do jc=2,loc_sz_x2-1
             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
                  (/1,1,1,1/) )
             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
                  sim%mesh2dx%delta_eta1
             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
                  sim%mesh2dx%delta_eta2
             jac_m=sim%tx%jacobian_matrix(x1,x2)
             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
             Ex=-(sim%phi_split(2,jc)-sim%buf1(jc))/2/ &
                  sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(1,jc+1)- &
                  sim%phi_split(1,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)

             if (abs(Ex).gt.100) then
                   write(*,*) 'ou',sim%my_rank,1, jc, Ex
             end if
             Ey=-(sim%phi_split(1,jc+1)-sim%phi_split(1,jc-1))/2/ &
                  sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,jc)- &
                  sim%buf1(jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)  

             det=sim%tx%jacobian(x1,x2)
             if(sim%test==2) then
                Ex=1.0_f64
                Ey=0.0_f64
             endif
             if(sim%test==3) then
                Ex=0.0_f64
                Ey=1.0_f64
             endif
             sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
                  sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
          end do
          !pour ic=loc_sz_x1
          do jc=2,loc_sz_x2-1
             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
                  (/1,1,1,1/) )
             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
                  sim%mesh2dx%delta_eta1
             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
                  sim%mesh2dx%delta_eta2
             jac_m=sim%tx%jacobian_matrix(x1,x2)
             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
             Ex=-(sim%buf2(jc)-sim%phi_split(loc_sz_x1-1,jc))/2/ &
                  sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(loc_sz_x1,jc+1)- &
                  sim%phi_split(loc_sz_x1,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
            
             Ey=-(sim%phi_split(loc_sz_x1,jc+1)-sim%phi_split(loc_sz_x1,jc-1))/2/ &
                  sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(jc)- &
                  sim%phi_split(loc_sz_x1-1,jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
             det=sim%tx%jacobian(x1,x2)
             if(sim%test==2) then
                Ex=1.0_f64
                Ey=0.0_f64
             endif
             if(sim%test==3) then
                Ex=0.0_f64
                Ey=1.0_f64
             endif
             sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
                  sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)

          end do
          !pour ic=1 et jc=1
             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
                  (/1,1,1,1/) )
             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
                  sim%mesh2dx%delta_eta1
             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
                  sim%mesh2dx%delta_eta2
             jac_m=sim%tx%jacobian_matrix(x1,x2)
             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
          Ex=-(sim%phi_split(2,1)-sim%buf1(1))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(1,2)- &
               sim%buf3(1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)

          Ey=-(sim%phi_split(1,2)-sim%buf3(1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,1)- &
               sim%buf1(1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
          
!!$          !******************REPARER
!!$          Ey=0
!!$          if(sim%my_rank==0) then
!!$             write(*,*) 'Ex at point (1,1) = ', Ex
!!$             write(*,*) 'phi_split (1,2) = ', sim%phi_split(1,2)
!!$             write(*,*) 'phi_split (2,1) = ', sim%phi_split(2,1)
!!$             write(*,*) 'phi_split (1,0) = ', sim%buf3(1)
!!$             write(*,*) 'phi_split (0,1) = ', sim%buf1(1)
!!$             write(*,*) ' sim%mesh2dx%delta_eta1',sim%mesh2dx%delta_eta1
!!$             write(*,*) ' sim%mesh2dx%delta_eta2',sim%mesh2dx%delta_eta2
!!$             write(*,*) 'inv_jac(1,1)',inv_jac(1,1)
!!$             write(*,*) 'inv_jac(2,1)',inv_jac(2,1)
!!$          end if

          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
!!$          if(sim%my_rank==0) then
!!$             write(*,*) 'Ex at point (1,1) = ', Ex
!!$             write(*,*) 'phi_split (1,2) = ', sim%phi_split(1,2)
!!$             write(*,*) 'phi_split (2,1) = ', sim%phi_split(2,1)
!!$             write(*,*) 'phi_split (1,0) = ', sim%buf3(1)
!!$             write(*,*) 'phi_split (0,1) = ', sim%buf1(1)
!!$             write(*,*) ' sim%mesh2dx%delta_eta1',sim%mesh2dx%delta_eta1
!!$             write(*,*) ' sim%mesh2dx%delta_eta2',sim%mesh2dx%delta_eta2
!!$             write(*,*) 'inv_jac(1,1)',inv_jac(1,1)
!!$             write(*,*) 'inv_jac(2,1)',inv_jac(2,1)
!!$             write(*,*) 'det',det
!!$          end if
          !write(*,*) 'Ey ', Ey, Ex
!!$          write(*,*) ' sim%mesh2dx%delta_eta1',sim%mesh2dx%delta_eta1
!!$          write(*,*) ' sim%mesh2dx%delta_eta2',sim%mesh2dx%delta_eta2
!!$          write(*,*) 'det',det
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
          ! stop
          !pour ic=1 et jc=loc_sz_x2

!!$             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$                  (/1,1,1,loc_sz_x2/) )
!!$             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta1
!!$             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta2
!!$             jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
          Ex=-(sim%phi_split(2,loc_sz_x2)-sim%buf1(loc_sz_x2))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf2(1)- &
               sim%phi_split(1,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
          !write(*,*) 'Ex at 5',  Ex
!!$          Ey=0
!!$!******************REPARER
          Ey=-(sim%buf4(1)-sim%phi_split(1,loc_sz_x2-1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,loc_sz_x2)- &
               sim%buf1(loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
          !write(*,*) 'pro mpi 1', sim%Enorm
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
!!$          if(sim%my_rank==0) then
!!$             write(*,*) 'check', 2, loc_sz_x2
!!$             write(*,*) 'Ex at point (2,1) = ', Ex, Ey
!!$             write(*,*) 'phi_split (2,2) = ', sim%phi_split(2,2)
!!$             write(*,*) 'phi_split (1,1) = ', sim%phi_split(1,1)
!!$             write(*,*) 'phi_split (1,0) = ', sim%buf2(1)
!!$             write(*,*) 'phi_split (0,2) = ', sim%buf1(2)
!!$             write(*,*) ' sim%mesh2dx%delta_eta1',sim%mesh2dx%delta_eta1
!!$             write(*,*) ' sim%mesh2dx%delta_eta2',sim%mesh2dx%delta_eta2
!!$             write(*,*) 'inv_jac(1,1)',inv_jac(1,1)
!!$             write(*,*) 'det',det
!!$             write(*,*) 'pro mpi 2', sim%Enorm
!!$          end if
          

!!$             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$                  (/1,1,1,1/) )
!!$             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta1
!!$             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta2
!!$             jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
          Ex=-(sim%buf2(1)-sim%phi_split(loc_sz_x1-1,1))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(loc_sz_x1,2)- &
               sim%buf3(loc_sz_x1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)

          Ey=-(sim%phi_split(loc_sz_x1,2)-sim%buf3(loc_sz_x1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(1)- &
               sim%phi_split(loc_sz_x1-1,1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
          !write(*,*) 'pro mpi 3', sim%Enorm
          !pour ic=loc_sz_x1 et jc=loc_sz_x2
!!$             global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$                  (/1,1,1,1/) )
!!$             x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta1
!!$             x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$                  sim%mesh2dx%delta_eta2
!!$             jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$             inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
          Ex=-(sim%buf2(loc_sz_x2)-sim%phi_split(loc_sz_x1-1,loc_sz_x2))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf4(loc_sz_x1)- &
               sim%phi_split(loc_sz_x1,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)

          Ey=-(sim%buf4(loc_sz_x1)-sim%phi_split(loc_sz_x1,loc_sz_x2-1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(loc_sz_x2)- &
               sim%phi_split(loc_sz_x1-1,loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
!!$          !******************REPARER
!!$          Ey=0
          det=sim%tx%jacobian(x1,x2)
          if(sim%test==2) then
             Ex=1.0_f64
             Ey=0.0_f64
          endif
          if(sim%test==3) then
             Ex=0.0_f64
             Ey=1.0_f64
          endif
          sim%Enorm=sim%Enorm + sim%mesh2dx%delta_eta1* &
               sim%mesh2dx%delta_eta2*det*(Ex**2+Ey**2)
          write(*,*) 'iter = ',itime, ' t = ', t ,' energy  = ', log(sqrt(sim%Enorm))
!!$          write(*,*) 'iter = ',itime, ' t = ', t ,' energy  = ', sim%Enorm

          buffer(buffer_counter) = sim%Enorm
          if(buffer_counter==BUFFER_SIZE) then
             call sll_collective_reduce_real64(sll_world_collective, &
                  buffer, &
                  BUFFER_SIZE, &
                  MPI_SUM, &
                  0, &
                  buffer_result )
             buffer_counter=1
             if (sim%my_rank==0) then
                open(file_id_4,file='energy',position='append')
                if(itime==BUFFER_SIZE) then 
                   rewind(file_id_4)
                endif

                buffer_result(:)=log(sqrt(buffer_result(:)))
                do i=1,BUFFER_SIZE
                   write(file_id_4,*) t, buffer_result(i)
                enddo
                close(file_id_4)
             end if
          else
             buffer_counter=buffer_counter+1
          end if
          !sim%Enorm = 0.0_f64
         end if
 end do

    write(*,*) 'number of iteration', itime
    write(*,*) 'final time ',t

!!$    if (sim%test==1) then
!!$       if (sim%my_rank==0) then
!!$          t=0.d0
!!$          open(699,file='asymptotic')
!!$          do while(t<sim%tmax)
!!$             call solexact(sim,t,E2norm_ex)
!!$             write(699,*) t,log(sqrt(E2norm_ex))
!!$             t=t+sim%dt
!!$          end do
!!$          close(699)
!!$       end if
!!$    end if



    call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
         loc_sz_v1, loc_sz_v2, loc_sz_x1, loc_sz_x2) 
    allocate (plotf2d_c1(loc_sz_x1,loc_sz_v1))
    do i = 1, loc_sz_x1
       do j = 1, loc_sz_v1
plotf2d_c1(i,j) = sim%fn_v1v2(j,1,i,1)
!write(*,*) 'check this ', sim%fn_v1v2(j,1,i,1)
          !plotf2d_c1(i,j) = max(0.0_f64,sim%fn_v1v2(j,1,i,1))
          if (plotf2d_c1(i,j).gt.100) then
             write(*,*) ' plotf2d_c1(i,j)',  plotf2d_c1(i,j)
          end if
       end do
    end do

    allocate (plotf2d_c2(loc_sz_x2,loc_sz_v2))
    do i = 1, loc_sz_x2
       do j = 1, loc_sz_v2
          plotf2d_c2(i,j) = max(0.0_f64,sim%fn_v1v2(1,j,1,i))
       end do
    end do
    !write(*,*) 'after plotf2d_c2'

!!$
!!$ allocate (xmil(loc_sz_x1))
!!$ allocate (node(loc_sz_v1))
!!$  call sll_new_file_id(file_id_3,ierr)
!!$ open(file_id_3,file='distribution')
!!$ do i=1,loc_sz_x1
!!$    xmil(i)=sim%mesh2dx%eta1_min+sim%mesh2dx%delta_eta1*i-sim%mesh2dx%delta_eta1/2
!!$ end do
!!$ do j=1,loc_sz_v1
!!$    node(j)=sim%mesh2dv%eta1_min+sim%mesh2dv%delta_eta1/sim%degree*(j-1)
!!$ end do
!!$ do i = 1, loc_sz_x1
!!$    do j = 1, loc_sz_v1
!!$       df = sim%fn_v1v2(j,1,i,1)
!!$       write(file_id_3,*) xmil(i), node(j), df 
!!$    end do
!!$ end do
!!$ close(file_id_3)
!!$ deallocate(xmil)
!!$ deallocate(node)
    !stop

!!$    write(*,*) 'plotf2d',plotf2d
!!$    stop

!!$    allocate (plotphi2d(loc_sz_x1,loc_sz_x2))
!!$
!!$          plotphi2d(i,j) = sim%phi_x1(i,j)
!!$       end do
!!$    end do

 global_indices(1:4) =  local_to_global_4D(sim%sequential_v1v2x1, (/1,1,1,1/) )
 write (*,*) 'Vxmax = ', sim%mesh2dv%eta1_max
 write (*,*) 'Vxmin = ', sim%mesh2dv%eta1_min
 call sll_gnuplot_rect_2d_parallel( &
      sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
      sim%mesh2dx%delta_eta1, &
      sim%mesh2dv%eta1_min+(global_indices(1)-1)*sim%mesh2dv%delta_eta1/sim%degree, &
      sim%mesh2dv%delta_eta1/sim%degree, &
      size(plotf2d_c1,1), &
      size(plotf2d_c1,2), &
      plotf2d_c1, &
      "plotf2d_c1", &
      0, &
      ierr)

 call sll_gnuplot_rect_2d_parallel( &
      sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
      sim%mesh2dx%delta_eta2, &
      sim%mesh2dv%eta2_min+(global_indices(2)-1)*sim%mesh2dv%delta_eta2/sim%degree, &
      sim%mesh2dv%delta_eta2/sim%degree, &
      size(plotf2d_c2,1), &
      size(plotf2d_c2,2), &
      plotf2d_c2, &
      "plotf2d_c2", &
      0, &
      ierr)

!!!!!!!!!!!!!!!!!!!simplify
!!$    allocate (f_exact(loc_sz_x1,loc_sz_v1))
!!$ do ix = 1, loc_sz_x1
!!$    do jvx = 1, loc_sz_v1
!!$       do iy = 1, loc_sz_x2
!!$          do jvy = 1, loc_sz_v2
!!$             x1=(ix-1)*sim%mesh2dx%delta_eta1+sim%mesh2dx%eta1_min
!!$             v1=sim%mesh2dv%eta1_min+(jvx-1)*sim%mesh2dv%delta_eta1/sim%degree
!!$             x2=(iy-1)*sim%mesh2dx%delta_eta2+sim%mesh2dx%eta2_min
!!$             v2=sim%mesh2dv%eta2_min+(jvy-1)*sim%mesh2dv%delta_eta2/sim%degree
!!$
!!$             sim%params(11)=t
!!$             f_exact(i,j)=sim%init_func(v1,v2,x1,x2,sim%params)
!!$          end do
!!$       end do
!!$    end do



    global_indices(1:4) =  local_to_global_4D(sim%sequential_v1v2_layout, (/1,1,1,1/) )
    !if(sim%test==1) then
       allocate (plotrho_split(loc_sz_x1,loc_sz_x2))
       do i = 1, loc_sz_x1
          do j=1,loc_sz_x2
             !plotrho_split(i,j) =  max(0.0_f64,sim%rho_split(i,j))
             plotrho_split(i,j) =  sim%rho_split(i,j)
          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
            sim%mesh2dx%delta_eta1, &
            sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
            sim%mesh2dx%delta_eta2, &
            plotrho_split, &
            "plotrho_split", &
            0, &
            ierr)
    !end if
    write (*,*) 'Vxmax = ', sim%mesh2dv%eta1_max
    write (*,*) 'Vxmin = ', sim%mesh2dv%eta1_min
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
         sim%mesh2dx%delta_eta1, &
         sim%mesh2dv%eta1_min+(global_indices(1)-1)*sim%mesh2dv%delta_eta1/sim%degree, &
         sim%mesh2dv%delta_eta1/sim%degree, &
         size(plotf2d_c1,1),size(plotf2d_c1,2), &
         plotf2d_c1, &
         "plotf2d_c1", &
         0, &
         ierr)
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
         sim%mesh2dx%delta_eta2, &
         sim%mesh2dv%eta2_min+(global_indices(2)-1)*sim%mesh2dv%delta_eta2/sim%degree, &
         sim%mesh2dv%delta_eta2/sim%degree, &
         size(plotf2d_c2,1),size(plotf2d_c2,2), &
         plotf2d_c2, &
         "plotf2d_c2", &
         0, &
         ierr)


    if((sim%test==0) .or. (sim%test==6))then
       allocate (f_x_exact(loc_sz_x1,loc_sz_v1))

       do i = 1, loc_sz_x1
          do j = 1, loc_sz_v1
!!$          global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2, &
!!$               (/1,1,1,1/) )
!!$          f_x_exact(i,j) = exp(-4*(modulo(((i-1)*sim%mesh2dx%delta_eta1 &
!!$               -(sim%mesh2dv%eta1_min+(j-1)*sim%mesh2dv%delta_eta1/sim%degree)*t),sim%mesh2dx%eta1_max-sim%mesh2dx%eta1_min)+sim%mesh2dx%eta1_min)**2)
!!$          xref=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$               sim%mesh2dx%delta_eta1
!!$          yref=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$               sim%mesh2dx%delta_eta2
!!$          x1=sim%tx%x1(xref,yref)
             x1=(i-1)*sim%mesh2dx%delta_eta1+sim%mesh2dx%eta1_min
             v1=sim%mesh2dv%eta1_min+(j-1)*sim%mesh2dv%delta_eta1/sim%degree

             sim%params(11)=t
             f_x_exact(i,j)=sim%init_func(v1,v2,x1,x2,sim%params)

          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
            sim%mesh2dx%delta_eta1, &
            sim%mesh2dv%eta1_min+(global_indices(1)-1)*sim%mesh2dv%delta_eta1/sim%degree, &
            sim%mesh2dv%delta_eta1/sim%degree, &
            f_x_exact, &
            "plotfxtransport", &
            0, &
            ierr)
    end if


    if((sim%test==4).or.(sim%test==7))then
       allocate (f_y_exact(loc_sz_x2,loc_sz_v2))
       do i = 1, loc_sz_x2
          do j = 1, loc_sz_v2
             x2=(i-1)*sim%mesh2dx%delta_eta2+sim%mesh2dx%eta2_min
             v2=sim%mesh2dv%eta2_min+(j-1)*sim%mesh2dv%delta_eta2/sim%degree

             sim%params(11)=t
             f_y_exact(i,j)=sim%init_func(v1,v2,x1,x2,sim%params)
          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
            sim%mesh2dx%delta_eta2, &
            sim%mesh2dv%eta2_min+(global_indices(2)-1)*sim%mesh2dv%delta_eta2/sim%degree, &
            sim%mesh2dv%delta_eta2/sim%degree, &
            f_y_exact, &
            "plotfytransport", &
            0, &
            ierr)
    end if

    if(sim%test==2)then
       allocate (f_vx_exact(loc_sz_x1,loc_sz_v1))
       do i = 1, loc_sz_x1
          do j = 1, loc_sz_v1
             v1=sim%mesh2dv%eta1_min+(j-1)*sim%mesh2dv%delta_eta1/sim%degree
             v2=0
             x1=0
             x2=0
             sim%params(11)=t
             f_vx_exact(i,j)=sim%init_func(v1,v2,x1,x2,sim%params)
          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
            sim%mesh2dx%delta_eta1, &
            sim%mesh2dv%eta1_min+(global_indices(1)-1)*sim%mesh2dv%delta_eta1/sim%degree, &
            sim%mesh2dv%delta_eta1/sim%degree, &
            f_vx_exact, &
            "plotfvxtransport", &
            0, &
            ierr)
    end if

    if(sim%test==3)then
       allocate (f_vy_exact(loc_sz_x2,loc_sz_v2))
       do i = 1, loc_sz_x2
          do j = 1, loc_sz_v2
             v2=sim%mesh2dv%eta2_min+(j-1)*sim%mesh2dv%delta_eta2/sim%degree
             v1=0
             x1=0
             x2=0
             sim%params(11)=t
             f_vy_exact(i,j)=sim%init_func(v1,v2,x1,x2,sim%params)
          end do
       end do
       call sll_gnuplot_rect_2d_parallel( &
            sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
            sim%mesh2dx%delta_eta2, &
            sim%mesh2dv%eta2_min+(global_indices(2)-1)*sim%mesh2dv%delta_eta2/sim%degree, &
            sim%mesh2dv%delta_eta2/sim%degree, &
            f_vy_exact, &
            "plotfvytransport", &
            0, &
            ierr)
    end if



!!$    call sll_gnuplot_rect_2d_parallel( &
!!$         sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
!!$         sim%mesh2dx%delta_eta1, &
!!$         sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
!!$         sim%mesh2dx%delta_eta2, &
!!$         plotphi2d, &
!!$         "plotphi2d", &
!!$         0, &
!!$         ierr)

!!$!when compute the normL2 of err with point center
!!$ if(sim%test==0)then
!!$    call normL2(sim,f_x_exact,plotf2d_c1,erreurL2)
!!$ end if
!!$ if(sim%test==2)then
!!$    call normL2(sim,f_vx_exact,plotf2d_c1,erreurL2) 
!!$ end if
!!$ if(sim%test==3)then
!!$    call normL2(sim,f_vy_exact,plotf2d_c2,erreurL2) 
!!$ end if
!!$ if(sim%test==4)then
!!$    call normL2(sim,f_y_exact,plotf2d_c2,erreurL2) 
!!$ end if
!!$  write(*,*) 'erreurL2 Nhung=',erreurL2
!!$  call sll_new_file_id(file_id_1,ierr)
!!$  inquire(file='log(err)', exist=exist)
!!$  if (exist) then
!!$     open(file_id_1,file='log(err)',status='old',position='append', action='write')
!!$  else
!!$     open(file_id_1, file='log(err)', status="new", action="write")
!!$  end if
!!$ if(sim%test==0)then
!!$    write(file_id_1,*) -log(sim%mesh2dx%delta_eta1),log(erreurL2)
!!$ elseif(sim%test==2)then
!!$    write(file_id_1,*) -log(sim%mesh2dv%delta_eta1/sim%degree),log(erreurL2)
!!$ elseif(sim%test==3)then
!!$    write(file_id_2,*) -log(sim%mesh2dv%delta_eta2/sim%degree), log(erreurL2)
!!$ elseif(sim%test==4)then
!!$    write(file_id_2,*) -log(sim%mesh2dx%delta_eta2), log(erreurL2)
!!$ end if
!!$  close(file_id_1)

!!$    !pour etudier le taux de convergence
!!$    !**********************************
!!$    sim%params(11)=t
!!$    call fn_L2_norm(sim,erreurL2_G)
!!$    if (sim%my_rank==0) then
!!$       call sll_new_file_id(file_id_2,ierr)
!!$       inquire(file='logerr', exist=exist)
!!$       if (exist) then
!!$          open(file_id_2,file='logerr',status='old',position='append', action='write')
!!$       else
!!$          open(file_id_2, file='logerr', status="new", action="write")
!!$       end if
!!$
!!$       if(sim%test==0)then
!!$          write(file_id_2,*) -log(sim%mesh2dx%delta_eta1), log(erreurL2_G)
!!$       elseif(sim%test==2)then
!!$          write(file_id_2,*) -log(sim%mesh2dv%delta_eta1/sim%degree), log(erreurL2_G)
!!$       elseif(sim%test==3)then
!!$          write(file_id_2,*) -log(sim%mesh2dv%delta_eta2/sim%degree), log(erreurL2_G)
!!$       elseif(sim%test==4)then
!!$          write(file_id_2,*) -log(sim%mesh2dx%delta_eta2), log(erreurL2_G)
!!$       elseif(sim%test==6)then
!!$          write(file_id_2,*) -log(sim%mesh2dv%delta_eta1), log(erreurL2_G)
!!$          !write(file_id_2,*) -log(sim%dt), log(erreurL2_G)
!!$       elseif(sim%test==7)then
!!$          write(file_id_2,*) -log(sim%mesh2dx%delta_eta2), log(erreurL2_G)
!!$       end if
!!$
!!$       close(file_id_2)
!!$    end if
!!$    !end etudier taux de convergence
!!$    !******************************

    if (sim%test .eq. 1) then
       write(*,*) 'we r using the Landau damping 1d xvx test case'
    else if (sim%test .eq. 9) then
       write(*,*) 'we r using the Landau damping 1d yvy test case'
    else if (sim%test .eq. 5) then
       write(*,*) 'we r using the Landau damping 2D  test case'
    else if (sim%test .eq. 0) then
       write(*,*) 'the x-transport test case'
    else if (sim%test .eq. 2) then
       write(*,*) 'the vx-transport test case'
    else if (sim%test .eq. 3) then
       write(*,*) 'the vy-transport test case'
    else if (sim%test .eq. 4) then
       write(*,*) 'the y-transport test case'
    else if (sim%test .eq. 4) then
       write(*,*) 'landau damping 2d test case'
    else if (sim%test .eq. 6) then
       write(*,*) 'the xvx-transport test case'
    else if (sim%test .eq. 10) then
       write(*,*) 'the two stream instability 1D test case'
    else if (sim%test .eq. 11) then
       write(*,*) 'the galaxy 1D test case'
    else if (sim%test .eq. 12) then
       write(*,*) 'the galaxy 2D test case'
    endif

    if (abs(sim%eps).gt.1.e-10_f64 ) then
       write(*,*) 'the decentered flux with esp = ', sim%eps
    else 
       write (*,*) 'the centered flux'
    endif

    if (sim%nsch == 0) then
       write(*,*) 'Euler scheme'
    else if(sim%nsch == 1) then
       write (*,*) 'R-K second order'
    else if(sim%nsch == 2) then
       write (*,*) 'R-K  order 4'
    endif




  end subroutine run_vp_cart

  subroutine velocity_mesh_connectivity(sim)

    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume) :: sim
    sll_int32 :: ierr
    sll_int32 :: i,j,k,ii,jj,ib,iel,iglob
    sll_int32 :: ic,jc,iploc,jploc,ino
    sll_real64  :: x,y,xref,yref
    sll_real64 :: det
    sll_real64 :: phi1,phi2,dphi1(2),dphi2(2),dphi1ref(2),dphi2ref(2)
    sll_real64 :: dtau(2),vnorm(2)
    sll_real64,dimension(2,2) :: jacob,invjacob
    sll_real64,dimension(:,:),allocatable :: lag,dlag,xynoe
    sll_real64,dimension(:),allocatable :: ploc
    sll_real64,dimension(:,:),allocatable :: mloc,av1loc,av2loc,bv1loc,bv2loc,&
         absbv1loc,absbv2loc
    sll_real64,dimension(:),allocatable :: gauss,weight
    sll_real64 :: void  ! only for a valid address
    sll_int32 :: ifac,isol,nsym,mp
    sll_int32 :: ll,ib1,ib2,jb1,jb2,counter

    SLL_ALLOCATE(sim%interp_pts_1D(sim%degree+1),ierr)

    ! for the moment: equally spaced points
    do i=0,sim%degree
       sim%interp_pts_1D(i+1)=real(i,f64)/sim%degree
    end do

    ! array of points coordinates

    SLL_ALLOCATE(sim%vcoords(2,sim%np_v1*sim%np_v2),ierr)

    ! connectivity
    sim%np_loc=(sim%degree+1)**2
    ! connectivity for surface elements
    SLL_ALLOCATE(sim%connec(sim%np_loc,sim%nc_v1*sim%nc_v2),ierr)
    ! connectivity for (bounadry) line elements
    SLL_ALLOCATE(sim%connecline(sim%degree+1,2*sim%nc_v1+2*sim%nc_v2),ierr)

    ! surface connectivity
    do ic=0,sim%nc_v1-1
       do jc=0,sim%nc_v2-1
          do iploc=0,sim%degree
             do jploc=0,sim%degree
                sim%connec(jploc*(sim%degree+1)+iploc+1, &
                     jc*sim%nc_v1+ic+1) = &
                     (sim%nc_v1 * sim%degree+1)* &
                     (jploc+jc*sim%degree)+ &
                     iploc+1+ic*sim%degree
             end do
          end do
       end do
    end do

    ! boundary connectivity
    counter=0

    !south boundary
    jc=0
    jploc=0
    do ic=0,sim%nc_v1-1
       counter=counter+1
       do iploc=0,sim%degree
          sim%connecline(iploc+1,counter)= &
               (sim%nc_v1 * sim%degree+1)* &
               (jploc+jc*sim%degree)+ &
               iploc+1+ic*sim%degree
       end do
    end do

    !east boundary
    ic=sim%nc_v1-1
    iploc=sim%degree
    do jc=0,sim%nc_v2-1
       counter=counter+1
       do jploc=0,sim%degree
          sim%connecline(jploc+1,counter)= &
               (sim%nc_v1 * sim%degree+1)* &
               (jploc+jc*sim%degree)+ &
               iploc+1+ic*sim%degree
       end do
    end do

    !north boundary
    jc=sim%nc_v2-1
    jploc=sim%degree
    do ic=sim%nc_v1-1,0,-1
       counter=counter+1
       do iploc=sim%degree,0,-1
          sim%connecline(sim%degree-iploc+1,counter)= &
               (sim%nc_v1 * sim%degree+1)* &
               (jploc+jc*sim%degree)+ &
               iploc+1+ic*sim%degree
       end do
    end do

    !west boundary
    ic=0
    iploc=0
    do jc=sim%nc_v2-1,0,-1
       counter=counter+1
       do jploc=sim%degree,0,-1
          sim%connecline(sim%degree-jploc+1,counter)= &
               (sim%nc_v1 * sim%degree+1)* &
               (jploc+jc*sim%degree)+ &
               iploc+1+ic*sim%degree
       end do
    end do

    sim%nel_boundary = counter

    ! geometry
    do ic=0,sim%nc_v1-1
       do jc=0,sim%nc_v2-1
          do iploc=0,sim%degree
             do jploc=0,sim%degree
                xref=sim%mesh2dv%eta1_min+ &
                     (ic+sim%interp_pts_1D(iploc+1))*sim%mesh2dv%delta_eta1
                yref=sim%mesh2dv%eta2_min+ &
                     (jc+sim%interp_pts_1D(jploc+1))*sim%mesh2dv%delta_eta2
                x=sim%tv%x1(xref,yref)
                y=sim%tv%x2(xref,yref)
                ino=sim%connec(jploc*(sim%degree+1)+iploc+1, &
                     jc*sim%nc_v1+ic+1)
                sim%vcoords(1,ino)=x
                sim%vcoords(2,ino)=y
             end do
          end do
       end do
    end do

!!$    do ino=1,sim%np_v1*sim%np_v2
!!$       write(*,*) ino,sim%vcoords(1,ino),sim%vcoords(2,ino)
!!$    end do
!!$    stop

    SLL_ALLOCATE(sim%prof(sim%np_v1*sim%np_v2),ierr)
    SLL_ALLOCATE(sim%mkld(sim%np_v1*sim%np_v2+1),ierr)

    sim%prof=0
    do k=1,sim%nc_v1*sim%nc_v2
       do ii=1,sim%np_loc
          do jj=1,sim%np_loc
             i=sim%connec(ii,k)
             j=sim%connec(jj,k)
             sim%prof(j)=max(sim%prof(j),j-i)
             sim%prof(i)=max(sim%prof(i),i-j)
          enddo
       enddo
    enddo

    sim%mkld(1)=1
    do i=1,sim%np_v1*sim%np_v2
       sim%mkld(i+1)=sim%mkld(i)+sim%prof(i)
    enddo
    !write (*,*) 'mkld = ', sim%mkld
    !write (*,*)
    sim%nsky=sim%mkld(sim%np_v1*sim%np_v2+1)-1
    !write(*,*) 'nsky = ', sim%nsky

    SLL_ALLOCATE(sim%M_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%M1_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M1_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M1_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Av1_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av1_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av1_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Av2_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av2_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av2_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv1p_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1p_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1p_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv1m_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1m_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1m_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv2p_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2p_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2p_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv2m_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2m_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2m_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%p(sim%np_v1*sim%np_v2),ierr)
    sim%p=0

    sim%M_low=0
    sim%M_sup=0
    sim%M_diag=0

    sim%Av1_low=0
    sim%Av1_sup=0
    sim%Av1_diag=0

    sim%Av2_low=0
    sim%Av2_sup=0
    sim%Av2_diag=0

    sim%Bv1p_low=0
    sim%Bv1p_sup=0
    sim%Bv1p_diag=0

    sim%Bv1m_low=0
    sim%Bv1m_sup=0
    sim%Bv1m_diag=0

    sim%Bv2p_low=0
    sim%Bv2p_sup=0
    sim%Bv2p_diag=0

    sim%Bv2m_low=0
    sim%Bv2m_sup=0
    sim%Bv2m_diag=0


    !write(*,*) 'sim%prof', sim%prof
!!$
!!$    stop

    ! init of Gauss points
    SLL_ALLOCATE(gauss(sim%degree+1),ierr)
    SLL_ALLOCATE(xynoe(2,sim%degree+1),ierr)
    SLL_ALLOCATE(weight(sim%degree+1),ierr)
    SLL_ALLOCATE(lag(sim%degree+1,sim%degree+1),ierr)
    SLL_ALLOCATE(dlag(sim%degree+1,sim%degree+1),ierr)
    SLL_ALLOCATE(ploc((sim%degree+1)**2),ierr)
    SLL_ALLOCATE(mloc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(av1loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(av2loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(bv1loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(bv2loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(absbv1loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(absbv2loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    call lag_gauss(sim%degree,gauss,weight,lag,dlag)

    ! matrix assembly
    ! surface assembly
    ! loop on the cells
    do ic=0,sim%nc_v1-1
       do jc=0,sim%nc_v2-1
          ploc=0
          mloc=0
          av1loc=0
          av2loc=0
          bv1loc=0
          bv2loc=0
          ! loop on the points
          ! loop on the Gauss points
          do iploc=0,sim%degree
             do jploc=0,sim%degree
                xref=sim%mesh2dv%eta1_min+ &
                     (ic+gauss(iploc+1))*sim%mesh2dv%delta_eta1
                yref=sim%mesh2dv%eta2_min+ &
                     (jc+gauss(jploc+1))*sim%mesh2dv%delta_eta2
                jacob=sim%tv%jacobian_matrix(xref,yref)
                invjacob=sim%tv%inverse_jacobian_matrix(xref,yref)
                det=sim%tv%jacobian(xref,yref)*sim%mesh2dv%delta_eta1*sim%mesh2dv%delta_eta2
                do ib1=1,sim%degree+1
                   do jb1=1,sim%degree+1
                      do ib2=1,sim%degree+1
                         do jb2=1,sim%degree+1
                            phi1=lag(ib1,iploc+1)*lag(jb1,jploc+1)
                            phi2=lag(ib2,iploc+1)*lag(jb2,jploc+1)
                            dphi1ref(1)=dlag(ib1,iploc+1)*lag(jb1,jploc+1)
                            dphi1ref(2)=lag(ib1,iploc+1)*dlag(jb1,jploc+1)
                            dphi2ref(1)=dlag(ib2,iploc+1)*lag(jb2,jploc+1)
                            dphi2ref(2)=lag(ib2,iploc+1)*dlag(jb2,jploc+1)
                            dphi1(1)=dphi1ref(1)*invjacob(1,1)/sim%mesh2dv%delta_eta1+ &
                                 dphi1ref(2)*invjacob(2,1)/sim%mesh2dv%delta_eta2
                            dphi1(2)=dphi1ref(1)*invjacob(1,2)/sim%mesh2dv%delta_eta1+ &
                                 dphi1ref(2)*invjacob(2,2)/sim%mesh2dv%delta_eta2
                            dphi2(1)=dphi2ref(1)*invjacob(1,1)/sim%mesh2dv%delta_eta1+ &
                                 dphi2ref(2)*invjacob(2,1)/sim%mesh2dv%delta_eta2
                            dphi2(2)=dphi2ref(1)*invjacob(1,2)/sim%mesh2dv%delta_eta1+ &
                                 dphi2ref(2)*invjacob(2,2)/sim%mesh2dv%delta_eta2
                            mloc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 mloc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                 phi1*phi2*det*weight(iploc+1)*weight(jploc+1)
                            av1loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 av1loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                 sim%tv%x1(xref,yref)*phi1*phi2*det*weight(iploc+1)*weight(jploc+1) 
                            av2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 av2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                 sim%tv%x2(xref,yref)*phi1*phi2*det*weight(iploc+1)*weight(jploc+1)
                            bv1loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 bv1loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                 dphi2(1)*phi1*det*weight(iploc+1)*weight(jploc+1) 
                            bv2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 bv2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                 dphi2(2)*phi1*det*weight(iploc+1)*weight(jploc+1)  
                            !write(*,*) 'phi1 =', phi1  
                         end do
                      end do
                      !write(*,*), 'phi1 =', phi1
                      ploc((jb1-1)*(sim%degree+1)+ib1)=ploc((jb1-1)*(sim%degree+1)+ib1)+ &
                           phi1*det*weight(iploc+1)*weight(jploc+1)
                   end do
                end do
             end do
          end do

          do ii=1,(sim%degree+1)**2
             i=sim%connec(ii,jc*sim%nc_v1+ic+1)
             sim%p(i)=sim%p(i)+ploc(ii)
             !write(*,*) ii, i
          end do

          do ii=1,(sim%degree+1)**2
             do jj=1,(sim%degree+1)**2
                i=sim%connec(ii,jc*sim%nc_v1+ic+1)
                j=sim%connec(jj,jc*sim%nc_v1+ic+1)
                if (i.eq.j) then
                   sim%M_diag(i)=sim%M_diag(i)+mloc(ii,jj)
                   sim%Av1_diag(i)=sim%Av1_diag(i)+av1loc(ii,jj)
                   sim%Av2_diag(i)=sim%Av2_diag(i)+av2loc(ii,jj)
                   sim%Bv1p_diag(i)=sim%Bv1p_diag(i)+bv1loc(ii,jj)
                   sim%Bv1m_diag(i)=sim%Bv1m_diag(i)+bv1loc(ii,jj)
                   sim%Bv2p_diag(i)=sim%Bv2p_diag(i)+bv2loc(ii,jj)
                   sim%Bv2m_diag(i)=sim%Bv2m_diag(i)+bv2loc(ii,jj)
                else if (j.gt.i) then
                   ll=sim%mkld(j+1)-j+i
                   sim%M_sup(ll)=sim%M_sup(ll)+mloc(ii,jj)
                   sim%Av1_sup(ll)=sim%Av1_sup(ll)+av1loc(ii,jj)
                   sim%Av2_sup(ll)=sim%Av2_sup(ll)+av2loc(ii,jj)
                   sim%Bv1p_sup(ll)=sim%Bv1p_sup(ll)+bv1loc(ii,jj)
                   sim%Bv2p_sup(ll)=sim%Bv2p_sup(ll)+bv2loc(ii,jj)
                   sim%Bv1m_sup(ll)=sim%Bv1m_sup(ll)+bv1loc(ii,jj)
                   sim%Bv2m_sup(ll)=sim%Bv2m_sup(ll)+bv2loc(ii,jj)
                else
                   ll=sim%mkld(i+1)-i+j
                   sim%M_low(ll)=sim%M_low(ll)+mloc(ii,jj)
                   sim%Av1_low(ll)=sim%Av1_low(ll)+av1loc(ii,jj)
                   sim%Av2_low(ll)=sim%Av2_low(ll)+av2loc(ii,jj)
                   sim%Bv1p_low(ll)=sim%Bv1p_low(ll)+bv1loc(ii,jj)
                   sim%Bv2p_low(ll)=sim%Bv2p_low(ll)+bv2loc(ii,jj)
                   sim%Bv1m_low(ll)=sim%Bv1m_low(ll)+bv1loc(ii,jj)
                   sim%Bv2m_low(ll)=sim%Bv2m_low(ll)+bv2loc(ii,jj)
                end if
             end do
          end do

          ! end loop on the cells
       end do
    end do

    ! -------------------------------------
    ! matrix assembly
    ! line assembly
    ! loop on the boundary element
    do iel=1,sim%nel_boundary
       bv1loc=0
       bv2loc=0
       absbv1loc=0
       absbv2loc=0

       ! first we get the coordinates of the nodes
       ! of element iel
       do iploc=0,sim%degree
          iglob=sim%connecline(iploc+1,iel)
          xynoe(1,iploc+1)=sim%vcoords(1,iglob)
          xynoe(2,iploc+1)=sim%vcoords(2,iglob)
          !          write(*,*) 'coords=',xynoe(1:2,iploc+1)
       end do


       ! loop on the Gauss points
       do iploc=0,sim%degree
          dtau=0
          do ib=0,sim%degree
             dtau(1)=dtau(1)+xynoe(1,ib+1)*dlag(ib+1,iploc+1)
             dtau(2)=dtau(2)+xynoe(2,ib+1)*dlag(ib+1,iploc+1)
             !write(*,*) 'dlag=',dlag(ib+1,iploc+1)
          end do
          vnorm(1)= dtau(2)
          vnorm(2)=-dtau(1)

          do ib1=1,sim%degree+1
             do ib2=1,sim%degree+1
                phi1=lag(ib1,iploc+1)
                phi2=lag(ib2,iploc+1)
                bv1loc(ib1,ib2)= bv1loc(ib1,ib2)-0.5d0*&
                     phi1*phi2*weight(iploc+1)*vnorm(1)
                bv2loc(ib1,ib2)= bv2loc(ib1,ib2)-0.5d0*&
                     phi1*phi2*weight(iploc+1)*vnorm(2)
                absbv1loc(ib1,ib2)= absbv1loc(ib1,ib2)+0.5d0*&
                     phi1*phi2*weight(iploc+1)*abs(vnorm(1))
                absbv2loc(ib1,ib2)= absbv2loc(ib1,ib2)+0.5d0*&
                     phi1*phi2*weight(iploc+1)*abs(vnorm(2))
             end do
          end do
       end do
       do ii=1,sim%degree+1
          do jj=1,sim%degree+1
             i=sim%connecline(ii,iel)
             j=sim%connecline(jj,iel)
             if (i.eq.j) then
                sim%Bv1m_diag(i)=sim%Bv1m_diag(i)+bv1loc(ii,jj)-absbv1loc(ii,jj)
                sim%Bv2m_diag(i)=sim%Bv2m_diag(i)+bv2loc(ii,jj)-absbv2loc(ii,jj)
                sim%Bv1p_diag(i)=sim%Bv1p_diag(i)+bv1loc(ii,jj)+absbv1loc(ii,jj)
                sim%Bv2p_diag(i)=sim%Bv2p_diag(i)+bv2loc(ii,jj)+absbv2loc(ii,jj)
             else if (j.gt.i) then
                ll=sim%mkld(j+1)-j+i
                sim%Bv1m_sup(ll)=sim%Bv1m_sup(ll)+bv1loc(ii,jj)-absbv1loc(ii,jj)
                sim%Bv2m_sup(ll)=sim%Bv2m_sup(ll)+bv2loc(ii,jj)-absbv2loc(ii,jj)
                sim%Bv1p_sup(ll)=sim%Bv1p_sup(ll)+bv1loc(ii,jj)+absbv1loc(ii,jj)
                sim%Bv2p_sup(ll)=sim%Bv2p_sup(ll)+bv2loc(ii,jj)+absbv2loc(ii,jj)
             else
                ll=sim%mkld(i+1)-i+j
                sim%Bv1m_low(ll)=sim%Bv1m_low(ll)+bv1loc(ii,jj)-absbv1loc(ii,jj)
                sim%Bv2m_low(ll)=sim%Bv2m_low(ll)+bv2loc(ii,jj)-absbv2loc(ii,jj)
                sim%Bv1p_low(ll)=sim%Bv1p_low(ll)+bv1loc(ii,jj)+absbv1loc(ii,jj)
                sim%Bv2p_low(ll)=sim%Bv2p_low(ll)+bv2loc(ii,jj)+absbv2loc(ii,jj)
             end if
          end do
       end do

       ! end loop on the cells
    end do

    !---------------------------------------

!!$write(*,*) 'sim%Bv1p_sup-sim%Bv1m_sup', sim%Bv1p_sup-sim%Bv1m_sup-sim%Bv1p_low+sim%Bv1m_low
!!$write(*,*) 'sim%Bv2p_sup-sim%Bv2m_sup', sim%Bv2p_sup-sim%Bv2m_sup-sim%Bv2p_low+sim%Bv2m_low
!!$write(*,*) 'sim%Bv1p_diag-sim%Bv1m_diag', sim%Bv1p_diag-sim%Bv1m_diag
!!$
!!$ stop
    !write(*,*) 'dimension of matrix : ', sim%np_v1*sim%np_v2
    !write(*,*) 'nsky = ', sim%nsky
!!$  write(*,*) 'matrix p : '
!!$   do i=1,sim%np_v2
!!$     write(*,*)  sim%p((i-1)*sim%np_v1+1:i*sim%np_v1)
!!$     write(*,*)
!!$   end do
    !stop

!!$    !write(*,*) 'M diag', sim%Av2_diag
!!$    write(*,*) 'M low', sim%M_low
!!$write(*,*) 'Bv1_diag', sim%Bv1_diag
!!$write(*,*) 'Bv1_low+Bv1_sup', sim%Bv1_low+sim%Bv1_sup
!!$stop
!!$write(*,*) 'Bv1_sup', sim%Bv1_sup



    ! LU decomposition of M
    ifac=1  ! we compute the LU decomposition
    isol=0  ! we do not solve the linear system
    nsym=1  ! we do not take into account the symetry of M
    mp=6    ! write the log on screen
    !copy the matrix M to matrix M1 to compute the energy 
    sim%M1_sup = sim%M_sup
    sim%M1_diag = sim%M_diag
    sim%M1_low = sim%M_low

    !stop
    call sol(sim%M_sup,sim%M_diag,sim%M_low,void,&
         sim%mkld,void,sim%np_v1*sim%np_v2,mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)

!!$    do i=1,sim%np_v1*sim%np_v2
!!$       write(*,*) sim%M_diag(i)
!!$    end do
!!$
!!$    write(*,*) 'fin LU'
!!$    stop



    SLL_DEALLOCATE_ARRAY(mloc,ierr)
    SLL_DEALLOCATE_ARRAY(av1loc,ierr)
    SLL_DEALLOCATE_ARRAY(av2loc,ierr)
    SLL_DEALLOCATE_ARRAY(bv1loc,ierr)
    SLL_DEALLOCATE_ARRAY(bv2loc,ierr)

!!$    deallocate(mloc)
!!$    deallocate(av1loc)
!!$    deallocate(av2loc)
!!$    deallocate(bv1loc)
!!$    deallocate(bv2loc)

  end subroutine velocity_mesh_connectivity

  subroutine delete_vp_cart( sim )
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%fn_v1v2, ierr )
    SLL_DEALLOCATE( sim%fn_star_v1v2, ierr )
    SLL_DEALLOCATE( sim%fnp1_v1v2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%dtfn_v1v2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_split, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
    call delete( sim%sequential_v1v2_layout)
    call delete( sim%phi_seq_x1_layout )
  end subroutine delete_vp_cart

  ! we put the reduction functions here for now, since we are only using
  ! simple data for the distribution function. This should go elsewhere.
  ! THIS SUBROUTINE IS JUST A PLACEHOLDER, IT IS NUMERICALLY INCORRECT.
  ! Change it later by something that uses some acceptable integrator in
  ! 1D.
  ! Design issues with this subroutine:
  ! 1. The distribution function needs to be preserved, thus this is an
  !    out-of-place operation.
  ! 2. There is probably a cleverer way to do this, but if the reduction
  !    happens in two steps a. reduction in x4 and b. reduction in x3, we
  !    need an array to store the intermediate result (after reducing in
  !    x4). This array should come as an argument.
  subroutine compute_charge_density()! mesh, numpts1, numpts2, f, partial, rho )
!!$    type(sll_logical_mesh_4d), pointer     :: mesh
!!$    sll_real64, intent(in),  dimension(:,:,:,:) :: f       ! local distr. func
!!$    sll_real64, intent(inout),  dimension(:,:,:):: partial ! intermediate res.
!!$    sll_real64, intent(inout), dimension(:,:)     :: rho     ! local rho
!!$    ! local sizes in the split directions have to be given by caller.
!!$    sll_int32, intent(in)                       :: numpts1
!!$    sll_int32, intent(in)                       :: numpts2
!!$    sll_real64                                  :: delta3
!!$    sll_real64                                  :: delta4
!!$    sll_int32                                   :: numpts3
!!$    sll_int32                                   :: numpts4
!!$    sll_int32 :: i, j, k, l
!!$    
!!$    delta4   = mesh%delta_x4
!!$    delta3   = mesh%delta_x3
!!$    partial(:,:,:) = 0.0
!!$    numpts3 = mesh%num_cells3
!!$    numpts4 = mesh%num_cells4
!!$    
!!$    ! This expects partial to be already initialized to zero!!!
!!$    do k=1,numpts3
!!$       do j=1,numpts2
!!$          do i=1,numpts1
!!$             ! This summation happens on a super-long stride... slow stuff
!!$             ! This loop should be substituted by a proper integration
!!$             ! function that we could use in the other directions as well...
!!$             do l=1,numpts4
!!$                partial(i,j,k) = partial(i,j,k) + f(i,j,k,l)*delta4
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$    ! Carry out the final reduction on x3. Note that rho is not initialized
!!$    ! to zero since it may already have the partial charge accumulation from
!!$    ! other species.
!!$    do j=1,numpts2
!!$       do i=1,numpts1
!!$          do k=1,numpts3
!!$             ! This summation happens on a very-long stride... slow stuff
!!$             ! This loop should be substituted by a proper integration
!!$             ! function that we could use in the other directions as well.
!!$             ! See above reduction function for same problem.
!!$             rho(i,j) = rho(i,j) + partial(i,j,k)*delta3
!!$          end do
!!$       end do
!!$    end do
  end subroutine compute_charge_density




!!$  subroutine plot_fields(itime, sim)
!!$    use sll_collective
!!$    use hdf5
!!$    use sll_hdf5_io_parallel
!!$    use sll_xml_io
!!$    sll_int32, intent(in) :: itime
!!$    character(len=4)      :: ctime
!!$    sll_int32             :: i_layout
!!$    character(len=1)      :: c_layout
!!$    class(sll_simulation_4d_vlasov_poisson_cart), intent(in) :: sim
!!$    type(layout_2D), pointer :: my_layout
!!$    character(len=7),  parameter :: hdf_file = "data.h5"  ! File name
!!$    sll_real64 :: tcpu1, tcpu2
!!$    sll_int32  :: my_rank
!!$    sll_int32  :: world_size
!!$    sll_int32  :: local_nx1
!!$    sll_int32  :: local_nx2
!!$    sll_int32  :: global_nx1
!!$    sll_int32  :: global_nx2
!!$    sll_int32  :: error
!!$    sll_int32  :: i
!!$    sll_int32  :: j
!!$    sll_int32  :: gi
!!$    sll_int32  :: gj
!!$    sll_int32,  dimension(2) :: global_indices
!!$    sll_real64, dimension(:,:), allocatable :: x1
!!$    sll_real64, dimension(:,:), allocatable :: x2
!!$    sll_real64 :: x1_min
!!$    sll_real64 :: x1_max
!!$    sll_real64 :: x2_min
!!$    sll_real64 :: x2_max
!!$    sll_real64 :: x3_min
!!$    sll_real64 :: x3_max
!!$    sll_real64 :: x4_min
!!$    sll_real64 :: x4_max
!!$    sll_real64 :: delta_x1
!!$    sll_real64 :: delta_x2
!!$    sll_real64 :: delta_x3
!!$    sll_real64 :: delta_x4 
!!$
!!$    integer(HID_T)                  :: hdf_file_id
!!$    sll_int32                       :: xml_file_id
!!$    integer(HSIZE_T), dimension(2)  :: array_dims 
!!$    integer(HSSIZE_T), dimension(2) :: offset 
!!$
!!$    array_dims(1) = sim%nc_x1
!!$    array_dims(2) = sim%nc_x2
!!$    world_size    = sll_get_collective_size(sll_world_collective)
!!$    my_rank       = sll_get_collective_rank(sll_world_collective)
!!$
!!$    tcpu1 = MPI_WTIME()
!!$
!!$    do i_layout = 1, 2
!!$
!!$       if (i_layout == 1) then
!!$          my_layout => sim%rho_seq_x1
!!$       else
!!$          my_layout => sim%rho_seq_x2
!!$       end if
!!$
!!$       call compute_local_sizes_2d( my_layout, local_nx1, local_nx2)        
!!$    
!!$       offset(1) =  get_layout_2D_i_min( my_layout, my_rank ) - 1
!!$       offset(2) =  get_layout_2D_j_min( my_layout, my_rank ) - 1
!!$
!!$       if (itime == 1) then
!!$
!!$          SLL_ALLOCATE(x1(local_nx1,local_nx2),error)
!!$          SLL_ALLOCATE(x2(local_nx1,local_nx2),error)
!!$       
!!$          x1_min = sim%mesh4d%x1_min
!!$          x1_max = sim%mesh4d%x1_max
!!$          x2_min = sim%mesh4d%x2_min
!!$          x2_max = sim%mesh4d%x2_max
!!$          x3_min = sim%mesh4d%x3_min
!!$          x3_max = sim%mesh4d%x3_max
!!$          x4_min = sim%mesh4d%x4_min
!!$          x4_max = sim%mesh4d%x4_max
!!$   
!!$          delta_x1 = sim%mesh4d%delta_x1
!!$          delta_x2 = sim%mesh4d%delta_x2
!!$          delta_x3 = sim%mesh4d%delta_x3
!!$          delta_x4 = sim%mesh4d%delta_x4
!!$   
!!$          do j = 1, local_nx2
!!$             do i = 1, local_nx1
!!$                global_indices =  local_to_global_2D( my_layout, (/i, j/) )
!!$                gi = global_indices(1)
!!$                gj = global_indices(2)
!!$                x1(i,j) = x1_min + (gi-1._f64)*delta_x1
!!$                x2(i,j) = x2_min + (gj-1._f64)*delta_x2
!!$             end do
!!$          end do
!!$       
!!$          call sll_hdf5_file_create("mesh_x"//c_layout//"_seq.h5",hdf_file_id,error)
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x1,"x1",error)
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,x2,"x2",error)
!!$          call sll_hdf5_file_close(hdf_file_id,error)
!!$
!!$          deallocate(x1)
!!$          deallocate(x2)
!!$
!!$       end if
!!$
!!$       call int2string(itime, ctime)
!!$       c_layout = char(i_layout+48)
!!$
!!$       call sll_hdf5_file_create("fields_x"//c_layout//"-"//ctime//".h5", &
!!$                                 hdf_file_id,error)
!!$
!!$       if (i_layout == 1) then
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%rho_x1, &
!!$                                    "rho_x"//c_layout,error)
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%phi_x1, &
!!$                                    "phi_x"//c_layout,error)
!!$       else
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%rho_x2, &
!!$                                    "rho_x"//c_layout,error)
!!$          call sll_hdf5_write_array(hdf_file_id,array_dims,offset,sim%phi_x2, &
!!$                                    "phi_x"//c_layout,error)
!!$       end if
!!$
!!$       call sll_hdf5_file_close(hdf_file_id,error)
!!$   
!!$       if (my_rank == 0) then
!!$          
!!$          !Conversion int64 -> int32
!!$          global_nx1 = transfer(array_dims(1),global_nx1)
!!$          global_nx2 = transfer(array_dims(2),global_nx2)
!!$       
!!$          call sll_xml_file_create("fields_x"//c_layout//"-"//ctime//".xmf", &
!!$                                   xml_file_id,error)
!!$          call sll_xml_grid_geometry(xml_file_id,          &
!!$                                  "mesh_x"//c_layout//"_seq.h5",global_nx1, &
!!$                                  "mesh_x"//c_layout//"_seq.h5",global_nx2, &
!!$                                  "x1", "x2" )
!!$          call sll_xml_field(xml_file_id,'rho_x'//c_layout,  &
!!$                             "fields_x"//c_layout//"-"//ctime//".h5:/rho_x"//c_layout, &
!!$                             global_nx1, global_nx2,'HDF','Node')
!!$          call sll_xml_field(xml_file_id,'phi_x'//c_layout,  &
!!$                             "fields_x"//c_layout//"-"//ctime//".h5:/phi_x"//c_layout, &
!!$                          global_nx1, global_nx2,'HDF','Node')
!!$          call sll_xml_file_close(xml_file_id,error)
!!$
!!$       end if
!!$
!!$   end do
!!$
!!$   tcpu2 = MPI_WTIME()
!!$   !if (my_rank == 0) &
!!$   !   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*world_size
!!$  
!!$  end subroutine plot_fields


  subroutine lag_gauss(degree,gauss,weight,lag,dlag)

    sll_int32 :: degree

    sll_real64,dimension(degree+1) :: gauss
    sll_real64,dimension(degree+1) :: weight
    sll_real64,dimension(degree+1,degree+1) :: lag
    sll_real64,dimension(degree+1,degree+1) :: dlag

    lag=0
    dlag=0
    gauss=0

    select case(degree)

    case(0)
       gauss(1) = 0_f64
       weight(1) = 2_f64
       lag(1,1) = 1_f64
       dlag(1,1) = 0_f64
    case(1)

       gauss(1) = -0.5773502691896257645091487805019574556476D0
       gauss(2) = 0.5773502691896257645091487805019574556476D0
       weight(1) = 0.1000000000000000000000000000000000000000D1
       weight(2) = 0.1000000000000000000000000000000000000000D1
       lag(1,1) = 0.7886751345948128822545743902509787278238D0
       lag(1,2) = 0.2113248654051871177454256097490212721762D0
       lag(2,1) = 0.2113248654051871177454256097490212721762D0
       lag(2,2) = 0.7886751345948128822545743902509787278238D0
       dlag(1,1) = -0.1D1
       dlag(1,2) = -0.1D1
       dlag(2,1) = 0.1D1
       dlag(2,2) = 0.1D1


    case(2)

       gauss(1) = -0.7745966692414833770358530799564799221666D0
       gauss(2) = 0.0D0
       gauss(3) = 0.7745966692414833770358530799564799221666D0
       weight(1) = 0.5555555555555555555555555555555555555560D0
       weight(2) = 0.8888888888888888888888888888888888888888D0
       weight(3) = 0.5555555555555555555555555555555555555560D0
       lag(1,1) = 0.6872983346207416885179265399782399610833D0
       lag(1,3) = -0.8729833462074168851792653997823996108333D-1
       lag(2,1) = 0.4000000000000000000000000000000000000000D0
       lag(2,2) = 0.1000000000000000000000000000000000000000D1
       lag(2,3) = 0.3999999999999999999999999999999999999992D0
       lag(3,1) = -0.8729833462074168851792653997823996108329D-1
       lag(3,3) = 0.6872983346207416885179265399782399610837D0
       dlag(1,1) = -0.2549193338482966754071706159912959844333D1
       dlag(1,2) = -0.1000000000000000000000000000000000000000D1
       dlag(1,3) = 0.549193338482966754071706159912959844333D0
       dlag(2,1) = 0.3098386676965933508143412319825919688666D1
       dlag(2,3) = -0.3098386676965933508143412319825919688666D1
       dlag(3,1) = -0.5491933384829667540717061599129598443332D0
       dlag(3,2) = 0.1000000000000000000000000000000000000000D1
       dlag(3,3) = 0.2549193338482966754071706159912959844333D1

    case(3)

       gauss(1) = -0.8611363115940525752239464888928095050957D0
       gauss(2) = -0.3399810435848562648026657591032446872006D0
       gauss(3) = 0.3399810435848562648026657591032446872006D0
       gauss(4) = 0.8611363115940525752239464888928095050957D0
       weight(1) = 0.3478548451374538573730639492219994072349D0
       weight(2) = 0.6521451548625461426269360507780005927648D0
       weight(3) = 0.6521451548625461426269360507780005927648D0
       weight(4) = 0.3478548451374538573730639492219994072349D0
       lag(1,1) = 0.6600056650728035304586090241410712591065D0
       lag(1,2) = 0.3373736432772532230452701535743525193096D-2
       lag(1,3) = 0.1661762313906394832923866663740817265688D-2
       lag(1,4) = 0.4924455046623182819230012194515868414855D-1
       lag(2,1) = 0.5209376877117036301110018070022304235785D0
       lag(2,2) = 0.1004885854825645727576017102247120797679D1
       lag(2,3) = -0.9921353572324654639393670446605140139555D-2
       lag(2,4) = -0.2301879032507389887619109530884603668331D0
       lag(3,1) = -0.2301879032507389887619109530884603668337D0
       lag(3,2) = -0.9921353572324654639393670446605140139348D-2
       lag(3,3) = 0.1004885854825645727576017102247120797680D1
       lag(3,4) = 0.5209376877117036301110018070022304235794D0
       lag(4,1) = 0.4924455046623182819230012194515868414852D-1
       lag(4,2) = 0.1661762313906394832923866663740817265822D-2
       lag(4,3) = 0.3373736432772532230452701535743525193161D-2
       lag(4,4) = 0.6600056650728035304586090241410712591064D0
       dlag(1,1) = -0.4315307347703724370206607038267511216231D1
       dlag(1,2) = -0.1030063844305963376996127662580753573578D1
       dlag(1,3) = 0.4998508518258898146158682533838475188246D0
       dlag(1,4) = -0.4401939455304877816988478382498684433006D0
       dlag(2,1) = 0.6070808640937936522112061914784890875756D1
       dlag(2,2) = -0.1439723163213963060623612928222340371667D1
       dlag(2,3) = -0.2969637859345816252235608844186941464072D1
       dlag(2,4) = 0.2195695238764699933604302714767248102833D1
       dlag(3,1) = -0.2195695238764699933604302714767248102829D1
       dlag(3,2) = 0.2969637859345816252235608844186941464070D1
       dlag(3,3) = 0.1439723163213963060623612928222340371668D1
       dlag(3,4) = -0.6070808640937936522112061914784890875754D1
       dlag(4,1) = 0.4401939455304877816988478382498684433001D0
       dlag(4,2) = -0.4998508518258898146158682533838475188246D0
       dlag(4,3) = 0.1030063844305963376996127662580753573578D1
       dlag(4,4) = 0.4315307347703724370206607038267511216230D1

    case(4)
       gauss(1) = -0.9061798459386639927976268782993929651254D0
       gauss(2) = -0.5384693101056830910363144207002088049674D0
       gauss(3) = 0.0D0
       gauss(4) = 0.5384693101056830910363144207002088049674D0
       gauss(5) = 0.9061798459386639927976268782993929651254D0
       weight(1) = 0.2369268850561890875142640407199173626416D0
       weight(2) = 0.4786286704993664680412915148356381929120D0
       weight(3) = 0.5688888888888888888888888888888888888888D0
       weight(4) = 0.4786286704993664680412915148356381929120D0
       weight(5) = 0.2369268850561890875142640407199173626416D0
       lag(1,1) = 0.6577278825775883905959416853145508163010D0
       lag(1,2) = 0.2206310329510026462660372235361697772277D-1
       lag(1,4) = -0.6618786099995526088371303768558919770633D-2
       lag(1,5) = -0.3237267008427455182670790754452363027991D-1
       lag(2,1) = 0.6076926946610149906004284764305419175436D0
       lag(2,2) = 0.1058797182171758427703178033099861242748D1
       lag(2,4) = 0.3922234075058382706049924394379542780080D-1
       lag(2,5) = 0.1755341081074128898504739055499048804561D0
       lag(3,1) = -0.4085820152617417192201361597504739840204D0
       lag(3,2) = -0.1134638401174469933019096956287147285025D0
       lag(3,3) = 0.1000000000000000000000000000000000000000D1
       lag(3,4) = -0.1134638401174469933019096956287147285029D0
       lag(3,5) = -0.4085820152617417192201361597504739840224D0
       lag(4,1) = 0.1755341081074128898504739055499048804554D0
       lag(4,2) = 0.3922234075058382706049924394379542780120D-1
       lag(4,4) = 0.1058797182171758427703178033099861242748D1
       lag(4,5) = 0.6076926946610149906004284764305419175430D0
       lag(5,1) = -0.3237267008427455182670790754452363027994D-1
       lag(5,2) = -0.6618786099995526088371303768558919770637D-2
       lag(5,4) = 0.2206310329510026462660372235361697772281D-1
       lag(5,5) = 0.6577278825775883905959416853145508163007D0
       dlag(1,1) = -0.6315836427348244333750187116696448638491D1
       dlag(1,2) = -0.1300170556020266432467383336339398362310D1
       dlag(1,3) = 0.3333333333333333333333333333333333333332D0
       dlag(1,4) = -0.3527563607185893186991638997617127757398D0
       dlag(1,5) = 0.4132077885315445293611787972420042209896D0
       dlag(2,1) = 0.1011127830306696558823072512093620641545D2
       dlag(2,2) = -0.2759999173503254749587529256345188128543D1
       dlag(2,3) = -0.2666666666666666666666666666666666666666D1
       dlag(2,4) = 0.2065853006980966251920623728547410404649D1
       dlag(2,5) = -0.2306021025433565979452708482027317580452D1
       dlag(3,1) = -0.5688255112620742704572067689025071136440D1
       dlag(3,2) = 0.5773266375785898115276372421470284119765D1
       dlag(3,4) = -0.5773266375785898115276372421470284119772D1
       dlag(3,5) = 0.5688255112620742704572067689025071136430D1
       dlag(4,1) = 0.2306021025433565979452708482027317580462D1
       dlag(4,2) = -0.2065853006980966251920623728547410404647D1
       dlag(4,3) = 0.2666666666666666666666666666666666666666D1
       dlag(4,4) = 0.2759999173503254749587529256345188128544D1
       dlag(4,5) = -0.1011127830306696558823072512093620641546D2
       dlag(5,1) = -0.4132077885315445293611787972420042209906D0
       dlag(5,2) = 0.3527563607185893186991638997617127757400D0
       dlag(5,3) = -0.3333333333333333333333333333333333333332D0
       dlag(5,4) = 0.1300170556020266432467383336339398362310D1
       dlag(5,5) = 0.6315836427348244333750187116696448638489D1

    case default
       write(*,*) 'degree ',degree,' not implemented !'
       stop

    end select


    ! remap to interval [0,1]
    gauss=(gauss+1)/2
    weight=weight/2


  end subroutine lag_gauss


  subroutine sourcenum(sim,Ex,Ey,w,source)

    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim   
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(in) :: w
    !sll_real64,intent(in)  :: Ex,Ey
    sll_real64  :: Ex,Ey
    !sll_real64,dimension(:,:),intent(out) :: source
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(out) :: source
    !sll_real64,dimension(:,:),allocatable :: source1,source2
    sll_real64,dimension(sim%np_v1*sim%np_v2) ::  source1,source2
    sll_int32 :: ierr,i

!!$    SLL_ALLOCATE(source1(sim%np_v1,sim%np_v2),ierr)
!!$    SLL_ALLOCATE(source2(sim%np_v1,sim%np_v2),ierr)
    sim%test=sim%params(8)
    !write(*,*) 'test =', sim%test
    if((sim%test==2) .or. (sim%test==6)) then
       Ex=1.0_f64
       Ey=0.0_f64
    endif
    if((sim%test==3).or. (sim%test==7)) then
       Ex=0.0_f64
       Ey=1.0_f64
    endif
    !correction the matrices B
!!$    do i=1,sim%np_v1*sim%np_v2
!!$       if (Ex.lt.0) then
!!$          if(sim%Bv1_diag(i).lt.0) then
!!$             Bv1_diag_corr(i)=0
!!$          end if
!!$       else
!!$          if(sim%Bv1_diag(i).gt.0) then
!!$             Bv1_diag_corr(i)=0
!!$          end if
!!$       endif
!!$       if (Ey.lt.0) then
!!$          if(sim%Bv2_diag(i).lt.0) then
!!$             Bv2_diag_corr(i)=0
!!$          end if
!!$       else
!!$          if(sim%Bv2_diag(i).gt.0) then
!!$             Bv2_diag_corr(i)=0
!!$          end if
!!$       endif
!!$    enddo
!!$
!!$
!!$
!!$    Bv1_diag_corr=0.0_f64
!!$    Bv2_diag_corr=0.0_f64

    source1=0
    source2=0

    if (Ex.lt.0) then
       call MULKU(sim%Bv1m_sup,sim%Bv1m_diag,sim%Bv1m_low, &
            sim%mkld,w,sim%np_v1*sim%np_v2,1,source1,sim%nsky)
    else
       call MULKU(sim%Bv1p_sup,sim%Bv1p_diag,sim%Bv1p_low, &
            sim%mkld,w,sim%np_v1*sim%np_v2,1,source1,sim%nsky)
    endif
!!$ if (Ex.lt.0) then
!!$    call MULKU(sim%Bv1m_low,sim%Bv1m_diag,sim%Bv1m_sup, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source1,sim%nsky)
!!$ else
!!$    call MULKU(sim%Bv1p_low,sim%Bv1p_diag,sim%Bv1p_sup, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source1,sim%nsky)
!!$ endif
!!$
!!$ if (Ey.lt.0) then
!!$    call MULKU(sim%Bv2m_low,sim%Bv2m_diag,sim%Bv2m_sup, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source2,sim%nsky)
!!$ else
!!$    call MULKU(sim%Bv2p_low,sim%Bv2p_diag,sim%Bv2p_sup, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source2,sim%nsky)
!!$ endif
    if (Ey.lt.0) then
       call MULKU(sim%Bv2m_sup,sim%Bv2m_diag,sim%Bv2m_low, &
            sim%mkld,w,sim%np_v1*sim%np_v2,1,source2,sim%nsky)
    else
       call MULKU(sim%Bv2p_sup,sim%Bv2p_diag,sim%Bv2p_low, &
            sim%mkld,w,sim%np_v1*sim%np_v2,1,source2,sim%nsky)
    endif

    source=Ex*source1+Ey*source2

    !source=-source
    !    write(*,*) 'source = ',  source
    !    stop
!!$    write(*,*) 'HELLO!!!!!!!!!!!!!!!!! COUCOU'
!!$    write(*,*) 'HELLO!!!!!!!!!!!!!!!!! COUCOU'
!!$    write(*,*) 'HELLO!!!!!!!!!!!!!!!!! COUCOU'
!!$    write(*,*) 'HELLO!!!!!!!!!!!!!!!!! COUCOU'
    !use this when we want to test the transport equation 
    if((sim%test==0).or.(sim%test==4)) then
       source=0.0_f64
    endif
      !source=0.0_f64

!!$    !can we do as following ? so we have to call only one 
!!$    !time the subroutine mulk
!!$    source=0
!!$    call MULKU(sim%Bv1_sup*Ex+Ey*sim%Bv2_sup,sim%Bv1_diag*Ex+Ey*sim%Bv2_diag, &
!!$         sim%Bv1_low*Ex+Ey*sim%Bv2_low, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source,sim%nsky)


!!$    SLL_DEALLOCATE_ARRAY(source1,ierr)
!!$    SLL_DEALLOCATE_ARRAY(source2,ierr)

  end subroutine sourcenum

  subroutine fluxnum(sim,wL,wR,vn,flux)

    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim   

    sll_real64,dimension(2),intent(in) :: vn
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(in) :: wL,wR
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(out) :: flux

    sll_real64,dimension(sim%np_v1*sim%np_v2) :: flux1,flux2
    sll_real64,dimension(sim%np_v1*sim%np_v2) :: wm
    sll_real64,dimension(sim%np_v1*sim%np_v2) :: temp
   ! sll_real64 :: time

    sim%test=sim%params(8)
    wm=(wL+wR)/2

    flux1=0
    flux2=0
   ! call set_time_mark(t0)
    call MULKU(sim%Av1_sup,sim%Av1_diag,sim%Av1_low, &
         sim%mkld,wm,sim%np_v1*sim%np_v2,1,flux1,sim%nsky)
    !time = time_elapsed_since(t0)
    !print *, 'time of one MULKU is  : ',time
    call MULKU(sim%Av2_sup,sim%Av2_diag,sim%Av2_low, &
         sim%mkld,wm,sim%np_v1*sim%np_v2,1,flux2,sim%nsky)

    flux=vn(1)*flux1+vn(2)*flux2-sim%eps/2*(wR-wL)
    if((sim%test==2).or.(sim%test==3)) then
       flux=0.0_f64
    end if

  end subroutine fluxnum

  ! time derivative of f
  subroutine dtf(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim 
  type(sll_time_mark)  :: t0,t2  
    sll_int32 :: i,j,k,beg
    sll_int32 :: count1,count2,count3,count4
    sll_int32  :: loc_sz_v1
    sll_int32  :: loc_sz_v2
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: ierr,ifac,isol,nsym,mp
    sll_real64 :: void,Ex,Ey,vn(2)
    sll_int32 :: ic,jc,icL,icR,jcL,jcR
    sll_real64 :: x1,x2,time,time2
    sll_real64,dimension(1:2,1:2) :: jac_m,inv_jac
    sll_int32,dimension(4)  :: global_indices

    sll_real64,dimension(:,:),allocatable  :: temp,flux,source

    SLL_ALLOCATE(flux(sim%np_v1,sim%np_v2),ierr)
    SLL_ALLOCATE(temp(sim%np_v1,sim%np_v2),ierr)
    SLL_ALLOCATE(source(sim%np_v1,sim%np_v2),ierr)

    ! solve  Mvx=vf
    ifac=0  ! we do not compute the LU decomposition
    isol=1  ! we do solve the linear system
    nsym=1  ! we do not take into account the symetry of M
    mp=6    ! write the log on screen

    !nsym=0    
    ! write(*,*) sim%my_rank, 'planter dans dtf1'
    call compute_local_sizes_2d( sim%phi_seq_x1_layout, loc_sz_x1, loc_sz_x2)
    ! init
    sim%dtfn_v1v2=0.0_f64
    !write(*,*) sim%my_rank, 'planter dans dtf2'
    !compute the fluxes in the x1 direction
    vn(1)=1*sim%surfx1(1,1) ! temporaire !!!!
    vn(2)=0
    ! Added the following because this is the proper layout that should be used
    ! to compute the local sizes for the following loop. ECG
    call compute_local_sizes_4d( sim%sequential_v1v2_layout, loc_sz_v1, loc_sz_v2, &
         loc_sz_x1, loc_sz_x2 )
    !write(*,*) 'entrer dans dtf1'
    !communication pour les processors
    sim%buf1 => get_buffer(sim%comm,1) 
    do k=1,loc_sz_x2
       do j=1,loc_sz_v2
          do i=1, loc_sz_v1
             sim%buf1(i+(j-1)*loc_sz_v1+(k-1)*loc_sz_v1*loc_sz_v2) &
                  = sim%fn_v1v2(i,j,1,k)
          end do
       end do
    end do
    call comm_send_real64(sim%comm,1,loc_sz_v1*loc_sz_v2*loc_sz_x2)
    sim%buf2 => get_buffer(sim%comm,2) 
    do k=1,loc_sz_x2
       do j=1,loc_sz_v2
          do i=1, loc_sz_v1
             sim%buf2(i+(j-1)*loc_sz_v1+(k-1)*loc_sz_v1*loc_sz_v2) &
                  = sim%fn_v1v2(i,j,loc_sz_x1,k)
          end do
       end do
    end do
    call comm_send_real64(sim%comm,2,loc_sz_v1*loc_sz_v2*loc_sz_x2)
    sim%buf3 => get_buffer(sim%comm,3) 
    do k=1,loc_sz_x1
       do j=1,loc_sz_v2
          do i=1, loc_sz_v1
             sim%buf3(i+(j-1)*loc_sz_v1+(k-1)*loc_sz_v1*loc_sz_v2) &
                  = sim%fn_v1v2(i,j,k,1)
          end do
       end do
    end do
    call comm_send_real64(sim%comm,3,loc_sz_v1*loc_sz_v2*loc_sz_x1)
    sim%buf4 => get_buffer(sim%comm,4) 
    do k=1,loc_sz_x1
       do j=1,loc_sz_v2
          do i=1, loc_sz_v1
             sim%buf4(i+(j-1)*loc_sz_v1+(k-1)*loc_sz_v1*loc_sz_v2) &
                  = sim%fn_v1v2(i,j,k,loc_sz_x2)
          end do
       end do
    end do
    call comm_send_real64(sim%comm,4,loc_sz_v1*loc_sz_v2*loc_sz_x1)

    do jc=1,loc_sz_x2
       do ic=1,loc_sz_x1-1
          icL=ic
          icR=ic+1
          call fluxnum(sim,sim%fn_v1v2(:,:,icL,jc), &
               sim%fn_v1v2(:,:,icR,jc),vn,flux)
          sim%dtfn_v1v2(:,:,icL,jc)=sim%dtfn_v1v2(:,:,icL,jc)-flux
          sim%dtfn_v1v2(:,:,icR,jc)=sim%dtfn_v1v2(:,:,icR,jc)+flux
       end do
    end do
    call comm_receive_real64(sim%comm,1,count1)
    call comm_receive_real64(sim%comm,2,count2)
    if ((count1.ne.loc_sz_v1*loc_sz_v2*loc_sz_x2).or. &
         (count2.ne.loc_sz_v1*loc_sz_v2*loc_sz_x2)) then
       write(*,*) 'problem avec send de mpi'
       stop
    endif

    !write(*,*) 'entrer dans dtf2'
    sim%buf1 => get_buffer(sim%comm,1)
    sim%buf2 => get_buffer(sim%comm,2)
    do jc=1,loc_sz_x2
       beg=(jc-1)*loc_sz_v1*loc_sz_v2
       call  fluxnum(sim,sim%buf1(beg+1:beg+loc_sz_v1*loc_sz_v2), &
            sim%fn_v1v2(:,:,1,jc),vn,flux)
       sim%dtfn_v1v2(:,:,1,jc)=sim%dtfn_v1v2(:,:,1,jc)+flux
       call  fluxnum(sim,sim%fn_v1v2(:,:,loc_sz_x1,jc), &
            sim%buf2(beg+1:beg+loc_sz_v1*loc_sz_v2),vn,flux)
       sim%dtfn_v1v2(:,:,loc_sz_x1,jc)=sim%dtfn_v1v2(:,:,loc_sz_x1,jc)-flux
    end do
    !write(*,*) sim%my_rank, 'planter dans dtf4'

    !write(*,*) 'sim%dtfn_v1v2',maxval(abs(sim%dtfn_v1v2(1,1,:,:)))

    !compute the fluxes in the x2 direction
    !write(*,*) 'ENTRER DANS dtf'
    vn(1)=0 ! temporaire !!!!
    vn(2)=1*sim%surfx2(1,1)
    !write(*,*) sim%my_rank, 'planter dans dtf1'
    do ic=1,loc_sz_x1
       do jc=1,loc_sz_x2-1
          jcL=jc
          jcR=jc+1
          call fluxnum(sim,sim%fn_v1v2(:,:,ic,jcL), &
               sim%fn_v1v2(:,:,ic,jcR),vn,flux)
          sim%dtfn_v1v2(:,:,ic,jcL)=sim%dtfn_v1v2(:,:,ic,jcL)-flux
          sim%dtfn_v1v2(:,:,ic,jcR)=sim%dtfn_v1v2(:,:,ic,jcR)+flux
       end do
    end do
    call comm_receive_real64(sim%comm,3,count3)
    call comm_receive_real64(sim%comm,4,count4)
    if ((count3.ne.loc_sz_v1*loc_sz_v2*loc_sz_x1).or. &
         (count4.ne.loc_sz_v1*loc_sz_v2*loc_sz_x1)) then
       write(*,*) 'problem avec send de mpi'
       stop
    endif
    sim%buf3 => get_buffer(sim%comm,3)
    sim%buf4 => get_buffer(sim%comm,4)
    do ic=1,loc_sz_x1
       beg=(ic-1)*loc_sz_v1*loc_sz_v2
       call  fluxnum(sim,sim%buf3(beg+1:beg+loc_sz_v1*loc_sz_v2), &
            sim%fn_v1v2(:,:,ic,1),vn,flux)
       sim%dtfn_v1v2(:,:,ic,1)=sim%dtfn_v1v2(:,:,ic,1)+flux
       call  fluxnum(sim,sim%fn_v1v2(:,:,ic,loc_sz_x2), &
            sim%buf4(beg+1:beg+loc_sz_v1*loc_sz_v2),vn,flux)
       sim%dtfn_v1v2(:,:,ic,loc_sz_x2)=sim%dtfn_v1v2(:,:,ic,loc_sz_x2)-flux
    end do
    !write(*,*) 'entrer dans dtf3'
    ! pour phi
    sim%buf1 => get_buffer(sim%comm,1) 
    do k=1,loc_sz_x2
       sim%buf1(k)=sim%phi_split(1,k)
    end do
    call comm_send_real64(sim%comm,1,loc_sz_x2)

    sim%buf2 => get_buffer(sim%comm,2) 
    do k=1,loc_sz_x2
       sim%buf2(k)=sim%phi_split(loc_sz_x1,k)
    end do
    call comm_send_real64(sim%comm,2,loc_sz_x2)

    sim%buf3 => get_buffer(sim%comm,3) 
    do k=1,loc_sz_x1
       sim%buf3(k)=sim%phi_split(k,1)
    end do
    call comm_send_real64(sim%comm,3,loc_sz_x1)

    sim%buf4 => get_buffer(sim%comm,4) 
    do k=1,loc_sz_x1
       sim%buf4(k)=sim%phi_split(k,loc_sz_x2)
    end do
    call comm_send_real64(sim%comm,4,loc_sz_x1)
    !write(*,*) 'entrer dans dtf4 pour phi'
    do ic=2,loc_sz_x1-1
       do jc=2,loc_sz_x2-1
          icL=ic-1
          icR=ic+1
          global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
               (/1,1,1,1/) )
          x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
               sim%mesh2dx%delta_eta1
          x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
               sim%mesh2dx%delta_eta2
          jac_m=sim%tx%jacobian_matrix(x1,x2)
          inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
          Ex=-(sim%phi_split(icR,jc)-sim%phi_split(icL,jc))/2/ &
               sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(ic,jc+1)- &
               sim%phi_split(ic,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
          Ey=-(sim%phi_split(ic,jc+1)-sim%phi_split(ic,jc-1))/2/ &
               sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,jc)- &
               sim%phi_split(icL,jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
  !call set_time_mark(t0)
          call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,ic,jc), &
               source)
 ! time = time_elapsed_since(t0)
  !print *, 'time of one sourcnum is  : ',time


          temp=sim%dtfn_v1v2(:,:,ic,jc)-sim%volume(1,1)*source
  !call set_time_mark(t2)
          call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
               sim%mkld, &
               sim%dtfn_v1v2(:,:,ic,jc),&
               sim%np_v1*sim%np_v2, &
               mp,ifac,isol,nsym,void,ierr,&
               sim%nsky)
  !time2 = time_elapsed_since(t2)
  !print *, 'time of one subr sol is : ',time2

       end do
    end do
    !write(*,*) 'entrer dans dtf4 pour phi 2'
    call comm_receive_real64(sim%comm,3,count3)
   !write(*,*) 'recevoir pour porte3'
    call comm_receive_real64(sim%comm,4,count4)
   !write(*,*) 'after receive1'
    if ((count3.ne.loc_sz_x1).or.(count4.ne.loc_sz_x1)) then
       write(*,*) 'problem avec send de mpi'
       stop
    endif
    !write(*,*) 'after receive2'
    sim%buf3 => get_buffer(sim%comm,3)
    sim%buf4 => get_buffer(sim%comm,4)
    !write(*,*) 'entrer dans dtf4 pour phi 21'
    !pour jc=1 et Loc_sz_x2
    do ic=2,loc_sz_x1-1
       icL=ic-1
       icR=ic+1
       global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
            (/1,1,1,1/) )
       x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
            sim%mesh2dx%delta_eta1
       x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
            sim%mesh2dx%delta_eta2
       jac_m=sim%tx%jacobian_matrix(x1,x2)
       inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)

       Ex=-(sim%phi_split(icR,1)-sim%phi_split(icL,1))/2/ &
            sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(ic,2)- &
            sim%buf3(ic))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
       Ey=-(sim%phi_split(ic,2)-sim%buf3(ic))/2/ &
            sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,1)- &
            sim%phi_split(icL,1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
       call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,ic,1), &
            source)
       temp=sim%dtfn_v1v2(:,:,ic,1)-sim%volume(1,1)*source
       call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
            sim%mkld, &
            sim%dtfn_v1v2(:,:,ic,1),&
            sim%np_v1*sim%np_v2, &
            mp,ifac,isol,nsym,void,ierr,&
            sim%nsky)
       global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
            (/1,1,1,1/) )
       x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
            sim%mesh2dx%delta_eta1
       x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
            sim%mesh2dx%delta_eta2
       jac_m=sim%tx%jacobian_matrix(x1,x2)
       inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)

       Ex=-(sim%phi_split(icR,loc_sz_x2)-sim%phi_split(icL,loc_sz_x2))/2/ &
            sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf4(ic)- &
            sim%phi_split(ic,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
       Ey=-(sim%buf4(ic)-sim%phi_split(ic,loc_sz_x2-1))/2/ &
            sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(icR,loc_sz_x2)- &
            sim%phi_split(icL,loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)
       call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,ic,loc_sz_x2), &
            source)
       temp=sim%dtfn_v1v2(:,:,ic,loc_sz_x2)-sim%volume(1,1)*source
       call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
            sim%mkld, &
            sim%dtfn_v1v2(:,:,ic,loc_sz_x2),&
            sim%np_v1*sim%np_v2, &
            mp,ifac,isol,nsym,void,ierr,&
            sim%nsky)
    end do
    !write(*,*) 'entrer dans dtf4 pour phi here apres 2'

    call comm_receive_real64(sim%comm,1,count1)
    call comm_receive_real64(sim%comm,2,count2)
    if ((count1.ne.loc_sz_x2).or.(count2.ne.loc_sz_x2)) then
       write(*,*) 'problem avec send de mpi'
       stop
    endif
    sim%buf1 => get_buffer(sim%comm,1)
    sim%buf2 => get_buffer(sim%comm,2)
    !write(*,*) 'entrer dans dtf4 pour phi 3'
    !pour ic=1
    do jc=2,loc_sz_x2-1
       global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
            (/1,1,1,1/) )
       x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
            sim%mesh2dx%delta_eta1
       x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
            sim%mesh2dx%delta_eta2
       jac_m=sim%tx%jacobian_matrix(x1,x2)
       inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
       Ex=-(sim%phi_split(2,jc)-sim%buf1(jc))/2/ &
            sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(1,jc+1)- &
            sim%phi_split(1,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
       Ey=-(sim%phi_split(1,jc+1)-sim%phi_split(1,jc-1))/2/ &
            sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,jc)- &
            sim%buf1(jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

       call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,1,jc), &
            source)


       temp=sim%dtfn_v1v2(:,:,1,jc)-sim%volume(1,1)*source

       call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
            sim%mkld, &
            sim%dtfn_v1v2(:,:,1,jc),&
            sim%np_v1*sim%np_v2, &
            mp,ifac,isol,nsym,void,ierr,&
            sim%nsky)

    end do
    !pour ic=loc_sz_x1
    do jc=2,loc_sz_x2-1
       global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
            (/1,1,1,1/) )
       x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
            sim%mesh2dx%delta_eta1
       x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
            sim%mesh2dx%delta_eta2
       jac_m=sim%tx%jacobian_matrix(x1,x2)
       inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
       Ex=-(sim%buf2(jc)-sim%phi_split(loc_sz_x1-1,jc))/2/ &
            sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(loc_sz_x1,jc+1)- &
            sim%phi_split(loc_sz_x1,jc-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
       Ey=-(sim%phi_split(loc_sz_x1,jc+1)-sim%phi_split(loc_sz_x1,jc-1))/2/ &
            sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(jc)- &
            sim%phi_split(loc_sz_x1-1,jc))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

       call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,loc_sz_x1,jc), &
            source)


       temp=sim%dtfn_v1v2(:,:,loc_sz_x1,jc)-sim%volume(1,1)*source

       call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
            sim%mkld, &
            sim%dtfn_v1v2(:,:,loc_sz_x1,jc),&
            sim%np_v1*sim%np_v2, &
            mp,ifac,isol,nsym,void,ierr,&
            sim%nsky)

    end do
    !write(*,*) 'entrer dans dtf4 pour phi 4'
    !pour ic=1 et jc=1
!!$    global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$         (/1,1,1,1/) )
!!$    x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$         sim%mesh2dx%delta_eta1
!!$    x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$         sim%mesh2dx%delta_eta2
!!$    jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$    inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
    Ex=-(sim%phi_split(2,1)-sim%buf1(1))/2/ &
         sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(1,2)- &
         sim%buf3(1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
    Ey=-(sim%phi_split(1,2)-sim%buf3(1))/2/ &
         sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,1)- &
         sim%buf1(1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

    call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,1,1), &
         source)


    temp=sim%dtfn_v1v2(:,:,1,1)-sim%volume(1,1)*source

    call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
         sim%mkld, &
         sim%dtfn_v1v2(:,:,1,1),&
         sim%np_v1*sim%np_v2, &
         mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)
    !pour ic=1 et jc=loc_sz_x2
!!$ global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$               (/1,1,1,1/) )
!!$          x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$               sim%mesh2dx%delta_eta1
!!$          x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$               sim%mesh2dx%delta_eta2
!!$          jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$          inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
    Ex=-(sim%phi_split(2,loc_sz_x2)-sim%buf1(loc_sz_x2))/2/ &
         sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf2(1)- &
         sim%phi_split(1,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
    Ey=-(sim%buf4(1)-sim%phi_split(1,loc_sz_x2-1))/2/ &
         sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%phi_split(2,loc_sz_x2)- &
         sim%buf1(loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

    !write(*,*) 'Ey =', Ey

    call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,1,loc_sz_x2), &
         source)


    temp=sim%dtfn_v1v2(:,:,1,loc_sz_x2)-sim%volume(1,1)*source

    call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
         sim%mkld, &
         sim%dtfn_v1v2(:,:,1,loc_sz_x2),&
         sim%np_v1*sim%np_v2, &
         mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)
    !pour ic=loc_sz_x1 et jc=1
!!$ global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$               (/1,1,1,1/) )
!!$          x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$               sim%mesh2dx%delta_eta1
!!$          x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$               sim%mesh2dx%delta_eta2
!!$          jac_m=sim%tx%jacobian_matrix(x1,x2)
          inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
    Ex=-(sim%buf2(1)-sim%phi_split(loc_sz_x1-1,1))/2/ &
         sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%phi_split(loc_sz_x1,2)- &
         sim%buf3(loc_sz_x1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
    Ey=-(sim%phi_split(loc_sz_x1,2)-sim%buf3(loc_sz_x1))/2/ &
         sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(1)- &
         sim%phi_split(loc_sz_x1-1,1))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

    call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,loc_sz_x1,1), &
         source)


    temp=sim%dtfn_v1v2(:,:,loc_sz_x1,1)-sim%volume(1,1)*source

    call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
         sim%mkld, &
         sim%dtfn_v1v2(:,:,loc_sz_x1,1),&
         sim%np_v1*sim%np_v2, &
         mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)
    !pour ic=loc_sz_x1 et jc=loc_sz_x2
!!$ global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, &
!!$               (/1,1,1,1/) )
!!$          x1=  sim%mesh2dx%eta1_min+real(global_indices(3)-1,f64)* &
!!$               sim%mesh2dx%delta_eta1
!!$          x2=  sim%mesh2dx%eta2_min+real(global_indices(4)-1,f64)* &
!!$               sim%mesh2dx%delta_eta2
!!$          jac_m=sim%tx%jacobian_matrix(x1,x2)
!!$          inv_jac=sim%tx%inverse_jacobian_matrix(x1,x2)
    Ex=-(sim%buf2(loc_sz_x2)-sim%phi_split(loc_sz_x1-1,loc_sz_x2))/2/ &
         sim%mesh2dx%delta_eta1*inv_jac(1,1)-(sim%buf4(loc_sz_x1)- &
         sim%phi_split(loc_sz_x1,loc_sz_x2-1))/2/sim%mesh2dx%delta_eta2*inv_jac(2,1)
    Ey=-(sim%buf4(loc_sz_x1)-sim%phi_split(loc_sz_x1,loc_sz_x2-1))/2/ &
         sim%mesh2dx%delta_eta2*inv_jac(2,2)-(sim%buf2(loc_sz_x2)- &
         sim%phi_split(loc_sz_x1-1,loc_sz_x2))/2/sim%mesh2dx%delta_eta1*inv_jac(1,2)

    call sourcenum(sim,Ex,Ey,sim%fn_v1v2(:,:,loc_sz_x1,loc_sz_x2), &
         source)


    temp=sim%dtfn_v1v2(:,:,loc_sz_x1,loc_sz_x2)-sim%volume(1,1)*source

    call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
         sim%mkld, &
         sim%dtfn_v1v2(:,:,loc_sz_x1,loc_sz_x2),&
         sim%np_v1*sim%np_v2, &
         mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)





    SLL_DEALLOCATE_ARRAY(flux,ierr)
    SLL_DEALLOCATE_ARRAY(temp,ierr)



  end subroutine dtf

  subroutine euler(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
  type(sll_time_mark)  :: t0 
    sll_int32 :: jc
    sll_int32  :: loc_sz_x2,loc_sz_x1
    sll_real64 :: time
    call compute_local_sizes_2d( sim%phi_seq_x1_layout, loc_sz_x1, loc_sz_x2)
    ! mpi communications
    !write(*,*) 'sim%dtfn_v1v2= ', sim%dtfn_v1v2
!!$    do jc=0,loc_sz_x2+1
!!$       write(*,*) 'avant',jc, sim%fn_v1v2(1,1,:,jc)
!!$    end do
!!$    do jc=0,loc_sz_x2+1
!!$       write(*,*) 'apres',jc, sim%fn_v1v2(1,1,:,jc)
!!$    end do
    !write(*,*) 'ici1'
  !call set_time_mark(t0)
    call dtf(sim)
  !time = time_elapsed_since(t0)
  !print *, 'time of simulation est  : ',time
    !write(*,*) 'ici2'
!!$    do jc=0,loc_sz_x2+1
!!$       write(*,*) 'apres',jc, sim%dtfn_v1v2(1,1,:,jc)
!!$    end do
    !stop
    !write(*,*) 'fist call dtf: sim%dtfn_v1v2= ', sim%dtfn_v1v2
!!$    sim%fn_v1v2 = sim%fn_v1v2 &
!!$         + sim%dt*sim%dtfn_v1v2/sim%volume(1,1)
    sim%fn_v1v2 = sim%fn_v1v2 &
         +sim%dtfn_v1v2/sim%volume(1,1)*sim%dt
    !write(*,*) 'volume = ',sim%volume(1,1)
    !stop
!!$    do jc=0,loc_sz_x2+1
!!$       write(*,*) 'apres',jc, sim%fn_v1v2(1,1,:,jc)
!!$    end do
    !write(*,*) 'second call dtf: sim%dtfn_v1v2= ', sim%dtfn_v1v2
  end subroutine euler

  subroutine RK2(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
    call dtf(sim)
    sim%fnp1_v1v2 = sim%fn_v1v2
    sim%fn_v1v2 = sim%fn_v1v2 &
         + sim%dt/2/sim%volume(1,1)*sim%dtfn_v1v2
    call dtf(sim)
    sim%fn_v1v2 = sim%fnp1_v1v2 &
         + sim%dt/sim%volume(1,1)*sim%dtfn_v1v2
  end subroutine RK2

  subroutine RK4(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
    call dtf(sim)
    sim%fn_star_v1v2 = sim%fn_v1v2
    sim%fn_v1v2 = sim%fn_star_v1v2 &
         + sim%dt/2/sim%volume(1,1)*sim%dtfn_v1v2
    sim%fnp1_v1v2 = sim%fn_star_v1v2+sim%dt/6/sim%volume(1,1) &
         *sim%dtfn_v1v2
    call dtf(sim)
    sim%fn_v1v2 = sim%fn_star_v1v2 &
         + sim%dt/2/sim%volume(1,1)*sim%dtfn_v1v2
    sim%fnp1_v1v2 = sim%fnp1_v1v2+sim%dt/3/sim%volume(1,1) &
         *sim%dtfn_v1v2
    call dtf(sim)
    sim%fn_v1v2 = sim%fn_star_v1v2 &
         + sim%dt/sim%volume(1,1)*sim%dtfn_v1v2
    sim%fnp1_v1v2 = sim%fnp1_v1v2+sim%dt/3/sim%volume(1,1) &
         *sim%dtfn_v1v2
    call dtf(sim)
    sim%fn_v1v2 = sim%fnp1_v1v2 &
         + sim%dt/6/sim%volume(1,1)*sim%dtfn_v1v2
  end subroutine RK4

!!$subroutine mpi_comm(sim)
!!$ class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim   
!!$ sll_int32 :: ierr
!!$ sll_int32  :: loc_sz_v1
!!$ sll_int32  :: loc_sz_v2
!!$ sll_int32  :: loc_sz_x1
!!$ sll_int32  :: loc_sz_x2
!!$ sll_int32  :: ranktop
!!$ sll_int32  :: rankbottom
!!$ sll_int32  :: rankright
!!$ sll_int32  :: rankleft
!!$ sll_int32  :: datasize,datasizephi
!!$
!!$ call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
!!$      loc_sz_v1, &
!!$      loc_sz_v2, &
!!$      loc_sz_x1, &
!!$      loc_sz_x2 )
!!$ ranktop=mod(sim%my_rank+1,sim%world_size)
!!$ rankbottom=sim%my_rank-1
!!$ ranktop=mod(sim%my_rank+1,sim%world_size)
!!$ rankbottom=sim%my_rank-1
!!$ if (rankbottom.lt.0) rankbottom=sim%world_size-1
!!$ datasize=loc_sz_v1*loc_sz_v2*loc_sz_x1
!!$ datasizephi=loc_sz_x1
!!$
!!$ Call mpi_SENDRECV(sim%fn_v1v2(1,1,1,loc_sz_x2),datasize, &
!!$      MPI_DOUBLE_PRECISION,ranktop,sim%my_rank,              &
!!$      sim%fn_v1v2(1,1,1,0),datasize,            &
!!$      MPI_DOUBLE_PRECISION,rankbottom,rankbottom,              &
!!$      MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)   
!!$
!!$ Call mpi_SENDRECV(sim%phi_split(1,loc_sz_x2),datasizephi, &
!!$      MPI_DOUBLE_PRECISION,ranktop,sim%my_rank,              &
!!$      sim%phi_split(1,0),datasizephi,            &
!!$      MPI_DOUBLE_PRECISION,rankbottom,rankbottom,              &
!!$      MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)   
!!$
!!$ ! bottom communications
!!$ Call mpi_SENDRECV(sim%fn_v1v2(1,1,1,1),datasize, &
!!$      MPI_DOUBLE_PRECISION,rankbottom,sim%my_rank,              &
!!$      sim%fn_v1v2(1,1,1,loc_sz_x2+1),datasize,            &
!!$      MPI_DOUBLE_PRECISION,ranktop,ranktop,              &
!!$      MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)       
!!$
!!$ Call mpi_SENDRECV(sim%phi_split(1,1),datasizephi, &
!!$      MPI_DOUBLE_PRECISION,rankbottom,sim%my_rank,              &
!!$      sim%phi_split(1,loc_sz_x2+1),datasizephi,            &
!!$      MPI_DOUBLE_PRECISION,ranktop,ranktop,              &
!!$      MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)       
!!$
!!$
!!$
!!$end subroutine mpi_comm

  !compute the L2 norm for the test cases in one dimension
  subroutine normL2(sim,w1,w2,res)
    !for instance, but after i want to code this sub with any dimension of vector
    !n in fact here i compute the norm L2 with the center points
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim 
    sll_real64,dimension(sim%nc_x1,sim%np_v1), intent(in) :: w1,w2
    sll_real64, intent(out) :: res
    sll_int32 :: i,j
    res=0.0_f64
    do i=1,sim%nc_x1
       do j=1,sim%np_v1
          res=res+(w1(i,j)-w2(i,j))*(w1(i,j)-w2(i,j))
       end do
    end do
    res=sqrt(res*sim%mesh2dx%delta_eta1*sim%mesh2dv%delta_eta1/sim%degree)
!!$
    !check dv/degree
  end subroutine normL2

  subroutine exact_x_transport()

  end subroutine exact_x_transport


  subroutine fn_L2_norm(sim,norml2_glob)

    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(in) :: sim

    sll_real64,dimension(:,:),allocatable :: lag,dlag
    sll_real64,dimension(:),allocatable :: gauss,weight
    sll_int32  :: ierr,loc_sz_x1,loc_sz_x2,loc_sz_v1,loc_sz_v2
    sll_int32,dimension(4)  :: global_indices
    sll_real64 :: x1,x2,v1,v2,f,norml2,det,vol_loc,vol_glob
    sll_int32  :: ib1,ib2,icv1,icv2,igv1,igv2,icx1,icx2,iv1,iv2
    sll_real64,intent(out) :: norml2_glob

    SLL_ALLOCATE(weight(sim%degree+1),ierr)
    SLL_ALLOCATE(gauss(sim%degree+1),ierr)
    SLL_ALLOCATE(lag(sim%degree+1,sim%degree+1),ierr)
    SLL_ALLOCATE(dlag(sim%degree+1,sim%degree+1),ierr)
    call lag_gauss(sim%degree,gauss,weight,lag,dlag)

    call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
         loc_sz_v1, &
         loc_sz_v2, &
         loc_sz_x1, &
         loc_sz_x2 )



    global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, (/1,1,1,1/) )

    ! loop on cells and elems
    normL2=0
    vol_loc=0
    do icx1=1,loc_sz_x1
       do icx2=1,loc_sz_x2
          do icv1=1,sim%nc_v1
             do icv2=1,sim%nc_v2
                ! loop on velocity gauss points
                do igv1=1,sim%degree+1
                   do igv2=1,sim%degree+1
                      f=0
                      x1=sim%mesh2dx%eta1_min+ &
                           (icx1-1)*sim%mesh2dx%delta_eta1
                      x2=sim%mesh2dx%eta2_min+ &
                           (icx2-1)*sim%mesh2dx%delta_eta2
                      ! we suppose a uniform mesh in x1,x2
                      v1=sim%mesh2dv%eta1_min+ &
                           (icv1-1+gauss(igv1))*sim%mesh2dv%delta_eta1
                      v2=sim%mesh2dv%eta2_min+ &
                           (icv2-1+gauss(igv2))*sim%mesh2dv%delta_eta2
!!$                    jacob=sim%tv%jacobian_matrix(v1,v2)
!!$                    invjacob=sim%tv%inverse_jacobian_matrix(v1,v2)
                      det=sim%tv%jacobian(v1,v2)* &
                           sim%mesh2dv%delta_eta1*sim%mesh2dv%delta_eta2 &
                           *sim%mesh2dx%delta_eta1*sim%mesh2dx%delta_eta2
                      do ib1=1,sim%degree+1
                         do ib2=1,sim%degree+1
                            iv1=(icv1-1)*sim%degree+ib1
                            iv2=(icv2-1)*sim%degree+ib2
                            f=f+sim%fn_v1v2(iv1,iv2,icx1,icx2)*lag(ib1,igv1)*lag(ib2,igv2)
                         end do
                      end do
                      f=f-sim%init_func(v1,v2,x1,x2,sim%params)
                      vol_loc=vol_loc+weight(igv1)*weight(igv2)*det
                      normL2=normL2+f*f*weight(igv1)*weight(igv2)*det
                   end do
                end do
             end do
          end do
       end do
    end do

    Call MPI_ALLREDUCE(normL2,normL2_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
         MPI_COMM_WORLD,ierr)
    Call MPI_ALLREDUCE(vol_loc,vol_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
         MPI_COMM_WORLD,ierr)
    normL2_glob=sqrt(normL2_glob/vol_glob)


  end subroutine fn_L2_norm

!!$  !calcul the first moment 
!!$  subroutine first_moment(sim,moment)
!!$
!!$    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(in) :: sim
!!$
!!$    sll_real64,dimension(:,:),allocatable :: lag,dlag
!!$    sll_real64,dimension(:),allocatable :: gauss,weight
!!$    sll_int32  :: ierr,loc_sz_x1,loc_sz_x2,loc_sz_v1,loc_sz_v2
!!$    sll_int32,dimension(4)  :: global_indices
!!$    sll_real64 :: x1,x2,v1,v2,f,mm,det,vol_loc,vol_glob
!!$    sll_int32  :: ib1,ib2,icv1,icv2,igv1,igv2,icx1,icx2,iv1,iv2
!!$    sll_real64,intent(out) :: moment
!!$
!!$    SLL_ALLOCATE(weight(sim%degree+1),ierr)
!!$    SLL_ALLOCATE(gauss(sim%degree+1),ierr)
!!$    SLL_ALLOCATE(lag(sim%degree+1,sim%degree+1),ierr)
!!$    SLL_ALLOCATE(dlag(sim%degree+1,sim%degree+1),ierr)
!!$    call lag_gauss(sim%degree,gauss,weight,lag,dlag)
!!$
!!$    call compute_local_sizes_4d( sim%sequential_v1v2_layout, &
!!$         loc_sz_v1, &
!!$         loc_sz_v2, &
!!$         loc_sz_x1, &
!!$         loc_sz_x2 )
!!$
!!$
!!$
!!$    global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2_layout, (/1,1,1,1/) )
!!$
!!$    ! loop on cells and elems
!!$    mm=0
!!$    vol_loc=0
!!$    do icx1=1,loc_sz_x1
!!$       do icx2=1,loc_sz_x2
!!$          do icv1=1,sim%nc_v1
!!$             do icv2=1,sim%nc_v2
!!$                ! loop on velocity gauss points
!!$                do igv1=1,sim%degree+1
!!$                   do igv2=1,sim%degree+1
!!$                      f=0
!!$                      x1=sim%mesh2dx%eta1_min+ &
!!$                           (icx1-1)*sim%mesh2dx%delta_eta1
!!$                      x2=sim%mesh2dx%eta2_min+ &
!!$                           (icx2-1)*sim%mesh2dx%delta_eta2
!!$                      ! we suppose a uniform mesh in x1,x2
!!$                      v1=sim%mesh2dv%eta1_min+ &
!!$                           (icv1-1+gauss(igv1))*sim%mesh2dv%delta_eta1
!!$                      v2=sim%mesh2dv%eta2_min+ &
!!$                           (icv2-1+gauss(igv2))*sim%mesh2dv%delta_eta2
!!$                      det=sim%tv%jacobian(v1,v2)* &
!!$                           sim%mesh2dv%delta_eta1*sim%mesh2dv%delta_eta2 &
!!$                           *sim%mesh2dx%delta_eta1*sim%mesh2dx%delta_eta2
!!$                      do ib1=1,sim%degree+1
!!$                         do ib2=1,sim%degree+1
!!$                            iv1=(icv1-1)*sim%degree+ib1
!!$                            iv2=(icv2-1)*sim%degree+ib2
!!$                            f=f+sim%fn_v1v2(iv1,iv2,icx1,icx2)*lag(ib1,igv1)*lag(ib2,igv2)
!!$                         end do
!!$                      end do
!!$                      f=f-sim%init_func(v1,v2,x1,x2,sim%params)
!!$                      vol_loc=vol_loc+weight(igv1)*weight(igv2)*det
!!$                      normL2=normL2+f*f*weight(igv1)*weight(igv2)*det
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    Call MPI_ALLREDUCE(normL2,normL2_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!!$         MPI_COMM_WORLD,ierr)
!!$    Call MPI_ALLREDUCE(vol_loc,vol_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!!$         MPI_COMM_WORLD,ierr)
!!$    normL2_glob=sqrt(normL2_glob/vol_glob)
!!$
!!$
!!$  end subroutine first_moment
  !Solution Exact of the test case Landau damping 1D x_vx
  !compute with the formula of Eric Sonnendrucker
  subroutine solexact(sim,t,E2norm_ex)

  class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(in) :: sim
sll_int32 :: i,nx,ierr
sll_real64,intent(in) ::t
sll_real64, intent(out) :: E2norm_ex
sll_real64 :: kx, eps,dx
sll_real64, dimension(:),pointer :: E,xmil
!sll_real64, dimension(nx) ::xmil
nx=sim%nc_x1
eps=sim%params(5)
!kx= 2.* sll_pi / (sim%params(2) - sim%params(1))
kx=0.2_f64
dx=sim%mesh2dx%delta_eta1
 SLL_ALLOCATE(E(0:nx-1),ierr)
 SLL_ALLOCATE(xmil(nx),ierr)
do i=1,nx
   xmil(i)=sim%mesh2dx%eta1_min+i*sim%mesh2dx%delta_eta1- &
        sim%mesh2dx%delta_eta1/2
end do
!write(*,*) 'entrer en solexact'
!write(*,*) 'kx =',kx, nx
    if (kx>0.19 .and. kx< 0.21) then
       do i=1,nx
          E(i-1)=4*eps*1.129664*exp(-0.0000551*t)*sin(kx*xmil(i))*cos(1.0640*t-0.00127377)
       enddo

    else if (kx>0.29 .and. kx< 0.31)then
       do i=1,nx

          E(i-1)=4*eps*0.63678*exp(-0.0126*t)*sin(kx*xmil(i))*cos(1.1598*t-0.114267)
       enddo
    else if (kx>0.39 .and. kx< 0.41) then
       do i=1,nx 
          E(i-1)=4*eps*0.424666*exp(-0.0661*t)*sin(kx*xmil(i))*cos(1.2850*t-0.3357725);
       enddo
    else if (kx>0.49 .and. kx< 0.51)then
       do i=1,nx
          E(i-1)=4*eps*0.3677*exp(-0.1533*t)*sin(kx*xmil(i))*cos(1.4156*t-0.536245);
       enddo
    endif

    ! L2-norm
    E2norm_ex=0.d0;
    do i=0,nx-1
       !E2norm_ex=E2norm_ex+E(i)*E(i)*dx*sim%mesh2dx%delta_eta2
       E2norm_ex=E2norm_ex+E(i)*E(i)*dx
    enddo
!write(*,*) 'nx = ', nx, dx
 SLL_DEALLOCATE_ARRAY(E,ierr)
 SLL_DEALLOCATE_ARRAY(xmil,ierr)
  end subroutine solexact


end module sll_simulation_4d_vp_eulerian_cartesian_finite_volume_module
