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
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_4d_vp_eulerian_cartesian_finite_volume
     ! Parallel environment parameters
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
     sll_int32  :: num_iterations

     ! Mesh parameters
     sll_int32  :: nc_v1  ! velocity cells
     sll_int32  :: nc_v2
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2

     sll_int32 :: degree ! polynomial degree

     sll_int32  :: np_v1  ! velocity nodes
     sll_int32  :: np_v2

     ! number of local nodes in each element
     sll_int32 :: np_loc

     ! for initializers
     type(sll_logical_mesh_2d), pointer    :: mesh2dx,mesh2dv
     class(sll_coordinate_transformation_2d_base),pointer     :: tx,tv
     type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan

     procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params

     ! 1D interpolation points repartition
     sll_real64, dimension(:),pointer  :: interp_pts_1D 
     ! finite element approximation in the velocity space
     ! connectivity array
     sll_int32,dimension(:,:),pointer  :: connec
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
     sll_real64, dimension(:),pointer  :: M_diag,M_low,M_sup
     sll_real64, dimension(:),pointer  :: Av1_diag,Av1_low,Av1_sup
     sll_real64, dimension(:),pointer  :: Av2_diag,Av2_low,Av2_sup
     sll_real64, dimension(:),pointer  :: Bv1_diag,Bv1_low,Bv1_sup
     sll_real64, dimension(:),pointer  :: Bv2_diag,Bv2_low,Bv2_sup
     ! distribution functions at time steps n, star and n+1 
     ! communications are needed only in the x3 direction
     sll_real64, dimension(:,:,:,:), pointer     :: fn_v1v2x1
     sll_real64, dimension(:,:,:,:), pointer     :: fnp1_v1v2x1
     sll_real64, dimension(:,:,:,:), pointer     :: dtfn_v1v2x1
     ! charge density
     sll_real64, dimension(:,:), allocatable     :: rho_x1,Ex_x1,Ey_x1
     ! potential 
     sll_real64, dimension(:,:), allocatable     :: phi_x1

     type(layout_4d),pointer :: sequential_v1v2x1
     type(layout_2d),pointer :: phi_seq_x1
     
   contains
     procedure, pass(sim) :: run => run_vp_cart
     procedure, pass(sim) :: init_from_file => init_vp_cart
    
  end type sll_simulation_4d_vp_eulerian_cartesian_finite_volume

  interface delete
     module procedure delete_vp_cart
  end interface delete

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
    sll_int32  :: loc_sz_v1
    sll_int32  :: loc_sz_v2
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_real64 :: vmin
    sll_real64 :: vmax
    sll_real64 :: dv
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_int32  :: istat
    sll_int32  :: ic
    sll_int32  :: jc
    
    sll_real64,dimension(:,:),allocatable :: plotf2d,plotphi2d
    sll_real64 :: t
    sll_int32,dimension(4)  :: global_indices

    sim%world_size = sll_get_collective_size(sll_world_collective)  
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)  

    ! allocate the layouts...
    sim%sequential_v1v2x1  => new_layout_4D( sll_world_collective )
    sim%phi_seq_x1       => new_layout_2D( sll_world_collective )

    sim%degree=sim%params(6)

    sim%nc_v1 = sim%mesh2dv%num_cells1
    sim%nc_v2 = sim%mesh2dv%num_cells2   
    sim%nc_x1 = sim%mesh2dx%num_cells1   
    sim%nc_x2 = sim%mesh2dx%num_cells2  

    sim%np_v1 = sim%degree * sim%nc_v1 + 1
    sim%np_v2 = sim%degree * sim%nc_v2 + 1

    sim%np_v1 = sim%degree * sim%nc_v1 + 1
    sim%np_v2 = sim%degree * sim%nc_v2 + 1

    sim%nproc_v1 = 1
    sim%nproc_v2 = 1
    sim%nproc_x1 = 1
    sim%nproc_x2 = sim%world_size
    
    ! init the layout for the distribution function
    ! the mesh is split only in the x3 direction
    call initialize_layout_with_distributed_4D_array( &
         sim%np_v1, &
         sim%np_v2, &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nproc_v1, &
         sim%nproc_v2, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%sequential_v1v2x1)


    
    ! potential layout
    call initialize_layout_with_distributed_2D_array( &
         sim%nc_x1, &
         sim%nc_x2, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%phi_seq_x1)

    sim%poisson_plan=>new_poisson_2d_periodic_plan_cartesian_par( &
    sim%phi_seq_x1, &
    sim%nc_x1, &
    sim%nc_x2, &
    sim%mesh2dx%eta1_max-sim%mesh2dx%eta1_min, &
    sim%mesh2dx%eta2_max-sim%mesh2dx%eta2_min )
    
    
    call compute_local_sizes_2d( sim%phi_seq_x1, loc_sz_x1, loc_sz_x2)
    ! iz=0 corresponds to the mean values of rho and phi 
    SLL_ALLOCATE(sim%rho_x1(loc_sz_x1,loc_sz_x2),ierr)
    SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,0:loc_sz_x2+1),ierr)
    !SLL_ALLOCATE(sim%phi_x1(loc_sz_x1,loc_sz_x2),ierr)
    print *, 'rank = ', sim%my_rank, 'local sizes of phi: ', loc_sz_x1, loc_sz_x2
       
    ! Allocate the array needed to store the local chunk of the distribution
    ! function data. First compute the local sizes.
    call compute_local_sizes_4d( sim%sequential_v1v2x1, &
         loc_sz_v1, &
         loc_sz_v2, &
         loc_sz_x1, &
         loc_sz_x2 )

    write(*,*) 'sim%np_v1',sim%np_v1,loc_sz_v1

    ! iz=0 and iz=loc_sz_x2+1 correspond to ghost cells.
    SLL_ALLOCATE(sim%fn_v1v2x1(loc_sz_v1,loc_sz_v2,loc_sz_x1,0:loc_sz_x2+1),ierr)
    SLL_ALLOCATE(sim%fnp1_v1v2x1(loc_sz_v1,loc_sz_v2,loc_sz_x1,0:loc_sz_x2+1),ierr)
    SLL_ALLOCATE(sim%dtfn_v1v2x1(loc_sz_v1,loc_sz_v2,loc_sz_x1,0:loc_sz_x2+1),ierr)
    
    
    
    

    ! initialize here the distribution function

    ! the function is passed by the user when the init_vp subroutine is called.
    ! The routine sll_4d_parallel_array_initializer_cartesian is in 
    ! src/parallal_array_initializers/sll_parallel_array_initializer_module.F90
    ! the particular initializer is in
    ! parallel_array_initializers/sll_common_array_initializers_module.F90

!!$    call sll_4d_parallel_array_initializer_cartesian( &
!!$         sim%sequential_v1v2x1, &
!!$         sim%mesh4d, &
!!$         sim%fn_v1v2x1(:,:,:,1:loc_sz_x2), &
!!$         sim%init_func, &
!!$         sim%params)
    call sll_4d_parallel_array_initializer( &
         sim%sequential_v1v2x1, &
         sim%mesh2dv, &
         sim%mesh2dx, &
         sim%fn_v1v2x1(:,:,:,1:loc_sz_x2), &
         sim%init_func , &
         sim%params, &
         sim%tv, &
         sim%tx, &
         sim%degree, &
         sim%degree)
    

    call velocity_mesh_connectivity(sim)


    !
    do i=1,loc_sz_x1
       do j=1,loc_sz_x2
          sim%rho_x1(i,j)=sum(sim%fn_v1v2x1(:,:,i,j))
          !sim%rho_x1(i,j)=1.d0
       enddo
    enddo

!!$    do i=1,loc_sz_x1
!!$       do j=1,loc_sz_x2
!!$          sim%rho_x1(i,j)=0_f64
!!$       enddo
!!$    enddo

    ! solve the poisson equation
    call solve_poisson_2d_periodic_cartesian_par(sim%poisson_plan, sim%rho_x1,&
         sim%phi_x1(1:loc_sz_x1,1:loc_sz_x2))


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
    

!!$    ! simple computation
!!$    ! to be generalized with generic transformation...
!!$    volume=sim%mesh2dx%delta_eta1 * &
!!$         sim%mesh2dx%delta_eta2 

    global_indices(1:4)=local_to_global_4D(sim%sequential_v1v2x1, (/1,1,1,1/) )
! cell volumes
    do j=1,loc_sz_x2
       do i=1,loc_sz_x1
          ic=i+global_indices(3)-1
          jc=j+global_indices(4)-1
          sim%volume(i,j)=cell_volume( sim%tx, ic, jc,3)          
       end do
    end do
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

    !write(*,*) 'vertical',sim%my_rank,sum(sim%surfx1(1,:))
    !write(*,*) 'horizontal',sim%my_rank,sum(sim%surfx2(1,:))


    ! time loop
    t=0
    !compute the time step
    sim%cfl=sim%params(7)
!!$    write(*,*) 'Vxmax = ', sim%mesh2dv%eta1_max
!!$    write(*,*) 'Vxmax = ', sim%mesh2dv%eta2_max
!!$    write(*,*) 'deltav = ', sim%mesh2dx%delta_eta1
!!$    write(*,*) 'surf/volum =  ', (sim%surfx1(1,1)+sim%surfx2(1,1))/ &
!!$         sim%volume(1,1)*sim%mesh2dx%delta_eta1
    sim%dt = sim%cfl*sim%volume(1,1)/(sim%surfx1(1,1)+sim%surfx2(1,1))/ &
         max(sim%mesh2dv%eta1_max,sim%mesh2dv%eta2_max)
    write(*,*) 'dt = ', sim%dt
    do while(t.lt.sim%tmax)
       t=t+sim%dt
       call RK2(sim)
    sim%fn_v1v2x1 = sim%fnp1_v1v2x1
   end do

    call compute_local_sizes_4d( sim%sequential_v1v2x1, &
         loc_sz_v1, loc_sz_v2, loc_sz_x1, loc_sz_x2) 
    
    allocate (plotf2d(loc_sz_x1,loc_sz_v1))

    do i = 1, loc_sz_x1
       do j = 1, loc_sz_v1
          plotf2d(i,j) = sim%fn_v1v2x1(j,1,i,1)
       end do
    end do
    
    allocate (plotphi2d(loc_sz_x1,loc_sz_x2))

    do i = 1, loc_sz_x1
       do j = 1, loc_sz_x2
          plotphi2d(i,j) = sim%phi_x1(i,j)
       end do
    end do

    global_indices(1:4) =  local_to_global_4D(sim%sequential_v1v2x1, (/1,1,1,1/) )
    
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
         sim%mesh2dx%delta_eta1, &
         sim%mesh2dv%eta1_min+(global_indices(1)-1)*sim%mesh2dv%delta_eta1/sim%degree, &
         sim%mesh2dv%delta_eta1/sim%degree, &
         plotf2d, &
         "plotf2d", &
         0, &
         ierr)
!!$
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh2dx%eta1_min+(global_indices(3)-1)*sim%mesh2dx%delta_eta1, &
         sim%mesh2dx%delta_eta1, &
         sim%mesh2dx%eta2_min+(global_indices(4)-1)*sim%mesh2dx%delta_eta2, &
         sim%mesh2dx%delta_eta2, &
         plotphi2d, &
         "plotphi2d", &
         0, &
         ierr)
    

  end subroutine run_vp_cart

  subroutine velocity_mesh_connectivity(sim)

    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume) :: sim
    sll_int32 :: ierr
    sll_int32 :: i,j,k,ii,jj
    sll_int32 :: ic,jc,iploc,jploc,ino
    sll_real64  :: x,y,xref,yref
    sll_real64 :: det
    sll_real64 :: phi1,phi2,dphi1(2),dphi2(2),dphi1ref(2),dphi2ref(2)
    sll_real64,dimension(2,2) :: jacob,invjacob
    sll_real64,dimension(:,:),allocatable :: lag,dlag
    sll_real64,dimension(:,:),allocatable :: mloc,av1loc,av2loc,bv1loc,bv2loc
    sll_real64,dimension(:),allocatable :: gauss,weight
    sll_real64 :: void  ! only for a valid address
    sll_int32 :: ifac,isol,nsym,mp
    sll_int32 :: ll,ib1,ib2,jb1,jb2

    SLL_ALLOCATE(sim%interp_pts_1D(sim%degree+1),ierr)

    ! for the moment: equally spaced points
    do i=0,sim%degree
       sim%interp_pts_1D(i+1)=real(i,f64)/sim%degree
    end do

    ! array of points coordinates
    SLL_ALLOCATE(sim%vcoords(2,sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%vcoords(2,sim%np_v1*sim%np_v2),ierr)
    
    ! connectivity
    sim%np_loc=(sim%degree+1)**2
    SLL_ALLOCATE(sim%connec(sim%np_loc,sim%nc_v1*sim%nc_v2),ierr)
    
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
    
!!$    do ic=1,sim%nc_v1*sim%nc_v2
!!$       write(*,*) sim%connec(:,ic) 
!!$    end do
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

    sim%nsky=sim%mkld(sim%np_v1*sim%np_v2+1)-1

    SLL_ALLOCATE(sim%M_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%M_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Av1_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av1_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av1_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Av2_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av2_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Av2_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv1_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv1_diag(sim%np_v1*sim%np_v2),ierr)

    SLL_ALLOCATE(sim%Bv2_low(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2_sup(sim%nsky),ierr)
    SLL_ALLOCATE(sim%Bv2_diag(sim%np_v1*sim%np_v2),ierr)

    sim%M_low=0
    sim%M_sup=0
    sim%M_diag=0

    sim%Av1_low=0
    sim%Av1_sup=0
    sim%Av1_diag=0

    sim%Av2_low=0
    sim%Av2_sup=0
    sim%Av2_diag=0

    sim%Bv1_low=0
    sim%Bv1_sup=0
    sim%Bv1_diag=0

    sim%Bv2_low=0
    sim%Bv2_sup=0
    sim%Bv2_diag=0
    


!!$    write(*,*) sim%prof
!!$
!!$    stop

    ! init of Gauss points
    SLL_ALLOCATE(gauss(sim%degree+1),ierr)
    SLL_ALLOCATE(weight(sim%degree+1),ierr)
    SLL_ALLOCATE(lag(sim%degree+1,sim%degree+1),ierr)
    SLL_ALLOCATE(dlag(sim%degree+1,sim%degree+1),ierr)

    SLL_ALLOCATE(mloc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(av1loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(av2loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(bv1loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    SLL_ALLOCATE(bv2loc((sim%degree+1)**2,(sim%degree+1)**2),ierr)
    call lag_gauss(sim%degree,gauss,weight,lag,dlag)


    ! matrix assembly

    ! loop on the cells
    do ic=0,sim%nc_v1-1
       do jc=0,sim%nc_v2-1
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
               ! write(*,*) 'xref,yref=',xref,yref
                jacob=sim%tv%jacobian_matrix(xref,yref)
                invjacob=sim%tv%inverse_jacobian_matrix(xref,yref)
                det=sim%tv%jacobian(xref,yref)*sim%mesh2dv%delta_eta1*sim%mesh2dv%delta_eta2

!!$                write(*,*) 'det=',det
!!$                write(*,*) 'jacob=',jacob
!!$                write(*,*) 'invjacob=',invjacob
!!$                write(*,*) 'gauss=',weight
!!$                write(*,*)

                do ib1=1,sim%degree+1
                   do jb1=1,sim%degree+1
                      do ib2=1,sim%degree+1
                         do jb2=1,sim%degree+1
                            phi1=lag(ib1,iploc+1)*lag(jb1,jploc+1)
                            phi2=lag(ib2,iploc+1)*lag(jb2,jploc+1)
!!$                            write(*,*) 'phi=',xref,yref,'i,j',(jb1-1)*(sim%degree+1)+ib1, &
!!$                                 (jb2-1)*(sim%degree+1)+ib2,phi1,phi2
                            dphi1ref(1)=dlag(ib1,iploc+1)*lag(jb1,jploc+1)
                            dphi1ref(2)=lag(ib1,iploc+1)*dlag(jb1,jploc+1)
                            dphi2ref(1)=dlag(ib2,iploc+1)*lag(jb2,jploc+1)
                            dphi2ref(2)=lag(ib2,iploc+1)*dlag(jb2,jploc+1)
                            dphi1(1)=dphi1ref(1)*invjacob(1,1)+dphi1ref(2)*invjacob(2,1)
                            dphi1(2)=dphi1ref(1)*invjacob(1,2)+dphi1ref(2)*invjacob(2,2)
                            dphi2(1)=dphi2ref(1)*invjacob(1,1)+dphi2ref(2)*invjacob(2,1)
                            dphi2(2)=dphi2ref(1)*invjacob(1,2)+dphi2ref(2)*invjacob(2,2)
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
                                dphi1(1)*phi2*det*weight(iploc+1)*weight(jploc+1)  
                            bv2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)=&
                                 bv2loc((jb1-1)*(sim%degree+1)+ib1,(jb2-1)*(sim%degree+1)+ib2)+&
                                dphi1(2)*phi2*det*weight(iploc+1)*weight(jploc+1)                    
                         end do
                      end do
                   end do
                end do

             end do
          end do

!!$          do ii=1,(sim%degree+1)**2
!!$             write(*,*) 'bv1loc=',bv1loc(ii,:)
!!$          end do
             
          do ii=1,(sim%degree+1)**2
             do jj=1,(sim%degree+1)**2
                i=sim%connec(ii,jc*sim%nc_v1+ic+1)
                j=sim%connec(jj,jc*sim%nc_v1+ic+1)
                !if (cal.eq.1) then
                if (i.eq.j) then
                   sim%M_diag(i)=sim%M_diag(i)+mloc(ii,jj)
                   sim%Av1_diag(i)=sim%Av1_diag(i)+av1loc(ii,jj)
                   sim%Av2_diag(i)=sim%Av2_diag(i)+av2loc(ii,jj)
                   sim%Bv1_diag(i)=sim%Bv1_diag(i)+bv1loc(ii,jj)
                   sim%Bv2_diag(i)=sim%Bv2_diag(i)+bv2loc(ii,jj)
                else if (j.gt.i) then
                   ll=sim%mkld(j+1)-j+i
                   sim%M_sup(ll)=sim%M_sup(ll)+mloc(ii,jj)
                   sim%Av1_sup(ll)=sim%Av1_sup(ll)+av1loc(ii,jj)
                   sim%Av2_sup(ll)=sim%Av2_sup(ll)+av2loc(ii,jj)
                   sim%Bv1_sup(ll)=sim%Bv1_sup(ll)+bv1loc(ii,jj)
                   sim%Bv2_sup(ll)=sim%Bv2_sup(ll)+bv2loc(ii,jj)
                else
                   ll=sim%mkld(i+1)-i+j
                   sim%M_low(ll)=sim%M_low(ll)+mloc(ii,jj)
                   sim%Av1_low(ll)=sim%Av1_low(ll)+av1loc(ii,jj)
                   sim%Av2_low(ll)=sim%Av2_low(ll)+av1loc(ii,jj)
                   sim%Bv1_low(ll)=sim%Bv1_low(ll)+bv1loc(ii,jj)
                   sim%Bv2_low(ll)=sim%Bv2_low(ll)+bv2loc(ii,jj)
                end if
             end do
          end do

          ! end loop on the cells
       end do
    end do


!!$    write(*,*) 'A1_sup=',sim%Av1_sup
!!$    write(*,*) 'A1_low=',sim%Av1_low
!!$    stop

    ! LU decomposition of M
    ifac=1  ! we compute the LU decomposition
    isol=0  ! we do not solve the linear system
    nsym=1  ! we do not take into account the symetry of M
    mp=6    ! write the log on screen

    call sol(sim%M_sup,sim%M_diag,sim%M_low,void,&
         sim%mkld,void,sim%np_v1*sim%np_v2,mp,ifac,isol,nsym,void,ierr,&
         sim%nsky)
    


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
    SLL_DEALLOCATE( sim%fn_v1v2x1, ierr )
    SLL_DEALLOCATE( sim%fnp1_v1v2x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%dtfn_v1v2x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1, ierr )
    call delete( sim%sequential_v1v2x1 )
    call delete( sim%phi_seq_x1 )
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

    case(1)
       gauss(1) = -0.5773502691896257645091487805019574556476D0
       gauss(2) = 0.5773502691896257645091487805019574556476D0
       weight(1) = 0.1000000000000000000000000000000000000000D1
       weight(2) = 0.1000000000000000000000000000000000000000D1
       lag(1,1) = 0.7886751345948128822545743902509787278238D0
       lag(1,2) = 0.2113248654051871177454256097490212721762D0
       lag(2,1) = 0.2113248654051871177454256097490212721762D0
       lag(2,2) = 0.7886751345948128822545743902509787278238D0
       dlag(1,1) = -0.5000000000000000000000000000000000000000D0
       dlag(1,2) = -0.5000000000000000000000000000000000000000D0
       dlag(2,1) = 0.5000000000000000000000000000000000000000D0
       dlag(2,2) = 0.5000000000000000000000000000000000000000D0
 
    case(2)
       gauss(1) = -0.7745966692414833770358530799564799221666D0
       gauss(2) = 0.0D0
       gauss(3) = 0.7745966692414833770358530799564799221666D0
       weight(1) = 0.5555555555555555555555555555555555555560D0
       weight(2) = 0.8888888888888888888888888888888888888888D0
       weight(3) = 0.5555555555555555555555555555555555555560D0
       lag(1,1) = 0.6872983346207416885179265399782399610833D0
       lag(1,3) = -0.8729833462074168851792653997823996108329D-1
       lag(2,1) = 0.4000000000000000000000000000000000000001D0
       lag(2,2) = 0.1D1
       lag(2,3) = 0.4000000000000000000000000000000000000001D0
       lag(3,1) = -0.8729833462074168851792653997823996108329D-1
       lag(3,3) = 0.6872983346207416885179265399782399610833D0
       dlag(1,1) = -0.1274596669241483377035853079956479922167D1
       dlag(1,2) = -0.5000000000000000000000000000000000000000D0
       dlag(1,3) = 0.2745966692414833770358530799564799221666D0
       dlag(2,1) = 0.1549193338482966754071706159912959844333D1
       dlag(2,3) = -0.1549193338482966754071706159912959844333D1
       dlag(3,1) = -0.2745966692414833770358530799564799221666D0
       dlag(3,2) = 0.5000000000000000000000000000000000000000D0
       dlag(3,3) = 0.1274596669241483377035853079956479922167D1

    case(3)
       gauss(1) = -0.8611363115940525752239464888928095050957D0
       gauss(2) = -0.3399810435848562648026657591032446872006D0
       gauss(3) = 0.3399810435848562648026657591032446872006D0
       gauss(4) = 0.8611363115940525752239464888928095050957D0
       weight(1) = 0.3478548451374538573730639492219994072349D0
       weight(2) = 0.6521451548625461426269360507780005927648D0
       weight(3) = 0.6521451548625461426269360507780005927648D0
       weight(4) = 0.3478548451374538573730639492219994072349D0
       lag(1,1) = 0.6600056650728035304586090241410712591072D0
       lag(1,2) = 0.3373736432772532230452701535743525193130D-2
       lag(1,3) = 0.1661762313906394832923866663740817265855D-2
       lag(1,4) = 0.4924455046623182819230012194515868414856D-1
       lag(2,1) = 0.5209376877117036301110018070022304235776D0
       lag(2,2) = 0.1004885854825645727576017102247120797681D1
       lag(2,3) = -0.9921353572324654639393670446605140139452D-2
       lag(2,4) = -0.2301879032507389887619109530884603668341D0
       lag(3,1) = -0.2301879032507389887619109530884603668341D0
       lag(3,2) = -0.9921353572324654639393670446605140139450D-2
       lag(3,3) = 0.1004885854825645727576017102247120797682D1
       lag(3,4) = 0.5209376877117036301110018070022304235776D0
       lag(4,1) = 0.4924455046623182819230012194515868414856D-1
       lag(4,2) = 0.1661762313906394832923866663740817265855D-2
       lag(4,3) = 0.3373736432772532230452701535743525193130D-2
       lag(4,4) = 0.6600056650728035304586090241410712591073D0
       dlag(1,1) = -0.2157653673851862185103303519133755608117D1
       dlag(1,2) = -0.5150319221529816884980638312903767867892D0
       dlag(1,3) = 0.2499254259129449073079341266919237594124D0
       dlag(1,4) = -0.2200969727652438908494239191249342216503D0
       dlag(2,1) = 0.3035404320468968261056030957392445437885D1
       dlag(2,2) = -0.7198615816069815303118064641111701858325D0
       dlag(2,3) = -0.1484818929672908126117804422093470732036D1
       dlag(2,4) = 0.1097847619382349966802151357383624051417D1
       dlag(3,1) = -0.1097847619382349966802151357383624051417D1
       dlag(3,2) = 0.1484818929672908126117804422093470732036D1
       dlag(3,3) = 0.719861581606981530311806464111170185832D0
       dlag(3,4) = -0.3035404320468968261056030957392445437884D1
       dlag(4,1) = 0.2200969727652438908494239191249342216502D0
       dlag(4,2) = -0.2499254259129449073079341266919237594124D0
       dlag(4,3) = 0.5150319221529816884980638312903767867892D0
       dlag(4,4) = 0.2157653673851862185103303519133755608117D1

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
      lag(1,1) = 0.6577278825775883905959416853145508163014D0
      lag(1,2) = 0.2206310329510026462660372235361697772283D-1
      lag(1,4) = -0.6618786099995526088371303768558919770637D-2
      lag(1,5) = -0.3237267008427455182670790754452363027997D-1
      lag(2,1) = 0.6076926946610149906004284764305419175436D0
      lag(2,2) = 0.1058797182171758427703178033099861242748D1
      lag(2,4) = 0.3922234075058382706049924394379542780120D-1
      lag(2,5) = 0.1755341081074128898504739055499048804553D0
      lag(3,1) = -0.4085820152617417192201361597504739840206D0
      lag(3,2) = -0.1134638401174469933019096956287147285027D0
      lag(3,3) = 0.1000000000000000000000000000000000000000D1
      lag(3,4) = -0.1134638401174469933019096956287147285028D0
      lag(3,5) = -0.4085820152617417192201361597504739840207D0
      lag(4,1) = 0.1755341081074128898504739055499048804554D0
      lag(4,2) = 0.3922234075058382706049924394379542780120D-1
      lag(4,4) = 0.1058797182171758427703178033099861242748D1
      lag(4,5) = 0.6076926946610149906004284764305419175436D0
      lag(5,1) = -0.3237267008427455182670790754452363027997D-1
      lag(5,2) = -0.6618786099995526088371303768558919770637D-2
      lag(5,4) = 0.2206310329510026462660372235361697772283D-1
      lag(5,5) = 0.6577278825775883905959416853145508163015D0
      dlag(1,1) = -0.3157918213674122166875093558348224319248D1
      dlag(1,2) = -0.6500852780101332162336916681696991811558D0
      dlag(1,3) = 0.1666666666666666666666666666666666666667D0
      dlag(1,4) = -0.1763781803592946593495819498808563878701D0
      dlag(1,5) = 0.2066038942657722646805893986210021104953D0
      dlag(2,1) = 0.5055639151533482794115362560468103207732D1
      dlag(2,2) = -0.1379999586751627374793764628172594064273D1
      dlag(2,3) = -0.1333333333333333333333333333333333333333D1
      dlag(2,4) = 0.1032926503490483125960311864273705202324D1
      dlag(2,5) = -0.1153010512716782989726354241013658790230D1
      dlag(3,1) = -0.2844127556310371352286033844512535568221D1
      dlag(3,2) = 0.2886633187892949057638186210735142059883D1
      dlag(3,4) = -0.2886633187892949057638186210735142059882D1
      dlag(3,5) = 0.2844127556310371352286033844512535568222D1
      dlag(4,1) = 0.1153010512716782989726354241013658790231D1
      dlag(4,2) = -0.1032926503490483125960311864273705202324D1
      dlag(4,3) = 0.1333333333333333333333333333333333333334D1
      dlag(4,4) = 0.1379999586751627374793764628172594064273D1
      dlag(4,5) = -0.5055639151533482794115362560468103207737D1
      dlag(5,1) = -0.2066038942657722646805893986210021104953D0
      dlag(5,2) = 0.1763781803592946593495819498808563878701D0
      dlag(5,3) = -0.1666666666666666666666666666666666666667D0
      dlag(5,4) = 0.6500852780101332162336916681696991811557D0
      dlag(5,5) = 0.3157918213674122166875093558348224319247D1

    case default
       write(*,*) 'degree ',degree,' not implemented !'
       stop

    end select


    ! remap to interval [0,1]
    gauss=(gauss+1)/2
    weight=weight/2
    dlag=dlag*2


  end subroutine lag_gauss


  subroutine sourcenum(sim,Ex,Ey,w,source)
    
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(in) :: sim   
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(in) :: w
    sll_real64,intent(in)  :: Ex,Ey
    sll_real64,dimension(:,:),intent(out) :: source

    sll_real64,dimension(:,:),allocatable :: source1,source2
    sll_int32 :: ierr

    SLL_ALLOCATE(source1(sim%np_v1,sim%np_v2),ierr)
    SLL_ALLOCATE(source2(sim%np_v1,sim%np_v2),ierr)


    source1=0
    source2=0
    
    call MULKU(sim%Bv1_sup,sim%Bv1_diag,sim%Bv1_low, &
         sim%mkld,w,sim%np_v1*sim%np_v2,1,source1,sim%nsky)
    call MULKU(sim%Bv2_sup,sim%Bv2_diag,sim%Bv2_low, &
         sim%mkld,w,sim%np_v1*sim%np_v2,1,source2,sim%nsky)

    source=Ex*source1+Ey*source2

!!$    !can we do as following ? so we have to call only one 
!!$    !time the subroutine mulk
!!$    source=0
!!$    call MULKU(sim%Bv1_sup*Ex+Ey*sim%Bv2_sup,sim%Bv1_diag*Ex+Ey*sim%Bv2_diag, &
!!$         sim%Bv1_low*Ex+Ey*sim%Bv2_low, &
!!$         sim%mkld,w,sim%np_v1*sim%np_v2,1,source,sim%nsky)


    SLL_DEALLOCATE_ARRAY(source1,ierr)
    SLL_DEALLOCATE_ARRAY(source2,ierr)

  end subroutine sourcenum

  subroutine fluxnum(sim,wL,wR,vn,flux)
    
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(in) :: sim   
    sll_real64,dimension(2),intent(in) :: vn
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(in) :: wL,wR
    sll_real64,dimension(sim%np_v1*sim%np_v2),intent(out) :: flux

    sll_real64,dimension(sim%np_v1*sim%np_v2) :: flux1,flux2
    sll_real64,dimension(sim%np_v1*sim%np_v2) :: wm
    sll_real64,dimension(sim%np_v1*sim%np_v2) :: temp

    sll_real64   :: eps=0.01

    wm=(wL+wR)/2

    flux1=0
    flux2=0
    
    call MULKU(sim%Av1_sup,sim%Av1_diag,sim%Av1_low, &
         sim%mkld,wm,sim%np_v1*sim%np_v2,1,flux1,sim%nsky)
    call MULKU(sim%Av2_sup,sim%Av2_diag,sim%Av2_low, &
         sim%mkld,wm,sim%np_v1*sim%np_v2,1,flux2,sim%nsky)

    flux=vn(1)*flux1+vn(2)*flux2-eps/2*(wR-wL)


  end subroutine fluxnum

  ! time derivative of f
  subroutine dtf(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim   
    
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: ierr,ifac,isol,nsym,mp
    sll_real64 :: void,Ex,Ey,vn(2)
    sll_int32 :: ic,jc,icL,icR,jcL,jcR

    sll_real64,dimension(:,:),allocatable  :: temp,flux,source

    SLL_ALLOCATE(flux(sim%np_v1,sim%np_v2),ierr)
    SLL_ALLOCATE(temp(sim%np_v1,sim%np_v2),ierr)
    SLL_ALLOCATE(source(sim%np_v1,sim%np_v2),ierr)

        ! solve  Mvx=vf
    ifac=0  ! we do not compute the LU decomposition
    isol=1  ! we do solve the linear system
    nsym=1  ! we do not take into account the symetry of M
    mp=6    ! write the log on screen

    
    call compute_local_sizes_2d( sim%phi_seq_x1, loc_sz_x1, loc_sz_x2)
    ! init
    sim%dtfn_v1v2x1=0

    !compute the fluxes in the x1 direction
    vn(1)=1*sim%surfx1(1,1) ! temporaire !!!!
    vn(2)=0
    do jc=1,loc_sz_x2
       do ic=0,loc_sz_x1-1
          icL=ic
          if (icL.le.0) icL=loc_sz_x1
          icR=ic+1
          call fluxnum(sim,sim%fn_v1v2x1(:,:,icL,jc), &
               sim%fn_v1v2x1(:,:,icR,jc),vn,flux)
!!$          sim%dtfn_v1v2x1(:,:,icL,jc)=sim%dtfn_v1v2x1(:,:,icL,jc) &
!!$               -flux/sim%mesh2dx%delta_eta1
!!$          sim%dtfn_v1v2x1(:,:,icR,jc)=sim%dtfn_v1v2x1(:,:,icR,jc) &
!!$               +flux/sim%mesh2dx%delta_eta1
          sim%dtfn_v1v2x1(:,:,icL,jc)=sim%dtfn_v1v2x1(:,:,icL,jc)-flux
          sim%dtfn_v1v2x1(:,:,icR,jc)=sim%dtfn_v1v2x1(:,:,icR,jc)+flux
       end do
    end do

    !compute the fluxes in the x2 direction
    vn(1)=0 ! temporaire !!!!
    vn(2)=1*sim%surfx2(1,1)
    do ic=1,loc_sz_x1
       do jc=0,loc_sz_x2-1
          jcL=jc
          jcR=jc+1
          call fluxnum(sim,sim%fn_v1v2x1(:,:,ic,jcL), &
               sim%fn_v1v2x1(:,:,ic,jcR),vn,flux)
!!$          sim%dtfn_v1v2x1(:,:,ic,jcL)=sim%dtfn_v1v2x1(:,:,ic,jcL) &
!!$               -flux/sim%mesh2dx%delta_eta2
!!$          sim%dtfn_v1v2x1(:,:,ic,jcR)=sim%dtfn_v1v2x1(:,:,ic,jcR) &
!!$               +flux/sim%mesh2dx%delta_eta2
          sim%dtfn_v1v2x1(:,:,ic,jcL)=sim%dtfn_v1v2x1(:,:,ic,jcL)-flux
          sim%dtfn_v1v2x1(:,:,ic,jcR)=sim%dtfn_v1v2x1(:,:,ic,jcR)+flux
 
       end do
    end do

    ! source terms
    do ic=1,loc_sz_x1
       do jc=1,loc_sz_x2
          icL=ic-1
          icR=ic+1
          if(ic.le.1) then
             icL=loc_sz_x1
          elseif(ic.ge.loc_sz_x1)then
             icR=1
          end if
          Ex=-(sim%phi_x1(icR,jc)-sim%phi_x1(icL,jc))/2/sim%mesh2dx%delta_eta1
          Ey=-(sim%phi_x1(ic,jc+1)-sim%phi_x1(ic,jc-1))/2/sim%mesh2dx%delta_eta2
          call sourcenum(sim,Ex,Ey,sim%fn_v1v2x1(:,:,ic,jc), &
               source)

          temp=sim%dtfn_v1v2x1(:,:,ic,jc)+sim%volume(1,1)*source
          call sol(sim%M_sup,sim%M_diag,sim%M_low,temp, &
               sim%mkld, &
               sim%dtfn_v1v2x1(:,:,ic,jc),&
               sim%np_v1*sim%np_v2, &
               mp,ifac,isol,nsym,void,ierr,&
               sim%nsky)
          !if we call the sol here, it means that we have to call many times
          !so if we call the sol out of the loop so we have to call it only 
          !one time and we have juste modify the temp=>vector,no?
       end do
    end do


    SLL_DEALLOCATE_ARRAY(flux,ierr)
    SLL_DEALLOCATE_ARRAY(temp,ierr)



  end subroutine dtf

  subroutine euler(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
    sim%fnp1_v1v2x1=0
    ! mpi communications
    call mpi_comm(sim)
    call dtf(sim)
    sim%fnp1_v1v2x1 = (sim%volume(1,1)*sim%fn_v1v2x1 &
         + sim%dt*sim%dtfn_v1v2x1)/sim%volume(1,1)
  end subroutine euler

  subroutine RK2(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim
    sim%fnp1_v1v2x1=0
    ! mpi communications
    call mpi_comm(sim)
    call dtf(sim)
    sim%fnp1_v1v2x1 = sim%fnp1_v1v2x1+sim%fn_v1v2x1
    sim%fn_v1v2x1 = (sim%volume(1,1)*sim%fn_v1v2x1 &
         + sim%dt*sim%dtfn_v1v2x1)/sim%volume(1,1)/2
    call mpi_comm(sim)
    call dtf(sim)
    sim%fnp1_v1v2x1 = (sim%volume(1,1)*sim%fnp1_v1v2x1 &
         + sim%dt*sim%dtfn_v1v2x1)/sim%volume(1,1)
  end subroutine RK2

  subroutine mpi_comm(sim)
    class(sll_simulation_4d_vp_eulerian_cartesian_finite_volume), intent(inout) :: sim   
    sll_int32 :: ierr
    sll_int32  :: loc_sz_v1
    sll_int32  :: loc_sz_v2
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: ranktop
    sll_int32  :: rankbottom
    sll_int32  :: datasize,datasizephi

    call compute_local_sizes_4d( sim%sequential_v1v2x1, &
         loc_sz_v1, &
         loc_sz_v2, &
         loc_sz_x1, &
         loc_sz_x2 )
    ranktop=mod(sim%my_rank+1,sim%world_size)
    rankbottom=sim%my_rank-1
    if (rankbottom.lt.0) rankbottom=sim%world_size-1
    datasize=loc_sz_v1*loc_sz_v2*loc_sz_x1
    datasizephi=loc_sz_x1

    Call mpi_SENDRECV(sim%fn_v1v2x1(1,1,1,loc_sz_x2),datasize, &
         MPI_DOUBLE_PRECISION,ranktop,sim%my_rank,              &
         sim%fn_v1v2x1(1,1,1,0),datasize,            &
         MPI_DOUBLE_PRECISION,rankbottom,rankbottom,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)   

    Call mpi_SENDRECV(sim%phi_x1(1,loc_sz_x2),datasizephi, &
         MPI_DOUBLE_PRECISION,ranktop,sim%my_rank,              &
         sim%phi_x1(1,0),datasizephi,            &
         MPI_DOUBLE_PRECISION,rankbottom,rankbottom,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)   

    ! bottom communications
    Call mpi_SENDRECV(sim%fn_v1v2x1(1,1,1,1),datasize, &
         MPI_DOUBLE_PRECISION,rankbottom,sim%my_rank,              &
         sim%fn_v1v2x1(1,1,1,loc_sz_x2+1),datasize,            &
         MPI_DOUBLE_PRECISION,ranktop,ranktop,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)       

    Call mpi_SENDRECV(sim%phi_x1(1,1),datasizephi, &
         MPI_DOUBLE_PRECISION,rankbottom,sim%my_rank,              &
         sim%phi_x1(1,loc_sz_x2+1),datasizephi,            &
         MPI_DOUBLE_PRECISION,ranktop,ranktop,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)       



  end subroutine mpi_comm
end module sll_simulation_4d_vp_eulerian_cartesian_finite_volume_module
