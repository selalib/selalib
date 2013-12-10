module sll_simulation_4d_drift_kinetic_cartesian_finite_volume
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
  use sll_gnuplot_parallel
  implicit none

  type, extends(sll_simulation_base_class) :: &
       sll_simulation_4d_drift_kinetic_cart_finite_volume
     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_v3  ! only processor vor the v direction
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     ! Mesh parameters
     sll_int32  :: nc_v3  ! velocity nodes
     sll_int32  :: nc_x1
     sll_int32  :: nc_x2
     sll_int32  :: nc_x3
     ! for initializers
     type(sll_logical_mesh_4d), pointer    :: mesh4d
     type(poisson_2d_periodic_plan_cartesian_par), pointer :: poisson_plan

     procedure(sll_scalar_initializer_4d), nopass, pointer :: init_func
     sll_real64, dimension(:), pointer :: params

     ! distribution functions at time steps n, star and n+1 
     ! communications are needed only in the x3 direction
     sll_real64, dimension(:,:,:,:), pointer     :: fn_v3x1x2
     sll_real64, dimension(:,:,:,:), pointer     :: fnp1_v3x1x2
     sll_real64, dimension(:,:,:,:), pointer     :: fstar_v3x1x2
     ! charge density
     sll_real64, dimension(:,:,:), allocatable     :: rho_x1x2
     ! potential 
     sll_real64, dimension(:,:,:), allocatable     :: phi_x1x2

     type(layout_4d),pointer :: sequential_v3x1x2
     type(layout_3d),pointer :: rho_seq_x1x2,phi_seq_x1x2
     
   contains
     procedure, pass(sim) :: run => run_dk_cart
     procedure, pass(sim) :: init_from_file => init_dk_cart
    
  end type sll_simulation_4d_drift_kinetic_cart_finite_volume

  interface delete
     module procedure delete_dk_cart
  end interface delete

contains

  subroutine init_dk_cart( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_drift_kinetic_cart_finite_volume), intent(inout) :: sim
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
       print *, 'init_dk_cart() failed to open file ', filename
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
  end subroutine init_dk_cart


  subroutine initialize_dk4d( &
   sim, &
   mesh4d, &
   init_func, &
   params )

   type(sll_simulation_4d_drift_kinetic_cart_finite_volume), intent(inout)     :: sim
   type(sll_logical_mesh_4d), pointer                    :: mesh4d
   procedure(sll_scalar_initializer_4d)                  :: init_func
   sll_real64, dimension(:), target                      :: params
   sim%mesh4d  => mesh4d
   sim%init_func => init_func
   sim%params    => params
 end subroutine initialize_dk4d


  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_dk_cart(sim)
    class(sll_simulation_4d_drift_kinetic_cart_finite_volume), intent(inout) :: sim
    sll_int32  :: loc_sz_v3
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: l
    sll_real64 :: vmin
    sll_real64 :: vmax
    sll_real64 :: dv
    sll_int32  :: ierr
    sll_int32  :: itime
    sll_int32  :: nc_v3
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: ranktop
    sll_int32  :: rankbottom
    sll_int32  :: message_id
    sll_int32  :: datasize
    sll_int32  :: istat
    sll_int32  :: tagtop,tagbottom

    
    sll_real64,dimension(:,:),allocatable :: plotf2d
    sll_int32,dimension(4)  :: global_indices

    sim%world_size = sll_get_collective_size(sll_world_collective)  
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)  

    ! allocate the layouts...
    sim%sequential_v3x1x2  => new_layout_4D( sll_world_collective )
    sim%rho_seq_x1x2       => new_layout_3D( sll_world_collective )
    sim%phi_seq_x1x2       => new_layout_3D( sll_world_collective )

    nc_v3 = sim%mesh4d%num_cells1
    nc_x1 = sim%mesh4d%num_cells2   
    nc_x2 = sim%mesh4d%num_cells3   
    nc_x3 = sim%mesh4d%num_cells4   

    sim%nproc_v3 = 1
    sim%nproc_x1 = 1
    sim%nproc_x2 = 1
    sim%nproc_x3 = sim%world_size
    
    ! init the layout for the distribution function
    ! the mesh is split only in the x3 direction
    call initialize_layout_with_distributed_4D_array( &
         nc_v3, &
         nc_x1, &
         nc_x2, &
         nc_x3, &
         sim%nproc_v3, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%sequential_v3x1x2)

    ! charge density layout
    call initialize_layout_with_distributed_3D_array( &
         nc_x1, &
         nc_x2, &
         nc_x3, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%rho_seq_x1x2)
    
    ! potential layout
    call initialize_layout_with_distributed_3D_array( &
         nc_x1, &
         nc_x2, &
         nc_x3, &
         sim%nproc_x1, &
         sim%nproc_x2, &
         sim%nproc_x3, &
         sim%phi_seq_x1x2)
    
    
    call compute_local_sizes_3d( sim%rho_seq_x1x2, loc_sz_x1, loc_sz_x2, &
         loc_sz_x3)

    ! iz=0 corresponds to the mean values of rho and phi 
    SLL_ALLOCATE(sim%rho_x1x2(loc_sz_x1,loc_sz_x2,0:loc_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi_x1x2(loc_sz_x1,loc_sz_x2,0:loc_sz_x3),ierr)

       
    ! Allocate the array needed to store the local chunk of the distribution
    ! function data. First compute the local sizes.
    call compute_local_sizes_4d( sim%sequential_v3x1x2, &
         loc_sz_v3, &
         loc_sz_x1, &
         loc_sz_x2, &
         loc_sz_x3 )

    ! iz=0 and iz=loc_sz_x3+1 correspond to ghost cells.
    SLL_ALLOCATE(sim%fn_v3x1x2(loc_sz_v3,loc_sz_x1,loc_sz_x2,0:loc_sz_x3+1),ierr)
    SLL_ALLOCATE(sim%fnp1_v3x1x2(loc_sz_v3,loc_sz_x1,loc_sz_x2,0:loc_sz_x3+1),ierr)
    SLL_ALLOCATE(sim%fstar_v3x1x2(loc_sz_v3,loc_sz_x1,loc_sz_x2,0:loc_sz_x3+1),ierr)
    
    
    
    

    ! initialize here the distribution function

    ! the function is passed by the user when the init_dk subroutine is called.
    ! The routine sll_4d_parallel_array_initializer_cartesian is in 
    ! src/parallal_array_initializers/sll_parallel_array_initializer_module.F90
    ! the particular initializer is in
    ! parallel_array_initializers/sll_common_array_initializers_module.F90

    call sll_4d_parallel_array_initializer_cartesian( &
         sim%sequential_v3x1x2, &
         sim%mesh4d, &
         sim%fn_v3x1x2(:,:,:,1:loc_sz_x3), &
         sim%init_func, &
         sim%params)

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



    ! mpi communications 

    ranktop=mod(sim%my_rank+1,sim%world_size)
    rankbottom=sim%my_rank-1
    if (rankbottom.lt.0) rankbottom=sim%world_size-1
    message_id=1
    datasize=loc_sz_v3*loc_sz_x1*loc_sz_x2
    tagtop=sim%my_rank
    tagbottom=ranktop

    ! top communications
    write(*,*) 'coucou1',sim%my_rank,' bottom:',rankbottom,' top:',ranktop,'datasize=',size(sim%fn_v3x1x2(:,:,:,loc_sz_x3+1))
    Call mpi_SENDRECV(sim%fn_v3x1x2(1,1,1,loc_sz_x3),datasize, &
         MPI_DOUBLE_PRECISION,ranktop,sim%my_rank,              &
         sim%fn_v3x1x2(1,1,1,0),datasize,            &
         MPI_DOUBLE_PRECISION,rankbottom,rankbottom,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)   

    ! bottom communications
    write(*,*) 'coucou2',sim%my_rank
    Call mpi_SENDRECV(sim%fn_v3x1x2(1,1,1,1),datasize, &
         MPI_DOUBLE_PRECISION,rankbottom,sim%my_rank,              &
         sim%fn_v3x1x2(1,1,1,loc_sz_x3+1),datasize,            &
         MPI_DOUBLE_PRECISION,ranktop,ranktop,              &
         MPI_COMM_WORLD,MPI_STATUS_IGNORE ,ierr)       

    write(*,*) 'end comm',sim%my_rank

    call compute_local_sizes_4d( sim%sequential_v3x1x2, &
         loc_sz_v3, loc_sz_x1, loc_sz_x2, loc_sz_x3) 
    
    allocate (plotf2d(loc_sz_v3,loc_sz_x1))

    do i = 1, loc_sz_x1
       do j = 1, loc_sz_v3
          plotf2d(j,i) = sim%fn_v3x1x2(j,i,1,0)
       end do
    end do
    
    global_indices(1:4) =  local_to_global_4D(sim%sequential_v3x1x2, (/1,1,1,1/) )
    
    call sll_gnuplot_rect_2d_parallel( &
         sim%mesh4d%eta1_min+(global_indices(1)-1)*sim%mesh4d%delta_eta1, &
         sim%mesh4d%delta_eta1, &
         sim%mesh4d%eta2_min+(global_indices(2)-1)*sim%mesh4d%delta_eta2, &
         sim%mesh4d%delta_eta2, &
         plotf2d, &
         "plotf2d", &
         0, &
         ierr)
    

  end subroutine run_dk_cart

  subroutine delete_dk_cart( sim )
    class(sll_simulation_4d_drift_kinetic_cart_finite_volume) :: sim
    sll_int32 :: ierr
    SLL_DEALLOCATE( sim%fn_v3x1x2, ierr )
    SLL_DEALLOCATE( sim%fnp1_v3x1x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%fstar_v3x1x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%rho_x1x2, ierr )
    SLL_DEALLOCATE_ARRAY( sim%phi_x1x2, ierr )
    call delete( sim%sequential_v3x1x2 )
    call delete( sim%rho_seq_x1x2 )
    call delete( sim%phi_seq_x1x2 )
  end subroutine delete_dk_cart

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
!!$    array_dims(1) = nc_x1
!!$    array_dims(2) = nc_x2
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

end module  sll_simulation_4d_drift_kinetic_cartesian_finite_volume
