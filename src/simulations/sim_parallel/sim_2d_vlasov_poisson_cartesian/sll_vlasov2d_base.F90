#define MPI_MASTER 0
module sll_m_vlasov2d_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"
  use sll_m_cartesian_meshes
  use sll_m_collective
  use sll_m_remapper
  use sll_m_xml_io
  use sll_m_init_functions, only: VA_VALIS, PSTD, METH_BSL_CUBIC_SPLINES, SPECTRAL

  implicit none

  public :: initialize_vlasov2d_base, free_vlasov2d_base
  public :: compute_charge, compute_current
  public :: write_energy
  public :: read_input_file
  public :: transposexv, transposevx, write_xmf_file

  type, public :: vlasov2d_base

     logical                                  :: transposed      
     type(sll_cartesian_mesh_1d), pointer     :: geomx
     type(sll_cartesian_mesh_1d), pointer     :: geomv
     sll_real64, dimension(:,:),  pointer     :: f
     sll_real64, dimension(:,:),  pointer     :: ft
     sll_real64, dimension(:),    pointer     :: ex
     sll_real64, dimension(:),    pointer     :: jx
     sll_real64, dimension(:),    pointer     :: rho
     sll_real64                               :: dt     
     sll_int32                                :: nbiter 
     sll_int32                                :: fdiag  
     sll_int32                                :: fthdiag
     sll_int32                                :: nc_eta1
     sll_int32                                :: nc_eta2
     sll_int32                                :: np_eta1
     sll_int32                                :: np_eta2
     sll_real64                               :: eta1_min
     sll_real64                               :: eta2_min
     sll_real64                               :: eta1_max
     sll_real64                               :: eta2_max
     sll_real64                               :: delta_eta1
     sll_real64                               :: delta_eta2
     sll_int32                                :: va 
     sll_int32                                :: num_case 
     sll_real64                               :: eps

     type(layout_2d), pointer                 :: layout_x
     type(layout_2d), pointer                 :: layout_v
     type(remap_plan_2d_real64), pointer      :: x_to_v 
     type(remap_plan_2d_real64), pointer      :: v_to_x

  end type vlasov2d_base

  sll_int32, public  :: poisson_type 
  sll_int32, public  :: maxwell_type

  private

  sll_int32  :: i, j
  sll_int32  :: ierr
  sll_int32  :: ithf  !< file unit to store energy time evolution

contains

  subroutine read_input_file(this)

    class(vlasov2d_base),intent(inout)    :: this

    sll_int32  :: idata            !< file unit for namelist

    sll_int32  :: nx               ! dimensions de l'espace physique
    sll_int32  :: nvx              ! dimensions de l'espace des vitesses
    sll_real64 :: x0               ! coordonnees debut du maillage espace physique
    sll_real64 :: vx0              ! coordonnees debut du maillage espace vitesses
    sll_real64 :: x1               ! coordonnees fin du maillage espace physique
    sll_real64 :: vx1              ! coordonnees fin du maillage espace vitesses
    sll_real64 :: dt               ! time step
    sll_int32  :: nbiter           ! number of loops over time
    sll_int32  :: fdiag            ! diagnostics frequency
    sll_int32  :: fthdiag          ! time history frequency
    sll_int32  :: va               ! algo charge type
    sll_int32  :: num_case         ! test case
    sll_real64 :: eps = 0.05_f64   ! perturbation amplitude
    sll_int32  :: meth             ! method
    sll_int32  :: psize
    sll_int32  :: prank
    sll_int32  :: comm

    namelist /time/          dt, nbiter
    namelist /diag/          fdiag, fthdiag
    namelist /phys_space/    x0,x1,nx
    namelist /vel_space/     vx0,vx1,nvx
    namelist /test_case/     num_case, eps
    namelist /algo_charge/   va, meth
    namelist /field_solvers/ poisson_type, maxwell_type

    va           = VA_VALIS
    meth         = METH_BSL_CUBIC_SPLINES
    poisson_type = SPECTRAL
    maxwell_type = PSTD

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    if (prank == MPI_MASTER) then

       call initialize_file(idata, ithf)
       read(idata,NML=time)
       read(idata,NML=diag)
       read(idata,NML=phys_space)
       read(idata,NML=vel_space)
       read(idata,NML=test_case)
       read(idata,NML=algo_charge)
       read(idata,NML=field_solvers)
       close(idata)
       write(*,*) 'physical space: nx, x0, x1, dx'
       write(*,"((i3,1x),3(g13.3,1x))")nx, x0, x1, this%delta_eta1
       write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
       write(*,"((i3,1x),3(g13.3,1x))")nvx, vx0, vx1, this%delta_eta2
       write(*,*) 'dt,nbiter,fdiag,fthdiag'
       write(*,"(g13.3,1x,3i5)") dt,nbiter,fdiag,fthdiag
       write(*,*) " Algo charge "
       select case(va) 
       case(0)
          print*, 'Valis' 
       case(1)
          print*, 'Vlasov-Poisson'
       case(2)
          print*, " diag charge "
       case(3)
          print*, "classic algorithm"
       end select

    end if

    call mpi_bcast(dt,           1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(nbiter,       1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(fdiag,        1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(fthdiag,      1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(x0,           1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(x1,           1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(nx,           1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(vx0,          1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(vx1,          1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(nvx,          1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(va,           1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(meth,         1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(num_case,     1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(eps,          1, MPI_REAL8  , MPI_MASTER, comm, ierr)
    call mpi_bcast(poisson_type, 1, MPI_INTEGER, MPI_MASTER, comm, ierr)
    call mpi_bcast(maxwell_type, 1, MPI_INTEGER, MPI_MASTER, comm, ierr)

    this%dt         = dt
    this%nbiter     = nbiter
    this%fdiag      = fdiag
    this%fthdiag    = fthdiag

    this%geomx      => new_cartesian_mesh_1d(nx,x0,x1)
    this%geomv      => new_cartesian_mesh_1d(nvx,vx0,vx1)

    this%nc_eta1    = this%geomx%num_cells
    this%nc_eta2    = this%geomv%num_cells

    this%np_eta1    = this%geomx%num_cells+1
    this%np_eta2    = this%geomv%num_cells+1

    this%eta1_min   = this%geomx%eta_min
    this%eta2_min   = this%geomv%eta_min

    this%eta1_max   = this%geomx%eta_max
    this%eta2_max   = this%geomv%eta_max

    this%delta_eta1 = this%geomx%delta_eta
    this%delta_eta2 = this%geomv%delta_eta

    this%va         = va
    this%num_case   = num_case
    this%eps        = eps


  end subroutine read_input_file

  subroutine initialize_vlasov2d_base(this)

    use sll_m_hdf5_io_serial

    class(vlasov2d_base),intent(inout)    :: this
    sll_int32                             :: error
    sll_int32                             :: file_id
    sll_real64, dimension(:), allocatable :: eta1
    sll_real64, dimension(:), allocatable :: eta2
    sll_int32                             :: prank
    sll_int32                             :: psize
    sll_int32                             :: comm
    sll_int32                             :: loc_sz_i
    sll_int32                             :: loc_sz_j

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    error = 0

    this%transposed = .false.

    this%layout_x => new_layout_2d( sll_world_collective )        
    call initialize_layout_with_distributed_array( &
         this%np_eta1, this%np_eta2, 1,int(psize,4),this%layout_x)

    if ( prank == MPI_MASTER ) call sll_view_lims( this%layout_x )
    call flush(6)

    call compute_local_sizes(this%layout_x, loc_sz_i,loc_sz_j)        
    SLL_CLEAR_ALLOCATE(this%f(1:loc_sz_i,1:loc_sz_j),ierr)


    this%layout_v => new_layout_2d( sll_world_collective )
    call initialize_layout_with_distributed_array( &
         this%np_eta1, this%np_eta2,int(psize,4),1,this%layout_v)

    if ( prank == MPI_MASTER ) call sll_view_lims( this%layout_v )
    call flush(6)

    call compute_local_sizes(this%layout_v, loc_sz_i,loc_sz_j)        
    SLL_CLEAR_ALLOCATE(this%ft(1:loc_sz_i,1:loc_sz_j),ierr)

    this%x_to_v => new_remap_plan( this%layout_x, this%layout_v, this%f)     
    this%v_to_x => new_remap_plan( this%layout_v, this%layout_x, this%ft)  

    if (prank == MPI_MASTER) then

       SLL_ALLOCATE(eta1(this%np_eta1),error)
       SLL_ALLOCATE(eta2(this%np_eta2),error)

       do i = 1, this%np_eta1
          eta1(i) = this%eta1_min + (i-1)*this%delta_eta1
       end do

       do j = 1, this%np_eta2
          eta2(j) = this%eta2_min + (j-1)*this%delta_eta2
       end do

       call sll_hdf5_file_create("mesh2d.h5",file_id,error)
       call sll_hdf5_write_array(file_id,eta1,"/x1",error)
       call sll_hdf5_write_array(file_id,eta2,"/x2",error)
       call sll_hdf5_file_close(file_id, error)
    end if

    nullify(this%ex)
    nullify(this%jx)
    nullify(this%rho)

  end subroutine initialize_vlasov2d_base

  subroutine free_vlasov2d_base(this)

    class(vlasov2d_base),intent(inout) :: this

    SLL_DEALLOCATE_ARRAY(this%f, ierr)

  end subroutine free_vlasov2d_base

  subroutine compute_charge(this)
    class(vlasov2d_base),intent(inout) :: this

    sll_int32                             :: comm
    sll_real64, dimension(this%np_eta2+1) :: locrho
    sll_int32                             :: loc_sz_i
    sll_int32                             :: loc_sz_j
    sll_int32                             :: gi
    sll_int32                             :: gj
    sll_int32                             :: global_indices(2)
    sll_int32                             :: error

    call compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j)        

    locrho = 0.0_f64
    do i=1,loc_sz_i
       global_indices = local_to_global(this%layout_v,(/i,1/)) 
       gi = global_indices(1)
       gj = global_indices(2)
       locrho(gi) = sum(this%ft(i,:))*this%delta_eta2 
    end do
    this%rho = 0.0_f64
    comm  = sll_world_collective%comm
    call mpi_barrier(comm,error)
    call mpi_allreduce(locrho,this%rho,this%np_eta1,MPI_REAL8,MPI_SUM,comm,error)

  end subroutine compute_charge

  subroutine compute_current(this)

    class(vlasov2d_base),intent(inout) :: this

    sll_int32                             :: comm
    sll_real64, dimension(this%np_eta2+1) :: locjx
    sll_int32                             :: loc_sz_i
    sll_int32                             :: loc_sz_j
    sll_int32                             :: gi
    sll_int32                             :: global_indices(2)
    sll_int32                             :: error
    sll_int32                             :: nc_x
    sll_real64                            :: v
    sll_real64                            :: v_min
    sll_real64                            :: delta_v

    nc_x    = this%nc_eta1
    v_min   = this%eta2_min
    delta_v = this%delta_eta2

    call compute_local_sizes(this%layout_v,loc_sz_i,loc_sz_j)        

    locjx = 0.0_f64
    do i=1,loc_sz_i
       global_indices = local_to_global(this%layout_v,(/i,1/)) 
       gi = global_indices(1)
       do j = 1, loc_sz_j
          v  = v_min +(j-1)*delta_v
          locjx(gi) = locjx(gi)+this%ft(i,j) * delta_v * v
       end do
    end do
    this%jx = 0.0_f64
    comm  = sll_world_collective%comm
    call mpi_barrier(comm,error)
    call mpi_allreduce(locjx,this%jx,this%np_eta1,MPI_REAL8,MPI_SUM,comm,error)

  end subroutine compute_current

  subroutine thdiag(this,nrj,t)

    class(vlasov2d_base),intent(inout) :: this
    sll_real64, intent(in)             :: t
    sll_real64, intent(in)             :: nrj  
    sll_real64,dimension(11)           :: auxloc
    sll_real64,dimension(13)           :: aux
    sll_real64                         :: cell_volume

    aux    = 0.0
    auxloc = 0.0
    cell_volume = this%delta_eta1 * this%delta_eta2
    auxloc(1) = cell_volume * sum(this%f) ! avg(f)
    auxloc(2) = cell_volume * sum(abs(this%f)) ! L1 norm
    auxloc(3) = cell_volume * sum(this%f*this%f) ! L2 norm

    aux(13)=t
    aux(12)=nrj
    write(*,"('time ', g12.3,' test nrj', 4f10.5)") t, nrj, aux(1:3)
    call time_history(ithf, "thf","(13(1x,e15.6))",aux(1:13))

  end subroutine thdiag

  subroutine transposexv(this)

    class(vlasov2d_base),intent(inout) :: this

    SLL_ASSERT(.not. this%transposed)
    call apply_remap_2D( this%x_to_v, this%f, this%ft )
    this%transposed = .true.

  end subroutine transposexv

  subroutine transposevx(this)

    class(vlasov2d_base),intent(inout) :: this

    SLL_ASSERT(this%transposed)
    call apply_remap_2D( this%v_to_x, this%ft, this%f )
    this%transposed = .false.

  end subroutine transposevx

  subroutine write_xmf_file(this, iplot)

    use hdf5
    use sll_m_hdf5_io_serial

    class(vlasov2d_base),intent(in) :: this
    sll_int32, intent(in)           :: iplot
    sll_int32                       :: error
    character(len=4)                :: cplot
    sll_int32                       :: prank
    sll_int32, parameter            :: one = 1
    sll_int32                       :: file_id
    sll_int32                       :: nx1, nx2

    call int2string(iplot,cplot)
    call write_fx1x2(this,cplot)

    prank = sll_get_collective_rank(sll_world_collective)
    if (prank == MPI_MASTER) then

       nx1 = this%np_eta1
       nx2 = this%np_eta2

       call sll_xml_file_create("fvalues_"//cplot//".xmf",file_id,error)
       call write_grid(file_id,nx1,nx2,"x1","x2",cplot)
       write(file_id,"(a)")"</Domain>"
       write(file_id,"(a)")"</Xdmf>"
       close(file_id)

    endif

  end subroutine write_xmf_file

  subroutine write_grid(file_id,nx,ny,xname,yname,cplot)

    sll_int32 :: file_id, nx, ny
    character(len=*) :: cplot, xname, yname

    write(file_id,"(a)")"<Grid Name='"//xname//yname//"' GridType='Uniform'>"
    write(file_id, &
         "(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='",ny,nx,"'/>"
    write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
    write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx, &
         "' NumberType='Float' Precision='8' Format='HDF'>"
    write(file_id,"(a)")"mesh2d.h5:/"//xname
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a,i5,a)")"<DataItem Dimensions='",ny, &
         "' NumberType='Float' Precision='8' Format='HDF'>"
    write(file_id,"(a)")"mesh2d.h5:/"//yname
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a)")"</Geometry>"
    call write_attribute(file_id,nx,ny,"df",cplot,xname,yname)
    write(file_id,"(a)")"</Grid>"

  end subroutine write_grid

  subroutine write_attribute(file_id,nx,ny,fname,cplot,xname,yname)

    sll_int32                    :: file_id
    sll_int32                    :: nx
    sll_int32                    :: ny
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: cplot
    character(len=*), optional   :: xname
    character(len=*), optional   :: yname

    write(file_id,"(a)") &
         "<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
    write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ny,nx, &
         "' NumberType='Float' Precision='8' Format='HDF'>"
    if( present(xname) .and. present(yname)) then
       write(file_id,"(a)")fname//xname//yname//"_"//cplot//".h5:/values"
    else
       write(file_id,"(a)")fname//"_"//cplot//".h5:/values"
    end if
    write(file_id,"(a)")"</DataItem>"
    write(file_id,"(a)")"</Attribute>"

  end subroutine write_attribute

  subroutine write_fx1x2(this,cplot)

    use sll_m_hdf5_io_serial
    class(vlasov2d_base),intent(in)     :: this
    character(len=*)                    :: cplot
    sll_int32                           :: error
    sll_int32                           :: file_id
    sll_int32                           :: prank

    prank = sll_get_collective_rank(sll_world_collective)
    if (prank == MPI_MASTER) then
       call sll_hdf5_file_create('dfx1x2_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%f,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if

  end subroutine write_fx1x2

  subroutine write_energy(this, time)

    class(vlasov2d_base) :: this
    sll_real64 :: time, nrj

    nrj=sum(this%ex*this%ex)*this%delta_eta1
    call thdiag(this,nrj,time)

  end subroutine write_energy

end module sll_m_vlasov2d_base
