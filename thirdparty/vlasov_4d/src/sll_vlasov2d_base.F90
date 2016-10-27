module sll_vlasov2d_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_cartesian_meshes
use sll_m_constants
use sll_m_remapper
use sll_m_xml_io
use init_functions

implicit none

public :: initialize_vlasov4d_base, free_vlasov4d_base
public :: compute_charge, compute_current
public :: write_energy
public :: read_input_file

type, public :: vlasov4d_base
  logical                                  :: transposed      
  type(sll_cartesian_mesh_1d), pointer     :: geomx
  type(sll_cartesian_mesh_1d), pointer     :: geomv
  sll_real64, dimension(:,:),  pointer     :: f
  sll_real64, dimension(:,:),  pointer     :: ft
  sll_real64, dimension(:,:), pointer      :: ex
  sll_real64, dimension(:,:), pointer      :: jx
  sll_real64, dimension(:,:), pointer      :: rho
  sll_real64                               :: dt     
  sll_int32                                :: nbiter 
  sll_int32                                :: fdiag  
  sll_int32                                :: fthdiag
  sll_int32                                :: nc_eta1
  sll_int32                                :: nc_eta1
  sll_int32                                :: np_eta2
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
end type vlasov4d_base

sll_int32, public  :: poisson_type 
sll_int32, public  :: maxwell_type

private

sll_int32  :: i, j, k, l
sll_int32  :: loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l
sll_int32  :: global_indices(4), gi, gj, gk, gl
sll_int32  :: ierr
sll_int32  :: ithf  !< file unit to store energy time evolution

contains

subroutine read_input_file(this)

 class(vlasov4d_base),intent(inout)    :: this
 sll_int32                             :: psize
 sll_int32                             :: prank
 sll_int32                             :: comm
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

 call initialize_file(idata, ithf)
 read(idata,NML=time)
 read(idata,NML=diag)
 read(idata,NML=phys_space)
 read(idata,NML=vel_space)
 read(idata,NML=test_case)
 read(idata,NML=algo_charge)
 read(idata,NML=field_solvers)
 close(idata)

 this%dt         = dt
 this%nbiter     = nbiter
 this%fdiag      = fdiag
 this%fthdiag    = fthdiag

 this%geomx      => new_cartesian_mesh_1d(nx,x0,x1)
 this%geomv      => new_cartesian_mesh_1d(nvx,vx0,vx1)

 this%nc_eta1    = this%geomx%num_cells1
 this%nc_eta2    = this%geomv%num_cells2

 this%np_eta1    = this%geomx%num_cells1+1
 this%np_eta2    = this%geomv%num_cells2+1

 this%eta1_min   = this%geomx%eta1_min
 this%eta2_min   = this%geomv%eta2_min

 this%eta1_max   = this%geomx%eta1_max
 this%eta2_max   = this%geomv%eta2_max

 this%delta_eta1 = this%geomx%delta_eta1
 this%delta_eta2 = this%geomv%delta_eta1

 this%va         = va
 this%num_case   = num_case
 this%eps        = eps

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

 endif

 !if (.not. is_power_of_two(int(psize,i64))) then     
 !   print *, 'This test needs to run in a number of processes which is ',&
 !        'a power of 2.'
 !   call sll_halt_collective()
 !   stop
 !end if

end subroutine read_input_file

subroutine initialize_vlasov2d_base(this)

 use sll_hdf5_io_serial

 class(vlasov4d_base),intent(inout)    :: this
 sll_int32                             :: error
 sll_int32                             :: file_id
 sll_real64, dimension(:), allocatable :: eta1
 sll_real64, dimension(:), allocatable :: eta2

 error = 0

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

end subroutine initialize_vlasov4d_base

subroutine free_vlasov4d_base(this)

 class(vlasov4d_base),intent(inout) :: this

 SLL_DEALLOCATE_ARRAY(this%f, ierr)

end subroutine free_vlasov4d_base

subroutine compute_charge(this)

  class(vlasov4d_base),intent(inout) :: this
  sll_int32                          :: error
  sll_real64                         :: dvx

  sll_real64, dimension(this%np_eta1,this%np_eta2) :: locrho

  dvx = this%delta_eta2

  this%rho = sum(this%f(:,:))*dvx 

end subroutine compute_charge

subroutine compute_current(this)

  class(vlasov4d_base),intent(inout)                  :: this
  sll_int32                                           :: error
  sll_real64                                          :: vx
  sll_real64                                          :: dvx

  jx = 0.0_f64
  do j = 1, this%np_eta2
  do i = 1, this%np_eta1
     vx = this%eta2_min+(j-1)*this%delta_eta2
     jx(i) = jx(i) + dvx * this%f(i,j) * vx
  end do

end subroutine compute_current

subroutine thdiag(this,nrj,t)

  class(vlasov4d_base),intent(inout) :: this
  sll_real64, intent(in) :: t,nrj  
  sll_real64,dimension(11) :: auxloc
  sll_real64,dimension(13) :: aux
  sll_int32 :: comm, error
  sll_real64 :: cell_volume

  aux    = 0.0
  auxloc = 0.0
  cell_volume = this%delta_eta1 * this%delta_eta2 &
              * this%delta_eta3 * this%delta_eta4
  auxloc(1) = cell_volume * sum(this%f) ! avg(f)
  auxloc(2) = cell_volume * sum(abs(this%f)) ! L1 norm
  auxloc(3) = cell_volume * sum(this%f*this%f) ! L2 norm
  
  aux(13)=t
  aux(12)=nrj
  write(*,"('time ', g12.3,' test nrj', 4f10.5)") t, nrj, aux(1:3)
  call time_history(ithf, "thf","(13(1x,e15.6))",aux(1:13))

end subroutine thdiag

subroutine write_xmf_file(this, iplot)

 use hdf5
 use sll_hdf5_io_serial

 class(vlasov2d_base),intent(in) :: this
 sll_int32, intent(in)           :: iplot
 sll_int32                       :: error
 character(len=4)                :: cplot
 sll_int32                       :: prank
 sll_int32, parameter            :: one = 1
 sll_int32                       :: file_id
 sll_int32                       :: nx1, nx2, nx3, nx4

 call int2string(iplot,cplot)
 call write_fx1x2(this,cplot)
 call write_fx1x3(this,cplot)
 call write_fx2x4(this,cplot)
 call write_fx3x4(this,cplot)

 prank = sll_get_collective_rank(sll_world_collective)
 if (prank == MPI_MASTER) then

    if(associated(this%ex)) then
       call sll_hdf5_file_create('ex_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%ex,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if
    if(associated(this%ey)) then
       call sll_hdf5_file_create('ey_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%ey,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if
    if(associated(this%rho)) then
       call sll_hdf5_file_create('rho_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%rho,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if
    if(associated(this%jx)) then
       call sll_hdf5_file_create('jx_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%jx,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if
    if(associated(this%jy)) then
       call sll_hdf5_file_create('jy_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%jy,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if
    if(associated(this%bz)) then
       call sll_hdf5_file_create('bz_'//cplot//".h5",file_id,error)
       call sll_hdf5_write_array(file_id,this%bz,"/values",error)
       call sll_hdf5_file_close(file_id, error)
    end if

    nx1 = this%np_eta1
    nx2 = this%np_eta2
    nx3 = this%np_eta3
    nx4 = this%np_eta4

    call sll_xml_file_create("fvalues_"//cplot//".xmf",file_id,error)
    call write_grid(this,file_id,nx1,nx2,"x1","x2",cplot)
    call write_grid(this,file_id,nx1,nx3,"x1","x3",cplot)
    call write_grid(this,file_id,nx2,nx4,"x2","x4",cplot)
    call write_grid(this,file_id,nx3,nx4,"x3","x4",cplot)
    write(file_id,"(a)")"</Domain>"
    write(file_id,"(a)")"</Xdmf>"
    close(file_id)

 endif

end subroutine write_xmf_file

subroutine write_grid(this,file_id,nx,ny,xname,yname,cplot)

 class(vlasov4d_base),intent(in) :: this
 sll_int32 :: file_id, nx, ny
 character(len=*) :: cplot, xname, yname

 write(file_id,"(a)")"<Grid Name='"//xname//yname//"' GridType='Uniform'>"
 write(file_id, &
  "(a,2i5,a)")"<Topology TopologyType='2DRectMesh' NumberOfElements='",ny,nx,"'/>"
 write(file_id,"(a)")"<Geometry GeometryType='VXVY'>"
 write(file_id,"(a,i5,a)")"<DataItem Dimensions='",nx, &
                          "' NumberType='Float' Precision='8' Format='HDF'>"
 write(file_id,"(a)")"mesh4d.h5:/"//xname
 write(file_id,"(a)")"</DataItem>"
 write(file_id,"(a,i5,a)")"<DataItem Dimensions='",ny, &
                          "' NumberType='Float' Precision='8' Format='HDF'>"
 write(file_id,"(a)")"mesh4d.h5:/"//yname
 write(file_id,"(a)")"</DataItem>"
 write(file_id,"(a)")"</Geometry>"
 if (xname == "x1" .and. yname == "x2") then
    if(associated(this%ex)) &
    call write_attribute(file_id,nx,ny,"ex",cplot)
    if(associated(this%ey)) &
    call write_attribute(file_id,nx,ny,"ey",cplot)
    if(associated(this%rho)) &
    call write_attribute(file_id,nx,ny,"rho",cplot)
    if(associated(this%bz)) &
    call write_attribute(file_id,nx,ny,"bz",cplot)
    if(associated(this%jx)) &
    call write_attribute(file_id,nx,ny,"jx",cplot)
    if(associated(this%jy)) &
    call write_attribute(file_id,nx,ny,"jy",cplot)
 end if
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
use sll_hdf5_io_serial
class(vlasov4d_base),intent(in)     :: this
character(len=*)                    :: cplot
sll_int32                           :: error
sll_int32                           :: file_id
sll_int32                           :: prank
sll_int32                           :: comm
sll_real64, dimension(:,:), pointer :: fij
sll_real64                          :: sumloc

prank = sll_get_collective_rank(sll_world_collective)
comm  = sll_world_collective%comm
call compute_local_sizes(this%layout_p,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
SLL_CLEAR_ALLOCATE(fij(1:loc_sz_i,1:loc_sz_j),error)
do j=1,loc_sz_j
   do i=1,loc_sz_i
      sumloc = sum(this%f(i,j,:,:))
      call mpi_reduce(sumloc,fij(i,j),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
   end do
end do
if (prank == MPI_MASTER) then
   call sll_hdf5_file_create('fx1x2_'//cplot//".h5",file_id,error)
   call sll_hdf5_write_array(file_id,fij,"/values",error)
   call sll_hdf5_file_close(file_id, error)
end if
end subroutine write_fx1x2

subroutine write_fx1x3(this,cplot)
use sll_hdf5_io_serial
class(vlasov4d_base),intent(in)     :: this
character(len=*)                    :: cplot
sll_int32                           :: error
sll_int32                           :: file_id
sll_real64, dimension(:,:), pointer :: fik
sll_int32                           :: prank
sll_int32                           :: comm
sll_real64                          :: sumloc

prank = sll_get_collective_rank(sll_world_collective)
comm  = sll_world_collective%comm
call compute_local_sizes(this%layout_p,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
SLL_CLEAR_ALLOCATE(fik(1:loc_sz_i,1:loc_sz_k),error)
do k=1,loc_sz_k
   do i=1,loc_sz_i
      sumloc= sum(this%f(i,:,k,:))
      call mpi_reduce(sumloc,fik(i,k),1,MPI_REAL8,MPI_SUM,MPI_MASTER,comm,error)
   end do
end do
if (prank == MPI_MASTER) then
   call sll_hdf5_file_create('fx1x3_'//cplot//".h5",file_id,error)
   call sll_hdf5_write_array(file_id,fik,"/values",error)
   call sll_hdf5_file_close(file_id, error)
end if
end subroutine write_fx1x3

subroutine write_fx2x4(this,cplot)
use hdf5
use sll_hdf5_io_parallel
class(vlasov4d_base),intent(in)     :: this
character(len=*)                    :: cplot
integer(HID_T)                      :: pfile_id
integer(HSSIZE_T)                   :: offset(2)
integer(HSIZE_T)                    :: global_dims(2)
sll_int32                           :: error
sll_int32                           :: prank
sll_real64, dimension(:,:), pointer :: fjl

prank = sll_get_collective_rank(sll_world_collective)
call compute_local_sizes(this%layout_p,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
SLL_CLEAR_ALLOCATE(fjl(1:loc_sz_j,1:loc_sz_l),error)
do l=1,loc_sz_l
   do j=1,loc_sz_j
      fjl(j,l) = sum(this%f(:,j,:,l))
   end do
end do
global_dims = (/this%nc_eta2,this%nc_eta4/)
offset(1) = get_layout_j_min(this%layout_p,prank)-1
offset(2) = get_layout_l_min(this%layout_p,prank)-1
call sll_hdf5_file_create('fx2x4_'//cplot//".h5",MPI_COMM_WORLD,pfile_id,error)
call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fjl,"/values",error)
call sll_hdf5_file_close(pfile_id, error)

end subroutine write_fx2x4

subroutine write_fx3x4(this,cplot)
use hdf5
use sll_hdf5_io_parallel
class(vlasov4d_base),intent(in)     :: this
character(len=*)                    :: cplot
integer(HID_T)                      :: pfile_id
integer(HSSIZE_T)                   :: offset(2)
integer(HSIZE_T)                    :: global_dims(2)
sll_int32                           :: error
sll_int32                           :: prank
sll_real64, dimension(:,:), pointer :: fkl

prank = sll_get_collective_rank(sll_world_collective)
call compute_local_sizes(this%layout_p,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
SLL_CLEAR_ALLOCATE(fkl(1:loc_sz_k,1:loc_sz_l),error)
do l=1,loc_sz_l
   do k=1,loc_sz_k
      fkl(k,l) = sum(this%f(:,:,k,l))
   end do
end do
global_dims = (/this%nc_eta3,this%nc_eta4/)
offset(1) = get_layout_k_min(this%layout_p,prank)-1
offset(2) = get_layout_l_min(this%layout_p,prank)-1
call sll_hdf5_file_create('fx3x4_'//cplot//".h5",MPI_COMM_WORLD,pfile_id,error)
call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fkl,"/values",error)
call sll_hdf5_file_close(pfile_id, error)

end subroutine write_fx3x4

subroutine write_energy(this, time)
  class(vlasov4d_base) :: this
  sll_real64 :: time, nrj
  
  nrj=sum(this%ex*this%ex+this%ey*this%ey)*this%delta_eta1*this%delta_eta2
  call thdiag(this,nrj,time)
  
end subroutine write_energy

end module sll_vlasov2d_base
