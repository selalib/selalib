#define MPI_MASTER 0
module sll_vlasov2d_base

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_cartesian_meshes.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
#include "sll_utilities.h"

use sll_remapper
use sll_xml_io
use init_functions
use sll_collective
use sll_remapper

implicit none

public :: initialize_vlasov2d_base, free_vlasov2d_base
public :: compute_charge, compute_current
public :: write_energy
public :: read_input_file

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
 
end subroutine read_input_file

subroutine initialize_vlasov2d_base(this)

  use sll_hdf5_io_serial
 
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
  sll_real64                         :: dvx

  dvx = this%delta_eta2

  this%rho = sum(this%f(:,:))*dvx 

end subroutine compute_charge

subroutine compute_current(this)

  class(vlasov2d_base),intent(inout) :: this
  sll_real64                         :: vx
  sll_real64                         :: dvx
  sll_int32                          :: i
  sll_int32                          :: j

  dvx = this%delta_eta2
  this%jx = 0.0_f64
  do j = 1, this%np_eta2
  do i = 1, this%np_eta1
     vx = this%eta2_min+(j-1)*this%delta_eta2
     this%jx(i) = this%jx(i) + dvx * this%f(i,j) * vx
  end do
  end do

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

!subroutine write_fx1x2(this,cplot)
!
!  use sll_hdf5_io_serial
!  class(vlasov2d_base),intent(in)     :: this
!  character(len=*)                    :: cplot
!  sll_int32                           :: error
!  sll_int32                           :: file_id
!  
!  call sll_hdf5_file_create('fx1x2_'//cplot//".h5",file_id,error)
!  call sll_hdf5_write_array(file_id,this%f,"/values",error)
!  call sll_hdf5_file_close(file_id, error)
!
!end subroutine write_fx1x2

subroutine write_energy(this, time)

  class(vlasov2d_base) :: this
  sll_real64 :: time, nrj
  
  nrj=sum(this%ex*this%ex)*this%delta_eta1
  call thdiag(this,nrj,time)
  
end subroutine write_energy

end module sll_vlasov2d_base
