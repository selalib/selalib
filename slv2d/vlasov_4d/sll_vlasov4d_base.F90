module sll_vlasov4d_base

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_logical_meshes.h"
#include "sll_constants.h"

 use sll_module_interpolators_1d_base
 use sll_module_interpolators_2d_base
 use sll_remapper
 use sll_xml_io

 implicit none
 
 public :: initialize_vlasov4d_base, free_vlasov4d_base
 public :: compute_charge, compute_current
 public :: transposexv, transposevx
 public :: write_xmf_file, write_energy
 public :: read_input_file

 type, public :: vlasov4d_base
   logical                                  :: transposed      
   sll_real64, dimension(:,:,:,:),  pointer :: f
   sll_real64, dimension(:,:,:,:),  pointer :: ft
   type(layout_4D), pointer                 :: layout_x
   type(layout_4D), pointer                 :: layout_v
   type(remap_plan_4D_real64), pointer      :: x_to_v 
   type(remap_plan_4D_real64), pointer      :: v_to_x
   sll_real64, dimension(:,:), pointer      :: ex
   sll_real64, dimension(:,:), pointer      :: ey
   sll_real64, dimension(:,:), pointer      :: jx
   sll_real64, dimension(:,:), pointer      :: jy
   sll_real64, dimension(:,:), pointer      :: bz
   sll_real64, dimension(:,:), pointer      :: rho
   type(sll_logical_mesh_2d), pointer       :: geomx
   type(sll_logical_mesh_2d), pointer       :: geomv
   sll_real64                               :: dt     
   sll_int32                                :: nbiter 
   sll_int32                                :: fdiag  
   sll_int32                                :: fthdiag
   sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
   sll_real64 :: eta1_min, eta2_min, eta3_min, eta4_min
   sll_real64 :: eta1_max, eta2_max, eta3_max, eta4_max
   sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_eta4
   sll_int32  :: va 
   sll_int32  :: num_case 
   sll_real64 :: eps
 end type vlasov4d_base


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
  sll_int32  :: idata !< file unit for namelist

  sll_int32  :: nx, ny           ! dimensions de l'espace physique
  sll_int32  :: nvx, nvy         ! dimensions de l'espace des vitesses
  sll_real64 :: x0, y0           ! coordonnees debut du maillage espace physique
  sll_real64 :: vx0, vy0         ! coordonnees debut du maillage espace vitesses
  sll_real64 :: x1, y1           ! coordonnees fin du maillage espace physique
  sll_real64 :: vx1, vy1         ! coordonnees fin du maillage espace vitesses
  sll_real64 :: dt               ! time step
  sll_int32  :: nbiter           ! number of loops over time
  sll_int32  :: fdiag            ! diagnostics frequency
  sll_int32  :: fthdiag          ! time history frequency
  sll_int32  :: va = 0           ! algo charge type
  sll_int32  :: num_case         ! test case
  sll_real64 :: eps = 0.05_f64   ! perturbation amplitude
  sll_int32  :: meth = 0         ! method

  namelist /time/ dt, nbiter
  namelist /diag/ fdiag, fthdiag
  namelist /phys_space/ x0,x1,y0,y1,nx,ny
  namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
  namelist /test_case/ num_case, eps
  namelist /algo_charge/ va, meth

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

  end if

  call mpi_bcast(dt,      1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(nbiter,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fdiag,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(fthdiag, 1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(x0,      1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(y0,      1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(x1,      1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(y1,      1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(nx,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(ny,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(vx0,     1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(vy0,     1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(vx1,     1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(vy1,     1,MPI_REAL8   ,MPI_MASTER,comm,ierr)
  call mpi_bcast(nvx,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(nvy,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(va,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(meth,    1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(num_case,1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
  call mpi_bcast(eps,     1,MPI_REAL8   ,MPI_MASTER,comm,ierr)

  this%dt         = dt
  this%nbiter     = nbiter
  this%fdiag      = fdiag
  this%fthdiag    = fthdiag

  this%geomx      => new_logical_mesh_2d(nx,ny,       &
                                         eta1_min=x0, &
                                         eta1_max=x1, & 
                                         eta2_min=y0, &
                                         eta2_max=y1)

  this%geomv      => new_logical_mesh_2d(nvx,nvy,vx0,vx1,vy0,vy1)

  this%nc_eta1    = this%geomx%num_cells1
  this%nc_eta2    = this%geomx%num_cells2
  this%nc_eta3    = this%geomv%num_cells1
  this%nc_eta4    = this%geomv%num_cells2

  this%eta1_min   = this%geomx%eta1_min
  this%eta2_min   = this%geomx%eta2_min
  this%eta3_min   = this%geomv%eta1_min
  this%eta4_min   = this%geomv%eta2_min

  this%eta1_max   = this%geomx%eta1_max
  this%eta2_max   = this%geomx%eta2_max
  this%eta3_max   = this%geomv%eta1_max
  this%eta4_max   = this%geomv%eta2_max

  this%delta_eta1 = this%geomx%delta_eta1
  this%delta_eta2 = this%geomx%delta_eta2
  this%delta_eta3 = this%geomv%delta_eta1
  this%delta_eta4 = this%geomv%delta_eta2

  this%va         = va
  this%num_case   = num_case
  this%eps        = eps

  if (prank == MPI_MASTER) then

       write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
       write(*,"(2(i3,1x),6(g13.3,1x))") &
        this%nc_eta1, this%nc_eta2, this%eta1_min, this%eta1_max, &
        this%eta2_min, this%eta2_max, this%delta_eta1, this%delta_eta2
       write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
       write(*,"(2(i3,1x),6(g13.3,1x))") &
          this%nc_eta3, this%nc_eta4, this%eta3_min, this%eta3_max, &
          this%eta4_min, this%eta4_max, this%delta_eta3, this%delta_eta4
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

  if (.not. is_power_of_two(int(psize,i64))) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'greater than 4 and a power of 2.'
     call sll_halt_collective()
     stop
  end if

 end subroutine read_input_file

 subroutine initialize_vlasov4d_base(this)

  use sll_hdf5_io

  class(vlasov4d_base),intent(inout)    :: this
  sll_int32                             :: error
  sll_int32                             :: psize
  sll_int32                             :: prank
  sll_int32                             :: comm
  sll_int32                             :: file_id
  sll_real64, dimension(:), allocatable :: eta1
  sll_real64, dimension(:), allocatable :: eta2
  sll_real64, dimension(:), allocatable :: eta3
  sll_real64, dimension(:), allocatable :: eta4

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  error = 0

  this%transposed = .false.

  this%layout_x => new_layout_4D( sll_world_collective )        
  call initialize_layout_with_distributed_4D_array( &
             this%nc_eta1, this%nc_eta2, this%nc_eta3, this%nc_eta4,    &
             1,1,1,int(psize,4),this%layout_x)

  if ( prank == MPI_MASTER ) call sll_view_lims_4D( this%layout_x )
  call flush(6)

  call compute_local_sizes_4d(this%layout_x, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%f(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

  this%layout_v => new_layout_4D( sll_world_collective )
  call initialize_layout_with_distributed_4D_array( &
              this%nc_eta1, this%nc_eta2, this%nc_eta3, this%nc_eta4,    &
              1,int(psize,4),1,1,this%layout_v)

  if ( prank == MPI_MASTER ) call sll_view_lims_4D( this%layout_v )
  call flush(6)

  call compute_local_sizes_4d(this%layout_v, &
                              loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
  SLL_CLEAR_ALLOCATE(this%ft(1:loc_sz_i,1:loc_sz_j,1:loc_sz_k,1:loc_sz_l),ierr)

  this%x_to_v => new_remap_plan( this%layout_x, this%layout_v, this%f)     
  this%v_to_x => new_remap_plan( this%layout_v, this%layout_x, this%ft)     
  
  if(prank == MPI_MASTER) then

     SLL_ALLOCATE(eta1(this%nc_eta1),error)
     SLL_ALLOCATE(eta2(this%nc_eta2),error)
     SLL_ALLOCATE(eta3(this%nc_eta3),error)
     SLL_ALLOCATE(eta4(this%nc_eta4),error)

     do i = 1, this%nc_eta1
        eta1(i) = this%eta1_min + (i-1)*this%delta_eta1
     end do

     do j = 1, this%nc_eta2
        eta2(j) = this%eta2_min + (j-1)*this%delta_eta2
     end do

     do k = 1, this%nc_eta3
        eta3(k) = this%eta1_min + (k-1)*this%delta_eta1
     end do

     do l = 1, this%nc_eta4
        eta4(l) = this%eta2_min + (l-1)*this%delta_eta2
     end do

     print *,'Printing layout x: '
     call sll_view_lims_4D(this%layout_x)
     print *,'Printing layout v: '
     call sll_view_lims_4D(this%layout_v)

     call sll_hdf5_file_create("mesh4d.h5",file_id,error)
     call sll_hdf5_write_array(file_id,eta1,"/x1",error)
     call sll_hdf5_write_array(file_id,eta2,"/x2",error)
     call sll_hdf5_write_array(file_id,eta3,"/x3",error)
     call sll_hdf5_write_array(file_id,eta4,"/x4",error)
     call sll_hdf5_file_close(file_id, error)

  end if

  nullify(this%ex)
  nullify(this%ey)
  nullify(this%bz)
  nullify(this%jx)
  nullify(this%jy)
  nullify(this%rho)

 end subroutine initialize_vlasov4d_base

 subroutine free_vlasov4d_base(this)

  class(vlasov4d_base),intent(inout) :: this

  call delete_layout_4D(this%layout_x)
  call delete_layout_4D(this%layout_v)
  SLL_DEALLOCATE_ARRAY(this%f, ierr)
  SLL_DEALLOCATE_ARRAY(this%ft, ierr)

 end subroutine free_vlasov4d_base

 subroutine compute_charge(this)

   class(vlasov4d_base),intent(inout) :: this
   sll_int32                          :: error
   sll_int32                          :: c
   sll_int32                          :: comm
   sll_real64                         :: dvxvy

   sll_real64, dimension(this%geomx%num_cells1,this%geomx%num_cells2) :: locrho

   dvxvy = this%geomv%delta_eta1*this%geomv%delta_eta2

   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v, &
        loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
   
   locrho(:,:) = 0.
   do j=1,loc_sz_j
   do i=1,loc_sz_i
      global_indices = local_to_global_4D(this%layout_v,(/i,j,1,1/)) 
      gi = global_indices(1)
      gj = global_indices(2)
      locrho(gi,gj) = sum(this%ft(i,j,:,:))*dvxvy 
   end do
   end do
   this%rho(:,:) = 0.
   comm  = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%num_cells1*this%geomx%num_cells2
   call mpi_allreduce(locrho,this%rho,c,MPI_REAL8,MPI_SUM,comm,error)

 end subroutine compute_charge

 subroutine compute_current(this)

   class(vlasov4d_base),intent(inout)                 :: this
   sll_int32                                          :: error
   sll_real64                                         :: vx
   sll_real64                                         :: vy 
   sll_real64, dimension(this%geomx%num_cells1,this%geomx%num_cells2) :: locjx
   sll_real64, dimension(this%geomx%num_cells1,this%geomx%num_cells2) :: locjy
   sll_int32                                          :: c
   sll_int32                                          :: comm
   sll_real64                                         :: dvxvy

   dvxvy = this%geomv%delta_eta1*this%geomv%delta_eta2
   SLL_ASSERT(this%transposed)
   call compute_local_sizes_4d(this%layout_v,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

   locjx(:,:) = 0.; locjy(:,:) = 0.
   do l=1,loc_sz_l
   do k=1,loc_sz_k
   do j=1,loc_sz_j
   do i=1,loc_sz_i
      global_indices = local_to_global_4D(this%layout_v,(/i,j,k,l/)) 
      gi = global_indices(1)
      gj = global_indices(2)
      gk = global_indices(3)
      gl = global_indices(4)
      vx = this%geomv%eta1_min+(gk-1)*this%geomv%delta_eta1
      vy = this%geomv%eta2_min+(gl-1)*this%geomv%delta_eta2
      locjx(gi,gj) = locjx(gi,gj) + dvxvy*this%ft(i,j,k,l) * vx
      locjy(gi,gj) = locjy(gi,gj) + dvxvy*this%ft(i,j,k,l) * vy
   end do
   end do
   end do
   end do

   this%jx(:,:) = 0.; this%jy(:,:) = 0.
   comm   = sll_world_collective%comm
   call mpi_barrier(comm,error)
   c=this%geomx%num_cells1*this%geomx%num_cells2
   call mpi_allreduce(locjx,this%jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,this%jy,c, MPI_REAL8,MPI_SUM,comm,error)

 end subroutine compute_current

 subroutine thdiag(this,nrj,t)

   class(vlasov4d_base),intent(inout) :: this
   sll_real64, intent(in) :: t,nrj  
   sll_real64,dimension(11) :: auxloc
   sll_int32 :: prank, psize
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: comm, error
   sll_real64 :: cell_volume

   comm  = sll_world_collective%comm
   prank = sll_get_collective_rank(sll_world_collective)
   psize = sll_get_collective_size(sll_world_collective)

   auxloc = 0.0
   cell_volume = this%geomx%delta_eta1 * this%geomx%delta_eta2 &
               * this%geomv%delta_eta1 * this%geomv%delta_eta2
   auxloc(1) = cell_volume * sum(this%f) ! avg(f)
   auxloc(2) = cell_volume * sum(abs(this%f)) ! L1 norm
   auxloc(3) = cell_volume * sum(this%f*this%f) ! L2 norm
   
   call mpi_reduce(auxloc,aux(1:11),11,MPI_REAL8,MPI_SUM,MPI_MASTER,comm, error)

   if (prank == MPI_MASTER) then
      diag=0.
!      aux=0.
      aux(13)=t
      aux(12)=nrj
      write(*,"('time ', g12.3,' test nrj',f10.5, f10.5)") t, nrj, aux(1)
      call time_history(ithf, "thf","(13(1x,e15.6))",aux(1:13))
   end if

 end subroutine thdiag

 subroutine transposexv(this)

   class(vlasov4d_base),intent(inout) :: this

   SLL_ASSERT(.not. this%transposed)
   call apply_remap_4D( this%x_to_v, this%f, this%ft )
   this%transposed = .true.

 end subroutine transposexv

 subroutine transposevx(this)

   class(vlasov4d_base),intent(inout) :: this

   SLL_ASSERT(this%transposed)
   call apply_remap_4D( this%v_to_x, this%ft, this%f )
   this%transposed = .false.

 end subroutine transposevx

 subroutine write_xmf_file(this, iplot)

  use hdf5
  use sll_hdf5_io

  class(vlasov4d_base),intent(in) :: this
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

     nx1 = this%geomx%num_cells1
     nx2 = this%geomv%num_cells1
     nx3 = this%geomx%num_cells2
     nx4 = this%geomv%num_cells2

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

  write(file_id,"(a)")"<Attribute Name='"//fname//"' AttributeType='Scalar' Center='Node'>"
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
 use sll_hdf5_io
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
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
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
 use sll_hdf5_io
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
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
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
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_CLEAR_ALLOCATE(fjl(1:loc_sz_j,1:loc_sz_l),error)
 do l=1,loc_sz_l
    do j=1,loc_sz_j
       fjl(j,l) = sum(this%f(:,j,:,l))
    end do
 end do
 global_dims = (/this%geomx%num_cells2,this%geomv%num_cells2/)
 offset(1) = get_layout_4D_j_min(this%layout_x,prank)-1
 offset(2) = get_layout_4D_l_min(this%layout_x,prank)-1
 call sll_hdf5_file_create('fx2x4_'//cplot//".h5",pfile_id,error)
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
 call compute_local_sizes_4d(this%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        
 SLL_CLEAR_ALLOCATE(fkl(1:loc_sz_k,1:loc_sz_l),error)
 do l=1,loc_sz_l
    do k=1,loc_sz_k
       fkl(k,l) = sum(this%f(:,:,k,l))
    end do
 end do
 global_dims = (/this%geomv%num_cells1,this%geomv%num_cells2/)
 offset(1) = get_layout_4D_k_min(this%layout_x,prank)-1
 offset(2) = get_layout_4D_l_min(this%layout_x,prank)-1
 call sll_hdf5_file_create('fx3x4_'//cplot//".h5",pfile_id,error)
 call sll_hdf5_write_array_2d(pfile_id,global_dims,offset,fkl,"/values",error)
 call sll_hdf5_file_close(pfile_id, error)

 end subroutine write_fx3x4

 subroutine write_energy(this, time)
   class(vlasov4d_base) :: this
   sll_real64 :: time, nrj
   
   nrj=sum(this%ex*this%ex+this%ey*this%ey)*this%delta_eta1*this%delta_eta2
   nrj=0.5_f64*log(nrj)
   call thdiag(this,nrj,time)
   
 end subroutine write_energy

end module sll_vlasov4d_base
