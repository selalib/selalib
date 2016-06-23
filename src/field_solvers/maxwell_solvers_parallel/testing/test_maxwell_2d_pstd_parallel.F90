!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
#define MPI_MASTER 0
program test_maxwell_2d_periodic_cart_par

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use hdf5, only: &
    hsize_t, hid_t

  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_hdf5_io_parallel, only: &
    sll_o_hdf5_file_close, &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

  use sll_m_maxwell_2d_periodic_cartesian_par, only: &
    sll_s_ampere_te, &
    sll_s_delete_maxwell_2d_periodic_plan_cartesian_par, &
    sll_s_faraday_te, &
    sll_t_maxwell_2d_periodic_plan_cartesian_par, &
    sll_f_new_maxwell_2d_periodic_plan_cartesian_par

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_i_min, &
    sll_o_get_layout_j_min, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_o_view_lims

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_m_xml_io, only: &
    sll_s_xml_file_close, &
    sll_s_xml_file_create

  use sll_mpi, only: &
    mpi_real8, &
    mpi_reduce, &
    mpi_sum, &
    mpi_wtime

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32   :: ncx
  sll_int32   :: ncy
  sll_int32   :: nx_loc
  sll_int32   :: ny_loc
  sll_int32   :: error
  sll_int32   :: i
  sll_int32   :: j
  sll_int32   :: gi
  sll_int32   :: gj
  sll_int32   :: prank
  sll_int64   :: psize 
  sll_int32   :: nprocx
  sll_int32   :: nprocy
  sll_int32   :: e
  sll_int32   :: istep
  sll_int32   :: nstep
  sll_int32   :: mode = 1
  sll_real32  :: ok 
  sll_real64  :: dt
  sll_real64  :: time = 0.0_f64
  sll_real64  :: omega
  sll_real64  :: err_l2
  sll_real64  :: err_glob
  sll_real64  :: Lx
  sll_real64  :: Ly
  sll_real64  :: dx
  sll_real64  :: dy
  sll_real64  :: x
  sll_real64  :: y
  sll_real64  :: tcpu1
  sll_real64  :: tcpu2

  type (sll_t_maxwell_2d_periodic_plan_cartesian_par), pointer :: plan

  type(sll_t_layout_2d), pointer            :: layout_x
  type(sll_t_layout_2d), pointer            :: layout_y
  sll_real64, dimension(:,:), pointer :: ex
  sll_real64, dimension(:,:), pointer :: ey
  sll_real64, dimension(:,:), pointer :: hz
  sll_int32, dimension(2)             :: global
  character(len=4)                    :: cstep

  ok = 1.0

  !Boot parallel environment
  call sll_s_boot_collective()

  ! Number of cells is equal to number of points in this case
  ncx = 64
  ncy = 64
  Lx  = 2.0*sll_p_pi
  Ly  = 2.0*sll_p_pi

  psize = int(sll_f_get_collective_size(sll_v_world_collective),i64)
  prank = sll_f_get_collective_rank(sll_v_world_collective)

  dx = Lx/ncx
  dy = Ly/ncy

  omega  = sqrt((mode*sll_p_pi/Lx)**2+(mode*sll_p_pi/Ly)**2)

  e      = int(log(real(psize))/log(2.))
  print *, 'running on ', 2**e, 'processes'

  ! Layout and local sizes for FFTs in x-direction
  layout_x => sll_f_new_layout_2d( sll_v_world_collective )
  layout_y => sll_f_new_layout_2d( sll_v_world_collective )
  nprocx = 2**e
  nprocy = 2**e

  call sll_o_initialize_layout_with_distributed_array( ncx, &
                                                    ncy, &
                                                      1, &
                                                 nprocy, &
                                                layout_x )

  call sll_o_initialize_layout_with_distributed_array( ncx, &
                                                    ncy, &
                                                 nprocx, &
                                                      1, &
                                                layout_y )
  flush( output_unit )
  call sll_o_view_lims(layout_x)
  flush( output_unit )
  call sll_o_view_lims(layout_y)
  flush( output_unit )

  plan => sll_f_new_maxwell_2d_periodic_plan_cartesian_par(layout_x, &
                                                     layout_y, &
                                                     ncx, ncy, &
                                                     Lx, Ly)
  plan%e_0  = 1.0_f64 
  plan%mu_0 = 1.0_f64 

  hz => plan%fz_x

  !Ex si sequential along y
  call sll_o_compute_local_sizes(plan%layout_y,nx_loc,ny_loc)
  SLL_CLEAR_ALLOCATE(ex(1:nx_loc,1:ny_loc), error)

  !Ey si sequential along x
  call sll_o_compute_local_sizes(plan%layout_x,nx_loc,ny_loc)
  SLL_CLEAR_ALLOCATE(ey(1:nx_loc,1:ny_loc), error)

  do j=1,ny_loc
     do i=1,nx_loc
        global = sll_o_local_to_global( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        hz(i,j) = - cos(mode*x) * cos(mode*y) * cos(omega*time)
     end do
  end do

  call write_fields_2d('fields-0000')

  tcpu1 = MPI_WTIME()

  nstep = 1000
  dt = .5 / sqrt (1./(dx*dx)+1./(dy*dy))
  print*, ' dt = ', dt
  call sll_s_faraday_te(plan,0.5*dt,ex,ey)
  time = 0.5*dt

  do istep = 1, nstep

     call sll_s_int2string(istep,cstep)
     call sll_s_ampere_te(plan,dt,ex,ey)
     time = time + dt
     call sll_s_faraday_te(plan,dt,ex,ey)
     call write_fields_2d('fields-'//cstep)

     err_l2 = 0.0_f64
     do j=1,ny_loc
     do i=1,nx_loc
        global = sll_o_local_to_global( layout_x, (/i, j/))
        gi = global(1)
        gj = global(2)
        x  = (gi-1)*dx
        y  = (gj-1)*dy
        err_l2 = err_l2 + abs(hz(i,j) + cos(mode*x)*cos(mode*y)*cos(omega*time))
     end do
     end do

     call mpi_reduce(err_l2,err_glob,1,mpi_real8,MPI_SUM,0, &
                     sll_v_world_collective%comm, error)

     if ( prank == MPI_MASTER ) then
        write(*,"(10x,' istep = ',I6)",advance="no") istep
        write(*,"(' time = ',g12.3,' sec')",advance="no") time
        write(*,"(' L2 error = ',g15.5)") err_glob / (ncx*ncy)
     end if

  end do
     
  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) then
     write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize
  end if

  SLL_DEALLOCATE_ARRAY(ex, error)
  SLL_DEALLOCATE_ARRAY(ey, error)

  call sll_s_delete_maxwell_2d_periodic_plan_cartesian_par(plan)

  call sll_s_halt_collective()

contains

  !> Experimental interface for an hdf5 array writer in parallel
  subroutine write_fields_2d( file_name )

    character(len=*), intent(in)     :: file_name
    integer(HSIZE_T), dimension(1:2) :: global_dims
    sll_int32                        :: error
    sll_int32                        :: file_id
    integer(HID_T)                   :: hfile_id
    integer(HSIZE_T), dimension(1:2) :: offset

    global_dims(:) = (/ int(ncx,HSIZE_T),int(ncy,HSIZE_T) /)

    call sll_o_hdf5_file_create(file_name//'.h5',sll_v_world_collective%comm, &
                              hfile_id,error)
    offset(1) = int(sll_o_get_layout_i_min(layout_y,prank) - 1,HSIZE_T)
    offset(2) = int(sll_o_get_layout_j_min(layout_y,prank) - 1,HSIZE_T)
    call sll_o_hdf5_write_array(hfile_id,global_dims,offset, &
                              ex,'ex_values',error)
    offset(1) = int(sll_o_get_layout_i_min(layout_x,prank) - 1,HSIZE_T)
    offset(2) = int(sll_o_get_layout_j_min(layout_x,prank) - 1,HSIZE_T)
    call sll_o_hdf5_write_array(hfile_id,global_dims,offset, &
                              ey,'ey_values',error)
    call sll_o_hdf5_write_array(hfile_id,global_dims,offset, &
                              hz,'hz_values',error)
    call sll_o_hdf5_file_close(hfile_id,error)

    if (prank == MPI_MASTER) then

       call sll_s_xml_file_create(file_name//".xmf",file_id,error)
       write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
       write(file_id,"(a,2i5,a)") &
          "<Topology TopologyType='2DCoRectMesh' NumberOfElements='", &
          ncy,ncx,"'/>"
       write(file_id,"(a)")  &
          "<Geometry GeometryType='ORIGIN_DXDY'>"
       write(file_id,"(a)")  &
          "<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
       write(file_id,"(2f12.5)") 0., 0.
       write(file_id,"(a)")"</DataItem>"
       write(file_id,"(a)")  &
          "<DataItem Dimensions='2' NumberType='Float' Format='XML'>"
       write(file_id,"(2f12.5)") dx, dy
       write(file_id,"(a)")"</DataItem>"
       write(file_id,"(a)")"</Geometry>"
       write(file_id,"(a)")  &
          "<Attribute Name='Ex' AttributeType='Scalar' Center='Node'>"
       write(file_id,"(a,2i5,a)") &
          "<DataItem Dimensions='",ncy,ncx, &
          "' NumberType='Float' Precision='8' Format='HDF'>"
       write(file_id,"(a)") file_name//".h5:/ex_values"
       write(file_id,"(a)")"</DataItem>"
       write(file_id,"(a)")"</Attribute>"
       write(file_id,"(a)")  &
          "<Attribute Name='Ey' AttributeType='Scalar' Center='Node'>"
       write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ncy,ncx, &
          "' NumberType='Float' Precision='8' Format='HDF'>"
       write(file_id,"(a)") file_name//".h5:/ey_values"
       write(file_id,"(a)")"</DataItem>"
       write(file_id,"(a)")"</Attribute>"
       write(file_id,"(a)")  &
          "<Attribute Name='Hz' AttributeType='Scalar' Center='Node'>"
       write(file_id,"(a,2i5,a)")"<DataItem Dimensions='",ncy,ncx, &
          "' NumberType='Float' Precision='8' Format='HDF'>"
       write(file_id,"(a)") file_name//".h5:/hz_values"
       write(file_id,"(a)")"</DataItem>"
       write(file_id,"(a)")"</Attribute>"
       call sll_s_xml_file_close(file_id,error)

    end if

  end subroutine write_fields_2d

end program test_maxwell_2d_periodic_cart_par
