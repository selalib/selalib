program vp4d_multigrid

#include "selalib-mpi.h"

  use sll_vlasov4d_base
  use sll_vlasov4d_poisson
  use init_functions
  use sll_multigrid_2d

  implicit none

  type(vlasov4d_poisson)    :: vlasov 
  type(multigrid_2d)        :: poisson

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
  type(cubic_spline_1d_interpolator), target :: spl_x3
  type(cubic_spline_1d_interpolator), target :: spl_x4

  type(layout_2D), pointer                :: layout_mg

  sll_real64, dimension(:,:), allocatable :: phi_local
  sll_real64, dimension(:,:), allocatable :: rho_local
  sll_real64, dimension(:,:), allocatable :: den
  sll_real64, dimension(:,:), allocatable :: ex_local
  sll_real64, dimension(:,:), allocatable :: ey_local
  sll_real64, dimension(:,:), allocatable :: ex_global
  sll_real64, dimension(:,:), allocatable :: ey_global

  sll_int32  :: iter = 0
  sll_int32  :: prank, comm, comm2d
  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
  sll_int64  :: psize
  sll_real64 :: tcpu1, tcpu2
  sll_int32  :: i,j,k,l,error
  sll_real64 :: vx,vy,v2,x,y
  sll_real64 :: kx, ky
  sll_int32  :: gi, gj, gk, gl
  sll_int32  :: global_indices(4)
  sll_real64 :: time
  sll_int32  :: nxprocs, nyprocs
  sll_int32  :: sx, ex, sy, ey
  sll_int32  :: block_sz

  sll_int32, dimension(8)     :: neighbor
  sll_int32, parameter        :: N =7, S =3, W =5, E =1
  sll_int32, parameter        :: NW=6, SW=4, NE=8, SE=2
  sll_int32, parameter        :: ndims = 2
  sll_int32, dimension(ndims) :: dims
  sll_int32, dimension(ndims) :: coords
  logical                     :: reorder
  logical,dimension(ndims)    :: periods

  sll_int32                   :: block_t
  integer                     :: statut(MPI_STATUS_SIZE)

  call sll_boot_collective()

  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  tcpu1 = MPI_WTIME()

  if (prank == MPI_MASTER) then
     print*,'MPI Version of slv2d running on ',psize, ' processors'
  end if

  call read_input_file(vlasov)

  call spl_x1%initialize(vlasov%nc_eta1+1,  &
                         vlasov%eta1_min,   &
                         vlasov%eta1_max,   &
                         SLL_PERIODIC)

  call spl_x2%initialize(vlasov%nc_eta2+1,  &
                         vlasov%eta2_min,   &
                         vlasov%eta2_max,   &
                         SLL_PERIODIC)

  call spl_x3%initialize(vlasov%nc_eta3+1,  &
                         vlasov%eta3_min,   &
                         vlasov%eta3_max,   &
                         SLL_PERIODIC)

  call spl_x4%initialize(vlasov%nc_eta4+1,  &
                         vlasov%eta4_min,   &
                         vlasov%eta4_max,   &
                         SLL_PERIODIC)

  call initialize(vlasov,spl_x1,spl_x2,spl_x3,spl_x4,error)

  call compute_local_sizes_4d(vlasov%layout_x, &
         loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

  kx  = 2_f64*sll_pi/(vlasov%eta1_max-vlasov%eta1_min)
  ky  = 2_f64*sll_pi/(vlasov%eta2_max-vlasov%eta2_min)
    
  do l=1,loc_sz_l 
  do k=1,loc_sz_k
  do j=1,loc_sz_j
  do i=1,loc_sz_i

     global_indices = local_to_global_4D(vlasov%layout_x,(/i,j,k,l/)) 
     gi = global_indices(1)
     gj = global_indices(2)
     gk = global_indices(3)
     gl = global_indices(4)

     x  = vlasov%eta1_min+(gi-1)*vlasov%delta_eta1
     y  = vlasov%eta2_min+(gj-1)*vlasov%delta_eta2
     vx = vlasov%eta3_min+(gk-1)*vlasov%delta_eta3
     vy = vlasov%eta4_min+(gl-1)*vlasov%delta_eta4

     v2 = vx*vx+vy*vy

     vlasov%f(i,j,k,l) = landau_cos_prod(vlasov%eps,kx, ky, x, y, v2)

  end do
  end do
  end do
  end do


  dims = 0
  CALL MPI_DIMS_CREATE(int(psize,4),ndims,dims,error)
  nxprocs = dims(1)
  nyprocs = dims(2)

  periods(1) = .true.
  periods(2) = .true.
  reorder    = .true.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods,reorder,comm2d,error)

  neighbor(:) = MPI_PROC_NULL

  CALL MPI_CART_SHIFT(comm2d,0,1,neighbor(N),neighbor(S),error)
  CALL MPI_CART_SHIFT(comm2d,1,1,neighbor(W),neighbor(E),error)
  CALL MPI_COMM_RANK(comm2d,prank,error)
  CALL MPI_CART_COORDS(comm2d,prank,ndims,coords,error)

  if (neighbor(N) /= MPI_PROC_NULL) then
     CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)+1/),neighbor(NW),error)
     CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)+1/),neighbor(NE),error)
  end if

  if (neighbor(S) /= MPI_PROC_NULL) then
     CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)-1/),neighbor(SW),error)
     CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)-1/),neighbor(SE),error)
  end if

  call flush(6)
  print"('proc: ',i2,', coords: ',2i3,', neighbors: ',8i3)",prank,coords,neighbor
  call flush(6)
  call MPI_BARRIER(comm,error)

  layout_mg => new_layout_2D( sll_world_collective )        

  call initialize_layout_with_distributed_2D_array( vlasov%nc_eta1, &
                                                    vlasov%nc_eta2, &
                                                    nxprocs,          &
                                                    nyprocs,          &
                                                    layout_mg)
  
  if ( prank == MPI_MASTER ) then
    print*, "Display layout for multigrid"
    call sll_view_lims_2D( layout_mg )
    call write_to_file(layout_mg, 'mesh')
  end if
  call flush(6)

  call compute_local_sizes_2d(layout_mg, loc_sz_i,loc_sz_j)        

  sx = get_layout_2D_i_min(layout_mg, prank)
  ex = get_layout_2D_i_max(layout_mg, prank)
  sy = get_layout_2D_j_min(layout_mg, prank)
  ey = get_layout_2D_j_max(layout_mg, prank)
  
  call initialize( poisson, sx, ex, sy, ey,          &
                   vlasov%nc_eta1*vlasov%delta_eta1, &
                   vlasov%nc_eta2*vlasov%delta_eta2, &
                   nxprocs,                          &
                   nyprocs,                          &
                   vlasov%nc_eta1,                   &
                   vlasov%nc_eta2,                   &
                   comm2d,                           &
                   neighbor  )

  SLL_CLEAR_ALLOCATE(phi_local(sx-1:ex+1,sy-1:ey+1), error)
  SLL_CLEAR_ALLOCATE(rho_local(sx-1:ex+1,sy-1:ey+1), error)
  SLL_CLEAR_ALLOCATE(den(sx-1:ex+1,sy-1:ey+1), error); den = 1.0_f64
  SLL_CLEAR_ALLOCATE(ex_local(sx-1:ex+1,sy-1:ey+1), error)
  SLL_CLEAR_ALLOCATE(ey_local(sx-1:ex+1,sy-1:ey+1), error)

  call transposexv(vlasov)
  call compute_charge(vlasov)
  call transposevx(vlasov)

  call global_plot(vlasov%rho, "rho")

  do j = sy, ey
     do i = sx, ex
        x = vlasov%eta1_min + (i-1) * vlasov%delta_eta1
        y = vlasov%eta2_min + (j-1) * vlasov%delta_eta2
        rho_local(i,j) = vlasov%rho(i,j)
     end do
  end do

  call global_scale( rho_local )

  call local_plot(rho_local, "rhs")

  call solve(poisson, phi_local, rho_local, den)
  !call global_scale( phi )

  call local_plot(phi_local, "phi")
  call local_plot(rho_local, "res")

  !compute Ex and Ey
  do j = sy,ey
     do i = sx,ex
        ex_local(i,j) = 0.5*(phi_local(i+1,j)-phi_local(i-1,j))/vlasov%delta_eta1
        ey_local(i,j) = 0.5*(phi_local(i,j+1)-phi_local(i,j-1))/vlasov%delta_eta2
     end do
  end do

  call local_plot(ex_local, "ex_local")
  call local_plot(ex_local, "ey_local")

  SLL_CLEAR_ALLOCATE(ex_global(1:vlasov%nc_eta1,1:vlasov%nc_eta2), error)
  SLL_CLEAR_ALLOCATE(ey_global(1:vlasov%nc_eta1,1:vlasov%nc_eta2), error)

  call MPI_TYPE_CREATE_SUBARRAY(2, &
                               (/vlasov%nc_eta1,vlasov%nc_eta2/), &
                               (/ex-sx+1,ey-sy+1/), &
                               (/sx, sy/), &
                               MPI_ORDER_FORTRAN, &
                               MPI_REAL8, block_t, error)

  block_sz = (ex-sx+1)*(ey-sy+1)
  call MPI_ALLGATHER(ex_local,block_sz,MPI_REAL8, &
                     ex_global,block_sz,MPI_REAL8, &
                     MPI_SUM,comm2d,error)
  call MPI_ALLGATHER(ey_local,block_sz,MPI_REAL8, &
                     ey_global,block_sz,MPI_REAL8, &
                     MPI_SUM,comm2d,error)

  vlasov%ex(1:vlasov%nc_eta1,1:vlasov%nc_eta2) = ex_global
  vlasov%ey(1:vlasov%nc_eta1,1:vlasov%nc_eta2) = ey_global

  vlasov%ex(vlasov%np_eta1,:) = vlasov%ex(1,:) 
  vlasov%ex(:,vlasov%np_eta2) = vlasov%ex(:,1) 
  vlasov%ey(vlasov%np_eta1,:) = vlasov%ey(1,:) 
  vlasov%ey(:,vlasov%np_eta2) = vlasov%ey(:,1) 

  call global_plot(vlasov%ex, 'ex')
  call global_plot(vlasov%ey, 'ey')

  time = 0.0_f64
  call advection_x1(vlasov,0.5*vlasov%dt)
  call advection_x2(vlasov,0.5*vlasov%dt)

  do iter=1, vlasov%nbiter

     if (iter == 1 .or. mod(iter, vlasov%fdiag) == 0) then 
        call write_xmf_file(vlasov,iter/ vlasov%fdiag)
     end if

     call transposexv(vlasov)

     call compute_charge(vlasov)

     time = time + 0.5*vlasov%dt
     if (mod(iter, vlasov%fthdiag).eq.0) then 
        call write_energy(vlasov,time)
     end if

     call advection_x3(vlasov, vlasov%dt)
     call advection_x4(vlasov, vlasov%dt)

     call transposevx(vlasov)

     time = time + 0.5*vlasov%dt

     call advection_x1(vlasov, vlasov%dt)
     call advection_x2(vlasov, vlasov%dt)


  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call sll_halt_collective()

contains

   subroutine local_plot(field, fieldname)

      character(len=*) :: fieldname
      sll_real64, intent(in) :: field(sx-1:ex+1,sy-1:ey+1)

      call sll_gnuplot_rect_2d_parallel((sx-1)*vlasov%delta_eta1, &
                                        vlasov%delta_eta1,        &
                                        (sy-1)*vlasov%delta_eta2, &
                                        vlasov%delta_eta2,        &
                                        ex-sx+1,                  &
                                        ey-sy+1,                  &
                                        field(sx:ex,sy:ey),       &
                                        fieldname,                &
                                        iter,                     &
                                        error)  

   end subroutine local_plot

   subroutine global_plot(field, fieldname)

      character(len=*) :: fieldname
      sll_real64, intent(in) :: field(:,:)
   
      if (prank == MPI_MASTER) &
      call sll_gnuplot_2d(vlasov%eta1_min, &
                          vlasov%eta1_max, &
                          vlasov%np_eta1,  &
                          vlasov%eta2_min, &
                          vlasov%eta2_max, &
                          vlasov%np_eta2,  &
                          field,           &
                          fieldname,       &
                          iter,            &
                          error)  

   end subroutine global_plot

   subroutine global_scale( field )

   sll_real64 :: field(sx-1:ex+1,sy-1:ey+1)
   sll_real64 :: avloc, av

   avloc = sum(field(sx:ex,sy:ey))

   call MPI_ALLREDUCE(avloc,av,1,MPI_REAL8,MPI_SUM,comm2d,error)

   av=av/float(vlasov%nc_eta1*vlasov%nc_eta2)

   field=field-av

   end subroutine global_scale


end program vp4d_multigrid
