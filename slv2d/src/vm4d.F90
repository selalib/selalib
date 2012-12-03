program vm4d

#include "selalib.h"
  use used_precision  
  use geometry_module
  use diagnostiques_module
  use sll_vlasov4d
  use sll_cubic_spline_interpolator_1d
  use maxwell2dfdtd_module
  use poisson2d_periodic
  use remapper

  implicit none

  type(geometry)   :: geomx 
  type(geometry)   :: geomv 
  type(vlasov4d)   :: vlas4d 
  type (maxwell2dfdtd) :: maxw2dfdtd 
  type(poisson2dpp):: poisson 

  type(cubic_spline_1d_interpolator), target :: spl_x1
  type(cubic_spline_1d_interpolator), target :: spl_x2
  type(cubic_spline_1d_interpolator), target :: spl_x3
  type(cubic_spline_1d_interpolator), target :: spl_x4

  sll_real64, dimension(:,:,:,:), pointer :: f4d
  sll_real64, dimension(:,:),     pointer :: bz
  sll_real64, dimension(:,:),     pointer :: ex
  sll_real64, dimension(:,:),     pointer :: ey 
  sll_real64, dimension(:,:),     pointer :: jx
  sll_real64, dimension(:,:),     pointer :: jy 
  sll_real64, dimension(:,:),     pointer :: rho 

  sll_int32  :: nbiter, iter 
  sll_int32  :: fdiag, fthdiag  
  sll_int32  :: jstartx, jendx, jstartv, jendv

  sll_real64 :: dt     
  sll_real64 :: nrj, tcpu1, tcpu2

  sll_int32  :: prank, comm
  sll_int64  :: psize

  sll_real64, dimension(:,:,:,:), allocatable :: local_array1
  sll_real64, dimension(:,:,:,:), allocatable :: local_array2
  sll_int32                                   :: loc_sz_i_init
  sll_int32                                   :: loc_sz_j_init
  sll_int32                                   :: loc_sz_k_init
  sll_int32                                   :: loc_sz_l_init
  sll_int32                                   :: loc_sz_i_final
  sll_int32                                   :: loc_sz_j_final
  sll_int32                                   :: loc_sz_k_final
  sll_int32                                   :: loc_sz_l_final
  sll_int32                                   :: npi
  sll_int32                                   :: npj
  sll_int32                                   :: npk
  sll_int32                                   :: npl
  sll_int32                                 :: gi, gj, gk, gl
  sll_int32                                   :: ierr
  type(layout_4D), pointer                  :: layout1
  type(layout_4D), pointer                  :: layout2
  type(remap_plan_4D_real64), pointer       :: rmp4
  sll_int32, dimension(4)                   :: global_indices, g


  call sll_boot_collective()
  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm   = sll_world_collective%comm

  tcpu1 = MPI_WTIME()
  if (.not. is_power_of_two(psize)) then     
     print *, 'This test needs to run in a number of processes which is ',&
          'a power of 2.'
     stop
  end if
  if (prank == MPI_MASTER) then
     print*,'MPI Version of slv2d running on ',psize, ' processors'
  end if

  call initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)

  if (prank == MPI_MASTER) then
     ! write some run data
     write(*,*) 'physical space: nx, ny, x0, x1, y0, y1, dx, dy'
     write(*,"(2(i3,1x),6(g13.3,1x))") geomx%nx, geomx%ny, geomx%x0, &
          geomx%x0+(geomx%nx)*geomx%dx, &
          geomx%y0, geomx%y0+(geomx%ny)*geomx%dy, &
          geomx%dx, geomx%dy   
     write(*,*) 'velocity space: nvx, nvy, vx0, vx1, vy0, vy1, dvx, dvy'
     write(*,"(2(i3,1x),6(g13.3,1x))") geomv%nx, geomv%ny, geomv%x0, &
          geomv%x0+(geomv%nx-1)*geomv%dx, &
          geomv%y0, geomv%y0+(geomv%ny-1)*geomv%dy, &
          geomv%dx, geomv%dy
     write(*,*) 'dt,nbiter,fdiag,fthdiag'
     write(*,"(g13.3,1x,3i3)") dt,nbiter,fdiag,fthdiag
  endif

  layout1  => new_layout_4D( sll_world_collective )        
  call factorize_in_random_2powers_4d(psize, npi, npj, npk, npl)
  if( prank .eq. 0 ) then
     print *, 'source configuration: ', npi, npj, npk, npl
  end if

  call initialize_layout_with_distributed_4D_array( &
          geomx%nx, &
          geomx%ny, &
          geomv%nx, &
          geomv%ny, &
          npi, &
          npj, &
          npk, &
          npl, &
          layout1 )

  call compute_local_sizes_4d(layout1,       &
                              loc_sz_i_init, &
                              loc_sz_j_init, &
                              loc_sz_k_init, &
                              loc_sz_l_init )        

  SLL_ALLOCATE(local_array1(loc_sz_i_init,loc_sz_j_init,loc_sz_k_init,loc_sz_l_init),ierr)
 
  layout2  => new_layout_4D( sll_world_collective )
  call factorize_in_random_2powers_4d(psize, npi, npj, npk, npl)

  if( prank .eq. 0 ) then
     print *, 'target configuration: ', npi, npj, npk, npl
  end if

  call initialize_layout_with_distributed_4D_array( &
          geomx%nx, &
          geomx%ny, &
          geomv%nx, &
          geomv%ny, &
          npi, &
          npj, &
          npk, &
          npl, &
          layout2 )

  call compute_local_sizes_4d(layout2, &
                              loc_sz_i_final, &
                              loc_sz_j_final, &
                              loc_sz_k_final, &
                              loc_sz_l_final )

  SLL_ALLOCATE(local_array2(loc_sz_i_final,loc_sz_j_final,loc_sz_k_final,loc_sz_l_final), ierr )
   
  rmp4 => new_remap_plan( layout1, layout2, local_array1)     

  call apply_remap_4D( rmp4, local_array1, local_array2 ) 

  print *, prank, 'Printing layout1: '
  call sll_view_lims_4D( layout1 )
  print *, prank, 'Printing layout2: '
  call sll_view_lims_4D( layout2 )

  call delete_layout_4D( layout1 )
  call delete_layout_4D( layout2 )
  SLL_DEALLOCATE_ARRAY(local_array1, ierr)
  SLL_DEALLOCATE_ARRAY(local_array2, ierr)



      


  call initlocal(jstartv,jendv,jstartx,jendx)

  call plot_mesh4d(geomx,geomv,jstartx,jendx,jstartv,jendv)

  call transposevx(vlas4d,f4d)
  call advection_x1(vlas4d,f4d,0.5*dt)
  call advection_x2(vlas4d,f4d,0.5*dt)

  do iter=1,nbiter

     if (mod(iter,fdiag) == 0) then 
        call plot_df(f4d,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, XY_VIEW)
        call plot_df(f4d,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, XVX_VIEW)
        call plot_df(f4d,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, YVY_VIEW)
        call plot_df(f4d,iter/fdiag,geomx,geomv,1,geomx%ny,jstartv,jendv, VXVY_VIEW)
     end if

     call transposexv(vlas4d,f4d)

     call densite_charge(vlas4d,rho)
     call densite_courant(vlas4d,jx,jy)
     
     !call solve(poisson,ex,ey,rho,nrj)
     if (iter == 1) then
     !   call solve_faraday(maxw2dfdtd,ex,ey,bz,0.5_f64*dt)
         call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,0.5_f64*dt)
     else
     !   call solve_faraday(maxw2dfdtd,ex,ey,bz,dt)
         call solve_ampere(maxw2dfdtd,ex,ey,bz,jx,jy,nrj,dt)
     end if
     
     call cl_periodiques(maxw2dfdtd,ex,ey,bz,jx,jy,dt)

     call advection_x3(vlas4d,ex,dt)
     call advection_x4(vlas4d,ey,dt)

     call transposevx(vlas4d,f4d)

     call advection_x1(vlas4d,f4d,dt)
     call advection_x2(vlas4d,f4d,dt)

     if (mod(iter,fthdiag).eq.0) then 
      call thdiag(vlas4d,f4d,nrj,iter*dt,jstartv)
     endif

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) &
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize

  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains

  subroutine initglobal(geomx,geomv,dt,nbiter,fdiag,fthdiag)

    type(geometry)  :: geomx       ! geometrie globale du probleme
    type(geometry)  :: geomv       ! geometrie globale du probleme
    sll_real64      :: dt          ! pas de temps
    sll_int32       :: nbiter      ! nombre d'iterations en temps
    sll_int32       :: fdiag       ! frequences des diagnostiques
    sll_int32       :: fthdiag     ! frequences des historiques en temps
    sll_int32       :: nx, ny      ! dimensions de l'espace physique
    sll_int32       :: nvx, nvy    ! dimensions de l'espace des vitesses
    sll_real64      :: x0, y0      ! coordonnees debut du maillage espace physique
    sll_real64      :: vx0, vy0    ! coordonnees debut du maillage espace vitesses
    sll_real64      :: x1, y1      ! coordonnees fin du maillage espace physique
    sll_real64      :: vx1, vy1    ! coordonnees fin du maillage espace vitesses
    sll_int32       :: iflag,ierr  ! indicateur d'erreur

    namelist /time/ dt, nbiter
    namelist /diag/ fdiag, fthdiag
    namelist /phys_space/ x0,x1,y0,y1,nx,ny
    namelist /vel_space/ vx0,vx1,vy0,vy1,nvx,nvy
    sll_int32 :: prank, psize, comm

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm   = sll_world_collective%comm

    if (prank == MPI_MASTER) then

       call fichinit()
       read(idata,NML=time)
       read(idata,NML=diag)
       read(idata,NML=phys_space)
       read(idata,NML=vel_space)

    end if

    call mpi_bcast(dt,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(nbiter,  1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(fdiag,   1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(fthdiag, 1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(x0,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(y0,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(x1,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(y1,      1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(nx,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(ny,      1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(vx0,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(vy0,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(vx1,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(vy1,     1,MPI_REAL8,MPI_MASTER,comm,ierr)
    call mpi_bcast(nvx,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)
    call mpi_bcast(nvy,     1,MPI_INTEGER ,MPI_MASTER,comm,ierr)

    call new(geomx,x0,y0,x1,y1,nx,ny,iflag,"perxy")
    call new(geomv,vx0,vy0,vx1,vy1,nvx,nvy,iflag,"natxy")

  end subroutine initglobal

  subroutine initlocal(jstartv,jendv,jstartx,jendx )

    sll_int32      :: jstartv
    sll_int32      :: jendv
    sll_int32      :: jstartx
    sll_int32      :: jendx

    sll_int32  :: ipiece_size_v
    sll_int32  :: ipiece_sizex

    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,iv,jv,iflag
    sll_int32  :: prank, psize, comm
    sll_real64 :: xi, eps, kx, ky



    prank      = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm        = sll_world_collective%comm

    ipiece_size_v = (geomv%ny + psize-1) / psize
    ipiece_sizex = (geomx%ny + psize-1) / psize

    jstartv = prank * ipiece_size_v + 1
    jendv   = min(jstartv - 1 + ipiece_size_v, geomv%ny)
    jstartx = prank * ipiece_sizex + 1
    jendx   = min(jstartx - 1 + ipiece_sizex, geomx%ny)

    SLL_ASSERT(jstartx<jendx)
    SLL_ASSERT(jstartv<jendv)
    print*,'init zone ',prank,jstartx,jendx,ipiece_sizex, &
         jstartv,jendv,ipiece_size_v

    SLL_ALLOCATE(f4d(geomx%nx,geomx%ny,geomv%nx,jstartv:jendv),iflag)

    SLL_ALLOCATE(ex(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(ey(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(jx(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(jy(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(bz(geomx%nx,geomx%ny),iflag)
    SLL_ALLOCATE(rho(geomx%nx,geomx%ny),iflag)

    xi  = 0.90_f64
    eps = 0.05_f64
    kx  = 2_f64*sll_pi/((geomx%nx)*geomx%dx)
    ky  = 2_f64*sll_pi/((geomx%ny)*geomx%dy)
    do jv=jstartv,jendv
       vy = geomv%y0+(jv-1)*geomv%dy
       do iv=1,geomv%nx
          vx = geomv%x0+(iv-1)*geomv%dx
          v2 = vx*vx+vy*vy
          do j=1,geomx%ny
             y=geomx%y0+(j-1)*geomx%dy
             do i=1,geomx%nx
                x=geomx%x0+(i-1)*geomx%dx
                f4d(i,j,iv,jv)=(1+eps*cos(kx*x))*1/(2*sll_pi)*exp(-.5*v2)
             end do
          end do
       end do
    end do

    call spl_x1%initialize(geomx%nx, geomx%x0, geomx%x1, PERIODIC_SPLINE)
    call spl_x2%initialize(geomx%ny, geomx%y0, geomx%y1, PERIODIC_SPLINE)
    call spl_x3%initialize(geomv%nx, geomv%x0, geomv%x1, PERIODIC_SPLINE)
    call spl_x4%initialize(geomv%ny, geomv%y0, geomv%y1, PERIODIC_SPLINE)
    call new(vlas4d,geomx,geomv,spl_x1,spl_x2,spl_x3,spl_x4,jstartx,jendx,jstartv,jendv,iflag)

    call new(poisson, ex, ey, geomx, iflag)
    call transposexv(vlas4d,f4d)
    call densite_charge(vlas4d,rho)
    call solve(poisson,ex,ey,rho,nrj)
    call new(maxw2dfdtd,geomx,iflag, jstartx, jendx)

!    call free(poisson)

  end subroutine initlocal

  subroutine factorize_in_random_2powers_4d(n, n1, n2, n3, n4)
    sll_int64, intent(in) :: n
    sll_int32, intent(out) ::n1, n2, n3, n4
    sll_int32   :: expo, expo1, expo2, expo3, expo4
    sll_real32 ::  rand_real
    if (.not.is_power_of_two(n)) then   
       print*, 'The number of processors must be a power of 2'
       stop
    endif 
    expo = int(log(real(n))/log(2.))  
    call random_number(rand_real)
    expo1 = int(rand_real*expo)
    call random_number(rand_real)
    expo2 = int(rand_real*(expo-expo1))
    call random_number(rand_real)
    expo3 = int(rand_real*(expo-(expo1+expo2)))
    expo4 = expo - (expo1 + expo2 + expo3)
    n1 = 2**expo1
    n2 = 2**expo2
    n3 = 2**expo3
    n4 = 2**expo4
  end subroutine factorize_in_random_2powers_4d

end program vm4d
