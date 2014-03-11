program vm4d_spectral

#define MPI_MASTER 0
#include "selalib-mpi.h"

  use geometry_module
  use diagnostiques_module
  use sll_vlasov4d_base
  use sll_vlasov4d_spectral

  implicit none

  type(geometry)            :: geomx 
  type(geometry)            :: geomv 
  type(maxwell_2d_pstd)     :: maxwell
  type(poisson_2d_periodic) :: poisson 
  type(vlasov4d_spectral)   :: vlasov4d 

  type(cubic_spline_2d_interpolator), target :: spl_x3x4

  sll_int32  :: nbiter, iter , fdiag, fthdiag  
  sll_real64 :: dt, nrj, tcpu1, tcpu2, mass0

  sll_int32  :: prank, comm
  sll_int64  :: psize

  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
  sll_int32  :: va

  va=0 !(va=0 --> Valis ; va=1 --> Vlasov-Poisson ; va=2 --> diag charge ; va=3 --> classic algorithm)

  call sll_boot_collective()
  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  tcpu1 = MPI_WTIME()
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
     write(*,"(g13.3,1x,3i5)") dt,nbiter,fdiag,fthdiag
  endif

  call initlocal()

  !f --> ft
  call transposexv(vlasov4d)

  call compute_charge(vlasov4d)
  call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
  vlasov4d%exn=vlasov4d%ex
  vlasov4d%eyn=vlasov4d%ey

  !ft --> f
  call transposevx(vlasov4d)

  mass0=sum(vlasov4d%rho(1:vlasov4d%geomx%nx,1:vlasov4d%geomx%ny))*vlasov4d%geomv%dx*vlasov4d%geomv%dy

  !###############
  !TIME LOOP
  !###############
  do iter=1,nbiter

     !print *,'iter',iter
     if (iter ==1 .or. mod(iter,fdiag) == 0) then 
        call write_xmf_file(vlasov4d,iter/fdiag)
     end if


     if ((va==0).or.(va==3)) then 
        !f --> ft, current (this%jx,this%jy), ft-->f
        call transposexv(vlasov4d)
        !compute this%jx, this%jy (zero average) at time tn
        call compute_current(vlasov4d)

        call transposevx(vlasov4d)

        !compute vlasov4d%bz=B^{n+1/2} from Ex^n, Ey^n, B^{n-1/2}  !!!!Attention initialisation B^{-1/2}
        vlasov4d%bzn=vlasov4d%bz
        call solve_faraday(vlasov4d,maxwell,dt)  
        
        !compute vlasov4d%bzn=B^n=0.5(B^{n+1/2}+B^{n-1/2})          
        vlasov4d%bzn=0.5_8*(vlasov4d%bz+vlasov4d%bzn)
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
        
        !compute (vlasov4d%ex,vlasov4d%ey)=E^{n+1/2} from vlasov4d%bzn=B^n
        call solve_ampere(vlasov4d,maxwell,0.5_8*dt) 

        if (va==3) then 
           vlasov4d%jx3=vlasov4d%jx
           vlasov4d%jy3=vlasov4d%jy
        endif
     endif

     !advec x + compute this%jx1
     call advection_x1(vlasov4d,0.5_8*dt)

     !advec y + compute this%jy1
     call advection_x2(vlasov4d,0.5_8*dt)


     if (va==1) then 
        call transposexv(vlasov4d)
        !compute rho^{n+1}
        call compute_charge(vlasov4d)
        call transposevx(vlasov4d)
        !compute E^{n+1} via Poisson
        call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
        print *,'nrj-1',nrj
     endif

     call transposexv(vlasov4d)

     call advection_x3x4(vlasov4d,dt)

     call transposevx(vlasov4d)

     !copy jy^{**}
     vlasov4d%jy=vlasov4d%jy1
     !advec y + compute this%jy1
     call advection_x2(vlasov4d,0.5*dt)
     
     !copy jx^*
     vlasov4d%jx=vlasov4d%jx1
     !advec x + compute this%jx1
     call advection_x1(vlasov4d,0.5_f64*dt)

     if (va==0) then 
        !compute the good jy current
        vlasov4d%jy=0.5_f64*(vlasov4d%jy+vlasov4d%jy1)
        !compute the good jx current
        vlasov4d%jx=0.5_f64*(vlasov4d%jx+vlasov4d%jx1)

        !compute E^{n+1} from B^{n+1/2}, vlasov4d%jx, vlasov4d%jy, E^n
        vlasov4d%ex=vlasov4d%exn
        vlasov4d%ey=vlasov4d%eyn
        vlasov4d%bzn=vlasov4d%bz     
        
        call solve_ampere(vlasov4d, maxwell, dt) 

        !copy ex and ey at t^n for the next loop
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
     endif

     if (va==3) then 
        !f --> ft, current (this%jx,this%jy), ft-->f
        call transposexv(vlasov4d)
        !compute this%jx, this%jy (zero average) at time tn
        call compute_current(vlasov4d)

        call transposevx(vlasov4d)

        !compute J^{n+1/2}=0.5*(J^n+J^{n+1})
        vlasov4d%jy=0.5_f64*(vlasov4d%jy+vlasov4d%jy3)
        vlasov4d%jx=0.5_f64*(vlasov4d%jx+vlasov4d%jx3)
        
        !compute E^{n+1} from B^{n+1/2}, vlasov4d%jx, vlasov4d%jy, E^n
        vlasov4d%ex=vlasov4d%exn
        vlasov4d%ey=vlasov4d%eyn
        vlasov4d%bzn=vlasov4d%bz     
        
        call solve_ampere(vlasov4d, maxwell, dt) 

        !copy ex and ey at t^n for the next loop
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
     endif


     if (va==0) then 
        call transposexv(vlasov4d)
        !compute rho^{n+1}
        call compute_charge(vlasov4d)
        call transposevx(vlasov4d)
        !compute E^{n+1} via Poisson
        call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
        print *,'nrj-2',nrj        

        print *,'verif charge conservation',maxval(vlasov4d%exn-vlasov4d%ex),maxval(vlasov4d%eyn-vlasov4d%ey)
        print *,'verif charge conservation',minval(vlasov4d%exn-vlasov4d%ex),minval(vlasov4d%eyn-vlasov4d%ey)
        print *,'mass',maxval(vlasov4d%exn),maxval(vlasov4d%ex),(sum(vlasov4d%rho(1:vlasov4d%geomx%nx,1:vlasov4d%geomx%ny))*vlasov4d%geomv%dx*vlasov4d%geomv%dy-mass0)/mass0
        print *,' ',minval(vlasov4d%eyn),minval(vlasov4d%ey),maxval(vlasov4d%bz),maxval(vlasov4d%bzn)
!        vlasov4d%ex=vlasov4d%exn
!        vlasov4d%ey=vlasov4d%eyn
     endif

     if (mod(iter,fthdiag).eq.0) then 
        nrj=sum(vlasov4d%ex*vlasov4d%ex+vlasov4d%ey*vlasov4d%ey) &
           *(vlasov4d%geomx%dx)*(vlasov4d%geomx%dy)
        nrj=0.5_wp*log(nrj)
        print *,'nrj-3',nrj
        call thdiag(vlasov4d,nrj,iter*dt)
     endif

  end do

  tcpu2 = MPI_WTIME()
  if (prank == MPI_MASTER) then
       write(*,"(//10x,' Wall time = ', G15.3, ' sec' )") (tcpu2-tcpu1)*psize
  end if

  call delete(poisson)
  call sll_halt_collective()

  print*,'PASSED'

!####################################################################################

contains

!####################################################################################

  subroutine initlocal()
    
    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,k,l,error
    sll_real64 :: xi, eps, kx, ky
    sll_int32  :: gi, gj, gk, gl
    sll_int32, dimension(4) :: global_indices
    sll_int32 :: psize

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    call spl_x3x4%initialize(geomv%nx, geomv%ny,                        &
    &                        geomv%x0, geomv%x1, geomv%y0, geomv%y1,    &
    &                        SLL_PERIODIC, SLL_PERIODIC)

    call new(vlasov4d,geomx,geomv,spl_x3x4,error)

    call compute_local_sizes_4d(vlasov4d%layout_x,loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

    xi  = 0.90_f64
    eps = 0.05_f64
    kx  = 2_f64*sll_pi/(geomx%nx*geomx%dx)
    ky  = 2_f64*sll_pi/(geomx%ny*geomx%dy)

    do l=1,loc_sz_l 
       do k=1,loc_sz_k
          do j=1,loc_sz_j
             do i=1,loc_sz_i
                
                global_indices = local_to_global_4D(vlasov4d%layout_x,(/i,j,k,l/)) 
                gi = global_indices(1)
                gj = global_indices(2)
                gk = global_indices(3)
                gl = global_indices(4)
                
                x  = geomx%x0+(gi-1)*geomx%dx
                y  = geomx%y0+(gj-1)*geomx%dy
                vx = geomv%x0+(gk-1)*geomv%dx
                vy = geomv%y0+(gl-1)*geomv%dy
                
                v2 = vx*vx+vy*vy
                vlasov4d%f(i,j,k,l)=(1._f64+eps*cos(kx*x))*1._f64/(2._f64*sll_pi)*exp(-0.5_f64*v2)
                
             end do
          end do
       end do
    end do

    print *,'init'


    call initialize(maxwell, geomx%x0, geomx%x1, geomx%nx, &
                    geomx%y0, geomx%y1, geomx%ny, TE_POLARIZATION)

    call initialize(poisson, geomx%x0, geomx%x1, geomx%nx, &
                    geomx%y0, geomx%y1, geomx%ny, error)

  end subroutine initlocal

  subroutine solve_ampere(vlasov4d, maxwell2d, dt)
    type(vlasov4d_spectral)   :: vlasov4d 
    type(maxwell_2d_pstd)     :: maxwell2d
    sll_real64, intent(in)    :: dt
    
    call ampere(maxwell2d, vlasov4d%ex, vlasov4d%ey, vlasov4d%bzn, dt, vlasov4d%jx, vlasov4d%jy)

  end subroutine solve_ampere
  
  subroutine solve_faraday(vlasov4d, maxwell2d, dt)
    type(vlasov4d_spectral)   :: vlasov4d 
    type(maxwell_2d_pstd)     :: maxwell2d
    sll_real64, intent(in)    :: dt
    
    call faraday(maxwell2d, vlasov4d%exn, vlasov4d%eyn, vlasov4d%bz, dt)
    
  end subroutine solve_faraday
  
end program vm4d_spectral
