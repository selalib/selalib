program vm4d_spectral_charge

#define MPI_MASTER 0
#include "selalib-mpi.h"

  use sll_vlasov4d_base
  use sll_vlasov4d_spectral_charge

  implicit none

  type(maxwell_2d_pstd)          :: maxwell
  type(poisson_2d_periodic)      :: poisson 
  type(vlasov4d_spectral_charge) :: vlasov4d 

  type(cubic_spline_2d_interpolator), target :: spl_x3x4

  sll_int32  :: iter 
  sll_real64 :: nrj, tcpu1, tcpu2, mass0

  sll_int32  :: prank, comm
  sll_int64  :: psize

  sll_int32  :: loc_sz_i, loc_sz_j, loc_sz_k, loc_sz_l
!  sll_int32  :: va

!  va=0 !(va=0 --> Valis ; va=1 --> Vlasov-Poisson ; va=2 --> diag charge ; va=3 --> classic algorithm)

  call sll_boot_collective()
  prank = sll_get_collective_rank(sll_world_collective)
  psize = sll_get_collective_size(sll_world_collective)
  comm  = sll_world_collective%comm

  tcpu1 = MPI_WTIME()
  if (prank == MPI_MASTER) then
     print*,'MPI Version of slv2d running on ',psize, ' processors'
  end if

  call initlocal()

  !f --> ft
  call transposexv(vlasov4d)

  call compute_charge(vlasov4d)
  call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
  vlasov4d%exn=vlasov4d%ex
  vlasov4d%eyn=vlasov4d%ey

  !ft --> f
  call transposevx(vlasov4d)

  mass0=sum(vlasov4d%rho(1:vlasov4d%nc_eta1,1:vlasov4d%nc_eta2))*vlasov4d%delta_eta1*vlasov4d%delta_eta2
  print *,'mass init',mass0

  !###############
  !TIME LOOP
  !###############
  do iter=1,vlasov4d%nbiter

     !print *,'iter',iter
     if (iter ==1 .or. mod(iter,vlasov4d%fdiag) == 0) then 
        call write_xmf_file(vlasov4d,iter/vlasov4d%fdiag)
     end if


     if ((vlasov4d%va==0).or.(vlasov4d%va==3)) then 
        !f --> ft, current (this%jx,this%jy), ft-->f
        call transposexv(vlasov4d)
        !compute this%jx, this%jy (zero average) at time tn
        call compute_current(vlasov4d)

        call transposevx(vlasov4d)

        !compute vlasov4d%bz=B^{n+1/2} from Ex^n, Ey^n, B^{n-1/2}  !!!!Attention initialisation B^{-1/2}
        vlasov4d%bzn=vlasov4d%bz
        call solve_faraday(vlasov4d,maxwell,vlasov4d%dt)  
        
        !compute vlasov4d%bzn=B^n=0.5(B^{n+1/2}+B^{n-1/2})          
        vlasov4d%bzn=0.5_f64*(vlasov4d%bz+vlasov4d%bzn)
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
        
        !compute (vlasov4d%ex,vlasov4d%ey)=E^{n+1/2} from vlasov4d%bzn=B^n
        call solve_ampere(vlasov4d,maxwell,0.5_f64*vlasov4d%dt) 

        if (vlasov4d%va==3) then 
           vlasov4d%jx3=vlasov4d%jx
           vlasov4d%jy3=vlasov4d%jy
        endif
     endif

     !advec x + compute this%jx1
     call advection_x1(vlasov4d,0.5_f64*vlasov4d%dt)

     !advec y + compute this%jy1
     call advection_x2(vlasov4d,0.5_f64*vlasov4d%dt)


     if (vlasov4d%va==1) then 
        call transposexv(vlasov4d)
        !compute rho^{n+1}
        call compute_charge(vlasov4d)
        call transposevx(vlasov4d)
        !compute E^{n+1} via Poisson
        call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
     endif

     call transposexv(vlasov4d)

     call advection_x3x4(vlasov4d,vlasov4d%dt)

     call transposevx(vlasov4d)

     !copy jy^{**}
     vlasov4d%jy=vlasov4d%jy1
     !advec y + compute this%jy1
     call advection_x2(vlasov4d,0.5_f64*vlasov4d%dt)
     
     !copy jx^*
     vlasov4d%jx=vlasov4d%jx1
     !advec x + compute this%jx1
     call advection_x1(vlasov4d,0.5_f64*vlasov4d%dt)

     if (vlasov4d%va==0) then 
        !compute the good jy current
        vlasov4d%jy=0.5_f64*(vlasov4d%jy+vlasov4d%jy1)
        !compute the good jx current
        vlasov4d%jx=0.5_f64*(vlasov4d%jx+vlasov4d%jx1)
        print *,'sum jx jy',sum(vlasov4d%jx)*vlasov4d%delta_eta1*vlasov4d%delta_eta2,&
             sum(vlasov4d%jy)*vlasov4d%delta_eta1*vlasov4d%delta_eta2, &
             maxval(vlasov4d%jx),maxval(vlasov4d%jy)

        !compute E^{n+1} from B^{n+1/2}, vlasov4d%jx, vlasov4d%jy, E^n
        vlasov4d%ex=vlasov4d%exn
        vlasov4d%ey=vlasov4d%eyn
        vlasov4d%bzn=vlasov4d%bz     
        
        call solve_ampere(vlasov4d, maxwell, vlasov4d%dt) 

        !copy ex and ey at t^n for the next loop
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
     endif

     if (vlasov4d%va==3) then 
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
        
        call solve_ampere(vlasov4d, maxwell, vlasov4d%dt) 

        !copy ex and ey at t^n for the next loop
        vlasov4d%exn=vlasov4d%ex
        vlasov4d%eyn=vlasov4d%ey
     endif


     if (vlasov4d%va==0) then 
        call transposexv(vlasov4d)
        !compute rho^{n+1}
        call compute_charge(vlasov4d)
        call transposevx(vlasov4d)

        !compute E^{n+1} via Poisson
        call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
!        print *,'verif charge conservation',maxval(vlasov4d%exn-vlasov4d%ex), &
!             maxval(vlasov4d%eyn-vlasov4d%ey)
     endif
     
     if (vlasov4d%va==1) then 
        !recompute the electric field at time (n+1) for diagnostics
        call transposexv(vlasov4d)
        !compute rho^{n+1}
        call compute_charge(vlasov4d)
        call transposevx(vlasov4d)
        !compute E^{n+1} via Poisson
        call solve(poisson,vlasov4d%ex,vlasov4d%ey,vlasov4d%rho,nrj)
     endif


     if (mod(iter,vlasov4d%fthdiag).eq.0) then 
        call write_energy(vlasov4d, iter*vlasov4d%dt)
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

    use init_functions
    
    sll_real64 :: vx,vy,v2,x,y
    sll_int32  :: i,j,k,l,error
    sll_real64 :: eps, kx, ky
    sll_int32  :: gi, gj, gk, gl
    sll_int32, dimension(4) :: global_indices
    sll_int32 :: psize

    prank = sll_get_collective_rank(sll_world_collective)
    psize = sll_get_collective_size(sll_world_collective)
    comm  = sll_world_collective%comm

    call read_input_file(vlasov4d)

    call spl_x3x4%initialize(vlasov4d%nc_eta3, vlasov4d%nc_eta4,   & 
                             vlasov4d%eta3_min, vlasov4d%eta3_max, &
                             vlasov4d%eta4_min, vlasov4d%eta4_max, &
    &                        SLL_PERIODIC, SLL_PERIODIC)


    call initialize(vlasov4d,spl_x3x4,error)

    call compute_local_sizes_4d(vlasov4d%layout_x, &
                                loc_sz_i,loc_sz_j,loc_sz_k,loc_sz_l)        

    eps = 0.05_f64
    kx  = 2_f64*sll_pi/(vlasov4d%nc_eta1*vlasov4d%delta_eta1)
    ky  = 2_f64*sll_pi/(vlasov4d%nc_eta2*vlasov4d%delta_eta2)

    do l=1,loc_sz_l 
       do k=1,loc_sz_k
          do j=1,loc_sz_j
             do i=1,loc_sz_i
                
                global_indices = local_to_global_4D(vlasov4d%layout_x,(/i,j,k,l/)) 
                gi = global_indices(1)
                gj = global_indices(2)
                gk = global_indices(3)
                gl = global_indices(4)
                
                x  = vlasov4d%eta1_min+(gi-1)*vlasov4d%delta_eta1
                y  = vlasov4d%eta2_min+(gj-1)*vlasov4d%delta_eta2
                vx = vlasov4d%eta3_min+(gk-1)*vlasov4d%delta_eta3
                vy = vlasov4d%eta4_min+(gl-1)*vlasov4d%delta_eta4
                
                v2 = vx*vx+vy*vy

                select case(vlasov4d%num_case)
                case(1)
                    vlasov4d%f(i,j,k,l)= landau_1d(eps, kx, x, v2)
                case(2)
                    vlasov4d%f(i,j,k,l)= landau_1d(eps, ky, y, v2)
                case(3)
                    vlasov4d%f(i,j,k,l)= landau_cos_prod(eps, kx, ky, x, y, v2)
                case(4)
                    vlasov4d%f(i,j,k,l)= landau_cos_sum(eps, kx, ky, x, y, v2)
                case(5)
                    vlasov4d%f(i,j,k,l)= tsi(eps, kx, x, vx, v2)
                end select
                
             end do
          end do
       end do
    end do
    
    print *,'init'
    
    
    call initialize(maxwell, &
         vlasov4d%eta1_min, vlasov4d%eta1_max, vlasov4d%nc_eta1, &
         vlasov4d%eta2_min, vlasov4d%eta2_max, vlasov4d%nc_eta2, TE_POLARIZATION)
    
    call initialize(poisson, &
         vlasov4d%eta1_min, vlasov4d%eta1_max, vlasov4d%nc_eta1, &
         vlasov4d%eta2_min, vlasov4d%eta2_max, vlasov4d%nc_eta2, error)
    
    
  end subroutine initlocal
  
  subroutine solve_ampere(vlasov4d, maxwell2d, dt)
    type(vlasov4d_spectral_charge)   :: vlasov4d 
    type(maxwell_2d_pstd)     :: maxwell2d
    sll_real64, intent(in)    :: dt
    
    call ampere(maxwell2d, vlasov4d%ex, vlasov4d%ey, &
                vlasov4d%bzn, dt, vlasov4d%jx, vlasov4d%jy)

  end subroutine solve_ampere
  
  subroutine solve_faraday(vlasov4d, maxwell2d, dt)
    type(vlasov4d_spectral_charge)   :: vlasov4d 
    type(maxwell_2d_pstd)     :: maxwell2d
    sll_real64, intent(in)    :: dt
    
    call faraday(maxwell2d, vlasov4d%exn, vlasov4d%eyn, vlasov4d%bz, dt)
    
  end subroutine solve_faraday
  
end program vm4d_spectral_charge
