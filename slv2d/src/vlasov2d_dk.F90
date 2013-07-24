module vlasov2d_dk_module

#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_poisson_solvers.h"
 use used_precision
 use splinenn_class
 use splinepp_class
 use geometry_module
 use diagnostiques_module
 use polar_operators
 use polar_advection
 use sll_cubic_splines
 use sll_cubic_spline_interpolator_1d

 !use clock

 implicit none
 private
 public :: new, dealloc,advection_x, advection_v,&
           densite_charge,densite_charge_dk,transposexv,transposevx,thdiag,densite_courant,&
           compute_profile,solve_quasi_neutral,compute_field_dk,advection_x_dk,&
           advection_x3_dk,advection_x4_dk

 type, public :: vlasov2d
   sll_real64, dimension(:,:,:,:), pointer :: ft
   type (splinepp) :: splinex ! spline periodique pour X
   type (splinenn) :: splinev ! spline naturel pour V
   type (geometry) :: geomx, geomv
   logical :: transposed       ! permet de definir si f ou ft est derniere 
                               ! fonction de distribution mise a jour
   sll_int32 :: jstartx, jendx
   sll_int32 :: jstartv, jendv  ! definition de la bande de calcul
   class(sll_interpolator_1d_base), pointer :: interp_x3
   class(sll_interpolator_1d_base), pointer :: interp_x4
 end type vlasov2d

 ! variables globales 
 sll_real64, dimension(:,:),allocatable :: P_x, P_y
 type(cubic_spline_1d_interpolator), target :: spl_x3
 type(cubic_spline_1d_interpolator), target :: spl_x4

 interface new
   module procedure new_vlasov2d
 end interface

 interface dealloc
   module procedure dealloc_vlasov2d
 end interface

contains

 subroutine new_vlasov2d(this,geomx,geomv,error, jstartx, jendx, &
                        jstartv, jendv)

  type(vlasov2d),intent(out)      :: this
  type(geometry),intent(in)       :: geomx
  type(geometry),intent(in)       :: geomv
  sll_int32, intent(out)          :: error
  sll_int32, intent(in), optional ::  jstartx
  sll_int32, intent(in), optional ::  jendx
  sll_int32, intent(in), optional ::  jstartv
  sll_int32, intent(in), optional ::  jendv

  error = 0

  ! on commence par utiliser la fonction f(x,y,vx,vy)
  this%transposed=.false.
  ! definition des bandes de calcul (en n'oubliant pas le cas sequentiel)
  if (.not.(present(jstartx))) then
     this%jstartx = 1
  else
     this%jstartx = jstartx
  end if
  if (.not.(present(jendx))) then
     this%jendx = geomx%ny
  else
     this%jendx = jendx
  end if
  if (.not.(present(jstartv))) then 
      this%jstartv = 1
  else
     this%jstartv = jstartv
  end if
  if (.not.(present(jendv))) then
     this%jendv = geomv%nx
  else
     this%jendv = jendv
  end if
  
  ! initialisation de la geometrie
  this%geomx=geomx
  this%geomv=geomv
  ! initialisation des splines de l'espace des vitesses
  call new(this%splinev,geomv,error)
  ! initialisation des splines de l'espace physique
  call new(this%splinex,geomx,error)  

  call spl_x3%initialize( geomv%nx+1, geomv%x0, geomv%x1, SLL_PERIODIC)
  call spl_x4%initialize( geomv%ny+1, geomv%y0, geomv%y1, SLL_PERIODIC)

  this%interp_x3 => spl_x3
  this%interp_x4 => spl_x4


  ! allocation memoire
  SLL_ALLOCATE(this%ft(geomv%nx,geomv%ny,geomx%nx,this%jstartx:this%jendx),error)

  SLL_ALLOCATE(P_x(this%geomv%nx,this%geomv%ny),error)
  SLL_ALLOCATE(P_y(this%geomv%nx,this%geomv%ny),error)

 end subroutine new_vlasov2d

 subroutine dealloc_vlasov2d(this)
  type(vlasov2d),intent(out) :: this
  sll_int32 :: error
  SLL_DEALLOCATE(this%ft,error)
 end subroutine dealloc_vlasov2d

 !>
 !> fait une advection en x sur un pas de temps dt
 !>
 subroutine advection_x(this,f,dt)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_real64, intent(in) :: dt

  sll_real64 :: depx, depy   ! deplacement par rapport au maillage
  sll_real64 :: vx, vy       ! vitesse du point courant
  sll_int32 :: iv, jv ! indices de boucle

  ! verifier que la transposition est a jour
  if (this%transposed) stop 'advection_x: on travaille sur f et pas ft'
  do jv=this%jstartv,this%jendv
     vy = this%geomv%y0+(jv-1)*this%geomv%dy
     depy = vy*dt
     do iv=1,this%geomv%nx
        vx = this%geomv%x0+(iv-1)*this%geomv%dx
        depx = vx*dt
        call interpole(this%splinex,f(:,:,iv,jv),depx,depy,jv==0)
     end do
  end do

 end subroutine advection_x


 !>
 !> fait une advection en x sur un pas de temps dt
 !>
 subroutine advection_x_dk(this,plan_adv,f,adv_field)
  type(vlasov2d),intent(inout) :: this
  type(sll_plan_adv_polar),pointer        :: plan_adv
  sll_real64, dimension(:,:,:,:) :: adv_field
  sll_real64, dimension(:,:,:,this%jstartv:) :: f
  sll_int32 :: iv, jv ! indices de boucle
  sll_int32 :: ierr
  sll_real64, dimension(:,:),allocatable :: fn
  sll_real64, dimension(:,:),allocatable :: fnp1
  
  SLL_ALLOCATE(fn(this%geomx%nx,this%geomx%ny+1),ierr)
  SLL_ALLOCATE(fnp1(this%geomx%nx,this%geomx%ny+1),ierr)
  ! verifier que la transposition est a jour
  if (this%transposed) stop 'advection_x: on travaille sur f et pas ft'
  do jv=this%jstartv,this%jendv
     do iv=1,this%geomv%nx
        plan_adv%field = adv_field(1:2,:,:,iv)
        fn(:,1:this%geomx%ny) = f(:,:,iv,jv)
        fn(:,this%geomx%ny+1) = fn(:,1)
        call advect_CG_polar(plan_adv,fn,fnp1)
        f(:,1:this%geomx%ny,iv,jv)=fnp1
        !call interpole(this%splinex,f(:,:,iv,jv),depx,depy,jv==0)
     end do
  end do
  
  SLL_DEALLOCATE_ARRAY(fn,ierr)
  SLL_DEALLOCATE_ARRAY(fnp1,ierr)
  
 end subroutine advection_x_dk



 subroutine advection_x3_dk(this,dt)

  type(vlasov2d),intent(inout) :: this

  !sll_real64, dimension(:,:,:), intent(in) :: ex
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x3, nc_x4,i,j,l
  sll_real64 :: alpha
  sll_real64,dimension(:),allocatable::buf
  sll_int32::ierr

  SLL_ASSERT(this%transposed) 
  

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny
  
  SLL_ALLOCATE(buf(nc_x3+1),ierr)
  
  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        do l=1,nc_x4
           alpha = (this%geomv%y0+real(l-1,f64)*this%geomv%dy)*dt
           buf(1:nc_x3)=this%ft(:,l,i,j)
           buf(nc_x3+1)=buf(1)
           buf = this%interp_x3%interpolate_array_disp( &
                              nc_x3+1, buf, alpha )
           this%ft(:,l,i,j)=buf(1:nc_x3)                   
       end do
    end do
 end do
  SLL_DEALLOCATE_ARRAY(buf,ierr)

 end subroutine advection_x3_dk


 subroutine advection_x4_dk(this,ey,dt)

  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:), intent(in) :: ey
  sll_real64, intent(in) :: dt
  sll_int32  :: nc_x1, nc_x3, nc_x4,i,j,k
  sll_real64 :: alpha
  sll_real64,dimension(:),allocatable::buf
  sll_int32::ierr

  SLL_ASSERT(this%transposed) 

  nc_x1    = this%geomx%nx
  nc_x3    = this%geomv%nx
  nc_x4    = this%geomv%ny

  SLL_ALLOCATE(buf(nc_x4+1),ierr)

  do j=this%jstartx,this%jendx
     do i=1,nc_x1
        do k=1,nc_x3
           alpha = ey(i,j,k)*dt
           buf(1:nc_x4)=this%ft(k,:,i,j)
           buf(nc_x4+1)=buf(1)
           buf = this%interp_x4%interpolate_array_disp( &
                              nc_x4+1, buf, alpha )
           this%ft(k,:,i,j)=buf(1:nc_x4)                   
       end do
    end do
 end do
  SLL_DEALLOCATE_ARRAY(buf,ierr)


 end subroutine advection_x4_dk







 !>
 !> fait une advection en v sur un pas de temps dt
 !>
 subroutine advection_v(this,fx,fy,dt,bz)

  type(vlasov2d),intent(inout) :: this

  sll_real64, dimension(:,:), intent(in) :: fx, fy
  sll_real64, dimension(:,:), optional, intent(in) :: bz
  sll_real64, intent(in) :: dt

  sll_real64 :: depvx, depvy   ! deplacement par rapport au maillage
  sll_int32 :: i, j ,iv, jv, im1! indices de boucle
  sll_real64 :: ctheta, stheta, px, py

  ! verifier que la transposition est a jour
  if (.not.(this%transposed)) stop 'advection_v: on travaille sur ft et pas f'

  if (present(bz)) then

   do j=this%jstartx,this%jendx
    do i=1,this%geomx%nx
     ctheta = cos(bz(i,j)*dt)
     stheta = sin(bz(i,j)*dt)
     depvx  = -0.5*dt*fx(i,j)
     depvy  = -0.5*dt*fy(i,j)
     do jv=1,this%geomv%ny-1
      py = this%geomv%y0+(jv-1)*this%geomv%dy
      do iv=1,this%geomv%nx-1 
       px = this%geomv%x0+(iv-1)*this%geomv%dx
       P_x(iv,jv) = depvx+(px+depvx)*ctheta-(py+depvy)*stheta
       P_y(iv,jv) = depvy+(px+depvx)*stheta+(py+depvy)*ctheta
      end do
     end do

     call interpole(this%splinev,this%ft(:,:,i,j),this%ft(:,:,i,j),P_x,P_y)

    end do
   end do

  else

   do j=this%jstartx,this%jendx
    do i=1,this%geomx%nx
     im1=mod(i-1+this%geomx%nx,this%geomx%nx)
     depvx = fx(i,j)*dt
     depvy = fy(i,j)*dt
     !depvx =  fx(i,j)*dt;depvy=0._wp
     call interpole(this%splinev,this%ft(:,:,i,j),depvx,depvy,j==0)
     !call interpole(this%splinev,this%ft(:,:,i,j),depvx,depvy,(j .eq. 3) .and. (i .eq. 3))
    end do
   end do

  end if
       
 end subroutine advection_v

subroutine solve_quasi_neutral(this,plan_poisson,rho,phi)
  type(vlasov2d),intent(inout) :: this
  type(sll_plan_poisson_polar), pointer :: plan_poisson
  sll_real64, dimension(:,:,:), intent(in)  :: rho
  sll_real64, dimension(:,:,:), intent(out)  :: phi
  sll_int32 :: i,ierr

  do i = 1,this%geomv%nx
    call solve_poisson_polar(plan_poisson,rho(:,:,i),phi(:,:,i))
  enddo


end subroutine solve_quasi_neutral

!call compute_grad_field(plan%grad,plan%phi,plan%adv%field)

subroutine compute_field_dk(this,grad,phi,adv_field)
  type(vlasov2d),intent(inout) :: this
  sll_real64, dimension(:,:,:,:), intent(out)  :: adv_field
  sll_real64, dimension(:,:,:), intent(inout)  :: phi
  type(plan_polar_op), pointer :: grad
  sll_int32 :: i,ierr
  
  do i=1,this%geomv%nx
    call compute_grad_field(grad,phi(:,:,i),adv_field(1:2,:,:,i))
  enddo
  
  !for the moment finite differences
  do i=1,this%geomv%nx
    adv_field(3,:,:,i)=(phi(:,:,modulo(i+1-1+this%geomv%nx,this%geomv%nx)+1)&
      &-phi(:,:,modulo(i-1-1+this%geomv%nx,this%geomv%nx)+1))/(2._f64*this%geomv%dx)
  enddo
  
end subroutine compute_field_dk


!>------------------------------------------------
!> calcule la densite de charge rho a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
!> Poisson n'est pas parallele on transmet donc rho
!> a tous les processeurs
subroutine densite_charge_dk(this, rho)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:,:), intent(out)  :: rho
   sll_real64, dimension(this%geomx%nx,this%geomx%ny+1,this%geomv%nx+1) :: locrho
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   SLL_ASSERT(this%transposed)
   
   rho(:,:,:) = 0._f64
   locrho(:,:,:) = 0._f64
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do iv=1,this%geomv%nx!-1 
            do jv=1,this%geomv%ny!-1
               locrho(i,j,iv) = locrho(i,j,iv) + this%geomv%dy* &
                    this%ft(iv,jv,i,j) 
            end do
         end do
      end do
   end do
   call mpi_barrier(comm,error)
   c=this%geomx%nx*(this%geomx%ny+1)*(this%geomv%nx+1)
   call mpi_allreduce(locrho,rho,c,MPI_REAL8,MPI_SUM,comm,error)
   rho(:,:,this%geomv%nx+1)=rho(:,:,1)
   rho(:,this%geomx%ny+1,:)=rho(:,1,:)
   
end subroutine densite_charge_dk


subroutine densite_charge(this, rho)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: rho
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locrho
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   SLL_ASSERT(this%transposed)
   
   rho(:,:) = 0.
   locrho(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            do iv=1,this%geomv%nx-1 
               locrho(i,j) = locrho(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) 
            end do
         end do
      end do
   end do
   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locrho,rho,c,MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_charge


!>------------------------------------------------
!> calcule la densite de courant jx et jy a partir de ft
!> en fait le moment d'ordre 0 de f. Les constantes
!> ne sont pas introduites
!>------------------------------------------------
subroutine densite_courant(this, jx, jy)

   type(vlasov2d),intent(inout) :: this
   sll_int32 :: error
   sll_real64, dimension(:,:), intent(out)  :: jx, jy
   sll_real64 :: vx, vy       ! vitesse du point courant
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjx
   sll_real64, dimension(this%geomx%nx,this%geomx%ny) :: locjy
   sll_int32 :: i,j,iv,jv,c
   sll_int32 :: comm

   comm   = sll_world_collective%comm

   ! verifier que la transposition est a jour
   if (.not.(this%transposed)) &
      stop 'densite_courant: on travaille sur ft et pas f'
   !    rho(:,this%jstartx:this%jendx)=0.
   jx(:,:) = 0.; jy(:,:) = 0.
   locjx(:,:) = 0.; locjy(:,:) = 0.
   do j=this%jstartx,this%jendx
      do i=1,this%geomx%nx
         do jv=1,this%geomv%ny-1
            vy = this%geomv%y0+(jv-1)*this%geomv%dy
            do iv=1,this%geomv%nx-1 
               vx = this%geomv%x0+(iv-1)*this%geomv%dx
               locjx(i,j) = locjx(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vx
               locjy(i,j) = locjy(i,j) + this%geomv%dx*this%geomv%dy* &
                    this%ft(iv,jv,i,j) * vy
            end do
         end do
      end do
   end do

   call mpi_barrier(comm,error)
   c=this%geomx%nx*this%geomx%ny
   call mpi_allreduce(locjx,jx,c, MPI_REAL8,MPI_SUM,comm,error)
   call mpi_allreduce(locjy,jy,c, MPI_REAL8,MPI_SUM,comm,error)

end subroutine densite_courant

!>---------------------------------------------
!> transpose la fonction de distribution
!> ATTENTION: cet fonction fait intervenir des 
!> communications entre les processeurs.
!>---------------------------------------------
subroutine transposexv(this,f)
   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(in) :: f
   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(f, this%ft, sizexy, sizevxvy, num_threads)

this%transposed=.true.

end subroutine transposexv

subroutine transposevx(this,f)

   type(vlasov2d),intent(inout) :: this
   sll_real64, dimension(:,:,:,:),intent(out) :: f

   sll_int32 :: sizexy, sizevxvy
   sll_int32 :: my_num, num_threads

   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   sizexy = this%geomx%nx * this%geomx%ny
   sizevxvy = this%geomv%nx * this%geomv%ny
   call transpose(this%ft,f, sizevxvy, sizexy, num_threads)

   this%transposed=.false.

end subroutine transposevx

subroutine compute_profile(this,prof,geom,rpeak,deltar,kappa,R0)!,n0_case)
  type(vlasov2d),intent(in) :: this
  type(geometry),intent(in) :: geom
  sll_real64,dimension(:),intent(out)   :: prof
  sll_real64,intent(in) :: rpeak,deltar,kappa,R0
  !sll_int32,intent(in)  :: n0_case


    sll_real64 :: r, dr_loc, rth, tmp
    sll_int32     :: ir, Nr_loc
    sll_real64 :: profnorm,profnorm_tmp,invL
    
    !*** local variable initialization ***
    Nr_loc = geom%Nx
    dr_loc = geom%dx
    
    invL = kappa/R0
    
    !*** compute prof solution of :                            ***
    !***  2/(prof(r)+prof(r-1))*(prof(r)-prof(r-1))/dr               ***
    !***                  = -invL*cosh^-2(r-rpeak/deltar) ***
    prof(1) = 10._f64**19
    do ir=1,Nr_loc-1
      r            = geom%xgrid(ir)
      rth          = r + dr_loc*0.5_f64
      tmp          = -invL / &
        cosh((rth-rpeak)/deltar)**2
      tmp          = 0.5_f64*dr_loc*tmp
      prof(ir+1) = (1._f64+tmp)/(1._f64-tmp)*prof(ir)
    enddo
    
    !*** normalisation of the density at int(prof(r)rdr)/int(rdr) ***
    ! -> computation of int(prof(r)rdr)
    profnorm_tmp = 0._f64
    do ir = 2,Nr_loc-1
      profnorm_tmp = profnorm_tmp + prof(ir)*geom%xgrid(ir)
    enddo
    profnorm_tmp = profnorm_tmp + 0.5_f64* &
      (prof(1)*geom%xgrid(1) + prof(Nr_loc)*geom%xgrid(Nr_loc))
    ! -> division by int(rdr)
    profnorm_tmp = profnorm_tmp*2._f64*dr_loc/ & 
      (geom%xgrid(Nr_loc)**2-geom%xgrid(1)**2)
      
    profnorm       = profnorm_tmp
    prof(1:Nr_loc) = prof(1:Nr_loc)/profnorm
    
    do ir = 1,Nr_loc
      r            = geom%xgrid(ir)
      rth          = r + dr_loc*0.5_f64    
     ! print *,'prof',ir,prof(ir)-1._f64,-invL /cosh((rth-rpeak)/deltar)**2
    enddo
    
  
end subroutine compute_profile

subroutine thdiag(this,f,phi,t,jstartv)

   type(vlasov2d),intent(inout) :: this
   sll_int32, intent(in) :: jstartv
   sll_real64, dimension(:,:,:,jstartv:),intent(in) :: f
   sll_real64, dimension(:,:,:),intent(in) :: phi
   sll_int32 :: error
   sll_real64, intent(in) :: t  ! current time
   ! variables locales
   sll_int32 :: i,iv, j,jv
   sll_real64 :: x, vx, y, vy,nrj
   !sll_real64,dimension(7) :: diagloc
   sll_real64,dimension(11) :: auxloc
   sll_real64,dimension(13) :: aux
   sll_real64,dimension(0:9) :: diag
   sll_int32 :: my_num, num_threads
   sll_int32 :: comm

   comm   = sll_world_collective%comm
   my_num = sll_get_collective_rank(sll_world_collective)
   num_threads = sll_get_collective_size(sll_world_collective)

   if (my_num == MPI_MASTER) then
      diag=0.
      aux=0.
   end if

   !diagloc = 0._f64
!   auxloc  = 0._f64
!   do i = 1,this%geomx%nx
!      x = this%geomx%x0+(i-1)*this%geomx%dx
!      do j = 1,this%geomx%ny
!         y= this%geomx%y0+(j-1)*this%geomx%dy
!         do iv=1,this%geomv%nx
!            vx = this%geomv%x0+(iv-1)*this%geomv%dx
!            do jv=this%jstartv,this%jendv
!               vy = this%geomv%y0+(jv-1)*this%geomv%dy
!               !diagloc(2) = diagloc(2) + f(i,j,iv,jv)*f(i,j,iv,jv)
!
!               auxloc(1) = auxloc(1) + f(i,j,iv,jv)         ! avg(f)
!               auxloc(2) = auxloc(2) + x*f(i,j,iv,jv)       ! avg(x)
!               auxloc(3) = auxloc(3) + vx*f(i,j,iv,jv)      ! avg(vx)
!               auxloc(4) = auxloc(4) + x*x*f(i,j,iv,jv)     ! avg(x^2)
!               auxloc(5) = auxloc(5) + vx*vx*f(i,j,iv,jv)   ! avg(vx^2)
!               auxloc(6) = auxloc(6) + x*vx*f(i,j,iv,jv)    ! avg(x*vx)
!               auxloc(7) = auxloc(7) + y*f(i,j,iv,jv)       ! avg(y)
!               auxloc(8) = auxloc(8) + vy*f(i,j,iv,jv)      ! avg(vy)
!               auxloc(9) = auxloc(9) + y*y*f(i,j,iv,jv)     ! avg(y^2)
!               auxloc(10) = auxloc(10) + vy*vy*f(i,j,iv,jv) ! avg(vy^2)
!               auxloc(11) = auxloc(11) + y*vy*f(i,j,iv,jv)  ! avg(y*vy)
!
!            end do
!         end do
!      end do
!   end do

   auxloc  = 0._f64
   do jv=this%jstartv,this%jendv
      vy = this%geomv%y0+(jv-1)*this%geomv%dy
      do iv=1,this%geomv%nx
         vx = this%geomv%x0+(iv-1)*this%geomv%dx
         do j = 1,this%geomx%ny
            y = this%geomx%y0+(j-1)*this%geomx%dy
               i=1
               x = this%geomx%x0+(i-1)*this%geomx%dx
               auxloc(1) = auxloc(1) + 0.5_f64*x*f(i,j,iv,jv)             ! mass
               auxloc(2) = auxloc(2) + 0.5_f64*x*abs(f(i,j,iv,jv))        ! L1
               auxloc(3) = auxloc(3) + 0.5_f64*x*abs(f(i,j,iv,jv))**2     ! L2
               auxloc(4) = auxloc(4) + 0.5_f64*x*vy*vy*f(i,j,iv,jv)       ! ekin
               auxloc(5) = auxloc(5) + 0.5_f64*x*phi(i,j,iv)*f(i,j,iv,jv) ! epot
            do i = 2,this%geomx%nx-1
               x = this%geomx%x0+(i-1)*this%geomx%dx
               auxloc(1) = auxloc(1) + x*f(i,j,iv,jv)             ! mass
               auxloc(2) = auxloc(2) + x*abs(f(i,j,iv,jv))        ! L1
               auxloc(3) = auxloc(3) + x*abs(f(i,j,iv,jv))**2     ! L2
               auxloc(4) = auxloc(4) + x*vy*vy*f(i,j,iv,jv)       ! ekin
               auxloc(5) = auxloc(5) + x*phi(i,j,iv)*f(i,j,iv,jv) ! epot
            end do
               i=this%geomx%nx
               x = this%geomx%x0+(i-1)*this%geomx%dx
               auxloc(1) = auxloc(1) + 0.5_f64*x*f(i,j,iv,jv)             ! mass
               auxloc(2) = auxloc(2) + 0.5_f64*x*abs(f(i,j,iv,jv))        ! L1
               auxloc(3) = auxloc(3) + 0.5_f64*x*abs(f(i,j,iv,jv))**2     ! L2
               auxloc(4) = auxloc(4) + 0.5_f64*x*vy*vy*f(i,j,iv,jv)       ! ekin
               auxloc(5) = auxloc(5) + 0.5_f64*x*phi(i,j,iv)*f(i,j,iv,jv) ! epot
         end do
      end do
   end do




   auxloc=auxloc*this%geomx%dx*this%geomx%dy*this%geomv%dx*this%geomv%dy

   call mpi_reduce(auxloc,aux,11,MPI_REAL8,MPI_SUM,0,  &
                   comm, error)
   !call mpi_reduce(diagloc(2),diag(2),1,MPI_REAL8,MPI_SUM,0, &
   !                comm, error)

  nrj = 0._f64
  nrj = sum(phi(this%geomx%nx/2,1:this%geomx%ny,1:this%geomv%nx)**2)*this%geomx%dy*this%geomv%dx

if (my_num==MPI_MASTER) then
   aux(13)=t
   aux(12)=nrj
   !write(*,"('time ', g8.3,' test nrj',f10.5)") t, nrj
   !call time_history("thf","(13(1x,e15.6))",aux(1:13))
   !print "(13(1x,e15.10))",aux(1:13)
   print *,t,nrj,aux(4)+aux(5),aux(1:5)
   !stop
end if

end subroutine thdiag

end module vlasov2d_dk_module


